//
// pomerol2dcore
//
// Copyright (C) 2019 Junya Otsuki
//

#include <pomerol.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <cassert>
#include <memory>
#include <sys/stat.h>
#include <sys/types.h>

#include "ReadWrite.h"
#include "Params.h"
#include "OperatorPresetsExtra.h"

using namespace Pomerol;

typedef std::pair<double, QuantumNumbers> EigenSystem;


void make_dir(std::string &dirname){
    struct stat st;
    if( stat(dirname.c_str(), &st) != 0 ) {
        mkdir(dirname.c_str(), 0755);
    }
}


// Small routine to make fancy screen output for text.
void print_section (const std::string& str)
{
    std::cout << std::string(60, '=') << std::endl;
    std::cout << str << std::endl;
    std::cout << std::string(60, '=') << std::endl;
}


void print_time(clock_t start, const char *str)
{
    std::cout << "#Time: "
              << (double)(clock()-start)/CLOCKS_PER_SEC << " sec"
              << " (" << str << ")"<< std::endl << std::flush;
}


void print_commutation(bool if_commute, const std::string &op, bool verbose)
{
    std::string str_commute(if_commute ? " = 0" : "!= 0");
    if(verbose){
        std::cout << "[H, " << op << "] " << str_commute << std::endl;
    }
}


Lattice::Term* OneBodyTerm ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    auto *T = new Lattice::Term(2);
    bool Operators[2]          = { true, false };
    std::string Labels[2]      = { Label, Label };
    unsigned short Spins[2]    = { spin1, spin2 };
    unsigned short Orbitals[2] = { orbital1, orbital2 };
    T->OperatorSequence.assign(Operators,Operators+2);
    T->SiteLabels.assign(Labels,Labels+2);
    T->Orbitals.assign(Orbitals,Orbitals+2);
    T->Spins.assign(Spins,Spins+2);
    T->Value = Value;
    return T;
}


Lattice::Term* TwoBodyTerm ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short orbital3, unsigned short orbital4, unsigned short spin1, unsigned short spin2, unsigned short spin3, unsigned short spin4 )
{
    auto *T = new Lattice::Term(4);
    bool Operators[4]          =      { true,  true,  false, false };
    std::string Labels[4]      =      { Label, Label, Label, Label };
    unsigned short Spins[4]    =      { spin1, spin2, spin4, spin3 };
    unsigned short Orbitals[4] =      { orbital1, orbital2, orbital4, orbital3 };
    T->OperatorSequence.assign(Operators,Operators+4);
    T->SiteLabels.assign(Labels,Labels+4);
    T->Orbitals.assign(Orbitals,Orbitals+4);
    T->Spins.assign(Spins,Spins+4);
    T->Value = Value;
    return T;
}


class ThreeFreq
{
public:
    ThreeFreq(int n1, bool negative1, int n2, bool negative2, int n3, bool negative3);
    unsigned int size(){
        return vec_freqs.size();
    }

    struct three_freqs{
        int w1, w2, w3;
    };
    struct three_freqs get_freqs(int i){
        return vec_freqs[i];
    }
private:
    std::vector<three_freqs> vec_freqs;
};

ThreeFreq::ThreeFreq(int n1, bool negative1, int n2, bool negative2, int n3, bool negative3)
{
    // make frequency list
    int size1 = negative1 ? 2*n1 : n1;
    int size2 = negative2 ? 2*n2 : n2;
    int size3 = negative3 ? 2*n3 : n3;
    for (int i1 = 0; i1 < size1; i1++) {
        for (int i2 = 0; i2 < size2; i2++) {
            for (int i3 = 0; i3 < size3; i3++) {
                three_freqs freqs = {
                    negative1 ? i1 - n1 : i1,
                    negative2 ? i2 - n2 : i2,
                    negative3 ? i3 - n3 : i3,
                };
                vec_freqs.push_back(freqs);
            }
        }
    }
}


/* Generic tips:
 * The calculation is done by computing a set of objects in the following order:
 * Lattice -> IndexClassification -> IndexHamiltonian -> Symmetrizer ->
 * -> StatesClassification -> Hamiltonian -> FieldOperator
 ;
 * (for thermal objects, such as GFs in Matsubara domain)
 * -> DensityMatrix -> Greens Function
 *                  -> TwoParticle GF   -> Vertex4
 * The detailed explanation of each class is given below.
 */

int main(int argc, char* argv[])
{
    clock_t time_start = clock(), time_temp;
    boost::mpi::environment MpiEnv(argc, argv);
    boost::mpi::communicator world = boost::mpi::communicator();

    bool verbose = world.rank() == 0;

    if(argc != 2){
        std::cerr << "Usage: " << argv[0] << " input_file_name" << std::endl;
        exit(1);
    }
    std::string filein(argv[1]);
    Params prms;
    prms.read(filein);
    if(verbose){
        prms.print();
    }

    // -----------------------------------------------------------------------

    Lattice L;
    L.addSite(new Lattice::Site("A", prms.n_orb, 2));

    // -----------------------------------------------------------------------
    if(verbose) print_section("Indices");

    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare();

    // Print which indices we have
    if(verbose){
        L.printSites();
        IndexInfo.printIndices();
    }
    // Save the total number of indices.
    ParticleIndex IndexSize = IndexInfo.getIndexSize();
    assert(IndexSize == 2*prms.n_orb);

    // -----------------------------------------------------------------------
    // converter

    struct converter_info{
        ParticleIndex index;
        unsigned short spn;
        unsigned short orb;
    };
    std::vector<converter_info> converter(IndexSize);

    for(ParticleIndex i=0; i<IndexSize; i++){
        IndexClassification::IndexInfo info = IndexInfo.getInfo(i);
        unsigned short spn = info.Spin;
        unsigned short orb = info.Orbital;

        int j;  // converter index
        if(prms.index_order==0){
            // (0, up), (1, up), ..., (0, down), (1, down), ...
            j = prms.n_orb * spn + orb;
        }
        else{
            // (0, up), (0, down), (1, up), (1, down), ...
            j = 2 * orb + spn;
        }
        converter[j].index = i;
        converter[j].spn = spn;
        converter[j].orb = orb;
    }

    if(verbose){
        for(int j=0; j<converter.size(); j++){
            std::cout << "Converter " << j << " : ";
            std::cout << "Index " << converter[j].index;
            std::cout << ", orb " << converter[j].orb;
            std::cout << ", spn " << converter[j].spn << std::endl;
        }
    }

    // check if converter is properly set (no overlap in the indices)
    {
        // remove overlap in converter using std::set
        std::set<ParticleIndex> set_converter;
        for(auto info : converter){
            assert(info.index >=0 && info.index<converter.size());
            set_converter.insert(info.index);
        }
        // size is not reduced
        assert(set_converter.size() == converter.size());
    }

    // -----------------------------------------------------------------------
    if(verbose) print_section("Terms");

    // set H_0
    {
        ReadDataFile rdf(prms.file_h0, 2, 1);
        while( rdf.read_line() ){
            unsigned short s1 = converter[rdf.get_index(0)].spn;
            unsigned short o1 = converter[rdf.get_index(0)].orb;
            unsigned short s2 = converter[rdf.get_index(1)].spn;
            unsigned short o2 = converter[rdf.get_index(1)].orb;
            double val = rdf.get_val(0);
            // c^+_{o1,s1} c_{o2,s2}
            L.addTerm(OneBodyTerm("A", val, o1, o2, s1, s2));
        }
    }

    // set U_{ijkl}
    {
        ReadDataFile rdf(prms.file_umat, 4, 1);
        while( rdf.read_line() ){
            unsigned short s1 = converter[rdf.get_index(0)].spn;
            unsigned short o1 = converter[rdf.get_index(0)].orb;
            unsigned short s2 = converter[rdf.get_index(1)].spn;
            unsigned short o2 = converter[rdf.get_index(1)].orb;
            unsigned short s3 = converter[rdf.get_index(2)].spn;
            unsigned short o3 = converter[rdf.get_index(2)].orb;
            unsigned short s4 = converter[rdf.get_index(3)].spn;
            unsigned short o4 = converter[rdf.get_index(3)].orb;
            double val = rdf.get_val(0);
            // (1/2) U c^+_{o1,s1} c^+_{o2,s2} c_{o4,s4} c_{o3,s3}
            L.addTerm(TwoBodyTerm("A", val/2., o1, o2, o3, o4, s1, s2, s3, s4));
        }
    }

    // -----------------------------------------------------------------------
    // Let us now print which sites and terms are defined.
    if(verbose){
        INFO("Terms with 2 operators");
        L.printTerms(2);
        INFO("Terms with 4 operators");
        L.printTerms(4);
    };

    // -----------------------------------------------------------------------
    if(verbose) print_section("Matrix element storage");

    IndexHamiltonian Storage(&L,IndexInfo);
    Storage.prepare();
    // Print out the Hamiltonian.
    if(verbose){
        INFO("H =\n" << Storage);
    }

    // define operators for checking symmetry
    OperatorPresets::N op_N(IndexSize);
    Operator op_Sz = OperatorPresets::SzSite(IndexInfo, "A");
    Operator op_Lz = OperatorPresets::LzSite(IndexInfo, "A");
    Operator op_Jz = op_Sz + op_Lz;
    if(verbose){
        INFO("N  = \n" << op_N);
        INFO("Sz = \n" << op_Sz);
        INFO("Lz = \n" << op_Lz);
        INFO("Jz = \n" << op_Jz);
    }

    // -----------------------------------------------------------------------
    if(verbose) print_section("Symmetry");

    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute();
    // [H, N] and [H, Sz] are checked in compute() by default

    // compute [H, N] for printing
    print_commutation(Storage.commutes(op_N), "N ", verbose);

    // compute [H, Sz] for printing
    bool flag_sz_commute = Storage.commutes(op_Sz);
    print_commutation(flag_sz_commute, "Sz", verbose);
    if(!flag_sz_commute && prms.flag_spin_conserve){
        std::cerr << "ERROR: Sz does not commute with H, but flag_spin_conserve=true" << std::endl;
        exit(3);
    }

    // Perform additional symmetry check with Lz and Jz
    print_commutation(Symm.checkSymmetry(op_Lz), "Lz", verbose);
    print_commutation(Symm.checkSymmetry(op_Jz), "Jz", verbose);

    if(verbose){
        INFO("Conserved quantum numbers " << Symm.getQuantumNumbers());
    }

    // -----------------------------------------------------------------------
    if(verbose) print_section("States classification");

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    // fileout
    if (!world.rank()) {
        unsigned long n_states = S.getNumberOfStates();
        BlockNumber n_blocks = S.NumberOfBlocks();
        if(verbose){
            INFO("Number of States is " << n_states);
            INFO("Number of Blocks is " << n_blocks);
        }
        std::ofstream fout(prms.file_states);
        fout << "Number of States is " << n_states << std::endl;
        fout << "Number of Blocks is " << n_blocks << std::endl;
        fout << "Block list:\n(# size quantum_numbers)" << std::endl;
        for(BlockNumber i=0; i<n_blocks; i++){
            // INFO(S.getQuantumNumbers(i));
            fout << i << "  " << S.getBlockSize(i) << "  " << S.getQuantumNumbers(i) << std::endl;
        }
        fout.close();
    }
    world.barrier();

    // -----------------------------------------------------------------------
    if(verbose) print_section("Hamiltonian");
    time_temp = clock();

    Hamiltonian H(IndexInfo, Storage, S);
    H.prepare();
    H.compute(world);

    if(verbose){
        INFO("Ground-state energy: " << H.getGroundEnergy());
    }

    // file out eigenvalues
    if (!world.rank()){
        // create a list of pairs of eigenvalue and quantum numers
        std::vector<EigenSystem> eigen;
        for(BlockNumber i=0; i<S.NumberOfBlocks(); i++){
            const HamiltonianPart &H_part = H.getPart(i);
            // INFO(H_part.getQuantumNumbers());
            // INFO(H_part.getEigenValues());
            QuantumNumbers q = H_part.getQuantumNumbers();
            for(InnerQuantumState j=0; j<H_part.getSize(); j++){
                eigen.emplace_back( std::make_pair( H_part.getEigenValue(j), q ) );
            }
        }
        // sort eigenvalues in ascending order
        std::sort(eigen.begin(), eigen.end());
        // print all eigenvalues and corresponding quantum numbers
        std::ofstream fout(prms.file_eigen);
        for( const auto &eig : eigen )
            fout << eig.first << "  " << eig.second << std::endl;
        fout.close();
        // fprint_eigen(eigen);
    }
    world.barrier();
    if(verbose) print_time(time_temp, "eigenstates");

    // Finish if no physical quantities will be computed
    if( !prms.flag_gf && !prms.flag_vx) {
        return 0;
    }

    // -----------------------------------------------------------------------
    if(verbose) print_section("Density Matrix");

    DensityMatrix rho(S, H, prms.beta);
    rho.prepare();
    rho.compute();

    /** Minimal magnitude of the weight (density matrix) to take it into account. */
    RealType DensityMatrixCutoff = 1e-10;

    // Truncate blocks that have only negligible contribute to GF and TwoParticleGF
    rho.truncateBlocks(DensityMatrixCutoff, verbose);

    // file out retained states
    if (!world.rank()) {
        std::ofstream fout(prms.file_retained);
        fout << "Block retained:\n(# size quantum_numbers)" << std::endl;
        for (BlockNumber i = 0; i < S.NumberOfBlocks(); i++) {
            // INFO(S.getQuantumNumbers(i));
            if (rho.isRetained(i))
                fout << i << "  " << S.getBlockSize(i) << "  " << S.getQuantumNumbers(i) << std::endl;
        }
        fout.close();
    }
    world.barrier();

    // -----------------------------------------------------------------------
    if(verbose) print_section("Creation/Annihilation operators");
    time_temp = clock();

//    std::vector <CreationOperator> CX;
    std::vector<std::unique_ptr<CreationOperator> > CX;
    for( auto info : converter ){
        CX.emplace_back(new CreationOperator(IndexInfo, S, H, info.index));
        CX.back()->prepare();
        CX.back()->compute();
    }

//    std::vector <AnnihilationOperator> C;
    std::vector<std::unique_ptr<AnnihilationOperator> > C;
    for( auto info : converter ){
        C.emplace_back(new AnnihilationOperator(IndexInfo, S, H, info.index));
        C.back()->prepare();
        C.back()->compute();
    }

    world.barrier();
    if(verbose) print_time(time_temp, "Creation/Annihilation op");

    // -----------------------------------------------------------------------
    if(prms.flag_gf){
        if(verbose) print_section("Single-particle Green function");
        time_temp = clock();

        std::unique_ptr<WriteDataFile> wdf;
        if (!world.rank()){
            wdf.reset(new WriteDataFile(prms.file_gf));
        }

        for(int i=0; i<IndexSize; i++){
            for(int j=0; j<IndexSize; j++){
                // skip if spin components of i and j are different
                if(prms.flag_spin_conserve && converter[i].spn != converter[j].spn){
                    continue;
                }

                GreensFunction GF(S,H,*C[i],*CX[j], rho);
                GF.prepare();
                GF.compute();

                std::vector<ComplexType> giw(prms.n_w);
                for(int iw=0; iw<prms.n_w; iw++){
                    giw[iw] = GF(iw);
                }

                if (!world.rank()) {
                    std::stringstream ss;
                    ss << "# " << i << " " << j;
                    wdf->write_str(ss.str());
                    wdf->write_vector(giw);
                }
            }
        }

        world.barrier();
        if(verbose) print_time(time_temp, "GF");
    }

    // -----------------------------------------------------------------------
    if(prms.flag_vx){
        if(verbose) print_section("Two-particle Green functions");

//        std::unique_ptr<WriteDataFile> wdf;
//        if (!world.rank()){
//            wdf.reset(new WriteDataFile(prms.file_gf));
//        }

        for(int i1=0; i1<IndexSize; i1++){
            for(int i2=0; i2<IndexSize; i2++){
                for(int i3=0; i3<IndexSize; i3++){
                    for(int i4=0; i4<IndexSize; i4++){
                        // check spin components
                        if(prms.flag_spin_conserve && converter[i1].spn + converter[i4].spn != converter[i2].spn + converter[i3].spn){
                            continue;
                        }

                        // TODO: check def of chi_{ijkl}
                        TwoParticleGF Chi(S, H, *C[i2], *C[i3], *CX[i1], *CX[i4], rho);
                        /** A difference in energies with magnitude less than this value is treated as zero. */
                        Chi.ReduceResonanceTolerance = 1e-8;
                        /** Minimal magnitude of the coefficient of a term to take it into account. */
                        Chi.CoefficientTolerance = 1e-16;
                        /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
                        Chi.MultiTermCoefficientTolerance = 1e-6;

                        time_temp = clock();
                        Chi.prepare();
                        std::vector <boost::tuple<ComplexType, ComplexType, ComplexType>> freqs_2pgf;
                        Chi.compute(false, freqs_2pgf, world);
                        world.barrier();
                        if(verbose) print_time(time_temp, "TwoParticleGF.compute");

                        time_temp = clock();

//                        ThreeFreq freq(prms.n_w2b, 2*prms.n_w2f, 2*prms.n_w2f);
                        ThreeFreq Freq(prms.n_w2b, false, prms.n_w2f, true, prms.n_w2f, true);
                        std::vector<ComplexType> x(Freq.size());

                        for(int i=0; i<Freq.size(); i++){
                            struct ThreeFreq::three_freqs freqs = Freq.get_freqs(i);
                            int wb = freqs.w1;
                            int wf1 = freqs.w2;
                            int wf2 = freqs.w3;
                            x[i] = -Chi(wb+wf1, wf2, wf1);
                        }

//                        if (!world.rank()) {
//                            std::stringstream ss;
//                            ss << "# " << i1 << " " << i2 << " " << i3 << " " << i4;
//                            wdf->write_str(ss.str());
//                            wdf->write_vector(x);
//                        }
                        if (!world.rank()) {
                            // make directory
                            make_dir(prms.dir_vx);
                            // filename
                            std::stringstream ss;
                            ss << prms.dir_vx << "/" << i1 << "_" << i2 << "_" << i3 << "_" << i4 << ".dat";
                            std::string filename(ss.str());
                            // fileout
                            WriteDataFile wdf(filename);
                            wdf.write_vector(x);
                        }

                        if(verbose) print_time(time_temp, "TwoParticleGF.getValues");
                    }
                }
            }
        }
    }

    if(verbose) print_time(time_start, "Total");
}
