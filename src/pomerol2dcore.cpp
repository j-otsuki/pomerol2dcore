#include <pomerol.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <cassert>
#include <memory>

#include "ReadWrite.h"
#include "Params.h"
#include "OperatorPresetsExtra.h"

using namespace Pomerol;

boost::mpi::communicator world;

typedef std::pair<double, QuantumNumbers> EigenSystem;


// Small routine to make fancy screen output for text.
void print_section (const std::string& str)
{
    if (!world.rank()) {
        std::cout << std::string(str.size(),'=') << std::endl;
        std::cout << str << std::endl;
        std::cout << std::string(str.size(),'=') << std::endl;
    }
}


void print_time(clock_t start, const char *str)
{
    if (!world.rank()) {
        std::cout << "#Time: "
                  << (double)(clock()-start)/CLOCKS_PER_SEC << " sec"
                  << " (" << str << ")"<< std::endl << std::flush;
    }
}


void print_commutation(bool if_commute, const std::string &op)
{
    std::string str_commute(if_commute ? " = 0" : "!= 0");
    if(!world.rank()){
        std::cout << "[H, " << op << "] " << str_commute << std::endl;
    }
}


Lattice::Term* OneBodyTerm ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    Lattice::Term *T = new Lattice::Term(2);
    bool Operators[2]          = { 1, 0 };
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
    Lattice::Term *T = new Lattice::Term(4);
    bool Operators[4]          =      { 1,     1,     0,     0     };
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
private:
    int n1, n2, n3;
public:
    ThreeFreq(int n1, int n2, int n3) : n1(n1), n2(n2), n3(n3) {}
    int size(){
        return n1*n2*n3;
    }
    int operator()(int i, int j, int k){
        return (i * n2 + j) * n3 + k;
    }
};


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
    world = boost::mpi::communicator();

    if(argc != 2){
        std::cerr << "Usage: " << argv[0] << " input_file_name" << std::endl;
        exit(1);
    }
    std::string filein(argv[1]);
    Params prms;
    prms.read(filein);
    prms.print();

    // -----------------------------------------------------------------------

    Lattice L;
    L.addSite(new Lattice::Site("A", prms.n_orb, 2));

    // -----------------------------------------------------------------------
    print_section("Indices");

    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare();

    // Print which indices we have
    if (!world.rank())
        L.printSites();
        IndexInfo.printIndices();
    // Save the total number of indices.
    ParticleIndex IndexSize = IndexInfo.getIndexSize();

    // -----------------------------------------------------------------------
    // i = (orbital, spin)
    // converter from i to (orbital, spin)
    std::vector<int> index2spn;
    std::vector<int> index2orb;
    for(ParticleIndex i=0; i<IndexSize; i++){
        IndexClassification::IndexInfo info = IndexInfo.getInfo(i);
        index2spn.push_back(info.Spin);
        index2orb.push_back(info.Orbital);
        std::cout << i << " " << index2orb[i] << " " << index2spn[i] << std::endl;
    }

    // -----------------------------------------------------------------------
    print_section("Terms");

    // set H_0
    {
        ReadDataFile rdf(prms.file_h0, 2, 1);
        while( rdf.read_line() ){
            int s1 = index2spn[rdf.get_index(0)];
            int o1 = index2orb[rdf.get_index(0)];
            int s2 = index2spn[rdf.get_index(1)];
            int o2 = index2orb[rdf.get_index(1)];
            double val = rdf.get_val(0);
            // c^+_{o1,s1} c_{o2,s2}
            L.addTerm(OneBodyTerm("A", val, o1, o2, s1, s2));
        }
    }

    // set U_{ijkl}
    {
        ReadDataFile rdf(prms.file_umat, 4, 1);
        while( rdf.read_line() ){
            int s1 = index2spn[rdf.get_index(0)];
            int o1 = index2orb[rdf.get_index(0)];
            int s2 = index2spn[rdf.get_index(1)];
            int o2 = index2orb[rdf.get_index(1)];
            int s3 = index2spn[rdf.get_index(2)];
            int o3 = index2orb[rdf.get_index(2)];
            int s4 = index2spn[rdf.get_index(3)];
            int o4 = index2orb[rdf.get_index(3)];
            double val = rdf.get_val(0);
            // c^+_{o1,s1} c^+_{o2,s2} c_{o4,s4} c_{o3,s3}
            L.addTerm(TwoBodyTerm("A", val, o1, o2, o3, o4, s1, s2, s3, s4));
        }
    }

    // -----------------------------------------------------------------------
    // Let us now print which sites and terms are defined.
    if (!world.rank()) {
        INFO("Terms with 2 operators");
        L.printTerms(2);
        INFO("Terms with 4 operators");
        L.printTerms(4);
    };

    // -----------------------------------------------------------------------
    print_section("Matrix element storage");

    IndexHamiltonian Storage(&L,IndexInfo);
    Storage.prepare();
    // Print out the Hamiltonian.
    if (!world.rank()) {
        INFO("H =\n" << Storage);
    }

    // define operators for checking symmetry
    OperatorPresets::N op_N(IndexSize);
    Operator op_Sz = OperatorPresets::SzSite(IndexInfo, "A");
    Operator op_Lz = OperatorPresets::LzSite(IndexInfo, "A");
    Operator op_Jz = op_Sz + op_Lz;
    if (!world.rank()) {
        INFO("N  = \n" << op_N);
        INFO("Sz = \n" << op_Sz);
        INFO("Lz = \n" << op_Lz);
        INFO("Jz = \n" << op_Jz);
    }

    // -----------------------------------------------------------------------
    print_section("Symmetry");

    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute();
    // [H, N] and [H, Sz] are checked in compute() by default

    // compute [H, N] for printing
    print_commutation(Storage.commutes(op_N), "N ");

    // compute [H, Sz] for printing
    bool flag_sz_commute = Storage.commutes(op_Sz);
    print_commutation(flag_sz_commute, "Sz");
    if(!flag_sz_commute && prms.flag_spin_conserve){
        std::cerr << "ERROR: Sz does not commute with H, but flag_spin_conserve=true" << std::endl;
        exit(3);
    }

    // Perform additional symmetry check with Lz and Jz
    print_commutation(Symm.checkSymmetry(op_Lz), "Lz");
    print_commutation(Symm.checkSymmetry(op_Jz), "Jz");

    INFO("Conserved quantum numbers " << Symm.getQuantumNumbers());

    // -----------------------------------------------------------------------
    print_section("States classification");

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    if (!world.rank()) {
        unsigned long n_states = S.getNumberOfStates();
        INFO("Number of States is " << n_states);
        BlockNumber n_blocks = S.NumberOfBlocks();
        INFO("Number of Blocks is " << n_blocks);
        // fileout
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
    print_section("Hamiltonian");
    time_temp = clock();

    Hamiltonian H(IndexInfo, Storage, S);
    H.prepare();
    H.compute(world);

    if (!world.rank())
        INFO("The value of ground energy is " << H.getGroundEnergy());

    if (!world.rank()){
        INFO("Eigenvalues");
        // create a list of pairs of eigenvalue and quantum numers
        std::vector<EigenSystem> eigen;
        for(BlockNumber i=0; i<S.NumberOfBlocks(); i++){
            HamiltonianPart H_part = H.getPart(i);
            // INFO(H_part.getQuantumNumbers());
            // INFO(H_part.getEigenValues());
            QuantumNumbers q = H_part.getQuantumNumbers();
            for(InnerQuantumState j=0; j<H_part.getSize(); j++){
                eigen.push_back( std::make_pair( H_part.getEigenValue(j), q ) );
            }
        }
        // sort eigenvalues in ascending order
        std::sort(eigen.begin(), eigen.end());
        // print all eigenvalues and corresponding quantum numbers
        std::ofstream fout(prms.file_eigen);
        for(int i=0; i<eigen.size(); i++)
            fout << eigen[i].first << "  " << eigen[i].second << std::endl;
        fout.close();
        // fprint_eigen(eigen);
    }
    world.barrier();
    print_time(time_temp, "eigenstates");

    // Finish if no physical quantities will be computed
    if( !prms.flag_gf && !prms.flag_vx) {
        return 0;
    }

    // -----------------------------------------------------------------------
    print_section("Density Matrix");

    DensityMatrix rho(S, H, prms.beta);
    rho.prepare();
    rho.compute();

    /** Minimal magnitude of the weight (density matrix) to take it into account. */
    RealType DensityMatrixCutoff = 1e-10;

    bool verbose = world.rank() == 0;
    // Truncate blocks that have only negligible contribute to GF and TwoParticleGF
    rho.truncateBlocks(DensityMatrixCutoff, verbose);

    if (!world.rank()) {
        // fileout
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
    print_section("Creation/Annihilation operators");
    time_temp = clock();

//    std::vector <CreationOperator> CX;
    std::vector<std::unique_ptr<CreationOperator> > CX;
    for (ParticleIndex i = 0; i < IndexSize; i++) {
//        CX.push_back(CreationOperator(IndexInfo, S, H, i));
        CX.emplace_back(new CreationOperator(IndexInfo, S, H, i));
        CX.back()->prepare();
        CX.back()->compute();
    }

//    std::vector <AnnihilationOperator> C;
    std::vector<std::unique_ptr<AnnihilationOperator> > C;
    for (ParticleIndex i = 0; i < IndexSize; i++) {
//        C.push_back(AnnihilationOperator(IndexInfo, S, H, i));
        C.emplace_back(new AnnihilationOperator(IndexInfo, S, H, i));
        C.back()->prepare();
        C.back()->compute();
    }

    world.barrier();
    print_time(time_temp, "Creation/Annihilation op");

    // -----------------------------------------------------------------------
    if(prms.flag_gf){
        print_section("Single-particle Green function");
        time_temp = clock();

        std::unique_ptr<WriteDataFile> wdf;
        if (!world.rank()){
            wdf.reset(new WriteDataFile(prms.file_gf));
        }

        for(ParticleIndex i=0; i<IndexSize; i++){
            for(ParticleIndex j=0; j<IndexSize; j++) {
                // skip if spin components of i and j are different
                if(prms.flag_spin_conserve && index2spn[i] != index2spn[j]){
                    continue;
                }

                GreensFunction GF(S,H,*C[i],*CX[j], rho);
                GF.prepare();
                GF.compute();

                std::vector<ComplexType> giw(prms.n_w);
                for(int i=0; i<prms.n_w; i++){
                    giw[i] = GF(i);
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
        print_time(time_temp, "GF");
    }

    // -----------------------------------------------------------------------
    if(prms.flag_vx){
        print_section("Two-particle Green functions");

        std::unique_ptr<WriteDataFile> wdf;
        if (!world.rank()){
            wdf.reset(new WriteDataFile(prms.file_gf));
        }

        for(ParticleIndex i=0; i<IndexSize; i++) {
            for (ParticleIndex j = 0; j < IndexSize; j++) {
                for (ParticleIndex k = 0; k < IndexSize; k++) {
                    for (ParticleIndex l = 0; l < IndexSize; l++) {
                        // TODO: check spin components

                        // TODO: check def of chi_{ijkl}
                        TwoParticleGF Chi(S, H, *C[i], *C[j], *CX[k], *CX[l], rho);
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
                        print_time(time_temp, "TwoParticleGF.compute");

                        time_temp = clock();

                        ThreeFreq freq(prms.n_w2b, 2*prms.n_w2f, 2*prms.n_w2f);
                        std::vector<ComplexType> x(freq.size());

                        for(int w1=0; w1<prms.n_w2b; w1++){
                            for(int w2=0; w2<2*prms.n_w2f; w2++) {
                                for (int w3=0; w3<2*prms.n_w2f; w3++) {
                                    x[freq(w1,w2,w3)] = Chi(w1, w2 - prms.n_w2f, w3 - prms.n_w2f);
                                }
                            }
                        }

                        if (!world.rank()) {
                            wdf->write_vector(x);
                        }

                        print_time(time_temp, "TwoParticleGF.getValues");
                    }
                }
            }
        }
    }

    print_time(time_start, "Total");
}
