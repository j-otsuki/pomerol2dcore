#include "OperatorPresetsExtra.h"

namespace Pomerol {
namespace OperatorPresets {

Operator SzSite(IndexClassification &indices, const std::string &Label, const std::vector<double> &SzValues)
{
    Operator op;
    for(ParticleIndex i=0; i<indices.getIndexSize(); i++){
        IndexClassification::IndexInfo info = indices.getInfo(i);
        if( info.SiteLabel == Label ){
            assert( info.Spin < SzValues.size() );
            op += OperatorPresets::n(i) * SzValues[info.Spin];
        }
    }
    return op;
}

Operator SzSite(IndexClassification &indices, const std::string &Label)
{
    // find # of spin components
    int Spins=1; // Spins = 2*S+1
    for(ParticleIndex i=0; i<indices.getIndexSize(); i++){
        IndexClassification::IndexInfo info = indices.getInfo(i);
        if( info.SiteLabel == Label )
            Spins = std::max(Spins, info.Spin+1);
    }
    // INFO("'Spins' is set at " << Spins);
    double S = (double)(Spins - 1) / 2.;

    std::vector<double> SzValues;  // {-S, -S+1, ..., S}
    for(int i=0; i<Spins; i++)  SzValues.push_back((double)i-S);

    return SzSite(indices, Label, SzValues);
}

Operator LzSite(IndexClassification &indices, const std::string &Label, const std::vector<double> &LzValues)
{
    Operator op;
    for(ParticleIndex i=0; i<indices.getIndexSize(); i++){
        IndexClassification::IndexInfo info = indices.getInfo(i);
        if( info.SiteLabel == Label ){
            assert( info.Orbital < LzValues.size() );
            op += OperatorPresets::n(i) * LzValues[info.Orbital];
        }
    }
    return op;
}

Operator LzSite(IndexClassification &indices, const std::string &Label)
{
    // find # of orbital components
    int Orbitals=1; // Orbitals = 2*L+1
    for(ParticleIndex i=0; i<indices.getIndexSize(); i++){
        IndexClassification::IndexInfo info = indices.getInfo(i);
        if( info.SiteLabel == Label )
            Orbitals = std::max(Orbitals, info.Orbital+1);
    }
    // INFO("'Orbitals' is set at " << Orbitals);
    double L = (double)(Orbitals - 1) / 2.;

    std::vector<double> LzValues;  // {-L, -L+1, ..., L}
    for(int i=0; i<Orbitals; i++)  LzValues.push_back((double)i-L);

    return LzSite(indices, Label, LzValues);
}


} // end of namespace OperatorPresets
} // end of namespace Pomerol
