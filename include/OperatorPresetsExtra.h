#ifndef POMEROL2DCORE_OPERATOR_PRESETS_EXTRA_H
#define POMEROL2DCORE_OPERATOR_PRESETS_EXTRA_H

// #include "Misc.h"
// #include "Operator.h"
#include <pomerol/OperatorPresets.h>

namespace Pomerol {
namespace OperatorPresets {

Operator SzSite(IndexClassification &indices, const std::string &Label, const std::vector<double> &SzValues);
Operator SzSite(IndexClassification &indices, const std::string &Label);
Operator LzSite(IndexClassification &indices, const std::string &Label, const std::vector<double> &LzValues);
Operator LzSite(IndexClassification &indices, const std::string &Label);

} // end of namespace OperatorPresets
} // end of namespace Pomerol

#endif // POMEROL2DCORE_OPERATOR_PRESETS_EXTRA_H
