/** \file OperatorPresets.h
**  \brief Declarations of the OperatorPresets class. This is some workaround to generate easy classes
**
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_OPERATOR_PRESETS_EXTRA_H__
#define __INCLUDE_OPERATOR_PRESETS_EXTRA_H__

// #include "Misc.h"
// #include "Operator.h"
#include <pomerol/OperatorPresets.h>

namespace Pomerol {
namespace OperatorPresets {

/*
Added by JO
 */
Operator SzSite(IndexClassification &indices, const std::string &Label, const std::vector<double> &SzValues);
Operator SzSite(IndexClassification &indices, const std::string &Label);
Operator LzSite(IndexClassification &indices, const std::string &Label, const std::vector<double> &LzValues);
Operator LzSite(IndexClassification &indices, const std::string &Label);


} // end of namespace OperatorPresets
} // end of namespace Pomerol

#endif // endif :: #ifndef __INCLUDE_OPERATOR_PRESETS_H__
