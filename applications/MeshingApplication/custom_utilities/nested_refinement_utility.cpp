//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "nested_refinement_utility.h"


namespace Kratos
{
/// Default constructor
NestedRefinementUtility::NestedRefinementUtility(ModelPart& rModelPart) :
    mrModelPart(rModelPart) {}

/// Destructor
NestedRefinementUtility::~NestedRefinementUtility() {}

/// Turn back information as a string.
std::string NestedRefinementUtility::Info() const {
    return "Nested refinement utility.";
}

/// Print information about this object.
void NestedRefinementUtility::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Nested refinement utility.";
}

/// Print object's data.
void NestedRefinementUtility::PrintData(std::ostream& rOStream) const {
    rOStream << "Nested refinement utility constructed with:\n";
    rOStream << "   Model part: " << mrModelPart.Info() << "\n";
}

///
void NestedRefinementUtility::Refine() {

}

}  // namespace Kratos.


