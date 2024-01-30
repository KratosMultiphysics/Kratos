//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Guillermo Casas
//

#ifndef KRATOS_MULTIAXIAL_CONTROL_MODULE_FEM_DEM_GENERALIZED_2D_UTILITIES
#define KRATOS_MULTIAXIAL_CONTROL_MODULE_FEM_DEM_GENERALIZED_2D_UTILITIES

// Project includes

// Application includes
#include "custom_utilities/multiaxial_control_module_generalized_2d_utilities.hpp"
#include "dem_structures_coupling_application_variables.h"

namespace Kratos
{
class KRATOS_API(DEM_STRUCTURES_COUPLING_APPLICATION) MultiaxialControlModuleFEMDEMGeneralized2DUtilities 
    : public MultiaxialControlModuleGeneralized2DUtilities
{
public:

KRATOS_CLASS_POINTER_DEFINITION(MultiaxialControlModuleFEMDEMGeneralized2DUtilities);

/// Defining a table with double argument and result type as table type.
typedef Table<double,double> TableType;

/// Default constructor.

MultiaxialControlModuleFEMDEMGeneralized2DUtilities(ModelPart& rDemModelPart,
                                ModelPart& rFemModelPart,
                                Parameters& rParameters
                                ) : MultiaxialControlModuleGeneralized2DUtilities(
                                    rDemModelPart,
                                    rFemModelPart,
                                    rParameters)
{

}

/// Destructor.

~MultiaxialControlModuleFEMDEMGeneralized2DUtilities() override {}

//***************************************************************************************************************
//***************************************************************************************************************

// Before FEM and DEM solution
void ExecuteInitialize() override;

// Before FEM and DEM solution
void ExecuteInitializeSolutionStep() override;

// After FEM and DEM solution
void ExecuteFinalizeSolutionStep() override;

//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member r_variables
///@{

///@}
///@name Protected member r_variables
///@{ template<class T, std::size_t dim>

///@}
///@name Protected Operators
///@{

///@}
///@name Protected Operations
///@{

Vector MeasureReactionStress(const Variable<array_1d<double,3>>& rVariable) override;

///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{

///@}
///@name Protected LifeCycle
///@{

///@}

private:

///@name Static Member r_variables
///@{

///@}
///@name Member r_variables
///@{
///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{

///@}
///@name Private  Access
///@{

///@}
///@name Private Inquiry
///@{

///@}
///@name Un accessible methods
///@{

/// Assignment operator.
MultiaxialControlModuleFEMDEMGeneralized2DUtilities & operator=(MultiaxialControlModuleFEMDEMGeneralized2DUtilities const& rOther);

///@}

}; // Class MultiaxialControlModuleFEMDEMGeneralized2DUtilities

}  // namespace Python.

#endif // KRATOS_MULTIAXIAL_CONTROL_MODULE_FEM_DEM_GENERALIZED_2D_UTILITIES
