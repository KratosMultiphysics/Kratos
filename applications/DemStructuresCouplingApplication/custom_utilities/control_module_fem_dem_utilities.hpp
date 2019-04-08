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
//


#ifndef KRATOS_CONTROL_MODULE_FEM_DEM_UTILITIES
#define KRATOS_CONTROL_MODULE_FEM_DEM_UTILITIES

// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

#include "dem_structures_coupling_application_variables.h"


namespace Kratos
{
class ControlModuleFemDemUtilities
{
public:
typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

KRATOS_CLASS_POINTER_DEFINITION(ControlModuleFemDemUtilities);

/// Default constructor.

ControlModuleFemDemUtilities(const double FaceArea)
{
    mFaceArea = FaceArea;
}

/// Destructor.

virtual ~ControlModuleFemDemUtilities(){}

//***************************************************************************************************************
//***************************************************************************************************************

void UpdateChargingVelocity(ModelPart& rFEMTopModelPart, ModelPart& rDEMTopModelPart,
                                ModelPart& rFEMBotModelPart, ModelPart& rDEMBotModelPart)
{
    // TODO

    const int NCons = static_cast<int>(mrModelPart.Conditions().size());
    ModelPart::ConditionsContainerType::iterator con_begin = mrModelPart.ConditionsBegin();
    double FaceArea = 0.0;

    #pragma omp parallel for reduction(+:FaceArea)
    for(int i = 0; i < NCons; i++)
    {
        ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;

        FaceArea += itCond->GetGeometry().Area();
    }

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;
    ComponentType ReactionVarComponent = KratosComponents< ComponentType >::Get(mReactionVariableName);
    const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
    double FaceReaction = 0.0;

    #pragma omp parallel for reduction(+:FaceReaction)
    for(int i = 0; i<NNodes; i++)
    {
        ModelPart::NodesContainerType::iterator it = it_begin + i;

        FaceReaction += it->FastGetSolutionStepValue(ReactionVarComponent);
    }

    const double ReactionStress = FaceReaction/FaceArea;
    TableType::Pointer pTargetStressTable = mrModelPart.pGetTable(mpTargetStressTableId);
    const double DeltaTime = mrModelPart.GetProcessInfo()[DELTA_TIME];
    const double CurrentTime = mrModelPart.GetProcessInfo()[TIME];
    const double NextTargetStress = pTargetStressTable->GetValue(CurrentTime+DeltaTime);

    const double DeltaVelocity = (NextTargetStress - ReactionStress)*mCompressionLength/(mYoungModulus*DeltaTime) - mVelocity;
    mVelocity += mVelocityFactor * DeltaVelocity;

    ComponentType TargetStressVarComponent = KratosComponents< ComponentType >::Get(mTargetStressVariableName);
    ComponentType ReactionStressVarComponent = KratosComponents< ComponentType >::Get(mReactionStressVariableName);
    #pragma omp parallel for
    for(int i = 0; i<NNodes; i++)
    {
        ModelPart::NodesContainerType::iterator it = it_begin + i;

        it->FastGetSolutionStepValue(TargetStressVarComponent) = pTargetStressTable->GetValue(CurrentTime);
        it->FastGetSolutionStepValue(ReactionStressVarComponent) = ReactionStress;
    }
}



//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const
{
    return "";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const
{
}


///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member r_variables
///@{

    double mFaceArea;

///@}
///@name Protected member r_variables
///@{ template<class T, std::size_t dim>


///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{


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
ControlModuleFemDemUtilities & operator=(ControlModuleFemDemUtilities const& rOther);


///@}

}; // Class ControlModuleFemDemUtilities

}  // namespace Python.

#endif // KRATOS_CONTROL_MODULE_FEM_DEM_UTILITIES
