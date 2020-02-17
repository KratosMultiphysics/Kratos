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

#include "includes/table.h"
#include "includes/kratos_parameters.h"

#include "dem_structures_coupling_application_variables.h"


namespace Kratos
{
class ControlModuleFemDemUtilities
{
public:

KRATOS_CLASS_POINTER_DEFINITION(ControlModuleFemDemUtilities);

/// Defining a table with double argument and result type as table type.
typedef Table<double,double> TableType;

/// Default constructor.

ControlModuleFemDemUtilities(ModelPart& rFemModelPart,
                            ModelPart& rDemModelPart,
                            Parameters& rParameters
                            ) :
                            mrFemModelPart(rFemModelPart),
                            mrDemModelPart(rDemModelPart)
{
    KRATOS_TRY

    Parameters default_parameters( R"(
        {
            "variable_name": "DISPLACEMENT",
            "reaction_variable_name": "REACTION",
            "dem_force_variable_name": "TOTAL_FORCES",
            "target_stress_variable_name": "TARGET_STRESS",
            "reaction_stress_variable_name": "REACTION_STRESS",
            "loading_velocity_variable_name": "LOADING_VELOCITY",
            "imposed_direction" : 2,
            "target_stress_table_id" : 0,
            "initial_velocity" : 0.0,
            "limit_velocity" : 1.0,
            "velocity_factor" : 1.0,
            "compression_length" : 1.0,
            "face_area": 1.0,
            "young_modulus" : 1.0e7,
            "stress_increment_tolerance": 1000.0,
            "update_stiffness": true,
            "start_time" : 0.0
        }  )" );

    // Now validate agains defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    const int ImposedDirection = rParameters["imposed_direction"].GetInt();
    if(ImposedDirection == 2)
    {
        mVariableName = rParameters["variable_name"].GetString() + std::string("_Z");
        mReactionVariableName = rParameters["reaction_variable_name"].GetString() + std::string("_Z");
        mDemForceVariableName = rParameters["dem_force_variable_name"].GetString() + std::string("_Z");
        mTargetStressVariableName = rParameters["target_stress_variable_name"].GetString() + std::string("_Z");
        mReactionStressVariableName = rParameters["reaction_stress_variable_name"].GetString() + std::string("_Z");
        mLoadingVelocityVariableName = rParameters["loading_velocity_variable_name"].GetString() + std::string("_Z");
        mDemVelocityVariableName = "IMPOSED_VELOCITY_Z_VALUE";
    }
    else if(ImposedDirection == 1)
    {
        mVariableName = rParameters["variable_name"].GetString() + std::string("_Y");
        mReactionVariableName = rParameters["reaction_variable_name"].GetString() + std::string("_Y");
        mDemForceVariableName = rParameters["dem_force_variable_name"].GetString() + std::string("_Y");
        mTargetStressVariableName = rParameters["target_stress_variable_name"].GetString() + std::string("_Y");
        mReactionStressVariableName = rParameters["reaction_stress_variable_name"].GetString() + std::string("_Y");
        mLoadingVelocityVariableName = rParameters["loading_velocity_variable_name"].GetString() + std::string("_Y");
        mDemVelocityVariableName = "IMPOSED_VELOCITY_Y_VALUE";
    }
    else
    {
        mVariableName = rParameters["variable_name"].GetString() + std::string("_X");
        mReactionVariableName = rParameters["reaction_variable_name"].GetString() + std::string("_X");
        mDemForceVariableName = rParameters["dem_force_variable_name"].GetString() + std::string("_X");
        mTargetStressVariableName = rParameters["target_stress_variable_name"].GetString() + std::string("_X");
        mReactionStressVariableName = rParameters["reaction_stress_variable_name"].GetString() + std::string("_X");
        mLoadingVelocityVariableName = rParameters["loading_velocity_variable_name"].GetString() + std::string("_X");
        mDemVelocityVariableName = "IMPOSED_VELOCITY_X_VALUE";
    }
    mTargetStressTableId = rParameters["target_stress_table_id"].GetInt();
    mVelocity = rParameters["initial_velocity"].GetDouble();
    Variable<double> DemVelocityVar = KratosComponents< Variable<double> >::Get( mDemVelocityVariableName );
    mrDemModelPart[DemVelocityVar] = mVelocity;
    mLimitVelocity = rParameters["limit_velocity"].GetDouble();
    mVelocityFactor = rParameters["velocity_factor"].GetDouble();
    mStartTime = rParameters["start_time"].GetDouble();
    mFaceArea = rParameters["face_area"].GetDouble();
    mStressIncrementTolerance = rParameters["stress_increment_tolerance"].GetDouble();
    mUpdateStiffness = rParameters["update_stiffness"].GetBool();
    mReactionStressOld = 0.0;
    mStiffness = rParameters["young_modulus"].GetDouble()*mFaceArea/rParameters["compression_length"].GetDouble();

    mrDemModelPart.GetProcessInfo()[TARGET_STRESS_Z] = 0.0;

    KRATOS_CATCH("");
}

/// Destructor.

virtual ~ControlModuleFemDemUtilities(){}

//***************************************************************************************************************
//***************************************************************************************************************

// Before FEM and DEM solution
void ExecuteInitialize()
{
    KRATOS_TRY;

    const int NNodes = static_cast<int>(mrFemModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator it_begin = mrFemModelPart.NodesBegin();
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;
    ComponentType VarComponent = KratosComponents< ComponentType >::Get(mVariableName);

    #pragma omp parallel for
    for(int i = 0; i<NNodes; i++)
    {
        ModelPart::NodesContainerType::iterator it = it_begin + i;

        it->Fix(VarComponent);
        it->FastGetSolutionStepValue(VarComponent) = 0.0;
    }

    KRATOS_CATCH("");
}

// Before FEM and DEM solution
void ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    const int NNodes = static_cast<int>(mrFemModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator it_begin = mrFemModelPart.NodesBegin();
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;
    ComponentType VarComponent = KratosComponents< ComponentType >::Get(mVariableName);
    const double DeltaTime = mrFemModelPart.GetProcessInfo()[DELTA_TIME];

    #pragma omp parallel for
    for(int i = 0; i<NNodes; i++)
    {
        ModelPart::NodesContainerType::iterator it = it_begin + i;

        it->FastGetSolutionStepValue(VarComponent) += mVelocity * DeltaTime;
    }

    KRATOS_CATCH("");
}

// After FEM and DEM solution
void ExecuteFinalizeSolutionStep()
{
    const double CurrentTime = mrFemModelPart.GetProcessInfo()[TIME];

    if(CurrentTime >= mStartTime && mTargetStressTableId > 0)
    {
        // Calculate ReactionStress
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;
        ComponentType ReactionVarComponent = KratosComponents< ComponentType >::Get(mReactionVariableName);
        const int NNodes = static_cast<int>(mrFemModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrFemModelPart.NodesBegin();
        double FaceReaction = 0.0;

        #pragma omp parallel for reduction(+:FaceReaction)
        for(int i = 0; i<NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;

            FaceReaction += it->FastGetSolutionStepValue(ReactionVarComponent);
        }

        ComponentType DemForceVarComponent = KratosComponents< ComponentType >::Get(mDemForceVariableName);
        const int NDemNodes = static_cast<int>(mrDemModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_dem_begin = mrDemModelPart.NodesBegin();

        // The sign of DEM TOTAL_FORCES is opposed to the sign of FEM REACTION
        #pragma omp parallel for reduction(-:FaceReaction)
        for(int i = 0; i<NDemNodes; i++)
        {
            ModelPart::NodesContainerType::iterator it = it_dem_begin + i;

            FaceReaction -= it->FastGetSolutionStepValue(DemForceVarComponent);
        }

        const double ReactionStress = FaceReaction/mFaceArea;

        // Update K if required
        const double DeltaTime = mrFemModelPart.GetProcessInfo()[DELTA_TIME];
        double K_estimated = mStiffness;
        if(mUpdateStiffness == true) {
            if(std::abs(mVelocity) > 1.0e-4*std::abs(mLimitVelocity) &&
                std::abs(ReactionStress-mReactionStressOld) > mStressIncrementTolerance) {

                K_estimated = std::abs((ReactionStress-mReactionStressOld)/(mVelocity * DeltaTime));
            }
            mReactionStressOld = ReactionStress;
            mStiffness = K_estimated;
        }

        // Update velocity
        TableType::Pointer pTargetStressTable = mrFemModelPart.pGetTable(mTargetStressTableId);
        const double NextTargetStress = pTargetStressTable->GetValue(CurrentTime+DeltaTime);
        const double df_target = NextTargetStress - ReactionStress;

        double delta_velocity = df_target/(K_estimated * DeltaTime) - mVelocity;

        if(std::abs(df_target) < mStressIncrementTolerance) {

            delta_velocity = -mVelocity;
        }

        mVelocity += mVelocityFactor * delta_velocity;

        if(std::abs(mVelocity) > std::abs(mLimitVelocity))
        {
            if(mVelocity >= 0.0)
            {
                mVelocity = std::abs(mLimitVelocity);
            }
            else
            {
                mVelocity = - std::abs(mLimitVelocity);
            }
        }

        Variable<double> DemVelocityVar = KratosComponents< Variable<double> >::Get( mDemVelocityVariableName );
        mrDemModelPart[DemVelocityVar] = mVelocity;

        // Output of TARGET_STRESS and REACTION_STRESS (only printed at the FEM nodes).
        // TODO: This output is updated after printing the current step.
        ComponentType TargetStressVarComponent = KratosComponents< ComponentType >::Get(mTargetStressVariableName);
        ComponentType ReactionStressVarComponent = KratosComponents< ComponentType >::Get(mReactionStressVariableName);
        ComponentType LoadingVelocityVarComponent = KratosComponents< ComponentType >::Get(mLoadingVelocityVariableName);
        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;

            it->FastGetSolutionStepValue(TargetStressVarComponent) = pTargetStressTable->GetValue(CurrentTime);
            it->FastGetSolutionStepValue(ReactionStressVarComponent) = ReactionStress;
            // it->FastGetSolutionStepValue(ReactionStressVarComponent) = pTargetStressTable->GetValue(CurrentTime)-ReactionStress;
            it->FastGetSolutionStepValue(LoadingVelocityVarComponent) = mVelocity;
        }

        mrDemModelPart.GetProcessInfo()[TARGET_STRESS_Z] = std::abs(pTargetStressTable->GetValue(CurrentTime));

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

    ModelPart& mrFemModelPart;
    ModelPart& mrDemModelPart;
    std::string mVariableName;
    std::string mDemVelocityVariableName;
    std::string mReactionVariableName;
    std::string mDemForceVariableName;
    std::string mTargetStressVariableName;
    std::string mReactionStressVariableName;
    std::string mLoadingVelocityVariableName;
    unsigned int mTargetStressTableId;
    double mVelocity;
    double mLimitVelocity;
    double mVelocityFactor;
    double mStartTime;
    double mFaceArea;
    double mReactionStressOld;
    double mStressIncrementTolerance;
    double mStiffness;
    bool mUpdateStiffness;

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
