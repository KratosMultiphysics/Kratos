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
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"

#include "includes/table.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_elements/spheric_continuum_particle.h"

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
            "imposed_direction" : 2,
            "alternate_axis_loading": false,
            "target_stress_table_id" : 0,
            "initial_velocity" : 0.0,
            "limit_velocity" : 1.0,
            "velocity_factor" : 1.0,
            "compression_length" : 1.0,
            "face_area": 1.0,
            "young_modulus" : 1.0e7,
            "stress_increment_tolerance": 100.0,
            "update_stiffness": true,
            "start_time" : 0.0,
            "stress_averaging_time": 1.0e-5
        }  )" );

    // Now validate agains defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mImposedDirection = rParameters["imposed_direction"].GetInt();
    mTargetStressTableId = rParameters["target_stress_table_id"].GetInt();
    mVelocity = rParameters["initial_velocity"].GetDouble();
    if (mImposedDirection == 0) {
        const Variable<double>& DemVelocityVar = KratosComponents< Variable<double> >::Get("IMPOSED_VELOCITY_X_VALUE");
        mrDemModelPart[DemVelocityVar] = mVelocity;
    } else if (mImposedDirection == 1) {
        const Variable<double>& DemVelocityVar = KratosComponents< Variable<double> >::Get("IMPOSED_VELOCITY_Y_VALUE");
        mrDemModelPart[DemVelocityVar] = mVelocity;
    } else {
        const Variable<double>& DemVelocityVar = KratosComponents< Variable<double> >::Get("IMPOSED_VELOCITY_Z_VALUE");
        mrDemModelPart[DemVelocityVar] = mVelocity;
    }
    mLimitVelocity = rParameters["limit_velocity"].GetDouble();
    mVelocityFactor = rParameters["velocity_factor"].GetDouble();
    mStartTime = rParameters["start_time"].GetDouble();
    mStressIncrementTolerance = rParameters["stress_increment_tolerance"].GetDouble();
    mUpdateStiffness = rParameters["update_stiffness"].GetBool();
    mReactionStressOld = 0.0;
    mStiffness = rParameters["young_modulus"].GetDouble()*rParameters["face_area"].GetDouble()/rParameters["compression_length"].GetDouble();
    mStressAveragingTime = rParameters["stress_averaging_time"].GetDouble();
    mVectorOfLastStresses.resize(0);

    mAlternateAxisLoading = rParameters["alternate_axis_loading"].GetBool();
    mXCounter = 1;
    mYCounter = 2;
    mZCounter = 3;

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

        if (mImposedDirection == 0) {
            // X direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Fix(DISPLACEMENT_X);
                it->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
            }
        } else if (mImposedDirection == 1) {
            // Y direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Fix(DISPLACEMENT_Y);
                it->FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
            }
        } else if (mImposedDirection == 2) {
            // Z direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Fix(DISPLACEMENT_Z);
                it->FastGetSolutionStepValue(DISPLACEMENT_Z) = 0.0;
            }
        }

    KRATOS_CATCH("");
}

// Before FEM and DEM solution
void ExecuteInitializeSolutionStep()
{
    // DEM variables
    ModelPart::ElementsContainerType& rElements = mrDemModelPart.GetCommunicator().LocalMesh().Elements();
    // FEM variables
    const double CurrentTime = mrFemModelPart.GetProcessInfo()[TIME];
    const double DeltaTime = mrFemModelPart.GetProcessInfo()[DELTA_TIME];
    const int NNodes = static_cast<int>(mrFemModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator it_begin = mrFemModelPart.NodesBegin();
    TableType::Pointer pTargetStressTable = mrFemModelPart.pGetTable(mTargetStressTableId);

    // Calculate face_area
    double face_area = 0.0;
    // DEM modelpart
    #pragma omp parallel for reduction(+:face_area)
    for (int i = 0; i < (int)rElements.size(); i++) {
        ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + i;
        Element* p_element = ptr_itElem->get();
        SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*>(p_element);
        const double radius = pDemElem->GetRadius();
        face_area += Globals::Pi*radius*radius;
    }
    // FEM modelpart
    const int NCons = static_cast<int>(mrFemModelPart.Conditions().size());
    ModelPart::ConditionsContainerType::iterator con_begin = mrFemModelPart.ConditionsBegin();
    #pragma omp parallel for reduction(+:face_area)
    for(int i = 0; i < NCons; i++)
    {
        ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
        face_area += itCond->GetGeometry().Area();
    }

    // Calculate ReactionStress
    double FaceReaction = 0.0;

    if (mImposedDirection == 0) { // X direction
        #pragma omp parallel for reduction(+:FaceReaction)
        for(int i = 0; i<NNodes; i++){
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            FaceReaction += it->FastGetSolutionStepValue(REACTION_X);
        }
    } else if (mImposedDirection == 1) { // Y direction
        #pragma omp parallel for reduction(+:FaceReaction)
        for(int i = 0; i<NNodes; i++){
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            FaceReaction += it->FastGetSolutionStepValue(REACTION_Y);
        }
    } else if (mImposedDirection == 2) { // Z direction
        #pragma omp parallel for reduction(+:FaceReaction)
        for(int i = 0; i<NNodes; i++){
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            FaceReaction += it->FastGetSolutionStepValue(REACTION_Z);
        }
    }

    const int NDemNodes = static_cast<int>(mrDemModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator it_dem_begin = mrDemModelPart.NodesBegin();

    // The sign of DEM TOTAL_FORCES is opposed to the sign of FEM REACTION
    if (mImposedDirection == 0) { // X direction
        #pragma omp parallel for reduction(-:FaceReaction)
        for(int i = 0; i<NDemNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_dem_begin + i;
            FaceReaction -= it->FastGetSolutionStepValue(TOTAL_FORCES_X);
        }
    } else if (mImposedDirection == 1) { // Y direction
        #pragma omp parallel for reduction(-:FaceReaction)
        for(int i = 0; i<NDemNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_dem_begin + i;
            FaceReaction -= it->FastGetSolutionStepValue(TOTAL_FORCES_Y);
        }
    } else if (mImposedDirection == 2) { // Z direction
        #pragma omp parallel for reduction(-:FaceReaction)
        for(int i = 0; i<NDemNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_dem_begin + i;
            FaceReaction -= it->FastGetSolutionStepValue(TOTAL_FORCES_Z);
        }
    }

    double ReactionStress = FaceReaction/face_area;

    ReactionStress = UpdateVectorOfHistoricalStressesAndComputeNewAverage(ReactionStress);

    // Update K if required
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

    const bool is_time_to_apply_cm = IsTimeToApplyCM();

    if (is_time_to_apply_cm == true)
    {
        // Update Imposed displacement
        if (mImposedDirection == 0) { // X direction
            const Variable<double>& DemVelocityVar = KratosComponents< Variable<double> >::Get("IMPOSED_VELOCITY_X_VALUE");
            mrDemModelPart[DemVelocityVar] = mVelocity;

            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(DISPLACEMENT_X) += mVelocity * DeltaTime;
            }
        } else if (mImposedDirection == 1) { // Y direction
            const Variable<double>& DemVelocityVar = KratosComponents< Variable<double> >::Get("IMPOSED_VELOCITY_Y_VALUE");
            mrDemModelPart[DemVelocityVar] = mVelocity;
            
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(DISPLACEMENT_Y) += mVelocity * DeltaTime;
            }
        } else if (mImposedDirection == 2) { // Z direction
            const Variable<double>& DemVelocityVar = KratosComponents< Variable<double> >::Get("IMPOSED_VELOCITY_Z_VALUE");
            mrDemModelPart[DemVelocityVar] = mVelocity;

            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(DISPLACEMENT_Z) += mVelocity * DeltaTime;
            }
        }
    }

    // Save calculated velocity and reaction for print (only at FEM nodes)
    if (mImposedDirection == 0) { // X direction
        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->FastGetSolutionStepValue(TARGET_STRESS_X) = pTargetStressTable->GetValue(CurrentTime);
            it->FastGetSolutionStepValue(REACTION_STRESS_X) = ReactionStress;
            it->FastGetSolutionStepValue(LOADING_VELOCITY_X) = mVelocity;
        }
    } else if (mImposedDirection == 1) { // Y direction
        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->FastGetSolutionStepValue(TARGET_STRESS_Y) = pTargetStressTable->GetValue(CurrentTime);
            it->FastGetSolutionStepValue(REACTION_STRESS_Y) = ReactionStress;
            it->FastGetSolutionStepValue(LOADING_VELOCITY_Y) = mVelocity;
        }
    } else if (mImposedDirection == 2) { // Z direction
        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->FastGetSolutionStepValue(TARGET_STRESS_Z) = pTargetStressTable->GetValue(CurrentTime);
            it->FastGetSolutionStepValue(REACTION_STRESS_Z) = ReactionStress;
            it->FastGetSolutionStepValue(LOADING_VELOCITY_Z) = mVelocity;

            mrDemModelPart.GetProcessInfo()[TARGET_STRESS_Z] = std::abs(pTargetStressTable->GetValue(CurrentTime));
        }
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
    unsigned int mImposedDirection;
    unsigned int mTargetStressTableId;
    double mVelocity;
    double mLimitVelocity;
    double mVelocityFactor;
    double mStartTime;
    double mReactionStressOld;
    double mStressIncrementTolerance;
    double mStiffness;
    bool mUpdateStiffness;
    std::vector<double> mVectorOfLastStresses;
    double mStressAveragingTime;
    bool mAlternateAxisLoading;
    unsigned int mXCounter;
    unsigned int mYCounter;
    unsigned int mZCounter;

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

double UpdateVectorOfHistoricalStressesAndComputeNewAverage(const double& last_reaction) {
    KRATOS_TRY;
    int length_of_vector = mVectorOfLastStresses.size();
    if (length_of_vector == 0) { //only the first time
        int number_of_steps_for_stress_averaging = (int) (mStressAveragingTime / mrFemModelPart.GetProcessInfo()[DELTA_TIME]);
        if(number_of_steps_for_stress_averaging < 1) number_of_steps_for_stress_averaging = 1;
        mVectorOfLastStresses.resize(number_of_steps_for_stress_averaging);
        KRATOS_INFO("DEM") << " 'number_of_steps_for_stress_averaging' is "<< number_of_steps_for_stress_averaging << std::endl;
    }

    length_of_vector = mVectorOfLastStresses.size();

    if(length_of_vector > 1) {
        for(int i=1; i<length_of_vector; i++) {
            mVectorOfLastStresses[i-1] = mVectorOfLastStresses[i];
        }
    }
    mVectorOfLastStresses[length_of_vector-1] = last_reaction;

    double average = 0.0;
    for(int i=0; i<length_of_vector; i++) {
        average += mVectorOfLastStresses[i];
    }
    average /= (double) length_of_vector;
    return average;

    KRATOS_CATCH("");
}

bool IsTimeToApplyCM(){
    const double current_time = mrFemModelPart.GetProcessInfo()[TIME];
    bool apply_cm = false;

    if(current_time >= mStartTime) {
        if (mAlternateAxisLoading == true) {
            const int step = mrFemModelPart.GetProcessInfo()[STEP];
            if (mImposedDirection == 0) {
                if(step == mXCounter){
                    apply_cm = true;
                    mXCounter += 3;
                }
            } else if (mImposedDirection == 1) {
                if(step == mYCounter){
                    apply_cm = true;
                    mYCounter += 3;
                }
            } else if (mImposedDirection == 2) {
                if(step == mZCounter){
                    apply_cm = true;
                    mZCounter += 3;
                }
            }
        } else {
            apply_cm = true;
        }
    }

    return apply_cm;
}

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
