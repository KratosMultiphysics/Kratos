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


#if !defined(KRATOS_CONTROL_MODULE_PROCESS )
#define  KRATOS_CONTROL_MODULE_PROCESS

#include "includes/table.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "dem_structures_coupling_application_variables.h"

namespace Kratos
{

class ControlModuleProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ControlModuleProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ControlModuleProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process() ,
            mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "imposed_direction" : 0,
                "alternate_axis_loading": false,
                "target_stress_table_id" : 0,
                "initial_velocity" : 0.0,
                "limit_velocity" : 1.0,
                "velocity_factor" : 1.0,
                "compression_length" : 1.0,
                "young_modulus" : 1.0e7,
                "stress_increment_tolerance": 100.0,
                "update_stiffness": true,
                "start_time" : 0.0,
                "stress_averaging_time": 1.0e-5
            }  )" );

        // Now validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mImposedDirection = rParameters["imposed_direction"].GetInt();
        mTargetStressTableId = rParameters["target_stress_table_id"].GetInt();
        mVelocity = rParameters["initial_velocity"].GetDouble();
        mLimitVelocity = rParameters["limit_velocity"].GetDouble();
        mVelocityFactor = rParameters["velocity_factor"].GetDouble();
        mStartTime = rParameters["start_time"].GetDouble();
        mStressIncrementTolerance = rParameters["stress_increment_tolerance"].GetDouble();
        mUpdateStiffness = rParameters["update_stiffness"].GetBool();
        mReactionStressOld = 0.0;
        mStiffness = rParameters["young_modulus"].GetDouble()/rParameters["compression_length"].GetDouble(); // mStiffness is actually a stiffness over an area
        mStressAveragingTime = rParameters["stress_averaging_time"].GetDouble();
        mVectorOfLastStresses.resize(0);
        // NOTE: Alternate axis loading only works for X,Y,Z loading. Radial loading is always applied.
        mAlternateAxisLoading = rParameters["alternate_axis_loading"].GetBool();
        mXCounter = 1;
        mYCounter = 2;
        mZCounter = 3;

        mApplyCM = false;

        // Initialize Variables
        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
        array_1d<double,3> zero_vector = ZeroVector(3);
        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->SetValue(TARGET_STRESS,zero_vector);
            it->SetValue(REACTION_STRESS,zero_vector);
            it->SetValue(LOADING_VELOCITY,zero_vector);
        }

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ControlModuleProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ControlModuleProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

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
        } else {
            // Radial direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Fix(DISPLACEMENT_X);
                it->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
                it->Fix(DISPLACEMENT_Y);
                it->FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
            }
        }

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const double CurrentTime = mrModelPart.GetProcessInfo()[TIME];
        TableType::Pointer pTargetStressTable = mrModelPart.pGetTable(mTargetStressTableId);
        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        double ReactionStress = CalculateReactionStress();
        ReactionStress = UpdateVectorOfHistoricalStressesAndComputeNewAverage(ReactionStress);

        // Check whether this is a loading step for the current axis
        IsTimeToApplyCM();

        if (mApplyCM == true) {

            // Update K if required
            const double DeltaTime = mrModelPart.GetProcessInfo()[DELTA_TIME];
            if (mAlternateAxisLoading == false) {
                if(mUpdateStiffness == true) {
                    mStiffness = EstimateStiffness(ReactionStress,DeltaTime);
                }
            }
            mReactionStressOld = ReactionStress;

            // Update velocity
            const double NextTargetStress = pTargetStressTable->GetValue(CurrentTime+DeltaTime);
            const double df_target = NextTargetStress - ReactionStress;
            double delta_velocity = df_target/(mStiffness * DeltaTime) - mVelocity;

            if(std::abs(df_target) < mStressIncrementTolerance) {
                delta_velocity = -mVelocity;
            }

            mVelocity += mVelocityFactor * delta_velocity;

            if(std::abs(mVelocity) > std::abs(mLimitVelocity)) {
                if(mVelocity >= 0.0) {
                    mVelocity = std::abs(mLimitVelocity);
                } else {
                    mVelocity = - std::abs(mLimitVelocity);
                }
            }

            // Update Imposed displacement
            if (mImposedDirection == 0) { // X direction
                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->FastGetSolutionStepValue(DISPLACEMENT_X) += mVelocity * DeltaTime;
                    // Save calculated velocity and reaction for print
                    it->GetValue(TARGET_STRESS_X) = pTargetStressTable->GetValue(CurrentTime);
                    it->GetValue(REACTION_STRESS_X) = ReactionStress;
                    it->GetValue(LOADING_VELOCITY_X) = mVelocity;
                }
            } else if (mImposedDirection == 1) { // Y direction
                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->FastGetSolutionStepValue(DISPLACEMENT_Y) += mVelocity * DeltaTime;
                    // Save calculated velocity and reaction for print
                    it->GetValue(TARGET_STRESS_Y) = pTargetStressTable->GetValue(CurrentTime);
                    it->GetValue(REACTION_STRESS_Y) = ReactionStress;
                    it->GetValue(LOADING_VELOCITY_Y) = mVelocity;
                }
            } else if (mImposedDirection == 2) { // Z direction
                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->FastGetSolutionStepValue(DISPLACEMENT_Z) += mVelocity * DeltaTime;
                    // Save calculated velocity and reaction for print
                    it->GetValue(TARGET_STRESS_Z) = pTargetStressTable->GetValue(CurrentTime);
                    it->GetValue(REACTION_STRESS_Z) = ReactionStress;
                    it->GetValue(LOADING_VELOCITY_Z) = mVelocity;
                }
            } else { // Radial direction
                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                    double cos_theta = it->X()/external_radius;
                    double sin_theta = it->Y()/external_radius;
                    it->FastGetSolutionStepValue(DISPLACEMENT_X) += mVelocity * cos_theta * DeltaTime;
                    it->FastGetSolutionStepValue(DISPLACEMENT_Y) += mVelocity * sin_theta * DeltaTime;
                    // Save calculated velocity and reaction for print
                    it->GetValue(TARGET_STRESS_X) = pTargetStressTable->GetValue(CurrentTime) * cos_theta;
                    it->GetValue(TARGET_STRESS_Y) = pTargetStressTable->GetValue(CurrentTime) * sin_theta;
                    it->GetValue(REACTION_STRESS_X) = ReactionStress * cos_theta;
                    it->GetValue(REACTION_STRESS_Y) = ReactionStress * sin_theta;
                    it->GetValue(LOADING_VELOCITY_X) = mVelocity * cos_theta;
                    it->GetValue(LOADING_VELOCITY_Y) = mVelocity * sin_theta;
                }
            }
        } else {
            // Save calculated velocity and reaction for print
            if (mImposedDirection == 0) { // X direction
                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->GetValue(TARGET_STRESS_X) = pTargetStressTable->GetValue(CurrentTime);
                    it->GetValue(REACTION_STRESS_X) = ReactionStress;
                    it->GetValue(LOADING_VELOCITY_X) = 0.0;
                }
            } else if (mImposedDirection == 1) { // Y direction
                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->GetValue(TARGET_STRESS_Y) = pTargetStressTable->GetValue(CurrentTime);
                    it->GetValue(REACTION_STRESS_Y) = ReactionStress;
                    it->GetValue(LOADING_VELOCITY_Y) = 0.0;
                }
            } else if (mImposedDirection == 2) { // Z direction
                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->GetValue(TARGET_STRESS_Z) = pTargetStressTable->GetValue(CurrentTime);
                    it->GetValue(REACTION_STRESS_Z) = ReactionStress;
                    it->GetValue(LOADING_VELOCITY_Z) = 0.0;
                }
            } else { // Radial direction
                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                    double cos_theta = it->X()/external_radius;
                    double sin_theta = it->Y()/external_radius;
                    it->GetValue(TARGET_STRESS_X) = pTargetStressTable->GetValue(CurrentTime) * cos_theta;
                    it->GetValue(TARGET_STRESS_Y) = pTargetStressTable->GetValue(CurrentTime) * sin_theta;
                    it->GetValue(REACTION_STRESS_X) = ReactionStress * cos_theta;
                    it->GetValue(REACTION_STRESS_Y) = ReactionStress * sin_theta;
                    it->GetValue(LOADING_VELOCITY_X) = 0.0;
                    it->GetValue(LOADING_VELOCITY_Y) = 0.0;
                }
            }
        }

        KRATOS_CATCH("");
    }


    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override
    {
        // Update K with latest ReactionStress after the axis has been loaded
        if (mApplyCM == true) {
            if (mAlternateAxisLoading == true) {
                const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];
                double ReactionStress = CalculateReactionStress();
                if(mUpdateStiffness == true) {
                    mStiffness = EstimateStiffness(ReactionStress,delta_time);
                }
            }
        }
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ControlModuleProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ControlModuleProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
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
    bool mApplyCM;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ControlModuleProcess& operator=(ControlModuleProcess const& rOther);

    /// Copy constructor.
    //ControlModuleProcess(ControlModuleProcess const& rOther);

    double UpdateVectorOfHistoricalStressesAndComputeNewAverage(const double& last_reaction) {
        KRATOS_TRY;
        int length_of_vector = mVectorOfLastStresses.size();
        if (length_of_vector == 0) { //only the first time
            int number_of_steps_for_stress_averaging = (int) (mStressAveragingTime / mrModelPart.GetProcessInfo()[DELTA_TIME]);
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

    void IsTimeToApplyCM(){
        const double current_time = mrModelPart.GetProcessInfo()[TIME];
        mApplyCM = false;

        if(current_time >= mStartTime) {
            if (mAlternateAxisLoading == true) {
                const unsigned int step = mrModelPart.GetProcessInfo()[STEP];
                if (mImposedDirection == 0) {
                    if(step == mXCounter){
                        mApplyCM = true;
                        mXCounter += 3;
                    }
                } else if (mImposedDirection == 1) {
                    if(step == mYCounter){
                        mApplyCM = true;
                        mYCounter += 3;
                    }
                } else if (mImposedDirection == 2) {
                    if(step == mZCounter){
                        mApplyCM = true;
                        mZCounter += 3;
                    }
                } else {
                    mApplyCM = true;
                }
            } else {
                mApplyCM = true;
            }
        }
    }

    double CalculateReactionStress() {
        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        // Calculate face_area
        const int NCons = static_cast<int>(mrModelPart.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = mrModelPart.ConditionsBegin();
        double FaceArea = 0.0;
        #pragma omp parallel for reduction(+:FaceArea)
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            FaceArea += itCond->GetGeometry().Area();
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
        } else { // Radial direction
            #pragma omp parallel for reduction(+:FaceReaction)
            for(int i = 0; i<NNodes; i++){
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                // Unit normal vector pointing outwards
                array_1d<double,2> n;
                n[0] = it->X();
                n[1] = it->Y();
                double inv_norm = 1.0/norm_2(n);
                n[0] *= inv_norm;
                n[1] *= inv_norm;
                // Scalar product between reaction and normal
                double n_dot_r = n[0] * it->FastGetSolutionStepValue(REACTION_X) +
                                    n[1] * it->FastGetSolutionStepValue(REACTION_Y);
                FaceReaction += n_dot_r;
            }
        }

        double ReactionStress;
        if (std::abs(FaceArea) > 1.0e-12) {
            ReactionStress = FaceReaction/FaceArea;
        } else {
            ReactionStress = 0.0;
        }

        return ReactionStress;
    }

    double EstimateStiffness(const double& rReactionStress, const double& rDeltaTime) {
        double K_estimated = mStiffness;
        if(std::abs(mVelocity) > 1.0e-12 && std::abs(rReactionStress-mReactionStressOld) > mStressIncrementTolerance) {
            K_estimated = std::abs((rReactionStress-mReactionStressOld)/(mVelocity * rDeltaTime));
        }
        return K_estimated;
    }

}; // Class ControlModuleProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ControlModuleProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ControlModuleProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_CONTROL_MODULE_PROCESS defined */
