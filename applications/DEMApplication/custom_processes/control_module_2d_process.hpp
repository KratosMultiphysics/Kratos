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


#if !defined(KRATOS_CONTROL_MODULE_2D_PROCESS )
#define  KRATOS_CONTROL_MODULE_2D_PROCESS

// Project includes
#include "includes/table.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application includes
#include "custom_elements/spheric_continuum_particle.h"

#include "DEM_application_variables.h"

namespace Kratos
{

class ControlModule2DProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ControlModule2DProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ControlModule2DProcess(
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
                "control_module_time_step": 1.0e-5,
                "target_stress_table_id" : 0,
                "initial_velocity" : 0.0,
                "limit_velocity" : 0.1,
                "velocity_factor" : 1.0,
                "compression_length" : 1.0,
                "young_modulus" : 1.0e9,
                "stress_increment_tolerance": 100.0,
                "update_stiffness": true,
                "stiffness_alpha": 1.0,
                "start_time" : 0.0,
                "stress_averaging_time": 1e-7
            }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mImposedDirection = rParameters["imposed_direction"].GetInt();
        mTargetStressTableId = rParameters["target_stress_table_id"].GetInt();
        mVelocity = rParameters["initial_velocity"].GetDouble();
        mLimitVelocity = rParameters["limit_velocity"].GetDouble();
        mVelocityFactor = rParameters["velocity_factor"].GetDouble();
        mCompressionLength = rParameters["compression_length"].GetDouble();
        mStartTime = rParameters["start_time"].GetDouble();
        mStressIncrementTolerance = rParameters["stress_increment_tolerance"].GetDouble();
        mUpdateStiffness = rParameters["update_stiffness"].GetBool();
        mReactionStressOld = 0.0;
        mStiffness = rParameters["young_modulus"].GetDouble()/mCompressionLength; // mStiffness is actually a stiffness over an area
        mStiffnessAlpha = rParameters["stiffness_alpha"].GetDouble();
        mStressAveragingTime = rParameters["stress_averaging_time"].GetDouble();
        mVectorOfLastStresses.resize(0);

        // NOTE: Alternate axis loading only works for X,Y,Z loading. Radial loading is always applied.
        mAlternateAxisLoading = rParameters["alternate_axis_loading"].GetBool();
        mCMTimeStep = rParameters["control_module_time_step"].GetDouble();
        mStep = 0;

        // Initialize Variables
        mrModelPart.GetProcessInfo().SetValue(TARGET_STRESS_Z,0.0);
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
    ~ControlModule2DProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ControlModule2DProcess algorithms.
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
                it->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_X) = 0.0;
                it->FastGetSolutionStepValue(VELOCITY_X) = mVelocity;
                it->FastGetSolutionStepValue(DISPLACEMENT_Z) = 0.0;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Z) = 0.0;
                it->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
            }
        } else if (mImposedDirection == 1) {
            // Y direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Y) = 0.0;
                it->FastGetSolutionStepValue(VELOCITY_Y) = mVelocity;
                it->FastGetSolutionStepValue(DISPLACEMENT_Z) = 0.0;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Z) = 0.0;
                it->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
            }
        } else if (mImposedDirection == 2) {
            // Z direction
            mrModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE] = 0.0;
        } else {
            // Radial direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                it->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_X) = 0.0;
                it->FastGetSolutionStepValue(VELOCITY_X) = mVelocity * cos_theta;
                it->FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Y) = 0.0;
                it->FastGetSolutionStepValue(VELOCITY_Y) = mVelocity * sin_theta;
                it->FastGetSolutionStepValue(DISPLACEMENT_Z) = 0.0;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Z) = 0.0;
                it->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
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

        mStep += 1;
        if (mStep == 1) {
            unsigned int cm_num_super_steps = (int)(mCMTimeStep / mrModelPart.GetProcessInfo()[DELTA_TIME]);
            if(cm_num_super_steps < 1) {cm_num_super_steps = 1;}
            mXStartStep = 1;
            mXEndStep = mXStartStep + cm_num_super_steps - 1;
            mYStartStep = mXEndStep + 1;
            mYEndStep = mYStartStep + cm_num_super_steps - 1;
            mZStartStep = mYEndStep + 1;
            mZEndStep = mZStartStep + cm_num_super_steps - 1;
        }

        double ReactionStress = CalculateReactionStress();
        ReactionStress = UpdateVectorOfHistoricalStressesAndComputeNewAverage(ReactionStress);

        // Check whether this is a loading step for the current axis
        IsTimeToApplyCM();

        if (mApplyCM == true) {

            // Update K if required
            double delta_time;
            if (mAlternateAxisLoading == false) {
                delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];
                if(mUpdateStiffness == true) {
                    double K_estimated = EstimateStiffness(ReactionStress,delta_time);
                    mStiffness = mStiffnessAlpha * K_estimated + (1.0 - mStiffnessAlpha) * mStiffness;
                }
            } else {
                delta_time = mCMTimeStep;
                // delta_time = 3.0*mCMTimeStep;
            }

            // Update velocity
            const double NextTargetStress = pTargetStressTable->GetValue(CurrentTime+delta_time);
            const double df_target = NextTargetStress - ReactionStress;
            double delta_velocity = df_target/(mStiffness * delta_time) - mVelocity;
            if(std::abs(df_target) < mStressIncrementTolerance) { delta_velocity = -mVelocity; }
            if (mAlternateAxisLoading == false) {
                mReactionStressOld = ReactionStress;
                // v_new = vf*v_new + (1-vf)*v_old
                mVelocity += mVelocityFactor * delta_velocity;
                if(std::abs(mVelocity) > std::abs(mLimitVelocity)) {
                    if(mVelocity >= 0.0) { mVelocity = std::abs(mLimitVelocity); }
                    else { mVelocity = - std::abs(mLimitVelocity); }
                }
            } else if (mIsStartStep == true) {
                mReactionStressOld = ReactionStress;
                mVelocity += mVelocityFactor * delta_velocity;
                if(std::abs(mVelocity) > std::abs(mLimitVelocity)) {
                    if(mVelocity >= 0.0) { mVelocity = std::abs(mLimitVelocity); }
                    else { mVelocity = - std::abs(mLimitVelocity); }
                }
            }

            UpdateMotionAndSaveVariablesToPrint(mVelocity, pTargetStressTable->GetValue(CurrentTime), ReactionStress);
        }
        else {
            UpdateMotionAndSaveVariablesToPrint(0.0, pTargetStressTable->GetValue(CurrentTime), ReactionStress);
        }

        // UpdateMotionAndSaveVariablesToPrint(mVelocity, pTargetStressTable->GetValue(CurrentTime), ReactionStress);

        KRATOS_CATCH("");
    }

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override
    {
        // Update K with latest ReactionStress after the axis has been loaded
        if (mAlternateAxisLoading == true && mApplyCM == true && mIsEndStep == true) {
            double ReactionStress = CalculateReactionStress();
            if(mUpdateStiffness == true) {
                double K_estimated = EstimateStiffness(ReactionStress,mCMTimeStep);
                mStiffness = mStiffnessAlpha * K_estimated + (1.0 - mStiffnessAlpha) * mStiffness;
            }
        }
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ControlModule2DProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ControlModule2DProcess";
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
    double mCompressionLength;
    double mStartTime;
    double mReactionStressOld;
    double mStressIncrementTolerance;
    double mStiffness;
    bool mUpdateStiffness;
    std::vector<double> mVectorOfLastStresses;
    double mStressAveragingTime;
    bool mAlternateAxisLoading;
    unsigned int mXStartStep;
    unsigned int mXEndStep;
    unsigned int mYStartStep;
    unsigned int mYEndStep;
    unsigned int mZStartStep;
    unsigned int mZEndStep;
    bool mApplyCM;
    bool mIsStartStep;
    bool mIsEndStep;
    unsigned int mStep;
    double mCMTimeStep;
    double mStiffnessAlpha;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ControlModule2DProcess& operator=(ControlModule2DProcess const& rOther);

    /// Copy constructor.
    //ControlModule2DProcess(ControlModule2DProcess const& rOther);

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
        mIsStartStep = false;
        mIsEndStep = false;
        // TODO
        // double r = (double)(rand()/(RAND_MAX)) - 0.5;
        // unsigned int cm_num_super_steps = (int)(r + mCMTimeStep / mrModelPart.GetProcessInfo()[DELTA_TIME]);
        unsigned int cm_num_super_steps = (int)(mCMTimeStep / mrModelPart.GetProcessInfo()[DELTA_TIME]);

        if(current_time >= mStartTime) {
            if (mAlternateAxisLoading == true) {
                if (mImposedDirection == 0) {
                    if (mStep >= mXStartStep && mStep <= mXEndStep) {
                        mApplyCM = true;
                        if(mStep == mXStartStep){
                            mIsStartStep = true;
                        } else if (mStep == mXEndStep) {
                            mIsEndStep = true;
                        }
                        if (mIsEndStep == true) {
                            mXStartStep += 3*cm_num_super_steps;
                            mXEndStep += 3*cm_num_super_steps;
                        }
                    }
                } else if (mImposedDirection == 1) {
                    if (mStep >= mYStartStep && mStep <= mYEndStep) {
                        mApplyCM = true;
                        if(mStep == mYStartStep){
                            mIsStartStep = true;
                        } else if (mStep == mYEndStep) {
                            mIsEndStep = true;
                        }
                        if (mIsEndStep == true) {
                            mYStartStep += 3*cm_num_super_steps;
                            mYEndStep += 3*cm_num_super_steps;
                        }
                    }
                } else if (mImposedDirection == 2) {
                    if (mStep >= mZStartStep && mStep <= mZEndStep) {
                        mApplyCM = true;
                        if(mStep == mZStartStep){
                            mIsStartStep = true;
                        } else if (mStep == mZEndStep) {
                            mIsEndStep = true;
                        }
                        if (mIsEndStep == true) {
                            mZStartStep += 3*cm_num_super_steps;
                            mZEndStep += 3*cm_num_super_steps;
                        }
                    }
                } else {
                    mApplyCM = true;
                    mIsStartStep = true;
                    mIsEndStep = true;
                }
            } else {
                mApplyCM = true;
                mIsStartStep = true;
                mIsEndStep = true;
            }
        }
    }

    double CalculateReactionStress() {
        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        // Calculate face_area
        double face_area = 0.0;
        if (mImposedDirection == 2) { // Z direction
            ModelPart::ElementsContainerType& rElements = mrModelPart.GetCommunicator().LocalMesh().Elements();

            #pragma omp parallel for reduction(+:face_area)
            for (int i = 0; i < (int)rElements.size(); i++) {
                ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + i;
                Element* p_element = ptr_itElem->get();
                SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*>(p_element);
                const double radius = pDemElem->GetRadius();
                face_area += Globals::Pi*radius*radius;
            }
        } else { // X and Y directions
            const int NCons = static_cast<int>(mrModelPart.Conditions().size());
            ModelPart::ConditionsContainerType::iterator con_begin = mrModelPart.ConditionsBegin();
            #pragma omp parallel for reduction(+:face_area)
            for(int i = 0; i < NCons; i++) {
                ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
                face_area += itCond->GetGeometry().Area();
            }
        }

        // Calculate ReactionStress
        double FaceReaction = 0.0;
        if (mImposedDirection == 0) { // X direction
            #pragma omp parallel for reduction(+:FaceReaction)
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                FaceReaction -= it->FastGetSolutionStepValue(CONTACT_FORCES_X);
            }
        } else if (mImposedDirection == 1) { // Y direction
            #pragma omp parallel for reduction(+:FaceReaction)
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                FaceReaction -= it->FastGetSolutionStepValue(CONTACT_FORCES_Y);
            }
        } else if (mImposedDirection == 2) { // Z direction
            ModelPart::ElementsContainerType& rElements = mrModelPart.GetCommunicator().LocalMesh().Elements();

            #pragma omp parallel for reduction(+:FaceReaction)
            for (int i = 0; i < (int)rElements.size(); i++) {
                ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + i;
                Element* p_element = ptr_itElem->get();
                SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*>(p_element);
                BoundedMatrix<double, 3, 3> stress_tensor = ZeroMatrix(3,3);
                noalias(stress_tensor) = (*(pDemElem->mSymmStressTensor));
                const double radius = pDemElem->GetRadius();
                FaceReaction += stress_tensor(2,2) * Globals::Pi*radius*radius;
            }
        } else { // Radial direction
            #pragma omp parallel for reduction(+:FaceReaction)
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                // Unit normal vector pointing outwards
                array_1d<double,2> n;
                n[0] = it->X();
                n[1] = it->Y();
                double inv_norm = 1.0/norm_2(n);
                n[0] *= inv_norm;
                n[1] *= inv_norm;
                // Scalar product between reaction and normal
                double n_dot_r = n[0] * it->FastGetSolutionStepValue(CONTACT_FORCES_X) + n[1] * it->FastGetSolutionStepValue(CONTACT_FORCES_Y);
                FaceReaction -= n_dot_r;
            }
            // FaceReaction = 0.0;
            // for(int i = 0; i<NNodes; i++) {
            //     ModelPart::NodesContainerType::iterator it = it_begin + i;
            //     // Unit normal vector pointing outwards
            //     array_1d<double,2> n;
            //     n[0] = it->X();
            //     n[1] = it->Y();
            //     double inv_norm = 1.0/norm_2(n);
            //     n[0] *= inv_norm;
            //     n[1] *= inv_norm;
            //     // Scalar product between reaction and normal
            //     double n_dot_r = n[0] * it->FastGetSolutionStepValue(CONTACT_FORCES_X) + n[1] * it->FastGetSolutionStepValue(CONTACT_FORCES_Y);
            //     FaceReaction -= n_dot_r;
            // }
        }

        double ReactionStress;
        if (std::abs(face_area) > 1.0e-12) {
            ReactionStress = FaceReaction/face_area;
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

    void UpdateMotionAndSaveVariablesToPrint(const double rVelocity, const double& rTargetStress, const double& rReactionStress) {

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        // Update Imposed displacement
        if (mImposedDirection == 0) { // X direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(VELOCITY_X) = rVelocity;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_X) = it->FastGetSolutionStepValue(VELOCITY_X) * mrModelPart.GetProcessInfo()[DELTA_TIME];
                it->FastGetSolutionStepValue(DISPLACEMENT_X) += it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_X);
                it->X() = it->X0() + it->FastGetSolutionStepValue(DISPLACEMENT_X);
                // Save calculated velocity and reaction for print
                // TODO: print errors:
                // max_relative
                // max_absolute
                // average_relative
                // average_absolute
                it->GetValue(TARGET_STRESS_X) = rTargetStress;
                it->GetValue(REACTION_STRESS_X) = rReactionStress;
                it->GetValue(LOADING_VELOCITY_X) = rVelocity;
            }
        } else if (mImposedDirection == 1) { // Y direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(VELOCITY_Y) = rVelocity;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Y) = it->FastGetSolutionStepValue(VELOCITY_Y) * mrModelPart.GetProcessInfo()[DELTA_TIME];
                it->FastGetSolutionStepValue(DISPLACEMENT_Y) += it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Y);
                it->Y() = it->Y0() + it->FastGetSolutionStepValue(DISPLACEMENT_Y);
                // Save calculated velocity and reaction for print
                it->GetValue(TARGET_STRESS_Y) = rTargetStress;
                it->GetValue(REACTION_STRESS_Y) = rReactionStress;
                it->GetValue(LOADING_VELOCITY_Y) = rVelocity;
            }
        } else if (mImposedDirection == 2) { // Z direction
            mrModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE] += rVelocity*mrModelPart.GetProcessInfo()[DELTA_TIME]/mCompressionLength;
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                // Save calculated velocity and reaction for print
                it->GetValue(TARGET_STRESS_Z) = rTargetStress;
                it->GetValue(REACTION_STRESS_Z) = rReactionStress;
                it->GetValue(LOADING_VELOCITY_Z) = rVelocity;
            }
            mrModelPart.GetProcessInfo()[TARGET_STRESS_Z] = std::abs(rTargetStress);
        } else { // Radial direction
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                it->FastGetSolutionStepValue(VELOCITY_X) = rVelocity * cos_theta;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_X) = it->FastGetSolutionStepValue(VELOCITY_X) * mrModelPart.GetProcessInfo()[DELTA_TIME];
                it->FastGetSolutionStepValue(DISPLACEMENT_X) += it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_X);
                it->X() = it->X0() + it->FastGetSolutionStepValue(DISPLACEMENT_X);
                it->FastGetSolutionStepValue(VELOCITY_Y) = rVelocity * sin_theta;
                it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Y) = it->FastGetSolutionStepValue(VELOCITY_Y) * mrModelPart.GetProcessInfo()[DELTA_TIME];
                it->FastGetSolutionStepValue(DISPLACEMENT_Y) += it->FastGetSolutionStepValue(DELTA_DISPLACEMENT_Y);
                it->Y() = it->Y0() + it->FastGetSolutionStepValue(DISPLACEMENT_Y);
                // Save calculated velocity and reaction for print
                it->GetValue(TARGET_STRESS_X) = rTargetStress * cos_theta;
                it->GetValue(TARGET_STRESS_Y) = rTargetStress * sin_theta;
                it->GetValue(REACTION_STRESS_X) = rReactionStress * cos_theta;
                it->GetValue(REACTION_STRESS_Y) = rReactionStress * sin_theta;
                it->GetValue(LOADING_VELOCITY_X) = rVelocity * cos_theta;
                it->GetValue(LOADING_VELOCITY_Y) = rVelocity * sin_theta;
            }
        }
    }


}; // Class ControlModule2DProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ControlModule2DProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ControlModule2DProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_CONTROL_MODULE_2D_PROCESS defined */
