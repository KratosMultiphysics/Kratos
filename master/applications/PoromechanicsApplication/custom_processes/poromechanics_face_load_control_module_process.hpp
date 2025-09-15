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


#if !defined(KRATOS_POROMECHANICS_FACE_LOAD_CONTROL_MODULE_PROCESS )
#define  KRATOS_POROMECHANICS_FACE_LOAD_CONTROL_MODULE_PROCESS

#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

/// FaceLoad control module for displacements

class PoromechanicsFaceLoadControlModuleProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PoromechanicsFaceLoadControlModuleProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    PoromechanicsFaceLoadControlModuleProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process(Flags()) , mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "imposed_direction" : 0,
                "face_load_table_id" : 0,
                "initial_velocity" : 0.0,
                "limit_velocity" : 1.0,
                "velocity_factor" : 1.0,
                "initial_stiffness" : 1.0e7,
                "force_increment_tolerance": 100.0,
                "update_stiffness": true,
                "force_averaging_time": 1.0e-5
            }  )" );

        rParameters.ValidateAndAssignDefaults(default_parameters);

        double ImposedDirection = rParameters["imposed_direction"].GetInt();

        if(ImposedDirection == 2)
        {
            mDisplacementName = "DISPLACEMENT" + std::string("_Z");
            mReactionName = "REACTION" + std::string("_Z");
            mAverageReactionName = "AVERAGE_REACTION" + std::string("_Z");
            mTargetReactionName = "TARGET_REACTION" + std::string("_Z");
            mLoadingVelocityName = "LOADING_VELOCITY" + std::string("_Z");
        }
        else if(ImposedDirection == 1)
        {
            mDisplacementName = "DISPLACEMENT" + std::string("_Y");
            mReactionName = "REACTION" + std::string("_Y");
            mAverageReactionName = "AVERAGE_REACTION" + std::string("_Y");
            mTargetReactionName = "TARGET_REACTION" + std::string("_Y");
            mLoadingVelocityName = "LOADING_VELOCITY" + std::string("_Y");
        }
        else
        {
            mDisplacementName = "DISPLACEMENT" + std::string("_X");
            mReactionName = "REACTION" + std::string("_X");
            mAverageReactionName = "AVERAGE_REACTION" + std::string("_X");
            mTargetReactionName = "TARGET_REACTION" + std::string("_X");
            mLoadingVelocityName = "LOADING_VELOCITY" + std::string("_X");
        }

        mFaceLoadTableId = rParameters["face_load_table_id"].GetInt();
        mVelocity = rParameters["initial_velocity"].GetDouble();
        mLimitVelocity = rParameters["limit_velocity"].GetDouble();
        mVelocityFactor = rParameters["velocity_factor"].GetDouble();
        mForceIncrementTolerance = rParameters["force_increment_tolerance"].GetDouble();
        mUpdateStiffness = rParameters["update_stiffness"].GetBool();
        mReactionForceOld = 0.0;
        mStiffness = rParameters["initial_stiffness"].GetDouble();
        mForceAveragingTime = rParameters["force_averaging_time"].GetDouble(); // If set equal to the time step, no average is performed
        mVectorOfLastForces.resize(0);

        // Initialize Variables
        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
        array_1d<double,3> zero_vector = ZeroVector(3);
        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->SetValue(AVERAGE_REACTION,zero_vector);
            it->SetValue(TARGET_REACTION,zero_vector);
            it->SetValue(LOADING_VELOCITY,zero_vector);
        }

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~PoromechanicsFaceLoadControlModuleProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the PoromechanicsFaceLoadControlModuleProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        typedef Variable<double> component_type;
        const component_type& disp_component = KratosComponents< component_type >::Get(mDisplacementName);
        const component_type& loading_vel_component = KratosComponents< component_type >::Get(mLoadingVelocityName);

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->Fix(disp_component);
            it->FastGetSolutionStepValue(disp_component) = 0.0;
            it->GetValue(loading_vel_component) = mVelocity;
        }

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const double CurrentTime = mrModelPart.GetProcessInfo()[TIME];

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        double ReactionForce = CalculateReactionForce();
        ReactionForce = UpdateVectorOfHistoricalForcesAndComputeNewAverage(ReactionForce);

        // Update K if required
        const double DeltaTime = mrModelPart.GetProcessInfo()[DELTA_TIME];
        if(mUpdateStiffness == true) {
            mStiffness = EstimateStiffness(ReactionForce,DeltaTime);
        }
        mReactionForceOld = ReactionForce;

        // Update velocity
        const double NextTargetReaction = this->GetTargetReaction(CurrentTime+DeltaTime);
        const double df_target = NextTargetReaction - ReactionForce;
        double delta_velocity = df_target/(mStiffness * DeltaTime) - mVelocity;

        if(std::abs(df_target) < mForceIncrementTolerance) {
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

        typedef Variable<double> component_type;
        const component_type& disp_component = KratosComponents< component_type >::Get(mDisplacementName);

        // Update Imposed displacement
        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->FastGetSolutionStepValue(disp_component) += mVelocity * DeltaTime;
        }

        KRATOS_CATCH("");
    }


    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override
    {
        const double CurrentTime = mrModelPart.GetProcessInfo()[TIME];

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        double ReactionForce = CalculateReactionForce();

        typedef Variable<double> component_type;
        const component_type& target_reaction_component = KratosComponents< component_type >::Get(mTargetReactionName);
        const component_type& average_reaction_component = KratosComponents< component_type >::Get(mAverageReactionName);
        const component_type& loading_vel_component = KratosComponents< component_type >::Get(mLoadingVelocityName);

        // Update Imposed displacement
        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            // Save calculated velocity and reaction for print
            it->GetValue(target_reaction_component) = this->GetTargetReaction(CurrentTime);
            it->GetValue(average_reaction_component) = ReactionForce;
            it->GetValue(loading_vel_component) = mVelocity;
        }
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PoromechanicsFaceLoadControlModuleProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PoromechanicsFaceLoadControlModuleProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    // unsigned int mImposedDirection;
    std::string mDisplacementName;
    std::string mReactionName;
    std::string mAverageReactionName;
    std::string mTargetReactionName;
    std::string mLoadingVelocityName;
    unsigned int mFaceLoadTableId;
    double mVelocity;
    double mLimitVelocity;
    double mVelocityFactor;
    double mReactionForceOld;
    double mForceIncrementTolerance;
    double mStiffness;
    bool mUpdateStiffness;
    std::vector<double> mVectorOfLastForces;
    double mForceAveragingTime;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual double GetTargetReaction(const double& time) {

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

        // Get Target Face Load from table
        TableType::Pointer pFaceLoadTable = mrModelPart.pGetTable(mFaceLoadTableId);
        const double face_load = pFaceLoadTable->GetValue(time);

        return face_load*FaceArea;
    }

    virtual double UpdateVectorOfHistoricalForcesAndComputeNewAverage(const double& last_reaction) {
        KRATOS_TRY;
        int length_of_vector = mVectorOfLastForces.size();
        if (length_of_vector == 0) { //only the first time
            int number_of_steps_for_force_averaging = (int) (mForceAveragingTime / mrModelPart.GetProcessInfo()[DELTA_TIME]);
            if(number_of_steps_for_force_averaging < 1) number_of_steps_for_force_averaging = 1;
            mVectorOfLastForces.resize(number_of_steps_for_force_averaging);
            KRATOS_INFO("Poro control module") << " 'number_of_steps_for_force_averaging' is "<< number_of_steps_for_force_averaging << std::endl;
        }

        length_of_vector = mVectorOfLastForces.size();

        if(length_of_vector > 1) {
            for(int i=1; i<length_of_vector; i++) {
                mVectorOfLastForces[i-1] = mVectorOfLastForces[i];
            }
        }
        mVectorOfLastForces[length_of_vector-1] = last_reaction;

        double average = 0.0;
        for(int i=0; i<length_of_vector; i++) {
            average += mVectorOfLastForces[i];
        }
        average /= (double) length_of_vector;
        return average;

        KRATOS_CATCH("");
    }

    virtual double CalculateReactionForce() {

        typedef Variable<double> component_type;
        const component_type& reaction_component = KratosComponents< component_type >::Get(mReactionName);

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        // Calculate Average ReactionForce
        double FaceReaction = 0.0;
        #pragma omp parallel for reduction(+:FaceReaction)
        for(int i = 0; i<NNodes; i++){
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            FaceReaction += it->FastGetSolutionStepValue(reaction_component);
        }

        return FaceReaction/NNodes;
    }

    virtual double EstimateStiffness(const double& rReactionForce, const double& rDeltaTime) {
        double K_estimated = mStiffness;
        if(std::abs(mVelocity) > 1.0e-12 && std::abs(rReactionForce-mReactionForceOld) > mForceIncrementTolerance) {
            K_estimated = std::abs((rReactionForce-mReactionForceOld)/(mVelocity * rDeltaTime));
        }
        return K_estimated;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    PoromechanicsFaceLoadControlModuleProcess& operator=(PoromechanicsFaceLoadControlModuleProcess const& rOther);

    /// Copy constructor.
    //PoromechanicsFaceLoadControlModuleProcess(PoromechanicsFaceLoadControlModuleProcess const& rOther);



}; // Class PoromechanicsFaceLoadControlModuleProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PoromechanicsFaceLoadControlModuleProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PoromechanicsFaceLoadControlModuleProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_POROMECHANICS_FACE_LOAD_CONTROL_MODULE_PROCESS defined */
