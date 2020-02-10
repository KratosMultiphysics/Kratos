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
                "variable_name": "VARIABLE_NAME",
                "reaction_variable_name": "REACTION_VARIABLE_NAME",
                "target_stress_variable_name": "TARGET_STRESS_VARIABLE_NAME",
                "reaction_stress_variable_name": "REACTION_STRESS_VARIABLE_NAME",
                "loading_velocity_variable_name": "LOADING_VELOCITY_VARIABLE_NAME",
                "radial_displacement" : false,
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

        mVariableName = rParameters["variable_name"].GetString();
        mReactionVariableName = rParameters["reaction_variable_name"].GetString();
        mTargetStressVariableName = rParameters["target_stress_variable_name"].GetString();
        mReactionStressVariableName = rParameters["reaction_stress_variable_name"].GetString();
        mLoadingVelocityVariableName = rParameters["loading_velocity_variable_name"].GetString();
        mTargetStressTableId = rParameters["target_stress_table_id"].GetInt();
        mVelocity = rParameters["initial_velocity"].GetDouble();
        mLimitVelocity = rParameters["limit_velocity"].GetDouble();
        mVelocityFactor = rParameters["velocity_factor"].GetDouble();
        mStartTime = rParameters["start_time"].GetDouble();
        mStressIncrementTolerance = rParameters["stress_increment_tolerance"].GetDouble();
        mUpdateStiffness = rParameters["update_stiffness"].GetBool();
        mReactionStressOld = 0.0;
        mStiffness = rParameters["young_modulus"].GetDouble()*rParameters["face_area"].GetDouble()/rParameters["compression_length"].GetDouble();

        mRadialDisplacement = rParameters["radial_displacement"].GetBool();
        if(mRadialDisplacement == true) {
            mVariableNameX = rParameters["variable_name"].GetString() + std::string("_X");
            mReactionVariableNameX = rParameters["reaction_variable_name"].GetString() + std::string("_X");
            mTargetStressVariableNameX = rParameters["target_stress_variable_name"].GetString() + std::string("_X");
            mReactionStressVariableNameX = rParameters["reaction_stress_variable_name"].GetString() + std::string("_X");
            mLoadingVelocityVariableNameX = rParameters["loading_velocity_variable_name"].GetString() + std::string("_X");
            mVariableNameY = rParameters["variable_name"].GetString() + std::string("_Y");
            mReactionVariableNameY = rParameters["reaction_variable_name"].GetString() + std::string("_Y");
            mTargetStressVariableNameY = rParameters["target_stress_variable_name"].GetString() + std::string("_Y");
            mReactionStressVariableNameY = rParameters["reaction_stress_variable_name"].GetString() + std::string("_Y");
            mLoadingVelocityVariableNameY = rParameters["loading_velocity_variable_name"].GetString() + std::string("_Y");
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
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;

        if(mRadialDisplacement == true) {
            ComponentType VarComponentX = KratosComponents< ComponentType >::Get(mVariableNameX);
            ComponentType VarComponentY = KratosComponents< ComponentType >::Get(mVariableNameY);

            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                it->Fix(VarComponentX);
                it->FastGetSolutionStepValue(VarComponentX) = 0.0;
                it->Fix(VarComponentY);
                it->FastGetSolutionStepValue(VarComponentY) = 0.0;
            }
        } else {
            ComponentType VarComponent = KratosComponents< ComponentType >::Get(mVariableName);

            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                it->Fix(VarComponent);
                it->FastGetSolutionStepValue(VarComponent) = 0.0;
            }
        }

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;
        const double DeltaTime = mrModelPart.GetProcessInfo()[DELTA_TIME];

        if(mRadialDisplacement == true) {
            ComponentType VarComponentX = KratosComponents< ComponentType >::Get(mVariableNameX);
            ComponentType VarComponentY = KratosComponents< ComponentType >::Get(mVariableNameY);

            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                double cos_theta = it->X()/external_radius;
                double sin_theta = it->Y()/external_radius;

                it->FastGetSolutionStepValue(VarComponentX) += mVelocity * cos_theta * DeltaTime;
                it->FastGetSolutionStepValue(VarComponentY) += mVelocity * sin_theta * DeltaTime;
            }
        } else {
            ComponentType VarComponent = KratosComponents< ComponentType >::Get(mVariableName);

            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                it->FastGetSolutionStepValue(VarComponent) += mVelocity * DeltaTime;
            }
        }

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;

        const double CurrentTime = mrModelPart.GetProcessInfo()[TIME];

        if(CurrentTime >= mStartTime && mTargetStressTableId > 0)
        {
            // Calculate ReactionStress
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
            const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
            double FaceReaction = 0.0;

            if(mRadialDisplacement == true) {
                ComponentType ReactionVarComponentX = KratosComponents< ComponentType >::Get(mReactionVariableNameX);
                ComponentType ReactionVarComponentY = KratosComponents< ComponentType >::Get(mReactionVariableNameY);

                #pragma omp parallel for reduction(+:FaceReaction)
                for(int i = 0; i<NNodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    // Unit normal vector pointing outwards
                    array_1d<double,2> n;
                    n[0] = it->X();
                    n[1] = it->Y();
                    double inv_norm = 1.0/norm_2(n);
                    n[0] *= inv_norm;
                    n[1] *= inv_norm;

                    // Scalar product between reaction and normal
                    double n_dot_r = n[0] * it->FastGetSolutionStepValue(ReactionVarComponentX) +
                                     n[1] * it->FastGetSolutionStepValue(ReactionVarComponentY);

                    FaceReaction += n_dot_r;
                }
            } else {
                ComponentType ReactionVarComponent = KratosComponents< ComponentType >::Get(mReactionVariableName);

                #pragma omp parallel for reduction(+:FaceReaction)
                for(int i = 0; i<NNodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    FaceReaction += it->FastGetSolutionStepValue(ReactionVarComponent);
                }
            }

            const double ReactionStress = FaceReaction/FaceArea;

            // Update K if required
            const double DeltaTime = mrModelPart.GetProcessInfo()[DELTA_TIME];
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
            TableType::Pointer pTargetStressTable = mrModelPart.pGetTable(mTargetStressTableId);
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

            if (mRadialDisplacement == true) {
                ComponentType TargetStressVarComponentX = KratosComponents< ComponentType >::Get(mTargetStressVariableNameX);
                ComponentType TargetStressVarComponentY = KratosComponents< ComponentType >::Get(mTargetStressVariableNameY);
                ComponentType ReactionStressVarComponentX = KratosComponents< ComponentType >::Get(mReactionStressVariableNameX);
                ComponentType ReactionStressVarComponentY = KratosComponents< ComponentType >::Get(mReactionStressVariableNameY);
                ComponentType LoadingVelocityVarComponentX = KratosComponents< ComponentType >::Get(mLoadingVelocityVariableNameX);
                ComponentType LoadingVelocityVarComponentY = KratosComponents< ComponentType >::Get(mLoadingVelocityVariableNameY);

                #pragma omp parallel for
                for(int i = 0; i<NNodes; i++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                    double cos_theta = it->X()/external_radius;
                    double sin_theta = it->Y()/external_radius;

                    it->FastGetSolutionStepValue(TargetStressVarComponentX) = pTargetStressTable->GetValue(CurrentTime) * cos_theta;
                    it->FastGetSolutionStepValue(TargetStressVarComponentY) = pTargetStressTable->GetValue(CurrentTime) * sin_theta;
                    it->FastGetSolutionStepValue(ReactionStressVarComponentX) = ReactionStress * cos_theta;
                    it->FastGetSolutionStepValue(ReactionStressVarComponentY) = ReactionStress * sin_theta;
                    it->FastGetSolutionStepValue(LoadingVelocityVarComponentX) = mVelocity * cos_theta;
                    it->FastGetSolutionStepValue(LoadingVelocityVarComponentY) = mVelocity * sin_theta;
                }
            } else {
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
            }
        }

        KRATOS_CATCH("");
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
    std::string mVariableName;
    std::string mReactionVariableName;
    std::string mTargetStressVariableName;
    std::string mReactionStressVariableName;
    std::string mLoadingVelocityVariableName;
    unsigned int mTargetStressTableId;
    double mVelocity;
    double mLimitVelocity;
    double mVelocityFactor;
    double mStartTime;
    double mReactionStressOld;
    double mStressIncrementTolerance;
    double mStiffness;
    bool mUpdateStiffness;

    bool mRadialDisplacement;
    std::string mVariableNameX;
    std::string mVariableNameY;
    std::string mReactionVariableNameX;
    std::string mReactionVariableNameY;
    std::string mTargetStressVariableNameX;
    std::string mTargetStressVariableNameY;
    std::string mReactionStressVariableNameX;
    std::string mReactionStressVariableNameY;
    std::string mLoadingVelocityVariableNameX;
    std::string mLoadingVelocityVariableNameY;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ControlModuleProcess& operator=(ControlModuleProcess const& rOther);

    /// Copy constructor.
    //ControlModuleProcess(ControlModuleProcess const& rOther);

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
