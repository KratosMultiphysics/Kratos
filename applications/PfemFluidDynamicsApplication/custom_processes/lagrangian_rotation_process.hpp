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

#if !defined(KRATOS_LAGRANGIAN_ROTATION_PROCESS )
#define  KRATOS_LAGRANGIAN_ROTATION_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/math_utils.h"

#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

/// Process used to rotate lagrangian model parts using Rodrigues' rotation formula

class LagrangianRotationProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(LagrangianRotationProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    LagrangianRotationProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "VELOCITY",
                "is_fixed": false,
                "angular_velocity": 0.0,
                "rotation_axis_initial_point": [0.0,0.0,0.0],
                "rotation_axis_final_point": [0.0,0.0,1.0],
                "initial_time": 0.0
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();
        mis_fixed = rParameters["is_fixed"].GetBool();

        mangular_velocity = rParameters["angular_velocity"].GetDouble();
        maxis_initial_point = rParameters["rotation_axis_initial_point"].GetVector();
        maxis_final_point = rParameters["rotation_axis_final_point"].GetVector();
        minitial_time = rParameters["initial_time"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~LagrangianRotationProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyComponentTableProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        // Constructing auxiliary Rodrigues matrices
        array_1d<double,3> = rotation_axis;
        noalias(rotation_axis) = maxis_final_point - maxis_initial_point;
        const double axis_norm = norm_2(rotation_axis);
        if (axis_norm > 1.0e-15){
            rotation_axis[0] *= 1.0/axis_norm;
            rotation_axis[1] *= 1.0/axis_norm;
            rotation_axis[2] *= 1.0/axis_norm;
        }

        midentity_matrix = ZeroMatrix(3,3);
        midentity_matrix(0,0) = 1.0;
        midentity_matrix(1,1) = 1.0;
        midentity_matrix(2,2) = 1.0;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                maxis_matrix(i,j) = rotation_axis[i] * rotation_axis[j];
            }
        }

        mantisym_axis_matrix = ZeroMatrix(3, 3);
        mantisym_axis_matrix(0, 1) = - rotation_axis[2];
        mantisym_axis_matrix(0, 2) =   rotation_axis[1];
        mantisym_axis_matrix(1, 0) =   rotation_axis[2];
        mantisym_axis_matrix(1, 2) = - rotation_axis[0];
        mantisym_axis_matrix(2, 0) = - rotation_axis[1];
        mantisym_axis_matrix(2, 1) =   rotation_axis[0];

        // Fixing DOFs and perform initial rotation if necessary

        //TODO: seguir

        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        component_type var_component = KratosComponents< component_type >::Get(mvariable_name);

        const int nnodes = static_cast<int>(mr_model_part.Nodes().size());

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                if(mis_fixed)
                {
                    it->Fix(var_component);
                }

                it->FastGetSolutionStepValue(var_component) = minitial_value;
            }
        }


    CalculateRodriguesMatrices(time);

    const int nnodes = r_model_part.Nodes().size();

    if (nnodes > 0){
        auto it_begin = r_model_part.NodesBegin();
        array_1d<double, 3> P0;
        array_1d<double, 3> P;

        #pragma omp parallel for private(P0, P)
        for (int i = 0; i < nnodes; ++i){
            auto it = it_begin + i;
            RotateNode(*(it.base()), P0, P);
            array_1d<double, 3>& displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
            noalias(displacement) = P - P0;
            array_1d<double, 3>& velocity = it->FastGetSolutionStepValue(MESH_VELOCITY);
            noalias(velocity) = prod(mRp, P0 - mAInit);
        }
    }


        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        component_type var_component = KratosComponents< component_type >::Get(mvariable_name);

        const double Time = mr_model_part.GetProcessInfo()[TIME]/mTimeUnitConverter;
        double value = mpTable->GetValue(Time);

        const int nnodes = static_cast<int>(mr_model_part.Nodes().size());

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                it->FastGetSolutionStepValue(var_component) = value;
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "LagrangianRotationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LagrangianRotationProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    bool mis_fixed;
    std::string mvariable_name;

    double mangular_velocity;
    double minitial_time;
    array_1d<double, 3> maxis_initial_point;
    array_1d<double, 3> maxis_final_point;
    BoundedMatrix<double, 3, 3> midentity_matrix;
    BoundedMatrix<double, 3, 3> maxis_matrix;
    BoundedMatrix<double, 3, 3> mantisym_axis_matrix;
    BoundedMatrix<double, 3, 3> mrotation_dt_matrix;
    // BoundedMatrix<double, 3, 3> mrotation_matrix;



///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    void CalculateRodriguesMatrices(const double current_time)
    {
        const double delta_time = current_time - minitial_time;
        const double sin_theta = std::sin(delta_time * mangular_velocity);
        const double cos_theta = std::cos(delta_time * mangular_velocity);

        // Rotation matrix
        // noalias(mrotation_matrix) = cos_theta * midentity_matrix
        //                             + sin_theta * mantisym_axis_matrix
        //                             + (1.0 - cos_theta) * maxis_matrix;

        // Rotation matrix derivative (derivative of R with respect to time)
        noalias(mrotation_dt_matrix) = - mangular_velocity * sin_theta * midentity_matrix
                                       + mangular_velocity * cos_theta * mantisym_axis_matrix
                                       + mangular_velocity * sin_theta * maxis_matrix;
    }

    /// Assignment operator.
    LagrangianRotationProcess& operator=(LagrangianRotationProcess const& rOther);

    /// Copy constructor.
    //LagrangianRotationProcess(LagrangianRotationProcess const& rOther);

}; // Class LagrangianRotationProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  LagrangianRotationProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const LagrangianRotationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_LAGRANGIAN_ROTATION_PROCESS defined */
