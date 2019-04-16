//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED)
#define K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED

// System includes
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_processes/scalar_co_solving_process.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "rans_modelling_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// The base class for all KEpsilonCoSolvingProcess in Kratos.
/** The KEpsilonCoSolvingProcess is the base class for all KEpsilonCoSolvingProcess and defines a simple interface for them.
    Execute method is used to execute the KEpsilonCoSolvingProcess algorithms. While the parameters of this method
  can be very different from one KEpsilonCoSolvingProcess to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all KEpsilonCoSolvingProcess parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other KEpsilonCoSolvingProcess or the base KEpsilonCoSolvingProcess class.
*/
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class KEpsilonCoSolvingProcess
    : public ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    typedef ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    /// Pointer definition of KEpsilonCoSolvingProcess
    KRATOS_CLASS_POINTER_DEFINITION(KEpsilonCoSolvingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    KEpsilonCoSolvingProcess(ModelPart& rModelPart, Parameters& rParameters, Process& rYPlusModelProcess)
        : BaseType(rModelPart, rParameters, TURBULENT_VISCOSITY),
          mrYPlusModelProcess(rYPlusModelProcess)
    {
    }

    /// Destructor.
    ~KEpsilonCoSolvingProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ExecuteFinalizeSolutionStep() override
    {
        UpdateEffectiveViscosity();
    }

    int Check() override
    {
        int value = BaseType::Check();

        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_VISCOSITY_MIN);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_VISCOSITY_MAX);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENCE_RANS_C_MU);
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);
        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);

        NodesContainerType& r_nodes = this->mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
        }

        return value;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KEpsilonCoSolvingProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KEpsilonCoSolvingProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Process& mrYPlusModelProcess;

    ///@}
    ///@name Operations
    ///@{

    void UpdateBeforeSolveSolutionStep() override
    {
        mrYPlusModelProcess.Execute();
    }

    void UpdateAfterSolveSolutionStep() override
    {
    }

    void UpdateConvergenceVariable() override
    {
        UpdateTurbulentViscosity();
    }

    void UpdateEffectiveViscosity()
    {
        NodesContainerType& r_nodes = this->mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            const double nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);

            r_node.FastGetSolutionStepValue(VISCOSITY) = nu_t + nu;
        }
    }

    void UpdateTurbulentViscosity()
    {
        const ProcessInfo& r_process_info = this->mrModelPart.GetProcessInfo();
        const double nu_t_min = r_process_info[TURBULENT_VISCOSITY_MIN];
        const double nu_t_max = r_process_info[TURBULENT_VISCOSITY_MAX];
        const double c_mu = r_process_info[TURBULENCE_RANS_C_MU];

        NodesContainerType& r_nodes = this->mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double epsilon =
                r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
            const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
                c_mu, tke, epsilon, f_mu);

            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t;
        }

        unsigned int lower_number_of_nodes, higher_number_of_nodes;
        RansVariableUtils().ClipScalarVariable(
            lower_number_of_nodes, higher_number_of_nodes, nu_t_min, nu_t_max,
            TURBULENT_VISCOSITY, this->mrModelPart.Nodes());

        KRATOS_WARNING_IF(this->Info(),
                          (lower_number_of_nodes > 0 || higher_number_of_nodes > 0) &&
                              this->mEchoLevel > 0)
            << "TURBULENT_VISCOSITY out of bounds. [ " << lower_number_of_nodes
            << " nodes < " << std::scientific << nu_t_min << ", "
            << higher_number_of_nodes << " nodes > " << std::scientific << nu_t_max
            << "  out of total number of " << number_of_nodes << " nodes. ]\n";
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(
        KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    // KEpsilonCoSolvingProcess(KEpsilonCoSolvingProcess const& rOther);

    ///@}

}; // Class KEpsilonCoSolvingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::istream& operator>>(std::istream& rIStream,
                                KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::ostream& operator<<(std::ostream& rOStream,
                                const KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED defined
