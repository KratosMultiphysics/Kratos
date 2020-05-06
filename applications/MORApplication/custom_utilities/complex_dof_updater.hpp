#if !defined(KRATOS_COMPLEX_DOF_UPDATER)
#define KRATOS_COMPLEX_DOF_UPDATER

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
namespace ComplexDofUpdater
{

    void AssignComplexDofValue(
        Variable<double>& var,
        Variable<double>& var_real,
        Variable<double>& var_imag,
        Node<3>& node,
        const std::complex<double>& Z)
    {
        const int sign = std::real(Z) >= 0 ? 1 : -1;
        node.FastGetSolutionStepValue(var) = sign * std::abs(Z);
        node.FastGetSolutionStepValue(var_real) = std::real(Z);
        node.FastGetSolutionStepValue(var_imag) = std::imag(Z);
    }

#ifdef KRATOS_MOR_ASSIGN_COMPLEX_DOF
#undef KRATOS_MOR_ASSIGN_COMPLEX_DOF
#endif
#define KRATOS_MOR_ASSIGN_COMPLEX_DOF(name, it_node) \
    if( it_node->HasDofFor(name) ) \
    { \
        const size_t eq_id = it_node->GetDof(name).EquationId(); \
        AssignComplexDofValue(name, REAL_##name, IMAG_##name, (*it_node), rX(eq_id)); \
    }

    /** @brief Assign new values for the problem's degrees of freedom using the complex vector rX.
     *  @details value = std::abs(rX[dof.EquationId()])
     *      REAL_value = std::real(rX[dof.EquationId()])
     *      IMAG_value = std::imag(rX[dof.EquationId()])
     *  @param[in/out] rModelPart The model part.
     *  @param[in] rX The solution vector.
     */
    template<typename ComplexSystemVectorType>
    void AssignDofs(
        ModelPart& rModelPart,
        const ComplexSystemVectorType& rX)
    {
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );

        #pragma omp parallel for
        for( int i=0; i<num_nodes; ++i )
        {
            auto it_node = std::begin(rModelPart.Nodes()) + i;

            KRATOS_MOR_ASSIGN_COMPLEX_DOF(DISPLACEMENT_X, it_node);
            KRATOS_MOR_ASSIGN_COMPLEX_DOF(DISPLACEMENT_Y, it_node);
            KRATOS_MOR_ASSIGN_COMPLEX_DOF(DISPLACEMENT_Z, it_node);
            KRATOS_MOR_ASSIGN_COMPLEX_DOF(PRESSURE, it_node);
        }
    }

} // namespace ComplexDofUpdater

} // namespace Kratos

#endif