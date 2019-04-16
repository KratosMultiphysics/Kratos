//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_GENERIC_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED)
#define KRATOS_GENERIC_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED

// System includes

// Project includes
#include "includes/model_part.h"
#include "includes/checks.h"
#include "solving_strategies/schemes/residual_based_bossak_velocity_scheme.h"

// Application includes
#include "rans_modelling_application_variables.h"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace>
class GenericResidualBasedBossakVelocityScalarScheme
    : public ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GenericResidualBasedBossakVelocityScalarScheme);

    typedef Node<3> NodeType;

    typedef ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::SystemMatrixType SystemMatrixType;

    typedef typename BaseType::SystemVectorType SystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    /// Constructor.

    GenericResidualBasedBossakVelocityScalarScheme(const double AlphaBossak,
                                            const Variable<double>& rScalarVariable,
                                            const Variable<double>& rScalarRateVariable,
                                            const Variable<double>& rRelaxedScalarRateVariable)
        : ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(AlphaBossak, {}, {&rScalarVariable}, {&rScalarRateVariable}, {}, {}, {}),
        mrScalarRateVariable(rScalarRateVariable),
        mrRelaxedScalarRateVariable(rRelaxedScalarRateVariable)
    {
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::Update(rModelPart, rDofSet, rA, rDx, rb);

        // Updating the auxiliary variables
        const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int iNode = 0; iNode < number_of_nodes; ++iNode)
        {
            NodeType& r_node = *(rModelPart.NodesBegin() + iNode);
            const double scalar_rate_dot_old = r_node.FastGetSolutionStepValue(
                mrScalarRateVariable, 1);
            const double scalar_rate_dot = r_node.FastGetSolutionStepValue(
                mrScalarRateVariable, 0);

            r_node.FastGetSolutionStepValue(mrRelaxedScalarRateVariable) =
                this->mAlphaBossak * scalar_rate_dot_old + (1.0 - this->mAlphaBossak) * scalar_rate_dot;
        }

        KRATOS_CATCH("");
    }

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        int value = BaseType::Check(rModelPart);

        KRATOS_CHECK_VARIABLE_KEY(mrScalarRateVariable);
        KRATOS_CHECK_VARIABLE_KEY(mrRelaxedScalarRateVariable);

        const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int iNode = 0; iNode < number_of_nodes; ++iNode)
        {
            NodeType& r_node = *(rModelPart.NodesBegin() + iNode);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrScalarRateVariable, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrRelaxedScalarRateVariable, r_node);
        }


        return value;
        KRATOS_CATCH("");
    }

private:
    Variable<double> const mrScalarRateVariable;
    Variable<double> const mrRelaxedScalarRateVariable;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_GENERIC_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED defined