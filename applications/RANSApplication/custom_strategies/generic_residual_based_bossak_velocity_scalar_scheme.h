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
#include "includes/checks.h"
#include "includes/model_part.h"
#include "residual_based_bossak_velocity_scheme.h"

// Application includes
#include "rans_application_variables.h"

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

    using BaseType = ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>;

    using NodeType = ModelPart::NodeType;

    using SystemMatrixType = typename BaseType::SystemMatrixType;

    using SystemVectorType = typename BaseType::SystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    /// Constructor.

    GenericResidualBasedBossakVelocityScalarScheme(const double AlphaBossak,
                                                   const double RelaxationFactor,
                                                   const Variable<double>& rScalarVariable,
                                                   const Variable<double>& rScalarRateVariable,
                                                   const Variable<double>& rRelaxedScalarRateVariable)
        : ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(
              AlphaBossak, RelaxationFactor, {}, {&rScalarVariable}, {&rScalarRateVariable}, {}, {}, {}),
          mpScalarRateVariable(&rScalarRateVariable),
          mpRelaxedScalarRateVariable(&rRelaxedScalarRateVariable)
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
            const double scalar_rate_dot_old =
                r_node.FastGetSolutionStepValue(*mpScalarRateVariable, 1);
            const double scalar_rate_dot =
                r_node.FastGetSolutionStepValue(*mpScalarRateVariable, 0);

            r_node.FastGetSolutionStepValue(*mpRelaxedScalarRateVariable) =
                this->mAlphaBossak * scalar_rate_dot_old +
                (1.0 - this->mAlphaBossak) * scalar_rate_dot;
        }

        KRATOS_CATCH("");
    }

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        int value = BaseType::Check(rModelPart);

        const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int iNode = 0; iNode < number_of_nodes; ++iNode)
        {
            NodeType& r_node = *(rModelPart.NodesBegin() + iNode);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*mpScalarRateVariable), r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*mpScalarRateVariable), r_node);
        }

        return value;
        KRATOS_CATCH("");
    }

private:
    Variable<double> const *mpScalarRateVariable;
    Variable<double> const *mpRelaxedScalarRateVariable;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_GENERIC_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED defined