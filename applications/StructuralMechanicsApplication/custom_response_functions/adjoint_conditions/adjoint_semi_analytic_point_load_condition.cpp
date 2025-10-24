// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"

#include "adjoint_semi_analytic_point_load_condition.h"
#include "structural_mechanics_application_variables.h"
#include "includes/variables.h"
#include "custom_conditions/point_load_condition.h"
#include "utilities/indirect_scalar.h"

namespace Kratos
{
    template <class TPrimalCondition>
    AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::ThisExtensions::ThisExtensions(Condition* pCondition)
        : mpCondition{pCondition}
    {
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::ThisExtensions::GetFirstDerivativesVector(
        std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
    {
        
        auto& r_node = mpCondition->GetGeometry()[NodeId];
        const SizeType dimension = mpCondition->GetGeometry().WorkingSpaceDimension();
        //const SizeType num_dofs = (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) ?  2 * dimension : dimension; // *2 for rotation
        rVector.resize(dimension );
        std::size_t index = 0;

        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_X, Step);
        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_Y, Step);
        if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_Z, Step);

    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::ThisExtensions::GetSecondDerivativesVector(
        std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
    {
        auto& r_node = mpCondition->GetGeometry()[NodeId];
        const SizeType dimension = mpCondition->GetGeometry().WorkingSpaceDimension();
        //const SizeType num_dofs = (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) ?  2 * dimension : dimension; // *2 for rotation
        rVector.resize(dimension );
        std::size_t index = 0;

        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_X, Step);
        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_Y, Step);
        if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_Z, Step);

    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::ThisExtensions::GetAuxiliaryVector(
        std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
    {
        auto& r_node = mpCondition->GetGeometry()[NodeId];
        const SizeType dimension = mpCondition->GetGeometry().WorkingSpaceDimension();
        //const SizeType num_dofs = (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) ?  2 * dimension : dimension; // *2 for rotation
        rVector.resize(dimension ); 
        std::size_t index = 0;

        rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_X, Step);
        rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_Y, Step);
        if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_Z, Step);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::ThisExtensions::GetFirstDerivativesVariables(
        std::vector<VariableData const*>& rVariables) const
    {
        rVariables.resize(1);
        rVariables[0] = &ADJOINT_VECTOR_2;
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::ThisExtensions::GetSecondDerivativesVariables(
        std::vector<VariableData const*>& rVariables) const
    {
        rVariables.resize(1);
        rVariables[0] = &ADJOINT_VECTOR_3;
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::ThisExtensions::GetAuxiliaryVariables(
        std::vector<VariableData const*>& rVariables) const
    {
        rVariables.resize(1);
        rVariables[0] = &AUX_ADJOINT_VECTOR_1;
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
        const SizeType mat_size = number_of_nodes * dimension;

        if( rDesignVariable == POINT_LOAD )
        {
            if ((rOutput.size1() != mat_size) || (rOutput.size2() != mat_size))
                rOutput.resize(mat_size, mat_size, false);

            noalias(rOutput) = ZeroMatrix(mat_size,mat_size);
            for(IndexType i = 0; i < mat_size; ++i)
                rOutput(i,i) = 1.0;
        }
        else if( rDesignVariable == SHAPE_SENSITIVITY )
        {
            rOutput = ZeroMatrix(mat_size, mat_size);
        }
        else
        {
            rOutput = ZeroMatrix(0, mat_size);
        }

        KRATOS_CATCH( "" )
    }

    // TODO find out what to do with KRATOS_API
    template class AdjointSemiAnalyticPointLoadCondition<PointLoadCondition>;

} // Namespace Kratos


