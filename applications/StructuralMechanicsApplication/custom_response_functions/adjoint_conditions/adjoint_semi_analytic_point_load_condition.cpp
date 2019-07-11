// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"

#include "adjoint_semi_analytic_point_load_condition.h"
#include "structural_mechanics_application_variables.h"
#include "custom_conditions/point_load_condition.h"

namespace Kratos
{

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dim = this->GetGeometry().WorkingSpaceDimension();
        if (rResult.size() != dim * number_of_nodes)
            rResult.resize(dim*number_of_nodes,false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

        if(dim == 2)
        {
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                const IndexType index = i * 2;
                rResult[index    ] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
            }
        }
        else
        {
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                const IndexType index = i * 3;
                rResult[index    ] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z,pos + 2).EquationId();
            }
        }
        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = number_of_nodes * dimension;

        if (rElementalDofList.size() != num_dofs)
            rElementalDofList.resize(num_dofs);

        if(dimension == 2)
        {
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                const IndexType index = i * 2;
                rElementalDofList[index    ] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
                rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
            }
        }
        else
        {
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                const IndexType index = i * 3;
                rElementalDofList[index    ] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
                rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
                rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z);
            }
        }

        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::GetValuesVector(Vector& rValues, int Step)
    {
        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = number_of_nodes * dimension;

        if (rValues.size() != num_dofs)
            rValues.resize(num_dofs, false);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            const array_1d<double, 3 > & Displacement = this->GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
            IndexType index = i * dimension;
            for(IndexType k = 0; k < dimension; ++k)
                rValues[index + k] = Displacement[k];
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = number_of_nodes * dimension;
        rOutput = ZeroMatrix(0, num_dofs);

        KRATOS_CATCH( "" )
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

    template <class TPrimalCondition>
    int AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // verify that the variables are correctly initialized
        KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);

        // Check dofs
        const GeometryType& r_geom = this->GetGeometry();
        for (IndexType i = 0; i < r_geom.size(); ++i)
        {
            const auto& r_node = r_geom[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
        }

        return 0;

        KRATOS_CATCH( "" )
    }

    // TODO find out what to do with KRATOS_API
    template class AdjointSemiAnalyticPointLoadCondition<PointLoadCondition>;

} // Namespace Kratos


