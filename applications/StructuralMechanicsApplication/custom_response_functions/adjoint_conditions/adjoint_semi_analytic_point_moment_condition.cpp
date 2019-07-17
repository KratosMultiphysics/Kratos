// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kevin Braun, https://github.com/MFusseder
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"

#include "adjoint_semi_analytic_point_moment_condition.h"
#include "structural_mechanics_application_variables.h"
#include "custom_conditions/point_moment_condition_3d.h"

namespace Kratos
{

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointMomentCondition<TPrimalCondition>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dim = 3;
        if (rResult.size() != dim * number_of_nodes)
            rResult.resize(dim*number_of_nodes,false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_ROTATION_X);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            const IndexType index = i * 3;
            rResult[index    ] = this->GetGeometry()[i].GetDof(ADJOINT_ROTATION_X,pos    ).EquationId();
            rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_ROTATION_Y,pos + 1).EquationId();
            rResult[index + 2] = this->GetGeometry()[i].GetDof(ADJOINT_ROTATION_Z,pos + 2).EquationId();
        }
        
        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointMomentCondition<TPrimalCondition>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  3;
        const SizeType num_dofs = number_of_nodes * dimension;

        if (rElementalDofList.size() != num_dofs)
            rElementalDofList.resize(num_dofs);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            const IndexType index = i * 3;
            rElementalDofList[index    ] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_X);
            rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Y);
            rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Z);
        }
        
        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointMomentCondition<TPrimalCondition>::GetValuesVector(Vector& rValues, int Step)
    {
        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = 3;
        const SizeType num_dofs = number_of_nodes * dimension;

        if (rValues.size() != num_dofs)
            rValues.resize(num_dofs, false);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            const array_1d<double, 3 > & ROTATION = this->GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_ROTATION, Step);
            IndexType index = i * dimension;
            for(IndexType k = 0; k < dimension; ++k)
                rValues[index + k] = ROTATION[k];
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointMomentCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  3;
        const SizeType num_dofs = number_of_nodes * dimension;
        rOutput = ZeroMatrix(number_of_nodes, num_dofs);

        KRATOS_CATCH( "" )
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointMomentCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        // not yet implemented
    }

    template <class TPrimalCondition>
    int AdjointSemiAnalyticPointMomentCondition<TPrimalCondition>::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // verify that the variables are correctly initialized
        KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);
        KRATOS_CHECK_VARIABLE_KEY(ROTATION);

        // Check dofs
        const GeometryType& r_geom = this->GetGeometry();
        for (IndexType i = 0; i < r_geom.size(); ++i)
        {
            const auto& r_node = r_geom[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION, r_node);

            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z, r_node);
        }

        return 0;

        KRATOS_CATCH( "" )
    }

    // TODO find out what to do with KRATOS_API
    template class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointSemiAnalyticPointMomentCondition<PointMomentCondition3D>;

} // Namespace Kratos


