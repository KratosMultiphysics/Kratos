// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Mahmoud Sesa
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"

#include "adjoint_semi_analytic_surface_load_condition.h"
#include "structural_mechanics_application_variables.h"
#include "custom_conditions/surface_load_condition_3d.h"

namespace Kratos
{

    template <class TPrimalCondition>
    void AdjointSemiAnalyticSurfaceLoadCondition<TPrimalCondition>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes,false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            const IndexType index = i * 3;
            rResult[index    ] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[index + 2] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z,pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticSurfaceLoadCondition<TPrimalCondition>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = number_of_nodes * dimension;

        if (rElementalDofList.size() != num_dofs)
            rElementalDofList.resize(num_dofs);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            const IndexType index = i * 3;
            rElementalDofList[index] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
            rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
            rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z);
        }
        KRATOS_CATCH("")
    }

    // TODO Mahmoud: Segmentation fault occasionally happens here, the reason is not clear
    template <class TPrimalCondition>
    void AdjointSemiAnalyticSurfaceLoadCondition<TPrimalCondition>::GetValuesVector(Vector& rValues, int Step)
    {
        KRATOS_TRY

        const GeometryType & geom = this->GetGeometry();
        const SizeType number_of_nodes = geom.PointsNumber();
        const SizeType num_dofs = number_of_nodes * 3;

        if (rValues.size() != num_dofs)
            rValues.resize(num_dofs, false);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& Displacement = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);

            const IndexType index = i * 3;
            rValues[index]     = Displacement[0];
            rValues[index + 1] = Displacement[1];
            rValues[index + 2] = Displacement[2];
        }

         KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticSurfaceLoadCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        ProcessInfo process_info = rCurrentProcessInfo;
        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
        const SizeType vec_size = number_of_nodes * dimension;
        double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
        rOutput = ZeroMatrix(vec_size, vec_size);
        Vector RHS;

        this->CalculateRightHandSide(RHS, process_info);

        int i_2 = 0;
        for (auto& node_i : this->GetGeometry())
        {
            Vector perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of x
            node_i.X0() += delta;
            node_i.X()  += delta;
            this->CalculateRightHandSide(perturbed_RHS, process_info);
            row(rOutput, i_2) = (perturbed_RHS - RHS) / delta;
            node_i.X0() -= delta;
            node_i.X()  -= delta;

            // Reset the pertubed vector
            perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of y
            node_i.Y0() += delta;
            node_i.Y()  += delta;
            this->CalculateRightHandSide(perturbed_RHS, process_info);
            row(rOutput, i_2 + 1) = (perturbed_RHS - RHS) / delta;
            node_i.Y0()-= delta;
            node_i.Y() -= delta;

            // Reset the pertubed vector
            perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of z
            node_i.Z0() += delta;
            node_i.Z()  += delta;
            this->CalculateRightHandSide(perturbed_RHS, process_info);
            row(rOutput, i_2 + 2) = (perturbed_RHS - RHS) / delta;
            node_i.Z0() -= delta;
            node_i.Z()  -= delta;

            i_2 += 3;
        }

        KRATOS_CATCH( "" )
    }

    template <class TPrimalCondition>
    int AdjointSemiAnalyticSurfaceLoadCondition<TPrimalCondition>::Check( const ProcessInfo& rCurrentProcessInfo )
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
    template class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointSemiAnalyticSurfaceLoadCondition<SurfaceLoadCondition3D>;

} // Namespace Kratos


