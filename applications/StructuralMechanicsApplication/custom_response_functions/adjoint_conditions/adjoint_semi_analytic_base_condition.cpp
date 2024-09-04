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
#include "adjoint_semi_analytic_base_condition.h"
#include "structural_mechanics_application_variables.h"
#include "custom_conditions/point_load_condition.h"
#include "custom_conditions/surface_load_condition_3d.h"
#include "custom_conditions/small_displacement_surface_load_condition_3d.h"
#include "custom_conditions/line_load_condition.h"
#include "custom_conditions/small_displacement_line_load_condition.h"
#include "custom_response_functions/response_utilities/finite_difference_utility.h"

namespace Kratos
{

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dim = this->GetGeometry().WorkingSpaceDimension();
        if (rResult.size() != dim * number_of_nodes) {
            rResult.resize(dim*number_of_nodes,false);
        }

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

        if(dim == 2) {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = i * 2;
                rResult[index    ] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
            }
        }
        else {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = i * 3;
                rResult[index    ] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z,pos + 2).EquationId();
            }
        }
        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = number_of_nodes * dimension;

        if (rElementalDofList.size() != num_dofs) {
            rElementalDofList.resize(num_dofs);
        }

        if(dimension == 2) {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = i * 2;
                rElementalDofList[index    ] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
                rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
            }
        }
        else {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = i * 3;
                rElementalDofList[index    ] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
                rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
                rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z);
            }
        }

        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::GetValuesVector(Vector& rValues, int Step) const
    {
        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = number_of_nodes * dimension;

        // const GeometryType & geom = this->GetGeometry();
        // const NodeType & iNode = geom[0];
        // const SizeType ndofs = iNode.GetDofs().size();
        // KRATOS_WATCH(ndofs);
        // const SizeType num_dofs = number_of_nodes * ndofs /2 ; // coz node has adjoint + primal dofs
        // KRATOS_WATCH(num_dofs);

        if (rValues.size() != num_dofs) {
            rValues.resize(num_dofs, false);
        }
        //rValues.clear();

        for (IndexType i = 0; i < number_of_nodes; ++i) {

            const array_1d<double, 3 > & Displacement = this->GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
            //const array_1d<double, 3 > & Rotation = this->GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_ROTATION, Step);
            IndexType index = i * dimension;
            for(IndexType k = 0; k < dimension; ++k) {
                rValues[index + k] = Displacement[k];
            }
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable, std::vector< array_1d<double, 3 > >& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if (this->Has(rVariable)) {
            // Get result value for output
            const auto& output_value = this->GetValue(rVariable);

            // Resize Output
            const SizeType gauss_points_number = this->GetGeometry()
                .IntegrationPointsNumber(this->GetIntegrationMethod());
            if (rOutput.size() != gauss_points_number) {
                rOutput.resize(gauss_points_number);
            }

            // Write result value on all Gauss-Points
            for(IndexType i = 0; i < gauss_points_number; ++i) {
                rOutput[i] = output_value;
            }

        }
        else {
            KRATOS_ERROR << "Unsupported output variable." << std::endl;
        }

        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if (this->Has(rVariable)) {
            // Get result value for output
            const auto& output_value = this->GetValue(rVariable);

            // Resize Output
            const SizeType gauss_points_number = this->GetGeometry()
                .IntegrationPointsNumber(this->GetIntegrationMethod());
            if (rOutput.size() != gauss_points_number) {
                rOutput.resize(gauss_points_number);
            }

            // Write result value on all Gauss-Points
            for(IndexType i = 0; i < gauss_points_number; ++i) {
                rOutput[i] = output_value;
            }

        }
        else {
            KRATOS_ERROR << "Unsupported output variable." << std::endl;
        }

        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType local_size = number_of_nodes * dimension;

        // this is to derive w.r.t. to variables like PRESSURE
        if(this->Has(rDesignVariable)) {
            if ((rOutput.size1() != 1) || (rOutput.size2() != local_size)) {
                rOutput.resize(1, local_size, false);
            }
            noalias(rOutput) = ZeroMatrix(1,local_size);

            const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);
            Vector RHS;
            Vector perturbed_RHS;
            Vector sensitivity;
            sensitivity.resize( 2*local_size, false);
            sensitivity.clear();

            this->CalculateRightHandSide(RHS, rCurrentProcessInfo);

            const auto design_variable_value = this->pGetPrimalCondition()->GetValue(rDesignVariable);

            // perturb design variable
            this->pGetPrimalCondition()->SetValue(rDesignVariable, (design_variable_value + delta));
            this->pGetPrimalCondition()->CalculateRightHandSide(perturbed_RHS, rCurrentProcessInfo);

            noalias(sensitivity) = (perturbed_RHS - RHS) / delta;
            
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                //const IndexType index = i * num_dofs_per_node;
                //const IndexType index2 = i * dimension;
                for (IndexType j = 0; j < dimension; ++j)
                {
                    rOutput(0,i*dimension + j) = sensitivity[i*local_size + j];
                }
            }  

            // unperturb design variable
            this->pGetPrimalCondition()->SetValue(rDesignVariable, design_variable_value);
        }
        else {
            if ((rOutput.size1() != 0) || (rOutput.size2() != local_size)) {
                rOutput.resize(0, local_size, false);
            }
            noalias(rOutput) = ZeroMatrix(0, local_size);
        }

        KRATOS_CATCH( "" )
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
        const SizeType local_size = number_of_nodes * dimension;
        Vector RHS;
        Vector perturbed_RHS = Vector(0);
        const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);

        this->CalculateRightHandSide(RHS, rCurrentProcessInfo);

        if (rDesignVariable == SHAPE_SENSITIVITY) {
            const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> coord_directions = {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y, &SHAPE_SENSITIVITY_Z};
            Vector derived_RHS;

            if ( (rOutput.size1() != dimension * number_of_nodes) || (rOutput.size2() != local_size ) ) {
                rOutput.resize(dimension * number_of_nodes, local_size, false);
            }

            IndexType index = 0;

            Vector RHS;
            pGetPrimalCondition()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
            for(auto& node_i : mpPrimalCondition->GetGeometry()) {
                for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i) {
                    // Get pseudo-load contribution from utility
                    FiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalCondition(), RHS, *coord_directions[coord_dir_i],
                                                                                node_i, delta, derived_RHS, rCurrentProcessInfo);

                    KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

                    for(IndexType i = 0; i < derived_RHS.size(); ++i) {
                        rOutput( (coord_dir_i + index*dimension), i) = derived_RHS[i];
                    }
                }
                index++;
            }
        }
        else if (this->Has(rDesignVariable)) {
            if ((rOutput.size1() != dimension) || (rOutput.size2() != local_size)) {
                rOutput.resize(dimension, local_size, false);
            }
            noalias(rOutput) = ZeroMatrix(dimension,local_size);

            const auto variable_value = this->pGetPrimalCondition()->GetValue(rDesignVariable);
            array_1d<double,3> disturbance = ZeroVector(3);
            noalias(disturbance) = variable_value;
            for(IndexType dir_i = 0; dir_i < dimension; ++dir_i) {
                disturbance[dir_i] += delta;
                this->pGetPrimalCondition()->SetValue(rDesignVariable, disturbance);

                this->pGetPrimalCondition()->CalculateRightHandSide(perturbed_RHS, rCurrentProcessInfo);
                row(rOutput, dir_i) = (perturbed_RHS - RHS) / delta;

                disturbance[dir_i] -= delta;
                perturbed_RHS = Vector(0);
            }
            // unperturb design variable
            this->pGetPrimalCondition()->SetValue(rDesignVariable, variable_value);
        }
        else {
            if ((rOutput.size1() != 0) || (rOutput.size2() != local_size)) {
                rOutput.resize(0, local_size, false);
            }
            noalias(rOutput) = ZeroMatrix(0, local_size);
        }

        KRATOS_CATCH( "" )
    }

    template <class TPrimalCondition>
    int AdjointSemiAnalyticBaseCondition<TPrimalCondition>::Check( const ProcessInfo& rCurrentProcessInfo ) const
    {
        KRATOS_TRY

        int return_value = Condition::Check(rCurrentProcessInfo);
        KRATOS_ERROR_IF_NOT(mpPrimalCondition) << "Primal conditions pointer is nullptr!" << std::endl;

        // Check dofs
        const GeometryType& r_geom = this->GetGeometry();
        for (IndexType i = 0; i < r_geom.size(); ++i) {
            const auto& r_node = r_geom[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
        }

        return return_value;

        KRATOS_CATCH( "" )
    }

    // private
    template <class TPrimalCondition>
    double AdjointSemiAnalyticBaseCondition<TPrimalCondition>::GetPerturbationSize(const Variable<double>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const
    {
        double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
        if (rCurrentProcessInfo[ADAPT_PERTURBATION_SIZE]) {
                delta *= this->GetPerturbationSizeModificationFactor(rDesignVariable);
        }

        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
        return delta;
    }

    template <class TPrimalCondition>
    double AdjointSemiAnalyticBaseCondition<TPrimalCondition>::GetPerturbationSize(const Variable<array_1d<double,3>>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const
    {
        double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
        if (rCurrentProcessInfo[ADAPT_PERTURBATION_SIZE]) {
                delta *= this->GetPerturbationSizeModificationFactor(rDesignVariable);
        }

        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
        return delta;
    }

    template <class TPrimalCondition>
    double AdjointSemiAnalyticBaseCondition<TPrimalCondition>::GetPerturbationSizeModificationFactor(const Variable<double>& rDesignVariable) const
    {
        KRATOS_TRY;

        if ( mpPrimalCondition->Has(rDesignVariable) ) {
            const double variable_value = mpPrimalCondition->GetValue(rDesignVariable);
            return variable_value;
        }
        else {
            return 1.0;
        }

        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    double AdjointSemiAnalyticBaseCondition<TPrimalCondition>::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) const
    {
        KRATOS_TRY;

        // For shape derivatives the size of the geometry (length, area, ...) is used as default perturbation size modification factor.
        // Later on this value is multiplied with a user defined factor. This product is then used as final perturbation size for computing
        // derivatives with finite differences.
        if(rDesignVariable == SHAPE_SENSITIVITY) {
            const double domain_size = mpPrimalCondition->GetGeometry().DomainSize();
            KRATOS_DEBUG_ERROR_IF(domain_size <= 0.0)
                << "Pertubation size for shape derivatives of condition" << this->Id() << "<= 0.0" << std::endl;
            return domain_size;
        }
        else {
            return 1.0;
        }

        KRATOS_CATCH("")
    }

    // TODO find out what to do with KRATOS_API
    template class AdjointSemiAnalyticBaseCondition<PointLoadCondition>;
    template class AdjointSemiAnalyticBaseCondition<SurfaceLoadCondition3D>;
    template class AdjointSemiAnalyticBaseCondition<SmallDisplacementSurfaceLoadCondition3D>;
    template class AdjointSemiAnalyticBaseCondition<LineLoadCondition<3>>;
    template class AdjointSemiAnalyticBaseCondition<SmallDisplacementLineLoadCondition<3>>;


} // Namespace Kratos


