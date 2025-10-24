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
#include "utilities/indirect_scalar.h"

namespace Kratos
{
    template <class TPrimalCondition>
    AdjointSemiAnalyticBaseCondition<TPrimalCondition>::ThisExtensions::ThisExtensions(Condition* pCondition)
        : mpCondition{pCondition}
    {
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::ThisExtensions::GetFirstDerivativesVector(
        std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
    {

        auto& r_node = mpCondition->GetGeometry()[NodeId];
        const SizeType dimension = mpCondition->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) ?  2 * dimension : dimension; // *2 for rotation
        rVector.resize(num_dofs );
        std::size_t index = 0;

        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_X, Step);
        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_Y, Step);
        if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_Z, Step);

        if (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) {
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_ROTATION_VECTOR_2_X, Step);
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_ROTATION_VECTOR_2_Y, Step);
            if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_ROTATION_VECTOR_2_Z, Step);    
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::ThisExtensions::GetSecondDerivativesVector(
        std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
    {
        auto& r_node = mpCondition->GetGeometry()[NodeId];
        const SizeType dimension = mpCondition->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) ?  2 * dimension : dimension; // *2 for rotation
        rVector.resize(num_dofs );
        std::size_t index = 0;

        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_X, Step);
        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_Y, Step);
        if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_Z, Step);

        if (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) {
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_ROTATION_VECTOR_3_X, Step);
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_ROTATION_VECTOR_3_Y, Step);
            if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_ROTATION_VECTOR_3_Z, Step);    
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::ThisExtensions::GetAuxiliaryVector(
        std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
    {
        auto& r_node = mpCondition->GetGeometry()[NodeId];
        const SizeType dimension = mpCondition->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) ?  2 * dimension : dimension; // *2 for rotation
        rVector.resize(num_dofs ); 
        std::size_t index = 0;

        rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_X, Step);
        rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_Y, Step);
        if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_Z, Step);

        //std::cout << "in auxiliary vector" << std::endl;
        if (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) {
            //std::cout << "auxiliary vector has rotation" << std::endl;
            rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_ROTATION_VECTOR_1_X, Step);
            rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_ROTATION_VECTOR_1_Y, Step);
            if (mpCondition->GetGeometry().WorkingSpaceDimension() == 3) rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_ROTATION_VECTOR_1_Z, Step);    
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::ThisExtensions::GetFirstDerivativesVariables(
        std::vector<VariableData const*>& rVariables) const
    {
        if (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) {
            rVariables.resize(2);
            rVariables[0] = &ADJOINT_VECTOR_2;
            rVariables[1] = &ADJOINT_ROTATION_VECTOR_2;
        } else {
            rVariables.resize(1);
            rVariables[0] = &ADJOINT_VECTOR_2;
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::ThisExtensions::GetSecondDerivativesVariables(
        std::vector<VariableData const*>& rVariables) const
    {
        if (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) {
            rVariables.resize(2);
            rVariables[0] = &ADJOINT_VECTOR_3;
            rVariables[1] = &ADJOINT_ROTATION_VECTOR_3;
        } else {
            rVariables.resize(1);
            rVariables[0] = &ADJOINT_VECTOR_3;
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::ThisExtensions::GetAuxiliaryVariables(
        std::vector<VariableData const*>& rVariables) const
    {
        if (mpCondition->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X)) {
            rVariables.resize(2);
            rVariables[0] = &AUX_ADJOINT_VECTOR_1;
            rVariables[1] = &AUX_ADJOINT_ROTATION_VECTOR_1;
        } else {
            rVariables.resize(1);
            rVariables[0] = &AUX_ADJOINT_VECTOR_1;
        }
    }
    namespace AdjointSemiAnalyticBaseConditionHelperUtils
    {
        template <class TData>
        void CalculateOnIntegrationPoints(
            Condition& rPrimalCondition,
            const Condition& rAdjointCondition,
            const Variable<TData>& rVariable,
            std::vector<TData>& rValues,
            const ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY

            if (rAdjointCondition.Has(rVariable)) {
                // Get result value for output
                const auto& output_value = rAdjointCondition.GetValue(rVariable);

                // Resize Output
                const auto gauss_points_number = rAdjointCondition.GetGeometry().IntegrationPointsNumber(rAdjointCondition.GetIntegrationMethod());
                if (rValues.size() != gauss_points_number) {
                    rValues.resize(gauss_points_number);
                }

                // Write scalar result value on all Gauss-Points
                for (IndexType i = 0; i < gauss_points_number; ++i) {
                    rValues[i] = output_value;
                }
            }
            else {
                rPrimalCondition.CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
            }

            KRATOS_CATCH("");
        }
    } // namespace AdjointSemiAnalyticBaseConditionHelperUtils

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dim = this->GetGeometry().WorkingSpaceDimension();
        const SizeType block_size = this->GetBlockSize();
        const SizeType num_dofs = number_of_nodes * block_size;

        rResult.resize(num_dofs, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

        if(dim == 2) {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = i * block_size;
                rResult[index    ] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
                if (this->HasRotDof())
                    rResult[index + 2] = GetGeometry()[i].GetDof(ADJOINT_ROTATION_Z,pos + 2).EquationId();
            }
        }
        else {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = i * block_size;
                rResult[index    ] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z,pos + 2).EquationId();
                if (this->HasRotDof()) {
                    rResult[index + 3] = GetGeometry()[i].GetDof(ADJOINT_ROTATION_X, pos + 3).EquationId();
                    rResult[index + 4] = GetGeometry()[i].GetDof(ADJOINT_ROTATION_Y, pos + 4).EquationId();
                    rResult[index + 5] = GetGeometry()[i].GetDof(ADJOINT_ROTATION_Z, pos + 5).EquationId();
                }
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
        const SizeType block_size = this->GetBlockSize();
        const SizeType num_dofs = number_of_nodes * block_size;

        rElementalDofList.resize(num_dofs);

        if(dimension == 2) {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = i * block_size;
                rElementalDofList[index    ] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
                rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
                if (this->HasRotDof())
                    rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Z);
            }
        }
        else {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = i * block_size;
                rElementalDofList[index    ] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
                rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
                rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z);
                if (this->HasRotDof()) {
                    rElementalDofList[index + 3] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_X);
                    rElementalDofList[index + 4] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Y);
                    rElementalDofList[index + 5] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Z);
                }
            }
        }

        KRATOS_CATCH("")
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::GetValuesVector(Vector& rValues, int Step) const
    {
        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType block_size = this->GetBlockSize();
        const SizeType num_dofs = number_of_nodes * block_size;

        rValues.resize(num_dofs, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const array_1d<double, 3 > & Displacement = this->GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
            IndexType index = i * block_size;
            for(IndexType k = 0; k < dimension; ++k) {
                rValues[index + k] = Displacement[k];
            }
            if (this->HasRotDof()) {
                const array_1d<double, 3 > & r_rotation = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_ROTATION, Step);
                for(SizeType k = 0; k < dimension; ++k) {
                    rValues[index + dimension + k] = r_rotation[k];
                }
            }
        }
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        AdjointSemiAnalyticBaseConditionHelperUtils::CalculateOnIntegrationPoints(*mpPrimalCondition, *this, rVariable, rOutput, rCurrentProcessInfo);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        AdjointSemiAnalyticBaseConditionHelperUtils::CalculateOnIntegrationPoints(*mpPrimalCondition, *this, rVariable, rOutput, rCurrentProcessInfo);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        AdjointSemiAnalyticBaseConditionHelperUtils::CalculateOnIntegrationPoints(*mpPrimalCondition, *this, rVariable, rOutput, rCurrentProcessInfo);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 4>>& rVariable,
        std::vector<array_1d<double, 4>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        AdjointSemiAnalyticBaseConditionHelperUtils::CalculateOnIntegrationPoints(*mpPrimalCondition, *this, rVariable, rOutput, rCurrentProcessInfo);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6>>& rVariable,
        std::vector<array_1d<double, 6>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        AdjointSemiAnalyticBaseConditionHelperUtils::CalculateOnIntegrationPoints(*mpPrimalCondition, *this, rVariable, rOutput, rCurrentProcessInfo);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 9>>& rVariable,
        std::vector<array_1d<double, 9>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        AdjointSemiAnalyticBaseConditionHelperUtils::CalculateOnIntegrationPoints(*mpPrimalCondition, *this, rVariable, rOutput, rCurrentProcessInfo);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        AdjointSemiAnalyticBaseConditionHelperUtils::CalculateOnIntegrationPoints(*mpPrimalCondition, *this, rVariable, rOutput, rCurrentProcessInfo);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        AdjointSemiAnalyticBaseConditionHelperUtils::CalculateOnIntegrationPoints(*mpPrimalCondition, *this, rVariable, rOutput, rCurrentProcessInfo);
    }

    template <class TPrimalCondition>
    void AdjointSemiAnalyticBaseCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        //const SizeType dimension =  this->GetGeometry().WorkingSpaceDimension();
        const SizeType block_size = this->GetBlockSize();
        const SizeType local_size = number_of_nodes * block_size;

        // this is to derive w.r.t. to variables like PRESSURE
        if(this->Has(rDesignVariable)) {
            if ((rOutput.size1() != 1) || (rOutput.size2() != local_size)) {
                rOutput.resize(1, local_size, false);
            }
            noalias(rOutput) = ZeroMatrix(1,local_size);

            const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);
            Vector RHS;
            Vector perturbed_RHS;
            this->CalculateRightHandSide(RHS, rCurrentProcessInfo);

            const auto design_variable_value = this->pGetPrimalCondition()->GetValue(rDesignVariable);

            // perturb design variable
            this->pGetPrimalCondition()->SetValue(rDesignVariable, (design_variable_value + delta));
            this->pGetPrimalCondition()->CalculateRightHandSide(perturbed_RHS, rCurrentProcessInfo);

            row(rOutput, 0) = (perturbed_RHS - RHS) / delta;

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
        const SizeType block_size = this->GetBlockSize();
        const SizeType local_size = number_of_nodes * block_size;
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

            if (this->HasRotDof()){
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION,r_node)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION,r_node)

                KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X, r_node)
                KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y, r_node)
                KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z, r_node)
            }
        }

        return return_value;

        KRATOS_CATCH( "" )
    }

    template <class TPrimalCondition>
    bool AdjointSemiAnalyticBaseCondition<TPrimalCondition>::HasRotDof() const
    {
        return (this->GetGeometry()[0].HasDofFor(ADJOINT_ROTATION_X));
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
            // the rDesignVariable value may be negative, therefore adding a std::abs
            const double variable_value = std::abs(mpPrimalCondition->GetValue(rDesignVariable));
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


