// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    
//

// System includes

// External includes

// Project includes
#include "adjoint_damping_energy_dissipation_response_function.h"

namespace Kratos
{
    AdjointDampingEnergyDissipationResponseFunction::AdjointDampingEnergyDissipationResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    { 
        mResponsePartName = ResponseSettings["response_part_name"].GetString();

    }

    AdjointDampingEnergyDissipationResponseFunction::~AdjointDampingEnergyDissipationResponseFunction(){}

    void AdjointDampingEnergyDissipationResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointDampingEnergyDissipationResponseFunction::CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        
        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();
        
        ModelPart& response_part = mrModelPart.GetSubModelPart(mResponsePartName);
        const ProcessInfo &r_current_process_info = response_part.GetProcessInfo();

        // Assuming spring damper element with number of dofs = 12
        Vector velocity_vector = ZeroVector(12);
        Matrix damping_matrix;

        for (auto& element_i : response_part.Elements())
        {
            if(element_i.Id() == rAdjointElement.Id()){
                
                // get velocity vector
                auto& r_geometry = rAdjointElement.GetGeometry();
                for(IndexType n_i = 0; n_i < 2; ++n_i){
                    for(IndexType dir_i = 0; dir_i < 3; ++dir_i){
                        velocity_vector[dir_i + 2*n_i   *3] = r_geometry[n_i].FastGetSolutionStepValue(VELOCITY)[dir_i];
                        velocity_vector[dir_i +(2*n_i+1)*3] = r_geometry[n_i].FastGetSolutionStepValue(ANGULAR_VELOCITY)[dir_i];
                    }
                }
                element_i.CalculateDampingMatrix(damping_matrix, r_current_process_info);
                rResponseGradient = -2*prod(damping_matrix, velocity_vector);
                break;
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointDampingEnergyDissipationResponseFunction::CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointDampingEnergyDissipationResponseFunction::CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointDampingEnergyDissipationResponseFunction::CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointDampingEnergyDissipationResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointDampingEnergyDissipationResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointDampingEnergyDissipationResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
        
        rSensitivityGradient.clear();

        if (rVariable == NODAL_ROTATIONAL_DAMPING_RATIO || rVariable == NODAL_DAMPING_RATIO){
            // get response model part
            ModelPart& response_part = mrModelPart.GetSubModelPart(mResponsePartName);
            for (auto& element_i : response_part.Elements()){
                if(element_i.Id() == rAdjointElement.Id()){
                    // get process info of response model part
                    const ProcessInfo& r_current_process_info = response_part.GetProcessInfo();

                    // save original damping parameters
                    const auto variable_value = element_i.GetValue(rVariable);

                    // reset original damping parameters before computing the derivatives
                    element_i.SetValue(rVariable, rVariable.Zero());

                    // allocate derivative component 
                    Matrix damping_matrix_deriv;

                    // get velocity vector
                    Vector velocity_vector = ZeroVector(12);
                    auto& r_geometry = rAdjointElement.GetGeometry();
                    for(IndexType n_i = 0; n_i < 2; ++n_i){
                        for(IndexType dir_i = 0; dir_i < 3; ++dir_i){
                            velocity_vector[dir_i + 2*n_i   *3] = r_geometry[n_i].FastGetSolutionStepValue(VELOCITY)[dir_i];
                            velocity_vector[dir_i +(2*n_i+1)*3] = r_geometry[n_i].FastGetSolutionStepValue(ANGULAR_VELOCITY)[dir_i];
                        }
                    }

                    // compute partial derivative
                    // The following approach assumes a linear dependency between damping matrix and damping ratio
                    for (IndexType dir_i = 0; dir_i < 3; ++dir_i){
                        array_1d<double, 3> perturbed_nodal_damping = ZeroVector(3);
                        perturbed_nodal_damping[dir_i] = 1.0;
                        element_i.SetValue(rVariable, perturbed_nodal_damping);
                        element_i.CalculateDampingMatrix(damping_matrix_deriv,r_current_process_info);

                        rSensitivityGradient[dir_i] = -inner_prod(velocity_vector, prod(damping_matrix_deriv, velocity_vector));
                    }
                    break;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointDampingEnergyDissipationResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    double AdjointDampingEnergyDissipationResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);
        const ProcessInfo &r_current_process_info = response_part.GetProcessInfo();

        double response_value = 0.0;
        Matrix damping_matrix;
        Vector velocity_vector  = ZeroVector(12);

        for(auto& element_i : response_part.Elements())
        {
            // get velocity vector
            auto& r_geometry = element_i.GetGeometry();
                    for(IndexType n_i = 0; n_i < 2; ++n_i){
                        for(IndexType dir_i = 0; dir_i < 3; ++dir_i){
                            velocity_vector[dir_i + 2*n_i   *3] = r_geometry[n_i].FastGetSolutionStepValue(VELOCITY)[dir_i];
                            velocity_vector[dir_i +(2*n_i+1)*3] = r_geometry[n_i].FastGetSolutionStepValue(ANGULAR_VELOCITY)[dir_i];
                        }
                    }
            element_i.CalculateDampingMatrix(damping_matrix, r_current_process_info);
            // Compute the dissipation work integrand -v*C*v
            response_value -= inner_prod(velocity_vector, prod(damping_matrix, velocity_vector)); 
        }

        return response_value;

        KRATOS_CATCH("");
    }
    
} // namespace Kratos.


