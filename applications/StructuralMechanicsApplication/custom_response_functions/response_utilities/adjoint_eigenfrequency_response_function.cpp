// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kateryna Pindak
//

//  General notes:
//  CalculatePartialSensitivity calculates the sensitivity of the EIGENVALUES (and NOT! eigenfrequencies).
//  Depending on the input parameters, sensitivity of a chosen number of eigenvalues can be combined using user-defined weighting (similar to eigenfrequency_response_function_utility.h)

// System includes

// External includes

// Project includes
#include "adjoint_eigenfrequency_response_function.h"

namespace Kratos
{
//public Operations

// Constructor
    AdjointEigenfrequencyResponseFunction::AdjointEigenfrequencyResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings),
        mrAdjointModelPart(rModelPart)
        {
            CheckSettingsForGradientAnalysis(ResponseSettings);
            DetermineTracedEigenfrequencies(ResponseSettings);
            if(AreSeveralEigenfrequenciesTraced())
            {
                KRATOS_ERROR_IF_NOT(ResponseSettings.Has("weighting_method")) << "Several eigenfrequencies are traced but no weighting method specified in the parameters!" << std::endl;
                if(ResponseSettings["weighting_method"].GetString() == "linear_scaling")
                    CalculateLinearWeights(ResponseSettings);
                else
                    KRATOS_ERROR  << "Specified weighting method for eigenfrequencies is not implemented! Available weighting methods are: 'linear_scaling'." << std::endl;
            }
            else
                UseDefaultWeight();
        }
// Destructor
    AdjointEigenfrequencyResponseFunction::~AdjointEigenfrequencyResponseFunction(){}
// Calculate the initial value of the eigenfrequency for the primal model
    double AdjointEigenfrequencyResponseFunction::CalculateValue(ModelPart& rPrimalModelPart)
    {
        KRATOS_TRY;
        CheckIfAllNecessaryEigenvaluesAreComputed(rPrimalModelPart.GetProcessInfo());
        GetWeightedSumOfEigenvalues(rPrimalModelPart.GetProcessInfo());
        GetWeightedSumOfEigenfrequencies(rPrimalModelPart.GetProcessInfo());
        return mWeightedSumOfEigenfrequencies;
        KRATOS_CATCH("");
    }
// Calculate the adjoint load = -eigenvector^T * (delta K / delta u - eigenvalue * delta M / delta u) * eigenvector
// Here rResidualGradient is the derivative of the residuum wrt state variables delta R / delta u (calculated in another script depending on the element type).
// Adjoint variables are after calculated in CalculateSystemContributions() function in ../kratos/solving_strategies/schemes/residual_based_adjoint_static_scheme.h
    void AdjointEigenfrequencyResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                                                    const Matrix& rResidualGradient,
                                                                    Vector& rResponseGradient,
                                                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
// Adjust the size of the output vector:
        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);
        rResponseGradient.clear();
// Apply the step size for calculating gradient (step_size in adjoint input parameters):
        double gradient_delta = mDelta;
        Element::DofsVectorType adjoint_element_dofs;
        rAdjointElement.GetDofList(adjoint_element_dofs,rProcessInfo);
        mpElementInAdjointPart = mrAdjointModelPart.pGetElement(rAdjointElement.Id());
// Define DoFs of the structural problem:
        const SizeType num_nodes = mpElementInAdjointPart->GetGeometry().PointsNumber();
        const SizeType dimension = mpElementInAdjointPart->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = adjoint_element_dofs.size();
        const SizeType num_dofs_per_node = num_dofs/num_nodes;
// Check for rotation dofs:
        bool HasRotationDofs = (num_dofs_per_node/dimension == 2) ? true :false;
// Create a vector to store the result, a vector of initial values of the displacements and a helping matrix:
        Vector result = ZeroVector(num_dofs);
        Vector initial_last_step_state_variables = ZeroVector(num_dofs);
        Matrix extra = ZeroMatrix(num_dofs, num_dofs);
// Create a sum of eigenvectors:
        Vector weighted_sum_of_eigenvectors = ZeroVector(num_dofs);
        GetWeightedSumOfEigenvectors(*mpElementInAdjointPart, weighted_sum_of_eigenvectors, rProcessInfo);
// Create LHS matrices for central differencing:
        Matrix LHS_initial = ZeroMatrix(num_dofs, num_dofs);
        Matrix LHS_peturbed_plus = ZeroMatrix(num_dofs, num_dofs);
        Matrix LHS_peturbed_minus = ZeroMatrix(num_dofs, num_dofs);
// Create Mass matrices for central differencing:
        Matrix Massmatrix_initial = ZeroMatrix(num_dofs, num_dofs);
        Matrix Massmatrix_peturbed_plus = ZeroMatrix(num_dofs, num_dofs);
        Matrix Massmatrix_peturbed_minus = ZeroMatrix(num_dofs, num_dofs);
// Create Derivative matrices for central differencing:
        Matrix LHS_derivative = ZeroMatrix(num_dofs, num_dofs);
        Matrix Massmatrix_derivative = ZeroMatrix(num_dofs, num_dofs);
// Build vector of variables containing the DOF-variables of the primal problem
        std::vector<const Variable<double>*> primal_solution_variable_list;
        primal_solution_variable_list.reserve(num_dofs_per_node);
        primal_solution_variable_list.push_back(&DISPLACEMENT_X);
        primal_solution_variable_list.push_back(&DISPLACEMENT_Y);
        primal_solution_variable_list.push_back(&DISPLACEMENT_Z);
        if(HasRotationDofs)
        {
            primal_solution_variable_list.push_back(&ROTATION_X);
            primal_solution_variable_list.push_back(&ROTATION_Y);
            primal_solution_variable_list.push_back(&ROTATION_Z);
        }
// Save initial state into initial_last_step_state_variables vector
        for(IndexType i=0; i< num_nodes; ++i)
        {
            const IndexType index = i * num_dofs_per_node;
            auto& rNode = mpElementInAdjointPart->GetGeometry()[i];

            for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
            {
                initial_last_step_state_variables[index + j] = rNode.FastGetSolutionStepValue(*primal_solution_variable_list[j]);
            }
        }
// ***Initial LHS and Mass matrices are created and calculated for calculation result verification. They are not needed and may be deleted to optimize the code.
// Calculate LHS_initial and Massmatrix_initial
        LHS_initial.clear();
        Massmatrix_initial.clear();
        mpElementInAdjointPart->Calculate(LHS_PRIMAL_ELEMENT, LHS_initial, rProcessInfo);
        mpElementInAdjointPart->Calculate(MASS_MATRIX_PRIMAL_ELEMENT, Massmatrix_initial, rProcessInfo);

// Central differencing: go over every node
        for(IndexType i=0; i< num_nodes; ++i)
        {
// Create index for accessing every DoF in the given node
            const IndexType index = i * num_dofs_per_node;
            auto& rNode = mpElementInAdjointPart->GetGeometry()[i];
// Go over dofs of the node
            for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
            {
                std::size_t variable_component_index = primal_solution_variable_list[j]->GetComponentIndex();
                if(adjoint_element_dofs[index+j]->IsFree())
                {
// Positive peturbation
                    rNode.FastGetSolutionStepValue(*primal_solution_variable_list[j]) += gradient_delta;
                    rNode.Coordinates()[variable_component_index] = rNode.GetInitialPosition()[variable_component_index] +rNode.FastGetSolutionStepValue(*primal_solution_variable_list[j]);
                    
                    mpElementInAdjointPart->FinalizeNonLinearIteration(rProcessInfo);
                    mpElementInAdjointPart->InitializeNonLinearIteration(rProcessInfo);
                    
                    LHS_peturbed_plus.clear();
                    Massmatrix_peturbed_plus.clear();
                    mpElementInAdjointPart->Calculate(LHS_PRIMAL_ELEMENT, LHS_peturbed_plus, rProcessInfo);
                    mpElementInAdjointPart->Calculate(MASS_MATRIX_PRIMAL_ELEMENT, Massmatrix_peturbed_plus, rProcessInfo);

// Negative peturbation
                    rNode.FastGetSolutionStepValue(*primal_solution_variable_list[j]) -= 2 * gradient_delta;
                    rNode.Coordinates()[variable_component_index] = rNode.GetInitialPosition()[variable_component_index] +rNode.FastGetSolutionStepValue(*primal_solution_variable_list[j]);
                    
                    mpElementInAdjointPart->FinalizeNonLinearIteration(rProcessInfo);
                    mpElementInAdjointPart->InitializeNonLinearIteration(rProcessInfo);
                    
                    LHS_peturbed_minus.clear();
                    Massmatrix_peturbed_minus.clear();
                    mpElementInAdjointPart->Calculate(LHS_PRIMAL_ELEMENT, LHS_peturbed_minus, rProcessInfo);
                    mpElementInAdjointPart->Calculate(MASS_MATRIX_PRIMAL_ELEMENT, Massmatrix_peturbed_minus, rProcessInfo);
                   
                    KRATOS_ERROR_IF_NOT(LHS_derivative.size1()==LHS_peturbed_plus.size1() || LHS_derivative.size1()==LHS_peturbed_minus.size1()) << "LHS Matrices calculated for central differencing have different sizes" << std::endl;
                    KRATOS_ERROR_IF_NOT(Massmatrix_derivative.size1()==Massmatrix_peturbed_plus.size1() || Massmatrix_derivative.size1()==Massmatrix_peturbed_minus.size1())<< "Mass Matrices calculated for central differencing have different sizes" << std::endl;
                    
                    noalias(LHS_derivative) = (LHS_peturbed_plus - LHS_peturbed_minus) / (2*gradient_delta);
                    noalias(Massmatrix_derivative) = (Massmatrix_peturbed_plus - Massmatrix_peturbed_minus) / (2*gradient_delta);
                    noalias(extra) = LHS_derivative - mWeightedSumOfEigenvalues * Massmatrix_derivative;
                    result[index + j] =  inner_prod(weighted_sum_of_eigenvectors, prod(extra, weighted_sum_of_eigenvectors));

// Reset peturbation
                    rNode.FastGetSolutionStepValue(*primal_solution_variable_list[j]) = initial_last_step_state_variables[index + j];
                    rNode.Coordinates()[variable_component_index] = rNode.GetInitialPosition()[variable_component_index] +initial_last_step_state_variables[index + j];
                    mpElementInAdjointPart->FinalizeNonLinearIteration(rProcessInfo);
                    mpElementInAdjointPart->InitializeNonLinearIteration(rProcessInfo);

                }
                else
                {
                    LHS_derivative = ZeroMatrix(num_dofs, num_dofs);
                    Massmatrix_derivative = ZeroMatrix(num_dofs, num_dofs);
                    result[index + j] = 0.0;
                }

            }
        }
// Resulting adjoint load:
        noalias(rResponseGradient) = -result;
        KRATOS_CATCH(""); 
    } 
// CalculatePartialSensitivity is yet implemented only for the element sensitivity calculation with respect to double variables (no shape optimizaion)
// Calculate the Partial Sensitivity = -eigenvector^T * (delta K / delta s - eigenvalue * delta M / delta s) * eigenvector
// Sensitivity is after calculated in CalculateLocalSensitivity() function in ../kratos/solving_strategies/schemes/sensitivity_builder_scheme.h

    void AdjointEigenfrequencyResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix, rSensitivityGradient, rProcessInfo);
        KRATOS_WATCH(rSensitivityMatrix);
        KRATOS_CATCH("");
    }

    void AdjointEigenfrequencyResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR  << "No calculation of the sensitivity for Condition wrt. a double variable is yet implemented" << std::endl;
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointEigenfrequencyResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR  << "No calculation of the sensitivity for Element wrt. a 3D variable is yet implemented" << std::endl;
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointEigenfrequencyResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR  << "No calculation of the sensitivity for Condition wrt. a 3D variable is yet implemented" << std::endl;
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

//private Operations

    void AdjointEigenfrequencyResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
                                            const std::string& rVariableName,
                                            const Matrix& rSensitivityMatrix,
                                            Vector& rSensitivityGradient,
                                            const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
// Check if adjoint element has the design variable wrt. which sensitivity calculation is concucted and define the output size     
        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);
        std::string& design_variable_name = rAdjointElement.GetValue( DESIGN_VARIABLE_NAME );
        if (KratosComponents<Variable<double>>::Has(design_variable_name))
        {
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
            Element::DofsVectorType adjoint_element_dofs;
            rAdjointElement.GetDofList(adjoint_element_dofs,rProcessInfo);
            const SizeType num_dofs = adjoint_element_dofs.size();

// Calculate LHS and Mass matrix derivatives wrt. design variable
            Matrix LHS_design_variable_derivarive = ZeroMatrix(num_dofs, num_dofs);
            rAdjointElement.Calculate(LHS_DESIGN_DERIVATIVE, LHS_design_variable_derivarive, rProcessInfo);
        
            Matrix Massmatrix_design_variable_derivarive = ZeroMatrix(num_dofs, num_dofs);
            rAdjointElement.Calculate(MASS_DESIGN_DERIVATIVE, Massmatrix_design_variable_derivarive, rProcessInfo);
            
            KRATOS_ERROR_IF_NOT(LHS_design_variable_derivarive.size1()==num_dofs) << "LHS Design Variable Derivarive Matrix has wrong size" << std::endl;
            KRATOS_ERROR_IF_NOT(Massmatrix_design_variable_derivarive.size1()==num_dofs) << "Mass Design Variable Derivarive Matrix has wrong size" << std::endl;

// Create a sum of eigenvectors
            Vector weighted_sum_of_eigenvectors = ZeroVector(num_dofs);
            GetWeightedSumOfEigenvectors(rAdjointElement, weighted_sum_of_eigenvectors, rProcessInfo);

//Calculate the output
            Matrix extra = ZeroMatrix(num_dofs, num_dofs);
            noalias(extra) = LHS_design_variable_derivarive - mWeightedSumOfEigenvalues * Massmatrix_design_variable_derivarive;
            rSensitivityGradient[0] =  inner_prod(weighted_sum_of_eigenvectors, prod(extra, weighted_sum_of_eigenvectors));
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_variable_name))
        {
// Here the partial sensitivity calculation can be implemented for a 3D design variable
            KRATOS_ERROR  << "No calculation for the 3D variable is yet implemented" << std::endl;
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        }
        else
        {
            KRATOS_ERROR  << "No calculation possible for the specified variable" << std::endl;
        }

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, "");
        KRATOS_CATCH("");
    }

    void AdjointEigenfrequencyResponseFunction::CheckSettingsForGradientAnalysis(Parameters rResponseSettings)
    {
        const std::string gradient_mode = rResponseSettings["gradient_mode"].GetString();
        if (gradient_mode == "semi_analytic")
            mDelta = rResponseSettings["step_size"].GetDouble();
        else
            KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;
    }

    void AdjointEigenfrequencyResponseFunction::DetermineTracedEigenfrequencies(Parameters rResponseSettings)
    {
// In the input file, one or a combination of the eigenfrequncies can be chosen as objective function
// For example: "traced_eigenfrequencies": [1, 5, 9] could be combined using defined weighting ("weighting_factors": [0.5, 0.3, 0.2]).
        mTracedEigenfrequencyIds.resize(rResponseSettings["traced_eigenfrequencies"].size(),false);
        for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
        {
            mTracedEigenfrequencyIds[i] = rResponseSettings["traced_eigenfrequencies"][i].GetInt();
        }
    }

    bool AdjointEigenfrequencyResponseFunction::AreSeveralEigenfrequenciesTraced()
    {
        return (mTracedEigenfrequencyIds.size()>1);
    }

    void AdjointEigenfrequencyResponseFunction::CalculateLinearWeights(Parameters rResponseSettings)
    {
        KRATOS_ERROR_IF_NOT(rResponseSettings.Has("weighting_factors")) << "No weighting factors defined for given eigenfrequency response!" << std::endl;
        KRATOS_ERROR_IF_NOT(rResponseSettings["weighting_factors"].size() == mTracedEigenfrequencyIds.size()) << "The number of chosen eigenvalues does not fit to the number of weighting factors!" << std::endl;

        mWeightingFactors.resize(rResponseSettings["weighting_factors"].size(),false);

// Read weighting factors and sum them up
        double test_sum_weight_facs = 0.0;
        for(std::size_t i = 0; i < mWeightingFactors.size(); i++)
        {
            mWeightingFactors[i] = rResponseSettings["weighting_factors"][i].GetDouble();
            test_sum_weight_facs += mWeightingFactors[i];
        }

// Check the weighting factors and modify them for the case that their sum is unequal to one
        if(std::abs(test_sum_weight_facs - 1.0) > 1e-12)
        {
            for(std::size_t i = 0; i < mWeightingFactors.size(); i++)
                mWeightingFactors[i] /= test_sum_weight_facs;

            KRATOS_INFO("\n> EigenfrequencyResponseFunctionUtility") << "The sum of the chosen weighting factors is unequal one. A corresponding scaling process was exected!" << std::endl;
        }
    }
    void AdjointEigenfrequencyResponseFunction::UseDefaultWeight()
    {
        mWeightingFactors.resize(1,false);
        mWeightingFactors[0] = 1.0;
    }

    void AdjointEigenfrequencyResponseFunction::CheckIfAllNecessaryEigenvaluesAreComputed(const ProcessInfo& rProcessInfo)
    {
        const int num_of_computed_eigenvalues = (rProcessInfo[EIGENVALUE_VECTOR]).size();
        const int max_required_eigenvalue = *(std::max_element(mTracedEigenfrequencyIds.begin(), mTracedEigenfrequencyIds.end()));
        KRATOS_ERROR_IF(max_required_eigenvalue > num_of_computed_eigenvalues) << "The following Eigenfrequency shall be traced but their corresponding eigenvalue was not computed by the eigenvalue analysies: " << max_required_eigenvalue << std::endl;
    }
    
    double AdjointEigenfrequencyResponseFunction::GetEigenvalue(const int eigenfrequency_id, const ProcessInfo& rProcessInfo)
    {
        return (rProcessInfo[EIGENVALUE_VECTOR])[eigenfrequency_id-1];
    }

    void AdjointEigenfrequencyResponseFunction::GetWeightedSumOfEigenfrequencies(const ProcessInfo& rProcessInfo)
    {     
        double resp_function_value = 0.0;
        KRATOS_INFO("EigenfrequencyResponseFunctionUtility") << "CalculateValue:" << std::endl
        << "    #    Eigenfrequency [Hz]    weighting factor" <<std::endl;
        for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
        {
            const double eigenfrequency = std::sqrt(GetEigenvalue(mTracedEigenfrequencyIds[i], rProcessInfo)) / (2.0 * Globals::Pi);
            resp_function_value += mWeightingFactors[i] * eigenfrequency;
            std::stringstream msg;
            msg << std::setw(5) << mTracedEigenfrequencyIds[i]
                << std::setw(23) << eigenfrequency
                << std::setw(20) << mWeightingFactors[i]<< std::endl;
            KRATOS_INFO("") << msg.str();
        }
        mWeightedSumOfEigenfrequencies = resp_function_value;
    }

    void AdjointEigenfrequencyResponseFunction::GetWeightedSumOfEigenvalues(const ProcessInfo& rProcessInfo)
    {   
        double resp_function_value = 0.0;
        KRATOS_INFO("EigenfrequencyResponseFunctionUtility") << "CalculateValue:" << std::endl
        << "    #    Eigenfrequency [Hz]    weighting factor" <<std::endl;
        for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
        {
            const double eigenvalue = GetEigenvalue(mTracedEigenfrequencyIds[i], rProcessInfo);
            resp_function_value += mWeightingFactors[i] * eigenvalue;
            std::stringstream msg;
            msg << std::setw(5) << mTracedEigenfrequencyIds[i]
                << std::setw(23) << eigenvalue
                << std::setw(20) << mWeightingFactors[i]<< std::endl;
            KRATOS_INFO("") << msg.str();
        }
        mWeightedSumOfEigenvalues = resp_function_value;
    }


    void AdjointEigenfrequencyResponseFunction::GetWeightedSumOfEigenvectors(Element& rAdjointElement, Vector& rWeightedSumOfEigenvectors, const ProcessInfo& rProcessInfo)
    {
// Determine DoFs size
        Element::DofsVectorType adjoint_element_dofs;
        rAdjointElement.GetDofList(adjoint_element_dofs,rProcessInfo);
        const SizeType num_dofs = adjoint_element_dofs.size();
        rWeightedSumOfEigenvectors = ZeroVector(num_dofs);
// Sum the vectors
        for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
        {
            Vector current_eigenvector = ZeroVector(num_dofs);
            DetermineEigenvectorOfElement(rAdjointElement, mTracedEigenfrequencyIds[i], current_eigenvector, rProcessInfo);
            rWeightedSumOfEigenvectors += mWeightingFactors[i] * current_eigenvector;
        }
    }

    void AdjointEigenfrequencyResponseFunction::DetermineEigenvectorOfElement(ModelPart::ElementType& rElement, const int eigenfrequency_id, Vector& rEigenvectorOfElement, const ProcessInfo& CurrentProcessInfo)
    {
        std::vector<std::size_t> eq_ids;
        rElement.EquationIdVector(eq_ids, CurrentProcessInfo);

        if (rEigenvectorOfElement.size() != eq_ids.size())
        {
            rEigenvectorOfElement.resize(eq_ids.size(), false);
        }
// Sort the values of the eigenvector into the rEigenvectorOfElement according to the dof ordering at the element
        for (auto& r_node_i : rElement.GetGeometry())
        {
            const auto& r_node_dofs = r_node_i.GetDofs();
            const Matrix& rNodeEigenvectors = r_node_i.GetValue(EIGENVECTOR_MATRIX);          

            for (std::size_t dof_index = 0; dof_index < r_node_dofs.size(); dof_index++)
            {
                const auto& current_dof = *(std::begin(r_node_dofs) + dof_index);
                const std::size_t dof_index_at_element = std::distance(eq_ids.begin(), std::find(eq_ids.begin(), eq_ids.end(), current_dof->EquationId()));
                rEigenvectorOfElement(dof_index_at_element) = rNodeEigenvectors((eigenfrequency_id-1), dof_index);
            }
        }
    }

} // namespace Kratos.