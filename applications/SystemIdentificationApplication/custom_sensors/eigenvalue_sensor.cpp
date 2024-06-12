//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Talhah Ansari
//

// System includes
#include <cmath>

#include <iostream>
#include <iomanip>

// External includes

// Project includes

// Application includes
#include "system_identification_application_variables.h"
#include "custom_utilities/finite_difference_utility.h"

// Include base h
#include "eigenvalue_sensor.h"

namespace Kratos
{

/// Constructor.
EigenvalueSensor::EigenvalueSensor(
    const std::string &rName,
    const Point &rLocation,
    const double Weight)
    : BaseType(rName, rLocation, Weight)

{}

const Parameters EigenvalueSensor::GetSensorParameters() const
{
    Parameters parameters = Parameters(R"(
{
    "type"      : "eigenvalue_sensor",
    "name"      : "",
    "value"     : 0.0,
    "location"  : [0.0, 0.0, 0.0],
    "weight"    : 0.0
})");
    parameters["name"].SetString(this->GetName());
    parameters["value"].SetDouble(this->GetSensorValue());
    parameters["location"].SetVector(this->GetLocation());
    parameters["weight"].SetDouble(this->GetWeight());
    return parameters;
}

Parameters EigenvalueSensor::GetDefaultParameters()
{
    return Parameters(R"(
{
    "type"         : "eigenvalue_sensor",
    "name"         : "",
    "value"        : 0,
    "location"     : [0.0, 0.0, 0.0],
    "weight"       : 1.0,
    "variable_data": {}
})");
}

double EigenvalueSensor::CalculateValue(ModelPart &rModelPart)
{
    const std::string sensor_name = this->GetName();
    const int index_number = std::stoi(sensor_name.substr(11, -1));
    const double eigen_value = rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR][index_number - 1];
    const double eigenfrequency = std::sqrt(eigen_value) / (2*Globals::Pi);
    return eigenfrequency;
}

void EigenvalueSensor::CalculateGradient(
    const Element &rPrimalElement,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvalueSensor::CalculateGradient(
    const Condition &rPrimalCondition,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvalueSensor::CalculateFirstDerivativesGradient(
    const Element &rPrimalElement,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvalueSensor::CalculateFirstDerivativesGradient(
    const Condition &rPrimalCondition,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvalueSensor::CalculateSecondDerivativesGradient(
    const Element &rPrimalElement,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvalueSensor::CalculateSecondDerivativesGradient(
    const Condition &rPrimalCondition,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvalueSensor::CalculatePartialSensitivity(
    Element &rPrimalElement,
    const Variable<double> &rVariable,
    const Matrix &rSensitivityMatrix,
    Vector &rSensitivityGradient,
    const ProcessInfo &rProcessInfo)
{
    KRATOS_TRY
    rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
    //std::cout << "hi" << std::endl;
    this->CalculateElementContributionToPartialSensitivity(rPrimalElement, rVariable.Name(), rSensitivityMatrix, 
                                                            rSensitivityGradient, rProcessInfo);
    
    //KRATOS_WATCH(rSensitivityMatrix);
    KRATOS_CATCH("");
    
    // KRATOS_TRY;

    // if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
    // {
    //     rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
    // }

    // rSensitivityGradient.clear();

    // //CheckIfAllNecessaryEigenvaluesAreComputed();
    // //PerformSemiAnalyticSensitivityAnalysis();


    // // Predetermine all necessary eigenvalues and prefactors for gradient calculation
    // //const std::size_t num_of_traced_eigenfrequencies = mTracedEigenfrequencyIds.size();
    // //Vector traced_eigenvalues(num_of_traced_eigenfrequencies,0.0);
    
    // std::string sensor_name = this->GetName();
    // int index_number = std::stoi(sensor_name.substr(11, -1));
    // double eigen_value = rProcessInfo[EIGENVALUE_VECTOR][index_number - 1];
    // double gradient_prefactor = 1.0 / (4.0 * Globals::Pi * std::sqrt(eigen_value));

    // // computation of gradients
    // Matrix LHS;
    // Matrix mass_matrix;
    // rPrimalElement.CalculateLeftHandSide(LHS, rProcessInfo);
    // rPrimalElement.CalculateMassMatrix(mass_matrix, rProcessInfo);

    // const std::size_t num_dofs_element = mass_matrix.size1();
    // Vector aux_vector = Vector(num_dofs_element);
    // Matrix aux_matrix = Matrix(num_dofs_element,num_dofs_element);

    // // Predetermine the necessary eigenvectors
    // std::vector<Vector> eigenvector_of_element(1,Vector(0));
    // DetermineEigenvectorOfElement(rPrimalElement, index_number, eigenvector_of_element[0], rProcessInfo);
    
    // // Computation of derivative of LHS w.r.t. YOUNG MODULUS
    // Matrix derived_LHS = Matrix(num_dofs_element,num_dofs_element);
    // Matrix derived_mass_matrix = Matrix(num_dofs_element,num_dofs_element);
        
    // const double mdelta = rProcessInfo[PERTURBATION_SIZE];
    // //std::cout << "perturbation size is " << mdelta << std::endl;
    // CalculateLeftHandSideDerivative(rPrimalElement, LHS, mdelta, derived_LHS, rProcessInfo);
    // CalculateMassMatrixDerivative(rPrimalElement, mass_matrix, mdelta, derived_mass_matrix, rProcessInfo);
    // //std::cout << derived_mass_matrix << std::endl;

    // aux_vector.clear();
    // aux_matrix.clear();
    // //std::cout << "derived lhs is " << derived_LHS << std::endl;
    // //std::cout << "eigenvector is " << eigenvector_of_element[0] << std::endl;
    // noalias(aux_matrix) = derived_LHS - derived_mass_matrix * eigen_value;
    // noalias(aux_vector) = prod(aux_matrix , eigenvector_of_element[0]);
    // rSensitivityGradient[0] = this->GetWeight() *gradient_prefactor * inner_prod(eigenvector_of_element[0] , aux_vector);
    // std::cout << "sensitivity gradient is " << rSensitivityGradient[0] << std::endl;
    // KRATOS_CATCH("");

}

void EigenvalueSensor::CalculatePartialSensitivity(
    Condition &rPrimalCondition,
    const Variable<double> &rVariable,
    const Matrix &rSensitivityMatrix,
    Vector &rSensitivityGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void EigenvalueSensor::CalculatePartialSensitivity(
    Element &rPrimalElement,
    const Variable<array_1d<double, 3>> &rVariable,
    const Matrix &rSensitivityMatrix,
    Vector &rSensitivityGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void EigenvalueSensor::CalculatePartialSensitivity(
    Condition &rPrimalCondition,
    const Variable<array_1d<double, 3>> &rVariable,
    const Matrix &rSensitivityMatrix,
    Vector &rSensitivityGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

std::string EigenvalueSensor::Info() const
{
    std::stringstream msg;
    msg << "EigenvalueSensor " << this->GetName();
    return msg.str();
}

void EigenvalueSensor::PrintInfo(std::ostream &rOStream) const
{
    rOStream << Info() << std::endl;
}

void EigenvalueSensor::PrintData(std::ostream &rOStream) const
{
    PrintInfo(rOStream);
    rOStream << "    Location: " << this->GetLocation() << std::endl;
    rOStream << "    Value: " << this->GetSensorValue() << std::endl;
    rOStream << "    Weight: " << this->GetWeight() << std::endl;
    DataValueContainer::PrintData(rOStream);
}

void EigenvalueSensor::SetVectorToZero(
    Vector &rVector,
    const IndexType Size)
{
    if (rVector.size() != Size)
    {
        rVector.resize(Size, false);
    }

    rVector.clear();
}

// --------------------------------------------------------------------------
void EigenvalueSensor::DetermineEigenvectorOfElement(ModelPart::ElementType& rPrimalElement, const int eigenfrequency_id, 
                                    Vector& rEigenvectorOfElement, const ProcessInfo& CurrentProcessInfo)
{
    std::vector<std::size_t> eq_ids;
    rPrimalElement.EquationIdVector(eq_ids, CurrentProcessInfo);

    if (rEigenvectorOfElement.size() != eq_ids.size())
    {
        rEigenvectorOfElement.resize(eq_ids.size(), false);
    }

    // sort the values of the eigenvector into the rEigenvectorOfElement according to the dof ordering at the element
    for (auto& r_node_i : rPrimalElement.GetGeometry())
    {
        const auto& r_node_dofs = r_node_i.GetDofs();

        const Matrix& rNodeEigenvectors = r_node_i.GetValue(EIGENVECTOR_MATRIX);
        std::cout << r_node_i.GetValue(EIGENVECTOR_MATRIX) << std::endl;
        for (std::size_t dof_index = 0; dof_index < r_node_dofs.size(); dof_index++)
        {
            const auto& current_dof = *(std::begin(r_node_dofs) + dof_index);
            const std::size_t dof_index_at_element = std::distance(eq_ids.begin(), std::find(eq_ids.begin(), eq_ids.end(), current_dof->EquationId()));
            std::cout << "dof index at element "<< dof_index_at_element << " with value: " << rNodeEigenvectors((eigenfrequency_id-1), dof_index) << std::endl;
            rEigenvectorOfElement(dof_index_at_element) = rNodeEigenvectors((eigenfrequency_id-1), dof_index);
        }
    }
}

void EigenvalueSensor::CalculateLeftHandSideDesignVariableDerivative(Element &rPrimalElement,
                                             Matrix &rLHS,
                                             const double &rPerturbationSize,
                                             Matrix &rOutput,
                                             const Variable<double> &rDesignVariable,
                                             const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    // Calculates delta LHS / delta s using FiniteDifferenceUtility
    //const double FinDiffStepSize = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);

    // const SizeType number_of_nodes = rPrimalElement.GetGeometry().size();
    // const SizeType dimension = rPrimalElement.GetGeometry().WorkingSpaceDimension();
    // const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    // const SizeType mat_size = number_of_nodes * num_dofs_per_node;
    //std::cout << "mat_size " << mat_size << std::endl;
 
    //Matrix LHS;
    rPrimalElement.CalculateLeftHandSide(rLHS, rCurrentProcessInfo);
    std::cout << "LHS " << rLHS  << std::endl;
    //KRATOS_ERROR_IF_NOT(rLHS.size1() == mat_size) << "LHS has different size than dofs" << std::endl;

    rOutput.resize(rLHS.size1(), rLHS.size2(), false);
    rOutput.clear();
    FiniteDifferenceUtility::CalculateLeftHandSideDerivative(rPrimalElement, rLHS, rDesignVariable, rPerturbationSize, rOutput, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void EigenvalueSensor::CalculateMassMatrixDesignVariableDerivative(Element &rPrimalElement,
                                             Matrix &rMassmatrix,
                                             const double &rPerturbationSize,
                                             Matrix &rOutput,
                                             const Variable<double> &rDesignVariable,
                                             const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    // Calculates delta M / delta s using FiniteDifferenceUtility
    //const double FinDiffStepSize = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);

    // const SizeType number_of_nodes = rPrimalElement.GetGeometry().size();
    // const SizeType dimension = rPrimalElement.GetGeometry().WorkingSpaceDimension();
    // const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    // const SizeType mat_size = number_of_nodes * num_dofs_per_node;
        
    //Matrix Mass;
    rPrimalElement.CalculateMassMatrix(rMassmatrix, rCurrentProcessInfo);
    //KRATOS_ERROR_IF_NOT(rMassmatrix.size1() == mat_size) << "Mass Matrix has different size than dofs" << std::endl;

    rOutput.resize(rMassmatrix.size1(), rMassmatrix.size2(), false);
    rOutput.clear();
    FiniteDifferenceUtility::CalculateMassMatrixDerivative(rPrimalElement, rMassmatrix, rDesignVariable, rPerturbationSize, rOutput, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

// void EigenvalueSensor::CalculateLeftHandSideDerivative(Element& rElement,
//                                             const Matrix& rLHS,
//                                             const double& rPerturbationSize,
//                                             Matrix& rOutput,
//                                             const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY;
    
//     #pragma omp critical
//     {
//         // define working variables
//         Matrix LHS_perturbed;
//         //Matrix LHS_perturbed2;
//         Vector dummy;

//         if ( (rOutput.size1() != rLHS.size1()) || (rOutput.size2() != rLHS.size2() ) )
//             rOutput.resize(rLHS.size1(), rLHS.size2(), false);

        
//         // Save property pointer
//         Properties::Pointer p_global_properties = rElement.pGetProperties();

//         // Create new property and assign it to the element
//         Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
//         rElement.SetProperties(p_local_property);
        
//         // perturb the design variable
//         const double current_property_value = rElement.GetProperties()[YOUNG_MODULUS];
//         p_local_property->SetValue(YOUNG_MODULUS, (current_property_value + rPerturbationSize));
//         std::cout << "perturbed E is" << rElement.GetProperties()[YOUNG_MODULUS] << std::endl;
//         // compute LHS after perturbation
//         rElement.CalculateLocalSystem(LHS_perturbed, dummy, rCurrentProcessInfo);
//         //rElement.CalculateLeftHandSide(LHS_perturbed2, rCurrentProcessInfo);
//         std::cout << "perturbation size is " << rPerturbationSize << std::endl;
//         //compute derivative of RHS w.r.t. design variable with finite differences
//         noalias(rOutput) = (LHS_perturbed - rLHS) / rPerturbationSize;
//         //std::cout << "LHS is " << rLHS << std::endl;
//         //std::cout << "perturbed lhs is " << LHS_perturbed << std::endl;
//         std::cout << "difference is" << LHS_perturbed - rLHS << std::endl;
//         //std::cout << "difference2 is" << LHS_perturbed2 - rLHS << std::endl;
//         // unperturb the design variable
//         // Give element original properties back
//         rElement.SetProperties(p_global_properties);

//         // // perturb the design variable
//         // rElement.GetProperties()[YOUNG_MODULUS] += rPerturbationSize;
//         // //std::cout << "LHS matrix is " << rLHS << std::endl;
//         // // compute LHS after perturbation
//         // rElement.CalculateLocalSystem(LHS_perturbed, dummy, rCurrentProcessInfo);
//         // //std::cout << "perturbed matrix is " << LHS_perturbed << std::endl;
//         // //compute derivative of RHS w.r.t. design variable with finite differences
//         // noalias(rOutput) = (LHS_perturbed - rLHS) / rPerturbationSize;
//         // //std::cout << "difference of matrices is " << LHS_perturbed - rLHS << std::endl;
//         //     // unperturb the design variable
//         // rElement.GetProperties()[YOUNG_MODULUS] -= rPerturbationSize;
        
//     }

//     KRATOS_CATCH("");
// }

// void EigenvalueSensor::CalculateMassMatrixDerivative(Element& rElement,
//                                             const Matrix& rMassMatrix,
//                                             const double& rPerturbationSize,
//                                             Matrix& rOutput,
//                                             const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY;
    
//     #pragma omp critical
//     {
//         // define working variables
//         Matrix mass_matrix_perturbed;

//         if ( (rOutput.size1() != rMassMatrix.size1()) || (rOutput.size2() != rMassMatrix.size2() ) )
//             rOutput.resize(rMassMatrix.size1(), rMassMatrix.size2(), false);

//         // Save property pointer
//         Properties::Pointer p_global_properties = rElement.pGetProperties();

//         // Create new property and assign it to the element
//         Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
//         rElement.SetProperties(p_local_property);
        
//         // perturb the design variable
//         const double current_property_value = rElement.GetProperties()[YOUNG_MODULUS];
//         p_local_property->SetValue(YOUNG_MODULUS, (current_property_value + rPerturbationSize));

//         // compute mass matrix after perturbation
//         rElement.CalculateMassMatrix(mass_matrix_perturbed, rCurrentProcessInfo);
        
//         //compute derivative of RHS w.r.t. design variable with finite differences
//         noalias(rOutput) = (mass_matrix_perturbed - rMassMatrix) / rPerturbationSize;
       
//         // unperturb the design variable
//         // Give element original properties back
//         rElement.SetProperties(p_global_properties);
        
//     }

//     KRATOS_CATCH("");
// }

void EigenvalueSensor::CalculateElementContributionToPartialSensitivity(Element& rPrimalElement,
                                    const std::string& rVariableName,
                                    const Matrix& rSensitivityMatrix,
                                    Vector& rSensitivityGradient,
                                    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    // Check if Primal element has the design variable wrt. which sensitivity calculation is concucted and define the output size     
    rPrimalElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);
    const std::string& design_variable_name = rPrimalElement.GetValue(DESIGN_VARIABLE_NAME);
    const auto& rDesignVariable = KratosComponents<Variable<double>>::Get(design_variable_name);
    //std::cout << design_variable_name << std::endl;

    const std::string sensor_name = this->GetName();
    const int index_number = std::stoi(sensor_name.substr(11, -1));
    const double eigen_value = rProcessInfo[EIGENVALUE_VECTOR][index_number - 1];
    const double gradient_prefactor = 1.0 / (4.0 * Globals::Pi * std::sqrt(eigen_value));


    if (KratosComponents<Variable<double>>::Has(design_variable_name))
    {
        std::cout<< "hi 2" << std::endl;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        Element::DofsVectorType primal_element_dofs;
        rPrimalElement.GetDofList(primal_element_dofs,rProcessInfo);
        const SizeType num_dofs = primal_element_dofs.size();
        std::cout << "ndof " << num_dofs << std::endl;

        const double rPerturbationSize = rProcessInfo[PERTURBATION_SIZE];

        // Calculate LHS and Mass matrix derivatives wrt. design variable
        Matrix rLHS_design_variable_derivative = ZeroMatrix(num_dofs, num_dofs);
        Matrix rLHS = ZeroMatrix(num_dofs, num_dofs);
        CalculateLeftHandSideDesignVariableDerivative(rPrimalElement, rLHS, rPerturbationSize,rLHS_design_variable_derivative,rDesignVariable,rProcessInfo);
        //rPrimalElement.Calculate(LHS_DESIGN_DERIVATIVE, LHS_design_variable_derivative, rProcessInfo);
        std::cout << "LHS derivative is" << rLHS_design_variable_derivative << std::endl;
        Matrix rMassmatrix_design_variable_derivative = ZeroMatrix(num_dofs, num_dofs);
        Matrix rMassmatrix = ZeroMatrix(num_dofs, num_dofs);
        CalculateMassMatrixDesignVariableDerivative(rPrimalElement,rMassmatrix,rPerturbationSize,rMassmatrix_design_variable_derivative,
                                            rDesignVariable,rProcessInfo);
        //rPrimalElement.Calculate(MASS_DESIGN_DERIVATIVE, Massmatrix_design_variable_derivative, rProcessInfo);
        
        // Predetermine the necessary eigenvectors
        std::vector<Vector> eigenvector_of_element(1,Vector(0));
        DetermineEigenvectorOfElement(rPrimalElement, index_number, eigenvector_of_element[0], rProcessInfo);
        std::cout << "eigen vector is " << eigenvector_of_element[0] << std::endl;
        
        //Calculate the output
        Matrix extra = ZeroMatrix(num_dofs, num_dofs);
        noalias(extra) = rLHS_design_variable_derivative - eigen_value * rMassmatrix_design_variable_derivative;
        rSensitivityGradient[0] = -1*this->GetWeight() * gradient_prefactor* inner_prod(eigenvector_of_element[0], prod(extra, eigenvector_of_element[0]));
        std::cout << "sensitivity gradient is " << rSensitivityGradient[0] << std::endl;
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

    rPrimalElement.SetValue(DESIGN_VARIABLE_NAME, "");
    KRATOS_CATCH("");
}



}; // namespace Kratos
