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

// Include base h
#include "eigenvector_sensor.h"

namespace Kratos
{

/// Constructor.
EigenvectorSensor::EigenvectorSensor(
    const std::string &rName,
    const Point &rLocation,
    const double Weight,
    const Vector& rSensorValueVector)
    : BaseType(rName, rLocation, Weight),
    mSensorValueVector(rSensorValueVector)

{}

const Parameters EigenvectorSensor::GetSensorParameters() const
{
    Parameters parameters = Parameters(R"(
{
    "type"      : "eigenvector_sensor",
    "name"      : "",
    "value"     : [0.0, 0.0, 0.0],
    "location"  : [0.0, 0.0, 0.0],
    "weight"    : 0.0
})");
    parameters["name"].SetString(this->GetName());
    parameters["value"].SetVector(this->GetSensorValueVector());
    parameters["location"].SetVector(this->GetLocation());
    parameters["weight"].SetDouble(this->GetWeight());
    return parameters;
}

Parameters EigenvectorSensor::GetDefaultParameters()
{
    return Parameters(R"(
{
    "type"         : "eigenvector_sensor",
    "name"         : "",
    "value"        : [0.0, 0.0, 0.0],
    "location"     : [0.0, 0.0, 0.0],
    "weight"       : 1.0,
    "variable_data": {}
})");
}

void EigenvectorSensor::SetSensorValueVector(const Vector Value)
{
    mSensorValueVector.resize(Value.size(), false);
    mSensorValueVector = Value;
    
}

Vector EigenvectorSensor::GetSensorValueVector() const
{
    return mSensorValueVector;
}

Vector EigenvectorSensor::CalculateValueVector(ModelPart &rModelPart)
{
    std::string sensor_name = this->GetName();
    int index_number = std::stoi(sensor_name.substr(12, -1));

    const SizeType number_of_nodes = rModelPart.NumberOfNodes();
    //const SizeType dimension = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
    //const SizeType num_dofs_per_node =  2 * dimension; //3D case with rotation dofs also
    //const SizeType global_size = number_of_nodes * num_dofs_per_node;
    
    Node& rnode = rModelPart.GetNode(1);
    Node::DofsContainerType& rnode_dofs = rnode.GetDofs();
    //Element::DofsVectorType adjoint_element_dofs;
    //Element& relement = rModelPart.GetElement(1);
    //relement.GetDofList(adjoint_element_dofs,rModelPart.GetProcessInfo());
    //const SizeType num_dofs_per_node = adjoint_element_dofs.size() / relement.GetGeometry().size(); // coz adjoint has double the number of dofs
    const SizeType num_dofs_per_node = rnode_dofs.size() / 2;
    const SizeType global_size = number_of_nodes * num_dofs_per_node;

    Vector eigenvector;
    SetVectorToZero(eigenvector, global_size);

    for (IndexType i_node = 0; i_node < rModelPart.NumberOfNodes(); i_node++) {
        
        auto& r_node = rModelPart.GetNode(i_node+1);
        Matrix nodal_eigenvector = r_node.GetValue(EIGENVECTOR_MATRIX);

        for (SizeType j = 0; j < num_dofs_per_node; j++) {
            eigenvector(i_node * num_dofs_per_node + j) = nodal_eigenvector((index_number-1), j+num_dofs_per_node); //accessing the adjoint solution eigensolve
        }
    }
    //std::cout << "eigenvector calc is " << eigenvector << std::endl;
    return eigenvector;
}

void EigenvectorSensor::CalculateGradient(
    const Element &rAdjointElement,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvectorSensor::CalculateGradient(
    const Condition &rAdjointCondition,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvectorSensor::CalculateFirstDerivativesGradient(
    const Element &rAdjointElement,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvectorSensor::CalculateFirstDerivativesGradient(
    const Condition &rAdjointCondition,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvectorSensor::CalculateSecondDerivativesGradient(
    const Element &rAdjointElement,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvectorSensor::CalculateSecondDerivativesGradient(
    const Condition &rAdjointCondition,
    const Matrix &rResidualGradient,
    Vector &rResponseGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void EigenvectorSensor::CalculateGlobalPartialSensitivity(
    Element &rAdjointElement,
    const Variable<double> &rVariable,
    const Matrix &rSensitivityMatrix,
    Vector &rSensitivityGradient,
    const ProcessInfo &rProcessInfo,
    ModelPart& rModelPart)
{
    KRATOS_TRY
    rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

    this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix, 
                                                            rSensitivityGradient, rProcessInfo, rModelPart);
    
    Vector c_der = rSensitivityGradient;
    Vector c = this->GetSensorValueVector();
    Vector m = this->GetValue(SENSOR_MEASURED_VALUE_VECTOR);
    //std::cout << "current eigen is " << this->GetName() << std::endl;
    //std::cout << "measured is m " << m << std::endl;

    double cc = inner_prod(c,c);
    double mm = inner_prod(m,m);
    double cm = inner_prod(c,m);
    double mc_der = inner_prod(m,c_der);
    double cc_der = inner_prod(c,c_der);

    double dMAC = 2* ((cm * cc * mc_der) - ( cm * cm * cc_der)) / (cc * cc * mm);

    rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
    rSensitivityGradient[0] = -1*this->GetWeight() * dMAC;
    //KRATOS_WATCH(rSensitivityMatrix);
    KRATOS_CATCH("");
}

void EigenvectorSensor::CalculatePartialSensitivity(
    Condition &rAdjointCondition,
    const Variable<double> &rVariable,
    const Matrix &rSensitivityMatrix,
    Vector &rSensitivityGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void EigenvectorSensor::CalculatePartialSensitivity(
    Element &rAdjointElement,
    const Variable<array_1d<double, 3>> &rVariable,
    const Matrix &rSensitivityMatrix,
    Vector &rSensitivityGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void EigenvectorSensor::CalculatePartialSensitivity(
    Condition &rAdjointCondition,
    const Variable<array_1d<double, 3>> &rVariable,
    const Matrix &rSensitivityMatrix,
    Vector &rSensitivityGradient,
    const ProcessInfo &rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

std::string EigenvectorSensor::Info() const
{
    std::stringstream msg;
    msg << "EigenvectorSensor " << this->GetName();
    return msg.str();
}

void EigenvectorSensor::PrintInfo(std::ostream &rOStream) const
{
    rOStream << Info() << std::endl;
}

void EigenvectorSensor::PrintData(std::ostream &rOStream) const
{
    PrintInfo(rOStream);
    rOStream << "    Location: " << this->GetLocation() << std::endl;
    rOStream << "    Value: " << this->GetSensorValue() << std::endl;
    rOStream << "    Weight: " << this->GetWeight() << std::endl;
    DataValueContainer::PrintData(rOStream);
}

void EigenvectorSensor::SetVectorToZero(
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
void EigenvectorSensor::DetermineEigenvectorOfElement(ModelPart::ElementType& rElement, const int eigenfrequency_id, 
                                    Vector& rEigenvectorOfElement, const ProcessInfo& CurrentProcessInfo)
{
    std::vector<std::size_t> eq_ids;
    rElement.EquationIdVector(eq_ids, CurrentProcessInfo);

    if (rEigenvectorOfElement.size() != eq_ids.size())
    {
        rEigenvectorOfElement.resize(eq_ids.size(), false);
    }

    // sort the values of the eigenvector into the rEigenvectorOfElement according to the dof ordering at the element
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

void EigenvectorSensor::CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
                                    const std::string& rVariableName,
                                    const Matrix& rSensitivityMatrix,
                                    Vector& rSensitivityGradient,
                                    const ProcessInfo& rProcessInfo,
                                    ModelPart& rModelPart)
{
    KRATOS_TRY;
    // Check if adjoint element has the design variable wrt. which sensitivity calculation is concucted and define the output size     
    rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);
    std::string& design_variable_name = rAdjointElement.GetValue(DESIGN_VARIABLE_NAME);
    //std::cout << design_variable_name << std::endl;

    std::string sensor_name = this->GetName();
    int index_number = std::stoi(sensor_name.substr(12, -1));
    double eigen_value_i = rProcessInfo[EIGENVALUE_VECTOR][index_number - 1];
    //double gradient_prefactor = 1.0 / (4.0 * Globals::Pi * std::sqrt(eigen_value_i));

    if (KratosComponents<Variable<double>>::Has(design_variable_name))
    {
        
        Element::DofsVectorType adjoint_element_dofs;
        rAdjointElement.GetDofList(adjoint_element_dofs,rProcessInfo);
        const SizeType num_dofs = adjoint_element_dofs.size();
        //std::cout << "ndof " << num_dofs << std::endl;

        // Calculate LHS and Mass matrix derivatives wrt. design variable
        Matrix LHS_design_variable_derivative = ZeroMatrix(num_dofs, num_dofs);
        rAdjointElement.Calculate(LHS_DESIGN_DERIVATIVE, LHS_design_variable_derivative, rProcessInfo);
        //std::cout << "LHS derivative is" << LHS_design_variable_derivative << std::endl;
        Matrix Massmatrix_design_variable_derivative = ZeroMatrix(num_dofs, num_dofs);
        rAdjointElement.Calculate(MASS_DESIGN_DERIVATIVE, Massmatrix_design_variable_derivative, rProcessInfo);
        
        // Predetermine the necessary eigenvectors phi_i
        std::vector<Vector> eigenvector_of_element(1,Vector(0));
        DetermineEigenvectorOfElement(rAdjointElement, index_number, eigenvector_of_element[0], rProcessInfo);

        Vector eigen_vector_i = this->CalculateValueVector(rModelPart);

        rSensitivityGradient = ZeroVector(eigen_vector_i.size());
        //std::cout << "eigen vector is " << eigenvector_of_element[0] << std::endl;
        //GetWeightedSumOfEigenvectors(rAdjointElement, weighted_sum_of_eigenvectors, rProcessInfo);
        

        //Calculate the output
        Matrix extra = ZeroMatrix(num_dofs, num_dofs);
        noalias(extra) = LHS_design_variable_derivative - eigen_value_i * Massmatrix_design_variable_derivative;
        Vector F = prod(extra, eigenvector_of_element[0]);

        for (int j = 0; j < 20; j++){
            if (j == (index_number-1)){
                double prefactor = (-1/(2*eigen_value_i)) * inner_prod(eigenvector_of_element[0], prod(LHS_design_variable_derivative, eigenvector_of_element[0]));
                
                rSensitivityGradient += prefactor * eigen_vector_i;
            }
            else{
                double eigen_value_j = rProcessInfo[EIGENVALUE_VECTOR][j];

                double phi_j_phi_j = 0.0;
                const int num_node_dofs = num_dofs / rAdjointElement.GetGeometry().size();
                for (const auto& r_node : rModelPart.Nodes()) {
                    //auto& r_node = rModelPart.GetNode(i_node+1);
                    Matrix r_nodal_eigenvector_matrix = r_node.GetValue(EIGENVECTOR_MATRIX);
                    //phi_j_phi_j += inner_prod(r_nodal_eigenvector_matrix(j,-1), r_nodal_eigenvector_matrix(j,-1));
                    for (int k = 0; k < num_node_dofs; k++) {
                        phi_j_phi_j += r_nodal_eigenvector_matrix(j, k+num_node_dofs) * r_nodal_eigenvector_matrix(j, k+num_node_dofs); //adding num_node_dofs coz accessing the adjoint solve
                    }
                }       
                double prefactor = (1/(eigen_value_i - eigen_value_j)) * phi_j_phi_j;

                int current_node = 0;
                for (auto& r_node_i : rAdjointElement.GetGeometry())
                {
                    const int num_node_dofs = num_dofs / rAdjointElement.GetGeometry().size();
                    int node_id = r_node_i.Id();
                    //const Matrix& rNodeEigenvectors = r_node_i.GetValue(EIGENVECTOR_MATRIX);
                    for (int dof_index = 0; dof_index < num_node_dofs; dof_index++)
                    {
                        rSensitivityGradient[node_id * num_node_dofs + dof_index] += prefactor * F[current_node * num_node_dofs + dof_index]; //correct to add at right place of global system
                        // const auto& current_dof = *(std::begin(r_node_dofs) + dof_index);
                        // const std::size_t dof_index_at_element = std::distance(eq_ids.begin(), std::find(eq_ids.begin(), eq_ids.end(), current_dof->EquationId()));
                        // rEigenvectorOfElement(dof_index_at_element) = rNodeEigenvectors((eigenfrequency_id-1), dof_index);
                    }
                    current_node +=1;
                }               
            }        
        } 
        //std::cout << "sensitivity gradient is " << rSensitivityGradient[0] << std::endl;
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



}; // namespace Kratos
