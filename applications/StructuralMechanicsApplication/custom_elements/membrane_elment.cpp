// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//


// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "custom_elements/membrane_element.hpp"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

// Constructor
MembraneElement::MembraneElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{

}

// Constructor
MembraneElement::MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element( NewId, pGeometry, pProperties )
{

}

//***********************************************************************************
//***********************************************************************************

Element::Pointer MembraneElement::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const

{
    return Kratos::make_intrusive< MembraneElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

//***********************************************************************************
//***********************************************************************************

Element::Pointer MembraneElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const

{
    return Kratos::make_intrusive< MembraneElement >(NewId, pGeom, pProperties);
}

//***********************************************************************************
//***********************************************************************************
// Destructor
MembraneElement::~MembraneElement()
{
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)

{
  KRATOS_TRY;

  unsigned int num_nodes, local_size;
  unsigned int local_index = 0;

  num_nodes = GetGeometry().size();
  local_size = num_nodes * 3;

  const unsigned int d_pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

  if (rResult.size() != local_size)
      rResult.resize(local_size, false);

  for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
  {
      rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_X, d_pos).EquationId();
      rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_Y, d_pos + 1).EquationId();
      rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_Z, d_pos + 2).EquationId();
  }

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)

{
    unsigned int num_nodes, local_size;
    num_nodes = GetGeometry().size();
    local_size = num_nodes * 3;

    if (rElementalDofList.size() != local_size)
        rElementalDofList.resize(local_size);

    unsigned int local_index = 0;

    for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
    {
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_X);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_Z);
    }
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::Initialize()

{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)

{
    const IntegrationMethod integration_method = GeometryData::GI_GAUSS_4;
    TotalStiffnessMatrix(rLeftHandSideMatrix,integration_method);
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int system_size = number_of_nodes * dimension;

    Vector internal_forces = ZeroVector(system_size);
    const IntegrationMethod integration_method = GeometryData::GI_GAUSS_4;
    InternalForces(internal_forces,integration_method);
    rRightHandSideVector -= internal_forces;
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)

{
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
}


//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetValuesVector(
    Vector& rValues,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const unsigned int index = i * 3;
        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];
    }
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetFirstDerivativesVector(
    Vector& rValues,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const unsigned int index = i * 3;
        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
    }

}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetSecondDerivativesVector(
    Vector& rValues,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const unsigned int index = i * 3;
        rValues[index] = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];
    }
}

void MembraneElement::VoigtNotation(const Matrix& rMetric, Vector& rOutputVector, const std::string StrainStressCheck)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    rOutputVector = ZeroVector(dimension);
    rOutputVector[0] = rMetric(0,0);
    rOutputVector[1] = rMetric(1,1);
    rOutputVector[2] = rMetric(0,1);
    if (StrainStressCheck=="strain") rOutputVector[2]*=2.0;
}

void MembraneElement::StressPk2(Vector& rStress,const Matrix& rShapeFunctionGradientValues)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix material_tangent_modulus = ZeroMatrix(dimension);
    MaterialTangentModulus(material_tangent_modulus,rShapeFunctionGradientValues);

    Vector strain_vector = ZeroVector(dimension);
    StrainGreenLagrange(strain_vector,rShapeFunctionGradientValues);

    rStress = prod(material_tangent_modulus,strain_vector);
}

void MembraneElement::MaterialTangentModulus(Matrix& rTangentModulus,const Matrix& rShapeFunctionGradientValues)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    rTangentModulus = ZeroMatrix(dimension);
    Matrix G = ZeroMatrix(dimension);
    ContravariantMetric(G,rShapeFunctionGradientValues,"reference");
    const double nu = GetProperties()[POISSON_RATIO];
    const double E  = GetProperties()[YOUNG_MODULUS];

    const double lam = (E*nu)/((1.0+nu)*(1-(2.0*nu)));
    const double m = E/(2.0*(1+nu));

    rTangentModulus(0,0) = (lam * G(0,0) * G(0,0)) + m * ((G(0,0) *  G(0,0)) + (G(0,0) *  G(0,0)));
    rTangentModulus(0,1) = (lam * G(0,0) * G(1,1)) + m * ((G(0,1) *  G(0,1)) + (G(0,1) *  G(0,1)));
    rTangentModulus(0,2) = (lam * G(0,0) * G(0,1)) + m * ((G(0,0) *  G(0,1)) + (G(0,1) *  G(0,0)));
    rTangentModulus(1,0) = (lam * G(1,1) * G(0,0)) + m * ((G(1,0) *  G(1,0)) + (G(1,0) *  G(1,0)));
    rTangentModulus(1,1) = (lam * G(1,1) * G(1,1)) + m * ((G(1,1) *  G(1,1)) + (G(1,1) *  G(1,1)));
    rTangentModulus(1,2) = (lam * G(1,1) * G(0,1)) + m * ((G(1,0) *  G(1,1)) + (G(1,1) *  G(1,0)));
    rTangentModulus(2,0) = (lam * G(0,1) * G(0,0)) + m * ((G(0,0) *  G(1,0)) + (G(0,0) *  G(1,0)));
    rTangentModulus(2,1) = (lam * G(0,1) * G(1,1)) + m * ((G(0,1) *  G(1,1)) + (G(0,1) *  G(1,1)));
    rTangentModulus(2,2) = (lam * G(0,1) * G(0,1)) + m * ((G(0,0) *  G(1,1)) + (G(0,1) *  G(1,0)));
}

void MembraneElement::StrainGreenLagrange(Vector& rStrain, const Matrix& rShapeFunctionGradientValues)
{
    Matrix current_covariant_metric = ZeroMatrix(2);
    Matrix reference_covariant_metric = ZeroMatrix(2);

    CovariantMetric(current_covariant_metric,rShapeFunctionGradientValues,"current");
    CovariantMetric(reference_covariant_metric,rShapeFunctionGradientValues,"reference");

    Matrix strain_matrix = 0.50 * (current_covariant_metric-reference_covariant_metric);
    VoigtNotation(strain_matrix,rStrain,"strain");
}

void MembraneElement::DerivativeStrainGreenLagrange(Vector& rStrain, const Matrix& rShapeFunctionGradientValues, const SizeType DofR)
{
    Matrix current_covariant_metric_derivative = ZeroMatrix(2);
    DerivativeCurrentCovariantMetric(current_covariant_metric_derivative,rShapeFunctionGradientValues,DofR);
    Matrix strain_matrix_derivative = 0.50 * current_covariant_metric_derivative;
    VoigtNotation(strain_matrix_derivative,rStrain,"strain");
}

void MembraneElement::Derivative2StrainGreenLagrange(Vector& rStrain,
 const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS)
{
    Matrix current_covariant_metric_derivative = ZeroMatrix(2);
    Derivative2CurrentCovariantMetric(current_covariant_metric_derivative,rShapeFunctionGradientValues,DofR,DofS);
    Matrix strain_matrix_derivative = 0.50 * current_covariant_metric_derivative;
    VoigtNotation(strain_matrix_derivative,rStrain,"strain");
}

void MembraneElement::JacobiDeterminante(double& rDetJacobi, const Matrix& rShapeFunctionGradientValues)
{
    array_1d<Vector,2> reference_base_vectors;
    CovariantBaseVectors(reference_base_vectors,rShapeFunctionGradientValues,"reference");
    Vector g3 = ZeroVector(3);
    MathUtils<double>::CrossProduct(g3, reference_base_vectors[0], reference_base_vectors[1]);
    rDetJacobi = MathUtils<double>::Norm(g3);
}

void MembraneElement::DerivativeCurrentCovariantMetric(Matrix& rMetric,
      const Matrix& rShapeFunctionGradientValues, const SizeType DofR)
{
    rMetric = ZeroMatrix(2);
    array_1d<Vector,2> covariant_base_vectors;
    CovariantBaseVectors(covariant_base_vectors,rShapeFunctionGradientValues,"current");
    array_1d<Vector,2> derivative_covariant_base_vectors;
    DeriveCurrentCovariantBaseVectors(derivative_covariant_base_vectors,rShapeFunctionGradientValues,DofR);

    for (SizeType i=0;i<2;++i){
        for (SizeType j=0;j<2;++j){
            rMetric(i,j) = inner_prod(derivative_covariant_base_vectors[i],covariant_base_vectors[j]);
            rMetric(i,j) += inner_prod(derivative_covariant_base_vectors[j],covariant_base_vectors[i]);
        }
    }
}

void MembraneElement::Derivative2CurrentCovariantMetric(Matrix& rMetric,
      const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS)
{
    rMetric = ZeroMatrix(2);
    array_1d<Vector,2> derivative_covariant_base_vectors_dur;
    DeriveCurrentCovariantBaseVectors(derivative_covariant_base_vectors_dur,rShapeFunctionGradientValues,DofR);
    array_1d<Vector,2> derivative_covariant_base_vectors_dus;
    DeriveCurrentCovariantBaseVectors(derivative_covariant_base_vectors_dus,rShapeFunctionGradientValues,DofS);

    for (SizeType i=0;i<2;++i){
        for (SizeType j=0;j<2;++j){
            rMetric(i,j) = inner_prod(derivative_covariant_base_vectors_dur[i],derivative_covariant_base_vectors_dus[j]);
            rMetric(i,j) += inner_prod(derivative_covariant_base_vectors_dus[j],derivative_covariant_base_vectors_dur[i]);
        }
    }
}

void MembraneElement::DeriveCurrentCovariantBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const Matrix& rShapeFunctionGradientValues, const SizeType DofR)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().size();
    Vector dg1 = ZeroVector(dimension);
    Vector dg2 = ZeroVector(dimension);
    Vector dudur = ZeroVector(dimension*number_of_nodes);
    dudur[DofR] = 1.0;

    for (SizeType i=0;i<number_of_nodes;++i){
        dg1[0] += (dudur[(i*dimension)+0]) * rShapeFunctionGradientValues(i, 0);
        dg1[1] += (dudur[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 0);
        dg1[2] += (dudur[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 0);

        dg2[0] += (dudur[(i*dimension)+0]) * rShapeFunctionGradientValues(i, 1);
        dg2[1] += (dudur[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 1);
        dg2[2] += (dudur[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 1);
    }
    rBaseVectors[0] = dg1;
    rBaseVectors[1] = dg2;
}

void MembraneElement::CovariantBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const Matrix& rShapeFunctionGradientValues, const std::string Configuration)
{
    // pass/call this ShapeFunctionsLocalGradients[pnt]
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().size();
    Vector g1 = ZeroVector(dimension);
    Vector g2 = ZeroVector(dimension);
    Vector current_displacement = ZeroVector(dimension*number_of_nodes);
    if (Configuration=="reference") GetValuesVector(current_displacement);
    else if (Configuration=="current");
    else KRATOS_ERROR << "configuration " << Configuration << " not known" << std::endl;

    for (SizeType i=0;i<GetGeometry().PointsNumber();++i){
        g1[0] += (GetGeometry().GetPoint( i ).X0()+current_displacement[i*dimension]) * rShapeFunctionGradientValues(i, 0);
        g1[1] += (GetGeometry().GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 0);
        g1[2] += (GetGeometry().GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 0);

        g2[0] += (GetGeometry().GetPoint( i ).X0()+current_displacement[i*dimension]) * rShapeFunctionGradientValues(i, 1);
        g2[1] += (GetGeometry().GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 1);
        g2[2] += (GetGeometry().GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 1);
    }
    rBaseVectors[0] = g1;
    rBaseVectors[1] = g2;
}

void MembraneElement::CovariantMetric(Matrix& rMetric,
      const Matrix& rShapeFunctionGradientValues, const std::string Configuration)
{
    rMetric = ZeroMatrix(2);
    array_1d<Vector,2> covariant_base_vectors;
    CovariantBaseVectors(covariant_base_vectors,rShapeFunctionGradientValues,Configuration);
    for (SizeType i=0;i<2;++i){
        for (SizeType j=0;j<2;++j){
            rMetric(i,j) = inner_prod(covariant_base_vectors[i],covariant_base_vectors[j]);
        }
    }
}

void MembraneElement::ContravariantMetric(Matrix& rMetric,
      const Matrix& rShapeFunctionGradientValues, const std::string Configuration)
{
    rMetric = ZeroMatrix(2);
    Matrix covariant_metric = ZeroMatrix(2);
    CovariantMetric(covariant_metric,rShapeFunctionGradientValues,Configuration);
    rMetric(0,0) = covariant_metric(1,1);
    rMetric(1,1) = covariant_metric(0,0);
    rMetric(0,1) = -1.0*covariant_metric(1,0);
    rMetric(1,0) = -1.0*covariant_metric(0,1);
    rMetric/=(covariant_metric(1,1)*covariant_metric(0,0)) - (covariant_metric(1,0)*covariant_metric(0,1));
}

void MembraneElement::ContraVariantBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const Matrix& rShapeFunctionGradientValues, const std::string Configuration)
{
    Matrix covariant_metric = ZeroMatrix(2);
    ContravariantMetric(covariant_metric,rShapeFunctionGradientValues,Configuration);
    array_1d<Vector,2> covariant_base_vectors;
    CovariantBaseVectors(covariant_base_vectors,rShapeFunctionGradientValues,Configuration);

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    rBaseVectors[0] = ZeroVector(dimension);
    rBaseVectors[1] = ZeroVector(dimension);

    rBaseVectors[0] = covariant_metric(0,0)*covariant_base_vectors[0] + covariant_metric(0,1)*covariant_base_vectors[1];
    rBaseVectors[1] = covariant_metric(1,0)*covariant_base_vectors[0] + covariant_metric(1,1)*covariant_base_vectors[1];
}

void MembraneElement::InternalForces(Vector& rInternalForces,const IntegrationMethod ThisMethod)
{
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType number_dofs = dimension*number_of_nodes;
    rInternalForces = ZeroVector(number_dofs);

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(ThisMethod);

    const double thickness = GetProperties()[THICKNESS];

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
        // getting information for integration
        const double integration_weight_i = r_integration_points[point_number].Weight();
        const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

        double detJ = 0.0;
        JacobiDeterminante(detJ,shape_functions_gradients_i);
        Vector stress = ZeroVector(3);
        StressPk2(stress,shape_functions_gradients_i);

        Vector derivative_strain = ZeroVector(3);
        for (SizeType dof_r=0;dof_r<number_dofs;++dof_r)
        {
            DerivativeStrainGreenLagrange(derivative_strain,shape_functions_gradients_i,dof_r);
            rInternalForces[dof_r] += inner_prod(stress,derivative_strain)*detJ*integration_weight_i*thickness;
        }
    }
}


void MembraneElement::MaterialStiffnessMatrixEntryIJ(double& rEntryIJ,
 const Matrix& rMaterialTangentModulus,const double& rDetJ, const double& rWeight,
 const SizeType& rPositionI, const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues)
 {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const double thickness = GetProperties()[THICKNESS];

    Vector strain_derivative = ZeroVector(dimension);
    DerivativeStrainGreenLagrange(strain_derivative,rShapeFunctionGradientValues,rPositionI);
    Vector stress_derivative = prod(rMaterialTangentModulus,strain_derivative);

    DerivativeStrainGreenLagrange(strain_derivative,rShapeFunctionGradientValues,rPositionJ);
    rEntryIJ += inner_prod(stress_derivative,strain_derivative)*rDetJ*rWeight*thickness;
 }

void MembraneElement::InitialStressStiffnessMatrixEntryIJ(double& rEntryIJ,
 const Vector& rStressVector,const double& rDetJ, const double& rWeight,
 const SizeType& rPositionI, const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues)
 {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const double thickness = GetProperties()[THICKNESS];

    Vector strain_derivative_2 = ZeroVector(dimension);
    Derivative2StrainGreenLagrange(strain_derivative_2,rShapeFunctionGradientValues,rPositionI,rPositionJ);
    rEntryIJ += inner_prod(rStressVector,strain_derivative_2)*rDetJ*rWeight*thickness;
 }

void MembraneElement::TotalStiffnessMatrix(Matrix& rStiffnessMatrix,const IntegrationMethod ThisMethod)
{
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType number_dofs = dimension*number_of_nodes;
    rStiffnessMatrix = ZeroMatrix(number_dofs);

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(ThisMethod);

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
        // getting information for integration
        const double integration_weight_i = r_integration_points[point_number].Weight();
        const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

        Vector stress = ZeroVector(3);
        StressPk2(stress,shape_functions_gradients_i);
        double detJ = 0.0;
        JacobiDeterminante(detJ,shape_functions_gradients_i);
        Matrix material_tangent_modulus = ZeroMatrix(dimension);
        MaterialTangentModulus(material_tangent_modulus,shape_functions_gradients_i);

        for (SizeType dof_s=0;dof_s<number_dofs;++dof_s){
            for (SizeType dof_r=0;dof_r<number_dofs;++dof_r){

                //do not calculate symmetric entries
                if(dof_s>dof_r){
                    if (point_number==(r_integration_points.size()-1)) rStiffnessMatrix(dof_s,dof_r) = rStiffnessMatrix(dof_r,dof_s);
                }
                else{
                    MaterialStiffnessMatrixEntryIJ(rStiffnessMatrix(dof_s,dof_r),
                        material_tangent_modulus,detJ,integration_weight_i,dof_s,dof_r,shape_functions_gradients_i);
                    InitialStressStiffnessMatrixEntryIJ(rStiffnessMatrix(dof_s,dof_r),
                        stress,detJ,integration_weight_i,dof_s,dof_r,shape_functions_gradients_i);
                }
            }
        }
    }
}



//***********************************************************************************
//***********************************************************************************
int MembraneElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = this->GetGeometry().size();
    // const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < number_of_nodes; i++ ) {
        const Node<3> &r_node = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW ))
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    KRATOS_ERROR_IF( strain_size != 3) << "Wrong constitutive law used. This is a membrane element! "
        << "Expected strain size is 3 (el id = " << this->Id() << ")" << std::endl;
    return 0;

    KRATOS_CATCH("");
}

void MembraneElement::save(Serializer& rSerializer) const
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
      rSerializer.save("mpConstitutiveLaw", mpConstitutiveLaw);
    }

    void MembraneElement::load(Serializer& rSerializer)
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
      rSerializer.save("mpConstitutiveLaw", mpConstitutiveLaw);
    }



//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos.
