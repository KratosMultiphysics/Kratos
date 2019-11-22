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
    KRATOS_TRY;
    if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << Id() << std::endl;
    }
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)

{
    TotalStiffnessMatrix(rLeftHandSideMatrix,GetGeometry().GetDefaultIntegrationMethod());
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
    InternalForces(internal_forces,GetGeometry().GetDefaultIntegrationMethod());
    rRightHandSideVector = ZeroVector(system_size);
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
    if (StrainStressCheck=="strain"){
        rOutputVector[2]*=2.0;
    }
}

template <class T>
void MembraneElement::InPlaneTransformationMatrix(Matrix& rTransformationMatrix, const array_1d<Vector,2>& rTransformedBaseVectors,
    const T& rLocalReferenceBaseVectors)
{
    const double e_g_11 = inner_prod(rTransformedBaseVectors[0],rLocalReferenceBaseVectors[0]);
    const double e_g_12 = inner_prod(rTransformedBaseVectors[0],rLocalReferenceBaseVectors[1]);
    const double e_g_21 = inner_prod(rTransformedBaseVectors[1],rLocalReferenceBaseVectors[0]);
    const double e_g_22 = inner_prod(rTransformedBaseVectors[1],rLocalReferenceBaseVectors[1]);
    rTransformationMatrix = ZeroMatrix(3);
    rTransformationMatrix(0,0) = e_g_11*e_g_11;
    rTransformationMatrix(0,1) = e_g_12*e_g_12;
    rTransformationMatrix(0,2) = 2.0*e_g_11*e_g_12;
    rTransformationMatrix(1,0) = e_g_21*e_g_21;
    rTransformationMatrix(1,1) = e_g_22*e_g_22;
    rTransformationMatrix(1,2) = 2.0*e_g_21*e_g_22;
    rTransformationMatrix(2,0) = e_g_11*e_g_21;
    rTransformationMatrix(2,1) = e_g_12*e_g_22;
    rTransformationMatrix(2,2) = (e_g_11*e_g_22) + (e_g_12*e_g_21);
}

void MembraneElement::TransformStrains(Vector& rStrains,
  Vector& rReferenceStrains, const Matrix& rTransformationMatrix)
{
    // use contravariant basevectors here
    // transform base vecs needs only G3 which is equal for co and contra if it is orthogonal
    // tranform strains needs contra-variant !
    rStrains = ZeroVector(3);
    rReferenceStrains[2]/=2.0; // extract E12 from voigt strain vector
    rStrains = prod(rTransformationMatrix,rReferenceStrains);
    rStrains[2]*=2.0; // include E12 and E21 for voigt strain vector
}

void MembraneElement::AddPreStressPk2(Vector& rStress, const array_1d<Vector,2>& rTransformedBaseVectors){

    Vector pre_stress = ZeroVector(3);
    if (GetProperties().Has(PRESTRESS_VECTOR)){
        pre_stress = GetProperties()(PRESTRESS_VECTOR);

        if (GetProperties().Has(LOCAL_PRESTRESS_AXIS_1) && GetProperties().Has(LOCAL_PRESTRESS_AXIS_2)){

            array_1d<array_1d<double,3>,2> local_prestress_axis;
            local_prestress_axis[0] = GetProperties()(LOCAL_PRESTRESS_AXIS_1)/MathUtils<double>::Norm(GetProperties()(LOCAL_PRESTRESS_AXIS_1));
            local_prestress_axis[1] = GetProperties()(LOCAL_PRESTRESS_AXIS_2)/MathUtils<double>::Norm(GetProperties()(LOCAL_PRESTRESS_AXIS_2));

            Matrix transformation_matrix = ZeroMatrix(3);
            InPlaneTransformationMatrix(transformation_matrix,rTransformedBaseVectors,local_prestress_axis);
            pre_stress = prod(transformation_matrix,pre_stress);
        }
    }

    rStress += pre_stress;
}

void MembraneElement::StressPk2(Vector& rStress,
    const Matrix& rReferenceContraVariantMetric,const Matrix& rReferenceCoVariantMetric,const Matrix& rCurrentCoVariantMetric,
    const array_1d<Vector,2>& rTransformedBaseVectors,const Matrix& rTransformationMatrix)
{
    Vector strain_vector = ZeroVector(3);
    rStress = ZeroVector(3);
    StrainGreenLagrange(strain_vector,rReferenceCoVariantMetric,
        rCurrentCoVariantMetric,rTransformationMatrix);

    ProcessInfo temp_process_information;
    ConstitutiveLaw::Parameters element_parameters(GetGeometry(),GetProperties(),temp_process_information);
    element_parameters.SetStrainVector(strain_vector);
    element_parameters.SetStressVector(rStress);
    element_parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    element_parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    mpConstitutiveLaw->CalculateMaterialResponse(element_parameters,ConstitutiveLaw::StressMeasure_PK2);

    AddPreStressPk2(rStress,rTransformedBaseVectors);
}

void MembraneElement::MaterialTangentModulus(Matrix& rTangentModulus,const Matrix& rReferenceContraVariantMetric,
    const Matrix& rReferenceCoVariantMetric,const Matrix& rCurrentCoVariantMetric, const Matrix& rTransformationMatrix)
{
    rTangentModulus = ZeroMatrix(3);
    Vector strain_vector = ZeroVector(3);
    StrainGreenLagrange(strain_vector,rReferenceCoVariantMetric,
        rCurrentCoVariantMetric,rTransformationMatrix);

    ProcessInfo temp_process_information;
    ConstitutiveLaw::Parameters element_parameters(GetGeometry(),GetProperties(),temp_process_information);
    element_parameters.SetStrainVector(strain_vector);
    element_parameters.SetConstitutiveMatrix(rTangentModulus);
    element_parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    element_parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    mpConstitutiveLaw->CalculateMaterialResponse(element_parameters,ConstitutiveLaw::StressMeasure_PK2);

}

void MembraneElement::StrainGreenLagrange(Vector& rStrain, const Matrix& rReferenceCoVariantMetric,const Matrix& rCurrentCoVariantMetric,
    const Matrix& rTransformationMatrix)
{
    Matrix strain_matrix = 0.50 * (rCurrentCoVariantMetric-rReferenceCoVariantMetric);
    Vector reference_strain = ZeroVector(3);
    VoigtNotation(strain_matrix,reference_strain,"strain");
    TransformStrains(rStrain,reference_strain,rTransformationMatrix);
}

void MembraneElement::DerivativeStrainGreenLagrange(Vector& rStrain, const Matrix& rShapeFunctionGradientValues, const SizeType DofR,
    const array_1d<Vector,2> rCurrentCovariantBaseVectors, const Matrix& rTransformationMatrix)
{
    Matrix current_covariant_metric_derivative = ZeroMatrix(2);
    DerivativeCurrentCovariantMetric(current_covariant_metric_derivative,rShapeFunctionGradientValues,DofR,rCurrentCovariantBaseVectors);
    Matrix strain_matrix_derivative = 0.50 * current_covariant_metric_derivative;

    Vector reference_strain = ZeroVector(3);
    VoigtNotation(strain_matrix_derivative,reference_strain,"strain");
    TransformStrains(rStrain,reference_strain,rTransformationMatrix);
}

void MembraneElement::Derivative2StrainGreenLagrange(Vector& rStrain,
 const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS,
 const Matrix& rTransformationMatrix)
{
    Matrix current_covariant_metric_derivative = ZeroMatrix(2);
    Derivative2CurrentCovariantMetric(current_covariant_metric_derivative,rShapeFunctionGradientValues,DofR,DofS);

    Matrix strain_matrix_derivative = 0.50 * current_covariant_metric_derivative;

    Vector reference_strain = ZeroVector(3);
    VoigtNotation(strain_matrix_derivative,reference_strain,"strain");
    TransformStrains(rStrain,reference_strain,rTransformationMatrix);
}

void MembraneElement::JacobiDeterminante(double& rDetJacobi, const array_1d<Vector,2>& rReferenceBaseVectors)
{
    Vector3 g3 = ZeroVector(3);
    MathUtils<double>::CrossProduct(g3, rReferenceBaseVectors[0], rReferenceBaseVectors[1]);
    rDetJacobi = MathUtils<double>::Norm(g3);
}

void MembraneElement::DerivativeCurrentCovariantMetric(Matrix& rMetric,
      const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const array_1d<Vector,2> rCurrentCovariantBaseVectors)
{
    rMetric = ZeroMatrix(2);
    array_1d<Vector,2> derivative_covariant_base_vectors;
    DeriveCurrentCovariantBaseVectors(derivative_covariant_base_vectors,rShapeFunctionGradientValues,DofR);


    for (SizeType i=0;i<2;++i){
        for (SizeType j=0;j<2;++j){
            rMetric(i,j) = inner_prod(derivative_covariant_base_vectors[i],rCurrentCovariantBaseVectors[j]);
            rMetric(i,j) += inner_prod(derivative_covariant_base_vectors[j],rCurrentCovariantBaseVectors[i]);
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
            rMetric(i,j) += inner_prod(derivative_covariant_base_vectors_dus[i],derivative_covariant_base_vectors_dur[j]);
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
    if (Configuration=="current") GetValuesVector(current_displacement);
    else if (Configuration=="reference");
    else KRATOS_ERROR << "configuration " << Configuration << " not known" << std::endl;

    for (SizeType i=0;i<number_of_nodes;++i){
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

void MembraneElement::CovariantMetric(Matrix& rMetric,const array_1d<Vector,2>& rBaseVectorCovariant)
{
    rMetric = ZeroMatrix(2);
    for (SizeType i=0;i<2;++i){
        for (SizeType j=0;j<2;++j){
            rMetric(i,j) = inner_prod(rBaseVectorCovariant[i],rBaseVectorCovariant[j]);
        }
    }
}

void MembraneElement::ContravariantMetric(Matrix& rMetric,const Matrix& rCovariantMetric)
{
    rMetric = ZeroMatrix(2);
    rMetric(0,0) = rCovariantMetric(1,1);
    rMetric(1,1) = rCovariantMetric(0,0);
    rMetric(0,1) = -1.0*rCovariantMetric(1,0);
    rMetric(1,0) = -1.0*rCovariantMetric(0,1);
    rMetric/=(rCovariantMetric(1,1)*rCovariantMetric(0,0)) - (rCovariantMetric(1,0)*rCovariantMetric(0,1));
}

void MembraneElement::ContraVariantBaseVectors(array_1d<Vector,2>& rBaseVectors,const Matrix& rContraVariantMetric,
    const array_1d<Vector,2> rCovariantBaseVectors)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    rBaseVectors[0] = ZeroVector(dimension);
    rBaseVectors[1] = ZeroVector(dimension);

    rBaseVectors[0] = rContraVariantMetric(0,0)*rCovariantBaseVectors[0] + rContraVariantMetric(0,1)*rCovariantBaseVectors[1];
    rBaseVectors[1] = rContraVariantMetric(1,0)*rCovariantBaseVectors[0] + rContraVariantMetric(1,1)*rCovariantBaseVectors[1];
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

    array_1d<Vector,2> current_covariant_base_vectors;
    array_1d<Vector,2> reference_covariant_base_vectors;
    array_1d<Vector,2> reference_contravariant_base_vectors;

    array_1d<Vector,2> transformed_base_vectors;

    Matrix covariant_metric_current = ZeroMatrix(3);
    Matrix covariant_metric_reference = ZeroMatrix(3);
    Matrix contravariant_metric_reference = ZeroMatrix(3);
    Matrix inplane_transformation_matrix_material = ZeroMatrix(3);
    double detJ = 0.0;
    Vector stress = ZeroVector(3);
    Vector derivative_strain = ZeroVector(3);

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
        // getting information for integration
        const double integration_weight_i = r_integration_points[point_number].Weight();
        const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

        CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,"current");
        CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,"reference");

        CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
        CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);
        ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);

        ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);

        TransformBaseVectors(transformed_base_vectors,reference_contravariant_base_vectors);

        InPlaneTransformationMatrix(inplane_transformation_matrix_material,transformed_base_vectors,reference_contravariant_base_vectors);


        JacobiDeterminante(detJ,reference_covariant_base_vectors);
        StressPk2(stress,contravariant_metric_reference,covariant_metric_reference,
            covariant_metric_current,transformed_base_vectors,inplane_transformation_matrix_material);

        for (SizeType dof_r=0;dof_r<number_dofs;++dof_r)
        {
            DerivativeStrainGreenLagrange(derivative_strain,shape_functions_gradients_i,
                dof_r,current_covariant_base_vectors,inplane_transformation_matrix_material);
            rInternalForces[dof_r] += inner_prod(stress,derivative_strain)*detJ*integration_weight_i*thickness;
        }
    }
}


void MembraneElement::MaterialStiffnessMatrixEntryIJ(double& rEntryIJ,
 const Matrix& rMaterialTangentModulus,const double& rDetJ, const double& rWeight,
 const SizeType& rPositionI, const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues,
 const array_1d<Vector,2>& rCurrentCovariantBaseVectors,const array_1d<Vector,2> rLocalContraVariantBaseVectorsReference,
 const array_1d<Vector,2>& rTransformedBaseVectors, const Matrix& rTransformationMatrix)
 {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const double thickness = GetProperties()[THICKNESS];

    Vector strain_derivative = ZeroVector(dimension);
    DerivativeStrainGreenLagrange(strain_derivative,rShapeFunctionGradientValues,rPositionI,
        rCurrentCovariantBaseVectors,rTransformationMatrix);

    Vector stress_derivative = prod(rMaterialTangentModulus,strain_derivative);

    DerivativeStrainGreenLagrange(strain_derivative,rShapeFunctionGradientValues,rPositionJ,
        rCurrentCovariantBaseVectors,rTransformationMatrix);

    rEntryIJ += inner_prod(stress_derivative,strain_derivative)*rDetJ*rWeight*thickness;
 }

void MembraneElement::InitialStressStiffnessMatrixEntryIJ(double& rEntryIJ,
 const Vector& rStressVector,const double& rDetJ, const double& rWeight,
 const SizeType& rPositionI, const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues,
 const array_1d<Vector,2>& rLocalContraVariantBaseVectorsReference, const array_1d<Vector,2>& rTransformedBaseVectors,
 const Matrix& rTransformationMatrix)
 {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const double thickness = GetProperties()[THICKNESS];

    Vector strain_derivative_2 = ZeroVector(dimension);
    Derivative2StrainGreenLagrange(strain_derivative_2,rShapeFunctionGradientValues,rPositionI,rPositionJ,
        rTransformationMatrix);
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

    array_1d<Vector,2> current_covariant_base_vectors;
    array_1d<Vector,2> reference_covariant_base_vectors;
    array_1d<Vector,2> reference_contravariant_base_vectors;
    array_1d<Vector,2> transformed_base_vectors;

    Matrix covariant_metric_current = ZeroMatrix(3);
    Matrix covariant_metric_reference = ZeroMatrix(3);
    Matrix contravariant_metric_reference = ZeroMatrix(3);
    Matrix inplane_transformation_matrix_material = ZeroMatrix(3);
    double detJ = 0.0;
    Vector stress = ZeroVector(3);
    Vector derivative_strain = ZeroVector(3);

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
        // getting information for integration
        const double integration_weight_i = r_integration_points[point_number].Weight();
        const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

        CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,"current");
        CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,"reference");

        CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
        CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);
        ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);

        ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);

        TransformBaseVectors(transformed_base_vectors,reference_contravariant_base_vectors);

        InPlaneTransformationMatrix(inplane_transformation_matrix_material,transformed_base_vectors,reference_contravariant_base_vectors);

        JacobiDeterminante(detJ,reference_covariant_base_vectors);
        StressPk2(stress,contravariant_metric_reference,covariant_metric_reference,covariant_metric_current,
            transformed_base_vectors,inplane_transformation_matrix_material);

        Matrix material_tangent_modulus = ZeroMatrix(dimension);
        MaterialTangentModulus(material_tangent_modulus,contravariant_metric_reference,covariant_metric_reference,covariant_metric_current,
            inplane_transformation_matrix_material);

        for (SizeType dof_s=0;dof_s<number_dofs;++dof_s){
            for (SizeType dof_r=0;dof_r<number_dofs;++dof_r){

                //do not calculate symmetric entries
                if(dof_s>dof_r){
                    if (point_number==(r_integration_points.size()-1)) rStiffnessMatrix(dof_s,dof_r) = rStiffnessMatrix(dof_r,dof_s);
                }
                else{
                    MaterialStiffnessMatrixEntryIJ(rStiffnessMatrix(dof_s,dof_r),
                        material_tangent_modulus,detJ,integration_weight_i,dof_s,dof_r,shape_functions_gradients_i,
                        current_covariant_base_vectors,reference_contravariant_base_vectors,transformed_base_vectors,
                        inplane_transformation_matrix_material);
                    InitialStressStiffnessMatrixEntryIJ(rStiffnessMatrix(dof_s,dof_r),
                        stress,detJ,integration_weight_i,dof_s,dof_r,shape_functions_gradients_i,
                        reference_contravariant_base_vectors,transformed_base_vectors,
                        inplane_transformation_matrix_material);
                }
            }
        }
    }
}

void MembraneElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    // element with two nodes can only represent results at one node
    const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
    const unsigned int& write_points_number =
        GetGeometry().IntegrationPointsNumber(integration_method);
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }

    else if (rVariable == LOCAL_AXIS_1) {
        array_1d<Vector,2> base_vectors_current_cov;
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);
        const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(integration_method);

        Vector base_1 = ZeroVector(3);
        Vector base_2 = ZeroVector(3);
        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
            const double integration_weight_i = r_integration_points[point_number].Weight();
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
            CovariantBaseVectors(base_vectors_current_cov,shape_functions_gradients_i,"reference");
            base_1 += base_vectors_current_cov[0]*integration_weight_i;
            base_2 += base_vectors_current_cov[1]*integration_weight_i;
        }


        array_1d<Vector,2> base_vectors_integrated;
        base_vectors_integrated[0] = base_1;
        base_vectors_integrated[1] = base_2;

        array_1d<Vector,2> base_vectors_integrated_transformed;
        TransformBaseVectors(base_vectors_integrated_transformed,base_vectors_integrated);

        for (SizeType i =0; i<3; ++i) {
            rOutput[0][i] = base_vectors_integrated_transformed[0][i]; // write integrated basevector to 1st GP
        }

    } else if (rVariable == LOCAL_AXIS_2) {
        array_1d<Vector,2> base_vectors_current_cov;
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);
        const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(integration_method);

        Vector base_1 = ZeroVector(3);
        Vector base_2 = ZeroVector(3);
        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
            const double integration_weight_i = r_integration_points[point_number].Weight();
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
            CovariantBaseVectors(base_vectors_current_cov,shape_functions_gradients_i,"reference");
            base_1 += base_vectors_current_cov[0]*integration_weight_i;
            base_2 += base_vectors_current_cov[1]*integration_weight_i;
        }

        array_1d<Vector,2> base_vectors_integrated;
        base_vectors_integrated[0] = base_1;
        base_vectors_integrated[1] = base_2;

        array_1d<Vector,2> base_vectors_integrated_transformed;
        TransformBaseVectors(base_vectors_integrated_transformed,base_vectors_integrated);

        for (SizeType i =0; i<3; ++i) {
            rOutput[0][i] = base_vectors_integrated_transformed[1][i]; // write integrated basevector to 1st GP
        }


    } else if (rVariable == LOCAL_AXIS_3) {
        array_1d<Vector,2> base_vectors_current_cov;
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);
        const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(integration_method);

        Vector base_1 = ZeroVector(3);
        Vector base_2 = ZeroVector(3);
        Vector base_3 = ZeroVector(3);
        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
            const double integration_weight_i = r_integration_points[point_number].Weight();
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
            CovariantBaseVectors(base_vectors_current_cov,shape_functions_gradients_i,"reference");
            base_1 += base_vectors_current_cov[0]*integration_weight_i;
            base_2 += base_vectors_current_cov[1]*integration_weight_i;
        }

        MathUtils<double>::CrossProduct(base_3, base_1, base_2);
        base_3 /= MathUtils<double>::Norm(base_3);

        for (SizeType i =0; i<3; ++i) {
            rOutput[0][i] = base_3[i]; // write integrated basevector to 1st GP
        }
    }


    KRATOS_CATCH("")
}

void MembraneElement::TransformBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const array_1d<Vector,2>& rLocalBaseVectors){

    // prepare for projection utility

    /* // create local cartesian coordinate system aligned to global material vectors (orthotropic)
    if (GetProperties().Has(LOCAL_MATERIAL_AXIS_1) && GetProperties().Has(LOCAL_MATERIAL_AXIS_2)){
        Vector local_material_axis_1 = ZeroVector(3);
        Vector local_material_axis_2 = ZeroVector(3);

        for(unsigned int i=0; i<3;i++){
            local_material_axis_1[i] = GetProperties()(LOCAL_MATERIAL_AXIS_1)[i];
            local_material_axis_2[i] = GetProperties()(LOCAL_MATERIAL_AXIS_2)[i];
        }

        rBaseVectors[0] = local_material_axis_1/MathUtils<double>::Norm(local_material_axis_1);
        rBaseVectors[1] = local_material_axis_2/MathUtils<double>::Norm(local_material_axis_2);

    } */

    // create local cartesian coordinate system aligned to global material vectors (orthotropic)
    if (GetProperties().Has(PRESTRESS_AXIS_1_GLOBAL) && GetProperties().Has(PRESTRESS_AXIS_2_GLOBAL)){
        array_1d<double,3> global_prestress_axis1, global_prestress_axis2;
                for(unsigned int i=0; i<3;i++){
                    global_prestress_axis1[i] = GetProperties()(PRESTRESS_AXIS_1_GLOBAL)[i];
                    global_prestress_axis2[i] = GetProperties()(PRESTRESS_AXIS_2_GLOBAL)[i];
                }

        Vector base_3 = ZeroVector(3);
        MathUtils<double>::CrossProduct(base_3, rLocalBaseVectors[0], rLocalBaseVectors[1]);
        base_3 /= MathUtils<double>::Norm(base_3);

        rBaseVectors[0] = ZeroVector(3);
        rBaseVectors[1] = ZeroVector(3);


        // simple projection
        MathUtils<double>::CrossProduct(rBaseVectors[1], base_3, global_prestress_axis1);
        rBaseVectors[1] /= MathUtils<double>::Norm(rBaseVectors[1]);

        MathUtils<double>::CrossProduct(rBaseVectors[0],rBaseVectors[1], base_3);
        rBaseVectors[0] /= MathUtils<double>::Norm(rBaseVectors[0]);
    }
    else {
        // create local cartesian coordinate system
        rBaseVectors[0] = ZeroVector(3);
        rBaseVectors[1] = ZeroVector(3);
        rBaseVectors[0] = rLocalBaseVectors[0] / MathUtils<double>::Norm(rLocalBaseVectors[0]);
        rBaseVectors[1] = rLocalBaseVectors[1] - (inner_prod(rLocalBaseVectors[1],rBaseVectors[0]) * rBaseVectors[0]);
        rBaseVectors[1] /= MathUtils<double>::Norm(rBaseVectors[1]);
    }
}

void MembraneElement::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    KRATOS_CATCH("")
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
