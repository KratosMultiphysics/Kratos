// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Klaus B. Sautter
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "membrane_element.hpp"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "utilities/integration_utilities.h"
#include "utilities/atomic_utilities.h"

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

void MembraneElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
  KRATOS_TRY;

  SizeType num_nodes, local_size;
  SizeType local_index = 0;

  num_nodes = GetGeometry().size();
  local_size = num_nodes * 3;

  const SizeType d_pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

  if (rResult.size() != local_size)
      rResult.resize(local_size, false);

  for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
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
    const ProcessInfo& rCurrentProcessInfo) const
{
    SizeType num_nodes, local_size;
    num_nodes = GetGeometry().size();
    local_size = num_nodes * 3;

    if (rElementalDofList.size() != local_size)
        rElementalDofList.resize(local_size);

    SizeType local_index = 0;

    for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
    {
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_X);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_Z);
    }
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());

        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != integration_points.size() )
            mConstitutiveLawVector.resize( integration_points.size() );

        if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
            const GeometryType& r_geometry = GetGeometry();
            const Properties& r_properties = GetProperties();
            const auto& N_values = r_geometry.ShapeFunctionsValues(GetIntegrationMethod());
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[point_number]->InitializeMaterial( r_properties, r_geometry, row(N_values , point_number ));
            }
        } else {
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
        }

    }
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)

{
    TotalStiffnessMatrix(rLeftHandSideMatrix,GetGeometry().GetDefaultIntegrationMethod(),rCurrentProcessInfo);
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)

{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = number_of_nodes * dimension;

    Vector internal_forces = ZeroVector(system_size);
    InternalForces(internal_forces,GetGeometry().GetDefaultIntegrationMethod(),rCurrentProcessInfo);
    rRightHandSideVector.resize(system_size);
    noalias(rRightHandSideVector) = ZeroVector(system_size);
    noalias(rRightHandSideVector) -= internal_forces;
    CalculateAndAddBodyForce(rRightHandSideVector,rCurrentProcessInfo);
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)

{
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

/**
 * ELEMENTS inherited from this class must implement this methods
 * if they need to add dynamic element contributions
 * note: first derivatives means the velocities if the displacements are the dof of the analysis
 * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
 * CalculateFirstDerivativesContributions,
 * CalculateFirstDerivativesLHS, CalculateFirstDerivativesRHS methods are : OPTIONAL
 */

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix the elemental left hand side matrix
 * @param rRightHandSideVector the elemental right hand side
 * @param rCurrentProcessInfo the current process info instance
 */
void MembraneElement::CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = number_of_nodes * dimension;

    if (rLeftHandSideMatrix.size1() != system_size || (rLeftHandSideMatrix.size2() != system_size)) {
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    }
    if (rRightHandSideVector.size() != system_size) {
        rRightHandSideVector.resize(system_size, false);
    }

    CalculateFirstDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateFirstDerivativesRHS(rRightHandSideVector, rCurrentProcessInfo);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the first derivatives contributions
 * @param rLeftHandSideMatrix the elemental left hand side matrix
 * @param rCurrentProcessInfo the current process info instance
 */
void MembraneElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = number_of_nodes * dimension;

    if (rLeftHandSideMatrix.size1() != system_size || (rLeftHandSideMatrix.size2() != system_size)) {
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    }
    CalculateDampingMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the first derivatives contributions
 * @param rRightHandSideVector the elemental right hand side vector
 * @param rCurrentProcessInfo the current process info instance
 */
void MembraneElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = number_of_nodes * dimension;

    if (rRightHandSideVector.size() != system_size) {
        rRightHandSideVector.resize(system_size, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(system_size);
}

/**
 * ELEMENTS inherited from this class must implement this methods
 * if they need to add dynamic element contributions
 * note: second derivatives means the accelerations if the displacements are the dof of the analysis
 * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
 * CalculateSecondDerivativesContributions,
 * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
 */


/**
 * this is called during the assembling process in order
 * to calculate the second derivative contributions for the LHS and RHS
 * @param rLeftHandSideMatrix the elemental left hand side matrix
 * @param rRightHandSideVector the elemental right hand side
 * @param rCurrentProcessInfo the current process info instance
 */
void MembraneElement::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = number_of_nodes * dimension;

    if (rLeftHandSideMatrix.size1() != system_size || (rLeftHandSideMatrix.size2() != system_size)) {
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    }
    if (rRightHandSideVector.size() != system_size) {
        rRightHandSideVector.resize(system_size, false);
    }
    CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateSecondDerivativesRHS(rRightHandSideVector, rCurrentProcessInfo);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the second derivatives contributions
 * @param rLeftHandSideMatrix the elemental left hand side matrix
 * @param rCurrentProcessInfo the current process info instance
 */
void MembraneElement::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = number_of_nodes * dimension;

    if (rLeftHandSideMatrix.size1() != system_size || (rLeftHandSideMatrix.size2() != system_size)) {
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    }
    CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the second derivatives contributions
 * @param rRightHandSideVector the elemental right hand side vector
 * @param rCurrentProcessInfo the current process info instance
 */
void MembraneElement::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = number_of_nodes * dimension;

    if (rRightHandSideVector.size() != system_size) {
        rRightHandSideVector.resize(system_size, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(system_size);
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetValuesVector(
    Vector& rValues,
    int Step) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (SizeType i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * 3;
        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];
    }
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetFirstDerivativesVector(
    Vector& rValues,
    int Step) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType mat_size = number_of_nodes * 3;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (SizeType i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * 3;
        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
    }

}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetSecondDerivativesVector(
    Vector& rValues,
    int Step) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType mat_size = number_of_nodes * 3;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (SizeType i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * 3;
        rValues[index] = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];
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
    noalias(rStrains) = prod(rTransformationMatrix,rReferenceStrains);
    rStrains[2]*=2.0; // include E12 and E21 for voigt strain vector
}

void MembraneElement::AddPreStressPk2(Vector& rStress, const array_1d<Vector,2>& rTransformedBaseVectors){

    Vector pre_stress = ZeroVector(3);
    if (GetProperties().Has(PRESTRESS_VECTOR)){
        pre_stress = GetProperties()(PRESTRESS_VECTOR);

        if (Has(LOCAL_PRESTRESS_AXIS_1) && Has(LOCAL_PRESTRESS_AXIS_2)){

            array_1d<array_1d<double,3>,2> local_prestress_axis;
            local_prestress_axis[0] = GetValue(LOCAL_PRESTRESS_AXIS_1)/MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_1));
            local_prestress_axis[1] = GetValue(LOCAL_PRESTRESS_AXIS_2)/MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_2));

            Matrix transformation_matrix = ZeroMatrix(3);
            InPlaneTransformationMatrix(transformation_matrix,rTransformedBaseVectors,local_prestress_axis);
            pre_stress = prod(transformation_matrix,pre_stress);

        } else if (Has(LOCAL_PRESTRESS_AXIS_1)) {

            Vector base_3 = ZeroVector(3);
            MathUtils<double>::UnitCrossProduct(base_3, rTransformedBaseVectors[0], rTransformedBaseVectors[1]);

            array_1d<array_1d<double,3>,2> local_prestress_axis;
            local_prestress_axis[0] = GetValue(LOCAL_PRESTRESS_AXIS_1)/MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_1));

            MathUtils<double>::UnitCrossProduct(local_prestress_axis[1], base_3, local_prestress_axis[0]);

            Matrix transformation_matrix = ZeroMatrix(3);
            InPlaneTransformationMatrix(transformation_matrix,rTransformedBaseVectors,local_prestress_axis);
            pre_stress = prod(transformation_matrix,pre_stress);
        }
    }
    noalias(rStress) += pre_stress;
}

void MembraneElement::MaterialResponse(Vector& rStress,
    const Matrix& rReferenceContraVariantMetric,const Matrix& rReferenceCoVariantMetric,const Matrix& rCurrentCoVariantMetric,
    const array_1d<Vector,2>& rTransformedBaseVectors,const Matrix& rTransformationMatrix,const SizeType& rIntegrationPointNumber,
    Matrix& rTangentModulus,const ProcessInfo& rCurrentProcessInfo)
{
    Vector strain_vector = ZeroVector(3);
    noalias(rStress) = ZeroVector(3);
    StrainGreenLagrange(strain_vector,rReferenceCoVariantMetric,
        rCurrentCoVariantMetric,rTransformationMatrix);

    // do this to consider the pre-stress influence in the check of the membrane state in the claw
    Vector initial_stress = ZeroVector(3);
    if (Has(MEMBRANE_PRESTRESS)){
        const Matrix& r_stress_input = GetValue(MEMBRANE_PRESTRESS);
        initial_stress += column(r_stress_input,rIntegrationPointNumber);
    } else {
        AddPreStressPk2(initial_stress,rTransformedBaseVectors);
    }
    rStress += initial_stress;

    ConstitutiveLaw::Parameters element_parameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    element_parameters.SetStrainVector(strain_vector);
    element_parameters.SetStressVector(rStress);
    element_parameters.SetConstitutiveMatrix(rTangentModulus);
    element_parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    element_parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    element_parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    mConstitutiveLawVector[rIntegrationPointNumber]->CalculateMaterialResponse(element_parameters,ConstitutiveLaw::StressMeasure_PK2);

    // do this to include the pre-stress in the actual stress state
    // rStress is reset in the claw and thus does not consider the initial stress anymore
    rStress += initial_stress;
}

void MembraneElement::StrainGreenLagrange(Vector& rStrain, const Matrix& rReferenceCoVariantMetric,const Matrix& rCurrentCoVariantMetric,
    const Matrix& rTransformationMatrix)
{
    Matrix strain_matrix = 0.50 * (rCurrentCoVariantMetric-rReferenceCoVariantMetric);
    Vector reference_strain = MathUtils<double>::StrainTensorToVector(strain_matrix,3);
    TransformStrains(rStrain,reference_strain,rTransformationMatrix);
}

void MembraneElement::DerivativeStrainGreenLagrange(Vector& rStrain, const Matrix& rShapeFunctionGradientValues, const SizeType DofR,
    const array_1d<Vector,2> rCurrentCovariantBaseVectors, const Matrix& rTransformationMatrix)
{
    Matrix current_covariant_metric_derivative = ZeroMatrix(2);
    DerivativeCurrentCovariantMetric(current_covariant_metric_derivative,rShapeFunctionGradientValues,DofR,rCurrentCovariantBaseVectors);
    Matrix strain_matrix_derivative = 0.50 * current_covariant_metric_derivative;
    Vector reference_strain = MathUtils<double>::StrainTensorToVector(strain_matrix_derivative,3);
    TransformStrains(rStrain,reference_strain,rTransformationMatrix);
}

void MembraneElement::Derivative2StrainGreenLagrange(Vector& rStrain,
 const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS,
 const Matrix& rTransformationMatrix)
{
    Matrix current_covariant_metric_derivative = ZeroMatrix(2);
    Derivative2CurrentCovariantMetric(current_covariant_metric_derivative,rShapeFunctionGradientValues,DofR,DofS);

    Matrix strain_matrix_derivative = 0.50 * current_covariant_metric_derivative;

    Vector reference_strain = MathUtils<double>::StrainTensorToVector(strain_matrix_derivative,3);
    TransformStrains(rStrain,reference_strain,rTransformationMatrix);
}

void MembraneElement::JacobiDeterminante(double& rDetJacobi, const array_1d<Vector,2>& rReferenceBaseVectors) const
{
    Vector3 g3 = ZeroVector(3);
    MathUtils<double>::CrossProduct(g3, rReferenceBaseVectors[0], rReferenceBaseVectors[1]);
    rDetJacobi = MathUtils<double>::Norm(g3);
    KRATOS_ERROR_IF(rDetJacobi<std::numeric_limits<double>::epsilon()) << "det of Jacobi smaller 0 for element with id" << Id() << std::endl;
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
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dof_nr = DofR%dimension;
    const SizeType node_nr = (DofR-dof_nr)/dimension;
    for (SizeType i=0;i<2;++i){
        rBaseVectors[i] = ZeroVector(dimension);
        rBaseVectors[i][dof_nr] = rShapeFunctionGradientValues(node_nr, i);
    }
}

void MembraneElement::CovariantBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const Matrix& rShapeFunctionGradientValues, const ConfigurationType& rConfiguration) const
{
    // pass/call this ShapeFunctionsLocalGradients[pnt]
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().size();
    Vector g1 = ZeroVector(dimension);
    Vector g2 = ZeroVector(dimension);

    Vector current_displacement = ZeroVector(dimension*number_of_nodes);
    if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);


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
    rMetric /= (rCovariantMetric(1,1)*rCovariantMetric(0,0)) - (rCovariantMetric(1,0)*rCovariantMetric(0,1));
}

void MembraneElement::ContraVariantBaseVectors(array_1d<Vector,2>& rBaseVectors,const Matrix& rContraVariantMetric,
    const array_1d<Vector,2> rCovariantBaseVectors)
{
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rBaseVectors[0] = ZeroVector(dimension);
    rBaseVectors[1] = ZeroVector(dimension);

    rBaseVectors[0] = rContraVariantMetric(0,0)*rCovariantBaseVectors[0] + rContraVariantMetric(0,1)*rCovariantBaseVectors[1];
    rBaseVectors[1] = rContraVariantMetric(1,0)*rCovariantBaseVectors[0] + rContraVariantMetric(1,1)*rCovariantBaseVectors[1];
}

void MembraneElement::InternalForces(Vector& rInternalForces,const IntegrationMethod& ThisMethod,const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = dimension*number_of_nodes;
    rInternalForces = ZeroVector(number_dofs);

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

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

        CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Current);
        CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);

        CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
        CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);
        ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);

        ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);

        TransformBaseVectors(transformed_base_vectors,reference_contravariant_base_vectors);

        InPlaneTransformationMatrix(inplane_transformation_matrix_material,transformed_base_vectors,reference_contravariant_base_vectors);


        JacobiDeterminante(detJ,reference_covariant_base_vectors);
        Matrix material_tangent_modulus = ZeroMatrix(dimension);
        MaterialResponse(stress,contravariant_metric_reference,covariant_metric_reference,covariant_metric_current,
            transformed_base_vectors,inplane_transformation_matrix_material,point_number,material_tangent_modulus,
            rCurrentProcessInfo);

        for (SizeType dof_r=0;dof_r<number_dofs;++dof_r)
        {
            DerivativeStrainGreenLagrange(derivative_strain,shape_functions_gradients_i,
                dof_r,current_covariant_base_vectors,inplane_transformation_matrix_material);
            rInternalForces[dof_r] += inner_prod(stress,derivative_strain)*detJ*integration_weight_i*thickness;
        }
    }
}

void MembraneElement::ReferenceLumpingFactors(Vector& rResult) const
{
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const IntegrationMethod integration_method = r_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& r_integrations_points = r_geom.IntegrationPoints( integration_method );
    const Matrix& r_Ncontainer = r_geom.ShapeFunctionsValues(integration_method);


    array_1d<Vector,2> reference_covariant_base_vectors;
    double detJ = 0.0;
    // Iterate over the integration points
    double domain_size = 0.0;
    for ( IndexType point_number = 0; point_number < r_integrations_points.size(); ++point_number ) {
        const Vector& rN = row(r_Ncontainer,point_number);
        const Matrix& shape_functions_gradients_i = r_geom.ShapeFunctionsLocalGradients(integration_method)[point_number];
        CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);

        JacobiDeterminante(detJ,reference_covariant_base_vectors);
        const double integration_weight = r_integrations_points[point_number].Weight() * detJ;

        // Computing domain size
        domain_size += integration_weight;

        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            rResult[i] += rN[i] * integration_weight;
        }
    }

    // Divide by the domain size
    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        rResult[i] /= domain_size;
    }
}


void MembraneElement::MaterialStiffnessMatrixEntryIJ(double& rEntryIJ,
    const Matrix& rMaterialTangentModulus,const SizeType& rPositionI,
    const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues,
    const array_1d<Vector,2>& rCurrentCovariantBaseVectors, const Matrix& rTransformationMatrix)
 {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    Vector strain_derivative = ZeroVector(dimension);
    DerivativeStrainGreenLagrange(strain_derivative,rShapeFunctionGradientValues,rPositionI,
        rCurrentCovariantBaseVectors,rTransformationMatrix);

    Vector stress_derivative = prod(rMaterialTangentModulus,strain_derivative);

    DerivativeStrainGreenLagrange(strain_derivative,rShapeFunctionGradientValues,rPositionJ,
        rCurrentCovariantBaseVectors,rTransformationMatrix);

    rEntryIJ += inner_prod(stress_derivative,strain_derivative);
 }

void MembraneElement::InitialStressStiffnessMatrixEntryIJ(double& rEntryIJ,
 const Vector& rStressVector, const SizeType& rPositionI,
 const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues,
 const Matrix& rTransformationMatrix)
 {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    Vector strain_derivative_2 = ZeroVector(dimension);
    Derivative2StrainGreenLagrange(strain_derivative_2,rShapeFunctionGradientValues,rPositionI,rPositionJ,
        rTransformationMatrix);
    rEntryIJ += inner_prod(rStressVector,strain_derivative_2);
 }


void MembraneElement::TotalStiffnessMatrix(Matrix& rStiffnessMatrix,const IntegrationMethod& ThisMethod,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = dimension*number_of_nodes;
    rStiffnessMatrix = ZeroMatrix(number_dofs);

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

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
    double temp_stiffness_entry;
    Vector stress = ZeroVector(3);

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
        // getting information for integration
        const double integration_weight_i = r_integration_points[point_number].Weight();
        const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

        CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Current);
        CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);

        CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
        CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);
        ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);

        ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);

        TransformBaseVectors(transformed_base_vectors,reference_contravariant_base_vectors);

        InPlaneTransformationMatrix(inplane_transformation_matrix_material,transformed_base_vectors,reference_contravariant_base_vectors);

        JacobiDeterminante(detJ,reference_covariant_base_vectors);

        Matrix material_tangent_modulus = ZeroMatrix(dimension);
        MaterialResponse(stress,contravariant_metric_reference,covariant_metric_reference,covariant_metric_current,
            transformed_base_vectors,inplane_transformation_matrix_material,point_number,material_tangent_modulus,
            rCurrentProcessInfo);


        for (SizeType dof_s=0;dof_s<number_dofs;++dof_s){
            for (SizeType dof_r=0;dof_r<number_dofs;++dof_r){

                //do not calculate symmetric entries
                if(dof_s>dof_r){
                    if (point_number==(r_integration_points.size()-1)) rStiffnessMatrix(dof_s,dof_r) = rStiffnessMatrix(dof_r,dof_s);
                }
                else{
                    temp_stiffness_entry = 0.0;
                    MaterialStiffnessMatrixEntryIJ(temp_stiffness_entry,
                        material_tangent_modulus,dof_s,dof_r,shape_functions_gradients_i,
                        current_covariant_base_vectors,inplane_transformation_matrix_material);
                    InitialStressStiffnessMatrixEntryIJ(temp_stiffness_entry,
                        stress,dof_s,dof_r,shape_functions_gradients_i,
                        inplane_transformation_matrix_material);
                    rStiffnessMatrix(dof_s,dof_r) += temp_stiffness_entry*detJ*integration_weight_i*thickness;
                }
            }
        }
    }
}

void MembraneElement::TransformBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const array_1d<Vector,2>& rLocalBaseVectors){

    // create local cartesian coordinate system aligned to global material vectors (orthotropic)
    if (Has(LOCAL_MATERIAL_AXIS_1) && Has(LOCAL_MATERIAL_AXIS_2)){
        rBaseVectors[0] = GetValue(LOCAL_MATERIAL_AXIS_1)/MathUtils<double>::Norm(GetValue(LOCAL_MATERIAL_AXIS_1));
        rBaseVectors[1] = GetValue(LOCAL_MATERIAL_AXIS_2)/MathUtils<double>::Norm(GetValue(LOCAL_MATERIAL_AXIS_2));
    } else if (Has(LOCAL_MATERIAL_AXIS_1)) {
        Vector base_3 = ZeroVector(3);
        MathUtils<double>::UnitCrossProduct(base_3, rLocalBaseVectors[0], rLocalBaseVectors[1]);

        rBaseVectors[0] = GetValue(LOCAL_MATERIAL_AXIS_1)/MathUtils<double>::Norm(GetValue(LOCAL_MATERIAL_AXIS_1));
        MathUtils<double>::UnitCrossProduct(rBaseVectors[1], base_3, rBaseVectors[0]);

    } else {
        // create local cartesian coordinate system
        rBaseVectors[0] = ZeroVector(3);
        rBaseVectors[1] = ZeroVector(3);
        rBaseVectors[0] = rLocalBaseVectors[0] / MathUtils<double>::Norm(rLocalBaseVectors[0]);
        rBaseVectors[1] = rLocalBaseVectors[1] - (inner_prod(rLocalBaseVectors[1],rBaseVectors[0]) * rBaseVectors[0]);
        rBaseVectors[1] /= MathUtils<double>::Norm(rBaseVectors[1]);
    }
}

void MembraneElement::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
                        std::vector< Vector >& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    // element with two nodes can only represent results at one node
    const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
    const SizeType& write_points_number =
        GetGeometry().IntegrationPointsNumber(integration_method);
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(integration_method);


    if (rVariable==PK2_STRESS_VECTOR || rVariable==PRINCIPAL_PK2_STRESS_VECTOR || rVariable == GREEN_LAGRANGE_STRAIN_VECTOR){
        Vector stress = ZeroVector(3);
        array_1d<Vector,2> current_covariant_base_vectors;
        array_1d<Vector,2> reference_covariant_base_vectors;
        array_1d<Vector,2> reference_contravariant_base_vectors;

        array_1d<Vector,2> transformed_base_vectors;

        Matrix covariant_metric_current = ZeroMatrix(3);
        Matrix covariant_metric_reference = ZeroMatrix(3);
        Matrix contravariant_metric_reference = ZeroMatrix(3);
        Matrix inplane_transformation_matrix_material = ZeroMatrix(3);

        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){

            if (rOutput[point_number].size() != 3) {
                rOutput[point_number].resize(3);
            }

            // getting information for integration
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

            CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Current);
            CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);

            CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
            CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);
            ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);

            ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);

            TransformBaseVectors(transformed_base_vectors,reference_contravariant_base_vectors);

            InPlaneTransformationMatrix(inplane_transformation_matrix_material,transformed_base_vectors,reference_contravariant_base_vectors);


            if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR){
                    Vector strain_vector = ZeroVector(3);
                    StrainGreenLagrange(strain_vector,covariant_metric_reference,
                    covariant_metric_current,inplane_transformation_matrix_material);
                    strain_vector[2] /= 2.0;
                    noalias(rOutput[point_number]) = strain_vector;
            }
            else {
                Matrix material_tangent_modulus = ZeroMatrix(dimension);
                MaterialResponse(stress,contravariant_metric_reference,covariant_metric_reference,covariant_metric_current,
                    transformed_base_vectors,inplane_transformation_matrix_material,point_number,material_tangent_modulus,
                    rCurrentProcessInfo);

                if (rVariable==PRINCIPAL_PK2_STRESS_VECTOR){
                    Vector principal_stresses = ZeroVector(2);

                    if (rOutput[point_number].size() != 2) {
                        rOutput[point_number].resize(2);
                    }


                    PrincipalVector(principal_stresses,stress);
                    noalias(rOutput[point_number]) = principal_stresses;
                }  else {
                    noalias(rOutput[point_number]) = stress;
                }
            }

        }

    }  else if (rVariable==CAUCHY_STRESS_VECTOR || rVariable==PRINCIPAL_CAUCHY_STRESS_VECTOR){

        Vector stress = ZeroVector(3);
        array_1d<Vector,2> current_covariant_base_vectors;
        array_1d<Vector,2> current_contravariant_base_vectors;
        array_1d<Vector,2> reference_covariant_base_vectors;
        array_1d<Vector,2> reference_contravariant_base_vectors;

        array_1d<Vector,2> transformed_base_vectors;

        Matrix covariant_metric_current = ZeroMatrix(3);
        Matrix covariant_metric_reference = ZeroMatrix(3);
        Matrix contravariant_metric_reference = ZeroMatrix(3);
        Matrix contravariant_metric_current = ZeroMatrix(3);
        Matrix inplane_transformation_matrix_material = ZeroMatrix(3);

        Matrix deformation_gradient = ZeroMatrix(3);
        double det_deformation_gradient = 0.0;

        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){


            if (rOutput[point_number].size() != 3) {
                rOutput[point_number].resize(3);
            }

            // getting information for integration
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

            CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Current);
            CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);

            CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
            CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);
            ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);
            ContravariantMetric(contravariant_metric_current,covariant_metric_current);

            ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);
            ContraVariantBaseVectors(current_contravariant_base_vectors,contravariant_metric_current,current_covariant_base_vectors);

            TransformBaseVectors(transformed_base_vectors,reference_contravariant_base_vectors);


            InPlaneTransformationMatrix(inplane_transformation_matrix_material,
                transformed_base_vectors,reference_contravariant_base_vectors);

            Matrix material_tangent_modulus = ZeroMatrix(dimension);
            MaterialResponse(stress,contravariant_metric_reference,covariant_metric_reference,covariant_metric_current,
                transformed_base_vectors,inplane_transformation_matrix_material,point_number,material_tangent_modulus,
                rCurrentProcessInfo);


            DeformationGradient(deformation_gradient,det_deformation_gradient,current_covariant_base_vectors,reference_contravariant_base_vectors);


            Matrix stress_matrix_local_cs = MathUtils<double>::StressVectorToTensor(stress);

            // transform stresses to original bases
            Matrix stress_matrix = ZeroMatrix(3);
            for (SizeType i=0;i<2;++i){
                for (SizeType j=0;j<2;++j){
                    stress_matrix += outer_prod(transformed_base_vectors[i],transformed_base_vectors[j]) * stress_matrix_local_cs(i,j);
                }
            }


            // calculate cauchy (this needs to be done in the original base)
            Matrix temp_stress_matrix = prod(deformation_gradient,stress_matrix);
            Matrix temp_stress_matrix_2 = prod(temp_stress_matrix,trans(deformation_gradient));
            Matrix cauchy_stress_matrix = temp_stress_matrix_2 / det_deformation_gradient;


            // transform stresses to local orthogonal base
            Matrix local_stress = ZeroMatrix(2);
            for (SizeType i=0;i<2;++i){
                for (SizeType j=0;j<2;++j){
                    local_stress(i,j) = inner_prod(transformed_base_vectors[i],prod(cauchy_stress_matrix,transformed_base_vectors[j]));
                }
            }
            stress = MathUtils<double>::StressTensorToVector(local_stress,3);


            if (rVariable==PRINCIPAL_CAUCHY_STRESS_VECTOR){
                Vector principal_stresses = ZeroVector(2);

                if (rOutput[point_number].size() != 2) {
                    rOutput[point_number].resize(2);
                }

                PrincipalVector(principal_stresses,stress);
                rOutput[point_number] = principal_stresses;
            }  else {
                rOutput[point_number] = stress;
            }
        }
    }
}


void MembraneElement::DeformationGradient(Matrix& rDeformationGradient, double& rDetDeformationGradient,
     const array_1d<Vector,2>& rCurrentCovariantBase, const array_1d<Vector,2>& rReferenceContraVariantBase)
{
    // attention: this is not in the local orthonogal coordinate system

    // calculate out of plane local vectors (membrane has no thickness change so g3 and G3 are normalized)
    Vector current_cov_3 = ZeroVector(3);
    MathUtils<double>::UnitCrossProduct(current_cov_3,rCurrentCovariantBase[0],rCurrentCovariantBase[1]);

    Vector reference_contra_3 = ZeroVector(3);
    MathUtils<double>::UnitCrossProduct(reference_contra_3,rReferenceContraVariantBase[0],rReferenceContraVariantBase[1]);

    // calculate deformation gradient
    rDeformationGradient = ZeroMatrix(3);
    for (SizeType i=0;i<2;++i){
        rDeformationGradient += outer_prod(rCurrentCovariantBase[i],rReferenceContraVariantBase[i]);
    }

    // add contribution of out of plane base vectors
    rDeformationGradient += outer_prod(current_cov_3,reference_contra_3);


    // calculate det(F)
    rDetDeformationGradient = MathUtils<double>::Det(rDeformationGradient);
}

void MembraneElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    // element with two nodes can only represent results at one node
    const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
    const SizeType& write_points_number =
        GetGeometry().IntegrationPointsNumber(integration_method);
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }

    if (rVariable == LOCAL_AXIS_1 || rVariable == LOCAL_AXIS_2 || rVariable == LOCAL_AXIS_3) {

        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);
        const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(integration_method);

        int which_axis = 0;
        if (rVariable == LOCAL_AXIS_2) which_axis = 1;

        array_1d<Vector,2> reference_covariant_base_vectors;
        array_1d<Vector,2> reference_contravariant_base_vectors;
        Matrix covariant_metric_reference = ZeroMatrix(3);
        Matrix contravariant_metric_reference = ZeroMatrix(3);
        array_1d<Vector,2> transformed_base_vectors;

        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){

            if (rOutput[point_number].size() != 3) {
                rOutput[point_number].resize(3);
            }

            // getting information for integration
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

            CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);
            CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);
            ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);
            ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);
            TransformBaseVectors(transformed_base_vectors,reference_contravariant_base_vectors);

            if (rVariable == LOCAL_AXIS_3){
                Vector base_vec_3 = ZeroVector(3);
                MathUtils<double>::UnitCrossProduct(base_vec_3,transformed_base_vectors[0],transformed_base_vectors[1]);

                for (SizeType i =0; i<3; ++i) {
                    rOutput[point_number][i] = base_vec_3[i];
                }
            }
            else {
                for (SizeType i =0; i<3; ++i) {
                    rOutput[point_number][i] = transformed_base_vectors[which_axis][i];
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void MembraneElement::Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == LOCAL_ELEMENT_ORIENTATION) {
        rOutput = ZeroMatrix(3);
        array_1d<Vector,2> base_vectors_current_cov;
        const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);
        const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(integration_method);

        Vector base_1 = ZeroVector(3);
        Vector base_2 = ZeroVector(3);
        Vector base_3 = ZeroVector(3);

        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
            const double integration_weight_i = r_integration_points[point_number].Weight();
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
            CovariantBaseVectors(base_vectors_current_cov,shape_functions_gradients_i,ConfigurationType::Reference);
            base_1 += base_vectors_current_cov[0]*integration_weight_i;
            base_2 += base_vectors_current_cov[1]*integration_weight_i;
        }

        MathUtils<double>::UnitCrossProduct(base_3, base_1, base_2);

        column(rOutput,0) = base_1;
        column(rOutput,1) = base_2;
        column(rOutput,2) = base_3;
    }
    else if (rVariable == MEMBRANE_PRESTRESS) {
        std::vector< Vector > prestress_matrix;
        CalculateOnIntegrationPoints(PK2_STRESS_VECTOR,prestress_matrix,rCurrentProcessInfo);
        const auto& r_integration_points = GetGeometry().IntegrationPoints(GetGeometry().GetDefaultIntegrationMethod());

        rOutput = ZeroMatrix(3,r_integration_points.size());

        // each column represents 1 GP
        for (SizeType i=0;i<r_integration_points.size();++i){
            column(rOutput,i) = prestress_matrix[i];
        }
    }
}

void MembraneElement::Calculate(const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == STRAIN_ENERGY) {
        const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
        const auto& r_geom = GetGeometry();

        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(integration_method);
        const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(integration_method);

        array_1d<Vector,2> current_covariant_base_vectors;
        array_1d<Vector,2> reference_covariant_base_vectors;
        array_1d<Vector,2> reference_contravariant_base_vectors;

        array_1d<Vector,2> transformed_base_vectors;

        Matrix covariant_metric_current = ZeroMatrix(3);
        Matrix covariant_metric_reference = ZeroMatrix(3);
        Matrix contravariant_metric_reference = ZeroMatrix(3);
        Matrix inplane_transformation_matrix_material = ZeroMatrix(3);
        double detJ = 0.0;
        rOutput = 0.0; // total strain energy
        Vector strain_vector = ZeroVector(3);
        Vector stress_vector = ZeroVector(3);


        ConstitutiveLaw::Parameters element_parameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        element_parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        element_parameters.SetStressVector(stress_vector);

        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
            // reset gauss point strain
            double strain_energy_gp = 0.0; // strain energy per gauss point

            // getting information for integration
            const double integration_weight_i = r_integration_points[point_number].Weight();
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];


            CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Current);
            CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);
            CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
            CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);

            ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);
            ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);

            TransformBaseVectors(transformed_base_vectors,reference_contravariant_base_vectors);
            InPlaneTransformationMatrix(inplane_transformation_matrix_material,transformed_base_vectors,reference_contravariant_base_vectors);
            JacobiDeterminante(detJ,reference_covariant_base_vectors);

            StrainGreenLagrange(strain_vector,covariant_metric_reference,covariant_metric_current,inplane_transformation_matrix_material);

            // add strain energy from material law
            element_parameters.SetStrainVector(strain_vector);
            mConstitutiveLawVector[point_number]->CalculateValue(element_parameters,STRAIN_ENERGY,strain_energy_gp);


            Vector pre_stress_vector = ZeroVector(3);
            AddPreStressPk2(pre_stress_vector,transformed_base_vectors);

            // add strain energy from pre_stress -> constant
            strain_energy_gp += inner_prod(strain_vector,pre_stress_vector);

            // integrate over reference domain
            strain_energy_gp *= detJ*integration_weight_i;

            // sum up the gauss point contributions
            rOutput += strain_energy_gp;
        }
        rOutput *= GetProperties()[THICKNESS];
    }
    else if (rVariable == KINETIC_ENERGY) {
        const auto& r_geom = GetGeometry();
        const SizeType number_dofs = r_geom.WorkingSpaceDimension()*r_geom.size();

        Matrix mass_matrix = ZeroMatrix(number_dofs,number_dofs);
        CalculateMassMatrix(mass_matrix,rCurrentProcessInfo);
        Vector current_nodal_velocities = ZeroVector(number_dofs);
        GetFirstDerivativesVector(current_nodal_velocities);
        rOutput = 0.50 * inner_prod(current_nodal_velocities,prod(mass_matrix,current_nodal_velocities));
    }
    else if (rVariable == ENERGY_DAMPING_DISSIPATION) {
        const auto& r_geom = GetGeometry();
        const SizeType number_dofs = r_geom.WorkingSpaceDimension()*r_geom.size();

        // Attention! this is only the current state and must be integrated over time (*dt)
        Matrix damping_matrix = ZeroMatrix(number_dofs,number_dofs);
        CalculateDampingMatrix(damping_matrix,rCurrentProcessInfo);
        Vector current_nodal_velocities = ZeroVector(number_dofs);
        GetFirstDerivativesVector(current_nodal_velocities);
        rOutput = inner_prod(current_nodal_velocities,prod(damping_matrix,current_nodal_velocities));
    }
    else if (rVariable == EXTERNAL_ENERGY) {
        // Dead Load contribution to external energy
        const auto& r_geom = GetGeometry();
        const SizeType number_dofs = r_geom.WorkingSpaceDimension()*r_geom.size();

        Vector dead_load_rhs = ZeroVector(number_dofs);
        CalculateAndAddBodyForce(dead_load_rhs, rCurrentProcessInfo);
        Vector current_nodal_displacements = ZeroVector(number_dofs);
        GetValuesVector(current_nodal_displacements, 0);
        rOutput = inner_prod(dead_load_rhs,current_nodal_displacements);
    }
}

void MembraneElement::CalculateConsistentMassMatrix(MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_dofs = dimension*number_of_nodes;

    if (number_of_nodes == 3){
        // consistent mass matrix for triangular element can be easily pre-computed
        const BoundedMatrix<double, 3, 3> fill_matrix = CalculateReferenceArea()*(IdentityMatrix(3)/12.0);
        for (SizeType i=0; i<3; ++i){
            for (SizeType j=0; j<3; ++j){
                    project(rMassMatrix, range((i*3),((i+1)*3)),range((j*3),((j+1)*3))) += fill_matrix;
                }
        }
        for (SizeType i=0; i<number_dofs; ++i) rMassMatrix(i,i) *= 2.0;
    }
    else {
        const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);

        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(integration_method);
        const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(integration_method);
        const Matrix& rNcontainer = r_geom.ShapeFunctionsValues(integration_method);
        array_1d<Vector,2> reference_covariant_base_vectors;

        double detJ = 0.0;



        for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){

            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
            CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);
            JacobiDeterminante(detJ,reference_covariant_base_vectors);

            const double integration_weight = r_integration_points[point_number].Weight();
            const Vector& rN = row(rNcontainer,point_number);


            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                    const SizeType index_i = i * dimension;

                    for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                        const SizeType index_j = j * dimension;
                        const double NiNj_weight = rN[i] * rN[j] * integration_weight * detJ;

                        for ( IndexType k = 0; k < dimension; ++k )
                            rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
                    }
                }
        }
    }

    rMassMatrix *= GetProperties()[THICKNESS]*StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    KRATOS_CATCH("");
}

void MembraneElement::CalculateMassMatrix(MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();

    // LUMPED MASS MATRIX
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != mat_size) {
        rMassMatrix.resize(mat_size, mat_size, false);
    }
    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);


    if (StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(GetProperties(), rCurrentProcessInfo) ){
        Vector lumped_mass_vector = ZeroVector(mat_size);
        CalculateLumpedMassVector(lumped_mass_vector,rCurrentProcessInfo);
        for (SizeType i=0;i<mat_size;++i) rMassMatrix(i,i) = lumped_mass_vector[i];
    }
    else {
        // CONSISTENT MASS MATRIX
        CalculateConsistentMassMatrix(rMassMatrix,rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

void MembraneElement::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType local_size = dimension*number_of_nodes;

    if (rLumpedMassVector.size() != local_size) {
        rLumpedMassVector.resize(local_size, false);
    }

    const double total_mass = CalculateReferenceArea() * GetProperties()[THICKNESS] * StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

    Vector lump_fact =  ZeroVector(number_of_nodes);
    ReferenceLumpingFactors(lump_fact);

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const double temp = lump_fact[i] * total_mass;

        for (SizeType j = 0; j < 3; ++j)
        {
            const SizeType index = i * 3 + j;
            rLumpedMassVector[index] = temp;
        }
    }
    KRATOS_CATCH("")
}

void MembraneElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType local_size = dimension*number_of_nodes;

    if (rDestinationVariable == NODAL_MASS) {
        VectorType element_mass_vector(local_size);
        CalculateLumpedMassVector(element_mass_vector, rCurrentProcessInfo);

        for (SizeType i = 0; i < number_of_nodes; ++i) {
            double& r_nodal_mass = r_geom[i].GetValue(NODAL_MASS);
            int index = i * dimension;

            AtomicAdd(r_nodal_mass, element_mass_vector(index));
        }
    }

    KRATOS_CATCH("")
}

void MembraneElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        GetGeometry().WorkingSpaceDimension()*GetGeometry().size());
}

const Parameters MembraneElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static","implicit","explicit"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["DISPLACEMENT","ROTATION","VELOCITY","ACCELERATION"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","ROTATION"],
        "required_dofs"              : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","ROTATION_X","ROTATION_Y","ROTATION_Z"],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle3D3", "Quadrilateral3D4"],
        "element_integrates_in_time" : false,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStress"],
            "dimension"   : ["3D"],
            "strain_size" : [3]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   : "This element implements a pre-stressed membrane formulation."
    })");

    return specifications;
}

void MembraneElement::AddExplicitContribution(
    const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType local_size = dimension*number_of_nodes;

    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {

        Vector damping_residual_contribution = ZeroVector(local_size);
        Vector current_nodal_velocities = ZeroVector(local_size);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix;
        CalculateDampingMatrix(damping_matrix, rCurrentProcessInfo);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);

        for (SizeType i = 0; i < number_of_nodes; ++i) {
            SizeType index = dimension * i;
            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < dimension; ++j) {
                AtomicAdd(r_force_residual[j], (rRHSVector[index + j] - damping_residual_contribution[index + j]));
            }
        }
    } else if (rDestinationVariable == NODAL_INERTIA) {

        // Getting the vector mass
        VectorType mass_vector(local_size);
        CalculateLumpedMassVector(mass_vector, rCurrentProcessInfo);

        for (SizeType i = 0; i < number_of_nodes; ++i) {
            double& r_nodal_mass = GetGeometry()[i].GetValue(NODAL_MASS);
            SizeType index = i * dimension;

            AtomicAdd(r_nodal_mass, mass_vector[index]);
        }
    }

    KRATOS_CATCH("")
}

void MembraneElement::CalculateAndAddBodyForce(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType local_size = dimension*number_of_nodes;

    if (r_geom[0].SolutionStepsDataHas(VOLUME_ACCELERATION)){

        Vector lumped_mass_vector = ZeroVector(local_size);
        CalculateLumpedMassVector(lumped_mass_vector,rCurrentProcessInfo);

        for (SizeType i = 0; i < number_of_nodes; ++i) {
            for (SizeType j = 0; j < 3; ++j)
            {
                const SizeType index = i * 3 + j;
                rRightHandSideVector[index] += lumped_mass_vector[index] * r_geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION)[j];
            }
        }
    }

    KRATOS_CATCH("")
}

double MembraneElement::CalculateReferenceArea() const
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    const IntegrationMethod integration_method = GetIntegrationMethod();

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(integration_method);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(integration_method);
    array_1d<Vector,2> reference_covariant_base_vectors;

    double detJ = 0.0;
    double ref_area = 0.0;

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){

        const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
        CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference);
        JacobiDeterminante(detJ,reference_covariant_base_vectors);
        const double integration_weight = r_integration_points[point_number].Weight();

        ref_area += integration_weight * detJ;
    }
    return ref_area;
    KRATOS_CATCH("");
}



void MembraneElement::PrincipalVector(Vector& rPrincipalVector, const Vector& rNonPrincipalVector)
{
    // make sure to divide rNonPrincipalVector[2]/2 if strains are passed
    rPrincipalVector = ZeroVector(2);
    rPrincipalVector[0] = 0.50 * (rNonPrincipalVector[0]+rNonPrincipalVector[1]) + std::sqrt(0.25*(std::pow(rNonPrincipalVector[0]-rNonPrincipalVector[1],2.0)) + std::pow(rNonPrincipalVector[2],2.0));
    rPrincipalVector[1] = 0.50 * (rNonPrincipalVector[0]+rNonPrincipalVector[1]) - std::sqrt(0.25*(std::pow(rNonPrincipalVector[0]-rNonPrincipalVector[1],2.0)) + std::pow(rNonPrincipalVector[2],2.0));
}

//***********************************************************************************
//***********************************************************************************
int MembraneElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo[DOMAIN_SIZE]==3) << "DOMAIN_SIZE in element " << Id() << " is not 3" << std::endl;
    KRATOS_ERROR_IF_NOT(dimension==3) << "dimension in element " << Id() << " is not 3" << std::endl;

    if (GetProperties().Has(THICKNESS) == false ||
            GetProperties()[THICKNESS] <= numerical_limit) {
        KRATOS_ERROR << "THICKNESS not provided for element " << Id()
                     << std::endl;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( SizeType i = 0; i < number_of_nodes; i++ ) {
        const Node &r_node = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW ))
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->Check(GetProperties(),GetGeometry(),rCurrentProcessInfo);
            const SizeType strain_size = mConstitutiveLawVector[point_number]->GetStrainSize();
            KRATOS_ERROR_IF( strain_size != 3) << "Wrong constitutive law used. This is a membrane element! "
                << "Expected strain size is 3 (el id = " << this->Id() << ")" << std::endl;
        }
    } else KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    // Verify that the constitutive law has the correct dimension


    return 0;
    KRATOS_CATCH("");
}

void MembraneElement::save(Serializer& rSerializer) const
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
      rSerializer.save("mConstitutiveLawVector", mConstitutiveLawVector);
    }

    void MembraneElement::load(Serializer& rSerializer)
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
      rSerializer.load("mConstitutiveLawVector", mConstitutiveLawVector);
    }



//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos.
