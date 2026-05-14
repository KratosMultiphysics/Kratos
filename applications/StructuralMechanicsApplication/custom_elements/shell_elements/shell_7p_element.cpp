// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//

// System includes

// External includes

// Project includes
#include "shell_7p_element.hpp"

//READ SOLVER SCRIPT*********************************************************************************************************************************
// https://github.com/KratosMultiphysics/Kratos/blob/master/applications/StructuralMechanicsApplication/python_scripts/structural_mechanics_solver.py
//***************************************************************************************************************************************************
namespace Kratos
{

Shell7pElement::Shell7pElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

Shell7pElement::Shell7pElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    KRATOS_WATCH("new element")
}

Element::Pointer Shell7pElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
) const
{
    return Kratos::make_intrusive<Shell7pElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer Shell7pElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
) const
{
    return Kratos::make_intrusive<Shell7pElement>(NewId, pGeom, pProperties);
}

void Shell7pElement::GetDofList(
DofsVectorType& rElementalDofList,
const ProcessInfo& rCurrentProcessInfo
) const
{
KRATOS_TRY

const auto& r_geom = GetGeometry();
const SizeType number_of_nodes = r_geom.size();
const SizeType local_size = number_of_nodes * 6;

if (rElementalDofList.size() != local_size) {
    rElementalDofList.resize(local_size);
}

for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node) {
    const SizeType index = 6 * i_node;
    rElementalDofList[index]     = r_geom[i_node].pGetDof(DISPLACEMENT_X);
    rElementalDofList[index + 1] = r_geom[i_node].pGetDof(DISPLACEMENT_Y);
    rElementalDofList[index + 2] = r_geom[i_node].pGetDof(DISPLACEMENT_Z);
    rElementalDofList[index + 3] = r_geom[i_node].pGetDof(ROTATION_X);
    rElementalDofList[index + 4] = r_geom[i_node].pGetDof(ROTATION_Y);
    rElementalDofList[index + 5] = r_geom[i_node].pGetDof(ROTATION_Z);
}

KRATOS_CATCH("");
}

void Shell7pElement::EquationIdVector(
EquationIdVectorType& rResult,
const ProcessInfo& rCurrentProcessInfo
) const
{
KRATOS_TRY

const auto& r_geom = GetGeometry();
const SizeType number_of_nodes = r_geom.size();
const SizeType local_size = number_of_nodes * 6;
const SizeType u_dof = r_geom[0].GetDofPosition(DISPLACEMENT_X);
const SizeType dir_dof = r_geom[0].GetDofPosition(ROTATION_X);

if (rResult.size() != local_size) {
    rResult.resize(local_size);
}

for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node) {
    const SizeType index = 6 * i_node;
    rResult[index]     = r_geom[i_node].GetDof(DISPLACEMENT_X, u_dof).EquationId();
    rResult[index + 1] = r_geom[i_node].GetDof(DISPLACEMENT_Y, u_dof + 1).EquationId();
    rResult[index + 2] = r_geom[i_node].GetDof(DISPLACEMENT_Z, u_dof + 2).EquationId();
    rResult[index + 3] = r_geom[i_node].GetDof(ROTATION_X, dir_dof).EquationId();
    rResult[index + 4] = r_geom[i_node].GetDof(ROTATION_Y, dir_dof + 1).EquationId();
    rResult[index + 5] = r_geom[i_node].GetDof(ROTATION_Z, dir_dof + 2).EquationId();
}

KRATOS_CATCH("");
}

int Shell7pElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
KRATOS_TRY;

const double numerical_limit = std::numeric_limits<double>::epsilon();

if (GetProperties().Has(THICKNESS) == false ||
        GetProperties()[THICKNESS] <= numerical_limit) {
    KRATOS_ERROR << "THICKNESS not provided for element " << Id()
                    << std::endl;
}

return 0;

KRATOS_CATCH("");
}

//***********************************************************************************
//***********************************************************************************

void Shell7pElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void Shell7pElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = 6 * number_of_nodes;

    if (rRightHandSideVector.size() != number_dofs) {
        rRightHandSideVector.resize(number_dofs, false);
    }

    MatrixType LHS;
    CalculateLeftHandSide(LHS, rCurrentProcessInfo);

    Vector current_values = ZeroVector(number_dofs);
    for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node) {
        const SizeType index = 6 * i_node;
        const auto& u_dof = r_geom[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const auto& dir_dof = r_geom[i_node].FastGetSolutionStepValue(ROTATION);

        current_values[index] = u_dof[0];
        current_values[index + 1] = u_dof[1];
        current_values[index + 2] = u_dof[2];
        current_values[index + 3] = dir_dof[0];
        current_values[index + 4] = dir_dof[1];
        current_values[index + 5] = dir_dof[2];
    }

    noalias(rRightHandSideVector) = -prod(LHS, current_values);
}

void Shell7pElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)

{
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = 6*number_of_nodes;
    rLeftHandSideMatrix = ZeroMatrix(number_dofs);

    const IntegrationMethod integration_method = r_geom.GetDefaultIntegrationMethod();
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);     // Rows: GP, Columns: element nodes. Ncontainer(k,i) = N_i evaluated at GP k. For Quad4: Ncontainer(k,i) = N_i evaluated at GP k, where i=0..3 and k=0..3
    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(integration_method); // Stack of matrices: one for each GP. [GP][N_node][derivative direction] Quad4: For GP k: node 0 [ dN0/dξ  dN0/dη ]
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(integration_method);                                                                                                            // node 1 [ dN1/dξ  dN1/dη ]
                                                                                                                                                                                                                            // node 2 [ dN2/dξ  dN2/dη ]
                                                                                                                                                                                                                            // node 3 [ dN3/dξ  dN3/dη ]
    const double thickness = GetProperties()[THICKNESS];                // GetProperties() returns the Properties object assigned to the element, not to individual nodes. All nodes of the element share the same thickness value.

    array_1d<Vector,3> current_covariant_base_vectors;
    // array_1d<Vector,2> a3kvp;    
    array_1d<Vector,3> reference_covariant_base_vectors;    // array_1d<array_1d<double,3>,3> reference_covariant_base_vectors;    // outer and inner sizes are compile-time fixed to 3 for better performance (no dynamic memory allocation). access: reference_covariant_base_vectors[i][j]: i = which base vector (0=g1, 1=g2, 2=g3)
    array_1d<Vector,3> reference_contravariant_base_vectors;                                                                                                                                                         // j = spatial component  (0=x,  1=y,  2=z)
    array_1d<Vector,3> transformed_base_vectors;

    Matrix covariant_metric_current = ZeroMatrix(3);        // BoundedMatrix<double, 3, 3> covariant_metric_current = ZeroMatrix(3, 3); // i cannot resize() it — size is fixed at compile-time
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
        const Vector& Nshape = row(Ncontainer,point_number);        // Node shape function values at the current integration point. Nshape[i] = N_i evaluated at the current GP

        // CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Current,thickness);
        CovariantBaseVectors(reference_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Reference,thickness);
        // DirectorDerivatives(a3kvp,reference_covariant_base_vectors,shape_functions_gradients_i);

        // CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
        CovariantMetric(covariant_metric_reference,reference_covariant_base_vectors);
        ContravariantMetric(contravariant_metric_reference,covariant_metric_reference);
        // later: implement partial derivatives of G3 with respect to inplane parametric coordinates

        ContraVariantBaseVectors(reference_contravariant_base_vectors,contravariant_metric_reference,reference_covariant_base_vectors);

        JacobiDeterminante(detJ,reference_covariant_base_vectors);

        Matrix material_tangent_modulus = ZeroMatrix(6);
        BoundedMatrix<double, 12, 12> Dmatrix=ZeroMatrix(12,12);
        Matrix Bop = ZeroMatrix(12,number_dofs);             // DOFs vary by geometry type
        CalculateMaterialLaw(Dmatrix,contravariant_metric_reference,thickness,ConstitutiveLawType::gStVenantKirchhoff);
        CalculatelinearBOperator(Bop,reference_covariant_base_vectors,shape_functions_gradients_i,Nshape,number_of_nodes);

        double weight = integration_weight_i * detJ * thickness*0.5; 
        Matrix DB = ZeroMatrix(12,number_dofs); 
        noalias(DB) = prod(Dmatrix, Bop);
        rLeftHandSideMatrix += prod(trans(Bop), DB) * weight;
    }
}

void Shell7pElement::CovariantBaseVectors(array_1d<Vector,3>& rBaseVectors,
     const Matrix& rShapeFunctionGradientValues, const ConfigurationType& rConfiguration, const double& thickness) const
{
    // pass/call this ShapeFunctionsLocalGradients[pnt]
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().size();
    Vector g1 = ZeroVector(dimension);
    Vector g2 = ZeroVector(dimension);
    Vector g3 = ZeroVector(dimension);
    
    // Vector current_displacement = ZeroVector(dimension*number_of_nodes);
    //if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);


    for (SizeType i=0;i<number_of_nodes;++i){
        g1[0] += GetGeometry().GetPoint( i ).X0() * rShapeFunctionGradientValues(i, 0);
        g1[1] += GetGeometry().GetPoint( i ).Y0() * rShapeFunctionGradientValues(i, 0);
        g1[2] += GetGeometry().GetPoint( i ).Z0() * rShapeFunctionGradientValues(i, 0);

        g2[0] += GetGeometry().GetPoint( i ).X0() * rShapeFunctionGradientValues(i, 1);
        g2[1] += GetGeometry().GetPoint( i ).Y0() * rShapeFunctionGradientValues(i, 1);
        g2[2] += GetGeometry().GetPoint( i ).Z0() * rShapeFunctionGradientValues(i, 1);
    }
    MathUtils<double>::UnitCrossProduct(g3,g1,g2);          // g3 is consistent only if all elements have the same node ordering convention (all counter-clockwise or all clockwise when viewed from one side). i can "trust" that the mesher produces consistently ordered elements, or add a check/correction in Initialize() by comparing g3 against a reference direction if needed
    rBaseVectors[0] = g1;
    rBaseVectors[1] = g2;
    rBaseVectors[2] = g3*thickness*0.5;   
}

// void Shell7pElement::DirectorDerivatives(array_1d<Vector,2>& rDirectorDerivatives,const array_1d<Vector,3>& rBaseVectorCovariant,const Matrix& rShapeFunctionGradientValues) const


void Shell7pElement::CovariantMetric(Matrix& rMetric,const array_1d<Vector,3>& rBaseVectorCovariant)
{
    rMetric = ZeroMatrix(3);
    for (SizeType i=0;i<2;++i){
        for (SizeType j=0;j<2;++j){
            rMetric(i,j) = inner_prod(rBaseVectorCovariant[i],rBaseVectorCovariant[j]);
        }
    }
    rMetric(2,2)=inner_prod(rBaseVectorCovariant[2],rBaseVectorCovariant[2]); // later change this!! since g3 could be non othogonal to g1 and g2
}

void Shell7pElement::ContravariantMetric(Matrix& rMetric,const Matrix& rCovariantMetric)
{
    rMetric = ZeroMatrix(3);
    rMetric(0,0) = rCovariantMetric(1,1);
    rMetric(1,1) = rCovariantMetric(0,0);
    rMetric(0,1) = -1.0*rCovariantMetric(1,0);
    rMetric(1,0) = -1.0*rCovariantMetric(0,1);
    rMetric /= (rCovariantMetric(1,1)*rCovariantMetric(0,0)) - (rCovariantMetric(1,0)*rCovariantMetric(0,1));
    rMetric(2,2)=1.0/rCovariantMetric(2,2); 
}
    // later change this to:
// void Shell7pElement::ContravariantMetric(Matrix& rMetric, const Matrix& rCovariantMetric)
// {
//     rMetric = ZeroMatrix(3);
//     double det = MathUtils<double>::InvertMatrix(rCovariantMetric, rMetric);         // 1.Uses the general 3×3 inversion 2.Checks determinant 3.Works even if orthogonality assumption fails 4.is clearer and maintainable

//     if (std::abs(det) < 1e-14) {                                                     // check the determinant for singularity (near-zero = bad matrix condition):
//         KRATOS_ERROR << "Singular covariant metric detected. det = " << det << std::endl;
//     }
// }

void Shell7pElement::ContraVariantBaseVectors(array_1d<Vector,3>& rBaseVectors,const Matrix& rContraVariantMetric,
    const array_1d<Vector,3> rCovariantBaseVectors)
{
    rBaseVectors[0] = ZeroVector(3);
    rBaseVectors[1] = ZeroVector(3);
    rBaseVectors[2] = ZeroVector(3);

    rBaseVectors[0] = rContraVariantMetric(0,0)*rCovariantBaseVectors[0] + rContraVariantMetric(0,1)*rCovariantBaseVectors[1];
    rBaseVectors[1] = rContraVariantMetric(1,0)*rCovariantBaseVectors[0] + rContraVariantMetric(1,1)*rCovariantBaseVectors[1];
    rBaseVectors[2] = rContraVariantMetric(2,2)*rCovariantBaseVectors[2]; // later change this!!
    // rBaseVectors[0] = rContraVariantMetric(0,0)*rCovariantBaseVectors[0] + rContraVariantMetric(0,1)*rCovariantBaseVectors[1] + rContraVariantMetric(0,2)*rCovariantBaseVectors[2];
    // rBaseVectors[1] = rContraVariantMetric(1,0)*rCovariantBaseVectors[0] + rContraVariantMetric(1,1)*rCovariantBaseVectors[1] + rContraVariantMetric(1,2)*rCovariantBaseVectors[2];
    // rBaseVectors[2] = rContraVariantMetric(2,0)*rCovariantBaseVectors[0] + rContraVariantMetric(2,1)*rCovariantBaseVectors[1] + rContraVariantMetric(2,2)*rCovariantBaseVectors[2]; 
}

void Shell7pElement::JacobiDeterminante(double& rDetJacobi, const array_1d<Vector,3>& rReferenceBaseVectors) const
{
    array_1d<double, 3> g3;
    MathUtils<double>::CrossProduct(g3, rReferenceBaseVectors[0], rReferenceBaseVectors[1]);
    rDetJacobi = MathUtils<double>::Norm(g3);
    KRATOS_ERROR_IF(rDetJacobi<std::numeric_limits<double>::epsilon()) << "det of Jacobi smaller 0 for element with id" << Id() << std::endl;
}

void Shell7pElement::CalculateMaterialLaw(BoundedMatrix<double, 12, 12>& CL, const Matrix& Gmkon, const double& thickness,
const ConstitutiveLawType& option)
{
    const auto& r_properties = GetProperties();
    const double E = r_properties[YOUNG_MODULUS];
    const double nu = r_properties[POISSON_RATIO];
    const double G = E/(2.0*(1.0+nu));
    const double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    const double mu = G;
    CL = ZeroMatrix(12);

    if (option == ConstitutiveLawType::gStVenantKirchhoff) {
        double C[3][3][3][3] = {};
        double Theta3[2] = {1.0 / sqrt(3.0), -1.0 / sqrt(3.0)};
        double gpweight[2] = {1.0, 1.0};
        BoundedMatrix<double, 6, 6> CC = ZeroMatrix(6);

        for (SizeType i=0; i<3; ++i){
            for (SizeType j=0; j<3; ++j){
                for (SizeType k=0; k<3; ++k){
                    for (SizeType l=0; l<3; ++l){
                        C[i][j][k][l] = lambda*Gmkon(i,j)*Gmkon(k,l) + mu*( Gmkon(i,k)*Gmkon(j,l) + Gmkon(i,l)*Gmkon(k,j) );
                    }
                }
            }
        }
        CC(0,0) = C[0][0][0][0];
        CC(0,1) = C[0][0][1][0];
        CC(0,2) = C[0][0][2][0];            // sigma 11
        CC(0,3) = C[0][0][1][1];
        CC(0,4) = C[0][0][2][1];
        CC(0,5) = C[0][0][2][2];

        CC(1,0) = C[1][0][0][0];
        CC(1,1) = C[1][0][1][0];
        CC(1,2) = C[1][0][2][0];            // sigma 12
        CC(1,3) = C[1][0][1][1];
        CC(1,4) = C[1][0][2][1];
        CC(1,5) = C[1][0][2][2];

        CC(2,0) = C[2][0][0][0];
        CC(2,1) = C[2][0][1][0];
        CC(2,2) = C[2][0][2][0]*5.0/6.0;     // sigma 13 with shear correction factor alpha=5/6 for E13 and E23
        CC(2,3) = C[2][0][1][1];
        CC(2,4) = C[2][0][2][1]*5.0/6.0;
        CC(2,5) = C[2][0][2][2];

        CC(3,0) = C[1][1][0][0];
        CC(3,1) = C[1][1][1][0];
        CC(3,2) = C[1][1][2][0];            // sigma 22
        CC(3,3) = C[1][1][1][1];
        CC(3,4) = C[1][1][2][1];
        CC(3,5) = C[1][1][2][2];

        CC(4,0) = C[2][1][0][0];
        CC(4,1) = C[2][1][1][0];
        CC(4,2) = C[2][1][2][0]*5.0/6.0;    // sigma 23 with shear correction factor alpha=5/6 for E13 and E23
        CC(4,3) = C[2][1][1][1];
        CC(4,4) = C[2][1][2][1]*5.0/6.0;
        CC(4,5) = C[2][1][2][2];

        CC(5,0) = C[2][2][0][0];
        CC(5,1) = C[2][2][1][0];
        CC(5,2) = C[2][2][2][0];            // sigma 33
        CC(5,3) = C[2][2][1][1];
        CC(5,4) = C[2][2][2][1];
        CC(5,5) = C[2][2][2][2];

        for (SizeType k=0; k<2; ++k){
            for (SizeType i=0; i<6; ++i){
                const SizeType i6 = i + 6;

                for (SizeType j=0; j<6; ++j){
                    const SizeType j6 = j + 6;

                    CL(i,j) += CC(i,j)*gpweight[k];
                    CL(i6,j) += CC(i,j)*Theta3[k]*gpweight[k];
                    CL(j,i6) += CC(j,i)*Theta3[k]*gpweight[k]; 
                    CL(i6,j6) += CC(i,j)*Theta3[k]*Theta3[k]*gpweight[k];
                 }
            }
        }
       // CL(2,2) *= 5.0/6.0;    // shear correction factor alpha=5/6 for n13,n23
       // CL(2,4) *= 5.0/6.0;
       // CL(4,2) *= 5.0/6.0;
       // CL(4,4) *= 5.0/6.0;
       // CL(8,8) *= 0.7;        // shear correction factor betta=0.7 for m13,m23
       // CL(8,10) *= 0.7;
       // CL(10,8) *= 0.7;
       // CL(10,10) *= 0.7;
    }
 else {
    const double Ebar = E*(1.0-nu)/((1.0+nu)*(1.0-2.0*nu));
    const double hbar = thickness*thickness*thickness/12.0;
    const double hq = thickness*5.0/6.0;
    const double hq_bar = 0.7*hbar;
        CL = ZeroMatrix(12);
        CL(0,0) = Ebar*thickness;
        CL(1,1) = Ebar*thickness;
        CL(2,2) = Ebar*thickness;
        CL(0,1) = lambda*thickness;
        CL(0,2) = lambda*thickness;
        CL(1,0) = lambda*thickness;
        CL(1,2) = lambda*thickness;             // if we use this, then we need to reorder this cartesian consitutive or B-matrix
        CL(2,0) = lambda*thickness;
        CL(2,1) = lambda*thickness;
        CL(3,3) = G*thickness;
        CL(4,4) = G*hq;
        CL(5,5) = G*hq;
        CL(6,6) = Ebar*hbar;
        CL(7,7) = Ebar*hbar;
        CL(8,8) = Ebar*hbar;
        CL(6,7) = lambda*hbar;
        CL(6,8) = lambda*hbar;
        CL(7,6) = lambda*hbar;
        CL(7,8) = lambda*hbar;
        CL(8,6) = lambda*hbar;
        CL(8,7) = lambda*hbar;
        CL(9,9) = G*hbar;
        CL(10,10) = G*hq_bar;
        CL(11,11) = G*hq_bar;
    }


}
                                                                                                           // [N_node][derivative direction]          // Nshape[i] = N_i evaluated at the current GP
void Shell7pElement::CalculatelinearBOperator(Matrix& bop, const array_1d<Vector,3>& CovariantBaseVectors, const Matrix& ShapeFunctionGradientValues, const Vector& Nshape, const SizeType& number_of_nodes)
{
const double a1x=CovariantBaseVectors[0][0];
const double a1y=CovariantBaseVectors[0][1];
const double a1z=CovariantBaseVectors[0][2];
const double a2x=CovariantBaseVectors[1][0];
const double a2y=CovariantBaseVectors[1][1];
const double a2z=CovariantBaseVectors[1][2];
const double a3x=CovariantBaseVectors[2][0];
const double a3y=CovariantBaseVectors[2][1];
const double a3z=CovariantBaseVectors[2][2];

    for (SizeType i=0;i<number_of_nodes;++i){
            const SizeType index = i*6;
            const double dNd1 = ShapeFunctionGradientValues(i,0);
            const double dNd2 = ShapeFunctionGradientValues(i,1);
            const double N = Nshape[i];

            bop(0,index)   = dNd1*a1x;
            bop(0,index+1) = dNd1*a1y;
            bop(0,index+2) = dNd1*a1z;                  // alpha 11
            bop(0,index+3) = 0.0;
            bop(0,index+4) = 0.0;
            bop(0,index+5) = 0.0;

            bop(1,index)   = dNd2*a1x + dNd1*a2x;
            bop(1,index+1) = dNd2*a1y + dNd1*a2y;
            bop(1,index+2) = dNd2*a1z + dNd1*a2z;       // alpha 12
            bop(1,index+3) = 0.0;
            bop(1,index+4) = 0.0;
            bop(1,index+5) = 0.0;

            bop(2,index)   = dNd1*a3x;
            bop(2,index+1) = dNd1*a3y;
            bop(2,index+2) = dNd1*a3z;                  // alpha 13
            bop(2,index+3) = N*a1x;
            bop(2,index+4) = N*a1y;
            bop(2,index+5) = N*a1z;

            bop(3,index)   = dNd2*a2x;
            bop(3,index+1) = dNd2*a2y;
            bop(3,index+2) = dNd2*a2z;                  // alpha 22
            bop(3,index+3) = 0.0;
            bop(3,index+4) = 0.0;
            bop(3,index+5) = 0.0;

            bop(4,index)   = dNd2*a3x;
            bop(4,index+1) = dNd2*a3y;
            bop(4,index+2) = dNd2*a3z;                  // alpha 23
            bop(4,index+3) = N*a2x;
            bop(4,index+4) = N*a2y;
            bop(4,index+5) = N*a2z;

            bop(5,index)   = 0.0;
            bop(5,index+1) = 0.0;
            bop(5,index+2) = 0.0;                       // alpha 33
            bop(5,index+3) = N*a3x;
            bop(5,index+4) = N*a3y;
            bop(5,index+5) = N*a3z;

            bop(6,index)   = 0.0;
            bop(6,index+1) = 0.0;
            bop(6,index+2) = 0.0;                       // betta 11
            bop(6,index+3) = dNd1*a1x;
            bop(6,index+4) = dNd1*a1y;
            bop(6,index+5) = dNd1*a1z;

            bop(7,index)   = 0.0;
            bop(7,index+1) = 0.0;
            bop(7,index+2) = 0.0;                       // betta 12
            bop(7,index+3) = dNd2*a1x + dNd1*a2x;
            bop(7,index+4) = dNd2*a1y + dNd1*a2y;
            bop(7,index+5) = dNd2*a1z + dNd1*a2z;

            bop(8,index)   = 0.0;
            bop(8,index+1) = 0.0;
            bop(8,index+2) = 0.0;                       // betta 13
            bop(8,index+3) = dNd1*a3x;
            bop(8,index+4) = dNd1*a3y;
            bop(8,index+5) = dNd1*a3z;

            bop(9,index)   = 0.0;
            bop(9,index+1) = 0.0;
            bop(9,index+2) = 0.0;                       // betta 22
            bop(9,index+3) = dNd2*a2x;
            bop(9,index+4) = dNd2*a2y;
            bop(9,index+5) = dNd2*a2z;

            bop(10,index)   = 0.0;
            bop(10,index+1) = 0.0;
            bop(10,index+2) = 0.0;                       // betta 23
            bop(10,index+3) = dNd2*a3x;
            bop(10,index+4) = dNd2*a3y;
            bop(10,index+5) = dNd2*a3z;

            bop(11,index)   = 0.0;
            bop(11,index+1) = 0.0;
            bop(11,index+2) = 0.0;                       // betta 33
            bop(11,index+3) = 0.0;
            bop(11,index+4) = 0.0;
            bop(11,index+5) = 0.0;
    }
}

} // namespace Kratos