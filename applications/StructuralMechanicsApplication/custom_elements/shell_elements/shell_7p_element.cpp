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
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"

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
    rElementalDofList[index + 3] = r_geom[i_node].pGetDof(DIFFERENCE_DIRECTOR_X);
    rElementalDofList[index + 4] = r_geom[i_node].pGetDof(DIFFERENCE_DIRECTOR_Y);
    rElementalDofList[index + 5] = r_geom[i_node].pGetDof(DIFFERENCE_DIRECTOR_Z);
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
const SizeType v_dof = r_geom[0].GetDofPosition(DISPLACEMENT_X);
const SizeType w_dof = r_geom[0].GetDofPosition(DIFFERENCE_DIRECTOR_X);

if (rResult.size() != local_size) {
    rResult.resize(local_size);
}

for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node) {
    const SizeType index = 6 * i_node;
    rResult[index]     = r_geom[i_node].GetDof(DISPLACEMENT_X, v_dof).EquationId();
    rResult[index + 1] = r_geom[i_node].GetDof(DISPLACEMENT_Y, v_dof + 1).EquationId();
    rResult[index + 2] = r_geom[i_node].GetDof(DISPLACEMENT_Z, v_dof + 2).EquationId();
    rResult[index + 3] = r_geom[i_node].GetDof(DIFFERENCE_DIRECTOR_X, w_dof).EquationId();
    rResult[index + 4] = r_geom[i_node].GetDof(DIFFERENCE_DIRECTOR_Y, w_dof + 1).EquationId();
    rResult[index + 5] = r_geom[i_node].GetDof(DIFFERENCE_DIRECTOR_Z, w_dof + 2).EquationId();
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
        const auto& v_dof = r_geom[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const auto& w_dof = r_geom[i_node].FastGetSolutionStepValue(DIFFERENCE_DIRECTOR); // reaction for "ROTATION" variable is a REACTION_MOMENT, and they are computed at the moment not for differential vector, so results for moment are different than should be. vtk output label difference director as "rotation" and plot it as such

        current_values[index] = v_dof[0];
        current_values[index + 1] = v_dof[1];
        current_values[index + 2] = v_dof[2];
        current_values[index + 3] = w_dof[0];
        current_values[index + 4] = w_dof[1];
        current_values[index + 5] = w_dof[2];
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
    const double thickness = GetProperties()[THICKNESS];                // All nodes of the element share the same thickness value, since GetProperties() returns the Properties object assigned to the element, not to individual nodes

    array_1d<Vector,3> current_covariant_base_vectors;
    array_1d<Vector,3> akovr;    // array_1d<array_1d<double,3>,3> akovr;    // outer and inner sizes are compile-time fixed to 3 for better performance (no dynamic memory allocation). access: reference_covariant_base_vectors[i][j]: i = which base vector (0=g1, 1=g2, 2=g3)
    array_1d<Vector,3> gkovr;
    array_1d<Vector,3> akonr;                                                                                                                                                         // j = spatial component  (0=x,  1=y,  2=z)
    array_1d<Vector,2> a3kvp;       

    Matrix covariant_metric_current = ZeroMatrix(3);        // BoundedMatrix<double, 3, 3> covariant_metric_current = ZeroMatrix(3, 3); // i cannot resize() it — size is fixed at compile-time
    Matrix amkovr = ZeroMatrix(3);
    Matrix amkonr = ZeroMatrix(3);
    Matrix gmkovr = ZeroMatrix(3);
    Matrix gmkonr = ZeroMatrix(3);
    double detJ_surface = 0.0;

    array_1d<double,2> gpcoord3;
    gpcoord3[0] =  1.0 / std::sqrt(3.0);
    gpcoord3[1] = -1.0 / std::sqrt(3.0);
    array_1d<double,2> gpweight;
    gpweight[0] = 1.0;
    gpweight[1] = 1.0;

    double amdet_body = 0.0;
    double gmdet_body = 0.0;

    ////////////////////////////////////////////////////////////////BEGIN ANS STUFF////////////////////////////////////////////////////////////////
    int ansq=1; // activation flag of Q-mode of ANS  
    SizeType n_ans_points=4; 
    // 4 ANS sampling points: (r, s)
    array_1d<array_1d<double,2>,4> ans_points;
    ans_points[0][0] =  0.0; ans_points[0][1] = -1.0; // (xr1[0], xs1[0])
    ans_points[1][0] =  0.0; ans_points[1][1] =  1.0; // (xr1[1], xs1[1])
    ans_points[2][0] = -1.0; ans_points[2][1] =  0.0; // (xr2[0], xs2[0])
    ans_points[3][0] =  1.0; ans_points[3][1] =  0.0; // (xr2[1], xs2[1])
    array_1d<double,2> frq;
    array_1d<double,2> fsq;

    //array_1d<Vector,4> funct_q;                        
    //array_1d<Matrix,4> deriv_q;

    Matrix N_ans = ZeroMatrix(n_ans_points, number_of_nodes);   // N_ans(p,i) = N_i at point p
    array_1d<Matrix,4> DN_ans;                                  // DN_ansp=dNi/dxi (DN_ans[p](i, 0)), DN_ansp=dNi/deta (DN_ans[p](i, 1))
    array_1d<array_1d<Vector,3>,4> akovr_ans;
    Vector Np(number_of_nodes);

    for (SizeType p = 0; p < n_ans_points; ++p) {
        //funct_q[p].resize(number_of_nodes);
        DN_ans[p].resize(number_of_nodes, 2, false);

        array_1d<double,3> local_coords;
        local_coords[0] = ans_points[p][0]; // xi
        local_coords[1] = ans_points[p][1]; // eta
        local_coords[2] = 0.0;

        r_geom.ShapeFunctionsValues(Np, local_coords);
        row(N_ans, p) = Np;

        r_geom.ShapeFunctionsLocalGradients(DN_ans[p], local_coords);

        CovariantBaseVectorsMidsurface(akovr_ans[p], DN_ans[p], row(N_ans, p), ConfigurationType::Reference, thickness);

    }

    ////////////////////////////////////////////////////////////////END ANS STUFF////////////////////////////////////////////////////////////////
    //array_1d<Matrix,2> akovr_1q;
    //array_1d<Matrix,2> akonr_1q;
    //array_1d<Matrix,2> amkovr_1q;
    //array_1d<Matrix,2> amkonr_1q;
    //array_1d<Matrix,2> a3kvpr_1q;
//
    //array_1d<Matrix,2> akovr_2q;
    //array_1d<Matrix,2> akonr_2q;
    //array_1d<Matrix,2> amkovr_2q;
    //array_1d<Matrix,2> amkonr_2q;
    //array_1d<Matrix,2> a3kvpr_2q;

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
        // getting information for integration
        const double integration_weight_i = r_integration_points[point_number].Weight();
        const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
        const Vector& Nshape = row(Ncontainer,point_number);        // Node shape function values at the current integration point. Nshape[i] = N_i evaluated at the current GP

        // CovariantBaseVectors(current_covariant_base_vectors,shape_functions_gradients_i,ConfigurationType::Current,thickness);
        CovariantBaseVectorsMidsurface(akovr,shape_functions_gradients_i,Nshape,ConfigurationType::Reference,thickness);
        DirectorDerivatives(a3kvp,shape_functions_gradients_i,thickness);

        // CovariantMetric(covariant_metric_current,current_covariant_base_vectors);
        CovariantMetric(amkovr,akovr);
        ContravariantMetric(amkonr,amkovr,amdet_body);

        ContraVariantBaseVectors(akonr,amkonr,akovr);

        JacobiDeterminante(detJ_surface,akovr);

        ////////////////////////////////////////////////////////////////BEGIN ANS STUFF////////////////////////////////////////////////////////////////
        const auto& r_gp = r_integration_points[point_number];

        const double r  = r_gp.X();
        const double s = r_gp.Y();

        s8_ansqshapefunctions(frq, fsq, r, s);

        ////////////////////////////////////////////////////////////////END ANS STUFF////////////////////////////////////////////////////////////////

        BoundedMatrix<double, 12, 12> Dmatrix=ZeroMatrix(12,12);
        Matrix Bop = ZeroMatrix(12,number_dofs);
        CalculatelinearBOperator(Bop,akovr,a3kvp,shape_functions_gradients_i,Nshape,number_of_nodes);
        //-------------------------------------- modifications due to ans 
        BOperatorANSmodification(Bop,frq,fsq,akovr_ans,a3kvp,DN_ans,N_ans,number_of_nodes);


        //-------------------------------------- loop over GP in thickness direction for preintegration of constitutive law
        for (SizeType k=0; k<2; ++k){           // separate function PreintegrateThroughThicknessConstitutive() ?
            double Theta3 = gpcoord3[k];
            double tweight = gpweight[k];
            CovariantBaseVectorsShellBody(gkovr,shape_functions_gradients_i,Nshape,ConfigurationType::Reference,Theta3,thickness);
            CovariantMetric(gmkovr,gkovr);
            ContravariantMetric(gmkonr,gmkovr,gmdet_body);
            //ContraVariantBaseVectors(gkonr,gmkonr,gkovr);

            double scalefactor= std::sqrt(gmdet_body)/detJ_surface * tweight;

            CalculateMaterialLaw(Dmatrix,gmkonr,thickness,ConstitutiveLawType::gStVenantKirchhoff, Theta3, scalefactor);
        }
        Dmatrix(2,2) *= 5.0/6.0; 
        Dmatrix(2,4) *= 5.0/6.0; 
        Dmatrix(4,2) *= 5.0/6.0;        // separate function ApplyShearCorrections()?
        Dmatrix(4,4) *= 5.0/6.0;
        Dmatrix(8,8) *= 0.7;     
        Dmatrix(8,10) *= 0.7;
        Dmatrix(10,8) *= 0.7;
        Dmatrix(10,10) *= 0.7;

        double weight = integration_weight_i * detJ_surface; // * thickness*0.5; 
        Matrix DB = ZeroMatrix(12,number_dofs); 
        noalias(DB) = prod(Dmatrix, Bop);
        rLeftHandSideMatrix += prod(trans(Bop), DB) * weight;
    }
}

void Shell7pElement::CovariantBaseVectorsMidsurface(array_1d<Vector,3>& akovr,
     const Matrix& rShapeFunctionGradientValues, const Vector& rNshape, const ConfigurationType& rConfiguration, const double& thickness) const
{
    // pass/call this ShapeFunctionsLocalGradients[pnt]
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    Vector a1 = ZeroVector(dimension);
    Vector a2 = ZeroVector(dimension);
    Vector a3 = ZeroVector(dimension);
    
    // Vector current_displacement = ZeroVector(dimension*number_of_nodes);
    //if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);


    for (SizeType i=0;i<number_of_nodes;++i){

        const Vector& nodal_normal = r_geom[i].GetValue(NORMAL);        // node.GetValue(NORMAL)

        a1[0] += r_geom.GetPoint( i ).X0() * rShapeFunctionGradientValues(i, 0);
        a1[1] += r_geom.GetPoint( i ).Y0() * rShapeFunctionGradientValues(i, 0);
        a1[2] += r_geom.GetPoint( i ).Z0() * rShapeFunctionGradientValues(i, 0);

        a2[0] += r_geom.GetPoint( i ).X0() * rShapeFunctionGradientValues(i, 1);
        a2[1] += r_geom.GetPoint( i ).Y0() * rShapeFunctionGradientValues(i, 1);
        a2[2] += r_geom.GetPoint( i ).Z0() * rShapeFunctionGradientValues(i, 1);

        a3[0] += nodal_normal[0] * rNshape[i];
        a3[1] += nodal_normal[1] * rNshape[i];
        a3[2] += nodal_normal[2] * rNshape[i];
    }
    akovr[0] = a1;
    akovr[1] = a2;
    akovr[2] = a3*thickness*0.5;
}

void Shell7pElement::DirectorDerivatives(array_1d<Vector,2>& a3kvp,const Matrix& rShapeFunctionGradientValues, const double& thickness) const

{
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    Vector a31 = ZeroVector(dimension);
    Vector a32 = ZeroVector(dimension);

    for (SizeType i=0;i<number_of_nodes;++i){

        const Vector& nodal_normal = r_geom[i].GetValue(NORMAL);
        a31[0] += nodal_normal[0]*rShapeFunctionGradientValues(i, 0);
        a31[1] += nodal_normal[1]*rShapeFunctionGradientValues(i, 0);
        a31[2] += nodal_normal[2]*rShapeFunctionGradientValues(i, 0);

        a32[0] += nodal_normal[0]*rShapeFunctionGradientValues(i, 1);
        a32[1] += nodal_normal[1]*rShapeFunctionGradientValues(i, 1);
        a32[2] += nodal_normal[2]*rShapeFunctionGradientValues(i, 1);
    }
    a3kvp[0] = a31*thickness*0.5;
    a3kvp[1] = a32*thickness*0.5;
}
void Shell7pElement::CovariantMetric(Matrix& rMetric,const array_1d<Vector,3>& rBaseVectorCovariant) const
{
    rMetric = ZeroMatrix(3);
    for (SizeType i=0;i<3;++i){
        for (SizeType j=0;j<3;++j){
            rMetric(i,j) = inner_prod(rBaseVectorCovariant[i],rBaseVectorCovariant[j]);
        }
    }
}
 
void Shell7pElement::ContravariantMetric(Matrix& rMetric, const Matrix& rCovariantMetric, double& detJ_body) const
{
    rMetric = ZeroMatrix(3);
    detJ_body = 0.0;
    MathUtils<double>::InvertMatrix3(rCovariantMetric, rMetric, detJ_body);         // 1.Uses the general 3×3 inversion 2.Checks determinant 3.Works even if orthogonality assumption fails   
    const double det_tol = 1.0e-14;

    KRATOS_ERROR_IF(detJ_body <= -det_tol) << "Negative covariant metric determinant detected. det = " << detJ_body << " in element " << Id() << std::endl;
    // check the determinant for singularity (near-zero = bad matrix condition) meaning large contravariant metric coeefs
    KRATOS_ERROR_IF(std::abs(detJ_body) < det_tol) << "Singular covariant metric detected. det = " << detJ_body << " in element " << Id() << std::endl;  
}

void Shell7pElement::ContraVariantBaseVectors(array_1d<Vector,3>& rBaseVectors,const Matrix& rContraVariantMetric,
    const array_1d<Vector,3> rCovariantBaseVectors) const
{
    rBaseVectors[0] = ZeroVector(3);
    rBaseVectors[1] = ZeroVector(3);
    rBaseVectors[2] = ZeroVector(3);

    rBaseVectors[0] = rContraVariantMetric(0,0)*rCovariantBaseVectors[0] + rContraVariantMetric(0,1)*rCovariantBaseVectors[1] + rContraVariantMetric(0,2)*rCovariantBaseVectors[2];
    rBaseVectors[1] = rContraVariantMetric(1,0)*rCovariantBaseVectors[0] + rContraVariantMetric(1,1)*rCovariantBaseVectors[1] + rContraVariantMetric(1,2)*rCovariantBaseVectors[2];
    rBaseVectors[2] = rContraVariantMetric(2,0)*rCovariantBaseVectors[0] + rContraVariantMetric(2,1)*rCovariantBaseVectors[1] + rContraVariantMetric(2,2)*rCovariantBaseVectors[2]; 
}

void Shell7pElement::JacobiDeterminante(double& DetJ, const array_1d<Vector,3>& akovr) const
{
    array_1d<double, 3> a3;
    MathUtils<double>::CrossProduct(a3, akovr[0], akovr[1]);
    DetJ = MathUtils<double>::Norm(a3);
    KRATOS_ERROR_IF(DetJ<std::numeric_limits<double>::epsilon()) << "det of Jacobi smaller 0 for element with id" << Id() << std::endl;
}

void Shell7pElement::CovariantBaseVectorsShellBody(array_1d<Vector,3>& gkovr,
     const Matrix& rShapeFunctionGradientValues, const Vector& rNshape, const ConfigurationType& rConfiguration, const double& Theta3, const double& thickness) const
{
    // pass/call this ShapeFunctionsLocalGradients[pnt]
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    Vector g1 = ZeroVector(dimension);
    Vector g2 = ZeroVector(dimension);
    Vector g3 = ZeroVector(dimension);
    
    // Vector current_displacement = ZeroVector(dimension*number_of_nodes);
    //if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);


    for (SizeType i=0;i<number_of_nodes;++i){

        const Vector& nodal_normal = r_geom[i].GetValue(NORMAL) * thickness*0.5;        // node.GetValue(NORMAL)

        g1[0] += (r_geom.GetPoint( i ).X0() + Theta3 * nodal_normal[0]) * rShapeFunctionGradientValues(i, 0);
        g1[1] += (r_geom.GetPoint( i ).Y0() + Theta3 * nodal_normal[1]) * rShapeFunctionGradientValues(i, 0);
        g1[2] += (r_geom.GetPoint( i ).Z0() + Theta3 * nodal_normal[2]) * rShapeFunctionGradientValues(i, 0);

        g2[0] += (r_geom.GetPoint( i ).X0() + Theta3 * nodal_normal[0]) * rShapeFunctionGradientValues(i, 1);
        g2[1] += (r_geom.GetPoint( i ).Y0() + Theta3 * nodal_normal[1]) * rShapeFunctionGradientValues(i, 1);
        g2[2] += (r_geom.GetPoint( i ).Z0() + Theta3 * nodal_normal[2]) * rShapeFunctionGradientValues(i, 1);

        g3[0] += nodal_normal[0] * rNshape[i];
        g3[1] += nodal_normal[1] * rNshape[i];
        g3[2] += nodal_normal[2] * rNshape[i];
    }
    gkovr[0] = g1;
    gkovr[1] = g2;
    gkovr[2] = g3;
}

void Shell7pElement::CalculateMaterialLaw(BoundedMatrix<double, 12, 12>& CL, const Matrix& gmkonr, const double& thickness,
const ConstitutiveLawType& option, const double& Theta3, const double& fact) const
{
    const auto& r_properties = GetProperties();
    const double E = r_properties[YOUNG_MODULUS];
    const double nu = r_properties[POISSON_RATIO];
    const double G = E/(2.0*(1.0+nu));
    const double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    const double mu = G;

    if (option == ConstitutiveLawType::gStVenantKirchhoff) {
        double C[3][3][3][3] = {};
        //double Theta3[2] = {1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)};
        //double gpweight[2] = {1.0, 1.0};
        BoundedMatrix<double, 6, 6> CC = ZeroMatrix(6);

        for (SizeType i=0; i<3; ++i){
            for (SizeType j=0; j<3; ++j){
                for (SizeType k=0; k<3; ++k){
                    for (SizeType l=0; l<3; ++l){
                        C[i][j][k][l] = lambda*gmkonr(i,j)*gmkonr(k,l) + mu*( gmkonr(i,k)*gmkonr(j,l) + gmkonr(i,l)*gmkonr(k,j) );
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
        CC(2,2) = C[2][0][2][0]; //*5.0/6.0;     // sigma 13 with shear correction factor alpha=5/6 for E13 and E23
        CC(2,3) = C[2][0][1][1];
        CC(2,4) = C[2][0][2][1]; //*5.0/6.0;
        CC(2,5) = C[2][0][2][2];

        CC(3,0) = C[1][1][0][0];
        CC(3,1) = C[1][1][1][0];
        CC(3,2) = C[1][1][2][0];            // sigma 22
        CC(3,3) = C[1][1][1][1];
        CC(3,4) = C[1][1][2][1];
        CC(3,5) = C[1][1][2][2];

        CC(4,0) = C[2][1][0][0];
        CC(4,1) = C[2][1][1][0];
        CC(4,2) = C[2][1][2][0]; //*5.0/6.0;    // sigma 23 with shear correction factor alpha=5/6 for E13 and E23
        CC(4,3) = C[2][1][1][1];
        CC(4,4) = C[2][1][2][1]; //*5.0/6.0;
        CC(4,5) = C[2][1][2][2];

        CC(5,0) = C[2][2][0][0];
        CC(5,1) = C[2][2][1][0];
        CC(5,2) = C[2][2][2][0];            // sigma 33
        CC(5,3) = C[2][2][1][1];
        CC(5,4) = C[2][2][2][1];
        CC(5,5) = C[2][2][2][2];

        
        for (SizeType i=0; i<6; ++i){
            const SizeType i6 = i + 6;
            for (SizeType j=0; j<6; ++j){
                const SizeType j6 = j + 6;
                CL(i,j) += CC(i,j)*fact;
                CL(i6,j) += CC(i,j)*Theta3*fact;
                CL(j,i6) += CC(j,i)*Theta3*fact; 
                CL(i6,j6) += CC(i,j)*Theta3*Theta3*fact;
             }
        }

       //CL(2,2) *= 5.0/6.0;    // shear correction factor alpha=5/6 for n13,n23
       //CL(2,4) *= 5.0/6.0;    // schould i use thickness_q=alpha*thickness instead of just thickness for shear terms?
       //CL(4,2) *= 5.0/6.0;
       //CL(4,4) *= 5.0/6.0;
       //CL(8,8) *= 0.7;        // shear correction factor betta=0.7 for m13,m23
       //CL(8,10) *= 0.7;
       //CL(10,8) *= 0.7;
       //CL(10,10) *= 0.7;
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
                                                                                                                                                                  // [N_node][derivative direction]   // Nshape[i] = N_i evaluated at the current GP
void Shell7pElement::CalculatelinearBOperator(Matrix& bop, const array_1d<Vector,3>& CovariantBaseVectors, const array_1d<Vector,2>& DirectorDerivatives, const Matrix& ShapeFunctionGradientValues, const Vector& Nshape, const SizeType& number_of_nodes) const
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
const double a31x=DirectorDerivatives[0][0];
const double a31y=DirectorDerivatives[0][1];
const double a31z=DirectorDerivatives[0][2];
const double a32x=DirectorDerivatives[1][0];
const double a32y=DirectorDerivatives[1][1];
const double a32z=DirectorDerivatives[1][2];

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

            bop(6,index)   = dNd1*a31x;
            bop(6,index+1) = dNd1*a31y;
            bop(6,index+2) = dNd1*a31z;                // betta 11
            bop(6,index+3) = dNd1*a1x;
            bop(6,index+4) = dNd1*a1y;
            bop(6,index+5) = dNd1*a1z;

            bop(7,index)   = dNd2*a31x + dNd1*a32x;
            bop(7,index+1) = dNd2*a31y + dNd1*a32y;
            bop(7,index+2) = dNd2*a31z + dNd1*a32z;    // betta 12
            bop(7,index+3) = dNd2*a1x + dNd1*a2x;
            bop(7,index+4) = dNd2*a1y + dNd1*a2y;
            bop(7,index+5) = dNd2*a1z + dNd1*a2z;

            bop(8,index)   = 0.0;
            bop(8,index+1) = 0.0;
            bop(8,index+2) = 0.0;                       // betta 13
            bop(8,index+3) = N*a31x + dNd1*a3x;
            bop(8,index+4) = N*a31y + dNd1*a3y;
            bop(8,index+5) = N*a31z + dNd1*a3z;

            bop(9,index)   = dNd2*a32x;
            bop(9,index+1) = dNd2*a32y;
            bop(9,index+2) = dNd2*a32z;                 // betta 22
            bop(9,index+3) = dNd2*a2x;
            bop(9,index+4) = dNd2*a2y;
            bop(9,index+5) = dNd2*a2z;

            bop(10,index)   = 0.0;
            bop(10,index+1) = 0.0;
            bop(10,index+2) = 0.0;                       // betta 23
            bop(10,index+3) = N*a32x + dNd2*a3x;
            bop(10,index+4) = N*a32y + dNd2*a3y;
            bop(10,index+5) = N*a32z + dNd2*a3z;

            bop(11,index)   = 0.0;
            bop(11,index+1) = 0.0;
            bop(11,index+2) = 0.0;                       // betta 33
            bop(11,index+3) = 0.0;
            bop(11,index+4) = 0.0;
            bop(11,index+5) = 0.0;
    }
}

void Shell7pElement::s8_ansqshapefunctions(array_1d<double,2>& frq, array_1d<double,2>& fsq,const double r, const double s) const
{

   frq[0] = 0.5 * (1.0 - s);
   frq[1] = 0.5 * (1.0 + s);
   fsq[0] = 0.5 * (1.0 - r);
   fsq[1] = 0.5 * (1.0 + r);
}

void Shell7pElement::BOperatorANSmodification(Matrix& Bop, const array_1d<double,2>& frq, const array_1d<double,2>& fsq,
    const array_1d<array_1d<Vector,3>,4>& akovr_ans,const array_1d<Vector,2>& a3kvp,const array_1d<Matrix,4>& DN_ans,
    const Matrix& N_ans, const SizeType& number_of_nodes) const
{
    for (SizeType inode = 0; inode < number_of_nodes; ++inode)
  {
    const SizeType node_start = inode*6;

    Bop(2,node_start+0)= 0.0;
    Bop(2,node_start+1)= 0.0;
    Bop(2,node_start+2)= 0.0;
    Bop(2,node_start+3)= 0.0;
    Bop(2,node_start+4)= 0.0;
    Bop(2,node_start+5)= 0.0;

    Bop(4,node_start+0)= 0.0;
    Bop(4,node_start+1)= 0.0;
    Bop(4,node_start+2)= 0.0;
    Bop(4,node_start+3)= 0.0;
    Bop(4,node_start+4)= 0.0;
    Bop(4,node_start+5)= 0.0;

    for (SizeType isamp = 0; isamp < 2; ++isamp)
    {
      const double a1x1 = akovr_ans[isamp][0][0];     // a1x
      const double a1y1 = akovr_ans[isamp][0][1];     // a1y
      const double a1z1 = akovr_ans[isamp][0][2];     // a1z
      const double a3x1 = akovr_ans[isamp][2][0];     // a3x
      const double a3y1 = akovr_ans[isamp][2][1];     // a3y
      const double a3z1 = akovr_ans[isamp][2][2];     // a3z

      const double a2x2 = akovr_ans[isamp+2][1][0];
      const double a2y2 = akovr_ans[isamp+2][1][1];
      const double a2z2 = akovr_ans[isamp+2][1][2];
      const double a3x2 = akovr_ans[isamp+2][2][0];
      const double a3y2 = akovr_ans[isamp+2][2][1];
      const double a3z2 = akovr_ans[isamp+2][2][2];

      //const double a31x1 = a3kvpc1q[isamp][0][0];
      //const double a31y1 = a3kvpc1q[isamp][1][0];
      //const double a31z1 = a3kvpc1q[isamp][2][0];

      //const double a32x2 = a3kvpc2q[isamp][0][1];
      //const double a32y2 = a3kvpc2q[isamp][1][1];
      //const double a32z2 = a3kvpc2q[isamp][2][1];

      const double N1 = N_ans(isamp, inode);
      const double N2 = N_ans(isamp+2, inode);

      const double dNd1 = DN_ans[isamp](inode, 0);
      const double dNd2 = DN_ans[isamp+2](inode, 1);

      const double N1_ans = frq[isamp];
      const double N2_ans = fsq[isamp];
/*--------------------------------------------------E13(CONST)-------- */
      Bop(2,node_start+0)+= dNd1*a3x1*N1_ans;
      Bop(2,node_start+1)+= dNd1*a3y1*N1_ans;
      Bop(2,node_start+2)+= dNd1*a3z1*N1_ans;
      Bop(2,node_start+3)+= N1*a1x1*N1_ans;
      Bop(2,node_start+4)+= N1*a1y1*N1_ans;
      Bop(2,node_start+5)+= N1*a1z1*N1_ans;
/*-------------------------------------------------E23(CONST)-------- */
      Bop(4,node_start+0)+= dNd2*a3x2*N2_ans;
      Bop(4,node_start+1)+= dNd2*a3y2*N2_ans;
      Bop(4,node_start+2)+= dNd2*a3z2*N2_ans;
      Bop(4,node_start+3)+= N2*a2x2*N2_ans;  
      Bop(4,node_start+4)+= N2*a2y2*N2_ans;
      Bop(4,node_start+5)+= N2*a2z2*N2_ans;
    }
  }


}

}   // namespace Kratos