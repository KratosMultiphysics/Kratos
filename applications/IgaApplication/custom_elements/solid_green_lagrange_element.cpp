//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:            BSD License
//                      Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:       Thomas Lillieskold
//  Co-author:          Andrea Gorgi
//
#include "solid_green_lagrange_element.h"

#include <algorithm>
#include <cmath>

#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "iga_application_variables.h"

namespace Kratos
{

//--------------------------------------------[1. Creación del elemento]--------------------------------------------

SolidGreenLagrangeElement::SolidGreenLagrangeElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

SolidGreenLagrangeElement::SolidGreenLagrangeElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer SolidGreenLagrangeElement::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

Element::Pointer SolidGreenLagrangeElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SolidGreenLagrangeElement>(NewId, pGeometry, pProperties);
}

//--------------------------------------------[2. Inicialización, cuadratura y grados de libertad]--------------------------------------------

void SolidGreenLagrangeElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    CheckGeometryDimension();

    // Calculate the integration weight on the reference configuration.
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    KRATOS_ERROR_IF(r_integration_points.size() != 1)
        << Info() << " requires exactly one integration point." << std::endl;

    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    KRATOS_ERROR_IF(r_DN_De.size() != 1)
        << Info() << " requires shape function derivatives at its integration point." << std::endl;

    KRATOS_ERROR_IF(r_DN_De[0].size1() != r_geometry.PointsNumber() ||
                    r_DN_De[0].size2() != r_geometry.LocalSpaceDimension())
        << Info() << " has inconsistent local shape function derivatives." << std::endl;

    const IndexType mDim = r_geometry.LocalSpaceDimension();
    Matrix InvJ0(mDim, mDim);
    GeometryType::JacobiansType J0;

    // Total Lagrangian integration requires the initial geometry.
    // Subtract the mesh motion to recover it when computing J0.
    Matrix delta_position(r_geometry.PointsNumber(), r_geometry.WorkingSpaceDimension());
    for (IndexType a = 0; a < delta_position.size1(); ++a) {
        const auto& r_node = r_geometry[a];
        for (IndexType i = 0; i < delta_position.size2(); ++i) {
            delta_position(a, i) = r_node.Coordinates()[i] - r_node.GetInitialPosition()[i];
        }
    }
    // Compute the Jacobians at the quadrature points using the initial coordinates
    // (current coordinates minus delta_position) and store them in J0.
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    // A validated 2D geometry embedded in 3D has a constant reference Z.
    // Its 3x2 Jacobian has a zero Z row, so the XY block is sufficient.
    const Matrix Jacobian = r_geometry.WorkingSpaceDimension() == r_geometry.LocalSpaceDimension()
        ? J0[0]
        : Matrix(project(J0[0], range(0, mDim), range(0, mDim)));

    double DetJ0 = MathUtils<double>::Det2(Jacobian);
    KRATOS_ERROR_IF_NOT(std::isfinite(DetJ0) && DetJ0 > 0.0)
        << Info() << " has invalid reference Jacobian determinant " << DetJ0 << std::endl;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix2(Jacobian, InvJ0, DetJ0);
    mDN_DX = prod(r_DN_De[0], InvJ0);


    const double Thickness = GetProperties()[THICKNESS];
    const double Weight = r_integration_points[0].Weight() * DetJ0 * Thickness;
    KRATOS_ERROR_IF_NOT(std::isfinite(Weight) && Weight > 0.0)
        << Info() << " has invalid reference integration weight " << Weight << std::endl;
    SetValue(INTEGRATION_WEIGHT, Weight);

    InitializeMaterial();

    KRATOS_CATCH("")
}

void SolidGreenLagrangeElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = GetGeometry();
    rResult.resize(2 * r_geometry.PointsNumber());
    for (IndexType a = 0; a < r_geometry.PointsNumber(); ++a) {
        rResult[2 * a] = r_geometry[a].GetDof(DISPLACEMENT_X).EquationId();
        rResult[2 * a + 1] = r_geometry[a].GetDof(DISPLACEMENT_Y).EquationId();
    }
}

void SolidGreenLagrangeElement::GetDofList(
    DofsVectorType& rDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = GetGeometry();
    rDofList.resize(2 * r_geometry.PointsNumber());

    for (IndexType a = 0; a < r_geometry.PointsNumber(); ++a) {
        rDofList[2 * a] = r_geometry[a].pGetDof(DISPLACEMENT_X);
        rDofList[2 * a + 1] = r_geometry[a].pGetDof(DISPLACEMENT_Y);
    }
}

Element::IntegrationMethod SolidGreenLagrangeElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

int SolidGreenLagrangeElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    CheckGeometryDimension();
    Element::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(pGetProperties()) << Info() << " requires properties." << std::endl;
    const auto& r_properties = GetProperties();

    KRATOS_ERROR_IF_NOT(r_properties.Has(THICKNESS))
        << Info() << " requires THICKNESS." << std::endl;

    const double Thickness = r_properties[THICKNESS];

    KRATOS_ERROR_IF_NOT(std::isfinite(Thickness) && Thickness > 0.0)
        << Info() << " requires a finite, positive THICKNESS." << std::endl;

    for (const auto& r_node : GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
    }

    KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
        << Info() << " requires CONSTITUTIVE_LAW." << std::endl;

    const auto& p_constitutivelaw = r_properties[CONSTITUTIVE_LAW];
    KRATOS_ERROR_IF_NOT(p_constitutivelaw)
        << Info() << " requires a valid constitutive law." << std::endl;

    KRATOS_ERROR_IF(
        p_constitutivelaw->WorkingSpaceDimension() != 2 ||
        p_constitutivelaw->GetStrainSize() != 3 ||
        p_constitutivelaw->GetStrainMeasure() != ConstitutiveLaw::StrainMeasure_GreenLagrange ||
        p_constitutivelaw->GetStressMeasure() != ConstitutiveLaw::StressMeasure_PK2)
        << Info() << " requires a 2D law with three components, "
        << "Green-Lagrange strain and PK2 stress." << std::endl;

    return p_constitutivelaw->Check(r_properties, GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void SolidGreenLagrangeElement::CheckGeometryDimension() const
{
    const auto& r_geometry = GetGeometry();
    //LocalSpaceDimension() is the dimension of the parametric space of the geometry,
    const SizeType local_dimension = r_geometry.LocalSpaceDimension();

    //WorkingSpaceDimension() is the dimension of the space in which the geometry is embedded.
    const SizeType working_dimension = r_geometry.WorkingSpaceDimension();

    KRATOS_ERROR_IF(local_dimension == 3)
        << Info() << " does not yet implement the 3D solid formulation." << std::endl;
    KRATOS_ERROR_IF(local_dimension != 2 || (working_dimension != 2 && working_dimension != 3))
        << Info() << " requires a two-dimensional or three-dimensional solid geometry."
        << std::endl;
    KRATOS_ERROR_IF(r_geometry.PointsNumber() == 0)
        << Info() << " requires control points." << std::endl;

    
    //The reference geometry must be parallel to the XY plane for 2D elements embedded in 3D.
    if (working_dimension == 3) {
        const auto& r_origin = r_geometry[0].GetInitialPosition();
        double xy_extent = 0.0;
        double z_extent = 0.0;

        for (const auto& r_node : r_geometry) {
            const auto& r_position = r_node.GetInitialPosition();

            xy_extent = std::max({xy_extent,std::abs(r_position[0] - r_origin[0]),
                                            std::abs(r_position[1] - r_origin[1])
                                }); 
            z_extent = std::max(z_extent, std::abs(r_position[2] - r_origin[2]));
        }
        constexpr double relative_tolerance = 1.0e-12;

        //Checks if reference geometry is parallel to XY plane
        KRATOS_ERROR_IF(z_extent > relative_tolerance * xy_extent)
            << Info() << " requires the 2D reference geometry to be parallel to the XY plane."
            << std::endl;
    }
}


std::string SolidGreenLagrangeElement::Info() const
{
    return "SolidGreenLagrangeElement #" + std::to_string(Id());
}

//--------------------------------------------[3. Cinemática de Green–Lagrange]--------------------------------------------

SolidGreenLagrangeElement::KinematicVariables SolidGreenLagrangeElement::CalculateKinematics() const
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    KRATOS_ERROR_IF(mDN_DX.size1() != r_geometry.PointsNumber() || mDN_DX.size2() != 2)
        << Info() << " must be initialized before calculating kinematics." << std::endl;

    KinematicVariables kinematics;

    // The deformation gradient F is computed as F_ij = delta_ij + sum_a(u_ai * dN_a/dX_j)
    for (IndexType a = 0; a < r_geometry.PointsNumber(); ++a) {
        const auto& r_displacement = r_geometry[a].FastGetSolutionStepValue(DISPLACEMENT); //displacement vector of node a
        for (IndexType i = 0; i < 2; ++i) {
            for (IndexType j = 0; j < 2; ++j) {
                kinematics.F(i, j) += r_displacement[i] * mDN_DX(a, j); 
        }
    }

    kinematics.DeterminantF = MathUtils<double>::Det2(kinematics.F);
    KRATOS_ERROR_IF_NOT(std::isfinite(kinematics.DeterminantF) && kinematics.DeterminantF > 0.0)
        << Info() << " has invalid deformation gradient determinant "
        << kinematics.DeterminantF << std::endl;

    // Cauchy-Green deformation tensor: C = F^T * F
    const Matrix right_cauchy_green = prod(trans(kinematics.F), kinematics.F);

    
    kinematics.Strain[0] = 0.5 * (right_cauchy_green(0, 0) - 1.0);  // Green-Lagrange strain: Exx= 0.5 * (C - I)
    kinematics.Strain[1] = 0.5 * (right_cauchy_green(1, 1) - 1.0);  // Green-Lagrange strain: Eyy= 0.5 * (C - I)
    kinematics.Strain[2] = right_cauchy_green(0, 1);                //Green-Lagrange Engineering strain: 2* Exy = Cxy

    return kinematics;

    KRATOS_CATCH("")
}

//--------------------------------------------[4. Respuesta constitutiva (pendiente)]--------------------------------------------

void SolidGreenLagrangeElement::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
}

void SolidGreenLagrangeElement::ApplyConstitutiveLaw(
    KinematicVariables& rKinematics,
    ConstitutiveLaw::Parameters& rValues,
    ConstitutiveVariables& rConstitutiveVariables)

{
    KRATOS_TRY

    auto& r_options = rValues.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    rValues.SetDeformationGradientF(rKinematics.F);
    rValues.SetDeterminantF(rKinematics.DeterminantF);
    rValues.SetStrainVector(rKinematics.Strain);

    rValues.SetStressVector(rConstitutiveVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiveVariables.D);

    mpConstitutiveLaw->CalculateMaterialResponse(
        rValues, ConstitutiveLaw::StressMeasure_PK2);

    KRATOS_CATCH("")
}


//--------------------------------------------[5. Ensamblado de volumen]--------------------------------------------

void SolidGreenLagrangeElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void SolidGreenLagrangeElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

    const SizeType mat_size = 2 * number_of_control_points;
    const double int_to_reference_weight = GetValue(INTEGRATION_WEIGHT);
    rLeftHandSideMatrix = ZeroMatrix(mat_size, mat_size);

    const Matrix& DN_DX = mDN_DX;

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();

    auto kinematics = CalculateKinematics();

    // The deformation gradient F is used to compute the B matrix, which relates strains to displacements.
    const Matrix& r_F = kinematics.F;

    // B maps nodal displacement variations to Green-Lagrange strain variations.
    Matrix B = ZeroMatrix(strain_size, mat_size);

    for (IndexType a = 0; a < number_of_control_points; ++a) {
        const double dN_dX = DN_DX(a, 0);
        const double dN_dY = DN_DX(a, 1);

        for (IndexType i = 0; i < 2; ++i) {
            const IndexType column = 2 * a + i;
            B(0, column) = r_F(i, 0) * dN_dX;                       // B11 = Fx * dN/dX
            B(1, column) = r_F(i, 1) * dN_dY;                       // B22 = Fy * dN/dY
            B(2, column) = r_F(i, 0) * dN_dY + r_F(i, 1) * dN_dX;   // B12 = Fx * dN/dY + Fy * dN/dX
        }
    }

    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    ConstitutiveVariables this_constitutive_variables(strain_size);

    ApplyConstitutiveLaw(kinematics, Values, this_constitutive_variables);

    const Matrix DB = prod(this_constitutive_variables.D, B);

    // The local material stiffness matrix is computed as K = B^T * D * B
    noalias(rLeftHandSideMatrix) += int_to_reference_weight * prod(trans(B), DB);

    const Matrix weighted_stress = int_to_reference_weight * MathUtils<double>::StressVectorToTensor(this_constitutive_variables.StressVector);

    //Local geometric stiffness matrix
    Matrix reduced_Kg(number_of_control_points, number_of_control_points);

    // Reduced geometric stiffness is w0 * DN_DX * S * trans(DN_DX).
    MathUtils<double>::BDBtProductOperation(reduced_Kg, weighted_stress, DN_DX);

    //Expand the reduced geometric stiffness matrix to the full size and add it to the left-hand side matrix
    MathUtils<double>::ExpandAndAddReducedMatrix(rLeftHandSideMatrix, reduced_Kg, 2);

    KRATOS_CATCH("")
}

void SolidGreenLagrangeElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const auto& r_N = r_geometry.ShapeFunctionsValues(GetIntegrationMethod());
    const auto& r_body_force = GetValue(BODY_FORCE);
    const double weight = GetValue(INTEGRATION_WEIGHT);

    auto kinematics = CalculateKinematics();
    ConstitutiveVariables constitutive_variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ApplyConstitutiveLaw(kinematics, values, constitutive_variables);

    // S is the second Piola-Kirchhoff stress tensor
    const Matrix S = MathUtils<double>::StressVectorToTensor(constitutive_variables.StressVector); 

    //P is the first Piola-Kirchhoff stress tensor -> P = F * S
    const Matrix P = prod(kinematics.F, S);  

    //Internal Forces = B^T * P
    const Matrix internal_forces = prod(mDN_DX, trans(P));

    rRightHandSideVector = ZeroVector(2 * r_geometry.size());
    for (IndexType a = 0; a < r_geometry.size(); ++a) {
        for (IndexType i = 0; i < 2; ++i) {
            rRightHandSideVector[2 * a + i] = weight * (r_N(0, a) * r_body_force[i] - internal_forces(a, i)); // RHS = weight * (N_a * b_i - B^T * P)
        }
    }

    KRATOS_CATCH("")
}

//--------------------------------------------[6. Ciclo de carga y resultados]--------------------------------------------

void SolidGreenLagrangeElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    auto kinematics = CalculateKinematics();
    ConstitutiveVariables constitutive_variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    ApplyConstitutiveLaw(kinematics, values, constitutive_variables);

    const Matrix S = MathUtils<double>::StressVectorToTensor(constitutive_variables.StressVector);
    const Matrix P = prod(kinematics.F, S);
    const Matrix sigma = prod(P, trans(kinematics.F)) / kinematics.DeterminantF;

    SetValue(CAUCHY_STRESS_XX, sigma(0, 0));
    SetValue(CAUCHY_STRESS_YY, sigma(1, 1));
    SetValue(CAUCHY_STRESS_XY, sigma(0, 1));

    KRATOS_CATCH("")
}

} // namespace Kratos
