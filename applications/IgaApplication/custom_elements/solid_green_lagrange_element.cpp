//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Thomas Lillieskold, Andrea Gorgi
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "custom_elements/solid_green_lagrange_element.h"

namespace Kratos
{

SolidGreenLagrangeElement::SolidGreenLagrangeElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : SolidElement(
        NewId,
        pGeometry)
{
}

SolidGreenLagrangeElement::SolidGreenLagrangeElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : SolidElement(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer SolidGreenLagrangeElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SolidGreenLagrangeElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer SolidGreenLagrangeElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SolidGreenLagrangeElement>(NewId, pGeom, pProperties);
}

Matrix SolidGreenLagrangeElement::CalculateDeformationGradient(const Matrix& rDN_DX) const
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_control_points = r_geometry.size();

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3)
        << "SolidGreenLagrangeElement only supports a 2D or 3D deformation gradient, but the current dimension is "
        << mDim << std::endl;
    KRATOS_ERROR_IF(rDN_DX.size1() != number_of_control_points || rDN_DX.size2() != mDim)
        << "SolidGreenLagrangeElement received shape-function gradients with incompatible dimensions." << std::endl;

    Matrix deformation_gradient = IdentityMatrix(mDim);
    for (IndexType a = 0; a < number_of_control_points; ++a) {
        const auto& r_displacement = r_geometry[a].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType i = 0; i < mDim; ++i) {
            for (IndexType j = 0; j < mDim; ++j) {
                deformation_gradient(i, j) += r_displacement[i] * rDN_DX(a, j); 
            }
        }
    }

    return deformation_gradient;
}

void SolidGreenLagrangeElement::CalculateBMatrix(
    const Matrix& rDeformationGradient,
    const Matrix& rDN_DX,
    Matrix& rB) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType strain_size = mDim == 2 ? 3 : 6;
    const SizeType mat_size = number_of_control_points * mDim;

    KRATOS_ERROR_IF(rDeformationGradient.size1() != mDim || rDeformationGradient.size2() != mDim)
        << "SolidGreenLagrangeElement received a deformation gradient with incompatible dimensions." << std::endl;
    KRATOS_ERROR_IF(rDN_DX.size1() != number_of_control_points || rDN_DX.size2() != mDim)
        << "SolidGreenLagrangeElement received shape-function gradients with incompatible dimensions." << std::endl;

    if (rB.size1() != strain_size || rB.size2() != mat_size) {
        rB.resize(strain_size, mat_size, false);
    }
    noalias(rB) = ZeroMatrix(strain_size, mat_size);

    for (IndexType a = 0; a < number_of_control_points; ++a) {
        const SizeType first_column = a * mDim;
        for (IndexType i = 0; i < mDim; ++i) {
            if (mDim == 2) {
                rB(0, first_column + i) = rDeformationGradient(i, 0) * rDN_DX(a, 0);
                rB(1, first_column + i) = rDeformationGradient(i, 1) * rDN_DX(a, 1);
                rB(2, first_column + i) = rDeformationGradient(i, 0) * rDN_DX(a, 1)
                    + rDeformationGradient(i, 1) * rDN_DX(a, 0);
            } else {
                rB(0, first_column + i) = rDeformationGradient(i, 0) * rDN_DX(a, 0);
                rB(1, first_column + i) = rDeformationGradient(i, 1) * rDN_DX(a, 1);
                rB(2, first_column + i) = rDeformationGradient(i, 2) * rDN_DX(a, 2);
                rB(3, first_column + i) = rDeformationGradient(i, 0) * rDN_DX(a, 1)
                    + rDeformationGradient(i, 1) * rDN_DX(a, 0);
                rB(4, first_column + i) = rDeformationGradient(i, 1) * rDN_DX(a, 2)
                    + rDeformationGradient(i, 2) * rDN_DX(a, 1);
                rB(5, first_column + i) = rDeformationGradient(i, 0) * rDN_DX(a, 2)
                    + rDeformationGradient(i, 2) * rDN_DX(a, 0);
            }
        }
    }
}

// Deconstructor

SolidGreenLagrangeElement::~SolidGreenLagrangeElement()
{
}

void SolidGreenLagrangeElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();

    // calculate the integration weight
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    mDim = r_DN_De[0].size2();

    Matrix InvJ0(mDim,mDim);
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    double DetJ0;

    Matrix Jacobian;
    CalculateInitialJacobian(r_geometry, Jacobian);
    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    switch (mDim) {
        case 2:
        {
            const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
            const double int_to_reference_weight = r_integration_points[0].Weight() * std::abs(DetJ0) * thickness;
            SetValue(INTEGRATION_WEIGHT, int_to_reference_weight);

            break;
        }
        case 3:
        {
            const double int_to_reference_weight = r_integration_points[0].Weight() * std::abs(DetJ0);
            SetValue(INTEGRATION_WEIGHT, int_to_reference_weight);
            break;
        }
        default:
            KRATOS_ERROR << "Unsupported dimension for SolidGreenLagrangeElement: " << mDim << std::endl;
    }

    
    
}


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

void SolidGreenLagrangeElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_control_points = r_geometry.size();
    const SizeType mat_size = number_of_control_points * mDim;
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    const double int_to_reference_weight = GetValue(INTEGRATION_WEIGHT);

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3)
        << "SolidGreenLagrangeElement only supports 2D and 3D elements, but the current dimension is "
        << mDim << std::endl;
    KRATOS_ERROR_IF(strain_size != (mDim == 2 ? 3 : 6))
        << "SolidGreenLagrangeElement has an incompatible constitutive strain size." << std::endl;

    if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size) {
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(GetIntegrationMethod());
    Matrix jacobian;
    CalculateInitialJacobian(r_geometry, jacobian);
    Matrix inverse_jacobian(mDim, mDim);
    double det_jacobian = 0.0;
    MathUtils<double>::InvertMatrix(jacobian, inverse_jacobian, det_jacobian);
    const Matrix DN_DX = prod(r_DN_De[0], inverse_jacobian);

    const Matrix deformation_gradient = CalculateDeformationGradient(DN_DX);
    const double det_deformation_gradient = MathUtils<double>::Det(deformation_gradient);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_deformation_gradient) && det_deformation_gradient > 0.0)
        << "SolidGreenLagrangeElement has an invalid deformation gradient determinant: "
        << det_deformation_gradient << std::endl;

    ConstitutiveVariables constitutive_variables(strain_size);

    Matrix B(strain_size, mat_size);
    CalculateBMatrix(deformation_gradient, DN_DX, B);

    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ApplyConstitutiveLaw(
        deformation_gradient,
        det_deformation_gradient,
        values,
        constitutive_variables,
        true,
        ConstitutiveLaw::StressMeasure_PK2);

    CalculateAndAddKm(
        rLeftHandSideMatrix,
        B,
        values.GetConstitutiveMatrix(),
        int_to_reference_weight);
    CalculateAndAddKg(
        rLeftHandSideMatrix,
        DN_DX,
        values.GetStressVector(),
        int_to_reference_weight);

    KRATOS_CATCH("")
}

void SolidGreenLagrangeElement::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& rB,
    const Matrix& rConstitutiveMatrix,
    const double IntegrationWeight) const
{
    noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(
        trans(rB), Matrix(prod(rConstitutiveMatrix, rB)));
}

void SolidGreenLagrangeElement::CalculateAndAddKg(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& rDN_DX,
    const Vector& rStressVector,
    const double IntegrationWeight) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const Matrix weighted_stress = IntegrationWeight *
        MathUtils<double>::StressVectorToTensor(rStressVector);
    Matrix reduced_geometric_stiffness(number_of_control_points, number_of_control_points);
    MathUtils<double>::BDBtProductOperation(reduced_geometric_stiffness, weighted_stress, rDN_DX);
    MathUtils<double>::ExpandAndAddReducedMatrix(
        rLeftHandSideMatrix, reduced_geometric_stiffness, mDim);
}

void SolidGreenLagrangeElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_control_points = r_geometry.size();

    // reading integration points and local gradients
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const SizeType mat_size = number_of_control_points * mDim;
    const double int_to_reference_weight = GetValue(INTEGRATION_WEIGHT);

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3)
        << "SolidGreenLagrangeElement only supports 2D and 3D elements, but the current dimension is "
        << mDim << std::endl;
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix jacobian;
    CalculateInitialJacobian(r_geometry, jacobian);
    Matrix inverse_jacobian(mDim, mDim);
    double det_jacobian = 0.0;
    MathUtils<double>::InvertMatrix(jacobian, inverse_jacobian, det_jacobian);
    const Matrix DN_DX = prod(r_DN_De[0], inverse_jacobian);

    auto N = row(N_gausspoint,0); // these are the N which correspond to the gauss point "i_point"

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    KRATOS_ERROR_IF(strain_size != (mDim == 2 ? 3 : 6))
        << "SolidGreenLagrangeElement has an incompatible constitutive strain size." << std::endl;

    const Matrix deformation_gradient = CalculateDeformationGradient(DN_DX);
    const double det_deformation_gradient = MathUtils<double>::Det(deformation_gradient);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_deformation_gradient) && det_deformation_gradient > 0.0)
        << "SolidGreenLagrangeElement has an invalid deformation gradient determinant: "
        << det_deformation_gradient << std::endl;

    Matrix B(strain_size, mat_size);
    CalculateBMatrix(deformation_gradient, DN_DX, B);

    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);

    ConstitutiveVariables constitutive_variables(strain_size);
    ApplyConstitutiveLaw(
        deformation_gradient,
        det_deformation_gradient,
        values,
        constitutive_variables,
        false,
        ConstitutiveLaw::StressMeasure_PK2);

    //-----------------------------------------------------------------------------------
    Vector volume_force_local = this->GetValue(BODY_FORCE);
    // // Calculating the local RHS
    for ( IndexType i = 0; i < number_of_control_points; ++i ) {
        const SizeType index = mDim * i;

        for (IndexType j = 0; j < mDim; ++j) {
            rRightHandSideVector[index + j] +=
                int_to_reference_weight * N[i] * volume_force_local[j];
        }
    }

    // RHS = ExtForces - K*temp;
    noalias(rRightHandSideVector) -= int_to_reference_weight * prod(
        trans(B), constitutive_variables.StressVector);

    KRATOS_CATCH("")
}

void SolidGreenLagrangeElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != mDim * number_of_control_points)
            rResult.resize(mDim * number_of_control_points, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * mDim;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            if (mDim == 3) {
                rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
            }
        }

        KRATOS_CATCH("")
    };

    void SolidGreenLagrangeElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(mDim * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            if (mDim == 3) {
                rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            }
        }

        KRATOS_CATCH("")
    };



int SolidGreenLagrangeElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Verify that the constitutive law exists
    if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
    {
        KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }
    else
    {
        // Verify that the constitutive law has the correct dimension
        if (mDim == 2 && this->GetProperties()[CONSTITUTIVE_LAW]->GetStrainSize() != 3)
        {
            KRATOS_ERROR << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
                         << this->Id() << std::endl;
        }
        else if (mDim == 3 && this->GetProperties()[CONSTITUTIVE_LAW]->GetStrainSize() != 6)
        {
            KRATOS_ERROR << "Wrong constitutive law used. This is a 3D element! Expected strain size is 6 (el id = ) "
                         << this->Id() << std::endl;
        }

        if (mDim == 2 && !this->GetProperties().Has(THICKNESS))
        {
            KRATOS_ERROR << "THICKNESS not provided for element " << this->Id() << std::endl;
        }
    }

    return Element::Check(rCurrentProcessInfo);
}


Element::IntegrationMethod SolidGreenLagrangeElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void SolidGreenLagrangeElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    Matrix jacobian;
    CalculateInitialJacobian(r_geometry, jacobian);
    Matrix inverse_jacobian(mDim, mDim);
    double det_jacobian = 0.0;
    MathUtils<double>::InvertMatrix(jacobian, inverse_jacobian, det_jacobian);
    const Matrix DN_DX = prod(r_DN_De[0], inverse_jacobian);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    KRATOS_ERROR_IF(mDim != 2 && mDim != 3)
        << "SolidGreenLagrangeElement only supports 2D and 3D elements, but the current dimension is "
        << mDim << std::endl;
    KRATOS_ERROR_IF(strain_size != (mDim == 2 ? 3 : 6))
        << "SolidGreenLagrangeElement has an incompatible constitutive strain size." << std::endl;

    const Matrix deformation_gradient = CalculateDeformationGradient(DN_DX);
    const double det_deformation_gradient = MathUtils<double>::Det(deformation_gradient);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_deformation_gradient) && det_deformation_gradient > 0.0)
        << "SolidGreenLagrangeElement has an invalid deformation gradient determinant: "
        << det_deformation_gradient << std::endl;

    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(strain_size);
    ApplyConstitutiveLaw(
        deformation_gradient,
        det_deformation_gradient,
        values,
        constitutive_variables,
        false,
        ConstitutiveLaw::StressMeasure_Cauchy);

    const Vector cauchy_stress = values.GetStressVector();
    const Matrix cauchy_stress_tensor = MathUtils<double>::StressVectorToTensor(cauchy_stress);
    ApplyConstitutiveLaw(
        deformation_gradient,
        det_deformation_gradient,
        values,
        constitutive_variables,
        false,
        ConstitutiveLaw::StressMeasure_PK2);
    const Vector& r_pk2_stress = values.GetStressVector();
    const Matrix pk2_stress_tensor = MathUtils<double>::StressVectorToTensor(r_pk2_stress);
    SetValue(PK2_STRESS_VECTOR, r_pk2_stress);
    SetValue(PK2_STRESS_TENSOR, pk2_stress_tensor);
    const Matrix right_cauchy_green = prod(trans(deformation_gradient), deformation_gradient);
    Matrix green_lagrange_strain_tensor(mDim, mDim);
    for (IndexType i = 0; i < mDim; ++i) {
        for (IndexType j = 0; j < mDim; ++j) {
            green_lagrange_strain_tensor(i, j) = 0.5 * (
                right_cauchy_green(i, j) - (i == j ? 1.0 : 0.0));
        }
    }
    SetValue(GREEN_LAGRANGE_STRAIN_TENSOR, green_lagrange_strain_tensor);

    if (mDim == 2) {
        SetValue(CAUCHY_STRESS_XX, cauchy_stress[0]);
        SetValue(CAUCHY_STRESS_YY, cauchy_stress[1]);
        SetValue(CAUCHY_STRESS_XY, cauchy_stress[2]);
        SetValue(PK2_STRESS_XX, r_pk2_stress[0]);
        SetValue(PK2_STRESS_YY, r_pk2_stress[1]);
        SetValue(PK2_STRESS_XY, r_pk2_stress[2]);
        const double sigma_xx = cauchy_stress[0];
        const double sigma_yy = cauchy_stress[1];
        const double sigma_xy = cauchy_stress[2];
        const double von_mises_argument = sigma_xx * sigma_xx - sigma_xx * sigma_yy + sigma_yy * sigma_yy + 3.0 * sigma_xy * sigma_xy;
        const double von_mises_stress = von_mises_argument > 0.0 ? std::sqrt(von_mises_argument) : 0.0;
        SetValue(VON_MISES_STRESS_IGA, von_mises_stress);
    }
    else if (mDim == 3) {
        SetValue(CAUCHY_STRESS_TENSOR, cauchy_stress_tensor);
    }

    if (mpConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
        mpConstitutiveLaw->FinalizeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

void SolidGreenLagrangeElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    KRATOS_ERROR_IF(mDim != 2 && mDim != 3)
        << "SolidGreenLagrangeElement only supports 2D and 3D elements, but the current dimension is "
        << mDim << std::endl;
    KRATOS_ERROR_IF(strain_size != (mDim == 2 ? 3 : 6))
        << "SolidGreenLagrangeElement has an incompatible constitutive strain size." << std::endl;

    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(GetIntegrationMethod());
    Matrix jacobian;
    CalculateInitialJacobian(r_geometry, jacobian);
    Matrix inverse_jacobian(mDim, mDim);
    double det_jacobian = 0.0;
    MathUtils<double>::InvertMatrix(jacobian, inverse_jacobian, det_jacobian);
    const Matrix DN_DX = prod(r_DN_De[0], inverse_jacobian);

    const Matrix deformation_gradient = CalculateDeformationGradient(DN_DX);
    const double det_deformation_gradient = MathUtils<double>::Det(deformation_gradient);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_deformation_gradient) && det_deformation_gradient > 0.0)
        << "SolidGreenLagrangeElement has an invalid deformation gradient determinant: "
        << det_deformation_gradient << std::endl;

    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(strain_size);
    Flags& r_options = values.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    values.SetDeformationGradientF(deformation_gradient);
    values.SetDeterminantF(det_deformation_gradient);
    values.SetStrainVector(constitutive_variables.StrainVector);
    values.SetStressVector(constitutive_variables.StressVector);
    values.SetConstitutiveMatrix(constitutive_variables.D);

    const Vector& r_knot_span_sizes = GetValue(KNOT_SPAN_SIZES);
    const double characteristic_geometry_length = 0.5 * std::sqrt(
        r_knot_span_sizes[0] * r_knot_span_sizes[0] +
        r_knot_span_sizes[1] * r_knot_span_sizes[1]);
    values.SetCharacteristicGeometryLength(characteristic_geometry_length);

    if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) 
        if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) {
            mpConstitutiveLaw->InitializeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
        }
}



void SolidGreenLagrangeElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();

    if (rOutput.size() != r_integration_points.size())
    {
        rOutput.resize(r_integration_points.size());
    }

    if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else if (rVariable.Name() == "EQUIVALENT_PLASTIC_STRAIN") {
        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        Matrix jacobian;
        CalculateInitialJacobian(r_geometry, jacobian);
        Matrix inverse_jacobian(mDim, mDim);
        double determinant;
        MathUtils<double>::InvertMatrix(jacobian, inverse_jacobian, determinant);
        const auto& r_gradients = r_geometry.ShapeFunctionsLocalGradients(GetIntegrationMethod());
        Matrix gradients = prod(r_gradients[0], inverse_jacobian);
        Matrix B = ZeroMatrix(strain_size, r_geometry.size() * mDim);
        CalculateB(B, gradients);
        Vector displacement(r_geometry.size() * mDim);
        GetSolutionCoefficientVector(displacement);
        Vector strain = prod(B, displacement);
        Vector stress = ZeroVector(strain_size);
        Matrix constitutive_matrix = ZeroMatrix(strain_size, strain_size);
        ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
        values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        values.SetStrainVector(strain);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(constitutive_matrix);
        const Vector& r_span_sizes = GetValue(KNOT_SPAN_SIZES);
        values.SetCharacteristicGeometryLength(
            0.5 * std::sqrt(r_span_sizes[0] * r_span_sizes[0] + r_span_sizes[1] * r_span_sizes[1]));
        mpConstitutiveLaw->CalculateValue(values, rVariable, rOutput[0]);
    } else {
        KRATOS_WATCH(rVariable);
        KRATOS_WARNING("VARIABLE PRINT STILL NOT IMPLEMENTED N THE IGA FRAMEWORK");
    }
}
void SolidGreenLagrangeElement::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 >>& rVariable,
        std::vector<array_1d<double, 3 >>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

    if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else {
            KRATOS_WATCH(rVariable);
            KRATOS_WARNING("VARIABLE PRINT STILL NOT IMPLEMENTED IN THE IGA FRAMEWORK");
    }
}


void SolidGreenLagrangeElement::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * mDim;
        const SizeType strain_size = (mDim == 2) ? 3 : 6;

        if (rB.size1() != strain_size || rB.size2() != mat_size)
            rB.resize(strain_size, mat_size, false);

        noalias(rB) = ZeroMatrix(strain_size, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const SizeType index = i * mDim;

            if (mDim == 2) {
                rB(0, index + 0) = r_DN_DX(i, 0); // exx
                rB(1, index + 1) = r_DN_DX(i, 1); // eyy

                rB(2, index + 0) = r_DN_DX(i, 1); // gamma_xy
                rB(2, index + 1) = r_DN_DX(i, 0);
            }
            else if (mDim == 3) {
                rB(0, index + 0) = r_DN_DX(i, 0); // exx
                rB(1, index + 1) = r_DN_DX(i, 1); // eyy
                rB(2, index + 2) = r_DN_DX(i, 2); // ezz

                rB(3, index + 0) = r_DN_DX(i, 1); // gamma_xy
                rB(3, index + 1) = r_DN_DX(i, 0);

                rB(4, index + 1) = r_DN_DX(i, 2); // gamma_yz
                rB(4, index + 2) = r_DN_DX(i, 1);

                rB(5, index + 0) = r_DN_DX(i, 2); // gamma_xz
                rB(5, index + 2) = r_DN_DX(i, 0);
            }
        }
    }



void SolidGreenLagrangeElement::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * mDim;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3>& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * mDim;

            for (IndexType j = 0; j < mDim; ++j) {
                rValues[index + j] = displacement[j];
            }
        }
    }

void SolidGreenLagrangeElement::ApplyConstitutiveLaw(
    const Matrix& rDeformationGradient,
    const double DeterminantF,
    ConstitutiveLaw::Parameters& rValues,
    ConstitutiveVariables& rConstitutiveVariables,
    const bool ComputeConstitutiveTensor,
    const ConstitutiveLaw::StressMeasure StressMeasure)
{
    Flags& r_options = rValues.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, ComputeConstitutiveTensor);

    rValues.SetDeformationGradientF(rDeformationGradient);
    rValues.SetDeterminantF(DeterminantF);
    rValues.SetStrainVector(rConstitutiveVariables.StrainVector);
    rValues.SetStressVector(rConstitutiveVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiveVariables.D);

    const Vector& r_knot_span_sizes = GetValue(KNOT_SPAN_SIZES);
    const double characteristic_geometry_length = 0.5 * std::sqrt(
        r_knot_span_sizes[0] * r_knot_span_sizes[0] +
        r_knot_span_sizes[1] * r_knot_span_sizes[1]);
    rValues.SetCharacteristicGeometryLength(characteristic_geometry_length);

    mpConstitutiveLaw->CalculateMaterialResponse(rValues, StressMeasure);
}

} // Namespace Kratos
