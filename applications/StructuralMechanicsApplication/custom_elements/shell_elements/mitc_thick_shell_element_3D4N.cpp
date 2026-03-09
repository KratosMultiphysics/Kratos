// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Massimo Petracca
//                   Alejandro Cornejo
//

#include "mitc_thick_shell_element_3D4N.hpp"

#include <string>
#include <iomanip>

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
MITCThickShellElement3D4N<TKinematics>::MITC4Params::MITC4Params(const ShellQ4_LocalCoordinateSystem& LCS)
    : Transformation(2, 2)
    , ShearStrains(4, 24, 0.0)
{

    const double x21 = LCS.X2() - LCS.X1();
    const double y21 = LCS.Y2() - LCS.Y1();
    const double x34 = LCS.X3() - LCS.X4();
    const double y34 = LCS.Y3() - LCS.Y4();
    const double x41 = LCS.X4() - LCS.X1();
    const double y41 = LCS.Y4() - LCS.Y1();
    const double x32 = LCS.X3() - LCS.X2();
    const double y32 = LCS.Y3() - LCS.Y2();

    Ax = - LCS.X1() + LCS.X2() + LCS.X3() - LCS.X4();
    Bx =   LCS.X1() - LCS.X2() + LCS.X3() - LCS.X4();
    Cx = - LCS.X1() - LCS.X2() + LCS.X3() + LCS.X4();
    Ay = - LCS.Y1() + LCS.Y2() + LCS.Y3() - LCS.Y4();
    By =   LCS.Y1() - LCS.Y2() + LCS.Y3() - LCS.Y4();
    Cy = - LCS.Y1() - LCS.Y2() + LCS.Y3() + LCS.Y4();

    double Alpha = std::atan(Ay / Ax);
    double Beta  = Globals::Pi * 0.5 - std::atan(Cx / Cy);

    Transformation(0, 0) =   std::sin(Beta);
    Transformation(0, 1) = - std::sin(Alpha);
    Transformation(1, 0) = - std::cos(Beta);
    Transformation(1, 1) =   std::cos(Alpha);

    ShearStrains(0, 2) = -0.5;
    ShearStrains(0, 3) = -y41 * 0.25;
    ShearStrains(0, 4) = x41 * 0.25;

    ShearStrains(0, 20) = 0.5;
    ShearStrains(0, 21) = ShearStrains(0, 3);
    ShearStrains(0, 22) = ShearStrains(0, 4);

    ShearStrains(1, 2) = -0.5;
    ShearStrains(1, 3) = -y21 * 0.25;
    ShearStrains(1, 4) = x21 * 0.25;

    ShearStrains(1, 8) = 0.5;
    ShearStrains(1, 9) = ShearStrains(1, 3);
    ShearStrains(1, 10) = ShearStrains(1, 4);

    ShearStrains(2, 8) = -0.5;
    ShearStrains(2, 9) = -y32 * 0.25;
    ShearStrains(2, 10) = x32 * 0.25;

    ShearStrains(2, 14) = 0.5;
    ShearStrains(2, 15) = ShearStrains(2, 9);
    ShearStrains(2, 16) = ShearStrains(2, 10);

    ShearStrains(3, 14) = 0.5;
    ShearStrains(3, 15) = -y34 * 0.25;
    ShearStrains(3, 16) = x34 * 0.25;

    ShearStrains(3, 20) = -0.5;
    ShearStrains(3, 21) = ShearStrains(3, 15);
    ShearStrains(3, 22) = ShearStrains(3, 16);
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
MITCThickShellElement3D4N<TKinematics>::EASOperatorStorage::EASOperatorStorage()
    : mInitialized(false)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperatorStorage::Initialize(const GeometryType& geom)
{
    if (!mInitialized) {
        noalias(alpha) = ZeroVector(5);
        noalias(alpha_converged) = ZeroVector(5);

        array_1d<double, 3> initial_displacement, initial_rotation;

        for (SizeType i = 0; i < 4; i++) {
            SizeType ii = i * 6;
            noalias(initial_displacement) = geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            noalias(initial_rotation) = geom[i].FastGetSolutionStepValue(ROTATION);

            displ(ii) = initial_displacement[0];
            displ(ii + 1) = initial_displacement[1];
            displ(ii + 2) = initial_displacement[2];
            displ_converged(ii) = initial_displacement[0];
            displ_converged(ii + 1) = initial_displacement[1];
            displ_converged(ii + 2) = initial_displacement[2];
            displ(ii + 3) = initial_rotation[0];
            displ(ii + 4) = initial_rotation[1];
            displ(ii + 5) = initial_rotation[2];
            displ_converged(ii + 3) = initial_rotation[0];
            displ_converged(ii + 4) = initial_rotation[1];
            displ_converged(ii + 5) = initial_rotation[2];
        }

        mInitialized = true;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperatorStorage::InitializeSolutionStep()
{
    noalias(displ) = displ_converged;
    noalias(alpha) = alpha_converged;
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperatorStorage::FinalizeSolutionStep()
{
    noalias(displ_converged) = displ;
    noalias(alpha_converged) = alpha;
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperatorStorage::FinalizeNonLinearIteration(const Vector& displacementVector)
{
    array_1d<double, 24> incrementalDispl;
    noalias(incrementalDispl) = displacementVector - displ;
    noalias(displ) = displacementVector;

    array_1d<double, 5> temp;
    noalias(temp) = prod(L, incrementalDispl);
    noalias(temp) -= residual;
    noalias(alpha) -= prod(Hinv, temp);
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperatorStorage::save(Serializer& rSerializer) const
{
    rSerializer.save("A0", alpha);
    rSerializer.save("A1", alpha_converged);
    rSerializer.save("U0", displ);
    rSerializer.save("U1", displ_converged);
    rSerializer.save("res", residual);
    rSerializer.save("Hinv", Hinv);
    rSerializer.save("mL", L);
    rSerializer.save("init", mInitialized);
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperatorStorage::load(Serializer& rSerializer)
{
    rSerializer.load("A0", alpha);
    rSerializer.load("A1", alpha_converged);
    rSerializer.load("U0", displ);
    rSerializer.load("U1", displ_converged);
    rSerializer.load("res", residual);
    rSerializer.load("Hinv", Hinv);
    rSerializer.load("mL", L);
    rSerializer.load("init", mInitialized);
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
MITCThickShellElement3D4N<TKinematics>::EASOperator::EASOperator(const ShellQ4_LocalCoordinateSystem& LCS, EASOperatorStorage& storage)
    : mF0inv(3, 3)
    , mEnhancedStrains(3)
    , mG(3, 5)
{

    // compute the jacobian at the element center

    double xi(0.0);
    double eta(0.0);

    Matrix dN(4, 2);
    ShellUtilities::ShapeFunc_NaturalDerivatives(xi, eta, dN);

    Matrix Jac0(2, 2);
    Jac0(0, 0) = dN(0, 0) * LCS.X1() + dN(1, 0) * LCS.X2() + dN(2, 0) * LCS.X3() + dN(3, 0) * LCS.X4();
    Jac0(0, 1) = dN(0, 0) * LCS.Y1() + dN(1, 0) * LCS.Y2() + dN(2, 0) * LCS.Y3() + dN(3, 0) * LCS.Y4();
    Jac0(1, 0) = dN(0, 1) * LCS.X1() + dN(1, 1) * LCS.X2() + dN(2, 1) * LCS.X3() + dN(3, 1) * LCS.X4();
    Jac0(1, 1) = dN(0, 1) * LCS.Y1() + dN(1, 1) * LCS.Y2() + dN(2, 1) * LCS.Y3() + dN(3, 1) * LCS.Y4();

    // save the jacobian determinant at center
    mJ0 = Jac0(0, 0) * Jac0(1, 1) - Jac0(1, 0) * Jac0(0, 1);

    // compute the transformation matrix used in the implementation of the EAS method
    // which operates in the natural coordinate system
    const double j11 = Jac0(0, 0);
    const double j22 = Jac0(1, 1);
    const double j12 = Jac0(0, 1);
    const double j21 = Jac0(1, 0);

    Matrix F0(3, 3);
    F0(0, 0) = j11 * j11;
    F0(0, 1) = j21 * j12;
    F0(0, 2) = 2.0 * j11 * j12;
    F0(1, 0) = j12 * j21;
    F0(1, 1) = j22 * j22;
    F0(1, 2) = 2.0 * j21 * j22;
    F0(2, 0) = j11 * j21;
    F0(2, 1) = j12 * j22;
    F0(2, 2) = j11 * j22 + j12 * j21;

    double dummyDet;
    MathUtils<double>::InvertMatrix3(F0, mF0inv, dummyDet);

    // Initialize these data to zero because they will
    // be integrated during the gauss loop
    storage.L.clear();
    storage.Hinv.clear();
    storage.residual.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperator::GaussPointComputation_Step1(
    double xi,
    double eta,
    const ShellUtilities::JacobianOperator& jac,
    Vector& generalizedStrains,
    EASOperatorStorage& storage)
{
    // construct the interpolation matrix in natural coordinate system
    Matrix E(3, 5, 0.0);
    E(0, 0) = xi;
    E(1, 1) = eta;
    E(2, 2) = xi;
    E(2, 3) = eta;
    E(0, 4) = xi * eta;
    E(1, 4) = -xi * eta;
    E(2, 4) = xi * xi - eta * eta;

    // construct the interpolation matrix in local coordinate system
    noalias(mG) = mJ0 / jac.Determinant() * prod(mF0inv, E);

    // Assuming that the input generalized strains has been already calculated on this
    // integration point, we just need to add the contribution of the enhanced strains
    noalias(mEnhancedStrains) = prod(mG, storage.alpha);

    generalizedStrains[0] += mEnhancedStrains[0]; // e.xx
    generalizedStrains[1] += mEnhancedStrains[1]; // e.yy
    generalizedStrains[2] += mEnhancedStrains[2]; // e.xy
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperator::GaussPointComputation_Step2(
    const Matrix& D,
    const Matrix& B,
    const Vector& S,
    EASOperatorStorage& storage)
{
    Matrix GTC(5, 3);
    noalias(GTC) = prod(trans(mG), project(D, range(0, 3), range(0, 3)));
    noalias(storage.Hinv) += prod(GTC, mG);
    noalias(storage.residual) -= prod(trans(mG), project(S, range(0, 3)));

    // compute L: [G'*Dmm, G'*Dmb, G'*Dms]*[Bm; Bm; Bs]
    int num_stress = D.size2(); // it can be 6 for thin shells or 8 for thick  shells
    Matrix GTD(5, num_stress, 0.0);
    noalias(project(GTD, range::all(), range(0, 3))) = GTC;
    noalias(project(GTD, range::all(), range(3, 6))) = prod(trans(mG), project(D, range(0, 3), range(3, 6)));
    if (num_stress == 8)
    {
        noalias(project(GTD, range::all(), range(6, 8))) = prod(trans(mG), project(D, range(0, 3), range(6, 8)));
    }
    noalias(storage.L) += prod(GTD, B);
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::EASOperator::ComputeModfiedTangentAndResidual(Matrix& rLeftHandSideMatrix,
        Vector& rRightHandSideVector,
        EASOperatorStorage& storage)
{
    // invert H
    Matrix Hcopy(storage.Hinv);
    permutation_matrix<Matrix::size_type> pm(5);
    lu_factorize(Hcopy, pm);
    noalias(storage.Hinv) = IdentityMatrix(5);
    lu_substitute(Hcopy, pm, storage.Hinv);

    // compute L' * H^-1
    Matrix LTHinv(24, 5);
    noalias(LTHinv) = prod(trans(storage.L), storage.Hinv);

    // modify the stiffness matrix and the residual vector
    // for the static condensation of the enhanced strain parameters
    noalias(rRightHandSideVector) -= prod(LTHinv, storage.residual); // R_mod = R - L' * H^-1 * residual
    noalias(rLeftHandSideMatrix) -= prod(LTHinv, storage.L);         // K_mod = K - L' * H^-1 * L
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
MITCThickShellElement3D4N<TKinematics>::MITCThickShellElement3D4N(IndexType NewId,
        GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
MITCThickShellElement3D4N<TKinematics>::MITCThickShellElement3D4N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
Element::Pointer MITCThickShellElement3D4N<TKinematics>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
    return Kratos::make_intrusive< MITCThickShellElement3D4N >(NewId, newGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
Element::Pointer MITCThickShellElement3D4N<TKinematics>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MITCThickShellElement3D4N >(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Initialize(rCurrentProcessInfo);

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        const auto& r_geometry = GetGeometry();
        mEASStorage.Initialize(r_geometry);
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());

        const auto& r_properties = GetProperties();
        if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
            auto N_values = Vector();
            for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
                mConstitutiveLawVector[point_number] = r_properties[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
            }
        } else {
            KRATOS_ERROR << "A constitutive law needs to be specified for the CS-DSG3 shell element with ID " << this->Id() << std::endl;
        }
    }

    KRATOS_CATCH("Initialize")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::CalculateMaterialResponse(
        ShellCrossSection::SectionParameters& rSectionParameters,
        const SizeType& rIntegrationPointNumber,
        const bool CalculateStress,
        const bool CalculateConstitutive,
        const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, CalculateStress);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, CalculateConstitutive);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    cl_values.SetStrainVector(rSectionParameters.GetGeneralizedStrainVector());
    cl_values.SetStressVector(rSectionParameters.GetGeneralizedStressVector());
    cl_values.SetConstitutiveMatrix(rSectionParameters.GetConstitutiveMatrix());

    mConstitutiveLawVector[rIntegrationPointNumber]->CalculateMaterialResponseCauchy(cl_values);

    KRATOS_CATCH("CalculateMaterialResponse")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::FinalizeMaterialResponse(
        ShellCrossSection::SectionParameters& rSectionParameters,
        const SizeType& rIntegrationPointNumber,
        const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    cl_values.SetStrainVector(rSectionParameters.GetGeneralizedStrainVector());
    cl_values.SetStressVector(rSectionParameters.GetGeneralizedStressVector());
    cl_values.SetConstitutiveMatrix(rSectionParameters.GetConstitutiveMatrix());

    mConstitutiveLawVector[rIntegrationPointNumber]->FinalizeMaterialResponseCauchy(cl_values);
    KRATOS_CATCH("FinalizeMaterialResponse")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::InitializeMaterialResponse(
        ShellCrossSection::SectionParameters& rSectionParameters,
        const SizeType& rIntegrationPointNumber,
        const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    cl_values.SetStrainVector(rSectionParameters.GetGeneralizedStrainVector());
    cl_values.SetStressVector(rSectionParameters.GetGeneralizedStressVector());
    cl_values.SetConstitutiveMatrix(rSectionParameters.GetConstitutiveMatrix());

    mConstitutiveLawVector[rIntegrationPointNumber]->InitializeMaterialResponseCauchy(cl_values);
    KRATOS_CATCH("FinalizeMaterialResponse")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::FinalizeNonLinearIteration(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    ShellQ4_LocalCoordinateSystem LCS(this->mpCoordinateTransformation->CreateLocalCoordinateSystem());
    Vector globalDisplacementVector(24);
    this->GetValuesVector(globalDisplacementVector);
    Vector localDisplacementVector(this->mpCoordinateTransformation->CalculateLocalDisplacements(LCS, globalDisplacementVector));

    mEASStorage.FinalizeNonLinearIteration(localDisplacementVector);
    KRATOS_CATCH("FinalizeNonLinearIteration")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    bool required = false;
    for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
        if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
            required = true;
            break;
        }
    }
    // Finalize material response if required
    if (required) {

        // Get some references.
        const auto& r_props = GetProperties();
        const auto& r_geom = GetGeometry();
        const Matrix& shapeFunctions = r_geom.ShapeFunctionsValues();
        const double thickness = r_props[THICKNESS];
        Vector iN(shapeFunctions.size2());

        // Compute the local coordinate system.
        ShellQ4_LocalCoordinateSystem localCoordinateSystem(
            this->mpCoordinateTransformation->CreateLocalCoordinateSystem());

        ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
            this->mpCoordinateTransformation->CreateReferenceCoordinateSystem());

        // Prepare all the parameters needed for the MITC formulation.
        // This is to be done here outside the Gauss Loop.
        MITC4Params shearParameters(referenceCoordinateSystem);

        ShellUtilities::JacobianOperator jacOp;
        array_1d<double, 4> dArea;

        // Instantiate all strain-displacement matrices.
        Matrix B(8, 24, 0.0);
        Vector Bdrilling(24, 0.0);

        // // Instantiate all section tangent matrices.
        Matrix D(8, 8, 0.0);

        // Instantiate strain and stress-resultant vectors
        Vector generalizedStrains(8);
        Vector generalizedStresses(8);

        // Get the current displacements in global coordinate system

        Vector globalDisplacements(24);
        this->GetValuesVector(globalDisplacements, 0);

        // Get the current displacements in local coordinate system
        Vector localDisplacements(this->mpCoordinateTransformation->CalculateLocalDisplacements(localCoordinateSystem, globalDisplacements));

        // Instantiate the EAS Operator.
        // This will apply the Enhanced Assumed Strain Method for the calculation
        // of the membrane contribution.
        EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

        ShellCrossSection::SectionParameters parameters(r_geom, r_props, rCurrentProcessInfo);
        parameters.SetGeneralizedStrainVector(generalizedStrains);
        parameters.SetGeneralizedStressVector(generalizedStresses);
        parameters.SetConstitutiveMatrix(D);
        Flags& r_options = parameters.GetOptions();
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Gauss Loop.
        for (SizeType integration_point = 0; integration_point < 4; integration_point++) {

            // Get a reference of the current integration point and shape functions
            const GeometryType::IntegrationPointType& ip = r_geom.IntegrationPoints()[integration_point];
            noalias(iN) = row(shapeFunctions, integration_point);

            // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
            // and Shape functions derivatives in the local coordinate system
            jacOp.Calculate(referenceCoordinateSystem, r_geom.ShapeFunctionLocalGradient(integration_point));

            // Compute all strain-displacement matrices
            CalculateBMatrix(ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling);

            // Calculate strain vectors in local coordinate system

            noalias(generalizedStrains) = prod(B, localDisplacements);

            // Apply the EAS method to modify the membrane part of the strains computed above.
            EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains, mEASStorage);

            // Calculate the response of the Cross Section
            parameters.SetShapeFunctionsValues(iN);
            parameters.SetShapeFunctionsDerivatives(jacOp.XYDerivatives());
            InitializeMaterialResponse(parameters, integration_point, rCurrentProcessInfo);
        }
    }
    BaseType::InitializeSolutionStep(rCurrentProcessInfo);  // TODO remove in the future and move mpCoordinateTransformation->FinalizeSolutionStep();
    mEASStorage.InitializeSolutionStep();

    KRATOS_CATCH("InitializeSolutionStep")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    bool required = false;
    for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
        if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
            required = true;
            break;
        }
    }

    // Finalize material response if required
    if (required) {

        // Get some references.
        const auto& r_props = GetProperties();
        const auto& r_geom = GetGeometry();
        const Matrix& shapeFunctions = r_geom.ShapeFunctionsValues();
        const double thickness = r_props[THICKNESS];
        Vector iN(shapeFunctions.size2());

        // Compute the local coordinate system.
        ShellQ4_LocalCoordinateSystem localCoordinateSystem(
            this->mpCoordinateTransformation->CreateLocalCoordinateSystem());

        ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
            this->mpCoordinateTransformation->CreateReferenceCoordinateSystem());

        // Prepare all the parameters needed for the MITC formulation.
        // This is to be done here outside the Gauss Loop.
        MITC4Params shearParameters(referenceCoordinateSystem);

        // Instantiate the Jacobian Operator.
        // This will store:
        // the jacobian matrix, its inverse, its determinant
        // and the derivatives of the shape functions in the local
        // coordinate system
        ShellUtilities::JacobianOperator jacOp;
        array_1d<double, 4> dArea;


        // Instantiate all strain-displacement matrices.
        Matrix B(8, 24, 0.0);
        Vector Bdrilling(24, 0.0);

        // // Instantiate all section tangent matrices.
        Matrix D(8, 8, 0.0);

        // Instantiate strain and stress-resultant vectors
        Vector generalizedStrains(8);
        Vector generalizedStresses(8);

        // Get the current displacements in global coordinate system

        Vector globalDisplacements(24);
        this->GetValuesVector(globalDisplacements, 0);

        // Get the current displacements in local coordinate system
        Vector localDisplacements(this->mpCoordinateTransformation->CalculateLocalDisplacements(localCoordinateSystem, globalDisplacements));

        // Instantiate the EAS Operator.
        // This will apply the Enhanced Assumed Strain Method for the calculation
        // of the membrane contribution.
        EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

        ShellCrossSection::SectionParameters parameters(r_geom, r_props, rCurrentProcessInfo);
        parameters.SetGeneralizedStrainVector(generalizedStrains);
        parameters.SetGeneralizedStressVector(generalizedStresses);
        parameters.SetConstitutiveMatrix(D);
        Flags& r_options = parameters.GetOptions();
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Gauss Loop.
        for (SizeType integration_point = 0; integration_point < 4; integration_point++) {

            // Get a reference of the current integration point and shape functions
            const GeometryType::IntegrationPointType& ip = r_geom.IntegrationPoints()[integration_point];
            noalias(iN) = row(shapeFunctions, integration_point);

            // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
            // and Shape functions derivatives in the local coordinate system
            jacOp.Calculate(referenceCoordinateSystem, r_geom.ShapeFunctionLocalGradient(integration_point));

            // Compute all strain-displacement matrices
            CalculateBMatrix(ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling);

            // Calculate strain vectors in local coordinate system
            noalias(generalizedStrains) = prod(B, localDisplacements);

            // Apply the EAS method to modify the membrane part of the strains computed above.
            EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains, mEASStorage);

            // Calculate the response of the Cross Section
            parameters.SetShapeFunctionsValues(iN);
            parameters.SetShapeFunctionsDerivatives(jacOp.XYDerivatives());
            FinalizeMaterialResponse(parameters, integration_point, rCurrentProcessInfo);
        }
    }

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo); // TODO remove in the future and move mpCoordinateTransformation->FinalizeSolutionStep();
    mEASStorage.FinalizeSolutionStep();

    KRATOS_CATCH("FinalizeSolutionStep")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_integration_points = GetGeometry().IntegrationPoints();
    const SizeType number_of_integration_points = r_integration_points.size();

    // Provide a default empty implementation: resize and set zeros
    rOutput.assign(number_of_integration_points, 0.0);

    if (mConstitutiveLawVector[0]->Has(rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("MITCThickShellElement3D4N::CalculateOnIntegrationPoints")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_integration_points = GetGeometry().IntegrationPoints();
    const SizeType number_of_integration_points = r_integration_points.size();

    // Provide a default empty implementation: resize and set zeros
    rOutput.assign(number_of_integration_points, ZeroVector(3)); // plane stress components

    if (mConstitutiveLawVector[0]->Has(rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("MITCThickShellElement3D4N::CalculateOnIntegrationPoints")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
int MITCThickShellElement3D4N<TKinematics>::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    mConstitutiveLawVector[0]->Check(GetProperties(), r_geom, rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT((r_geom.IntegrationPoints(this->GetIntegrationMethod())).size() == 4) << "MITCThickShellElement3D4N - needs a full integration scheme" << std::endl;

    const int points_number = r_geom.PointsNumber();
    KRATOS_ERROR_IF_NOT(points_number == 4) << "MITCThickShellElement3D4N - Wrong number of nodes" << points_number << std::endl;

    return 0;

    KRATOS_CATCH("Check")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::CalculateBMatrix(
    double xi,
    double eta,
    const ShellUtilities::JacobianOperator& Jac,
    const MITC4Params& mitc_params,
    const Vector& N,
    Matrix& B,
    Vector& Bdrill)
{
    KRATOS_TRY

    const Matrix& dNxy = Jac.XYDerivatives();

    // Membrane
    B(0,  0) =  dNxy(0, 0);
    B(0,  6) =  dNxy(1, 0);
    B(0, 12) =  dNxy(2, 0);
    B(0, 18) =  dNxy(3, 0);
    B(1,  1) =  dNxy(0, 1);
    B(1,  7) =  dNxy(1, 1);
    B(1, 13) =  dNxy(2, 1);
    B(1, 19) =  dNxy(3, 1);
    B(2,  0) =  dNxy(0, 1);
    B(2,  6) =  dNxy(1, 1);
    B(2, 12) =  dNxy(2, 1);
    B(2, 18) =  dNxy(3, 1);
    B(2,  1) =  dNxy(0, 0);
    B(2,  7) =  dNxy(1, 0);
    B(2, 13) =  dNxy(2, 0);
    B(2, 19) =  dNxy(3, 0);

    // Bending
    B(3,  4) =  dNxy(0, 0);
    B(3, 10) =  dNxy(1, 0);
    B(3, 16) =  dNxy(2, 0);
    B(3, 22) =  dNxy(3, 0);
    B(4,  3) = -dNxy(0, 1);
    B(4,  9) = -dNxy(1, 1);
    B(4, 15) = -dNxy(2, 1);
    B(4, 21) = -dNxy(3, 1);
    B(5,  3) = -dNxy(0, 0);
    B(5,  9) = -dNxy(1, 0);
    B(5, 15) = -dNxy(2, 0);
    B(5, 21) = -dNxy(3, 0);
    B(5,  4) =  dNxy(0, 1);
    B(5, 10) =  dNxy(1, 1);
    B(5, 16) =  dNxy(2, 1);
    B(5, 22) =  dNxy(3, 1);

    // Drilling
    Bdrill[0] = -0.5 * dNxy(0, 1);
    Bdrill[1] =  0.5 * dNxy(0, 0);
    Bdrill[5] = -N[0];
    Bdrill[6] = -0.5 * dNxy(1, 1);
    Bdrill[7] =  0.5 * dNxy(1, 0);
    Bdrill[11] = -N[1];
    Bdrill[12] = -0.5 * dNxy(2, 1);
    Bdrill[13] =  0.5 * dNxy(2, 0);
    Bdrill[17] = -N[2];
    Bdrill[18] = -0.5 * dNxy(3, 1);
    Bdrill[19] =  0.5 * dNxy(3, 0);
    Bdrill[23] = -N[3];

    // Shear
    // MITC modified shape functions
    Matrix MITCShapeFunctions(2, 4, 0.0);
    MITCShapeFunctions(1, 0) = 1.0 - xi;
    MITCShapeFunctions(0, 1) = 1.0 - eta;
    MITCShapeFunctions(1, 2) = 1.0 + xi;
    MITCShapeFunctions(0, 3) = 1.0 + eta;

    // Strain displacement matrix in natural coordinate system.
    // Interpolate the shear strains given in MITC4Params
    // Using the modified shape function
    Matrix BN = prod(MITCShapeFunctions, mitc_params.ShearStrains);

    // Modify the shear strain intensity in the tying points
    // to match the values that would be obtained using standard
    // interpolations
    double Temp1, Temp2, Temp3;
    Temp1 = mitc_params.Cx + xi * mitc_params.Bx;
    Temp3 = mitc_params.Cy + xi * mitc_params.By;
    Temp1 = Temp1 * Temp1 + Temp3 * Temp3;
    Temp1 = std::sqrt(Temp1) / (8.0 * Jac.Determinant());
    Temp2 = mitc_params.Ax + eta * mitc_params.Bx;
    Temp3 = mitc_params.Ay + eta * mitc_params.By;
    Temp2 = Temp2 * Temp2 + Temp3 * Temp3;
    Temp2 = std::sqrt(Temp2) / (8.0 * Jac.Determinant());

    row(BN, 0) *= Temp1;
    row(BN, 1) *= Temp2;

    // Transform the strain-displacement matrix from natural
    // to local coordinate system taking into account the element distortion
    project(B, slice(7, -1, 2), slice::all()) = prod(mitc_params.Transformation, BN);
    // Explanation of the 'slice':
    // The MITC4Params class and the first part of this method were coded
    // assuming the following notations for strain vector : [e.xx, e.yy, e.zz, 2 e.xy, 2 e.XZ, 2 e.YZ].
    // Since in Kratos the strain vector is assumed to be : [e.xx, e.yy, e.zz, 2 e.xy, 2 e.YZ, 2 e.XZ],
    // I had to "virtually" swap the 2 rows of the matrix B before the assignment from the product on the Right-Hand-Side.
    // In this way the matrix B is consistent with the notation used in Kratos without recoding all the MITC stuff.

    KRATOS_CATCH("CalculateBMatrix")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    if ((rLeftHandSideMatrix.size1() != 24) || (rLeftHandSideMatrix.size2() != 24)) {
        rLeftHandSideMatrix.resize(24, 24, false);
    }
    rLeftHandSideMatrix.clear();

    if (rRightHandSideVector.size() != 24) {
        rRightHandSideVector.resize(24, false);
    }
    rRightHandSideVector.clear();

    // Get some references.
    const auto& r_props = GetProperties();
    const auto& r_geom = GetGeometry();
    const Matrix& shapeFunctions = r_geom.ShapeFunctionsValues();
    const double thickness = r_props[THICKNESS];
    Vector iN(shapeFunctions.size2());

    // Compute the local coordinate system.
    ShellQ4_LocalCoordinateSystem localCoordinateSystem(
        this->mpCoordinateTransformation->CreateLocalCoordinateSystem());

    ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
        this->mpCoordinateTransformation->CreateReferenceCoordinateSystem());

    // Prepare all the parameters needed for the MITC formulation.
    // This is to be done here outside the Gauss Loop.
    MITC4Params shearParameters(referenceCoordinateSystem);

    // Instantiate the Jacobian Operator.
    // This will store:
    // the jacobian matrix, its inverse, its determinant
    // and the derivatives of the shape functions in the local
    // coordinate system
    ShellUtilities::JacobianOperator jacOp;
    array_1d<double, 4> dArea;

    // Instantiate all strain-displacement matrices.

    Matrix B(8, 24, 0.0);
    Vector Bdrilling(24, 0.0);

    // Instantiate all section tangent matrices.
    Matrix D(8, 8, 0.0);

    // Instantiate strain and stress-resultant vectors
    Vector generalizedStrains(8);
    Vector generalizedStresses(8);

    // Get the current displacements in global coordinate system

    Vector globalDisplacements(24);
    this->GetValuesVector(globalDisplacements, 0);

    // Get the current displacements in local coordinate system
    Vector localDisplacements(this->mpCoordinateTransformation->CalculateLocalDisplacements(localCoordinateSystem, globalDisplacements));

    // Instantiate the EAS Operator.
    // This will apply the Enhanced Assumed Strain Method for the calculation
    // of the membrane contribution.
    EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

    // auxiliary data

    Matrix BTD(24, 8); // auxiliary matrix to store the product B'*D
    // Initialize parameters for the cross section calculation

    ShellCrossSection::SectionParameters parameters(r_geom, r_props, rCurrentProcessInfo);
    parameters.SetGeneralizedStrainVector(generalizedStrains);
    parameters.SetGeneralizedStressVector(generalizedStresses);
    parameters.SetConstitutiveMatrix(D);
    Flags& r_options = parameters.GetOptions();
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, CalculateResidualVectorFlag);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, CalculateStiffnessMatrixFlag);

    // Gauss Loop.
    for (SizeType integration_point = 0; integration_point < 4; integration_point++) {

        // get a reference of the current integration point and shape functions
        const GeometryType::IntegrationPointType& ip = r_geom.IntegrationPoints()[integration_point];
        noalias(iN) = row(shapeFunctions, integration_point);

        // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
        // and Shape functions derivatives in the local coordinate system
        jacOp.Calculate(referenceCoordinateSystem, r_geom.ShapeFunctionLocalGradient(integration_point));

        // compute the 'area' of the current integration point
        const double dA = ip.Weight() * jacOp.Determinant();
        dArea[integration_point] = dA;

        // Compute all strain-displacement matrices

        CalculateBMatrix(ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling);

        // Calculate strain vectors in local coordinate system
        noalias(generalizedStrains) = prod(B, localDisplacements);
        const double drillingStrain = inner_prod(Bdrilling, localDisplacements);

        // Apply the EAS method to modify the membrane part of the strains computed above.
        EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains, mEASStorage);

        // Calculate the response of the Cross Section
        parameters.SetShapeFunctionsValues(iN);
        parameters.SetShapeFunctionsDerivatives(jacOp.XYDerivatives());
        CalculateMaterialResponse(parameters, integration_point, CalculateResidualVectorFlag, CalculateStiffnessMatrixFlag, rCurrentProcessInfo);

        // multiply the section tangent matrices and stress resultants by 'dA'
        D *= dA;
        double Ddrilling = D(2, 2) * (r_props.Has(STABILIZATION_FACTOR) ? r_props[STABILIZATION_FACTOR] : drilling_factor); // drilling stiffness
        generalizedStresses *= dA;
        const double drillingStress = Ddrilling * drillingStrain; // already multiplied by 'dA'

        // Add all contributions to the Stiffness Matrix
        noalias(BTD) = prod(trans(B), D);
        noalias(rLeftHandSideMatrix) += prod(BTD, B);
        noalias(rLeftHandSideMatrix) += outer_prod(Ddrilling * Bdrilling, Bdrilling);

        // Add all contributions to the residual vector
        noalias(rRightHandSideVector) -= prod(trans(B), generalizedStresses);
        noalias(rRightHandSideVector) -= Bdrilling * drillingStress;

        // Continue the calculation of the EAS method now that the contitutive response
        // has been computed
        EASOp.GaussPointComputation_Step2(D, B, generalizedStresses, mEASStorage);
    }

    // Now that the gauss integration is over, let the EAS operator modify
    // the local stiffness matrix and residual vector.
    // It will perform a static condensation to remove the enhanced strain parameters
    // at the element level
    EASOp.ComputeModfiedTangentAndResidual(rLeftHandSideMatrix, rRightHandSideVector, mEASStorage);

    // Let the CoordinateTransformation finalize the calculation.
    // This will handle the transformation of the local matrices/vectors to
    // the global coordinate system.
    this->mpCoordinateTransformation->FinalizeCalculations(
        localCoordinateSystem,
        globalDisplacements,
        localDisplacements,
        rLeftHandSideMatrix,
        rRightHandSideVector,
        CalculateResidualVectorFlag,
        CalculateStiffnessMatrixFlag);
    // Add body forces contributions. This doesn't depend on the coordinate system
    AddBodyForces(dArea, rRightHandSideVector);

    KRATOS_CATCH("CalculateAll")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::AddBodyForces(const array_1d<double,4>& dA, VectorType& rRightHandSideVector)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();

    // Get shape functions
    const Matrix& N = r_geometry.ShapeFunctionsValues();

    // auxiliary
    array_1d<double, 3> bf;

    double density = 0.0;
    if (r_props.Has(DENSITY)) {
        density = r_props[DENSITY];
    } else {
        // In composite shells the density is to be retrieved from the CL
        mConstitutiveLawVector[0]->GetValue(DENSITY, density);
        KRATOS_ERROR_IF(density <= 0.0) << "DENSITY is null as far as the shell CL is concerned... Please implement the GetValue(DENSITY) " <<  std::endl;
    }
    const double mass_per_unit_area = density * r_props[THICKNESS];

    // gauss loop to integrate the external force vector
    for (SizeType igauss = 0; igauss < 4; igauss++) {

        // interpolate nodal volume accelerations to this gauss point
        // and obtain the body force vector
        bf.clear();
        for (SizeType inode = 0; inode < 4; inode++) {
            if (r_geometry[inode].SolutionStepsDataHas(VOLUME_ACCELERATION)) { //temporary, will be checked once at the beginning only
                bf += N(igauss,inode) * r_geometry[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
            }
        }
        bf *= (mass_per_unit_area * dA[igauss]);

        // add it to the RHS vector
        for (SizeType inode = 0; inode < 4; inode++) {
            SizeType index = inode*6;
            double iN = N(igauss,inode);
            rRightHandSideVector[index + 0] += iN * bf[0];
            rRightHandSideVector[index + 1] += iN * bf[1];
            rRightHandSideVector[index + 2] += iN * bf[2];
        }
    }
    KRATOS_CATCH("MITCThickShellElement3D4N::AddBodyForces")
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
ShellCrossSection::SectionBehaviorType MITCThickShellElement3D4N<TKinematics>::GetSectionBehavior() const
{
    return ShellCrossSection::Thick;
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    rSerializer.save("EAS", mEASStorage);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template <ShellKinematics TKinematics>
void MITCThickShellElement3D4N<TKinematics>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    rSerializer.load("EAS", mEASStorage);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template class MITCThickShellElement3D4N<ShellKinematics::LINEAR>;
template class MITCThickShellElement3D4N<ShellKinematics::NONLINEAR_COROTATIONAL>;

/***********************************************************************************/
/***********************************************************************************/

}
