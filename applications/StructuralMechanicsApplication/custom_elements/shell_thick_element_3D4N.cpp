// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Massimo Petracca
//

#include "shell_thick_element_3D4N.hpp"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"
#include "custom_utilities/shell_utilities.h"

#include <string>
#include <iomanip>

namespace Kratos
{
// =====================================================================================
//
// Class JacobianOperator
//
// =====================================================================================

ShellThickElement3D4N::JacobianOperator::JacobianOperator()
    : mJac(2, 2, 0.0)
    , mInv(2, 2, 0.0)
    , mXYDeriv(4, 2, 0.0)
    , mDet(0.0)
{
}

void ShellThickElement3D4N::JacobianOperator::Calculate(const ShellQ4_LocalCoordinateSystem & CS, const Matrix & dN)
{
    mJac(0, 0) = dN(0, 0) * CS.X1() + dN(1, 0) * CS.X2() + dN(2, 0) * CS.X3() + dN(3, 0) * CS.X4();
    mJac(0, 1) = dN(0, 0) * CS.Y1() + dN(1, 0) * CS.Y2() + dN(2, 0) * CS.Y3() + dN(3, 0) * CS.Y4();
    mJac(1, 0) = dN(0, 1) * CS.X1() + dN(1, 1) * CS.X2() + dN(2, 1) * CS.X3() + dN(3, 1) * CS.X4();
    mJac(1, 1) = dN(0, 1) * CS.Y1() + dN(1, 1) * CS.Y2() + dN(2, 1) * CS.Y3() + dN(3, 1) * CS.Y4();

    mDet = mJac(0, 0) * mJac(1, 1) - mJac(1, 0) * mJac(0, 1);
    double mult = 1.0 / mDet;

    mInv(0, 0) =   mJac(1, 1) * mult;
    mInv(0, 1) = - mJac(0, 1) * mult;
    mInv(1, 0) = - mJac(1, 0) * mult;
    mInv(1, 1) =   mJac(0, 0) * mult;

    noalias( mXYDeriv ) = prod( dN, trans( mInv ) );
}

// =====================================================================================
//
// Class MITC4Params
//
// =====================================================================================

ShellThickElement3D4N::MITC4Params::MITC4Params(const ShellQ4_LocalCoordinateSystem & LCS)
    : Transformation(2, 2)
    , ShearStrains(4, 24, 0.0)
{

    double x21 = LCS.X2() - LCS.X1();
    double y21 = LCS.Y2() - LCS.Y1();
    double x34 = LCS.X3() - LCS.X4();
    double y34 = LCS.Y3() - LCS.Y4();
    double x41 = LCS.X4() - LCS.X1();
    double y41 = LCS.Y4() - LCS.Y1();
    double x32 = LCS.X3() - LCS.X2();
    double y32 = LCS.Y3() - LCS.Y2();

    Ax = - LCS.X1() + LCS.X2() + LCS.X3() - LCS.X4();
    Bx =   LCS.X1() - LCS.X2() + LCS.X3() - LCS.X4();
    Cx = - LCS.X1() - LCS.X2() + LCS.X3() + LCS.X4();
    Ay = - LCS.Y1() + LCS.Y2() + LCS.Y3() - LCS.Y4();
    By =   LCS.Y1() - LCS.Y2() + LCS.Y3() - LCS.Y4();
    Cy = - LCS.Y1() - LCS.Y2() + LCS.Y3() + LCS.Y4();

    double Alpha = std::atan( Ay / Ax );
    double Beta  = Globals::Pi * 0.5 - std::atan( Cx / Cy );

    Transformation(0, 0) =   std::sin(Beta);
    Transformation(0, 1) = - std::sin(Alpha);
    Transformation(1, 0) = - std::cos(Beta);
    Transformation(1, 1) =   std::cos(Alpha);

    ShearStrains(0,  2)  = - 0.5;
    ShearStrains(0,  3)  = - y41 * 0.25;
    ShearStrains(0,  4)  =   x41 * 0.25;

    ShearStrains(0, 20) =   0.5;
    ShearStrains(0, 21) = - y41 * 0.25;
    ShearStrains(0, 22) =   x41 * 0.25;

    ShearStrains(1,  2)  = - 0.5;
    ShearStrains(1,  3)  = - y21 * 0.25;
    ShearStrains(1,  4)  =   x21 * 0.25;

    ShearStrains(1,  8)  =   0.5;
    ShearStrains(1,  9) = - y21 * 0.25;
    ShearStrains(1, 10) =   x21 * 0.25;

    ShearStrains(2,  8)  = - 0.5;
    ShearStrains(2,  9) = - y32 * 0.25;
    ShearStrains(2, 10) =   x32 * 0.25;

    ShearStrains(2, 14) =   0.5;
    ShearStrains(2, 15) = - y32 * 0.25;
    ShearStrains(2, 16) =   x32 * 0.25;

    ShearStrains(3, 14) =   0.5;
    ShearStrains(3, 15) = - y34 * 0.25;
    ShearStrains(3, 16) =   x34 * 0.25;

    ShearStrains(3, 20) = - 0.5;
    ShearStrains(3, 21) = - y34 * 0.25;
    ShearStrains(3, 22) =   x34 * 0.25;
}

// =====================================================================================
//
// Class EASOperatorStorage
//
// =====================================================================================

ShellThickElement3D4N::EASOperatorStorage::EASOperatorStorage()
    : mInitialized(false)
{
}

void ShellThickElement3D4N::EASOperatorStorage::Initialize(const GeometryType& geom)
{
    if(!mInitialized)
    {
        noalias(alpha) = ZeroVector(5);
        noalias(alpha_converged) = ZeroVector(5);

        for(SizeType i = 0; i < 4; i++)
        {
            SizeType ii = i * 6;
            const array_1d<double, 3>& initialDispl = geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& initialRot = geom[i].FastGetSolutionStepValue(ROTATION);

            displ(ii    ) = initialDispl(0);
            displ(ii + 1) = initialDispl(1);
            displ(ii + 2) = initialDispl(2);
            displ_converged(ii    ) = initialDispl(0);
            displ_converged(ii + 1) = initialDispl(1);
            displ_converged(ii + 2) = initialDispl(2);

            displ(ii + 3) = initialRot(0);
            displ(ii + 4) = initialRot(1);
            displ(ii + 5) = initialRot(2);
            displ_converged(ii + 3) = initialRot(0);
            displ_converged(ii + 4) = initialRot(1);
            displ_converged(ii + 5) = initialRot(2);
        }

        mInitialized = true;
    }
}

void ShellThickElement3D4N::EASOperatorStorage::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    displ = displ_converged;
    alpha = alpha_converged;
}

void ShellThickElement3D4N::EASOperatorStorage::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    displ_converged = displ;
    alpha_converged = alpha;
}

void ShellThickElement3D4N::EASOperatorStorage::FinalizeNonLinearIteration(const Vector& displacementVector, ProcessInfo& CurrentProcessInfo)
{
    Vector incrementalDispl(24);
    noalias(incrementalDispl) = displacementVector - displ;
    noalias(displ) = displacementVector;

    array_1d<double, 5> temp;
    noalias(temp) = prod(L, incrementalDispl);
    noalias(temp) -= residual;
    noalias(alpha) -= prod(Hinv, temp);
}

void ShellThickElement3D4N::EASOperatorStorage::save(Serializer& rSerializer) const
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

void ShellThickElement3D4N::EASOperatorStorage::load(Serializer& rSerializer)
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

// =====================================================================================
//
// Class EASOperator
//
// =====================================================================================

ShellThickElement3D4N::EASOperator::EASOperator(const ShellQ4_LocalCoordinateSystem& LCS, EASOperatorStorage& storage)
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

    double j11 = Jac0(0,0);
    double j22 = Jac0(1,1);
    double j12 = Jac0(0,1);
    double j21 = Jac0(1,0);

    Matrix F0(3,3);
    F0(0,0) = j11*j11;
    F0(0,1) = j21*j12;
    F0(0,2) = 2.0*j11*j12;
    F0(1,0) = j12*j21;
    F0(1,1) = j22*j22;
    F0(1,2) = 2.0*j21*j22;
    F0(2,0) = j11*j21;
    F0(2,1) = j12*j22;
    F0(2,2) = j11*j22 + j12*j21;

    double dummyDet;
    MathUtils<double>::InvertMatrix3(F0, mF0inv, dummyDet);

    // initialize these data to zero because they will
    // be integrated during the gauss loop
    storage.L.clear();
    storage.Hinv.clear();
    storage.residual.clear();
}

void ShellThickElement3D4N::EASOperator::GaussPointComputation_Step1(double xi, double eta, const JacobianOperator& jac,
        Vector& generalizedStrains,
        EASOperatorStorage& storage)
{
    // construct the interpolation matrix in natural coordinate system
    Matrix E(3, 5, 0.0);
    E(0,0) = xi;
    E(1,1) = eta;
    E(2,2) = xi;
    E(2,3) = eta;
    E(0,4) = xi*eta;
    E(1,4) = -xi*eta;
    E(2,4) = xi*xi - eta*eta;

    // construct the interpolation matrix in local coordinate system
    noalias( mG ) = mJ0 / jac.Determinant() * prod( mF0inv, E );

    // assuming that the input generalized strains has been already calculated on this
    // integration point, we just need to add the contribution of the enhanced strains
    noalias( mEnhancedStrains ) = prod( mG, storage.alpha );

    generalizedStrains(0) += mEnhancedStrains(0); // e.xx
    generalizedStrains(1) += mEnhancedStrains(1); // e.yy
    generalizedStrains(2) += mEnhancedStrains(2); // e.xy
}

void ShellThickElement3D4N::EASOperator::GaussPointComputation_Step2(const Matrix& D,
        const Matrix& B,
        const Vector& S,
        EASOperatorStorage& storage)
{
    Matrix GTC(5, 3);
    noalias( GTC )                = prod( trans( mG ), project(D, range(0,3), range(0,3)) );
    noalias( storage.Hinv )      += prod( GTC, mG );
    noalias( storage.residual )  -= prod( trans( mG ), project(S, range(0,3)) );

    // compute L: [G'*Dmm, G'*Dmb, G'*Dms]*[Bm; Bm; Bs]
    int num_stress = D.size2(); // it can be 6 for thin shells or 8 for thick  shells
    Matrix GTD(5, num_stress, 0.0);
    noalias( project(GTD, range::all(), range(0,3)) ) = GTC;
    noalias( project(GTD, range::all(), range(3,6)) ) = prod( trans(mG), project(D, range(0,3), range(3,6)) );
    if(num_stress == 8)
        noalias( project(GTD, range::all(), range(6,8)) ) = prod( trans(mG), project(D, range(0,3), range(6,8)) );
    noalias( storage.L ) += prod( GTD, B );
}

void ShellThickElement3D4N::EASOperator::ComputeModfiedTangentAndResidual(Matrix& rLeftHandSideMatrix,
        Vector& rRightHandSideVector,
        EASOperatorStorage& storage)
{
    // invert H
    Matrix Hcopy( storage.Hinv );
    permutation_matrix<Matrix::size_type> pm( 5 );
    lu_factorize( Hcopy, pm );
    noalias( storage.Hinv ) = IdentityMatrix( 5 );
    lu_substitute( Hcopy, pm, storage.Hinv );

    // compute L' * H^-1
    Matrix LTHinv(24, 5);
    noalias( LTHinv ) = prod( trans( storage.L ), storage.Hinv );

    // modify the stiffness matrix and the residual vector
    // for the static condensation of the enhanced strain parameters
    noalias( rRightHandSideVector ) -= prod( LTHinv, storage.residual );   // R_mod = R - L' * H^-1 * residual
    noalias( rLeftHandSideMatrix  ) -= prod( LTHinv, storage.L );          // K_mod = K - L' * H^-1 * L
}

// =====================================================================================
//
// Class ShellThickElement3D4N
//
// =====================================================================================

ShellThickElement3D4N::ShellThickElement3D4N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        bool NLGeom)
    : BaseShellElement(NewId, pGeometry)
    , mpCoordinateTransformation( NLGeom ?
                                  new ShellQ4_CorotationalCoordinateTransformation(pGeometry) :
                                  new ShellQ4_CoordinateTransformation(pGeometry))
{
}

ShellThickElement3D4N::ShellThickElement3D4N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        bool NLGeom)
    : BaseShellElement(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation( NLGeom ?
                                  new ShellQ4_CorotationalCoordinateTransformation(pGeometry) :
                                  new ShellQ4_CoordinateTransformation(pGeometry))
{
}

ShellThickElement3D4N::ShellThickElement3D4N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        CoordinateTransformationBasePointerType pCoordinateTransformation)
    : BaseShellElement(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation(pCoordinateTransformation)
{
}

ShellThickElement3D4N::~ShellThickElement3D4N()
{
}

// TODO are the GetIntegrationMethod methods needed (implemented in the other 3 shells)

Element::Pointer ShellThickElement3D4N::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );

    return Kratos::make_shared< ShellThickElement3D4N >(NewId, newGeom, pProperties, mpCoordinateTransformation->Create(newGeom) );
}

Element::Pointer ShellThickElement3D4N::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared< ShellThickElement3D4N >(NewId, pGeom, pProperties, mpCoordinateTransformation->Create(pGeom) );
}

void ShellThickElement3D4N::Initialize()
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry();

    KRATOS_ERROR_IF_NOT((r_geom.IntegrationPoints(GetIntegrationMethod())).size() == 4)
        << "ShellThickElement3D4N - needs a full integration scheme" << std::endl;

    const int points_number = r_geom.PointsNumber();
    KRATOS_ERROR_IF_NOT(points_number == 4) << "ShellThickElement3D4N - Wrong number of nodes"
        << points_number << std::endl;

    BaseShellElement::Initialize();

    mpCoordinateTransformation->Initialize();

    this->SetupOrientationAngles();

    mEASStorage.Initialize(r_geom);

    KRATOS_CATCH("")
}

void ShellThickElement3D4N::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->InitializeNonLinearIteration(rCurrentProcessInfo);

    BaseInitializeNonLinearIteration(rCurrentProcessInfo);
}

void ShellThickElement3D4N::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->FinalizeNonLinearIteration(rCurrentProcessInfo);

    ShellQ4_LocalCoordinateSystem LCS( mpCoordinateTransformation->CreateLocalCoordinateSystem() );
    Vector globalDisplacementVector(24);
    GetValuesVector(globalDisplacementVector);
    Vector localDisplacementVector( mpCoordinateTransformation->CalculateLocalDisplacements( LCS, globalDisplacementVector ) );

    mEASStorage.FinalizeNonLinearIteration(localDisplacementVector, rCurrentProcessInfo);

    BaseFinalizeNonLinearIteration(rCurrentProcessInfo);
}

void ShellThickElement3D4N::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    BaseInitializeSolutionStep(rCurrentProcessInfo);

    mpCoordinateTransformation->InitializeSolutionStep(rCurrentProcessInfo);

    mEASStorage.InitializeSolutionStep(rCurrentProcessInfo);
}

void ShellThickElement3D4N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    BaseFinalizeSolutionStep(rCurrentProcessInfo);

    mpCoordinateTransformation->FinalizeSolutionStep(rCurrentProcessInfo);

    mEASStorage.FinalizeSolutionStep(rCurrentProcessInfo);
}

void ShellThickElement3D4N::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if((rMassMatrix.size1() != 24) || (rMassMatrix.size2() != 24))
        rMassMatrix.resize(24, 24, false);
    noalias(rMassMatrix) = ZeroMatrix(24, 24);

    // Compute the local coordinate system.

    ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem() );

    // lumped area

    double lump_area = referenceCoordinateSystem.Area() / 4.0;

    // Calculate avarage mass per unit area
    double av_mass_per_unit_area = 0.0;
    for(SizeType i = 0; i < 4; i++)
        av_mass_per_unit_area += mSections[i]->CalculateMassPerUnitArea(GetProperties());
    av_mass_per_unit_area /= 4.0;

    // Gauss Loop

    for(SizeType i = 0; i < 4; i++)
    {
        SizeType index = i * 6;

        double nodal_mass = av_mass_per_unit_area * lump_area;

        // translational mass
        rMassMatrix(index, index)            = nodal_mass;
        rMassMatrix(index + 1, index + 1)    = nodal_mass;
        rMassMatrix(index + 2, index + 2)    = nodal_mass;

        // rotational mass - neglected for the moment...
    }
}

// =====================================================================================
//
// Class ShellThickElement3D4N - Results on Gauss Points
//
// =====================================================================================

void ShellThickElement3D4N::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    SizeType size = GetGeometry().size();
    if (rValues.size() != size)
        rValues.resize(size);

    // The membrane formulation needs to iterate to find the correct
    // mid-surface strain values.
    // Check if we are doing a non-linear analysis type. If not, print warning

    if (!rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
    {
        KRATOS_WARNING_ONCE("ShellThickElement3D4N") << "Warning: Gauss point results have "
            << "been requested for a linear analysis.\nThe membrane formulation used "
            << "requires iteration to accurately determine recovered "
            << "quantities (strain, stress, etc...)." << std::endl;
    }

	if (rVariable == VON_MISES_STRESS ||
		rVariable == VON_MISES_STRESS_TOP_SURFACE ||
		rVariable == VON_MISES_STRESS_MIDDLE_SURFACE ||
		rVariable == VON_MISES_STRESS_BOTTOM_SURFACE)
	{
		// Von mises calcs

		// Get some references.
		const PropertiesType & props = GetProperties();
		const GeometryType & geom = GetGeometry();
		const Matrix & shapeFunctions = geom.ShapeFunctionsValues();
		Vector iN(shapeFunctions.size2());

		// Compute the local coordinate system.
		ShellQ4_LocalCoordinateSystem localCoordinateSystem(
			mpCoordinateTransformation->CreateLocalCoordinateSystem());

		ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
			mpCoordinateTransformation->CreateReferenceCoordinateSystem());

		// Prepare all the parameters needed for the MITC formulation.
		// This is to be done here outside the Gauss Loop.
		MITC4Params shearParameters(referenceCoordinateSystem);

		// Instantiate the Jacobian Operator.
		// This will store:
		// the jacobian matrix, its inverse, its determinant
		// and the derivatives of the shape functions in the local
		// coordinate system
		JacobianOperator jacOp;

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
		GetValuesVector(globalDisplacements, 0);

		// Get the current displacements in local coordinate system
		Vector localDisplacements(
			mpCoordinateTransformation->CalculateLocalDisplacements(localCoordinateSystem, globalDisplacements));

		// Instantiate the EAS Operator.
		// This will apply the Enhanced Assumed Strain Method for the calculation
		// of the membrane contribution.
		EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

		// Just to store the rotation matrix for visualization purposes
		Matrix R(8, 8);
		Matrix aux33(3, 3);

		// Initialize parameters for the cross section calculation
		ShellCrossSection::SectionParameters parameters(geom, props, rCurrentProcessInfo);
		parameters.SetGeneralizedStrainVector(generalizedStrains);
		parameters.SetGeneralizedStressVector(generalizedStresses);
		parameters.SetConstitutiveMatrix(D);
		Flags& options = parameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

		// Gauss Loop
		for (unsigned int i = 0; i < size; i++)
		{
			// get a reference of the current integration point and shape functions
			const GeometryType::IntegrationPointType & ip = geom.IntegrationPoints()[i];

			noalias(iN) = row(shapeFunctions, i);

			// Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
			// and Shape functions derivatives in the local coordinate system
			jacOp.Calculate(referenceCoordinateSystem, geom.ShapeFunctionLocalGradient(i));

			// Compute all strain-displacement matrices
			CalculateBMatrix(ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling);

			// Calculate strain vectors in local coordinate system
			noalias(generalizedStrains) = prod(B, localDisplacements);

			// Apply the EAS method to modify the membrane part of the strains computed above.
			EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains, mEASStorage);

			// Calculate the response of the Cross Section
			ShellCrossSection::Pointer & section = mSections[i];

			//add in shear stabilization
			double shearStabilisation = CalculateStenbergShearStabilization(referenceCoordinateSystem, section->GetThickness(GetProperties()));
			parameters.SetStenbergShearStabilization(shearStabilisation);
			//double shearStabilisation = (hMean*hMean) / (hMean*hMean + 0.1*h_e*h_e);

			// calculate force resultants
			parameters.SetShapeFunctionsValues(iN);
			parameters.SetShapeFunctionsDerivatives(jacOp.XYDerivatives());
			section->CalculateSectionResponse(parameters, ConstitutiveLaw::StressMeasure_PK2);

			// Compute stresses
			CalculateStressesFromForceResultants(generalizedStresses,
				section->GetThickness(GetProperties()));

			// Calculate von mises results
			CalculateVonMisesStress(generalizedStresses, rVariable, rValues[i]);

		} // end gauss loop
	}
	else if (rVariable == TSAI_WU_RESERVE_FACTOR)
	{
		// resize output
		SizeType size = 4;
		if (rValues.size() != size)
			rValues.resize(size);

		// Get some references.
		const PropertiesType & props = GetProperties();
		const GeometryType & geom = GetGeometry();
		const Matrix & shapeFunctions = geom.ShapeFunctionsValues();
		Vector iN(shapeFunctions.size2());

		// Compute the local coordinate system.
		ShellQ4_LocalCoordinateSystem localCoordinateSystem(
			mpCoordinateTransformation->CreateLocalCoordinateSystem());
		ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
			mpCoordinateTransformation->CreateReferenceCoordinateSystem());

		// Prepare all the parameters needed for the MITC formulation.
		// This is to be done here outside the Gauss Loop.
		MITC4Params shearParameters(referenceCoordinateSystem);

		// Instantiate the Jacobian Operator.
		JacobianOperator jacOp;

		// Instantiate all strain-displacement matrices.
		Matrix B(8, 24, 0.0);
		Vector Bdrilling(24, 0.0);

		// Instantiate all section tangent matrices.
		Matrix D(8, 8, 0.0);

		// Instantiate strain and stress-resultant vectors
		Vector generalizedStrains(8);
		Vector generalizedStresses(8);
		std::vector<VectorType> rlaminateStrains;
		std::vector<VectorType> rlaminateStresses;

		// Get the current displacements in global coordinate system
		Vector globalDisplacements(24);
		GetValuesVector(globalDisplacements, 0);

		// Get the current displacements in local coordinate system
		Vector localDisplacements(
			mpCoordinateTransformation->CalculateLocalDisplacements(localCoordinateSystem, globalDisplacements));

		// Instantiate the EAS Operator.
		EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

		// Initialize parameters for the cross section calculation
		ShellCrossSection::SectionParameters parameters(geom, props, rCurrentProcessInfo);
		parameters.SetGeneralizedStrainVector(generalizedStrains);
		parameters.SetGeneralizedStressVector(generalizedStresses);
		parameters.SetConstitutiveMatrix(D);
		Flags& options = parameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

		// Get all laminae strengths
		ShellCrossSection::Pointer & section = mSections[0];
		std::vector<Matrix> Laminae_Strengths =
			std::vector<Matrix>(section->NumberOfPlies());
		for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
		{
			Laminae_Strengths[ply].resize(3, 3, 0.0);
			Laminae_Strengths[ply].clear();
		}
		section->GetLaminaeStrengths(Laminae_Strengths, props);

		// Define variables
		Matrix R(8, 8);
		double total_rotation = 0.0;

		// Gauss Loop
		for (unsigned int i = 0; i < size; i++)
		{
			// get a reference of the current integration point and shape functions
			const GeometryType::IntegrationPointType & ip = geom.IntegrationPoints()[i];
			noalias(iN) = row(shapeFunctions, i);

			// Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
			// and Shape functions derivatives in the local coordinate system
			jacOp.Calculate(referenceCoordinateSystem, geom.ShapeFunctionLocalGradient(i));

			// Compute all strain-displacement matrices
			CalculateBMatrix(ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling);

			// Calculate strain vectors in local coordinate system
			noalias(generalizedStrains) = prod(B, localDisplacements);

			// Apply the EAS method to modify the membrane part of the strains computed above.
			EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains, mEASStorage);

			// Calculate the response of the Cross Section
			ShellCrossSection::Pointer & section = mSections[i];

			//Calculate lamina stresses
			CalculateLaminaStrains(section, generalizedStrains, rlaminateStrains);
			CalculateLaminaStresses(section, parameters, rlaminateStrains, rlaminateStresses);

			// Retrieve ply orientations
			section = mSections[i];
			Vector ply_orientation(section->NumberOfPlies());
			section->GetLaminaeOrientation(GetProperties(), ply_orientation);

			// Rotate lamina stress from element CS to section CS, and then
			// to lamina angle to lamina material principal directions
			for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
			{
				total_rotation = -ply_orientation[ply] - (section->GetOrientationAngle());
				section->GetRotationMatrixForGeneralizedStresses(total_rotation, R);
				//top surface of current ply
				rlaminateStresses[2 * ply] = prod(R, rlaminateStresses[2 * ply]);
				//bottom surface of current ply
				rlaminateStresses[2 * ply + 1] = prod(R, rlaminateStresses[2 * ply + 1]);
			}

			// Calculate Tsai-Wu criterion for each ply, take min of all plies
			double min_tsai_wu = 0.0;
			double temp_tsai_wu = 0.0;
			for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
			{
				temp_tsai_wu = CalculateTsaiWuPlaneStress(rlaminateStresses, Laminae_Strengths[ply], ply);
				if (ply == 0)
				{
					min_tsai_wu = temp_tsai_wu;
				}
				else if (temp_tsai_wu < min_tsai_wu)
				{
					min_tsai_wu = temp_tsai_wu;
				}
			}

			// Output min Tsai-Wu result
			rValues[i] = min_tsai_wu;

		}//Gauss loop
	}
	else
	{
		std::vector<double> temp(size);

		for (SizeType i = 0; i < size; i++)
			mSections[i]->GetValue(rVariable, GetProperties(), temp[i]);

		const Matrix & shapeFunctions = GetGeometry().ShapeFunctionsValues();
		Vector N(size);

		for (SizeType i = 0; i < size; i++)
		{
			noalias(N) = row(shapeFunctions, i);
			double& ival = rValues[i];
			ival = 0.0;
			for (SizeType j = 0; j < size; j++)
			{
				ival += N(j) * temp[j];
			}
		}
	}
}

void ShellThickElement3D4N::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    // The membrane formulation needs to iterate to find the correct
    // mid-surface strain values.
    // Check if we are doing a non-linear analysis type. If not, print warning

    if (!rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
    {
        KRATOS_WARNING_ONCE("ShellThickElement3D4N") << "Warning: Gauss point results have "
            << "been requested for a linear analysis.\nThe membrane formulation used "
            << "requires iteration to accurately determine recovered "
            << "quantities (strain, stress, etc...)." << std::endl;
    }


    if(TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(rVariable, rValues, rCurrentProcessInfo)) return;
}

void ShellThickElement3D4N::CalculateOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
        std::vector<array_1d<double,3> >& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == LOCAL_AXIS_1 ||
        rVariable == LOCAL_AXIS_2 ||
        rVariable == LOCAL_AXIS_3) {
        BaseShellElement::ComputeLocalAxis(rVariable, rOutput, mpCoordinateTransformation);
    }
    else if (rVariable == LOCAL_MATERIAL_AXIS_1 ||
             rVariable == LOCAL_MATERIAL_AXIS_2 ||
             rVariable == LOCAL_MATERIAL_AXIS_3) {
        BaseShellElement::ComputeLocalMaterialAxis(rVariable, rOutput, mpCoordinateTransformation);
    }
}

void ShellThickElement3D4N::Calculate(const Variable<Matrix>& rVariable, Matrix & Output, const ProcessInfo & rCurrentProcessInfo)
{
	if (rVariable == LOCAL_ELEMENT_ORIENTATION)
	{
		Output.resize(3, 3, false);

		// Compute the local coordinate system.
		ShellQ4_LocalCoordinateSystem localCoordinateSystem(
			mpCoordinateTransformation->CreateReferenceCoordinateSystem());
		Output = localCoordinateSystem.Orientation();
	}
}

// =====================================================================================
//
// Class ShellThickElement3D4N - Private methods
//
// =====================================================================================

void ShellThickElement3D4N::CalculateStressesFromForceResultants(VectorType & rstresses, const double & rthickness)
{
	// Refer http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch20.d/AFEM.Ch20.pdf

	// membrane forces -> in-plane stresses (av. across whole thickness)
	rstresses[0] /= rthickness;
	rstresses[1] /= rthickness;
	rstresses[2] /= rthickness;

	// bending moments -> peak in-plane stresses (@ top and bottom surface)
	rstresses[3] *= 6.0 / (rthickness*rthickness);
	rstresses[4] *= 6.0 / (rthickness*rthickness);
	rstresses[5] *= 6.0 / (rthickness*rthickness);

	// shear forces -> peak shear stresses (@ midsurface)
	rstresses[6] *= 1.5 / rthickness;
	rstresses[7] *= 1.5 / rthickness;
}

void ShellThickElement3D4N::CalculateLaminaStrains(ShellCrossSection::Pointer & section, const Vector & generalizedStrains, std::vector<VectorType>& rlaminateStrains)
{
	// Get laminate properties
	double thickness = section->GetThickness(GetProperties());
	double z_current = thickness / -2.0; // start from the top of the 1st layer

	// Establish current strains at the midplane
	// (element coordinate system)
	double e_x = generalizedStrains[0];
	double e_y = generalizedStrains[1];
	double e_xy = generalizedStrains[2];	//this is still engineering
											//strain (2xtensorial shear)
	double kap_x = generalizedStrains[3];
	double kap_y = generalizedStrains[4];
	double kap_xy = generalizedStrains[5];	//this is still engineering

	// Get ply thicknesses
	Vector ply_thicknesses = Vector(section->NumberOfPlies(), 0.0);
	section->GetPlyThicknesses(GetProperties(), ply_thicknesses);

	// Resize output vector. 2 Surfaces for each ply
	rlaminateStrains.resize(2 * section->NumberOfPlies());
	for (unsigned int i = 0; i < 2 * section->NumberOfPlies(); i++)
	{
		rlaminateStrains[i].resize(8, false);
		rlaminateStrains[i].clear();
	}

	// Loop over all plies - start from bottom ply, bottom surface
	for (unsigned int plyNumber = 0;
		plyNumber < section->NumberOfPlies(); ++plyNumber)
	{
		// Calculate strains at top surface, arranged in columns.
		// (element coordinate system)
		rlaminateStrains[2 * plyNumber][0] = e_x + z_current*kap_x;
		rlaminateStrains[2 * plyNumber][1] = e_y + z_current*kap_y;
		rlaminateStrains[2 * plyNumber][2] = e_xy + z_current*kap_xy;

		// constant transverse shear strain dist
		rlaminateStrains[2 * plyNumber][6] = generalizedStrains[6];
		rlaminateStrains[2 * plyNumber][7] = generalizedStrains[7];

		// Move to bottom surface of current layer
		z_current += ply_thicknesses[plyNumber];

		// Calculate strains at bottom surface, arranged in columns
		// (element coordinate system)
		rlaminateStrains[2 * plyNumber + 1][0] = e_x + z_current*kap_x;
		rlaminateStrains[2 * plyNumber + 1][1] = e_y + z_current*kap_y;
		rlaminateStrains[2 * plyNumber + 1][2] = e_xy + z_current*kap_xy;

		// constant transverse shear strain dist
		rlaminateStrains[2 * plyNumber + 1][6] = generalizedStrains[6];
		rlaminateStrains[2 * plyNumber + 1][7] = generalizedStrains[7];
	}
}

void ShellThickElement3D4N::CalculateLaminaStresses(ShellCrossSection::Pointer & section, ShellCrossSection::SectionParameters parameters, const std::vector<VectorType>& rlaminateStrains, std::vector<VectorType>& rlaminateStresses)
{
	// Setup flag to compute ply constitutive matrices
	// (units [Pa] and rotated to element orientation)
	section->SetupGetPlyConstitutiveMatrices();
	Flags& options = parameters.GetOptions();
	options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,true);
	section->CalculateSectionResponse(parameters,
		ConstitutiveLaw::StressMeasure_PK2);

	// Resize output vector. 2 Surfaces for each ply
	rlaminateStresses.resize(2 * section->NumberOfPlies());
	for (unsigned int i = 0; i < 2 * section->NumberOfPlies(); i++)
	{
		rlaminateStresses[i].resize(8, false);
		rlaminateStresses[i].clear();
	}

	// Loop over all plies - start from top ply, top surface
	for (unsigned int plyNumber = 0;
		plyNumber < section->NumberOfPlies(); ++plyNumber)
	{
		// determine stresses at currrent ply, top surface
		// (element coordinate system)
		rlaminateStresses[2 * plyNumber] = prod(
			section->GetPlyConstitutiveMatrix(plyNumber),
			rlaminateStrains[2 * plyNumber]);

		// determine stresses at currrent ply, bottom surface
		// (element coordinate system)
		rlaminateStresses[2 * plyNumber + 1] = prod(
			section->GetPlyConstitutiveMatrix(plyNumber),
			rlaminateStrains[2 * plyNumber + 1]);
	}
}

double ShellThickElement3D4N::CalculateTsaiWuPlaneStress(const std::vector<VectorType>& rlaminateStresses, const Matrix & rLamina_Strengths, const unsigned int & rPly)
{
	// Incoming lamina strengths are organized as follows:
	// Refer to 'shell_cross_section.cpp' for details.
	//
	//	|	T1,		C1,		T2	|
	//	|	C2,		S12,	S13	|
	//	|   S23		0		0	|

	// Convert raw lamina strengths into tsai strengths F_i and F_ij.
	// Refer Reddy (2003) Section 10.9.4 (re-ordered for kratos DOFs).
	// All F_i3 components ignored - thin shell theory.
	//

	// Controls F_12
	// Should be FALSE unless testing against other programs that ignore it.
	bool disable_in_plane_interaction = false;

	// First, F_i
	Vector F_i = Vector(3, 0.0);
	F_i[0] = 1.0 / rLamina_Strengths(0, 0) - 1.0 / rLamina_Strengths(0, 1);
	F_i[1] = 1.0 / rLamina_Strengths(0, 2) - 1.0 / rLamina_Strengths(1, 0);
	F_i[2] = 0.0;

	// Second, F_ij
	Matrix F_ij = Matrix(5, 5, 0.0);
	F_ij.clear();
	F_ij(0, 0) = 1.0 / rLamina_Strengths(0, 0) / rLamina_Strengths(0, 1);	// 11
	F_ij(1, 1) = 1.0 / rLamina_Strengths(0, 2) / rLamina_Strengths(1, 0);	// 22
	F_ij(2, 2) = 1.0 / rLamina_Strengths(1, 1) / rLamina_Strengths(1, 1);	// 12
	F_ij(0, 1) = F_ij(1, 0) = -0.5 / std::sqrt(rLamina_Strengths(0, 0)*rLamina_Strengths(0, 1)*rLamina_Strengths(0, 2)*rLamina_Strengths(1, 0));

	if (disable_in_plane_interaction)
	{
		F_ij(0, 1) = F_ij(1, 0) = 0.0;
	}

	// Third, addditional transverse shear terms
	F_ij(3, 3) = 1.0 / rLamina_Strengths(1, 2) / rLamina_Strengths(1, 2);	// 13
	F_ij(4, 4) = 1.0 / rLamina_Strengths(2, 0) / rLamina_Strengths(2, 0);	// 23

																			// Evaluate Tsai-Wu @ top surface of current layer
	double var_a = 0.0;
	double var_b = 0.0;
	for (SizeType i = 0; i < 3; i++)
	{
		var_b += F_i[i] * rlaminateStresses[2 * rPly][i];
		for (SizeType j = 0; j < 3; j++)
		{
			var_a += F_ij(i, j)*rlaminateStresses[2 * rPly][i] * rlaminateStresses[2 * rPly][j];
		}
	}
	var_a += F_ij(3, 3)*rlaminateStresses[2 * rPly][6] * rlaminateStresses[2 * rPly][6]; // Transverse shear 13
	var_a += F_ij(4, 4)*rlaminateStresses[2 * rPly][7] * rlaminateStresses[2 * rPly][7]; // Transverse shear 23

	double tsai_reserve_factor_top = (-1.0*var_b + std::sqrt(var_b*var_b + 4.0 * var_a)) / 2.0 / var_a;

	// Evaluate Tsai-Wu @ bottom surface of current layer
	var_a = 0.0;
	var_b = 0.0;
	for (SizeType i = 0; i < 3; i++)
	{
		var_b += F_i[i] * rlaminateStresses[2 * rPly + 1][i];
		for (SizeType j = 0; j < 3; j++)
		{
			var_a += F_ij(i, j)*rlaminateStresses[2 * rPly + 1][i] * rlaminateStresses[2 * rPly + 1][j];
		}
	}
	var_a += F_ij(3, 3)*rlaminateStresses[2 * rPly + 1][6] * rlaminateStresses[2 * rPly + 1][6]; // Transverse shear 13
	var_a += F_ij(4, 4)*rlaminateStresses[2 * rPly + 1][7] * rlaminateStresses[2 * rPly + 1][7]; // Transverse shear 23

	double tsai_reserve_factor_bottom = (-1.0*var_b + std::sqrt(var_b*var_b + 4.0 * var_a)) / 2.0 / var_a;

	// Return min of both surfaces as the result for the whole ply
	return std::min(tsai_reserve_factor_bottom, tsai_reserve_factor_top);
}

void ShellThickElement3D4N::CalculateVonMisesStress(const Vector & generalizedStresses, const Variable<double>& rVariable, double & rVon_Mises_Result)
{
	// calc von mises stresses at top, mid and bottom surfaces for
	// thick shell
	double sxx, syy, sxy, sxz, syz;
	double von_mises_top, von_mises_mid, von_mises_bottom;
	// top surface: membrane and +bending contributions
	//				(no transverse shear)
	sxx = generalizedStresses[0] + generalizedStresses[3];
	syy = generalizedStresses[1] + generalizedStresses[4];
	sxy = generalizedStresses[2] + generalizedStresses[5];
	von_mises_top = sxx*sxx - sxx*syy + syy*syy + 3.0*sxy*sxy;

	// mid surface: membrane and transverse shear contributions
	//				(no bending)
	sxx = generalizedStresses[0];
	syy = generalizedStresses[1];
	sxy = generalizedStresses[2];
	sxz = generalizedStresses[6];
	syz = generalizedStresses[7];
	von_mises_mid = sxx*sxx - sxx*syy + syy*syy +
		3.0*(sxy*sxy + sxz*sxz + syz*syz);

	// bottom surface:	membrane and -bending contributions
	//					(no transverse shear)
	sxx = generalizedStresses[0] - generalizedStresses[3];
	syy = generalizedStresses[1] - generalizedStresses[4];
	sxy = generalizedStresses[2] - generalizedStresses[5];
	von_mises_bottom = sxx*sxx - sxx*syy + syy*syy + 3.0*sxy*sxy;

	// Output requested quantity
	if (rVariable == VON_MISES_STRESS_TOP_SURFACE)
	{
		rVon_Mises_Result = std::sqrt(von_mises_top);
	}
	else if (rVariable == VON_MISES_STRESS_MIDDLE_SURFACE)
	{
		rVon_Mises_Result = std::sqrt(von_mises_mid);
	}
	else if (rVariable == VON_MISES_STRESS_BOTTOM_SURFACE)
	{
		rVon_Mises_Result = std::sqrt(von_mises_bottom);
	}
	else if (rVariable == VON_MISES_STRESS)
	{
		// take the greatest value and output
		rVon_Mises_Result =
        std::sqrt(std::max(von_mises_top,
				std::max(von_mises_mid, von_mises_bottom)));
	}
}

void ShellThickElement3D4N::CheckGeneralizedStressOrStrainOutput(const Variable<Matrix>& rVariable, int & ijob, bool & bGlobal)
{
	if (rVariable == SHELL_STRAIN)
	{
		ijob = 1;
	}
	else if (rVariable == SHELL_STRAIN_GLOBAL)
	{
		ijob = 1;
		bGlobal = true;
	}
	else if (rVariable == SHELL_CURVATURE)
	{
		ijob = 2;
	}
	else if (rVariable == SHELL_CURVATURE_GLOBAL)
	{
		ijob = 2;
		bGlobal = true;
	}
	else if (rVariable == SHELL_FORCE)
	{
		ijob = 3;
	}
	else if (rVariable == SHELL_FORCE_GLOBAL)
	{
		ijob = 3;
		bGlobal = true;
	}
	else if (rVariable == SHELL_MOMENT)
	{
		ijob = 4;
	}
	else if (rVariable == SHELL_MOMENT_GLOBAL)
	{
		ijob = 4;
		bGlobal = true;
	}
	else if (rVariable == SHELL_STRESS_TOP_SURFACE)
	{
		ijob = 5;
	}
	else if (rVariable == SHELL_STRESS_TOP_SURFACE_GLOBAL)
	{
		ijob = 5;
		bGlobal = true;
	}
	else if (rVariable == SHELL_STRESS_MIDDLE_SURFACE)
	{
		ijob = 6;
	}
	else if (rVariable == SHELL_STRESS_MIDDLE_SURFACE_GLOBAL)
	{
		ijob = 6;
		bGlobal = true;
	}
	else if (rVariable == SHELL_STRESS_BOTTOM_SURFACE)
	{
		ijob = 7;
	}
	else if (rVariable == SHELL_STRESS_BOTTOM_SURFACE_GLOBAL)
	{
		ijob = 7;
		bGlobal = true;
	}
	else if (rVariable == SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE)
	{
		ijob = 8;
	}
	else if (rVariable == SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE_GLOBAL)
	{
		ijob = 8;
		bGlobal = true;
	}
	else if (rVariable == SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE)
	{
		ijob = 9;
	}
	else if (rVariable == SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE_GLOBAL)
	{
		ijob = 9;
		bGlobal = true;
	}
}

double ShellThickElement3D4N::CalculateStenbergShearStabilization(const ShellQ4_LocalCoordinateSystem & referenceCoordinateSystem, const double & hMean)
{
	// Calculate Stenberg shear stabilisation as per
	// https://doi.org/10.1016/j.cma.2003.12.036 section 3.1

	// Determine longest element edge
	Vector edge_1 = Vector(referenceCoordinateSystem.P1() - referenceCoordinateSystem.P2());
	Vector edge_2 = Vector(referenceCoordinateSystem.P2() - referenceCoordinateSystem.P3());
	Vector edge_3 = Vector(referenceCoordinateSystem.P3() - referenceCoordinateSystem.P4());
	Vector edge_4 = Vector(referenceCoordinateSystem.P4() - referenceCoordinateSystem.P1());
	double h_e = inner_prod(edge_1, edge_1);
	if (inner_prod(edge_2, edge_2) > h_e)
	{
		h_e = inner_prod(edge_2, edge_2);
	}
	if (inner_prod(edge_3, edge_3) > h_e)
	{
		h_e = inner_prod(edge_3, edge_3);
	}
	if (inner_prod(edge_4, edge_4) > h_e)
	{
		h_e = inner_prod(edge_4, edge_4);
	}
	h_e = std::sqrt(h_e);

	return ((hMean*hMean) / (hMean*hMean + 0.1*h_e*h_e));
}

void ShellThickElement3D4N::DecimalCorrection(Vector& a)
{
    double norm = norm_2(a);
    double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
    for(SizeType i = 0; i < a.size(); i++)
        if(std::abs(a(i)) < tolerance)
            a(i) = 0.0;
}

void ShellThickElement3D4N::SetupOrientationAngles()
{
    if (this->Has(MATERIAL_ORIENTATION_ANGLE))
    {
        for (CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it)
        (*it)->SetOrientationAngle(this->GetValue(MATERIAL_ORIENTATION_ANGLE));
    }
    else
    {
        ShellQ4_LocalCoordinateSystem lcs( mpCoordinateTransformation->CreateReferenceCoordinateSystem() );

        Vector3Type normal;
        noalias( normal ) = lcs.Vz();

        Vector3Type dZ;
        dZ(0) = 0.0;
        dZ(1) = 0.0;
        dZ(2) = 1.0; // for the moment let's take this. But the user can specify its own triad! TODO

        Vector3Type dirX;
        MathUtils<double>::CrossProduct(dirX,   dZ, normal);

        // try to normalize the x vector. if it is near zero it means that we need
        // to choose a default one.
        double dirX_norm = dirX(0)*dirX(0) + dirX(1)*dirX(1) + dirX(2)*dirX(2);
        if(dirX_norm < 1.0E-12)
        {
            dirX(0) = 1.0;
            dirX(1) = 0.0;
            dirX(2) = 0.0;
        }
        else if(dirX_norm != 1.0)
        {
            dirX_norm = std::sqrt(dirX_norm);
            dirX /= dirX_norm;
        }

        Vector3Type elem_dirX = lcs.Vx();

        // now calculate the angle between the element x direction and the material x direction.
        Vector3Type& a = elem_dirX;
        Vector3Type& b = dirX;
        double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
        if(a_dot_b < -1.0) a_dot_b = -1.0;
        if(a_dot_b >  1.0) a_dot_b =  1.0;
        double angle = std::acos( a_dot_b );

        // if they are not counter-clock-wise, let's change the sign of the angle
        if(angle != 0.0)
        {
            const MatrixType& R = lcs.Orientation();
            if( dirX(0)*R(1, 0) + dirX(1)*R(1, 1) + dirX(2)*R(1, 2) < 0.0 )
                angle = -angle;
        }

        for(CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it)
            (*it)->SetOrientationAngle(angle);
    }
}

void ShellThickElement3D4N::CalculateBMatrix(double xi, double eta,
        const JacobianOperator& Jac, const MITC4Params& mitc_params,
        const Vector& N,
        Matrix& B, Vector& Bdrill)
{
    // some data

    const Matrix& dNxy = Jac.XYDerivatives();

    // membrane ************************************************************************************************

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

    // bending *************************************************************************************************

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

    // drilling ************************************************************************************************

    Bdrill( 0) = -0.5 * dNxy(0, 1);
    Bdrill( 1) =  0.5 * dNxy(0, 0);
    Bdrill( 5) = -N(0);
    Bdrill( 6) = -0.5 * dNxy(1, 1);
    Bdrill( 7) =  0.5 * dNxy(1, 0);
    Bdrill(11) = -N(1);
    Bdrill(12) = -0.5 * dNxy(2, 1);
    Bdrill(13) =  0.5 * dNxy(2, 0);
    Bdrill(17) = -N(2);
    Bdrill(18) = -0.5 * dNxy(3, 1);
    Bdrill(19) =  0.5 * dNxy(3, 0);
    Bdrill(23) = -N(3);

    // shear ***************************************************************************************************

    // MITC modified shape functions
    Matrix MITCShapeFunctions(2, 4, 0.0);
    MITCShapeFunctions(1, 0) = 1.0 - xi;
    MITCShapeFunctions(0, 1) = 1.0 - eta;
    MITCShapeFunctions(1, 2) = 1.0 + xi;
    MITCShapeFunctions(0, 3) = 1.0 + eta;

    // strain displacement matrix in natural coordinate system.
    // interpolate the shear strains given in MITC4Params
    // using the modified shape function
    Matrix BN = prod( MITCShapeFunctions, mitc_params.ShearStrains );

    // modify the shear strain intensity in the tying points
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

    row( BN, 0 ) *= Temp1;
    row( BN, 1 ) *= Temp2;

    // transform the strain-displacement matrix from natural
    // to local coordinate system taking into account the element distorsion

    project( B, slice(7, -1, 2), slice::all() ) = prod( mitc_params.Transformation, BN );

    // Explanation of the 'slice':
    // The MITC4Params class and the first part of this method were coded
    // assuming the following notations for strain vector : [e.xx, e.yy, e.zz, 2 e.xy, 2 e.XZ, 2 e.YZ].
    // Since in Kratos the strain vector is assumed to be : [e.xx, e.yy, e.zz, 2 e.xy, 2 e.YZ, 2 e.XZ],
    // I had to "virtually" swap the 2 rows of the matrix B before the assignment from the product on the Right-Hand-Side.
    // In this way the matrix B is consistent with the notation used in Kratos without recoding all the MITC stuff.
}

void ShellThickElement3D4N::CalculateAll(MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    // Resize the Left Hand Side if necessary,
    // and initialize it to Zero

    if((rLeftHandSideMatrix.size1() != 24) || (rLeftHandSideMatrix.size2() != 24))
        rLeftHandSideMatrix.resize(24, 24, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(24, 24);

    // Resize the Right Hand Side if necessary,
    // and initialize it to Zero

    if(rRightHandSideVector.size() != 24)
        rRightHandSideVector.resize(24, false);
    noalias(rRightHandSideVector) = ZeroVector(24);

    // Get some references.

    const PropertiesType & props = GetProperties();
    const GeometryType & geom = GetGeometry();
    const Matrix & shapeFunctions = geom.ShapeFunctionsValues();
    Vector iN(shapeFunctions.size2());

    // Compute the local coordinate system.

    ShellQ4_LocalCoordinateSystem localCoordinateSystem(
        mpCoordinateTransformation->CreateLocalCoordinateSystem() );

    ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem() );

    // Prepare all the parameters needed for the MITC formulation.
    // This is to be done here outside the Gauss Loop.

    MITC4Params shearParameters(referenceCoordinateSystem);

    // Instantiate the Jacobian Operator.
    // This will store:
    // the jacobian matrix, its inverse, its determinant
    // and the derivatives of the shape functions in the local
    // coordinate system

    JacobianOperator jacOp;
    array_1d<double, 4> dArea;

    // Instantiate all strain-displacement matrices.

    Matrix B(8, 24, 0.0);
    Vector Bdrilling(24, 0.0);

    // Instantiate all section tangent matrices.

    Matrix D(8, 8, 0.0);
    double Ddrilling(0.0);

    // Instantiate strain and stress-resultant vectors

    Vector generalizedStrains(8);
    Vector generalizedStresses(8);

    // Get the current displacements in global coordinate system

    Vector globalDisplacements(24);
    GetValuesVector(globalDisplacements, 0);

    // Get the current displacements in local coordinate system

    Vector localDisplacements(
        mpCoordinateTransformation->CalculateLocalDisplacements(localCoordinateSystem, globalDisplacements) );

    // Instantiate the EAS Operator.
    // This will apply the Enhanced Assumed Strain Method for the calculation
    // of the membrane contribution.

    EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

    // auxiliary data

    Matrix BTD(24, 8); // auxiliary matrix to store the product B'*D

    // Initialize parameters for the cross section calculation

    ShellCrossSection::SectionParameters parameters(geom, props, rCurrentProcessInfo);
    parameters.SetGeneralizedStrainVector( generalizedStrains );
    parameters.SetGeneralizedStressVector( generalizedStresses );
    parameters.SetConstitutiveMatrix( D );
    Flags& options = parameters.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, CalculateResidualVectorFlag);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, CalculateStiffnessMatrixFlag);

    // Gauss Loop.
    for(int i = 0; i < 4; i++)
    {

        // get a reference of the current integration point and shape functions

        const GeometryType::IntegrationPointType & ip = geom.IntegrationPoints()[i];

        noalias( iN ) = row( shapeFunctions, i );

        // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
        // and Shape functions derivatives in the local coordinate system

        jacOp.Calculate(referenceCoordinateSystem, geom.ShapeFunctionLocalGradient(i));

        // compute the 'area' of the current integration point

        double dA = ip.Weight() * jacOp.Determinant();
        dArea[i] = dA;

        // Compute all strain-displacement matrices

        CalculateBMatrix( ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling );

        // Calculate strain vectors in local coordinate system

        noalias( generalizedStrains ) = prod( B, localDisplacements );
        double drillingStrain = inner_prod( Bdrilling, localDisplacements );

        // Apply the EAS method to modify the membrane part of the strains computed above.

        EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains, mEASStorage);

        // Calculate the response of the Cross Section

        ShellCrossSection::Pointer & section = mSections[i];

        parameters.SetShapeFunctionsValues( iN );
        parameters.SetShapeFunctionsDerivatives( jacOp.XYDerivatives() );
		//add in shear stabilization
		double shearStabilisation = CalculateStenbergShearStabilization(referenceCoordinateSystem, section->GetThickness(GetProperties()));
		parameters.SetStenbergShearStabilization(shearStabilisation);
        section->CalculateSectionResponse( parameters, ConstitutiveLaw::StressMeasure_PK2 );
        Ddrilling = section->GetDrillingStiffness();

        // multiply the section tangent matrices and stress resultants by 'dA'

        D *= dA;
        Ddrilling *= dA;
        generalizedStresses *= dA;
        double drillingStress = Ddrilling * drillingStrain; // already multiplied by 'dA'

        // Add all contributions to the Stiffness Matrix

        noalias( BTD ) = prod( trans( B ), D );
        noalias( rLeftHandSideMatrix ) += prod( BTD, B );
        noalias( rLeftHandSideMatrix ) += outer_prod( Ddrilling * Bdrilling, Bdrilling );

        // Add all contributions to the residual vector

        noalias( rRightHandSideVector ) -= prod( trans( B ), generalizedStresses );
        noalias( rRightHandSideVector ) -= Bdrilling * drillingStress;

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

    mpCoordinateTransformation->FinalizeCalculations(localCoordinateSystem,
            globalDisplacements,
            localDisplacements,
            rLeftHandSideMatrix,
            rRightHandSideVector,
            CalculateResidualVectorFlag,
            CalculateStiffnessMatrixFlag);

    // Add body forces contributions. This doesn't depend on the coordinate system
    AddBodyForces(dArea, rRightHandSideVector);
}

void ShellThickElement3D4N::AddBodyForces(const array_1d<double,4> & dA, VectorType& rRightHandSideVector)
{
    const GeometryType& geom = GetGeometry();

    // Get shape functions
    const Matrix & N = geom.ShapeFunctionsValues();

    // auxiliary
    array_1d<double, 3> bf;

    // gauss loop to integrate the external force vector
    for(unsigned int igauss = 0; igauss < 4; igauss++)
    {
        // get mass per unit area
        double mass_per_unit_area = mSections[igauss]->CalculateMassPerUnitArea(GetProperties());

        // interpolate nodal volume accelerations to this gauss point
        // and obtain the body force vector
        bf.clear();
        if (GetProperties().Has( VOLUME_ACCELERATION ))
            noalias(bf) = GetProperties()[VOLUME_ACCELERATION];
        else if (this->Has( VOLUME_ACCELERATION ))
            noalias(bf) =  this->GetValue(VOLUME_ACCELERATION);
        else {
            for(unsigned int inode = 0; inode < 4; inode++)
            {
                if( geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
                    noalias(bf) += N(igauss,inode) * geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
            }
        }
        bf *= (mass_per_unit_area * dA[igauss]);

        // add it to the RHS vector
        for(unsigned int inode = 0; inode < 4; inode++)
        {
            unsigned int index = inode*6;
            double iN = N(igauss,inode);
            rRightHandSideVector[index + 0] += iN * bf[0];
            rRightHandSideVector[index + 1] += iN * bf[1];
            rRightHandSideVector[index + 2] += iN * bf[2];
        }
    }
}

bool ShellThickElement3D4N::TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    // Check the required output

    int ijob = 0;
    bool bGlobal = false;
	CheckGeneralizedStressOrStrainOutput(rVariable, ijob, bGlobal);

    // quick return

    if(ijob == 0) return false;

    // resize output

    SizeType size = 4;
    if(rValues.size() != size)
        rValues.resize(size);

    // Get some references.

    const PropertiesType & props = GetProperties();
    const GeometryType & geom = GetGeometry();
    const Matrix & shapeFunctions = geom.ShapeFunctionsValues();
    Vector iN(shapeFunctions.size2());

    // Compute the local coordinate system.

    ShellQ4_LocalCoordinateSystem localCoordinateSystem(
        mpCoordinateTransformation->CreateLocalCoordinateSystem() );

    ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem() );

    // Prepare all the parameters needed for the MITC formulation.
    // This is to be done here outside the Gauss Loop.

    MITC4Params shearParameters(referenceCoordinateSystem);

    // Instantiate the Jacobian Operator.
    // This will store:
    // the jacobian matrix, its inverse, its determinant
    // and the derivatives of the shape functions in the local
    // coordinate system

    JacobianOperator jacOp;

    // Instantiate all strain-displacement matrices.

    Matrix B(8, 24, 0.0);
    Vector Bdrilling(24, 0.0);

    // Instantiate all section tangent matrices.

    Matrix D(8, 8, 0.0);

    // Instantiate strain and stress-resultant vectors

    Vector generalizedStrains(8);
    Vector generalizedStresses(8);
	std::vector<VectorType> rlaminateStrains;
	std::vector<VectorType> rlaminateStresses;

    // Get the current displacements in global coordinate system

    Vector globalDisplacements(24);
    GetValuesVector(globalDisplacements, 0);

    // Get the current displacements in local coordinate system

    Vector localDisplacements(
        mpCoordinateTransformation->CalculateLocalDisplacements(localCoordinateSystem, globalDisplacements) );

    // Instantiate the EAS Operator.
    // This will apply the Enhanced Assumed Strain Method for the calculation
    // of the membrane contribution.

    EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

    // Just to store the rotation matrix for visualization purposes
    Matrix R(8, 8);
    Matrix aux33(3, 3);

    // Initialize parameters for the cross section calculation

    ShellCrossSection::SectionParameters parameters(geom, props, rCurrentProcessInfo);
    parameters.SetGeneralizedStrainVector( generalizedStrains );
    parameters.SetGeneralizedStressVector( generalizedStresses );
    parameters.SetConstitutiveMatrix( D );
    Flags& options = parameters.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    // Gauss Loop

    for(unsigned int i = 0; i < size; i++)
    {

        // get a reference of the current integration point and shape functions

        const GeometryType::IntegrationPointType & ip = geom.IntegrationPoints()[i];

        noalias( iN ) = row( shapeFunctions, i );

        // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
        // and Shape functions derivatives in the local coordinate system

        jacOp.Calculate(referenceCoordinateSystem, geom.ShapeFunctionLocalGradient(i));

        // Compute all strain-displacement matrices

        CalculateBMatrix( ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling );

        // Calculate strain vectors in local coordinate system

        noalias( generalizedStrains ) = prod( B, localDisplacements );

        // Apply the EAS method to modify the membrane part of the strains computed above.

        EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains, mEASStorage);

        // Calculate the response of the Cross Section
        ShellCrossSection::Pointer & section = mSections[i];

		//Add in shear stabilization
		double shearStabilisation = CalculateStenbergShearStabilization(referenceCoordinateSystem, section->GetThickness(GetProperties()));
		parameters.SetStenbergShearStabilization(shearStabilisation);

        if(ijob > 2)
        {
            if (ijob > 7)
			{
				//Calculate lamina stresses
				CalculateLaminaStrains(section,generalizedStrains,rlaminateStrains);
				CalculateLaminaStresses(section,parameters,rlaminateStrains,rlaminateStresses);
			}
			else
			{
				// calculate force resultants
				parameters.SetShapeFunctionsValues(iN);
				parameters.SetShapeFunctionsDerivatives(jacOp.XYDerivatives());
				section->CalculateSectionResponse(parameters, ConstitutiveLaw::StressMeasure_PK2);

				if (ijob > 4)
				{
					// Compute stresses
					CalculateStressesFromForceResultants(generalizedStresses,
						section->GetThickness(GetProperties()));
				}
			}
        }

        // save the results

        DecimalCorrection( generalizedStrains );
        DecimalCorrection( generalizedStresses );

        // now the results are in the element coordinate system
        // if necessary, rotate the results in the section (local) coordinate system
        if(section->GetOrientationAngle() != 0.0 && !bGlobal)
        {
			if (ijob > 7)
			{
				section->GetRotationMatrixForGeneralizedStresses(-(section->GetOrientationAngle()), R);
				for (unsigned int i = 0; i < rlaminateStresses.size(); i++)
				{
					rlaminateStresses[i] = prod(R, rlaminateStresses[i]);
				}

				section->GetRotationMatrixForGeneralizedStrains(-(section->GetOrientationAngle()), R);
				for (unsigned int i = 0; i < rlaminateStrains.size(); i++)
				{
					rlaminateStrains[i] = prod(R, rlaminateStrains[i]);
				}
			}
			else if (ijob > 2)
            {
                section->GetRotationMatrixForGeneralizedStresses( -(section->GetOrientationAngle()), R );
                generalizedStresses = prod( R, generalizedStresses );
            }
            else
            {
                section->GetRotationMatrixForGeneralizedStrains( -(section->GetOrientationAngle()), R );
                generalizedStrains = prod( R, generalizedStrains );
            }
        }

        Matrix & iValue = rValues[i];
        if(iValue.size1() != 3 || iValue.size2() != 3)
            iValue.resize(3, 3, false);

        if(ijob == 1) // strains
        {
            iValue(0, 0) = generalizedStrains(0);
            iValue(1, 1) = generalizedStrains(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = 0.5 * generalizedStrains(2);
            iValue(0, 2) = iValue(2, 0) = 0.5 * generalizedStrains(7);
            iValue(1, 2) = iValue(2, 1) = 0.5 * generalizedStrains(6);
        }
        else if(ijob == 2) // curvatures
        {
            iValue(0, 0) = generalizedStrains(3);
            iValue(1, 1) = generalizedStrains(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = 0.5 * generalizedStrains(5);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        }
        else if(ijob == 3) // forces
        {
            iValue(0, 0) = generalizedStresses(0);
            iValue(1, 1) = generalizedStresses(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = generalizedStresses(2);
            iValue(0, 2) = iValue(2, 0) = generalizedStresses(7);
            iValue(1, 2) = iValue(2, 1) = generalizedStresses(6);
        }
        else if(ijob == 4) // moments
        {
            iValue(0, 0) = generalizedStresses(3);
            iValue(1, 1) = generalizedStresses(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = generalizedStresses(5);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        }
		else if (ijob == 5) // SHELL_STRESS_TOP_SURFACE
		{
			iValue(0, 0) = generalizedStresses(0) +
				generalizedStresses(3);
			iValue(1, 1) = generalizedStresses(1) +
				generalizedStresses(4);
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) = generalizedStresses[2] +
				generalizedStresses[5];
			iValue(0, 2) = iValue(2, 0) = 0.0;
			iValue(1, 2) = iValue(2, 1) = 0.0;
		}
		else if (ijob == 6) // SHELL_STRESS_MIDDLE_SURFACE
		{
			iValue(0, 0) = generalizedStresses(0);
			iValue(1, 1) = generalizedStresses(1);
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) = generalizedStresses[2];
			iValue(0, 2) = iValue(2, 0) = generalizedStresses[6];
			iValue(1, 2) = iValue(2, 1) = generalizedStresses[7];
		}
		else if (ijob == 7) // SHELL_STRESS_BOTTOM_SURFACE
		{
			iValue(0, 0) = generalizedStresses(0) -
				generalizedStresses(3);
			iValue(1, 1) = generalizedStresses(1) -
				generalizedStresses(4);
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) = generalizedStresses[2] -
				generalizedStresses[5];
			iValue(0, 2) = iValue(2, 0) = 0.0;
			iValue(1, 2) = iValue(2, 1) = 0.0;
		}
		else if (ijob == 8) // SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE
		{
			iValue(0, 0) =
				rlaminateStresses[rlaminateStresses.size() - 1][0];
			iValue(1, 1) =
				rlaminateStresses[rlaminateStresses.size() - 1][1];
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) =
				rlaminateStresses[rlaminateStresses.size() - 1][2];
			iValue(0, 2) = iValue(2, 0) = rlaminateStresses[rlaminateStresses.size() - 1][6];
			iValue(1, 2) = iValue(2, 1) = rlaminateStresses[rlaminateStresses.size() - 1][7];
		}
		else if (ijob == 9) // SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE
		{
			iValue(0, 0) = rlaminateStresses[0][0];
			iValue(1, 1) = rlaminateStresses[0][1];
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) = rlaminateStresses[0][2];
			iValue(0, 2) = iValue(2, 0) = rlaminateStresses[0][6];
			iValue(1, 2) = iValue(2, 1) = rlaminateStresses[0][7];
		}

        // if requested, rotate the results in the global coordinate system
        if(bGlobal)
        {
            const Matrix& RG = localCoordinateSystem.Orientation();
            noalias( aux33 ) = prod( trans( RG ), iValue );
            noalias( iValue ) = prod( aux33, RG );
        }

    } // Gauss Loop

    return true;
}

ShellCrossSection::SectionBehaviorType ShellThickElement3D4N::GetSectionBehavior()
{
    return ShellCrossSection::Thick;
}


// =====================================================================================
//
// Class ShellThickElement3D4N - Serialization
//
// =====================================================================================

void ShellThickElement3D4N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  BaseShellElement );
    rSerializer.save("CTr", mpCoordinateTransformation);
    rSerializer.save("EAS", mEASStorage);
}

void ShellThickElement3D4N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  BaseShellElement );
    rSerializer.load("CTr", mpCoordinateTransformation);
    rSerializer.load("EAS", mEASStorage);
}

}
