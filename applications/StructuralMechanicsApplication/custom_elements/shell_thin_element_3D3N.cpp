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

#include "shell_thin_element_3D3N.hpp"
#include "custom_utilities/shellt3_corotational_coordinate_transformation.hpp"
#include "custom_utilities/shell_utilities.h"

#include <string>
#include <iomanip>

//----------------------------------------
// preprocessors for the integration
// method used by this element.

//#define OPT_1_POINT_INTEGRATION

//----------------------------------------
// preprocessors to handle the output
// in case of 3 integration points

//#define OPT_USES_INTERIOR_GAUSS_POINTS

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X) ShellUtilities::InterpToStandardGaussPoints(X)
#endif // OPT_USES_INTERIOR_GAUSS_POINTS
#endif // OPT_1_POINT_INTEGRATION

//#define OPT_AVARAGE_RESULTS

namespace Kratos
{
// =====================================================================================
//
// CalculationData
//
// =====================================================================================

ShellThinElement3D3N::CalculationData::CalculationData(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
        const ProcessInfo& rCurrentProcessInfo)
    : LCS0( pCoordinateTransformation->CreateReferenceCoordinateSystem() )
    , LCS( pCoordinateTransformation->CreateLocalCoordinateSystem() )
    , CurrentProcessInfo(rCurrentProcessInfo)

{
}

// =====================================================================================
//
// Class ShellThinElement3D3N
//
// =====================================================================================

ShellThinElement3D3N::ShellThinElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        bool NLGeom)
    : BaseShellElement(NewId, pGeometry)
    , mpCoordinateTransformation( NLGeom ?
                                  new ShellT3_CorotationalCoordinateTransformation(pGeometry) :
                                  new ShellT3_CoordinateTransformation(pGeometry))
{
}

ShellThinElement3D3N::ShellThinElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        bool NLGeom)
    : BaseShellElement(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation( NLGeom ?
                                  new ShellT3_CorotationalCoordinateTransformation(pGeometry) :
                                  new ShellT3_CoordinateTransformation(pGeometry))
{
}

ShellThinElement3D3N::ShellThinElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        CoordinateTransformationBasePointerType pCoordinateTransformation)
    : BaseShellElement(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation(pCoordinateTransformation)
{
}

ShellThinElement3D3N::~ShellThinElement3D3N()
{
}

Element::Pointer ShellThinElement3D3N::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
    return Kratos::make_shared< ShellThinElement3D3N >(NewId, newGeom, pProperties, mpCoordinateTransformation->Create(newGeom) );
}

Element::Pointer ShellThinElement3D3N::Create(IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared< ShellThinElement3D3N >(NewId, pGeom, pProperties, mpCoordinateTransformation->Create(pGeom) );
}

void ShellThinElement3D3N::Initialize()
{
    KRATOS_TRY

    const int points_number = GetGeometry().PointsNumber();

    KRATOS_ERROR_IF_NOT(points_number == 3) <<"ShellThinElement3D3N - Wrong number of nodes"
        << points_number << std::endl;

    BaseShellElement::Initialize();

    mpCoordinateTransformation->Initialize();

    this->SetupOrientationAngles();

    KRATOS_CATCH("")
}

void ShellThinElement3D3N::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->InitializeNonLinearIteration(rCurrentProcessInfo);

    BaseInitializeNonLinearIteration(rCurrentProcessInfo);
}

void ShellThinElement3D3N::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->FinalizeNonLinearIteration(rCurrentProcessInfo);

    BaseFinalizeNonLinearIteration(rCurrentProcessInfo);
}

void ShellThinElement3D3N::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    BaseInitializeSolutionStep(rCurrentProcessInfo);

    mpCoordinateTransformation->InitializeSolutionStep(rCurrentProcessInfo);
}

void ShellThinElement3D3N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    BaseFinalizeSolutionStep(rCurrentProcessInfo);

    mpCoordinateTransformation->FinalizeSolutionStep(rCurrentProcessInfo);
}

void ShellThinElement3D3N::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if((rMassMatrix.size1() != 18) || (rMassMatrix.size2() != 18))
        rMassMatrix.resize(18, 18, false);
    noalias(rMassMatrix) = ZeroMatrix(18, 18);

    // Compute the local coordinate system.

    ShellT3_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem() );

    // lumped area

    double lump_area = referenceCoordinateSystem.Area() / 3.0;

    // Calculate avarage mass per unit area

    const SizeType num_gps = GetNumberOfGPs();

    double av_mass_per_unit_area = 0.0;
    for(SizeType i = 0; i < num_gps; i++)
        av_mass_per_unit_area += mSections[i]->CalculateMassPerUnitArea(GetProperties());
    av_mass_per_unit_area /= double(num_gps);

    // loop on nodes
    for(SizeType i = 0; i < 3; i++)
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
// Class ShellThinElement3D3N - Results on Gauss Points
//
// =====================================================================================

void ShellThinElement3D3N::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_gps = GetNumberOfGPs();

    if (rValues.size() != num_gps)
        rValues.resize(num_gps);

    // The membrane formulation needs to iterate to find the correct
    // mid-surface strain values.
    // Check if we are doing a non-linear analysis type. If not, print warning

    if (!rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
    {
        KRATOS_WARNING_ONCE("ShellThinElement3D3N") << "Warning: Gauss point results have "
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
		// Initialize common calculation variables
		CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
		data.CalculateLHS = false;
		data.CalculateRHS = true;
		InitializeCalculationData(data);

		// Gauss Loop.
		for (SizeType i = 0; i < num_gps; i++)
		{
			// set the current integration point index
			data.gpIndex = i;
			ShellCrossSection::Pointer& section = mSections[i];

			// calculate beta0
			CalculateBeta0(data);

			// calculate the total strain displ. matrix
			CalculateBMatrix(data);

			// compute generalized strains
			noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

			// calculate section response
			CalculateSectionResponse(data);

			// Compute stresses
			CalculateStressesFromForceResultants(data.generalizedStresses,
				section->GetThickness(GetProperties()));

			// calculate von mises stress
			CalculateVonMisesStress(data, rVariable, rValues[i]);

		} // end gauss loop
	}
	else if (rVariable == TSAI_WU_RESERVE_FACTOR)
	{
		// resize output
		if (rValues.size() != num_gps)
			rValues.resize(num_gps);

		// Just to store the rotation matrix for visualization purposes
		Matrix R(mStrainSize, mStrainSize);

		// Initialize common calculation variables
		CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
		data.CalculateLHS = false;
		data.CalculateRHS = true;
		InitializeCalculationData(data);

		// Get all laminae strengths
		const PropertiesType & props = GetProperties();
		ShellCrossSection::Pointer & section = mSections[0];
		std::vector<Matrix> Laminae_Strengths =
			std::vector<Matrix>(section->NumberOfPlies());
		for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
		{
			Laminae_Strengths[ply].resize(3, 3, 0.0);
			Laminae_Strengths[ply].clear();
		}
		section->GetLaminaeStrengths(Laminae_Strengths, props);

		// Retrieve ply orientations
		Vector ply_orientation(section->NumberOfPlies());
		section->GetLaminaeOrientation(props, ply_orientation);
		double total_rotation = 0.0;

		// Gauss Loop.
		for (SizeType i = 0; i < num_gps; i++)
		{
			// set the current integration point index
			data.gpIndex = i;
			ShellCrossSection::Pointer& section = mSections[i];

			// calculate beta0
			CalculateBeta0(data);

			// calculate the total strain displ. matrix
			CalculateBMatrix(data);

			// compute generalized strains
			noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

			//Calculate lamina stresses
			CalculateLaminaStrains(data);
			CalculateLaminaStresses(data);

			// Rotate lamina stress to lamina material principal directions
			for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
			{
				total_rotation = -ply_orientation[ply] - (section->GetOrientationAngle());
				section->GetRotationMatrixForGeneralizedStresses(total_rotation, R);
				//top surface of current ply
				data.rlaminateStresses[2 * ply] = prod(R, data.rlaminateStresses[2 * ply]);
				//bottom surface of current ply
				data.rlaminateStresses[2 * ply + 1] = prod(R, data.rlaminateStresses[2 * ply + 1]);
			}

			// Calculate Tsai-Wu criterion for each ply, take min of all plies
			double min_tsai_wu = 0.0;
			double temp_tsai_wu = 0.0;
			for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
			{
				temp_tsai_wu = CalculateTsaiWuPlaneStress(data, Laminae_Strengths[ply], ply);
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

		}// Gauss loop
	}
	else
	{
		for (SizeType i = 0; i < num_gps; i++)
			mSections[i]->GetValue(rVariable, GetProperties(), rValues[i]);
	}

    if(this->Has(rVariable))
    {
        // Get result value for output
        const auto& output_value = this->GetValue(rVariable);

        // Write the same result on all Gauss-Points
        for(IndexType i = 0; i < num_gps; ++i)
            rValues[i] = output_value;
    }

    OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
}

void ShellThinElement3D3N::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_gps = GetNumberOfGPs();

	if (rVariable == LOCAL_AXIS_1)
	{
		// LOCAL_AXIS_1 output DOES NOT include the effect of section
		// orientation, which rotates the entrire element section in-plane
		// and is used in element stiffness calculation.

        if (rValues.size() != num_gps) rValues.resize(num_gps);

		for (SizeType i = 0; i < num_gps; ++i) rValues[i] = ZeroVector(3);
		// Initialize common calculation variables
		ShellT3_LocalCoordinateSystem localCoordinateSystem(mpCoordinateTransformation->CreateReferenceCoordinateSystem());

		for (SizeType GP = 0; GP < 1; GP++)
		{
			rValues[GP] = localCoordinateSystem.Vx();
		}
	}
	else if (rVariable == LOCAL_MATERIAL_ORIENTATION_VECTOR_1)
	{
		// LOCAL_MATERIAL_ORIENTATION_VECTOR_1 output DOES include the effect of
		// section orientation, which rotates the entrire element section
		// in-plane and is used in the element stiffness calculation.

		// Resize output
		if (rValues.size() != num_gps) rValues.resize(num_gps);
		for (SizeType i = 0; i < num_gps; ++i) rValues[i] = ZeroVector(3);

		// Initialize common calculation variables
		// Compute the local coordinate system.
		ShellT3_LocalCoordinateSystem localCoordinateSystem(mpCoordinateTransformation->CreateReferenceCoordinateSystem());

		// Get local axis 1 in flattened LCS space
		Vector3 localAxis1 = localCoordinateSystem.P2() - localCoordinateSystem.P1();

		// Perform rotation of local axis 1 to fiber1 in flattened LCS space
		Matrix localToFiberRotation = Matrix(3, 3, 0.0);
		double fiberSectionRotation = mSections[0]->GetOrientationAngle();
		double c = std::cos(fiberSectionRotation);
		double s = std::sin(fiberSectionRotation);
		localToFiberRotation(0, 0) = c;
		localToFiberRotation(0, 1) = -s;
		localToFiberRotation(1, 0) = s;
		localToFiberRotation(1, 1) = c;
		localToFiberRotation(2, 2) = 1.0;

		Vector3 temp = prod(localToFiberRotation, localAxis1);

		// Transform result back to global cartesian coords and normalize
		Matrix localToGlobalSmall = localCoordinateSystem.Orientation();
		Vector3 fiberAxis1 = prod(trans(localToGlobalSmall), temp);
		fiberAxis1 /= std::sqrt(inner_prod(fiberAxis1, fiberAxis1));

		//write results
		for (SizeType dir = 0; dir < 1; dir++)
		{
			rValues[dir] = fiberAxis1;
		}
	}
}

void ShellThinElement3D3N::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    // The membrane formulation needs to iterate to find the correct
    // mid-surface strain values.
    // Check if we are doing a non-linear analysis type. If not, print warning

    if(this->Has(rVariable))
    {
        // Resize Output
        const SizeType  write_points_number = GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != write_points_number)
            rValues.resize(write_points_number);

        // Get result value for output
        const auto& output_value = this->GetValue(rVariable);

        // Write the same result on all Gauss-Points
        for(IndexType i = 0; i < write_points_number; ++i)
            rValues[i] = output_value;
    }

	if (this->Id() == 1)
	{
		if (!rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
		{
			std::cout << "\nWARNING:\nGauss point results have been requested for a linear analysis."
				<< "\nThe membrane formulation used in the specified shell element"
				<< "(ShellThinElement3D3N) requires iteration to accurately determine "
				<< "recovered quantities (strain, stress, etc...).\n"
				<< "Please switch to 'analysis_type = Non-Linear' in your json file for accurate recovered quantities."
				<< std::endl;
		}
	}

    if(TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(rVariable, rValues, rCurrentProcessInfo)) return;

}

void ShellThinElement3D3N::CalculateOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
        std::vector<array_1d<double,3> >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    if(TryCalculateOnIntegrationPoints_MaterialOrientation(rVariable, rValues, rCurrentProcessInfo)) return;
}

void ShellThinElement3D3N::Calculate(const Variable<Matrix>& rVariable, Matrix & Output, const ProcessInfo & rCurrentProcessInfo)
{
	if (rVariable == LOCAL_ELEMENT_ORIENTATION)
	{
		Output.resize(3, 3, false);

		// Compute the local coordinate system.
		ShellT3_LocalCoordinateSystem localCoordinateSystem(mpCoordinateTransformation->CreateReferenceCoordinateSystem());
		Output = localCoordinateSystem.Orientation();
	}
}

// =====================================================================================
//
// Class ShellThinElement3D3N - Private methods
//
// =====================================================================================

void ShellThinElement3D3N::CheckGeneralizedStressOrStrainOutput(const Variable<Matrix>& rVariable, int & ijob, bool & bGlobal)
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

void ShellThinElement3D3N::CalculateStressesFromForceResultants(VectorType & rstresses, const double & rthickness)
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
}

void ShellThinElement3D3N::CalculateLaminaStrains(CalculationData & data)
{
	ShellCrossSection::Pointer& section = mSections[data.gpIndex];

	// Get laminate properties
	double thickness = section->GetThickness(GetProperties());
	double z_current = thickness / -2.0; // start from the top of the 1st layer

										 // Establish current strains at the midplane
										 // (element coordinate system)
	double e_x = data.generalizedStrains[0];
	double e_y = data.generalizedStrains[1];
	double e_xy = data.generalizedStrains[2];	//this is still engineering
												//strain (2xtensorial shear)
	double kap_x = data.generalizedStrains[3];
	double kap_y = data.generalizedStrains[4];
	double kap_xy = data.generalizedStrains[5];	//this is still engineering

												// Get ply thicknesses
	Vector ply_thicknesses = Vector(section->NumberOfPlies(), 0.0);
	section->GetPlyThicknesses(GetProperties(), ply_thicknesses);

	// Resize output vector. 2 Surfaces for each ply
	data.rlaminateStrains.resize(2 * section->NumberOfPlies());
	for (unsigned int i = 0; i < 2 * section->NumberOfPlies(); i++)
	{
		data.rlaminateStrains[i].resize(6, false);
		data.rlaminateStrains[i].clear();
	}

	// Loop over all plies - start from bottom ply, bottom surface
	for (unsigned int plyNumber = 0;
		plyNumber < section->NumberOfPlies(); ++plyNumber)
	{
		// Calculate strains at top surface, arranged in columns.
		// (element coordinate system)
		data.rlaminateStrains[2 * plyNumber][0] = e_x + z_current*kap_x;
		data.rlaminateStrains[2 * plyNumber][1] = e_y + z_current*kap_y;
		data.rlaminateStrains[2 * plyNumber][2] = e_xy + z_current*kap_xy;

		// Move to bottom surface of current layer
		z_current += ply_thicknesses[plyNumber];

		// Calculate strains at bottom surface, arranged in columns
		// (element coordinate system)
		data.rlaminateStrains[2 * plyNumber + 1][0] = e_x + z_current*kap_x;
		data.rlaminateStrains[2 * plyNumber + 1][1] = e_y + z_current*kap_y;
		data.rlaminateStrains[2 * plyNumber + 1][2] = e_xy + z_current*kap_xy;
	}
}

void ShellThinElement3D3N::CalculateLaminaStresses(CalculationData & data)
{
	ShellCrossSection::Pointer& section = mSections[data.gpIndex];

	// Setup flag to compute ply constitutive matrices
	// (units [Pa] and rotated to element orientation)
	section->SetupGetPlyConstitutiveMatrices();
	Flags& options = data.SectionParameters.GetOptions();
	options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
	section->CalculateSectionResponse(data.SectionParameters,
		ConstitutiveLaw::StressMeasure_PK2);
	CalculateSectionResponse(data);

	// Resize output vector. 2 Surfaces for each ply
	data.rlaminateStresses.resize(2 * section->NumberOfPlies());
	for (unsigned int i = 0; i < 2 * section->NumberOfPlies(); i++)
	{
		data.rlaminateStresses[i].resize(6, false);
		data.rlaminateStresses[i].clear();
	}

	// Loop over all plies - start from top ply, top surface
	for (unsigned int plyNumber = 0;
		plyNumber < section->NumberOfPlies(); ++plyNumber)
	{
		// determine stresses at currrent ply, top surface
		// (element coordinate system)
		data.rlaminateStresses[2 * plyNumber] = prod(
			section->GetPlyConstitutiveMatrix(plyNumber),
			data.rlaminateStrains[2 * plyNumber]);

		// determine stresses at currrent ply, bottom surface
		// (element coordinate system)
		data.rlaminateStresses[2 * plyNumber + 1] = prod(
			section->GetPlyConstitutiveMatrix(plyNumber),
			data.rlaminateStrains[2 * plyNumber + 1]);
	}
}

double ShellThinElement3D3N::CalculateTsaiWuPlaneStress(const CalculationData & data, const Matrix & rLamina_Strengths, const unsigned int & rPly)
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

	// Should be FALSE unless testing against other programs that ignore it.
	bool disable_in_plane_interaction = false;

	// First, F_i
	Vector F_i = Vector(3, 0.0);
	F_i[0] = 1.0 / rLamina_Strengths(0, 0) - 1.0 / rLamina_Strengths(0, 1);
	F_i[1] = 1.0 / rLamina_Strengths(0, 2) - 1.0 / rLamina_Strengths(1, 0);
	F_i[2] = 0.0;

	// Second, F_ij
	Matrix F_ij = Matrix(3, 3, 0.0);
	F_ij.clear();
	F_ij(0, 0) = 1.0 / rLamina_Strengths(0, 0) / rLamina_Strengths(0, 1);	// 11
	F_ij(1, 1) = 1.0 / rLamina_Strengths(0, 2) / rLamina_Strengths(1, 0);	// 22
	F_ij(2, 2) = 1.0 / rLamina_Strengths(1, 1) / rLamina_Strengths(1, 1);	// 12
	F_ij(0, 1) = F_ij(1, 0) = -0.5 / std::sqrt(rLamina_Strengths(0, 0)*rLamina_Strengths(0, 1)*rLamina_Strengths(0, 2)*rLamina_Strengths(1, 0));

	if (disable_in_plane_interaction)
	{
		F_ij(0, 1) = F_ij(1, 0) = 0.0;
	}


	// Evaluate Tsai-Wu @ top surface of current layer
	double var_a = 0.0;
	double var_b = 0.0;
	for (SizeType i = 0; i < 3; i++)
	{
		var_b += F_i[i] * data.rlaminateStresses[2 * rPly][i];
		for (SizeType j = 0; j < 3; j++)
		{
			var_a += F_ij(i, j)*data.rlaminateStresses[2 * rPly][i] * data.rlaminateStresses[2 * rPly][j];
		}
	}
	double tsai_reserve_factor_top = (-1.0*var_b + std::sqrt(var_b*var_b + 4.0 * var_a)) / 2.0 / var_a;

	// Evaluate Tsai-Wu @ bottom surface of current layer
	var_a = 0.0;
	var_b = 0.0;
	for (SizeType i = 0; i < 3; i++)
	{
		var_b += F_i[i] * data.rlaminateStresses[2 * rPly + 1][i];
		for (SizeType j = 0; j < 3; j++)
		{
			var_a += F_ij(i, j)*data.rlaminateStresses[2 * rPly + 1][i] * data.rlaminateStresses[2 * rPly + 1][j];
		}
	}
	double tsai_reserve_factor_bottom = (-1.0*var_b + std::sqrt(var_b*var_b + 4.0 * var_a)) / 2.0 / var_a;

	// Return min of both surfaces as the result for the whole ply
	return std::min(tsai_reserve_factor_bottom, tsai_reserve_factor_top);
}

void ShellThinElement3D3N::CalculateVonMisesStress(const CalculationData & data, const Variable<double>& rVariable, double & rVon_Mises_Result)
{
	// calc von mises stresses at top mid and bottom surfaces for
	// thin shell
	double von_mises_top, von_mises_mid, von_mises_bottom;
	double sxx, syy, sxy;

	// top surface: membrane and +bending contributions
	//				(no transverse shear)
	sxx = data.generalizedStresses[0] + data.generalizedStresses[3];
	syy = data.generalizedStresses[1] + data.generalizedStresses[4];
	sxy = data.generalizedStresses[2] + data.generalizedStresses[5];
	von_mises_top = sxx*sxx - sxx*syy + syy*syy + 3.0*sxy*sxy;

	// mid surface: membrane only contributions
	//				(no bending or transverse shear)
	sxx = data.generalizedStresses[0];
	syy = data.generalizedStresses[1];
	sxy = data.generalizedStresses[2];
	von_mises_mid = sxx*sxx - sxx*syy + syy*syy +
		3.0*(sxy*sxy);

	// bottom surface:	membrane and bending contributions
	//					(no transverse shear)
	sxx = data.generalizedStresses[0] - data.generalizedStresses[3];
	syy = data.generalizedStresses[1] - data.generalizedStresses[4];
	sxy = data.generalizedStresses[2] - data.generalizedStresses[5];
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

void ShellThinElement3D3N::DecimalCorrection(Vector& a)
{
    double norm = norm_2(a);
    double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
    for(SizeType i = 0; i < a.size(); i++)
        if(std::abs(a(i)) < tolerance)
            a(i) = 0.0;
}

void ShellThinElement3D3N::SetupOrientationAngles()
{
    if (this->Has(MATERIAL_ORIENTATION_ANGLE))
    {
        for (CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it)
        (*it)->SetOrientationAngle(this->GetValue(MATERIAL_ORIENTATION_ANGLE));
    }
    else
    {
        ShellT3_LocalCoordinateSystem lcs( mpCoordinateTransformation->CreateReferenceCoordinateSystem() );

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

void ShellThinElement3D3N::InitializeCalculationData(CalculationData& data)
{
    //-------------------------------------
    // Computation of all stuff that remain
    // constant throughout the calculations

    //-------------------------------------
    // geometry data

    const double x12 = data.LCS0.X1() - data.LCS0.X2();
    const double x23 = data.LCS0.X2() - data.LCS0.X3();
    const double x31 = data.LCS0.X3() - data.LCS0.X1();
    const double x21 = -x12;
    const double x32 = -x23;
    const double x13 = -x31;

    const double y12 = data.LCS0.Y1() - data.LCS0.Y2();
    const double y23 = data.LCS0.Y2() - data.LCS0.Y3();
    const double y31 = data.LCS0.Y3() - data.LCS0.Y1();
    const double y21 = -y12;
    const double y32 = -y23;
    const double y13 = -y31;

    const double A  = 0.5*(y21*x13 - x21*y13);
    const double A2 = 2.0*A;
    const double A4 = 4.0*A;
    const double AA4 = A * A4;

    const double LL21 = x21*x21 + y21*y21;
    const double LL32 = x32*x32 + y32*y32;
    const double LL13 = x13*x13 + y13*y13;

    // Note: here we compute the avarage thickness,
    // since L is constant over the element.
    // Now it is not necessary to compute the avarage
    // because the current implementation of the cross section
    // doesn't have a variable thickness
    // (for example as a function of the spatial coordinates...).
    // This is just a place-holder for future
    // implementation of a variable thickness

    double h = 0.0;
    for(unsigned int i = 0; i < mSections.size(); i++)
        h += mSections[i]->GetThickness(GetProperties());
    h /= (double)mSections.size();

    data.hMean = h;
    data.TotalArea = A;
    data.TotalVolume = A * h;

    // this is the integration weight
    // used during the gauss loop.
    // it is dArea because it will
    // multiply section stress resultants
    // and section constitutive matrices
    // that already take into accout the
    // thickness

    const SizeType num_gps = GetNumberOfGPs();

    data.dA = A / (double)num_gps;

    // crete the integration point locations
    if(data.gpLocations.size() != 0) data.gpLocations.clear();
    data.gpLocations.resize( num_gps );
#ifdef OPT_1_POINT_INTEGRATION
    array_1d<double,3>& gp0 = data.gpLocations[0];
    gp0[0] = 1.0/3.0;
    gp0[1] = 1.0/3.0;
    gp0[2] = 1.0/3.0;
#else
    array_1d<double,3>& gp0 = data.gpLocations[0];
    array_1d<double,3>& gp1 = data.gpLocations[1];
    array_1d<double,3>& gp2 = data.gpLocations[2];
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
    gp0[0] = 1.0/6.0;
    gp0[1] = 1.0/6.0;
    gp0[2] = 2.0/3.0;
    gp1[0] = 2.0/3.0;
    gp1[1] = 1.0/6.0;
    gp1[2] = 1.0/6.0;
    gp2[0] = 1.0/6.0;
    gp2[1] = 2.0/3.0;
    gp2[2] = 1.0/6.0;
#else
    gp0[0] = 0.5;
    gp0[1] = 0.5;
    gp0[2] = 0.0;
    gp1[0] = 0.0;
    gp1[1] = 0.5;
    gp1[2] = 0.5;
    gp2[0] = 0.5;
    gp2[1] = 0.0;
    gp2[2] = 0.5;
#endif // OPT_USES_INTERIOR_GAUSS_POINTS

#endif // OPT_1_POINT_INTEGRATION

    // cartesian derivatives
    data.dNxy.resize(3, 2, false);
    data.dNxy(0, 0) = (y13 - y12)/A2;
    data.dNxy(0, 1) = (x12 - x13)/A2;
    data.dNxy(1, 0) =        -y13/A2;
    data.dNxy(1, 1) =         x13/A2;
    data.dNxy(2, 0) =         y12/A2;
    data.dNxy(2, 1) =        -x12/A2;

    //-------------------------------------
    // template parameters

    const double b1    =  1.0;
    const double b2    =  2.0;
    const double b3    =  1.0;
    const double b4    =  0.0;
    const double b5    =  1.0;
    const double b6    = -1.0;
    const double b7    = -1.0;
    const double b8    = -1.0;
    const double b9    = -2.0;

    const double alpha   =  1.5;
    const double alpha_6 =  alpha / 6.0;

    //--------------------------------------
    // calculate L - Lumping matrix
    // for the construction of the basic
    // stiffness

    double L_mult = 0.5/A;

    data.L.resize(3, 9, false);

    data.L(0, 0) = L_mult * y23;
    data.L(1, 0) = 0.00;
    data.L(2, 0) = L_mult * x32;
    data.L(0, 1) = 0.00;
    data.L(1, 1) = L_mult * x32;
    data.L(2, 1) = L_mult * y23;
    data.L(0, 2) = L_mult * y23*(y13-y21)*alpha_6;
    data.L(1, 2) = L_mult * x32*(x31-x12)*alpha_6;
    data.L(2, 2) = L_mult * 2.00*(x31*y13-x12*y21)*alpha_6;
    data.L(0, 3) = L_mult * y31;
    data.L(1, 3) = 0.00;
    data.L(2, 3) = L_mult * x13;
    data.L(0, 4) = 0.00;
    data.L(1, 4) = L_mult * x13;
    data.L(2, 4) = L_mult * y31;
    data.L(0, 5) = L_mult * y31*(y21-y32)*alpha_6;
    data.L(1, 5) = L_mult * x13*(x12-x23)*alpha_6;
    data.L(2, 5) = L_mult * 2.00*(x12*y21-x23*y32)*alpha_6;
    data.L(0, 6) = L_mult * y12;
    data.L(1, 6) = 0.00;
    data.L(2, 6) = L_mult * x21;
    data.L(0, 7) = 0.00;
    data.L(1, 7) = L_mult * x21;
    data.L(2, 7) = L_mult * y12;
    data.L(0, 8) = L_mult * y12*(y32-y13)*alpha_6;
    data.L(1, 8) = L_mult * x21*(x23-x31)*alpha_6;
    data.L(2, 8) = L_mult * 2.00*(x23*y32-x31*y13)*alpha_6;

    //--------------------------------------
    // calculate Q1,Q2,Q3 - matrices
    // for the construction of the
    // higher order stiffness

    data.Q1.resize(3, 3, false);

    data.Q1(0,0) = b1*A2/(LL21*3.00);
    data.Q1(0,1) = b2*A2/(LL21*3.00);
    data.Q1(0,2) = b3*A2/(LL21*3.00);
    data.Q1(1,0) = b4*A2/(LL32*3.00);
    data.Q1(1,1) = b5*A2/(LL32*3.00);
    data.Q1(1,2) = b6*A2/(LL32*3.00);
    data.Q1(2,0) = b7*A2/(LL13*3.00);
    data.Q1(2,1) = b8*A2/(LL13*3.00);
    data.Q1(2,2) = b9*A2/(LL13*3.00);

    data.Q2.resize(3, 3, false);

    data.Q2(0,0) = b9*A2/(LL21*3.00);
    data.Q2(0,1) = b7*A2/(LL21*3.00);
    data.Q2(0,2) = b8*A2/(LL21*3.00);
    data.Q2(1,0) = b3*A2/(LL32*3.00);
    data.Q2(1,1) = b1*A2/(LL32*3.00);
    data.Q2(1,2) = b2*A2/(LL32*3.00);
    data.Q2(2,0) = b6*A2/(LL13*3.00);
    data.Q2(2,1) = b4*A2/(LL13*3.00);
    data.Q2(2,2) = b5*A2/(LL13*3.00);

    data.Q3.resize(3, 3, false);

    data.Q3(0,0) = b5*A2/(LL21*3.00);
    data.Q3(0,1) = b6*A2/(LL21*3.00);
    data.Q3(0,2) = b4*A2/(LL21*3.00);
    data.Q3(1,0) = b8*A2/(LL32*3.00);
    data.Q3(1,1) = b9*A2/(LL32*3.00);
    data.Q3(1,2) = b7*A2/(LL32*3.00);
    data.Q3(2,0) = b2*A2/(LL13*3.00);
    data.Q3(2,1) = b3*A2/(LL13*3.00);
    data.Q3(2,2) = b1*A2/(LL13*3.00);

    //--------------------------------------
    // calculate Te, TTu -
    // transformation matrices
    // for the construction of the
    // higher order stiffness

    data.Te.resize(3, 3, false);

    data.Te(0,0) = 1.0/AA4 * y23*y13*LL21;
    data.Te(0,1) = 1.0/AA4 * y31*y21*LL32;
    data.Te(0,2) = 1.0/AA4 * y12*y32*LL13;
    data.Te(1,0) = 1.0/AA4 * x23*x13*LL21;
    data.Te(1,1) = 1.0/AA4 * x31*x21*LL32;
    data.Te(1,2) = 1.0/AA4 * x12*x32*LL13;
    data.Te(2,0) = 1.0/AA4 * (y23*x31+x32*y13)*LL21;
    data.Te(2,1) = 1.0/AA4 * (y31*x12+x13*y21)*LL32;
    data.Te(2,2) = 1.0/AA4 * (y12*x23+x21*y32)*LL13;

    data.TTu.resize(3, 9, false);

    for(unsigned int i=0; i<3; i++)
    {
        data.TTu(i, 0) = 1.0/A4 * x32;
        data.TTu(i, 1) = 1.0/A4 * y32;
        data.TTu(i, 2) = 0.0;
        data.TTu(i, 3) = 1.0/A4 * x13;
        data.TTu(i, 4) = 1.0/A4 * y13;
        data.TTu(i, 5) = 0.0;
        data.TTu(i, 6) = 1.0/A4 * x21;
        data.TTu(i, 7) = 1.0/A4 * y21;
        data.TTu(i, 8) = 0.0;
    }
    data.TTu(0, 2) = 1.0;
    data.TTu(1, 5) = 1.0;
    data.TTu(2, 8) = 1.0;

    //--------------------------------------
    // calculate the displacement vector
    // in global and local coordinate systems

    data.globalDisplacements.resize(18, false);
    GetValuesVector( data.globalDisplacements );

    data.localDisplacements =
        mpCoordinateTransformation->CalculateLocalDisplacements(
            data.LCS, data.globalDisplacements);

    //--------------------------------------
    // Finally allocate all auxiliary
    // matrices to be used later on
    // during the element integration.
    // Just to avoid re-allocations

    data.B.resize(mStrainSize, 18, false);
    data.D.resize(mStrainSize, mStrainSize, false);
    data.BTD.resize(18, mStrainSize, false);

    data.generalizedStrains.resize(mStrainSize, false);
    data.generalizedStresses.resize(mStrainSize, false);

    data.N.resize(3, false);

    data.Q.resize(3, 3, false);
    data.Qh.resize(3, 9, false);
    data.TeQ.resize(3, 3, false);

    data.H1.resize(9, false);
    data.H2.resize(9, false);
    data.H3.resize(9, false);
    data.H4.resize(9, false);
    data.Bb.resize(3, 9, false);

    //--------------------------------------
    // Initialize the section parameters

    data.SectionParameters.SetElementGeometry( GetGeometry() );
    data.SectionParameters.SetMaterialProperties( GetProperties() );
    data.SectionParameters.SetProcessInfo( data.CurrentProcessInfo );

    data.SectionParameters.SetGeneralizedStrainVector( data.generalizedStrains );
    data.SectionParameters.SetGeneralizedStressVector( data.generalizedStresses );
    data.SectionParameters.SetConstitutiveMatrix( data.D );

    data.SectionParameters.SetShapeFunctionsDerivatives( data.dNxy );

    Flags& options = data.SectionParameters.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, data.CalculateRHS);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, data.CalculateLHS);
}

void ShellThinElement3D3N::CalculateBMatrix(CalculationData& data)
{
    //---------------------------------------------
    // geom data
    array_1d<double, 3>& gpLoc = data.gpLocations[data.gpIndex];
    double loc1 = gpLoc[0];
    double loc2 = gpLoc[1];
    double loc3 = gpLoc[2];

    const double x12 = data.LCS0.X1() - data.LCS0.X2();
    const double x23 = data.LCS0.X2() - data.LCS0.X3();
    const double x31 = data.LCS0.X3() - data.LCS0.X1();
    const double y12 = data.LCS0.Y1() - data.LCS0.Y2();
    const double y23 = data.LCS0.Y2() - data.LCS0.Y3();
    const double y31 = data.LCS0.Y3() - data.LCS0.Y1();

    //---------------------------------------------
    // set to zero the total B matrix
    data.B.clear();

    //---------------------------------------------
    // membrane basic part L
    // already computed.it is constant over the
    // element

    //---------------------------------------------
    // membrane higher order part Qh
    noalias( data.Q )  = loc1 * data.Q1;
    noalias( data.Q ) += loc2 * data.Q2;
    noalias( data.Q ) += loc3 * data.Q3;
    noalias( data.TeQ ) = prod( data.Te, data.Q );
    noalias( data.Qh  ) = /*1.5**/std::sqrt(data.beta0) * prod( data.TeQ, data.TTu );

    //---------------------------------------------
    // compute the bending B matrix (DKT)

    const double LL12 = x12*x12 + y12*y12;
    const double LL23 = x23*x23 + y23*y23;
    const double LL31 = x31*x31 + y31*y31;

    const double p4 = -6.0*x23/(LL23);
    const double p5 = -6.0*x31/(LL31);
    const double p6 = -6.0*x12/(LL12);

    const double t4 = -6.0*y23/(LL23);
    const double t5 = -6.0*y31/(LL31);
    const double t6 = -6.0*y12/(LL12);

    const double q4 = 3.0*x23*y23/(LL23);
    const double q5 = 3.0*x31*y31/(LL31);
    const double q6 = 3.0*x12*y12/(LL12);

    const double r4 = 3.0*y23*y23/(LL23);
    const double r5 = 3.0*y31*y31/(LL31);
    const double r6 = 3.0*y12*y12/(LL12);

    const double area = 0.5*(y12*x31 - x12*y31);

    data.H1[0] = p6*(1.0-2.0*loc2) + (p5-p6)*loc3;
    data.H1[1] = q6*(1.0-2.0*loc2) - (q5+q6)*loc3;
    data.H1[2]= -4.0+6.0*(loc2+loc3) + r6*(1.0-2*loc2) - loc3*(r5+r6);
    data.H1[3] = -p6*(1.0-2.0*loc2) + (p4+p6)*loc3;
    data.H1[4] = q6*(1.0-2.0*loc2) - (q6-q4)*loc3;
    data.H1[5]= -2.0+6.0*loc2  + r6*(1-2.0*loc2) + loc3*(r4-r6);
    data.H1[6] = -(p4+p5)*loc3;
    data.H1[7] = (q4-q5)*loc3;
    data.H1[8]= -loc3*(r5-r4);

    data.H2[0] = t6*(1.0-2.0*loc2) + (t5-t6)*loc3;
    data.H2[1] = 1.0 + r6*(1.0-2.0*loc2) - loc3*(r5+r6);
    data.H2[2]= -q6*(1.0-2.0*loc2)+loc3*(q5+q6);
    data.H2[3] = -t6*(1.0-2.0*loc2) + (t4+t6)*loc3;
    data.H2[4] = -1.0 + r6*(1-2*loc2) + loc3*(r4-r6);
    data.H2[5]= -q6*(1.0-2.0*loc2)-loc3*(q4-q6);
    data.H2[6] = -(t5+t4)*loc3;
    data.H2[7] = (r4-r5)*loc3;
    data.H2[8]= -loc3*(q4-q5);

    data.H3[0] = -p5*(1.0-2.0*loc3) - (p6-p5)*loc2;
    data.H3[1] = q5*(1.0-2.0*loc3) - (q5+q6)*loc2;
    data.H3[2]= -4.0 + 6.0*(loc2+loc3) + r5*(1.0-2.0*loc3) - loc2*(r5+r6);
    data.H3[3] = loc2*(p4+p6);
    data.H3[4] = loc2*(q4-q6);
    data.H3[5]= -loc2*(r6-r4);
    data.H3[6] = p5*(1.0-2.0*loc3) - (p4+p5)*loc2;
    data.H3[7] = q5*(1.0-2.0*loc3)+loc2*(q4-q5);
    data.H3[8]= -2.0+6.0*loc3+r5*(1.0-2.0*loc3)+loc2*(r4-r5);

    data.H4[0] = -t5*(1.0-2.0*loc3) - (t6-t5)*loc2;
    data.H4[1] = 1.0+r5*(1.0-2.0*loc3)-loc2*(r5+r6);
    data.H4[2]= -q5*(1.0-2.0*loc3)+loc2*(q5+q6);
    data.H4[3] = (t4+t6)*loc2;
    data.H4[4] = (r4-r6)*loc2;
    data.H4[5]= -loc2*(q4-q6);
    data.H4[6] = t5*(1.0-2.0*loc3) - (t5+t4)*loc2;
    data.H4[7] = -1.0+r5*(1.0-2.0*loc3)+loc2*(r4-r5);
    data.H4[8]= -q5*(1.0-2.0*loc3)-loc2*(q4-q5);

    double temp = 0.5 / area;
    for(int i =0; i<9; i++)
    {
        data.Bb(0, i) = temp * ( y31*data.H1[i] + y12*data.H3[i]);
        data.Bb(1, i) = temp * (-x31*data.H2[i] - x12*data.H4[i]);
        data.Bb(2, i) = temp * (-x31*data.H1[i] - x12*data.H3[i] + y31*data.H2[i] + y12*data.H4[i] );
    }

    //---------------------------------------------
    // assemble the 2 mamebrane contributions
    // and the bending contribution
    // into the combined B matrix
    for(int nodeid = 0; nodeid < 3; nodeid++)
    {
        int i = nodeid * 3;
        int j = nodeid * 6;

        // membrane map: [0,1,5] <- [0,1,2]

        data.B(0, j  ) = data.L(0, i  ) + data.Qh(0, i  );
        data.B(0, j+1) = data.L(0, i+1) + data.Qh(0, i+1);
        data.B(0, j+5) = data.L(0, i+2) + data.Qh(0, i+2);

        data.B(1, j  ) = data.L(1, i  ) + data.Qh(1, i  );
        data.B(1, j+1) = data.L(1, i+1) + data.Qh(1, i+1);
        data.B(1, j+5) = data.L(1, i+2) + data.Qh(1, i+2);

        data.B(2, j  ) = data.L(2, i  ) + data.Qh(2, i  );
        data.B(2, j+1) = data.L(2, i+1) + data.Qh(2, i+1);
        data.B(2, j+5) = data.L(2, i+2) + data.Qh(2, i+2);

        // bending map: [2,3,4] <- [0,1,2]

        data.B(3, j+2) = data.Bb(0, i  );
        data.B(3, j+3) = data.Bb(0, i+1);
        data.B(3, j+4) = data.Bb(0, i+2);

        data.B(4, j+2) = data.Bb(1, i  );
        data.B(4, j+3) = data.Bb(1, i+1);
        data.B(4, j+4) = data.Bb(1, i+2);

        data.B(5, j+2) = data.Bb(2, i  );
        data.B(5, j+3) = data.Bb(2, i+1);
        data.B(5, j+4) = data.Bb(2, i+2);
    }
}

void ShellThinElement3D3N::CalculateBeta0(CalculationData& data)
{
    data.beta0 = 1.0; // to be changed!
}

void ShellThinElement3D3N::CalculateSectionResponse(CalculationData& data)
{
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
    const Matrix & shapeFunctions = GetGeometry().ShapeFunctionsValues(mIntegrationMethod);
    for(int nodeid = 0; nodeid < GetGeometry().PointsNumber(); nodeid++)
        data.N(nodeid) = shapeFunctions(data.gpIndex, nodeid);
#else
    const array_1d<double,3>& loc = data.gpLocations[data.gpIndex];
    data.N(0) = 1.0 - loc[1] - loc[2];
    data.N(1) = loc[1];
    data.N(2) = loc[2];
#endif // !OPT_USES_INTERIOR_GAUSS_POINTS

    ShellCrossSection::Pointer& section = mSections[data.gpIndex];
    data.SectionParameters.SetShapeFunctionsValues( data.N );
    data.SectionParameters.SetMaterialProperties(GetProperties());
    section->CalculateSectionResponse( data.SectionParameters, ConstitutiveLaw::StressMeasure_PK2 );
}

void ShellThinElement3D3N::CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS)
{
    // calculate beta0
    CalculateBeta0( data );

    // calculate the total strain displ. matrix
    CalculateBMatrix( data );

    // compute generalized strains
    noalias( data.generalizedStrains ) = prod( data.B, data.localDisplacements );

    // calculate section response
    CalculateSectionResponse( data );

    // multiply the section tangent matrices and stress resultants by 'dA'
    data.D *= data.dA;
    //******
    Vector3Type& iSig = data.Sig[data.gpIndex];
    iSig(0) = data.generalizedStresses(0);
    iSig(1) = data.generalizedStresses(1);
    iSig(2) = data.generalizedStresses(2);
    //******
    data.generalizedStresses *= data.dA;

    // Add all contributions to the Stiffness Matrix
    noalias( data.BTD ) = prod( trans( data.B ), data.D );
    noalias( LHS ) += prod( data.BTD, data.B );

    // Add all contributions to the residual vector
    noalias( RHS ) -= prod( trans( data.B ), data.generalizedStresses );
}

void ShellThinElement3D3N::ApplyCorrectionToRHS(CalculationData& data, VectorType& RHS)
{
    Vector3Type meanS;
    meanS.clear();
    for(int i = 0; i < 3; i++)
        noalias( meanS ) += data.Sig[i];
    meanS /= 3.0;

    for(int i = 0; i < 3; i++)
    {
        int i1 = i;
        int i2 = i == 2 ? 0 : i+1;

        const Vector3Type& p1 = data.LCS0.Nodes()[i1];
        const Vector3Type& p2 = data.LCS0.Nodes()[i2];
        /*const Vector3Type& s1 = data.Sig[i1];
        const Vector3Type& s2 = data.Sig[i2];*/
        const Vector3Type& s1 = meanS;
        const Vector3Type& s2 = meanS;
        /*const Vector3Type& s1 = data.Sig[i];
        const Vector3Type& s2 = data.Sig[i];*/

        Vector3Type t = p2 - p1;
        Vector3Type z;
        z(0) = 0.0;
        z(1) = 0.0;
        z(2) = 1.0;
        Vector3Type n;
        MathUtils<double>::CrossProduct(n,  t,z);
        n /= MathUtils<double>::Norm3(n);

        double sx, sy;
        sx = s1(0)*n(0) + s1(2)*n(1);
        sy = s1(2)*n(0) + s1(1)*n(1);
        double q1 = std::sqrt(sx*sx + sy*sy);
        sx = s2(0)*n(0) + s2(2)*n(1);
        sy = s2(2)*n(0) + s2(1)*n(1);
        double q2 = std::sqrt(sx*sx + sy*sy);

        double q0 = (q1+q2)/2.0;

        double L = std::sqrt(t(0)*t(0) + t(1)*t(1));

        double m = 1.0/8.0*L*L*q0;

        RHS(i1*6 + 5) -= m;
        RHS(i2*6 + 5) += m;
    }
}

void ShellThinElement3D3N::AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector)
{
    const GeometryType& geom = GetGeometry();
    const SizeType num_gps = GetNumberOfGPs();

    // Get shape functions
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
    const Matrix & N = GetGeometry().ShapeFunctionsValues(mIntegrationMethod);
#else
    Matrix N(3,3);


    for(unsigned int igauss = 0; igauss < num_gps; igauss++)
    {
        const array_1d<double,3>& loc = data.gpLocations[igauss];
        N(igauss,0) = 1.0 - loc[1] - loc[2];
        N(igauss,1) = loc[1];
        N(igauss,2) = loc[2];
    }
#endif // !OPT_USES_INTERIOR_GAUSS_POINTS

    // auxiliary
    array_1d<double, 3> bf;

    // gauss loop to integrate the external force vector
    for(unsigned int igauss = 0; igauss < num_gps; igauss++)
    {
        // get mass per unit area
        double mass_per_unit_area = mSections[igauss]->CalculateMassPerUnitArea(GetProperties());

        // interpolate nodal volume accelerations to this gauss point
        // and obtain the body force vector
        bf.clear();
        for(unsigned int inode = 0; inode < 3; inode++)
        {
            if( geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
                bf += N(igauss,inode) * geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
        }
        bf *= (mass_per_unit_area * data.dA);

        // add it to the RHS vector
        for(unsigned int inode = 0; inode < 3; inode++)
        {
            unsigned int index = inode*6;
            double iN = N(igauss,inode);
            rRightHandSideVector[index + 0] += iN * bf[0];
            rRightHandSideVector[index + 1] += iN * bf[1];
            rRightHandSideVector[index + 2] += iN * bf[2];
        }
    }
}

void ShellThinElement3D3N::CalculateAll(MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    // Resize the Left Hand Side if necessary,
    // and initialize it to Zero

    if((rLeftHandSideMatrix.size1() != 18) || (rLeftHandSideMatrix.size2() != 18))
        rLeftHandSideMatrix.resize(18, 18, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(18, 18);

    // Resize the Right Hand Side if necessary,
    // and initialize it to Zero

    if(rRightHandSideVector.size() != 18)
        rRightHandSideVector.resize(18, false);
    noalias(rRightHandSideVector) = ZeroVector(18);

    // Initialize common calculation variables

    CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
    data.CalculateLHS = CalculateStiffnessMatrixFlag;
    data.CalculateRHS = CalculateResidualVectorFlag;
    InitializeCalculationData(data);

    // Gauss Loop.
    for(SizeType i = 0; i < GetNumberOfGPs(); ++i)
    {
        data.gpIndex = i;
        CalculateGaussPointContribution(data, rLeftHandSideMatrix, rRightHandSideVector);
    }

    //ApplyCorrectionToRHS(data, rRightHandSideVector);

    // Let the CoordinateTransformation finalize the calculation.
    // This will handle the transformation of the local matrices/vectors to
    // the global coordinate system.

    mpCoordinateTransformation->FinalizeCalculations(data.LCS,
            data.globalDisplacements,
            data.localDisplacements,
            rLeftHandSideMatrix,
            rRightHandSideVector,
            CalculateResidualVectorFlag,
            CalculateStiffnessMatrixFlag);

    // Add body forces contributions. This doesn't depend on the coordinate system

    AddBodyForces(data, rRightHandSideVector);
}

bool ShellThinElement3D3N::TryCalculateOnIntegrationPoints_MaterialOrientation(const Variable<array_1d<double,3> >& rVariable,
        std::vector<array_1d<double,3> >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    // Check the required output

    int ijob = 0;
    if(rVariable == MATERIAL_ORIENTATION_DX)
        ijob = 1;
    else if(rVariable == MATERIAL_ORIENTATION_DY)
        ijob = 2;
    else if(rVariable == MATERIAL_ORIENTATION_DZ)
        ijob = 3;

    // quick return

    if(ijob == 0) return false;

    const SizeType num_gps = GetNumberOfGPs();

    // resize output
    if(rValues.size() != num_gps)
        rValues.resize(num_gps);

    // Compute the local coordinate system.

    ShellT3_LocalCoordinateSystem localCoordinateSystem(
        mpCoordinateTransformation->CreateLocalCoordinateSystem() );

    Vector3Type eZ = localCoordinateSystem.Vz();

    // Gauss Loop
    if(ijob == 1)
    {
        Vector3Type eX = localCoordinateSystem.Vx();
        for(SizeType i = 0; i < num_gps; i++)
        {
            QuaternionType q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mSections[i]->GetOrientationAngle());
            q.RotateVector3(eX, rValues[i]);
        }
    }
    else if(ijob == 2)
    {
        Vector3Type eY = localCoordinateSystem.Vy();
        for(SizeType i = 0; i < num_gps; i++)
        {
            QuaternionType q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mSections[i]->GetOrientationAngle());
            q.RotateVector3(eY, rValues[i]);
        }
    }
    else if(ijob == 3)
    {
        for(SizeType i = 0; i < num_gps; i++)
        {
            noalias( rValues[i] ) = eZ;
        }
    }

    return true;
}

bool ShellThinElement3D3N::TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    // Check the required output

    int ijob = 0;
    bool bGlobal = false;
    CheckGeneralizedStressOrStrainOutput(rVariable, ijob, bGlobal);

    // quick return

    if(ijob == 0) return false;

    const SizeType num_gps = GetNumberOfGPs();

    // resize output
    if(rValues.size() != num_gps)
        rValues.resize(num_gps);

    // Just to store the rotation matrix for visualization purposes

    Matrix R(mStrainSize, mStrainSize);
    Matrix aux33(3, 3);

    // Initialize common calculation variables

    CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
    data.CalculateLHS = false;
    data.CalculateRHS = true;
    InitializeCalculationData(data);

    // Gauss Loop.

    for(SizeType i = 0; i < num_gps; i++)
    {
        // set the current integration point index
        data.gpIndex = i;
        ShellCrossSection::Pointer& section = mSections[i];

        // calculate beta0
        CalculateBeta0( data );

        // calculate the total strain displ. matrix
        CalculateBMatrix( data );

        // compute generalized strains
        noalias( data.generalizedStrains ) = prod( data.B, data.localDisplacements );

        // calculate section response
        if (ijob > 2)
		{
			if (ijob > 7)
			{
				//Calculate lamina stresses
				CalculateLaminaStrains(data);
				CalculateLaminaStresses(data);
			}
			else
			{
				// calculate force resultants
				CalculateSectionResponse(data);

				if (ijob > 4)
				{
					// Compute stresses
					CalculateStressesFromForceResultants(data.generalizedStresses,
						section->GetThickness(GetProperties()));
				}
			}
		}

        // adjust output
        DecimalCorrection( data.generalizedStrains );
        DecimalCorrection( data.generalizedStresses );

        // store the results, but first rotate them back to the section coordinate system.
        // we want to visualize the results in that system not in the element one!
        if(section->GetOrientationAngle() != 0.0 && !bGlobal)
        {
			if (ijob > 7)
			{
				section->GetRotationMatrixForGeneralizedStresses(-(section->GetOrientationAngle()), R);
				for (unsigned int i = 0; i < data.rlaminateStresses.size(); i++)
				{
					data.rlaminateStresses[i] = prod(R, data.rlaminateStresses[i]);
				}

				section->GetRotationMatrixForGeneralizedStrains(-(section->GetOrientationAngle()), R);
				for (unsigned int i = 0; i < data.rlaminateStrains.size(); i++)
				{
					data.rlaminateStrains[i] = prod(R, data.rlaminateStrains[i]);
				}
			}
			else if (ijob > 2)
            {
                section->GetRotationMatrixForGeneralizedStresses( -(section->GetOrientationAngle()), R );
                data.generalizedStresses = prod( R, data.generalizedStresses );
            }
            else
            {
                section->GetRotationMatrixForGeneralizedStrains( -(section->GetOrientationAngle()), R );
                data.generalizedStrains = prod( R, data.generalizedStrains );
            }
        }

        // save results
        Matrix & iValue = rValues[i];
        if(iValue.size1() != 3 || iValue.size2() != 3)
            iValue.resize(3, 3, false);

        if(ijob == 1) // strains
        {
            iValue(0, 0) = data.generalizedStrains(0);
            iValue(1, 1) = data.generalizedStrains(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(2);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        }
        else if(ijob == 2) // curvatures
        {
            iValue(0, 0) = data.generalizedStrains(3);
            iValue(1, 1) = data.generalizedStrains(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(5);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        }
        else if(ijob == 3) // forces
        {
            iValue(0, 0) = data.generalizedStresses(0);
            iValue(1, 1) = data.generalizedStresses(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(2);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        }
        else if(ijob == 4) // moments
        {
            iValue(0, 0) = data.generalizedStresses(3);
            iValue(1, 1) = data.generalizedStresses(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(5);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        }
		else if (ijob == 5) // SHELL_STRESS_TOP_SURFACE
		{
			iValue(0, 0) = data.generalizedStresses(0) +
				data.generalizedStresses(3);
			iValue(1, 1) = data.generalizedStresses(1) +
				data.generalizedStresses(4);
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) = data.generalizedStresses[2] +
				data.generalizedStresses[5];
			iValue(0, 2) = iValue(2, 0) = 0.0;
			iValue(1, 2) = iValue(2, 1) = 0.0;
		}
		else if (ijob == 6) // SHELL_STRESS_MIDDLE_SURFACE
		{
			iValue(0, 0) = data.generalizedStresses(0);
			iValue(1, 1) = data.generalizedStresses(1);
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) = data.generalizedStresses[2];
			iValue(0, 2) = iValue(2, 0) = 0.0;
			iValue(1, 2) = iValue(2, 1) = 0.0;
		}
		else if (ijob == 7) // SHELL_STRESS_BOTTOM_SURFACE
		{
			iValue(0, 0) = data.generalizedStresses(0) -
				data.generalizedStresses(3);
			iValue(1, 1) = data.generalizedStresses(1) -
				data.generalizedStresses(4);
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) = data.generalizedStresses[2] -
				data.generalizedStresses[5];
			iValue(0, 2) = iValue(2, 0) = 0.0;
			iValue(1, 2) = iValue(2, 1) = 0.0;
		}
		else if (ijob == 8) // SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE
		{
			iValue(0, 0) =
				data.rlaminateStresses[data.rlaminateStresses.size() - 1][0];
			iValue(1, 1) =
				data.rlaminateStresses[data.rlaminateStresses.size() - 1][1];
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) =
				data.rlaminateStresses[data.rlaminateStresses.size() - 1][2];
			iValue(0, 2) = iValue(2, 0) = 0.0;
			iValue(1, 2) = iValue(2, 1) = 0.0;
		}
		else if (ijob == 9) // SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE
		{
			iValue(0, 0) = data.rlaminateStresses[0][0];
			iValue(1, 1) = data.rlaminateStresses[0][1];
			iValue(2, 2) = 0.0;
			iValue(0, 1) = iValue(1, 0) = data.rlaminateStresses[0][2];
			iValue(0, 2) = iValue(2, 0) = 0.0;
			iValue(1, 2) = iValue(2, 1) = 0.0;
		}

        // if requested, rotate the results in the global coordinate system
        if(bGlobal)
        {
            const Matrix& RG = data.LCS.Orientation();
            noalias( aux33 ) = prod( trans( RG ), iValue );
            noalias( iValue ) = prod( aux33, RG );
        }
    }

    OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);

    return true;
}

ShellCrossSection::SectionBehaviorType ShellThinElement3D3N::GetSectionBehavior()
{
    return ShellCrossSection::Thin;
}

// =====================================================================================
//
// Class ShellThinElement3D3N - Serialization
//
// =====================================================================================

void ShellThinElement3D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  BaseShellElement );
    rSerializer.save("CTr", mpCoordinateTransformation);
}

void ShellThinElement3D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  BaseShellElement );
    rSerializer.load("CTr", mpCoordinateTransformation);
}

}
