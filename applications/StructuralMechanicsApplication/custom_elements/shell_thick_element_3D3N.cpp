// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//       contact:    A.Winterstein [at] tum.de
//

#include "shell_thick_element_3D3N.hpp"
#include "custom_utilities/shellt3_corotational_coordinate_transformation.hpp"

#include "custom_utilities/shell_utilities.h"

#include <string>
#include <iomanip>


/*
Element overview:---------------------------------------------------------------
This element represents a 3-node Shell element
based on the Discrete Shear Gap theory (DSG) by Bletzinger.
This element is formulated for small strains,
but can be used in Geometrically nonlinear problems
involving large displacements and rotations
using a Corotational Coordinate Transformation.
Material nonlinearity is handled by means of the cross section object.

Shell formulation references:---------------------------------------------------
1.    Bletzinger, K.U., Bischoff, M. and Ramm, E., 2000. A unified approach for
    shear-locking-free triangular and rectangular shell finite elements.
    Computers & Structures, 75(3), pp.321-334.
2.    Rama, G.,  Marinkovic, D.,  Zehn, M., 2016. Efficient co-rotational
    3-node shell element. American Journal of Engineering and Applied Sciences,
    Volume 9, Issue 2, Pages 420-431.
*/

namespace Kratos
{
// =========================================================================
//
// Definitions
//
// =========================================================================

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

// =========================================================================
//
// CalculationData
//
// =========================================================================

ShellThickElement3D3N::CalculationData::CalculationData
(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
    const ProcessInfo& rCurrentProcessInfo)
    : LCS0(pCoordinateTransformation->CreateReferenceCoordinateSystem())
    , LCS(pCoordinateTransformation->CreateLocalCoordinateSystem())
    , CurrentProcessInfo(rCurrentProcessInfo)

{
}

// =========================================================================
//
// Class ShellThickElement3D3N
//
// =========================================================================

ShellThickElement3D3N::ShellThickElement3D3N(IndexType NewId,
    GeometryType::Pointer pGeometry,
    bool NLGeom)
    : BaseShellElement(NewId, pGeometry)
    , mpCoordinateTransformation(NLGeom ?
        new ShellT3_CorotationalCoordinateTransformation(pGeometry) :
        new ShellT3_CoordinateTransformation(pGeometry))
{
}

ShellThickElement3D3N::ShellThickElement3D3N(IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties,
    bool NLGeom)
    : BaseShellElement(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation(NLGeom ?
        new ShellT3_CorotationalCoordinateTransformation(pGeometry) :
        new ShellT3_CoordinateTransformation(pGeometry))
{
}

ShellThickElement3D3N::ShellThickElement3D3N(IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties,
    CoordinateTransformationBasePointerType pCoordinateTransformation)
    : BaseShellElement(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation(pCoordinateTransformation)
{
}

ShellThickElement3D3N::~ShellThickElement3D3N()
{
}

Element::Pointer ShellThickElement3D3N::Create(IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
    return Kratos::make_shared< ShellThickElement3D3N >(NewId, newGeom,
        pProperties, mpCoordinateTransformation->Create(newGeom));
}

Element::Pointer ShellThickElement3D3N::Create(IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared< ShellThickElement3D3N >(NewId, pGeom,
        pProperties, mpCoordinateTransformation->Create(pGeom));
}

void ShellThickElement3D3N::Initialize()
{
    KRATOS_TRY

    const int points_number = GetGeometry().PointsNumber();

    KRATOS_ERROR_IF_NOT(points_number == 3) <<"ShellThickElement3D3N - Wrong number of nodes"
        << points_number << std::endl;

    BaseShellElement::Initialize();

    mpCoordinateTransformation->Initialize();

    this->SetupOrientationAngles();

    KRATOS_CATCH("")
}

void ShellThickElement3D3N::InitializeNonLinearIteration
(ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->InitializeNonLinearIteration(rCurrentProcessInfo);

    BaseInitializeNonLinearIteration(rCurrentProcessInfo);
}

void ShellThickElement3D3N::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->FinalizeNonLinearIteration(rCurrentProcessInfo);

    BaseFinalizeNonLinearIteration(rCurrentProcessInfo);
}

void ShellThickElement3D3N::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    BaseInitializeSolutionStep(rCurrentProcessInfo);

    mpCoordinateTransformation->InitializeSolutionStep(rCurrentProcessInfo);
}

void ShellThickElement3D3N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    BaseFinalizeSolutionStep(rCurrentProcessInfo);

    mpCoordinateTransformation->FinalizeSolutionStep(rCurrentProcessInfo);
}

void ShellThickElement3D3N::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if ((rMassMatrix.size1() != 18) || (rMassMatrix.size2() != 18))
        rMassMatrix.resize(18, 18, false);
    noalias(rMassMatrix) = ZeroMatrix(18, 18);

    // Compute the local coordinate system.
    ShellT3_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem());

    // Average mass per unit area over the whole element
    double av_mass_per_unit_area = 0.0;

    const SizeType num_gps = GetNumberOfGPs();

    for (SizeType i = 0; i < num_gps; i++)
        av_mass_per_unit_area += mSections[i]->CalculateMassPerUnitArea(GetProperties());
    av_mass_per_unit_area /= double(num_gps);

    // Flag for consistent or lumped mass matrix
    bool bconsistent_matrix = false;

    // Consistent mass matrix
    if (bconsistent_matrix)
    {
        // General matrix form as per Felippa plane stress CST eqn 31.27:
        // http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf
        //
        // Density and thickness are averaged over element.

        // Average thickness over the whole element
        double thickness = 0.0;
        for (SizeType i = 0; i < num_gps; i++)
            thickness += mSections[i]->GetThickness(GetProperties());
        thickness /= double(num_gps);

        // Populate mass matrix with integation results
        for (SizeType row = 0; row < 18; row++)
        {
            if (row % 6 < 3)
            {
                // translational entry
                for (SizeType col = 0; col < 3; col++)
                {
                    rMassMatrix(row, 6 * col + row % 6) = 1.0;
                }
            }
            else
            {
                // rotational entry
                for (SizeType col = 0; col < 3; col++)
                {
                    rMassMatrix(row, 6 * col + row % 6) =
                        thickness*thickness / 12.0;
                }
            }

            // Diagonal entry
            rMassMatrix(row, row) *= 2.0;
        }

        rMassMatrix *=
            av_mass_per_unit_area*referenceCoordinateSystem.Area() / 12.0;
    }// Consistent mass matrix
    else
    {
        // Lumped mass matrix

            // lumped area
        double lump_area = referenceCoordinateSystem.Area() / 3.0;

        // loop on nodes
        for (SizeType i = 0; i < 3; i++)
        {
            SizeType index = i * 6;

            double nodal_mass = av_mass_per_unit_area * lump_area;

            // translational mass
            rMassMatrix(index, index) = nodal_mass;
            rMassMatrix(index + 1, index + 1) = nodal_mass;
            rMassMatrix(index + 2, index + 2) = nodal_mass;

            // rotational mass - neglected for the moment...
        }
    }// Lumped mass matrix
}

// =====================================================================================
//
// Class ShellThickElement3D3N - Results on Gauss Points
//
// =====================================================================================

void ShellThickElement3D3N::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_gps = GetNumberOfGPs();
    // resize output
    if (rValues.size() != num_gps)
        rValues.resize(num_gps);

    int caseId = -1;
    if (rVariable == TSAI_WU_RESERVE_FACTOR)
    {
        caseId = 10;
    }
    else if (rVariable == VON_MISES_STRESS ||
        rVariable == VON_MISES_STRESS_TOP_SURFACE ||
        rVariable == VON_MISES_STRESS_MIDDLE_SURFACE ||
        rVariable == VON_MISES_STRESS_BOTTOM_SURFACE)
    {
        caseId = 20;
    }
    else if (rVariable == SHELL_ELEMENT_MEMBRANE_ENERGY ||
        SHELL_ELEMENT_MEMBRANE_ENERGY_FRACTION ||
        SHELL_ELEMENT_BENDING_ENERGY ||
        SHELL_ELEMENT_BENDING_ENERGY_FRACTION ||
        SHELL_ELEMENT_SHEAR_ENERGY ||
        SHELL_ELEMENT_SHEAR_ENERGY_FRACTION)
    {
        caseId = 30;
    }

    if (caseId > 19) // calculate stresses
    {
        // Initialize common calculation variables
        CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
        data.CalculateLHS = true;
        data.CalculateRHS = true;
        InitializeCalculationData(data);

        // Get the current displacements in global coordinate system and
        // transform to reference local system
        ShellT3_LocalCoordinateSystem referenceCoordinateSystem(
            mpCoordinateTransformation->CreateReferenceCoordinateSystem());
        MatrixType Rdisp(18, 18);
        referenceCoordinateSystem.ComputeTotalRotationMatrix(Rdisp);
        data.localDisplacements = prod(Rdisp, data.globalDisplacements);

        // Get strains
        noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

        // Calculate the response of the Cross Section - setup for 1 GP
        data.gpIndex = 0;
        ShellCrossSection::Pointer & section = mSections[0];
        CalculateSectionResponse(data);


        double resultDouble = 0.0;

        if (caseId == 30)
        {
            // Energy calcs - these haven't been verified or tested yet.
            CalculateShellElementEnergy(data, rVariable, resultDouble);
        }
        else if (caseId == 20)
        {
            //Von mises calcs

            // recover stresses
            CalculateStressesFromForceResultants(data.generalizedStresses,
                section->GetThickness(GetProperties()));

            // account for orientation
            if (section->GetOrientationAngle() != 0.0)
            {
                Matrix R(8, 8);
                section->GetRotationMatrixForGeneralizedStresses
                (-(section->GetOrientationAngle()), R);
                data.generalizedStresses = prod(R, data.generalizedStresses);
            }

            CalculateVonMisesStress(data, rVariable, resultDouble);
        }
        else
        {
            KRATOS_ERROR <<
                "Error: ELEMENT ShellThinElement3D4N, METHOD CalculateOnIntegrationPoints(double)"
                << std::endl;
        }

        // loop over gauss points - for output only
        for (unsigned int gauss_point = 0; gauss_point < num_gps; ++gauss_point)
        {
            // store the result calculated
            rValues[gauss_point] = resultDouble;
        }
    }
    else if (rVariable == TSAI_WU_RESERVE_FACTOR)
    {
        // resize output
        if (rValues.size() != num_gps)
            rValues.resize(num_gps);

        // Initialize common calculation variables
        CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
        data.CalculateLHS = true;
        data.CalculateRHS = true;
        InitializeCalculationData(data);
        data.gpIndex = 0;

        // Get the current displacements in global coordinate system and
        // transform to reference local system
        ShellT3_LocalCoordinateSystem referenceCoordinateSystem(
            mpCoordinateTransformation->CreateReferenceCoordinateSystem());
        MatrixType Rdisp(18, 18);
        referenceCoordinateSystem.ComputeTotalRotationMatrix(Rdisp);
        data.localDisplacements = prod(Rdisp, data.globalDisplacements);

        // Get strains
        noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

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

        // Define variables
        Matrix R(8, 8);

        // Retrieve ply orientations
        Vector ply_orientation(section->NumberOfPlies());
        section->GetLaminaeOrientation(props, ply_orientation);
        double total_rotation = 0.0;

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

        // Gauss Loop
        for (unsigned int gauss_point = 0; gauss_point < num_gps; ++gauss_point)
        {
            // Output min Tsai-Wu result
            rValues[gauss_point] = min_tsai_wu;

        }// Gauss loop

    } // Tsai wu
    else
    {
        for (SizeType i = 0; i < num_gps; i++)
            mSections[i]->GetValue(rVariable, GetProperties(), rValues[i]);
    }

    OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
}

void ShellThickElement3D3N::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(rVariable, rValues, rCurrentProcessInfo)) return;
}

void ShellThickElement3D3N::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable,
    std::vector<array_1d<double, 3> >& rOutput,
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

void ShellThickElement3D3N::Calculate(const Variable<Matrix>& rVariable, Matrix & Output, const ProcessInfo & rCurrentProcessInfo)
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
// Class ShellThickElement3D3N - Private methods
//
// =====================================================================================

void ShellThickElement3D3N::CalculateStressesFromForceResultants
(VectorType& rstresses, const double& rthickness)
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

void ShellThickElement3D3N::CalculateLaminaStrains(CalculationData& data)
{
    ShellCrossSection::Pointer& section = mSections[data.gpIndex];

    // Get laminate properties
    double thickness = section->GetThickness(GetProperties());
    double z_current = thickness / -2.0; // start from the top of the 1st layer

    // Establish current strains at the midplane
    // (element coordinate system)
    double e_x = data.generalizedStrains[0];
    double e_y = data.generalizedStrains[1];
    double e_xy = data.generalizedStrains[2];    //this is still engineering
                                                //strain (2xtensorial shear)
    double kap_x = data.generalizedStrains[3];
    double kap_y = data.generalizedStrains[4];
    double kap_xy = data.generalizedStrains[5];    //this is still engineering
                                                //strain (2xtensorial shear)

    // Get ply thicknesses
    Vector ply_thicknesses = Vector(section->NumberOfPlies(), 0.0);
    section->GetPlyThicknesses(GetProperties(), ply_thicknesses);

    // Resize output vector. 2 Surfaces for each ply
    data.rlaminateStrains.resize(2 * section->NumberOfPlies());
    for (unsigned int i = 0; i < 2 * section->NumberOfPlies(); i++)
    {
        data.rlaminateStrains[i].resize(8, false);
        data.rlaminateStrains[i].clear();
    }

    // Loop over all plies - start from top ply, top surface
    for (unsigned int plyNumber = 0;
        plyNumber < section->NumberOfPlies(); ++plyNumber)
    {
        // Calculate strains at top surface, arranged in columns.
        // (element coordinate system)
        data.rlaminateStrains[2 * plyNumber][0] = e_x + z_current*kap_x;
        data.rlaminateStrains[2 * plyNumber][1] = e_y + z_current*kap_y;
        data.rlaminateStrains[2 * plyNumber][2] = e_xy + z_current*kap_xy;

        if (data.parabolic_composite_transverse_shear_strains)
        {
            // assume parabolic transverse shear strain dist across whole
            // laminate
            data.rlaminateStrains[2 * plyNumber][6] =
                1.5*(1.0 - 4.0 * z_current*z_current / thickness / thickness)*
                data.generalizedStrains[6];
            data.rlaminateStrains[2 * plyNumber][7] =
                1.5*(1.0 - 4.0 * z_current*z_current / thickness / thickness)*
                data.generalizedStrains[7];
        }
        else
        {
            // constant transverse shear strain dist
            data.rlaminateStrains[2 * plyNumber][6] = data.generalizedStrains[6];
            data.rlaminateStrains[2 * plyNumber][7] = data.generalizedStrains[7];
        }


        // Move to bottom surface of current layer
        z_current += ply_thicknesses[plyNumber];

        // Calculate strains at bottom surface, arranged in columns
        // (element coordinate system)
        data.rlaminateStrains[2 * plyNumber + 1][0] = e_x + z_current*kap_x;
        data.rlaminateStrains[2 * plyNumber + 1][1] = e_y + z_current*kap_y;
        data.rlaminateStrains[2 * plyNumber + 1][2] = e_xy + z_current*kap_xy;

        if (data.parabolic_composite_transverse_shear_strains)
        {
            data.rlaminateStrains[2 * plyNumber + 1][6] =
                1.5*(1.0 - 4.0 * z_current*z_current / thickness / thickness)*
                data.generalizedStrains[6];
            data.rlaminateStrains[2 * plyNumber + 1][7] =
                1.5*(1.0 - 4.0 * z_current*z_current / thickness / thickness)*
                data.generalizedStrains[7];
        }
        else
        {
            data.rlaminateStrains[2 * plyNumber+1][6] = data.generalizedStrains[6];
            data.rlaminateStrains[2 * plyNumber+1][7] = data.generalizedStrains[7];
        }
    }
}

void ShellThickElement3D3N::CalculateLaminaStresses(CalculationData& data)
{
    ShellCrossSection::Pointer& section = mSections[data.gpIndex];

    // Setup flag to compute ply constitutive matrices
    // (units [Pa] and rotated to element orientation)
    section->SetupGetPlyConstitutiveMatrices();
    CalculateSectionResponse(data);

    // Resize output vector. 2 Surfaces for each ply
    data.rlaminateStresses.resize(2 * section->NumberOfPlies());
    for (unsigned int i = 0; i < 2 * section->NumberOfPlies(); i++)
    {
        data.rlaminateStresses[i].resize(8, false);
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

double ShellThickElement3D3N::CalculateTsaiWuPlaneStress(const CalculationData & data, const Matrix& rLamina_Strengths, const unsigned int& rPly)
{
    // Incoming lamina strengths are organized as follows:
    // Refer to 'shell_cross_section.cpp' for details.
    //
    //    |    T1,        C1,        T2    |
    //    |    C2,        S12,    S13    |
    //    |   S23        0        0    |

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
    Matrix F_ij = Matrix(5, 5, 0.0);
    F_ij.clear();
    F_ij(0, 0) = 1.0 / rLamina_Strengths(0, 0) / rLamina_Strengths(0, 1);    // 11
    F_ij(1, 1) = 1.0 / rLamina_Strengths(0, 2) / rLamina_Strengths(1, 0);    // 22
    F_ij(2, 2) = 1.0 / rLamina_Strengths(1, 1) / rLamina_Strengths(1, 1);    // 12
    F_ij(0, 1) = F_ij(1, 0) = -0.5 / std::sqrt(rLamina_Strengths(0, 0)*rLamina_Strengths(0, 1)*rLamina_Strengths(0, 2)*rLamina_Strengths(1, 0));

    if (disable_in_plane_interaction)
    {
        F_ij(0, 1) = F_ij(1, 0) = 0.0;
    }

    // Third, addditional transverse shear terms
    F_ij(3, 3) = 1.0 / rLamina_Strengths(1, 2) / rLamina_Strengths(1, 2);    // 13
    F_ij(4, 4) = 1.0 / rLamina_Strengths(2, 0) / rLamina_Strengths(2, 0);    // 23

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
    var_a += F_ij(3, 3)*data.rlaminateStresses[2 * rPly][6] * data.rlaminateStresses[2 * rPly][6]; // Transverse shear 13
    var_a += F_ij(4, 4)*data.rlaminateStresses[2 * rPly][7] * data.rlaminateStresses[2 * rPly][7]; // Transverse shear 23

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
    var_a += F_ij(3, 3)*data.rlaminateStresses[2 * rPly + 1][6] * data.rlaminateStresses[2 * rPly + 1][6]; // Transverse shear 13
    var_a += F_ij(4, 4)*data.rlaminateStresses[2 * rPly + 1][7] * data.rlaminateStresses[2 * rPly + 1][7]; // Transverse shear 23

    double tsai_reserve_factor_bottom = (-1.0*var_b + std::sqrt(var_b*var_b + 4.0 * var_a)) / 2.0 / var_a;

    // Return min of both surfaces as the result for the whole ply
    return std::min(tsai_reserve_factor_bottom, tsai_reserve_factor_top);
}

void ShellThickElement3D3N::CalculateVonMisesStress(const CalculationData & data, const Variable<double>& rVariable, double & rVon_Mises_Result)
{
    // calc von mises stresses at top, mid and bottom surfaces for
    // thick shell
    double sxx, syy, sxy, sxz, syz;
    double von_mises_top, von_mises_mid, von_mises_bottom;
    // top surface: membrane and +bending contributions
    //                (no transverse shear)
    sxx = data.generalizedStresses[0] + data.generalizedStresses[3];
    syy = data.generalizedStresses[1] + data.generalizedStresses[4];
    sxy = data.generalizedStresses[2] + data.generalizedStresses[5];
    von_mises_top = sxx*sxx - sxx*syy + syy*syy + 3.0*sxy*sxy;

    // mid surface: membrane and transverse shear contributions
    //                (no bending)
    sxx = data.generalizedStresses[0];
    syy = data.generalizedStresses[1];
    sxy = data.generalizedStresses[2];
    sxz = data.generalizedStresses[6];
    syz = data.generalizedStresses[7];
    von_mises_mid = sxx*sxx - sxx*syy + syy*syy +
        3.0*(sxy*sxy + sxz*sxz + syz*syz);

    // bottom surface:    membrane and -bending contributions
    //                    (no transverse shear)
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

void ShellThickElement3D3N::CalculateShellElementEnergy(const CalculationData & data, const Variable<double>& rVariable, double & rEnergy_Result)
{
    // Energy calcs - these haven't been verified or tested yet.

    // At each Gauss Point the energy of that Gauss Point's weighted area
    // dA*w_i is output. This means that the total energy of the element is
    // the sum of the Gauss Point energies. Accordingly, the total energy of
    // the system is the sum of Gauss Point energies over all elements.

    bool is_fraction_calc = false;
    double totalEnergy = 1.0;

    if (rVariable == SHELL_ELEMENT_MEMBRANE_ENERGY_FRACTION ||
        rVariable == SHELL_ELEMENT_BENDING_ENERGY_FRACTION ||
        rVariable == SHELL_ELEMENT_SHEAR_ENERGY_FRACTION)
    {
        // need to calculate total energy over current dA first
        totalEnergy = inner_prod(data.generalizedStresses, data.generalizedStrains)*data.TotalArea/3.0;
        is_fraction_calc = true;
    }

    if (rVariable == SHELL_ELEMENT_MEMBRANE_ENERGY || rVariable == SHELL_ELEMENT_MEMBRANE_ENERGY_FRACTION)
    {
        for (SizeType i = 0; i < 3; i++)
        {
            rEnergy_Result += data.generalizedStresses[i] * data.generalizedStrains[i] * data.TotalArea / 3.0;
        }

        if (is_fraction_calc)
        {
            rEnergy_Result /= totalEnergy;
        }
    }
    else if (rVariable == SHELL_ELEMENT_BENDING_ENERGY || rVariable == SHELL_ELEMENT_BENDING_ENERGY_FRACTION)
    {
        for (SizeType i = 3; i < 6; i++)
        {
            rEnergy_Result += data.generalizedStresses[i] * data.generalizedStrains[i] * data.TotalArea / 3.0;
        }

        if (is_fraction_calc)
        {
            rEnergy_Result /= totalEnergy;
        }
    }
    else if (rVariable == SHELL_ELEMENT_SHEAR_ENERGY || rVariable == SHELL_ELEMENT_SHEAR_ENERGY_FRACTION)
    {
        for (SizeType i = 6; i < 8; i++)
        {
            rEnergy_Result += data.generalizedStresses[i] * data.generalizedStrains[i] * data.TotalArea / 3.0;
        }

        if (is_fraction_calc)
        {
            rEnergy_Result /= totalEnergy;
        }
    }
}

void ShellThickElement3D3N::CheckGeneralizedStressOrStrainOutput(const Variable<Matrix>& rVariable, int & ijob, bool & bGlobal)
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
    else if (rVariable == SHELL_ORTHOTROPIC_4PLY_THROUGH_THICKNESS)
    {
        // TESTING VARIABLE
        ijob = 99;
    }
}

void ShellThickElement3D3N::DecimalCorrection(Vector& a)
{
    double norm = norm_2(a);
    double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
    for (SizeType i = 0; i < a.size(); i++)
        if (std::abs(a(i)) < tolerance)
            a(i) = 0.0;
}

void ShellThickElement3D3N::SetupOrientationAngles()
{
    if (this->Has(MATERIAL_ORIENTATION_ANGLE))
    {
        for (CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it)
        (*it)->SetOrientationAngle(this->GetValue(MATERIAL_ORIENTATION_ANGLE));
    }
    else
    {
        ShellT3_LocalCoordinateSystem lcs(mpCoordinateTransformation->CreateReferenceCoordinateSystem());

        Vector3Type normal;
        noalias(normal) = lcs.Vz();

        Vector3Type dZ;
        dZ(0) = 0.0;
        dZ(1) = 0.0;
        dZ(2) = 1.0; // for the moment let's take this. But the user can specify its own triad! TODO

        Vector3Type dirX;
        MathUtils<double>::CrossProduct(dirX, dZ, normal);

        // try to normalize the x vector. if it is near zero it means that we need
        // to choose a default one.
        double dirX_norm = dirX(0)*dirX(0) + dirX(1)*dirX(1) + dirX(2)*dirX(2);
        if (dirX_norm < 1.0E-12)
        {
            dirX(0) = 1.0;
            dirX(1) = 0.0;
            dirX(2) = 0.0;
        }
        else if (dirX_norm != 1.0)
        {
            dirX_norm = std::sqrt(dirX_norm);
            dirX /= dirX_norm;
        }

        Vector3Type elem_dirX = lcs.Vx();

        // now calculate the angle between the element x direction and the material x direction.
        Vector3Type& a = elem_dirX;
        Vector3Type& b = dirX;
        double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
        if (a_dot_b < -1.0) a_dot_b = -1.0;
        if (a_dot_b > 1.0) a_dot_b = 1.0;
        double angle = std::acos(a_dot_b);

        // if they are not counter-clock-wise, let's change the sign of the angle
        if (angle != 0.0)
        {
            const MatrixType& R = lcs.Orientation();
            if (dirX(0)*R(1, 0) + dirX(1)*R(1, 1) + dirX(2)*R(1, 2) < 0.0)
                angle = -angle;
        }

        for (CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it)
            (*it)->SetOrientationAngle(angle);
    }
}

void ShellThickElement3D3N::CalculateSectionResponse(CalculationData& data)
{
    const array_1d<double, 3>& loc = data.gpLocations[0];
    data.N(0) = 1.0 - loc[1] - loc[2];
    data.N(1) = loc[1];
    data.N(2) = loc[2];

    ShellCrossSection::Pointer& section = mSections[0];
    data.SectionParameters.SetShapeFunctionsValues(data.N);
    data.SectionParameters.SetMaterialProperties(GetProperties());

    if (data.ignore_shear_stabilization || data.basicTriCST)
    {
        //remove already added shear stabilization
        data.shearStabilisation = 1.0;
        data.SectionParameters.SetStenbergShearStabilization(data.shearStabilisation);
        std::cout << "Not applying shear stabilisation to shear part of material matrix!" << std::endl;
    }

    section->CalculateSectionResponse(data.SectionParameters, ConstitutiveLaw::StressMeasure_PK2);
}

void ShellThickElement3D3N::InitializeCalculationData(CalculationData& data)
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

    const double y13 = -y31;

    const double A = 0.5*(y21*x13 - x21*y13);
    const double A2 = 2.0*A;

    double h = 0.0;
    for (unsigned int i = 0; i < mSections.size(); i++)
        h += mSections[i]->GetThickness(GetProperties());
    h /= (double)mSections.size();

    data.hMean = h;
    data.TotalArea = A;
    data.dA = A;

    // create the integration point locations
    if (data.gpLocations.size() != 0) data.gpLocations.clear();

    data.gpLocations.resize(GetNumberOfGPs());
#ifdef OPT_1_POINT_INTEGRATION
    array_1d<double, 3>& gp0 = data.gpLocations[0];
    gp0[0] = 1.0 / 3.0;
    gp0[1] = 1.0 / 3.0;
    gp0[2] = 1.0 / 3.0;
#else
    array_1d<double, 3>& gp0 = data.gpLocations[0];
    array_1d<double, 3>& gp1 = data.gpLocations[1];
    array_1d<double, 3>& gp2 = data.gpLocations[2];
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
    gp0[0] = 1.0 / 6.0;
    gp0[1] = 1.0 / 6.0;
    gp0[2] = 2.0 / 3.0;
    gp1[0] = 2.0 / 3.0;
    gp1[1] = 1.0 / 6.0;
    gp1[2] = 1.0 / 6.0;
    gp2[0] = 1.0 / 6.0;
    gp2[1] = 2.0 / 3.0;
    gp2[2] = 1.0 / 6.0;
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
    data.dNxy(0, 0) = (y13 - y12) / A2;
    data.dNxy(0, 1) = (x12 - x13) / A2;
    data.dNxy(1, 0) = -y13 / A2;
    data.dNxy(1, 1) = x13 / A2;
    data.dNxy(2, 0) = y12 / A2;
    data.dNxy(2, 1) = -x12 / A2;

    // ---------------------------------------------------------------------
    // Total Formulation:
    //        as per Efficient Co-Rotational 3-Node Shell Element paper (2016)
    // ---------------------------------------------------------------------

    data.B.clear();

    //Membrane components
    //
    //node 1
    data.B(0, 0) = y23;
    data.B(1, 1) = x32;
    data.B(2, 0) = x32;
    data.B(2, 1) = y23;

    //node 2
    data.B(0, 6) = y31;
    data.B(1, 7) = x13;
    data.B(2, 6) = x13;
    data.B(2, 7) = y31;

    //node 3
    data.B(0, 12) = y12;
    data.B(1, 13) = x21;
    data.B(2, 12) = x21;
    data.B(2, 13) = y12;

    //Bending components
    //
    //node 1
    data.B(3, 4) = y23;
    data.B(4, 3) = -1.0 * x32;
    data.B(5, 3) = -1.0 * y23;
    data.B(5, 4) = x32;

    //node 2
    data.B(3, 10) = y31;
    data.B(4, 9) = -1.0 * x13;
    data.B(5, 9) = -1.0 * y31;
    data.B(5, 10) = x13;

    //node 3
    data.B(3, 16) = y12;
    data.B(4, 15) = -1.0 * x21;
    data.B(5, 15) = -1.0 * y12;
    data.B(5, 16) = x21;

    //Shear components
    //
    if (data.specialDSGc3)
    {
        // Don't do anything here.
        // The shear contribution will be added with 3GPs later on.
    }
    else if (data.smoothedDSG)
    {
        // Use smoothed DSG formulation according to [Nguyen-Thoi et al., 2013]
        std::cout << "Using smoothed DSG" << std::endl;
        CalculateSmoothedDSGBMatrix(data);
    }
    else if (data.basicTriCST == false)
    {
        // Use DSG method as per Bletzinger (2000)
        // Entries denoted "altered by -1" are multiplied by -1

        const double a = x21;
        const double b = y21;
        const double c = y31;
        const double d = x31;

        //node 1
        data.B(6, 2) = b - c;
        data.B(6, 4) = A;

        data.B(7, 2) = d - a;
        data.B(7, 3) = -1.0*A;


        //node 2
        data.B(6, 8) = c;
        data.B(6, 9) = -1.0 * b*c / 2.0;
        data.B(6, 10) = a*c / 2.0;

        data.B(7, 8) = -1.0 * d;
        data.B(7, 9) = b*d / 2.0;
        data.B(7, 10) = -1.0 * a*d / 2.0;

        //node 3
        data.B(6, 14) = -1.0 * b;
        data.B(6, 15) = b*c / 2.0;
        data.B(6, 16) = -b*d / 2.0;

        data.B(7, 14) = a;
        data.B(7, 15) = -1.0 * a*c / 2.0;
        data.B(7, 16) = a*d / 2.0;
    }
    else
    {
        // Basic CST displacement derived shear
        // strain displacement matrix.
        // Only for testing!
        std::cout << "Using basic CST shear formulation!" << std::endl;
        const Matrix & shapeFunctions =
            GetGeometry().ShapeFunctionsValues(mIntegrationMethod);

        //node 1
        data.B(6, 2) = y23;
        data.B(6, 3) = shapeFunctions(0, 0);

        data.B(7, 2) = x32;
        data.B(7, 4) = shapeFunctions(0, 0);

        //node 2
        data.B(6, 8) = y31;
        data.B(6, 9) = shapeFunctions(0, 1);

        data.B(7, 8) = x13;
        data.B(7, 10) = shapeFunctions(0, 1);

        //node 3
        data.B(6, 14) = y12;
        data.B(6, 15) = shapeFunctions(0, 2);

        data.B(7, 14) = x21;
        data.B(7, 16) = shapeFunctions(0, 2);
    }

    //Final multiplication
    data.B /= (A2);

    //determine longest side length (eqn 22) - for shear correction
    Vector P12 = Vector(data.LCS0.P1() - data.LCS0.P2());
    Vector P13 = Vector(data.LCS0.P1() - data.LCS0.P3());
    Vector P23 = Vector(data.LCS0.P2() - data.LCS0.P3());
    data.h_e = std::sqrt(inner_prod(P12, P12));
    double edge_length = std::sqrt(inner_prod(P13, P13));
    if (edge_length > data.h_e) { data.h_e = edge_length; }
    edge_length = std::sqrt(inner_prod(P23, P23));
    if (edge_length > data.h_e) { data.h_e = edge_length; }

    // Write Stenberg shear stabilisation coefficient
    data.shearStabilisation = (data.hMean*data.hMean)
        / (data.hMean*data.hMean + data.alpha*data.h_e*data.h_e);

    //--------------------------------------
    // Calculate material matrices
    //

    //allocate and setup
    data.SectionParameters.SetElementGeometry(GetGeometry());
    data.SectionParameters.SetMaterialProperties(GetProperties());
    data.SectionParameters.SetProcessInfo(data.CurrentProcessInfo);
    data.SectionParameters.SetGeneralizedStrainVector(data.generalizedStrains);
    data.SectionParameters.SetGeneralizedStressVector(data.generalizedStresses);
    data.SectionParameters.SetConstitutiveMatrix(data.D);
    data.SectionParameters.SetShapeFunctionsDerivatives(data.dNxy);
    data.SectionParameters.SetStenbergShearStabilization(data.shearStabilisation);
    Flags& options = data.SectionParameters.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, data.CalculateRHS);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,
        data.CalculateLHS);

    //--------------------------------------
    // calculate the displacement vector
    // in global and local coordinate systems

    data.globalDisplacements.clear();
    GetValuesVector(data.globalDisplacements);

    data.localDisplacements =
        mpCoordinateTransformation->CalculateLocalDisplacements(
            data.LCS, data.globalDisplacements);
}

void ShellThickElement3D3N::CalculateDSGc3Contribution(CalculationData & data, MatrixType & rLeftHandSideMatrix)
{
    std::cout << "DSGc3" << std::endl;

    // choose formulation options
    bool use_pure_bubble_mode = false;
    bool use_no_bubble_mode = true;
    bool use_reconstructed_shear_gaps = false;
    bool use_original_dsg = false;
    bool use_original_dsg_with_bubble = false;

    // Material matrix data.D is already multiplied with A!!!
    data.D *= 2.0; //now data.D is D*detJ

    // Modified B-matrix
    Matrix BSuper = Matrix(2, 9, 0.0);

    // Get geometry data
    const double x1 = data.LCS0.X1();
    const double y1 = data.LCS0.Y1();
    const double x2 = data.LCS0.X2();
    const double y2 = data.LCS0.Y2();
    const double x3 = data.LCS0.X3();
    const double y3 = data.LCS0.Y3();

    // Geometric quantities
    const double a = x2 - x1;
    const double b = y2 - y1;
    const double c = y3 - y1;
    const double d = x3 - x1;

    // Numeric integration
    int numberGPs = 3;
    double integration_weight = 1.0 / 6.0;
    if (use_original_dsg_with_bubble)
    {
        // use quartic integration
        numberGPs = 7;
    }

    std::vector< array_1d<double, 3> > quarticGPLocations;
    quarticGPLocations.resize(7);
    quarticGPLocations.clear();
    for (SizeType i = 0; i < 7; i++)
    {
        quarticGPLocations[i].clear();
        quarticGPLocations[i][0] = 0.0;
        quarticGPLocations[i][1] = 0.0;
    }
    Vector quarticGPweights = Vector(7, 0.0);

    quarticGPweights[0] = 1.0 / 40.0;

    quarticGPLocations[1][0] = 0.5;
    quarticGPweights[1] = 1.0 / 15.0;

    quarticGPLocations[2][0] = 1.0;
    quarticGPweights[2] = 1.0 / 40.0;

    quarticGPLocations[3][0] = 0.5;
    quarticGPLocations[3][1] = 0.5;
    quarticGPweights[3] = 1.0 / 15.0;

    quarticGPLocations[4][1] = 1.0;
    quarticGPweights[4] = 1.0 / 40.0;

    quarticGPLocations[5][1] = 0.5;
    quarticGPweights[5] = 1.0 / 15.0;

    quarticGPLocations[6][0] = 0.3;
    quarticGPLocations[6][1] = 0.3;
    quarticGPweights[6] = 9.0 / 40.0;


    //  Integration loop
    double loc1, loc2;
    for (int gauss_point = 0; gauss_point < numberGPs; gauss_point++)
    {
        if (use_original_dsg_with_bubble)
        {
            loc1 = quarticGPLocations[gauss_point][0];
            loc2 = quarticGPLocations[gauss_point][1];
        }
        else
        {
            loc1 = data.gpLocations[gauss_point][0];
            loc2 = data.gpLocations[gauss_point][1];
        }

        BSuper.clear();

        if (use_pure_bubble_mode)
        {
            std::cout << "use_pure_bubble_mode" << std::endl;
            // pure bubble mode formulation below with no alterations
            // DOFs are [w1,w2,w3,Rx1,...]
            BSuper(0, 0) = -2.0*b*loc1 + 1.0*b + 2.0*c*loc2 - 1.0*c;
            BSuper(0, 1) = 1.0*b*loc1 - 1.0*c*loc2 + 1.0*c;
            BSuper(0, 2) = 1.0*b*loc1 - 1.0*b - 1.0*c*loc2;
            BSuper(0, 3) = -1.0*b*(1.0*loc1 - 0.5);
            BSuper(0, 4) = 0.5*b*loc1 + 0.5*c*loc2;
            BSuper(0, 5) = -0.5*b*loc1 + 0.5*b + 0.5*c*loc2;
            BSuper(0, 6) = -1.0*c*(1.0*loc2 - 0.5);
            BSuper(0, 7) = 0.5*b*loc1 - 0.5*c*loc2 + 0.5*c;
            BSuper(0, 8) = 0.5*b*loc1 + 0.5*c*loc2;
            BSuper(1, 0) = 2.0*a*loc1 - 1.0*a - 2.0*d*loc2 + 1.0*d;
            BSuper(1, 1) = -1.0*a*loc1 + 1.0*d*loc2 - 1.0*d;
            BSuper(1, 2) = -1.0*a*loc1 + 1.0*a + 1.0*d*loc2;
            BSuper(1, 3) = 1.0*a*(1.0*loc1 - 0.5);
            BSuper(1, 4) = -0.5*a*loc1 - 0.5*d*loc2;
            BSuper(1, 5) = 0.5*a*loc1 - 0.5*a - 0.5*d*loc2;
            BSuper(1, 6) = 1.0*d*(1.0*loc2 - 0.5);
            BSuper(1, 7) = -0.5*a*loc1 + 0.5*d*loc2 - 0.5*d;
            BSuper(1, 8) = -0.5*a*loc1 - 0.5*d*loc2;

            BSuper /= (2.0*data.TotalArea);
        }
        else if (use_no_bubble_mode)
        {
            std::cout << "use_no_bubble_mode" << std::endl;
            // pure no bubble mode formulation below with no alterations
            // DOFs are [w1,w2,w3,Rx1,...]

            BSuper(0, 0) = b - c;
            BSuper(0, 1) = c;
            BSuper(0, 2) = -b;
            BSuper(0, 3) = 0.5*(b - c)*(b*loc1 + c*loc2);
            BSuper(0, 4) = -0.5*b*b * loc1 + 0.5*b*c*loc1 - 0.5*b*c*loc2 - 0.5*b*c + 0.5*c*c * loc2;
            BSuper(0, 5) = 0.5*b*b * loc1 - 0.5*b*c*loc1 + 0.5*b*c*loc2 + 0.5*b*c - 0.5*c*c * loc2;
            BSuper(0, 6) = -0.5*a*b*loc1 - 0.5*a*c*loc2 + 0.5*a*c + 0.5*b*d*loc1 - 0.5*b*d + 0.5*c*d*loc2;
            BSuper(0, 7) = -0.5*a*b*loc1 - 0.5*a*c*loc2 + 0.5*a*c + 0.5*b*d*loc1 + 0.5*c*d*loc2;
            BSuper(0, 8) = 0.5*a*b*loc1 + 0.5*a*c*loc2 - 0.5*b*d*loc1 - 0.5*b*d - 0.5*c*d*loc2;
            BSuper(1, 0) = -a + d;
            BSuper(1, 1) = -d;
            BSuper(1, 2) = a;
            BSuper(1, 3) = -0.5*a*b*loc1 + 0.5*a*c*loc1 - 0.5*a*c - 0.5*b*d*loc2 + 0.5*b*d + 0.5*c*d*loc2;
            BSuper(1, 4) = 0.5*a*b*loc1 - 0.5*a*c*loc1 + 0.5*b*d*loc2 + 0.5*b*d - 0.5*c*d*loc2;
            BSuper(1, 5) = -0.5*a*b*loc1 + 0.5*a*c*loc1 - 0.5*a*c - 0.5*b*d*loc2 + 0.5*c*d*loc2;
            BSuper(1, 6) = 0.5*(a - d)*(a*loc1 + d*loc2);
            BSuper(1, 7) = 0.5*a*a * loc1 - 0.5*a*d*loc1 + 0.5*a*d*loc2 - 0.5*a*d - 0.5*d*d * loc2;
            BSuper(1, 8) = -0.5*a*a * loc1 + 0.5*a*d*loc1 - 0.5*a*d*loc2 + 0.5*a*d + 0.5*d*d * loc2;

            BSuper /= (2.0*data.TotalArea);
        }
        else if (use_reconstructed_shear_gaps)
        {
            std::cout << "use_reconstructed_shear_gaps" << std::endl;
            // The reconstructed shear gap field is the same for
            // bubble and no bubble modes

            BSuper(0, 0) = b - c;
            BSuper(0, 1) = c;
            BSuper(0, 2) = -b;
            BSuper(0, 3) = -0.5*b*loc1 - 0.5*c*loc2;
            BSuper(0, 4) = -0.5*b*c + 0.5*b*loc1 + 0.5*c*loc2;
            BSuper(0, 5) = 0.5*b*c;
            BSuper(0, 6) = 0.5*a*c - 0.5*b*d - 0.5*b*loc1 - 0.5*c*loc2;
            BSuper(0, 7) = 0.5*a*c;
            BSuper(0, 8) = -0.5*b*d + 0.5*b*loc1 + 0.5*c*loc2;
            BSuper(1, 0) = -a + d;
            BSuper(1, 1) = -d;
            BSuper(1, 2) = a;
            BSuper(1, 3) = -0.5*a*c + 0.5*a*loc1 + 0.5*b*d + 0.5*d*loc2;
            BSuper(1, 4) = -0.5*a*loc1 + 0.5*b*d - 0.5*d*loc2;
            BSuper(1, 5) = -0.5*a*c;
            BSuper(1, 6) = 0.5*a*loc1 + 0.5*d*loc2;
            BSuper(1, 7) = -0.5*a*d;
            BSuper(1, 8) = 0.5*a*d - 0.5*a*loc1 - 0.5*d*loc2;

            BSuper /= (2.0*data.TotalArea);
        }
        else if (use_original_dsg)
        {
            std::cout << "use_original_dsg" << std::endl;
            BSuper(0, 0) = 1.0*b - 1.0*c;
            BSuper(0, 1) = c;
            BSuper(0, 2) = -b;
            BSuper(0, 3) = 0;
            BSuper(0, 4) = -0.5*b*c;
            BSuper(0, 5) = 0.5*b*c;
            BSuper(0, 6) = 0.5*a*c - 0.5*b*d;
            BSuper(0, 7) = 0.5*a*c;
            BSuper(0, 8) = -0.5*b*d;
            BSuper(1, 0) = -1.0*a + 1.0*d;
            BSuper(1, 1) = -d;
            BSuper(1, 2) = a;
            BSuper(1, 3) = -0.5*a*c + 0.5*b*d;
            BSuper(1, 4) = 0.5*b*d;
            BSuper(1, 5) = -0.5*a*c;
            BSuper(1, 6) = 0;
            BSuper(1, 7) = -0.5*a*d;
            BSuper(1, 8) = 0.5*a*d;

            BSuper /= (2.0*data.TotalArea);
        }
        else if (use_original_dsg_with_bubble)
        {
            std::cout << "use_original_dsg with bubble" << std::endl;

            BSuper(0, 0) = 1.0*b - 1.0*c;
            BSuper(0, 1) = 1.0*c*(1.0);
            BSuper(0, 2) = -1.0*b*(1.0);
            BSuper(0, 3) = 0.0;
            BSuper(0, 4) = -1.0*b*c*(3.0*loc1*loc1 + 12.0*loc1*loc2 - 3.0*loc1 + 3.0*loc2*loc2 - 3.0*loc2 + 0.5);
            BSuper(0, 5) = 1.0*b*c*(3.0*loc1*loc1 + 12.0*loc1*loc2 - 3.0*loc1 + 3.0*loc2*loc2 - 3.0*loc2 + 0.5);
            BSuper(0, 6) = 0.5*a*c - 0.5*b*d;
            BSuper(0, 7) = 6.0*a*c*loc1*loc2 + 3.0*a*c*loc2*loc2 - 3.0*a*c*loc2 + 0.5*a*c + 3.0*b*d*loc1*loc1 + 6.0*b*d*loc1*loc2 - 3.0*b*d*loc1;
            BSuper(0, 8) = -6.0*a*c*loc1*loc2 - 3.0*a*c*loc2*loc2 + 3.0*a*c*loc2 - 3.0*b*d*loc1*loc1 - 6.0*b*d*loc1*loc2 + 3.0*b*d*loc1 - 0.5*b*d;
            BSuper(1, 0) = - 1.0*a + 1.0*d;
            BSuper(1, 1) = -1.0*d*(1.0);
            BSuper(1, 2) = 1.0*a*(1.0);
            BSuper(1, 3) = -0.5*a*c + 0.5*b*d;
            BSuper(1, 4) = 3.0*a*c*loc1*loc1 + 6.0*a*c*loc1*loc2 - 3.0*a*c*loc1 + 6.0*b*d*loc1*loc2 + 3.0*b*d*loc2*loc2 - 3.0*b*d*loc2 + 0.5*b*d;
            BSuper(1, 5) = -3.0*a*c*loc1*loc1 - 6.0*a*c*loc1*loc2 + 3.0*a*c*loc1 - 0.5*a*c - 6.0*b*d*loc1*loc2 - 3.0*b*d*loc2*loc2 + 3.0*b*d*loc2;
            BSuper(1, 6) = 0.0;
            BSuper(1, 7) = -1.0*a*d*(3.0*loc1*loc1 + 12.0*loc1*loc2 - 3.0*loc1 + 3.0*loc2*loc2 - 3.0*loc2 + 0.5);
            BSuper(1, 8) = 1.0*a*d*(3.0*loc1*loc1 + 12.0*loc1*loc2 - 3.0*loc1 + 3.0*loc2*loc2 - 3.0*loc2 + 0.5);

            BSuper /= (2.0*data.TotalArea);

            integration_weight = quarticGPweights[gauss_point];
        }

        //comparison K matrix, in Bletzinger DOF ordering
        //Matrix temp2 = Matrix(prod(trans(BSuper), shearD));
        //shearK += prod(temp2, BSuper);

        data.B.clear();

        // Transfer from Bletzinger B matrix to Kratos B matrix
        // Dofs from [w1, w2, w3, px1, ...] to [w1, px1, py1, w2, ...]
        for (SizeType node = 0; node < 3; node++)
        {
            data.B(6, 2 + 6 * node) = BSuper(0, node);        // w
            data.B(6, 3 + 6 * node) = BSuper(0, 3 + node);    // phix
            data.B(6, 4 + 6 * node) = BSuper(0, 6 + node);    // phiy

            data.B(7, 2 + 6 * node) = BSuper(1, node);        // w
            data.B(7, 3 + 6 * node) = BSuper(1, 3 + node);    // phix
            data.B(7, 4 + 6 * node) = BSuper(1, 6 + node);    // phiy
        }
        // Add to stiffness matrix
        Matrix temp = Matrix(prod(trans(data.B), data.D*integration_weight));
        rLeftHandSideMatrix += prod(temp, data.B);
    }
}

void ShellThickElement3D3N::CalculateSmoothedDSGBMatrix(CalculationData & data)
{
    // Use smoothed DSG formulation according to [Nguyen-Thoi et al., 2013]
    //

    //meta-triangle centre coords
    const double x0 = (data.LCS0.X1() + data.LCS0.X2() + data.LCS0.X3()) / 3.0;
    const double y0 = (data.LCS0.Y1() + data.LCS0.Y2() + data.LCS0.Y3()) / 3.0;


    // Assemble sub triangle coords
    std::vector<Vector3> subTriangleXCoords = std::vector<Vector3>(3);
    std::vector<Vector3> subTriangleYCoords = std::vector<Vector3>(3);
    for (SizeType i = 0; i < 3; i++)
    {
        subTriangleXCoords[i].clear();
        subTriangleYCoords[i].clear();
        subTriangleXCoords[i][0] = x0;
        subTriangleYCoords[i][0] = y0;
    }
    subTriangleXCoords[0][1] = data.LCS0.X1(); //subtri_1 = 0, 1, 2
    subTriangleXCoords[0][2] = data.LCS0.X2();
    subTriangleYCoords[0][1] = data.LCS0.Y1();
    subTriangleYCoords[0][2] = data.LCS0.Y2();

    subTriangleXCoords[1][1] = data.LCS0.X2(); //subtri_2 = 0, 2, 3
    subTriangleXCoords[1][2] = data.LCS0.X3();
    subTriangleYCoords[1][1] = data.LCS0.Y2();
    subTriangleYCoords[1][2] = data.LCS0.Y3();

    subTriangleXCoords[2][1] = data.LCS0.X3(); //subtri_3 = 0, 3, 1
    subTriangleXCoords[2][2] = data.LCS0.X1();
    subTriangleYCoords[2][1] = data.LCS0.Y3();
    subTriangleYCoords[2][2] = data.LCS0.Y1();


    // The mapping controls how the nodally grouped entries of the virgin
    // sub-triangle B matrices are added to the meta-triangle shear B matrix
    std::vector<Vector3> matrixMapping = std::vector<Vector3>(3);
    matrixMapping[0][0] = 2; // node number, not index (number = index + 1)
    matrixMapping[0][1] = 3;
    matrixMapping[0][2] = 9; // signify no addition contribution needed

    matrixMapping[1][0] = 9;
    matrixMapping[1][1] = 2;
    matrixMapping[1][2] = 3;

    matrixMapping[2][0] = 3;
    matrixMapping[2][1] = 9;
    matrixMapping[2][2] = 2;


    // Setup variables
    Matrix smoothedShearMatrix = Matrix(2, 18, 0.0);
    Matrix virginSubTriangleShearMatrix = Matrix(2, 18, 0.0);
    Matrix convertedSubTriangleShearMatrix = Matrix(2, 18, 0.0);
    double a, b, c, d, subTriangleArea;


    // Loop over all sub triangles
    for (SizeType subTriangle = 0; subTriangle < 3; subTriangle++)
    {
        a = subTriangleXCoords[subTriangle][1] - subTriangleXCoords[subTriangle][0]; //x21
        b = subTriangleYCoords[subTriangle][1] - subTriangleYCoords[subTriangle][0]; //y21
        c = subTriangleYCoords[subTriangle][2] - subTriangleYCoords[subTriangle][0]; //y31
        d = subTriangleXCoords[subTriangle][2] - subTriangleXCoords[subTriangle][0]; //x31
        subTriangleArea = 0.5*(a*c - b*d);


        // Calculate the DSG shear B matrix for only the current subtriangle
        virginSubTriangleShearMatrix.clear();
        CalculateDSGShearBMatrix(virginSubTriangleShearMatrix, a, b, c, d, subTriangleArea);


        // Setup matrix to store shear B converted from subtriangle DOFs to
        // meta-triangle DOFs
        convertedSubTriangleShearMatrix.clear();


        // Express the subtriangle B matrix in terms of the 3 meta-triangle
        // nodes
        for (SizeType metaNode = 0; metaNode < 3; metaNode++)
        {
            for (SizeType row = 0; row < 2; row++)
            {
                for (SizeType col = 0; col < 6; col++)
                {
                    // add in B_0 / 3.0 for all nodes in B matrix
                    convertedSubTriangleShearMatrix(row, metaNode * 6 + col) += virginSubTriangleShearMatrix(row, col) / 3.0;
                }
            }

            if (matrixMapping[subTriangle][metaNode] == 9)
            {
                // centre point entry. do nothing, entries already covered
                // by operation above
            }
            else
            {
                // add in new entries
                for (SizeType row = 0; row < 2; row++)
                {
                    for (SizeType col = 0; col < 6; col++)
                    {
                        convertedSubTriangleShearMatrix(row, metaNode * 6 + col) += virginSubTriangleShearMatrix(row, 6*(matrixMapping[subTriangle][metaNode] - 1) + col);
                    }
                }
            }
        }

        // Add subtriangle contribution to overall meta-triangle B matrix
        smoothedShearMatrix += (convertedSubTriangleShearMatrix * subTriangleArea);
    }


    // Smooth by averaging over area
    smoothedShearMatrix /= data.TotalArea;
    smoothedShearMatrix *= (2.0*data.TotalArea); // to nullify data.B/=2A in main pipeline

    // copy over entries to main B matrix
    for (SizeType row = 0; row < 2; row++)
    {
        for (SizeType col = 0; col < 18; col++)
        {
            data.B(row + 6, col) = smoothedShearMatrix(row, col);
        }
    }
}

void ShellThickElement3D3N::CalculateDSGShearBMatrix(Matrix& shearBMatrix,const double & a, const double & b, const double & c, const double & d, const double & A)
{
    //node 1
    shearBMatrix(0, 2) = b - c;
    shearBMatrix(0, 4) = A;

    shearBMatrix(1, 2) = d - a;
    shearBMatrix(1, 3) = -1.0*A;

    //node 2
    shearBMatrix(0, 8) = c;
    shearBMatrix(0, 9) = -1.0 * b*c / 2.0;
    shearBMatrix(0, 10) = a*c / 2.0;

    shearBMatrix(1, 8) = -1.0 * d;
    shearBMatrix(1, 9) = b*d / 2.0;
    shearBMatrix(1, 10) = -1.0 * a*d / 2.0;

    //node 3
    shearBMatrix(0, 14) = -1.0 * b;
    shearBMatrix(0, 15) = b*c / 2.0;
    shearBMatrix(0, 16) = -b*d / 2.0;

    shearBMatrix(1, 14) = a;
    shearBMatrix(1, 15) = -1.0 * a*c / 2.0;
    shearBMatrix(1, 16) = a*d / 2.0;

    shearBMatrix /= (2.0*A);
}

void ShellThickElement3D3N::AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector)
{
    // This is hardcoded to use 1 gauss point, despite the declared def
    // using 3 gps.
    const GeometryType& geom = GetGeometry();

    // Get shape functions
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
    const Matrix & N = GetGeometry().ShapeFunctionsValues(mIntegrationMethod);
#else
    // Disabled to use 1 gp below
    /*
    Matrix N(3, 3);
    for (unsigned int igauss = 0; igauss < mNumGPs; igauss++)
    {
        const array_1d<double, 3>& loc = data.gpLocations[igauss];
        N(igauss, 0) = 1.0 - loc[1] - loc[2];
        N(igauss, 1) = loc[1];
        N(igauss, 2) = loc[2];
    }
    */
    // 1 gp used!
    Matrix N(1, 3);
    N(0, 0) = 1.0 / 3.0;
    N(0, 1) = 1.0 / 3.0;
    N(0, 2) = 1.0 / 3.0;

#endif // !OPT_USES_INTERIOR_GAUSS_POINTS

    // auxiliary
    array_1d<double, 3> bf;

    // gauss loop to integrate the external force vector
    //for (unsigned int igauss = 0; igauss < mNumGPs; igauss++)
    for (unsigned int igauss = 0; igauss < 1; igauss++)
    {
        // get mass per unit area
        double mass_per_unit_area = mSections[igauss]->CalculateMassPerUnitArea(GetProperties());

        // interpolate nodal volume accelerations to this gauss point
        // and obtain the body force vector
        bf.clear();
        if (GetProperties().Has( VOLUME_ACCELERATION ))
            noalias(bf) = GetProperties()[VOLUME_ACCELERATION];
        else if (this->Has( VOLUME_ACCELERATION ))
            noalias(bf) = this->GetValue(VOLUME_ACCELERATION);
        else {
            for (unsigned int inode = 0; inode < 3; inode++)
            {
                if (geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION))
                {
                    noalias(bf) += N(igauss, inode) * geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
                }
            }
        }
        bf *= (mass_per_unit_area * data.dA);

        // add it to the RHS vector
        for (unsigned int inode = 0; inode < 3; inode++)
        {
            unsigned int index = inode * 6;
            double iN = N(igauss, inode);
            rRightHandSideVector[index + 0] += iN * bf[0];
            rRightHandSideVector[index + 1] += iN * bf[1];
            rRightHandSideVector[index + 2] += iN * bf[2];
        }
    }
}

void ShellThickElement3D3N::CalculateAll(MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
        // Resize the Left Hand Side if necessary,
        // and initialize it to Zero

    if ((rLeftHandSideMatrix.size1() != 18) || (rLeftHandSideMatrix.size2() != 18))
        rLeftHandSideMatrix.resize(18, 18, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(18, 18);

    // Resize the Right Hand Side if necessary,
    // and initialize it to Zero

    if (rRightHandSideVector.size() != 18)
        rRightHandSideVector.resize(18, false);
    noalias(rRightHandSideVector) = ZeroVector(18);

    // Initialize common calculation variables
    CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
    data.CalculateLHS = CalculateStiffnessMatrixFlag;
    data.CalculateRHS = CalculateResidualVectorFlag;
    InitializeCalculationData(data);
    CalculateSectionResponse(data);

    // Calculate element stiffness
    Matrix BTD = Matrix(18, 8, 0.0);
    data.D *= data.TotalArea;
    BTD = prod(trans(data.B), data.D);
    noalias(rLeftHandSideMatrix) += prod(BTD, data.B);
    if (data.specialDSGc3)
    {
        CalculateDSGc3Contribution(data, rLeftHandSideMatrix);
    }

    //add in z_rot artificial stiffness
    double z_rot_multiplier = 0.001;
    double max_stiff = 0.0; //max diagonal stiffness
    for (int i = 0; i < 18; i++)
    {
        if (rLeftHandSideMatrix(i, i) > max_stiff)
        {
            max_stiff = rLeftHandSideMatrix(i, i);
        }
    }
    for (int i = 0; i < 3; i++)
    {
        rLeftHandSideMatrix(6 * i + 5, 6 * i + 5) =
            z_rot_multiplier*max_stiff;
    }

    // Add RHS term
    rRightHandSideVector -= prod(rLeftHandSideMatrix, data.localDisplacements);

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
    KRATOS_CATCH("")
}

bool ShellThickElement3D3N::TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Check the required output
    int ijob = 0;
    bool bGlobal = false;
    CheckGeneralizedStressOrStrainOutput(rVariable, ijob, bGlobal);

    // quick return
    if (ijob == 0) return false;

    const SizeType num_gps = GetNumberOfGPs();

    // resize output
    if (rValues.size() != num_gps)
        rValues.resize(num_gps);

    // Just to store the rotation matrix for visualization purposes
    Matrix R(8, 8);
    Matrix aux33(3, 3);

    // Initialize common calculation variables
    CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
    if (ijob > 2)
    {
        data.CalculateLHS = true; // calc constitutive mat for forces/moments
    }
    else
    {
        data.CalculateLHS = false;
    }
    data.CalculateRHS = true;
    InitializeCalculationData(data);

    // Get the current displacements in global coordinate system and
    // transform to reference local system
    ShellT3_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem());
    MatrixType Rdisp(18, 18);
    referenceCoordinateSystem.ComputeTotalRotationMatrix(Rdisp);
    data.localDisplacements = prod(Rdisp, data.globalDisplacements);

    // set the current integration point index
    data.gpIndex = 0;
    ShellCrossSection::Pointer& section = mSections[0];

    // compute strains
    noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

    //compute forces
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
            CalculateSectionResponse(data);

            if (ijob > 4)
            {
                // Compute stresses for isotropic materials
                CalculateStressesFromForceResultants(data.generalizedStresses,
                    section->GetThickness(GetProperties()));
            }
        }
        DecimalCorrection(data.generalizedStresses);
    }

    // adjust output
    DecimalCorrection(data.generalizedStrains);

    // store the results, but first rotate them back to the section
    // coordinate system. we want to visualize the results in that system not
    // in the element one!
    if (section->GetOrientationAngle() != 0.0 && !bGlobal)
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
            section->GetRotationMatrixForGeneralizedStresses(-(section->GetOrientationAngle()), R);
            data.generalizedStresses = prod(R, data.generalizedStresses);
        }
        else
        {
            section->GetRotationMatrixForGeneralizedStrains(-(section->GetOrientationAngle()), R);
            data.generalizedStrains = prod(R, data.generalizedStrains);
        }
    }

    // Gauss Loop.

    for (SizeType i = 0; i < num_gps; i++)
    {
        // save results
        Matrix & iValue = rValues[i];
        if (iValue.size1() != 3 || iValue.size2() != 3)
            iValue.resize(3, 3, false);

        if (ijob == 1) // strains
        {
            iValue(0, 0) = data.generalizedStrains(0);
            iValue(1, 1) = data.generalizedStrains(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(2);
            iValue(0, 2) = iValue(2, 0) = 0.5 * data.generalizedStrains(7);
            iValue(1, 2) = iValue(2, 1) = 0.5 * data.generalizedStrains(6);
        }
        else if (ijob == 2) // curvatures
        {
            iValue(0, 0) = data.generalizedStrains(3);
            iValue(1, 1) = data.generalizedStrains(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(5);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        }
        else if (ijob == 3) // forces
        {
            iValue(0, 0) = data.generalizedStresses(0);
            iValue(1, 1) = data.generalizedStresses(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(2);
            iValue(0, 2) = iValue(2, 0) = data.generalizedStresses(7);
            iValue(1, 2) = iValue(2, 1) = data.generalizedStresses(6);
        }
        else if (ijob == 4) // moments
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
            iValue(0, 2) = iValue(2, 0) = data.generalizedStresses[6];
            iValue(1, 2) = iValue(2, 1) = data.generalizedStresses[7];
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
            iValue(0, 2) = iValue(2, 0) =
                data.rlaminateStresses[data.rlaminateStresses.size() - 1][6];
            iValue(1, 2) = iValue(2, 1) =
                data.rlaminateStresses[data.rlaminateStresses.size() - 1][7];
        }
        else if (ijob == 9) // SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE
        {
            iValue(0, 0) = data.rlaminateStresses[0][0];
            iValue(1, 1) = data.rlaminateStresses[0][1];
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.rlaminateStresses[0][2];
            iValue(0, 2) = iValue(2, 0) = data.rlaminateStresses[0][6];
            iValue(1, 2) = iValue(2, 1) = data.rlaminateStresses[0][7];
        }
        else if (ijob == 99) // SHELL_ORTHOTROPIC_4PLY_THROUGH_THICKNESS
        {
            // Testing variable to get lamina stress/strain values
            // on each surface of a 4 ply laminate

            int surface = 0; // start from top ply top surface
                                // Output global results sequentially
            for (SizeType row = 0; row < 3; row++)
            {
                for (SizeType col = 0; col < 3; col++)
                {
                    if (surface > 7)
                    {
                        iValue(row, col) = 0.0;
                    }
                    else
                    {
                        iValue(row, col) = data.rlaminateStrains[surface][6];
                    }
                    surface++;
                }
            }


            bool tsai_wu_thru_output = false;
            if (tsai_wu_thru_output)
            {
                std::vector<Matrix> Laminae_Strengths =
                    std::vector<Matrix>(section->NumberOfPlies());
                for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
                {
                    Laminae_Strengths[ply].resize(3, 3, 0.0);
                    Laminae_Strengths[ply].clear();
                }
                const PropertiesType & props = GetProperties();
                section->GetLaminaeStrengths(Laminae_Strengths, props);
                Vector ply_orientation(section->NumberOfPlies());
                section->GetLaminaeOrientation(props, ply_orientation);

                CalculateLaminaStrains(data);
                CalculateLaminaStresses(data);

                // Rotate lamina stress from section CS
                // to lamina angle to lamina material principal directions
                for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
                {
                    double total_rotation = -ply_orientation[ply] - (section->GetOrientationAngle()); // already rotated to section CS
                    section->GetRotationMatrixForGeneralizedStresses(total_rotation, R);
                    //top surface of current ply
                    data.rlaminateStresses[2 * ply] = prod(R, data.rlaminateStresses[2 * ply]);
                    //bottom surface of current ply
                    data.rlaminateStresses[2 * ply + 1] = prod(R, data.rlaminateStresses[2 * ply + 1]);
                }

                // Calculate Tsai-Wu criterion for each ply
                Vector tsai_output = Vector(9, 0.0);
                tsai_output.clear();
                double temp_tsai_wu = 0.0;
                for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++)
                {
                    Vector lamina_stress_top = Vector(data.rlaminateStresses[2 * ply]);
                    Vector lamina_stress_bottom = Vector(data.rlaminateStresses[2 * ply + 1]);

                    // top surface
                    data.rlaminateStresses[2 * ply + 1] = Vector(lamina_stress_top);
                    temp_tsai_wu = CalculateTsaiWuPlaneStress(data, Laminae_Strengths[ply], ply);
                    tsai_output[2 * ply] = temp_tsai_wu;

                    // bottom surface
                    data.rlaminateStresses[2 * ply] = Vector(lamina_stress_bottom);
                    data.rlaminateStresses[2 * ply + 1] = Vector(lamina_stress_bottom);
                    temp_tsai_wu = CalculateTsaiWuPlaneStress(data, Laminae_Strengths[ply], ply);
                    tsai_output[2 * ply + 1] = temp_tsai_wu;
                }

                //dump into results

                int surface1 = 0; // start from top ply top surface
                                    // Output global results sequentially
                for (SizeType row = 0; row < 3; row++)
                {
                    for (SizeType col = 0; col < 3; col++)
                    {
                        if (surface1 > 7)
                        {
                            iValue(row, col) = 0.0;
                        }
                        else
                        {
                            iValue(row, col) = tsai_output[surface1];
                        }
                        surface1++;
                    }
                }

            } // tsai wu output
        }

        // if requested, rotate the results in the global coordinate system
        if (bGlobal)
        {
            const Matrix& RG = data.LCS.Orientation();
            noalias(aux33) = prod(trans(RG), iValue);
            noalias(iValue) = prod(aux33, RG);
        }
    }

    OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
    return true;
}

ShellCrossSection::SectionBehaviorType ShellThickElement3D3N::GetSectionBehavior()
{
    return ShellCrossSection::Thick;
}

// =====================================================================================
//
// Class ShellThickElement3D3N - Serialization
//
// =====================================================================================

void ShellThickElement3D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseShellElement);
    rSerializer.save("CTr", mpCoordinateTransformation);
}

void ShellThickElement3D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseShellElement);
    rSerializer.load("CTr", mpCoordinateTransformation);
}

} // namespace Kratos
