// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//       Contact:    A.Winterstein@tum.de
//

#include "shell_thin_element_3D4N.hpp"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"

#include "custom_utilities/shell_utilities.h"

#include <string>
#include <iomanip>

/*
Element overview:---------------------------------------------------------------
This element represents a 4-node Shell element.
The membrane part is Felippa's assumed Natural DEviatoric Strain (ANDES)
formulation, while the bending part is the Discrete Kirchhoff Quadrilateral.
This element is formulated for small strains,
but can be used in Geometrically nonlinear problems
involving large displacements and rotations
using a Corotational Coordinate Transformation.
Material nonlinearity is handled by means of the cross section object.

Shell formulation references:---------------------------------------------------
ANDES formulation:
Bjorn Haugen. "Buckling and Stability Problems for Thin Shell Structures
Using High Performance Finite Elements". Dissertation. Colorado: University
of Colorado, 1994.

ANDES filter matrix H modification as per:
Carlos A Felippa. "Supernatural QUAD4: a template formulation".    In: Computer
methods in applied mechanics and engineering 195.41 (2006), pp. 5316-5342.

DKQ original formulation:
Jean-Louis Batoz and Mabrouk Ben Tahar. "Evaluation of a new quadrilateral
thin plate bending element". In: International Journal for Numerical Methods
in Engineering 18.11 (1982), pp. 1655-1677.

Clearly presented DKQ formulation:
Fabian Rojas Barrales. "Development of a nonlinear quadrilateral layered
membrane element with drilling degrees of freedom and a nonlinear
quadrilateral thin flat layered shell element for the modeling of reinforced
concrete walls". Dissertation. Los Angeles, California: University of
Southern California, 2012.
*/

namespace Kratos
{
// =========================================================================
//
// Class JacobianOperator
//
// =========================================================================

ShellThinElement3D4N::JacobianOperator::JacobianOperator()
    : mJac(2, 2, 0.0)
    , mInv(2, 2, 0.0)
    , mXYDeriv(4, 2, 0.0)
    , mDet(0.0)
{
}

void ShellThinElement3D4N::JacobianOperator::Calculate
(const ShellQ4_LocalCoordinateSystem& CS, const Matrix& dN)
{
    mJac(0, 0) = dN(0, 0) * CS.X1() + dN(1, 0) * CS.X2() +
                 dN(2, 0) * CS.X3() + dN(3, 0) * CS.X4();
    mJac(0, 1) = dN(0, 0) * CS.Y1() + dN(1, 0) * CS.Y2() +
                 dN(2, 0) * CS.Y3() + dN(3, 0) * CS.Y4();
    mJac(1, 0) = dN(0, 1) * CS.X1() + dN(1, 1) * CS.X2() +
                 dN(2, 1) * CS.X3() + dN(3, 1) * CS.X4();
    mJac(1, 1) = dN(0, 1) * CS.Y1() + dN(1, 1) * CS.Y2() +
                 dN(2, 1) * CS.Y3() + dN(3, 1) * CS.Y4();

    mDet = mJac(0, 0) * mJac(1, 1) - mJac(1, 0) * mJac(0, 1);
    double mult = 1.0 / mDet;

    mInv(0, 0) = mJac(1, 1) * mult;
    mInv(0, 1) = -mJac(0, 1) * mult;
    mInv(1, 0) = -mJac(1, 0) * mult;
    mInv(1, 1) = mJac(0, 0) * mult;

    noalias(mXYDeriv) = prod(dN, trans(mInv));
}

// =========================================================================
//
// Class ShellThinElement3D4N
//
// =========================================================================

ShellThinElement3D4N::ShellThinElement3D4N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        bool NLGeom)
    : BaseShellElement(NewId, pGeometry)
    , mpCoordinateTransformation(NLGeom ?
                                 new ShellQ4_CorotationalCoordinateTransformation(pGeometry) :
                                 new ShellQ4_CoordinateTransformation(pGeometry))
{
}

ShellThinElement3D4N::ShellThinElement3D4N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        bool NLGeom)
    : BaseShellElement(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation(NLGeom ?
                                 new ShellQ4_CorotationalCoordinateTransformation(pGeometry) :
                                 new ShellQ4_CoordinateTransformation(pGeometry))
{
}

ShellThinElement3D4N::ShellThinElement3D4N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        CoordinateTransformationBasePointerType pCoordinateTransformation)
    : BaseShellElement(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation(pCoordinateTransformation)
{
}

ShellThinElement3D4N::~ShellThinElement3D4N()
{
}

//Basic methods

Element::Pointer ShellThinElement3D4N::Create(IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
    GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
    return Kratos::make_intrusive< ShellThinElement3D4N >(NewId, newGeom,
            pProperties, mpCoordinateTransformation->Create(newGeom));
}

Element::Pointer ShellThinElement3D4N::Create(IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< ShellThinElement3D4N >(NewId, pGeom,
            pProperties, mpCoordinateTransformation->Create(pGeom));
}

void ShellThinElement3D4N::Initialize()
{
    KRATOS_TRY

    const int points_number = GetGeometry().PointsNumber();

    KRATOS_ERROR_IF_NOT(points_number == 4) <<"ShellThinElement3D4N - Wrong number of nodes"
                                            << points_number << std::endl;

    BaseShellElement::Initialize();

    mpCoordinateTransformation->Initialize();

    this->SetupOrientationAngles();

    KRATOS_CATCH("")
}

void ShellThinElement3D4N::InitializeNonLinearIteration
(ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->InitializeNonLinearIteration();

    BaseInitializeNonLinearIteration(rCurrentProcessInfo);
}

void ShellThinElement3D4N::FinalizeNonLinearIteration
(ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->FinalizeNonLinearIteration();

    BaseFinalizeNonLinearIteration(rCurrentProcessInfo);
}

void ShellThinElement3D4N::InitializeSolutionStep
(ProcessInfo& rCurrentProcessInfo)
{
    BaseInitializeSolutionStep(rCurrentProcessInfo);

    mpCoordinateTransformation->InitializeSolutionStep();
}

void ShellThinElement3D4N::FinalizeSolutionStep
(ProcessInfo& rCurrentProcessInfo)
{
    BaseFinalizeSolutionStep(rCurrentProcessInfo);

    mpCoordinateTransformation->FinalizeSolutionStep();
}

void ShellThinElement3D4N::CalculateMassMatrix(MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo)
{
    if ((rMassMatrix.size1() != 24) || (rMassMatrix.size2() != 24)) {
        rMassMatrix.resize(24, 24, false);
    }
    noalias(rMassMatrix) = ZeroMatrix(24, 24);

    // Compute the local coordinate system.
    ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem());

    // Average mass per unit area over the whole element
    double av_mass_per_unit_area = 0.0;

    // Flag for consistent or lumped mass matrix
    bool bconsistent_matrix = false;

    // Consistent mass matrix
    if (bconsistent_matrix) {
        // Get shape function values and setup jacobian
        const GeometryType& geom = GetGeometry();
        const Matrix& shapeFunctions = geom.ShapeFunctionsValues();
        JacobianOperator jacOp;

        // Get integration points
        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints(mIntegrationMethod);

        // Setup matrix of shape functions
        Matrix N = Matrix(6, 24, 0.0);

        // Other variables
        double dA = 0.0;
        double thickness = 0.0;
        double drilling_factor = 1.0;    // sqrt of the actual factor applied,
        // 1.0 is no reduction.

        // Gauss loop
        for (SizeType gauss_point = 0; gauss_point < 4; gauss_point++) {
            // Calculate average mass per unit area and thickness at the
            // current GP
            av_mass_per_unit_area =
                mSections[gauss_point]->CalculateMassPerUnitArea(GetProperties());
            thickness = mSections[gauss_point]->GetThickness(GetProperties());

            // Calc jacobian and weighted dA at current GP
            jacOp.Calculate(referenceCoordinateSystem,
                            geom.ShapeFunctionLocalGradient(gauss_point));
            dA = integration_points[gauss_point].Weight() *
                 jacOp.Determinant();

            // Assemble shape function matrix over nodes
            for (SizeType node = 0; node < 4; node++) {
                // translational entries - dofs 1, 2, 3
                for (SizeType dof = 0; dof < 3; dof++) {
                    N(dof, 6 * node + dof) =
                        shapeFunctions(gauss_point, node);
                }

                // rotational inertia entries - dofs 4, 5
                for (SizeType dof = 0; dof < 2; dof++) {
                    N(dof + 3, 6 * node + dof + 3) =
                        thickness / std::sqrt(12.0) *
                        shapeFunctions(gauss_point, node);
                }

                // drilling rotational entry - artifical factor included
                N(5, 6 * node + 5) = thickness / std::sqrt(12.0) *
                                     shapeFunctions(gauss_point, node) /
                                     drilling_factor;
            }

            // Add contribution to total mass matrix
            rMassMatrix += prod(trans(N), N)*dA*av_mass_per_unit_area;
        }

    }// Consistent mass matrix
    else {
        // Lumped mass matrix

        // Calculate average mass per unit area over the whole element
        for (SizeType i = 0; i < 4; i++) {
            av_mass_per_unit_area += mSections[i]->CalculateMassPerUnitArea(GetProperties());
        }
        av_mass_per_unit_area /= 4.0;

        // lumped area
        double lump_area = referenceCoordinateSystem.Area() / 4.0;

        // Gauss Loop
        for (SizeType i = 0; i < 4; i++) {
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

void ShellThinElement3D4N::AddBodyForces(CalculationData& data,
        VectorType& rRightHandSideVector)
{
    const GeometryType& geom = GetGeometry();

    // Get shape functions
    const Matrix& N = geom.ShapeFunctionsValues();

    // auxiliary
    array_1d<double, 3> bf;

    // gauss loop to integrate the external force vector
    for (unsigned int igauss = 0; igauss < 4; igauss++) {
        // get mass per unit area
        double mass_per_unit_area =
            mSections[igauss]->CalculateMassPerUnitArea(GetProperties());

        // interpolate nodal volume accelerations to this gauss point
        // and obtain the body force vector
        bf.clear();
        for (unsigned int inode = 0; inode < 4; inode++) {
            if (geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION))
                //temporary, will be checked once at the beginning only
                bf += N(igauss, inode) *
                      geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
        }
        bf *= (mass_per_unit_area * data.dA[igauss]);

        // add it to the RHS vector
        for (unsigned int inode = 0; inode < 4; inode++) {
            unsigned int index = inode * 6;
            double iN = N(igauss, inode);
            rRightHandSideVector[index + 0] += iN * bf[0];
            rRightHandSideVector[index + 1] += iN * bf[1];
            rRightHandSideVector[index + 2] += iN * bf[2];
        }
    }
}

// =========================================================================
//
// Class ShellThinElement3D4N - Results on Gauss Points
//
// =========================================================================

void ShellThinElement3D4N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int caseId = -1;
    if (rVariable == TSAI_WU_RESERVE_FACTOR) {
        caseId = 10;
    } else if (rVariable == SHEAR_ANGLE) {
        caseId = 40;
    } else if (rVariable == VON_MISES_STRESS ||
               rVariable == VON_MISES_STRESS_TOP_SURFACE ||
               rVariable == VON_MISES_STRESS_MIDDLE_SURFACE ||
               rVariable == VON_MISES_STRESS_BOTTOM_SURFACE) {
        caseId = 20;
    } else if (rVariable == SHELL_ELEMENT_MEMBRANE_ENERGY ||
               SHELL_ELEMENT_MEMBRANE_ENERGY_FRACTION ||
               SHELL_ELEMENT_BENDING_ENERGY ||
               SHELL_ELEMENT_BENDING_ENERGY_FRACTION ||
               SHELL_ELEMENT_SHEAR_ENERGY ||
               SHELL_ELEMENT_SHEAR_ENERGY_FRACTION) {
        caseId = 30;
    }

    if (caseId > 19) {
        // resize output
        SizeType size = 4;
        if (rValues.size() != size) {
            rValues.resize(size);
        }

        // Compute the local coordinate system.
        ShellQ4_LocalCoordinateSystem localCoordinateSystem(
            mpCoordinateTransformation->CreateLocalCoordinateSystem());
        ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
            mpCoordinateTransformation->CreateReferenceCoordinateSystem());

        // Initialize common calculation variables
        CalculationData data(localCoordinateSystem,
                             referenceCoordinateSystem, rCurrentProcessInfo);
        data.CalculateLHS = false;
        data.CalculateRHS = true;
        InitializeCalculationData(data);

        // Get the current displacements in global coordinate system and
        // transform to reference local system
        MatrixType Rdisp(24, 24);
        referenceCoordinateSystem.ComputeTotalRotationMatrix(Rdisp);
        if (referenceCoordinateSystem.IsWarped()) {
            MatrixType W(24, 24);
            referenceCoordinateSystem.ComputeTotalWarpageMatrix(W);
            Rdisp = prod(W, Rdisp);
        }
        data.localDisplacements = prod(Rdisp, data.globalDisplacements);


        // loop over gauss points
        for (unsigned int gauss_point = 0; gauss_point < size; ++gauss_point) {
            // Compute all strain-displacement matrices
            data.gpIndex = gauss_point;
            CalculateBMatrix(data);

            // Calculate strain vectors in local coordinate system
            noalias(data.generalizedStrains) =
                prod(data.B, data.localDisplacements);

            // Calculate the response of the Cross Section
            ShellCrossSection::Pointer& section = mSections[gauss_point];
            CalculateSectionResponse(data);

            double resultDouble = 0.0;

            if (caseId == 30) {
                // Energy calcs - these haven't been verified or tested yet.
                CalculateShellElementEnergy(data, rVariable, resultDouble);
            } else if (caseId == 20) {
                //Von mises calcs

                // recover stresses
                CalculateStressesFromForceResultants(data.generalizedStresses,
                                                     section->GetThickness(GetProperties()));

                // account for orientation
                if (section->GetOrientationAngle() != 0.0) {
                    Matrix R(6, 6);
                    section->GetRotationMatrixForGeneralizedStresses
                    (-(section->GetOrientationAngle()), R);
                    data.generalizedStresses = prod(R, data.generalizedStresses);
                }

                CalculateVonMisesStress(data, rVariable, resultDouble);
            } else if (caseId == 40) {
                resultDouble = std::atan(0.5*data.generalizedStrains[2]) * 180 / Kratos::Globals::Pi;
            } else {
                KRATOS_ERROR <<
                             "Error: ELEMENT ShellThinElement3D4N, METHOD CalculateOnIntegrationPoints(double)"
                             << std::endl;
            }

            // store the result calculated
            rValues[gauss_point] = resultDouble;
        }
    } else if (rVariable == TSAI_WU_RESERVE_FACTOR) {
        // resize output
        SizeType size = 4;
        if (rValues.size() != size) {
            rValues.resize(size);
        }

        //CalculationData data = SetupStressOrStrainCalculation(rCurrentProcessInfo);
        // Compute the local coordinate system.
        ShellQ4_LocalCoordinateSystem localCoordinateSystem(
            mpCoordinateTransformation->CreateLocalCoordinateSystem());
        ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
            mpCoordinateTransformation->CreateReferenceCoordinateSystem());

        // Initialize common calculation variables
        CalculationData data(localCoordinateSystem,
                             referenceCoordinateSystem, rCurrentProcessInfo);
        data.CalculateLHS = true;
        data.CalculateRHS = true;
        InitializeCalculationData(data);

        // Get the current displacements in global coordinate system and
        // transform to reference local system
        MatrixType Rdisp(24, 24);
        referenceCoordinateSystem.ComputeTotalRotationMatrix(Rdisp);
        if (referenceCoordinateSystem.IsWarped()) {
            MatrixType W(24, 24);
            referenceCoordinateSystem.ComputeTotalWarpageMatrix(W);
            Rdisp = prod(W, Rdisp);
        }
        data.localDisplacements = prod(Rdisp, data.globalDisplacements);


        // Get all laminae strengths
        const PropertiesType& props = GetProperties();
        ShellCrossSection::Pointer& section = mSections[0];
        std::vector<Matrix> Laminae_Strengths =
            std::vector<Matrix>(section->NumberOfPlies());
        for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++) {
            Laminae_Strengths[ply].resize(3, 3, 0.0);
            Laminae_Strengths[ply].clear();
        }
        section->GetLaminaeStrengths(Laminae_Strengths,props);

        // Define variables
        Matrix R(6, 6);
        double total_rotation = 0.0;

        // Gauss Loop
        for (unsigned int gauss_point = 0; gauss_point < size; gauss_point++) {
            // Compute all strain-displacement matrices
            data.gpIndex = gauss_point;
            CalculateBMatrix(data);

            // Calculate strain vectors in local coordinate system
            noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

            // Retrieve ply orientations
            section = mSections[gauss_point];
            Vector ply_orientation(section->NumberOfPlies());
            section->GetLaminaeOrientation(props, ply_orientation);

            //Calculate lamina stresses
            CalculateLaminaStrains(data);
            CalculateLaminaStresses(data);

            // Rotate lamina stress from element CS to section CS, and then
            // to lamina angle to lamina material principal directions
            for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++) {
                total_rotation = -ply_orientation[ply] - (section->GetOrientationAngle());
                section->GetRotationMatrixForGeneralizedStresses(total_rotation, R);
                //top surface of current ply
                data.rlaminateStresses[2*ply] = prod(R, data.rlaminateStresses[2*ply]);
                //bottom surface of current ply
                data.rlaminateStresses[2 * ply +1] = prod(R, data.rlaminateStresses[2 * ply +1]);
            }

            // Calculate Tsai-Wu criterion for each ply, take min of all plies
            double min_tsai_wu = 0.0;
            double temp_tsai_wu = 0.0;
            for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++) {
                temp_tsai_wu = CalculateTsaiWuPlaneStress(data, Laminae_Strengths[ply],ply);
                if (ply == 0) {
                    min_tsai_wu = temp_tsai_wu;
                } else if (temp_tsai_wu < min_tsai_wu) {
                    min_tsai_wu = temp_tsai_wu;
                }
            }

            // Output min Tsai-Wu result
            rValues[gauss_point] = min_tsai_wu;

        }// Gauss loop
    } else {
        SizeType size = GetGeometry().size();
        if (rValues.size() != size) {
            rValues.resize(size);
        }

        std::vector<double> temp(size);

        for (SizeType i = 0; i < size; i++) {
            mSections[i]->GetValue(rVariable, GetProperties(), temp[i]);
        }

        const Matrix& shapeFunctions = GetGeometry().ShapeFunctionsValues();
        Vector N(size);

        for (SizeType i = 0; i < size; i++) {
            noalias(N) = row(shapeFunctions, i);
            double& ival = rValues[i];
            ival = 0.0;
            for (SizeType j = 0; j < size; j++) {
                ival += N(j) * temp[j];
            }
        }
    }

    KRATOS_CATCH("");
}

void ShellThinElement3D4N::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(
                rVariable, rValues, rCurrentProcessInfo)) {
        return;
    }
}

void ShellThinElement3D4N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3> >& rVariable,
    std::vector<array_1d<double, 3> >& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == LOCAL_AXIS_1 ||
            rVariable == LOCAL_AXIS_2 ||
            rVariable == LOCAL_AXIS_3) {
        BaseShellElement::ComputeLocalAxis(rVariable, rOutput, mpCoordinateTransformation);
    } else if (rVariable == LOCAL_MATERIAL_AXIS_1 ||
               rVariable == LOCAL_MATERIAL_AXIS_2 ||
               rVariable == LOCAL_MATERIAL_AXIS_3) {
        BaseShellElement::ComputeLocalMaterialAxis(rVariable, rOutput, mpCoordinateTransformation);
    }
}

void ShellThinElement3D4N::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == LOCAL_ELEMENT_ORIENTATION) {
        Output.resize(3, 3, false);

        // Compute the local coordinate system.
        ShellQ4_LocalCoordinateSystem localCoordinateSystem(
            mpCoordinateTransformation->CreateReferenceCoordinateSystem());
        Output = trans(localCoordinateSystem.Orientation());
    }
}

// =========================================================================
//
// Class ShellThinElement3D4N - Private methods
//
// =========================================================================

void ShellThinElement3D4N::CalculateStressesFromForceResultants
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
}

void ShellThinElement3D4N::CalculateLaminaStrains(CalculationData& data)
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

    // Get ply thicknesses
    Vector ply_thicknesses = Vector(section->NumberOfPlies(), 0.0);
    section->GetPlyThicknesses(GetProperties(), ply_thicknesses);

    // Resize output vector. 2 Surfaces for each ply
    data.rlaminateStrains.resize(2 * section->NumberOfPlies());
    for (unsigned int i = 0; i < 2 * section->NumberOfPlies(); i++) {
        data.rlaminateStrains[i].resize(6, false);
        data.rlaminateStrains[i].clear();
    }

    // Loop over all plies - start from bottom ply, bottom surface
    for (unsigned int plyNumber = 0;
            plyNumber < section->NumberOfPlies(); ++plyNumber) {
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

void ShellThinElement3D4N::CalculateLaminaStresses(CalculationData& data)
{
    ShellCrossSection::Pointer& section = mSections[data.gpIndex];

    // Setup flag to compute ply constitutive matrices
    // (units [Pa] and rotated to element orientation)
    section->SetupGetPlyConstitutiveMatrices();
    CalculateSectionResponse(data);

    // Resize output vector. 2 Surfaces for each ply
    data.rlaminateStresses.resize(2 * section->NumberOfPlies());
    for (unsigned int i = 0; i < 2 * section->NumberOfPlies(); i++) {
        data.rlaminateStresses[i].resize(6, false);
        data.rlaminateStresses[i].clear();
    }

    // Loop over all plies - start from top ply, top surface
    for (unsigned int plyNumber = 0;
            plyNumber < section->NumberOfPlies(); ++plyNumber) {
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

double ShellThinElement3D4N::CalculateTsaiWuPlaneStress(const CalculationData& data, const Matrix& rLamina_Strengths, const unsigned int& rPly)
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
    Matrix F_ij = Matrix(3, 3, 0.0);
    F_ij.clear();
    F_ij(0, 0) = 1.0 / rLamina_Strengths(0, 0) / rLamina_Strengths(0, 1);    // 11
    F_ij(1,1) = 1.0 / rLamina_Strengths(0, 2) / rLamina_Strengths(1, 0);    // 22
    F_ij(2, 2) = 1.0 / rLamina_Strengths(1, 1) / rLamina_Strengths(1, 1);    // 12
    F_ij(0, 1) = F_ij(1, 0) = -0.5 / std::sqrt(rLamina_Strengths(0, 0)*rLamina_Strengths(0, 1)*rLamina_Strengths(0, 2)*rLamina_Strengths(1, 0));

    if (disable_in_plane_interaction) {
        F_ij(0, 1) = F_ij(1, 0) = 0.0;
    }


    // Evaluate Tsai-Wu @ top surface of current layer
    double var_a = 0.0;
    double var_b = 0.0;
    for (SizeType i = 0; i < 3; i++) {
        var_b += F_i[i] * data.rlaminateStresses[2 * rPly][i];
        for (SizeType j = 0; j < 3; j++) {
            var_a += F_ij(i, j)*data.rlaminateStresses[2 * rPly][i] * data.rlaminateStresses[2 * rPly][j];
        }
    }
    double tsai_reserve_factor_top = (-1.0*var_b + std::sqrt(var_b*var_b + 4.0 * var_a)) / 2.0 / var_a;

    // Evaluate Tsai-Wu @ bottom surface of current layer
    var_a = 0.0;
    var_b = 0.0;
    for (SizeType i = 0; i < 3; i++) {
        var_b += F_i[i] * data.rlaminateStresses[2 * rPly + 1][i];
        for (SizeType j = 0; j < 3; j++) {
            var_a += F_ij(i, j)*data.rlaminateStresses[2 * rPly + 1][i] * data.rlaminateStresses[2 * rPly + 1][j];
        }
    }
    double tsai_reserve_factor_bottom = (-1.0*var_b + std::sqrt(var_b*var_b + 4.0 * var_a)) / 2.0 / var_a;

    // Return min of both surfaces as the result for the whole ply
    return std::min(tsai_reserve_factor_bottom, tsai_reserve_factor_top);
}

void ShellThinElement3D4N::CalculateVonMisesStress(const CalculationData& data, const Variable<double>& rVariable, double& rVon_Mises_Result)
{
    // calc von mises stresses at top mid and bottom surfaces for
    // thin shell
    double von_mises_top, von_mises_mid, von_mises_bottom;
    double sxx, syy, sxy;

    // top surface: membrane and +bending contributions
    //                (no transverse shear)
    sxx = data.generalizedStresses[0] + data.generalizedStresses[3];
    syy = data.generalizedStresses[1] + data.generalizedStresses[4];
    sxy = data.generalizedStresses[2] + data.generalizedStresses[5];
    von_mises_top = sxx*sxx - sxx*syy + syy*syy + 3.0*sxy*sxy;

    // mid surface: membrane only contributions
    //                (no bending or transverse shear)
    sxx = data.generalizedStresses[0];
    syy = data.generalizedStresses[1];
    sxy = data.generalizedStresses[2];
    von_mises_mid = sxx*sxx - sxx*syy + syy*syy +
                    3.0*(sxy*sxy);

    // bottom surface:    membrane and bending contributions
    //                    (no transverse shear)
    sxx = data.generalizedStresses[0] - data.generalizedStresses[3];
    syy = data.generalizedStresses[1] - data.generalizedStresses[4];
    sxy = data.generalizedStresses[2] - data.generalizedStresses[5];
    von_mises_bottom = sxx*sxx - sxx*syy + syy*syy + 3.0*sxy*sxy;

    // Output requested quantity
    if (rVariable == VON_MISES_STRESS_TOP_SURFACE) {
        rVon_Mises_Result = std::sqrt(von_mises_top);
    } else if (rVariable == VON_MISES_STRESS_MIDDLE_SURFACE) {
        rVon_Mises_Result = std::sqrt(von_mises_mid);
    } else if (rVariable == VON_MISES_STRESS_BOTTOM_SURFACE) {
        rVon_Mises_Result = std::sqrt(von_mises_bottom);
    } else if (rVariable == VON_MISES_STRESS) {
        // take the greatest value and output
        rVon_Mises_Result =
            std::sqrt(std::max(von_mises_top,
                               std::max(von_mises_mid, von_mises_bottom)));
    }
}

void ShellThinElement3D4N::CalculateShellElementEnergy(const CalculationData& data, const Variable<double>& rVariable, double& rEnergy_Result)
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
            rVariable == SHELL_ELEMENT_SHEAR_ENERGY_FRACTION) {
        // need to calculate total energy over current dA first
        totalEnergy = inner_prod(data.generalizedStresses, data.generalizedStrains)*data.dA[data.gpIndex];
        is_fraction_calc = true;
    }

    if (rVariable == SHELL_ELEMENT_MEMBRANE_ENERGY || rVariable == SHELL_ELEMENT_MEMBRANE_ENERGY_FRACTION) {
        for (SizeType i = 0; i < 3; i++) {
            rEnergy_Result += data.generalizedStresses[i] * data.generalizedStrains[i]* data.dA[data.gpIndex];
        }

        if (is_fraction_calc) {
            rEnergy_Result /= totalEnergy;
        }
    } else if (rVariable == SHELL_ELEMENT_BENDING_ENERGY || rVariable == SHELL_ELEMENT_BENDING_ENERGY_FRACTION) {
        for (SizeType i = 3; i < 6; i++) {
            rEnergy_Result += data.generalizedStresses[i] * data.generalizedStrains[i]* data.dA[data.gpIndex];
        }

        if (is_fraction_calc) {
            rEnergy_Result /= totalEnergy;
        }
    } else if (rVariable == SHELL_ELEMENT_SHEAR_ENERGY || rVariable == SHELL_ELEMENT_SHEAR_ENERGY_FRACTION) {
        rEnergy_Result = 0.0;
    }
}

void ShellThinElement3D4N::CheckGeneralizedStressOrStrainOutput(const Variable<Matrix>& rVariable, int& ijob, bool& bGlobal)
{
    if (rVariable == SHELL_STRAIN) {
        ijob = 1;
    } else if (rVariable == SHELL_STRAIN_GLOBAL) {
        ijob = 1;
        bGlobal = true;
    } else if (rVariable == SHELL_CURVATURE) {
        ijob = 2;
    } else if (rVariable == SHELL_CURVATURE_GLOBAL) {
        ijob = 2;
        bGlobal = true;
    } else if (rVariable == SHELL_FORCE) {
        ijob = 3;
    } else if (rVariable == SHELL_FORCE_GLOBAL) {
        ijob = 3;
        bGlobal = true;
    } else if (rVariable == SHELL_MOMENT) {
        ijob = 4;
    } else if (rVariable == SHELL_MOMENT_GLOBAL) {
        ijob = 4;
        bGlobal = true;
    } else if (rVariable == SHELL_STRESS_TOP_SURFACE) {
        ijob = 5;
    } else if (rVariable == SHELL_STRESS_TOP_SURFACE_GLOBAL) {
        ijob = 5;
        bGlobal = true;
    } else if (rVariable == SHELL_STRESS_MIDDLE_SURFACE) {
        ijob = 6;
    } else if (rVariable == SHELL_STRESS_MIDDLE_SURFACE_GLOBAL) {
        ijob = 6;
        bGlobal = true;
    } else if (rVariable == SHELL_STRESS_BOTTOM_SURFACE) {
        ijob = 7;
    } else if (rVariable == SHELL_STRESS_BOTTOM_SURFACE_GLOBAL) {
        ijob = 7;
        bGlobal = true;
    } else if (rVariable == SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE) {
        ijob = 8;
    } else if (rVariable == SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE_GLOBAL) {
        ijob = 8;
        bGlobal = true;
    } else if (rVariable == SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE) {
        ijob = 9;
    } else if (rVariable == SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE_GLOBAL) {
        ijob = 9;
        bGlobal = true;
    } else if (rVariable == SHELL_ORTHOTROPIC_4PLY_THROUGH_THICKNESS) {
        // TESTING VARIABLE
        ijob = 99;
    }
}

void ShellThinElement3D4N::DecimalCorrection(Vector& a)
{
    double norm = norm_2(a);
    double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
    for (SizeType i = 0; i < a.size(); i++)
        if (std::abs(a(i)) < tolerance) {
            a(i) = 0.0;
        }
}

void ShellThinElement3D4N::SetupOrientationAngles()
{
    if (this->Has(MATERIAL_ORIENTATION_ANGLE)) {
        for (CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it) {
            (*it)->SetOrientationAngle(this->GetValue(MATERIAL_ORIENTATION_ANGLE));
        }
    } else {
        ShellQ4_LocalCoordinateSystem lcs(mpCoordinateTransformation->
                                          CreateReferenceCoordinateSystem());

        Vector3Type normal;
        noalias(normal) = lcs.Vz();

        Vector3Type dZ;
        dZ(0) = 0.0;
        dZ(1) = 0.0;
        dZ(2) = 1.0;

        Vector3Type dirX;
        MathUtils<double>::CrossProduct(dirX, dZ, normal);

        // try to normalize the x vector. if it is near zero it means that we
        // need to choose a default one.
        double dirX_norm = dirX(0)*dirX(0) + dirX(1)*dirX(1) + dirX(2)*dirX(2);
        if (dirX_norm < 1.0E-12) {
            dirX(0) = 1.0;
            dirX(1) = 0.0;
            dirX(2) = 0.0;
        } else if (dirX_norm != 1.0) {
            dirX_norm = std::sqrt(dirX_norm);
            dirX /= dirX_norm;
        }

        Vector3Type elem_dirX = lcs.Vx();

        // now calculate the angle between the element x direction and the
        // material x direction.
        Vector3Type& a = elem_dirX;
        Vector3Type& b = dirX;
        double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
        if (a_dot_b < -1.0) {
            a_dot_b = -1.0;
        }
        if (a_dot_b > 1.0) {
            a_dot_b = 1.0;
        }
        double angle = std::acos(a_dot_b);

        // if they are not counter-clock-wise,
        // let's change the sign of the angle
        if (angle != 0.0) {
            const MatrixType& R = lcs.Orientation();
            if (dirX(0)*R(1, 0) + dirX(1)*R(1, 1) + dirX(2)*R(1, 2) < 0.0) {
                angle = -angle;
            }
        }

        for (CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it) {
            (*it)->SetOrientationAngle(angle);
        }
    }
}

void ShellThinElement3D4N::InitializeCalculationData(CalculationData& data)
{
    KRATOS_TRY
    //-------------------------------------
    // Computation of all stuff that remain
    // constant throughout the calculations

    //-------------------------------------
    // geometry data
    const double x12 = data.LCS0.X1() - data.LCS0.X2();
    const double x13 = data.LCS0.X1() - data.LCS0.X3();
    const double x23 = data.LCS0.X2() - data.LCS0.X3();
    const double x24 = data.LCS0.X2() - data.LCS0.X4();
    const double x34 = data.LCS0.X3() - data.LCS0.X4();
    const double x41 = data.LCS0.X4() - data.LCS0.X1();

    const double x21 = -x12;
    const double x31 = -x13;
    const double x32 = -x23;
    const double x42 = -x24;
    const double x43 = -x34;
    const double x14 = -x41;

    const double y12 = data.LCS0.Y1() - data.LCS0.Y2();
    const double y13 = data.LCS0.Y1() - data.LCS0.Y3();
    const double y23 = data.LCS0.Y2() - data.LCS0.Y3();
    const double y24 = data.LCS0.Y2() - data.LCS0.Y4();
    const double y34 = data.LCS0.Y3() - data.LCS0.Y4();
    const double y41 = data.LCS0.Y4() - data.LCS0.Y1();

    const double y21 = -y12;
    const double y31 = -y13;
    const double y32 = -y23;
    const double y42 = -y24;
    const double y43 = -y34;
    const double y14 = -y41;

    const double A = data.LCS0.Area();

    for (int i = 0; i < 4; i++) {
        data.r_cartesian[i] = Vector(3, 0.0);
    }
    data.r_cartesian[0] = data.LCS0.P1();
    data.r_cartesian[1] = data.LCS0.P2();
    data.r_cartesian[2] = data.LCS0.P3();
    data.r_cartesian[3] = data.LCS0.P4();

    //Precalculate dA to be multiplied with material matrix
    const GeometryType& geom = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points =
        geom.IntegrationPoints(mIntegrationMethod);
    data.dA.clear();
    for (int gp = 0; gp < 4; gp++) {
        //getting informations for integration
        double IntegrationWeight = integration_points[gp].Weight();

        // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
        // and Shape functions derivatives in the local coordinate system
        data.jacOp.Calculate(data.LCS0,
                             geom.ShapeFunctionLocalGradient(gp));

        // compute the 'area' of the current integration point
        data.dA[gp] = IntegrationWeight * data.jacOp.Determinant();
    }

    if (data.basicQuad == false) {
        // ---------------------------------------------------------------------
        //
        //                     ANDES MEMBRANE FORMULATION
        //
        // ---------------------------------------------------------------------


        // Unit vectors s_xi and s_eta (eqn 5.2.25)
        // eqns 5.2.29 -> 5.2.31
        data.s_xi.clear();
        data.s_eta.clear();

        //set values of SFs to xi = 1 and eta = 0
        ShellUtilities::ShapeFunc(1, 0, data.N);
        for (int i = 0; i < 4; i++) {
            data.s_xi(0) += data.r_cartesian[i][0] * data.N(i);
            data.s_xi(1) += data.r_cartesian[i][1] * data.N(i);
        }
        double s_xi_mag = std::sqrt(inner_prod(data.s_xi, data.s_xi));
        data.s_xi = data.s_xi / s_xi_mag;

        //set values of SFs to xi = 0 and eta = 1
        data.N.clear();
        ShellUtilities::ShapeFunc(0, 1, data.N);
        for (int i = 0; i < 4; i++) {
            data.s_eta(0) += data.r_cartesian[i][0] * data.N(i);
            data.s_eta(1) += data.r_cartesian[i][1] * data.N(i);
        }
        double s_eta_mag = std::sqrt(inner_prod(data.s_eta, data.s_eta));
        data.s_eta = data.s_eta / s_eta_mag;

        // calculate L - Lumping matrix (ref eqn 5.2.4)
        // for the construction of the basic
        // stiffness. Transpose from the presented version
        // to allow combination of B matrices later

        //Template constants
        double L_mult = 0.5 / A;
        const double alpha_6 = data.alpha / 6.0;
        const double alpha_3 = data.alpha / 3.0;
        data.L_mem.clear();

        //j = 1, i=4, k=2
        //ki = 24, ij = 41
        data.L_mem(0, 0) = L_mult * y24;
        data.L_mem(1, 0) = 0.00;
        data.L_mem(2, 0) = L_mult * x42;
        data.L_mem(0, 1) = 0.00;
        data.L_mem(1, 1) = L_mult * x42;
        data.L_mem(2, 1) = L_mult * y24;
        data.L_mem(0, 2) = L_mult * alpha_6*(y41*y41 - y21*y21);
        data.L_mem(1, 2) = L_mult *  alpha_6*(x41*x41 - x21*x21);
        data.L_mem(2, 2) = L_mult * alpha_3*(x21*y21 - x41*y41);

        //j = 2, i=1, k=3
        //ki = 31, ij = 12
        data.L_mem(0, 3) = L_mult * y31;
        data.L_mem(1, 3) = 0.00;
        data.L_mem(2, 3) = L_mult * x13;
        data.L_mem(0, 4) = 0.00;
        data.L_mem(1, 4) = L_mult * x13;
        data.L_mem(2, 4) = L_mult * y31;
        data.L_mem(0, 5) = L_mult * alpha_6*(y12*y12 - y32*y32);
        data.L_mem(1, 5) = L_mult * alpha_6*(x12*x12 - x32*x32);
        data.L_mem(2, 5) = L_mult * alpha_3*(x32*y32 - x12*y12);

        //j = 3, i=2, k=4
        //ki = 42, ij = 23
        data.L_mem(0, 6) = L_mult * y42;
        data.L_mem(1, 6) = 0.00;
        data.L_mem(2, 6) = L_mult * x24;
        data.L_mem(0, 7) = 0.00;
        data.L_mem(1, 7) = L_mult * x24;
        data.L_mem(2, 7) = L_mult * y42;
        data.L_mem(0, 8) = L_mult * alpha_6*(y23*y23 - y43*y43);
        data.L_mem(1, 8) = L_mult * alpha_6*(x23*x23 - x43*x43);
        data.L_mem(2, 8) = L_mult * alpha_3*(x43*y43 - x23*y23);

        //j = 4, i=3, k=1
        //ki = 13, ij = 34
        data.L_mem(0, 9) = L_mult * y13;
        data.L_mem(1, 9) = 0.00;
        data.L_mem(2, 9) = L_mult * x31;
        data.L_mem(0, 10) = 0.00;
        data.L_mem(1, 10) = L_mult * x31;
        data.L_mem(2, 10) = L_mult * y13;
        data.L_mem(0, 11) = L_mult * alpha_6*(y34*y34 - y14*y14);
        data.L_mem(1, 11) = L_mult * alpha_6*(x34*x34 - x14*x14);
        data.L_mem(2, 11) = L_mult * alpha_3*(x14*y14 - x34*y34);

        //--------------------------------------
        // calculate H - matrix
        // for the construction of the
        // higher order stiffness

        //H_mem_mod transformation matrix 'H' (ref eqn 5.2.26)

        //eqn 5.2.12
        double detJ = 0.0;
        for (int i = 0; i < 4; i++) {
            int j = i + 1;
            if (j == 4) {
                j = 0;
            }
            detJ += (data.r_cartesian[i][0] * data.r_cartesian[j][1] -
                     data.r_cartesian[j][0] * data.r_cartesian[i][1]);
        }
        detJ /= 8.0;
        const double f = 16.0 * detJ;

        //filter matrix for higher order rotation strain field
        //eqn 5.2.14
        Matrix H_theta = Matrix(5, 12, 0.0);
        H_theta(4, 0) = x42 / f;
        H_theta(4, 1) = x13 / f;
        H_theta(4, 2) = x24 / f;
        H_theta(4, 3) = x31 / f;
        H_theta(4, 4) = y42 / f;
        H_theta(4, 5) = y13 / f;
        H_theta(4, 6) = y24 / f;
        H_theta(4, 7) = y31 / f;

        H_theta(4, 8) = 0.25;
        H_theta(4, 9) = 0.25;
        H_theta(4, 10) = 0.25;
        H_theta(4, 11) = 0.25;

        for (int row = 0; row < 4; row++) {
            for (int col = 0; col < 4; col++) {
                if (row == col) {
                    H_theta(row, col + 8) = 0.75;
                } else {
                    H_theta(row, col + 8) = -0.25;
                }
            }
        }

        //xi and eta unit vectors
        Vector s_xi = Vector(data.r_cartesian[1] + data.r_cartesian[2]);
        s_xi /= std::sqrt(inner_prod(s_xi, s_xi));

        Vector s_eta = Vector(data.r_cartesian[2] + data.r_cartesian[3]);
        s_eta /= std::sqrt(inner_prod(s_eta, s_eta));

        //eqn 5.2.21 - modified as per Felippa supernatural quad paper
        //eqns 63, 65
        array_1d<double, 4> v_h;
        double A_0 = data.LCS0.Area();        //element area
        double A_1 = 0.5*(x34*y12 - x12*y34);
        double A_2 = 0.5*(x23*y14 - x14*y23);
        v_h[0] = (A_0 + A_1 + A_2) / 2.0 / A_0;
        v_h[1] = (-1.0*A_0 + A_1 - A_2) / 2.0 / A_0;
        v_h[2] = (A_0 - A_1 - A_2) / 2.0 / A_0;
        v_h[3] = (-1.0*A_0 - A_1 + A_2) / 2.0 / A_0;

        //filter matrix for higher order translations strain field
        Matrix H_tv = Matrix(2, 12, 0.0);
        for (int i = 0; i < 4; i++) {
            // x-components
            H_tv(0, i) = v_h[i] * s_xi[0];
            H_tv(1, i) = v_h[i] * s_eta[0];

            //y components
            H_tv(0, 4 + i) = v_h[i] * s_xi[1];
            H_tv(1, 4 + i) = v_h[i] * s_eta[1];
        }

        //eqn 5.2.26
        Matrix H = Matrix(7, 12, 0.0);
        for (int col = 0; col < 12; col++) {
            for (int row = 0; row < 5; row++) {
                H(row, col) = H_theta(row, col);
            }
            for (int row = 0; row < 2; row++) {
                H(row + 5, col) = H_tv(row, col);
            }
        }

        //Q matrices
        array_1d<double, 4> d_xi_i;
        array_1d<double, 4> d_eta_i;
        array_1d<double, 4> chi_xi_i;
        array_1d<double, 4> chi_eta_i;

        Vector r_xi = Vector(data.r_cartesian[1] + data.r_cartesian[2] -
                             data.r_cartesian[0] - data.r_cartesian[3]);
        r_xi /= 2.0;
        Vector r_eta = Vector(data.r_cartesian[2] + data.r_cartesian[3] -
                              data.r_cartesian[0] - data.r_cartesian[1]);
        r_eta /= 2.0;
        double l_xi = std::sqrt(inner_prod(r_xi, r_xi));
        double l_eta = std::sqrt(inner_prod(r_eta, r_eta));

        for (int i = 0; i < 4; i++) {
            //eqn 5.2.29
            const Vector vec1 = MathUtils<double>::CrossProduct(data.r_cartesian[i], s_xi);
            d_xi_i[i] = std::sqrt(inner_prod(vec1, vec1));
            chi_xi_i[i] = d_xi_i[i] / l_xi;

            const Vector vec2 = MathUtils<double>::CrossProduct(data.r_cartesian[i], s_eta);
            d_eta_i[i] = std::sqrt(inner_prod(vec2, vec2));
            chi_eta_i[i] = d_eta_i[i] / l_eta;
        }

        Vector r_24 = Vector(data.r_cartesian[1] - data.r_cartesian[3]);
        Vector r_13 = Vector(data.r_cartesian[0] - data.r_cartesian[2]);
        double l_24 = std::sqrt(inner_prod(r_24, r_24));
        double l_13 = std::sqrt(inner_prod(r_13, r_13));

        Vector e_24 = Vector(r_24 / l_24);
        const Vector vec1 = Vector(MathUtils<double>::CrossProduct(r_13, e_24));

        const double d_24 = std::sqrt(inner_prod(vec1, vec1));
        const double d_13 = d_24;
        const double chi_24 = d_24 / 2.0 / l_24;
        const double chi_13 = d_13 / 2.0 / l_13;

        const double chi_xi_t = l_eta / l_xi;
        const double chi_eta_t = l_xi / l_eta;

        double chi_xi_hat = 0.0;
        double chi_eta_hat = 0.0;
        for (int i = 0; i < 4; i++) {
            chi_xi_hat += chi_xi_i[i];
            chi_eta_hat += chi_eta_i[i];
        }
        chi_xi_hat /= 4.0;
        chi_eta_hat /= 4.0;

        // Template constants defined in eqn 5.2.41
        const double rho1 = 0.1;
        const double rho2 = -0.1;
        const double rho3 = -0.1;
        const double rho4 = 0.1;
        const double rho5 = 0.0;
        const double rho6 = 0.5;
        const double rho7 = 0.0;
        const double rho8 = -0.5;
        const double beta1 = 0.6;
        //double beta2 = 0.0; - entries disabled to save effort

        //s_13 and s_24 unit vectors
        Vector s_13 = Vector(data.r_cartesian[2] - data.r_cartesian[0]);
        s_13 /= std::sqrt(inner_prod(s_13, s_13));
        Vector s_24 = Vector(data.r_cartesian[3] - data.r_cartesian[1]);
        s_24 /= std::sqrt(inner_prod(s_24, s_24));

        Matrix Q1 = Matrix(3, 7, 0.0);
        Matrix Q2 = Matrix(3, 7, 0.0);
        Matrix Q3 = Matrix(3, 7, 0.0);
        Matrix Q4 = Matrix(3, 7, 0.0);

        Q1(0, 0) = rho1*chi_xi_i[0];
        Q1(0, 1) = rho2*chi_xi_i[0];
        Q1(0, 2) = rho3*chi_xi_i[0];
        Q1(0, 3) = rho4*chi_xi_i[0];

        Q1(0, 4) = data.alpha*chi_xi_t;
        Q1(0, 5) = -1.0 * beta1 * chi_xi_i[0] / chi_xi_hat / l_xi;

        Q1(1, 0) = -1.0* rho1*chi_eta_i[0];
        Q1(1, 1) = -1.0* rho4*chi_eta_i[0];
        Q1(1, 2) = -1.0* rho3*chi_eta_i[0];
        Q1(1, 3) = -1.0* rho2*chi_eta_i[0];

        Q1(1, 4) = -1.0 * data.alpha*chi_eta_t;
        Q1(1, 6) = -1.0 * beta1 * chi_eta_i[0] / chi_eta_hat / l_eta;

        Q1(2, 0) = rho5*chi_24;
        Q1(2, 1) = rho6*chi_24;
        Q1(2, 2) = rho7*chi_24;
        Q1(2, 3) = rho8*chi_24;

        //Q1(2, 5) = beta2 * c_24_xi / l_24;    - beta2 = 0!!!
        //Q1(2, 6) = -1.0 * beta2 * c_24_eta / l_24;    - beta2 = 0!!!

        Q2(0, 0) = -1.0*rho2*chi_xi_i[1];
        Q2(0, 1) = -1.0*rho1*chi_xi_i[1];
        Q2(0, 2) = -1.0*rho4*chi_xi_i[1];
        Q2(0, 3) = -1.0*rho3*chi_xi_i[1];

        Q2(0, 4) = -1.0*data.alpha*chi_xi_t;
        Q2(0, 5) = -1.0 * beta1 * chi_xi_i[1] / chi_xi_hat / l_xi;

        Q2(1, 0) = rho4*chi_eta_i[1];
        Q2(1, 1) = rho1*chi_eta_i[1];
        Q2(1, 2) = rho2*chi_eta_i[1];
        Q2(1, 3) = rho3*chi_eta_i[1];

        Q2(1, 4) = data.alpha*chi_eta_t;
        Q2(1, 6) = beta1 * chi_eta_i[1] / chi_eta_hat / l_eta;

        Q2(2, 0) = rho8*chi_13;
        Q2(2, 1) = rho5*chi_13;
        Q2(2, 2) = rho6*chi_13;
        Q2(2, 3) = rho7*chi_13;

        //Q2(2, 5) = -1.0* beta2 * c_13_xi / l_13;     beta2 = 0!!!
        //Q2(2, 6) = beta2 * c_13_eta / l_13;    beta2=0!!!


        Q3(0, 0) = rho3*chi_xi_i[2];
        Q3(0, 1) = rho4*chi_xi_i[2];
        Q3(0, 2) = rho1*chi_xi_i[2];
        Q3(0, 3) = rho2*chi_xi_i[2];

        Q3(0, 4) = data.alpha*chi_xi_t;
        Q3(0, 5) = beta1 * chi_xi_i[2] / chi_xi_hat / l_xi;

        Q3(1, 0) = -1.0* rho3*chi_eta_i[2];
        Q3(1, 1) = -1.0* rho2*chi_eta_i[2];
        Q3(1, 2) = -1.0* rho1*chi_eta_i[2];
        Q3(1, 3) = -1.0* rho4*chi_eta_i[2];

        Q3(1, 4) = -1.0 * data.alpha*chi_eta_t;
        Q3(1, 6) = beta1 * chi_eta_i[2] / chi_eta_hat / l_eta;

        Q3(2, 0) = rho7*chi_13;
        Q3(2, 1) = rho8*chi_13;
        Q3(2, 2) = rho5*chi_13;
        Q3(2, 3) = rho6*chi_13;

        //Q3(2, 5) = -1.0*beta2 * c_13_xi / l_13;    beta2 = 0!!!
        //Q3(2, 6) = beta2 * c_13_eta / l_13;    beta2 = 0!!!


        Q4(0, 0) = -1.0*rho4*chi_xi_i[3];
        Q4(0, 1) = -1.0*rho3*chi_xi_i[3];
        Q4(0, 2) = -1.0*rho2*chi_xi_i[3];
        Q4(0, 3) = -1.0*rho1*chi_xi_i[3];

        Q4(0, 4) = -1.0*data.alpha*chi_xi_t;
        Q4(0, 5) = beta1 * chi_xi_i[3] / chi_xi_hat / l_xi;

        Q4(1, 0) = rho2*chi_eta_i[3];
        Q4(1, 1) = rho3*chi_eta_i[3];
        Q4(1, 2) = rho4*chi_eta_i[3];
        Q4(1, 3) = rho1*chi_eta_i[3];

        Q4(1, 4) = data.alpha*chi_eta_t;
        Q4(1, 6) = -1.0* beta1 * chi_eta_i[3] / chi_eta_hat / l_eta;

        Q4(2, 0) = rho6*chi_13;
        Q4(2, 1) = rho7*chi_13;
        Q4(2, 2) = rho8*chi_13;
        Q4(2, 3) = rho5*chi_13;

        //Q4(2, 5) = beta2 * c_13_xi / l_13;    beta2 = 0
        //Q4(2, 6) = -1.0*beta2 * c_13_eta / l_13;    beta2 = 0

        Matrix T_13_inv = Matrix(3, 3, 0.0);
        T_13_inv(0, 0) = s_xi[0] * s_xi[0];
        T_13_inv(0, 1) = s_xi[1] * s_xi[1];
        T_13_inv(0, 2) = s_xi[0] * s_xi[1];
        T_13_inv(1, 0) = s_eta[0] * s_eta[0];
        T_13_inv(1, 1) = s_eta[1] * s_eta[1];
        T_13_inv(1, 2) = s_eta[0] * s_eta[1];
        T_13_inv(2, 0) = s_24[0] * s_24[0];
        T_13_inv(2, 1) = s_24[1] * s_24[1];
        T_13_inv(2, 2) = s_24[0] * s_24[1];

        Matrix T_24_inv = Matrix(3, 3, 0.0);
        T_24_inv(0, 0) = s_xi[0] * s_xi[0];
        T_24_inv(0, 1) = s_xi[1] * s_xi[1];
        T_24_inv(0, 2) = s_xi[0] * s_xi[1];
        T_24_inv(1, 0) = s_eta[0] * s_eta[0];
        T_24_inv(1, 1) = s_eta[1] * s_eta[1];
        T_24_inv(1, 2) = s_eta[0] * s_eta[1];
        T_24_inv(2, 0) = s_13[0] * s_13[0];
        T_24_inv(2, 1) = s_13[1] * s_13[1];
        T_24_inv(2, 2) = s_13[0] * s_13[1];

        Matrix T_13 = Matrix(3, 3, 0.0);
        Matrix T_24 = Matrix(3, 3, 0.0);
        double t13invdet = MathUtils<double>::Det(T_13_inv);
        MathUtils<double>::InvertMatrix(T_13_inv, T_13, t13invdet);
        double t24invdet = MathUtils<double>::Det(T_24_inv);
        MathUtils<double>::InvertMatrix(T_24_inv, T_24, t24invdet);

        data.B_h_1 = prod(T_13, Q1);
        data.B_h_2 = prod(T_24, Q2);
        data.B_h_3 = prod(T_13, Q3);
        data.B_h_4 = prod(T_24, Q4);

        //transform DOFs from Haugen to Kratos
        data.Z.clear();
        for (int i = 0; i < 3; i++) {
            data.Z(4 * i, i) = 1.0;
            data.Z(4 * i + 1, i + 3) = 1.0;
            data.Z(4 * i + 2, i + 6) = 1.0;
            data.Z(4 * i + 3, i + 9) = 1.0;
        }

        data.H_mem_mod.clear();
        data.H_mem_mod = prod(H, data.Z);

        //calculate Bh bar
        data.B_h_bar.clear();
        const Matrix& shapeFunctionsValues =
            geom.ShapeFunctionsValues(GetIntegrationMethod());
        for (int i = 0; i < 4; i++) {
            data.B_h_bar += shapeFunctionsValues(i, 0) * data.B_h_1;
            data.B_h_bar += shapeFunctionsValues(i, 1) * data.B_h_2;
            data.B_h_bar += shapeFunctionsValues(i, 2) * data.B_h_3;
            data.B_h_bar += shapeFunctionsValues(i, 3) * data.B_h_4;
        }
    }

    // ---------------------------------------------------------------------
    //
    //                       DKQ BENDING FORMULATION
    //
    // ---------------------------------------------------------------------

    //Calculate edge normal vectors for eqn 5.3.15
    //ref eqn 5.1.9 for calc of normal vector from edge vec
    Vector s_12 = Vector(data.LCS0.P1() - data.LCS0.P2());
    const double l_12 = std::sqrt(inner_prod(s_12, s_12));

    Vector s_23 = Vector(data.LCS0.P2() - data.LCS0.P3());
    const double l_23 = std::sqrt(inner_prod(s_23, s_23));

    Vector s_34 = Vector(data.LCS0.P3() - data.LCS0.P4());
    const double l_34 = std::sqrt(inner_prod(s_34, s_34));

    Vector s_41 = Vector(data.LCS0.P4() - data.LCS0.P1());
    const double l_41 = std::sqrt(inner_prod(s_41, s_41));

    Vector s_13 = Vector(data.LCS0.P1() - data.LCS0.P3());

    Vector s_24 = Vector(data.LCS0.P2() - data.LCS0.P4());

    //--------------------------------------
    // calculate DKQ bending stiffness

    data.DKQ_a.clear();
    data.DKQ_b.clear();
    data.DKQ_c.clear();
    data.DKQ_d.clear();
    data.DKQ_e.clear();

    //assemble a_k - eqn 3.86a : a[0] = a_5
    data.DKQ_a[0] = -1.0 * x12 / l_12 / l_12;
    data.DKQ_a[1] = -1.0 * x23 / l_23 / l_23;
    data.DKQ_a[2] = -1.0 * x34 / l_34 / l_34;
    data.DKQ_a[3] = -1.0 * x41 / l_41 / l_41;

    //assemble b_k - eqn 3.86b : b[0] = b_5
    data.DKQ_b[0] = 3.0 / 4.0 * x12 * y12 / l_12 / l_12;
    data.DKQ_b[1] = 3.0 / 4.0 * x23 * y23 / l_23 / l_23;
    data.DKQ_b[2] = 3.0 / 4.0 * x34 * y34 / l_34 / l_34;
    data.DKQ_b[3] = 3.0 / 4.0 * x41 * y41 / l_41 / l_41;

    //assemble c_k - eqn 3.86c : c[0] = c_5
    data.DKQ_c[0] = (x12 * x12 / 4.0 - y12 * y12 / 2.0) / l_12 / l_12;
    data.DKQ_c[1] = (x23 * x23 / 4.0 - y23 * y23 / 2.0) / l_23 / l_23;
    data.DKQ_c[2] = (x34 * x34 / 4.0 - y34 * y34 / 2.0) / l_34 / l_34;
    data.DKQ_c[3] = (x41 * x41 / 4.0 - y41 * y41 / 2.0) / l_41 / l_41;

    //assemble d_k - eqn 3.86d : d[0] = d_5
    data.DKQ_d[0] = -1.0 * y12 / l_12 / l_12;
    data.DKQ_d[1] = -1.0 * y23 / l_23 / l_23;
    data.DKQ_d[2] = -1.0 * y34 / l_34 / l_34;
    data.DKQ_d[3] = -1.0 * y41 / l_41 / l_41;

    //assemble e_k - eqn 3.86e : e[0] = e_5
    data.DKQ_e[0] = (-1.0 * x12 * x12 / 2.0 + y12 * y12 / 4.0) / l_12 /
                    l_12;
    data.DKQ_e[1] = (-1.0 * x23 * x23 / 2.0 + y23 * y23 / 4.0) / l_23 /
                    l_23;
    data.DKQ_e[2] = (-1.0 * x34 * x34 / 2.0 + y34 * y34 / 4.0) / l_34 /
                    l_34;
    data.DKQ_e[3] = (-1.0 * x41 * x41 / 2.0 + y41 * y41 / 4.0) / l_41 /
                    l_41;

    //prepare DKT assembly indices - for eqn 3.92 a->f
    data.DKQ_indices.clear();
    // 1st col = r, 2nd col = s
    data.DKQ_indices(0, 0) = 5;    //actual node number, not index!
    data.DKQ_indices(0, 1) = 8;
    data.DKQ_indices(1, 0) = 6;
    data.DKQ_indices(1, 1) = 5;
    data.DKQ_indices(2, 0) = 7;
    data.DKQ_indices(2, 1) = 6;
    data.DKQ_indices(3, 0) = 8;
    data.DKQ_indices(3, 1) = 7;

    //custom jacobian as per eqn 14
    const GeometryType::IntegrationPointsArrayType& integrationPoints =
        geom.IntegrationPoints(GetIntegrationMethod());

    for (int i = 0; i < 4; i++) {
        //set to current parametric integration location
        const array_1d<double,3>& r_coordinates = integrationPoints[i].Coordinates();
        double xi = r_coordinates[0];
        double eta = r_coordinates[1];

        Matrix DKQ_temp = Matrix(2, 2, 0.0);
        DKQ_temp(0, 0) = x21 + x34 + eta*(x12 + x34);
        DKQ_temp(0, 1) = y21 + y34 + eta*(y12 + y34);
        DKQ_temp(1, 0) = x32 + x41 + xi*(x12 + x34);
        DKQ_temp(1, 1) = y32 + y41 + xi*(y12 + y34);
        DKQ_temp = DKQ_temp / 4;
        double det = MathUtils<double>::Det(DKQ_temp);
        Matrix DKQ_temp_inv = Matrix(2, 2, 0.0);
        DKQ_temp_inv(0, 0) = DKQ_temp(1, 1);
        DKQ_temp_inv(1, 1) = DKQ_temp(0, 0);
        DKQ_temp_inv(0, 1) = -1.0 * DKQ_temp(0, 1);
        DKQ_temp_inv(1, 0) = -1.0 *DKQ_temp(1, 0);
        DKQ_temp_inv = DKQ_temp_inv / det;
        data.DKQ_invJac[i] = Matrix(DKQ_temp_inv);
    }

    //--------------------------------------
    // calculate the displacement vector
    // in global and local coordinate systems

    GetValuesVector(data.globalDisplacements);
    data.localDisplacements =
        mpCoordinateTransformation->CalculateLocalDisplacements(
            data.LCS, data.globalDisplacements);

    //--------------------------------------
    // Clear all auxiliary
    // matrices to be used later on
    // during the element integration.

    data.B.clear();
    data.D.clear();
    data.BTD.clear();
    data.generalizedStrains.clear();
    data.generalizedStresses.clear();


    //--------------------------------------
    // Initialize the section parameters

    data.SectionParameters.SetElementGeometry(GetGeometry());
    data.SectionParameters.SetMaterialProperties(GetProperties());
    data.SectionParameters.SetProcessInfo(data.CurrentProcessInfo);

    data.SectionParameters.SetGeneralizedStrainVector
    (data.generalizedStrains);
    data.SectionParameters.SetGeneralizedStressVector
    (data.generalizedStresses);
    data.SectionParameters.SetConstitutiveMatrix(data.D);

    Flags& options = data.SectionParameters.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, data.CalculateRHS);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,
                data.CalculateLHS);

    KRATOS_CATCH("")
}

void ShellThinElement3D4N::CalculateBMatrix(CalculationData& data)
{
    //---------------------------------------------
    // geom data
    const GeometryType& geom = GetGeometry();

    //---------------------------------------------
    // set to zero the total B matrix
    data.B.clear();
    Matrix B_mem = Matrix(3, 12, 0.0);

    // --------------------------------------------
    //     ANDES MEMBRANE FORMULATION
    //---------------------------------------------

    if (data.basicQuad == true) {
        // This section replaces the ANDES membrane B mat
        // with a basic unenhanced pure displacement-based
        // B-matrix. Only for comparison purposes!
        //
        // ref: ToP-1516-tutorial1.pdf

        std::cout << "Using basic membrane formulation!" << std::endl;

        const GeometryType::IntegrationPointsArrayType& integrationPoints =
            geom.IntegrationPoints(GetIntegrationMethod());

        //set to current parametric integration location
        const array_1d<double,3>& r_coordinates = integrationPoints[data.gpIndex].Coordinates();
        double xi = r_coordinates[0];
        double eta = r_coordinates[1];

        Matrix dN(4, 2, 0.0);

        // dN/dxi and dN/deta
        ShellUtilities::ShapeFunc_NaturalDerivatives(xi, eta, dN);
        Matrix temp = Matrix(prod(data.DKQ_invJac[data.gpIndex], trans(dN)));

        // dN/dx and dN/dy
        dN = trans(temp);
        for (int node = 0; node < 4; node++) {
            B_mem(0, 3 * node) = dN(node, 0);        //dx
            B_mem(1, 3 * node + 1) = dN(node, 1);    //dy
            B_mem(2, 3 * node) = B_mem(1, 3 * node + 1);
            B_mem(2, 3 * node + 1) = B_mem(0, 3 * node);
        }
    } else {
        //---------------------------------------------
        // membrane basic part L
        // already computed.it is constant over the
        // element

        //---------------------------------------------
        // membrane higher order part B_h
        // already calculated in calculate B_h_mats func
        Matrix B_h = Matrix(3, 7, 0.0);
        const Matrix& shapeFunctionsValues =
            geom.ShapeFunctionsValues(GetIntegrationMethod());
        B_h += shapeFunctionsValues(data.gpIndex, 0) * data.B_h_1;
        B_h += shapeFunctionsValues(data.gpIndex, 1) * data.B_h_2;
        B_h += shapeFunctionsValues(data.gpIndex, 2) * data.B_h_3;
        B_h += shapeFunctionsValues(data.gpIndex, 3) * data.B_h_4;
        B_h -= data.B_h_bar;

        //---------------------------------------------
        // combine membrane entries by transforming
        // B_h

        Matrix B_hH_mem = Matrix(prod(B_h, data.H_mem_mod));
        B_mem += data.L_mem + B_hH_mem;
    }

    // --------------------------------------------
    //     DKQ BENDING FORMULATION
    //---------------------------------------------
    Matrix B_bend_total = Matrix(3, 12, 0.0);

    Vector dpsiX_dxi = Vector(12, 0.0);
    Vector dpsiX_deta = Vector(12, 0.0);
    Vector dpsiY_dxi = Vector(12, 0.0);
    Vector dpsiY_deta = Vector(12, 0.0);
    const GeometryType::IntegrationPointsArrayType& integrationPoints =
        geom.IntegrationPoints(GetIntegrationMethod());

    //set to current parametric integration location
    const array_1d<double,3>& r_coordinates = integrationPoints[data.gpIndex].Coordinates();
    double xi = r_coordinates[0];
    double eta = r_coordinates[1];

    double ar, br, cr, dr, er;
    double as, bs, cs, ds, es;
    int r, s;
    for (int node = 0; node < 4; node++) {
        r = data.DKQ_indices(node, 0);    //actual node number r retrieved
        s = data.DKQ_indices(node, 1);    //actual node number s retrieved

        ar = data.DKQ_a[r - 5];    //converted to index, where node 5 = index 0
        br = data.DKQ_b[r - 5];
        cr = data.DKQ_c[r - 5];
        dr = data.DKQ_d[r - 5];
        er = data.DKQ_e[r - 5];

        as = data.DKQ_a[s - 5];
        bs = data.DKQ_b[s - 5];
        cs = data.DKQ_c[s - 5];
        ds = data.DKQ_d[s - 5];
        es = data.DKQ_e[s - 5];

        //eqn 3.92 a->f
        // d( ) / dxi
        // Compute vector of dPsi_x/dxi
        dpsiX_dxi[3 * node] = 1.5 *
                              (ar*ShellUtilities::dN_seren_dxi(r, xi, eta) -
                               as*ShellUtilities::dN_seren_dxi(s, xi, eta));

        dpsiX_dxi[3 * node + 1] = br*ShellUtilities::dN_seren_dxi(r, xi, eta) +
                                  bs*ShellUtilities::dN_seren_dxi(s, xi, eta);

        dpsiX_dxi[3 * node + 2] =
            ShellUtilities::dN_seren_dxi(node + 1, xi, eta) -
            cr*ShellUtilities::dN_seren_dxi(r, xi, eta) -
            cs*ShellUtilities::dN_seren_dxi(s, xi, eta);

        // Compute vector of dPsi_y/dxi
        dpsiY_dxi[3 * node] = 1.5 *
                              (dr*ShellUtilities::dN_seren_dxi(r, xi, eta) -
                               ds*ShellUtilities::dN_seren_dxi(s, xi, eta));

        dpsiY_dxi[3 * node + 1] = -1.0 *
                                  ShellUtilities::dN_seren_dxi(node + 1, xi, eta) +
                                  er*ShellUtilities::dN_seren_dxi(r, xi, eta) +
                                  es*ShellUtilities::dN_seren_dxi(s, xi, eta);

        dpsiY_dxi[3 * node + 2] = -1.0 * br*
                                  ShellUtilities::dN_seren_dxi(r, xi, eta) -
                                  bs*ShellUtilities::dN_seren_dxi(s, xi, eta);

        // d( ) / deta
        // Compute vector of dPsi_x/deta
        dpsiX_deta[3 * node] = 1.5 *
                               (ar*ShellUtilities::dN_seren_deta(r, xi, eta) -
                                as*ShellUtilities::dN_seren_deta(s, xi, eta));

        dpsiX_deta[3 * node + 1] = br*
                                   ShellUtilities::dN_seren_deta(r, xi, eta) +
                                   bs*ShellUtilities::dN_seren_deta(s, xi, eta);

        dpsiX_deta[3 * node + 2] =
            ShellUtilities::dN_seren_deta(node + 1, xi, eta) -
            cr*ShellUtilities::dN_seren_deta(r, xi, eta) -
            cs*ShellUtilities::dN_seren_deta(s, xi, eta);

        // Compute vector of dPsi_y/deta
        dpsiY_deta[3 * node] = 1.5 *
                               (dr*ShellUtilities::dN_seren_deta(r, xi, eta) -
                                ds*ShellUtilities::dN_seren_deta(s, xi, eta));

        dpsiY_deta[3 * node + 1] = -1.0 *
                                   ShellUtilities::dN_seren_deta(node + 1, xi, eta) +
                                   er*ShellUtilities::dN_seren_deta(r, xi, eta) +
                                   es*ShellUtilities::dN_seren_deta(s, xi, eta);

        dpsiY_deta[3 * node + 2] = -1.0 * br*
                                   ShellUtilities::dN_seren_deta(r, xi, eta) -
                                   bs*ShellUtilities::dN_seren_deta(s, xi, eta);
    }

    double j11, j12, j21, j22;
    j11 = data.DKQ_invJac[data.gpIndex](0, 0);
    j12 = data.DKQ_invJac[data.gpIndex](0, 1);
    j21 = data.DKQ_invJac[data.gpIndex](1, 0);
    j22 = data.DKQ_invJac[data.gpIndex](1, 1);

    //Assemble into B matrix
    Matrix B_bend_DKQ = Matrix(3, 12, 0.0);
    for (int col = 0; col < 12; col++) {
        B_bend_DKQ(0, col) += j11*dpsiX_dxi(col) + j12*dpsiX_deta(col);

        B_bend_DKQ(1, col) += j21*dpsiY_dxi(col) + j22*dpsiY_deta(col);

        B_bend_DKQ(2, col) += j11*dpsiY_dxi(col) + j12*dpsiY_deta(col) +
                              j21*dpsiX_dxi(col) + j22*dpsiX_deta(col);
    }

    B_bend_total.clear();
    B_bend_total += B_bend_DKQ;

    //---------------------------------------------
    // assemble the membrane contribution
    // and the bending contribution
    // into the combined B matrix

    for (int nodeid = 0; nodeid < 4; nodeid++) {
        int i = nodeid * 3;
        int j = nodeid * 6;

        // membrane map: [0,1,5] <- [0,1,2]

        data.B(0, j) = B_mem(0, i);
        data.B(0, j + 1) = B_mem(0, i + 1);
        data.B(0, j + 5) = B_mem(0, i + 2);

        data.B(1, j) = B_mem(1, i);
        data.B(1, j + 1) = B_mem(1, i + 1);
        data.B(1, j + 5) = B_mem(1, i + 2);

        data.B(2, j) = B_mem(2, i);
        data.B(2, j + 1) = B_mem(2, i + 1);
        data.B(2, j + 5) = B_mem(2, i + 2);

        // bending map: [2,3,4] <- [0,1,2]

        data.B(3, j + 2) = B_bend_total(0, i);
        data.B(3, j + 3) = B_bend_total(0, i + 1);
        data.B(3, j + 4) = B_bend_total(0, i + 2);

        data.B(4, j + 2) = B_bend_total(1, i);
        data.B(4, j + 3) = B_bend_total(1, i + 1);
        data.B(4, j + 4) = B_bend_total(1, i + 2);

        data.B(5, j + 2) = B_bend_total(2, i);
        data.B(5, j + 3) = B_bend_total(2, i + 1);
        data.B(5, j + 4) = B_bend_total(2, i + 2);
    }
}

void ShellThinElement3D4N::CalculateSectionResponse(CalculationData& data)
{
    const GeometryType& geom = GetGeometry();
    const Matrix& shapeFunctions = geom.ShapeFunctionsValues();
    Vector iN(shapeFunctions.size2());
    noalias(iN) = row(shapeFunctions, data.gpIndex);

    data.jacOp.Calculate(data.LCS0,
                         geom.ShapeFunctionLocalGradient(data.gpIndex));
    data.SectionParameters.SetShapeFunctionsDerivatives
    (data.jacOp.XYDerivatives());

    ShellCrossSection::Pointer& section = mSections[data.gpIndex];
    data.SectionParameters.SetShapeFunctionsValues(iN);
    data.SectionParameters.SetMaterialProperties(GetProperties());
    data.D.clear();
    section->CalculateSectionResponse(data.SectionParameters,
                                      ConstitutiveLaw::StressMeasure_PK2);
}

void ShellThinElement3D4N::CalculateGaussPointContribution
(CalculationData& data, MatrixType& LHS, VectorType& RHS)
{
    // calculate the total strain displ. matrix
    CalculateBMatrix(data);

    // calculate section response
    CalculateSectionResponse(data);

    // multiply the section tangent matrices and stress resultants by 'dA'
    data.D *= data.dA[data.gpIndex];

    // Add all contributions to the Stiffness Matrix
    data.BTD.clear();
    noalias(data.BTD) = prod(trans(data.B), data.D);
    noalias(LHS) += prod(data.BTD, data.B);
}

void ShellThinElement3D4N::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo,
                                        const bool CalculateStiffnessMatrixFlag,
                                        const bool CalculateResidualVectorFlag)
{
    // Resize the Left Hand Side if necessary,
    // and initialize it to Zero

    if ((rLeftHandSideMatrix.size1() != 24) ||
            (rLeftHandSideMatrix.size2() != 24)) {
        rLeftHandSideMatrix.resize(24, 24, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(24, 24);

    // Resize the Right Hand Side if necessary,
    // and initialize it to Zero

    if (rRightHandSideVector.size() != 24) {
        rRightHandSideVector.resize(24, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(24);

    // Compute the local coordinate system.
    ShellQ4_LocalCoordinateSystem localCoordinateSystem(
        mpCoordinateTransformation->CreateLocalCoordinateSystem());

    ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem());

    // Initialize common calculation variables
    CalculationData data(localCoordinateSystem, referenceCoordinateSystem,
                         rCurrentProcessInfo);
    data.CalculateLHS = CalculateStiffnessMatrixFlag;
    data.CalculateRHS = CalculateResidualVectorFlag;
    InitializeCalculationData(data);

    // Gauss Loop.
    for (SizeType i = 0; i < GetNumberOfGPs(); i++) {
        data.gpIndex = i;
        CalculateGaussPointContribution(data, rLeftHandSideMatrix,
                                        rRightHandSideVector);
    }

    // If basic membrane formulation is enabled - add drilling stiffness
    if (data.basicQuad == true) {
        double max_stiff = 0.0;
        for (int dof = 0; dof < 24; dof++) {
            if (rLeftHandSideMatrix(dof, dof) > max_stiff) {
                max_stiff = rLeftHandSideMatrix(dof, dof);
            }
        }
        for (int node = 0; node < 4; node++) {
            rLeftHandSideMatrix(6 * node + 5, 6 * node + 5) =
                max_stiff / 1000.0;
        }
    }

    // Add all contributions to the residual vector
    rRightHandSideVector -= prod(rLeftHandSideMatrix,
                                 data.localDisplacements);

    // Placeholders for extension into stability analysis
    //bool extractKm = false;
    //bool extractKg = false;

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

    // Add body forces contributions. This doesn't depend on the coordinate
    // system
    AddBodyForces(data, rRightHandSideVector);
}

bool ShellThinElement3D4N::
TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses
(const Variable<Matrix>& rVariable,
 std::vector<Matrix>& rValues,
 const ProcessInfo& rCurrentProcessInfo)
{
    // Check the required output
    int ijob = 0;
    bool bGlobal = false;
    CheckGeneralizedStressOrStrainOutput(rVariable, ijob, bGlobal);

    // quick return
    if (ijob == 0) {
        return false;
    }

    // resize output
    SizeType size = 4;
    if (rValues.size() != size) {
        rValues.resize(size);
    }

    // Compute the local coordinate system.
    ShellQ4_LocalCoordinateSystem localCoordinateSystem(
        mpCoordinateTransformation->CreateLocalCoordinateSystem());
    ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
        mpCoordinateTransformation->CreateReferenceCoordinateSystem());

    // Just to store the rotation matrix for visualization purposes
    Matrix R(6, 6);
    Matrix aux33(3, 3);

    // Initialize common calculation variables
    CalculationData data(localCoordinateSystem, referenceCoordinateSystem,
                         rCurrentProcessInfo);
    if (ijob > 2) {
        data.CalculateLHS = true; // calc constitutive mat for composites
    } else {
        data.CalculateLHS = false;
    }
    data.CalculateRHS = true;
    InitializeCalculationData(data);

    // Get the current displacements in global coordinate system and
    // transform to reference local system
    MatrixType Rdisp(24, 24);
    referenceCoordinateSystem.ComputeTotalRotationMatrix(Rdisp);
    if (referenceCoordinateSystem.IsWarped()) {
        MatrixType W(24, 24);
        referenceCoordinateSystem.ComputeTotalWarpageMatrix(W);
        Rdisp = prod(W, Rdisp);
    }
    data.localDisplacements = prod(Rdisp, data.globalDisplacements);

    // Gauss Loop
    for (unsigned int i = 0; i < size; i++) {
        // Compute all strain-displacement matrices
        data.gpIndex = i;
        CalculateBMatrix(data);

        // Calculate strain vectors in local coordinate system
        noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

        // Calculate the response of the Cross Section
        ShellCrossSection::Pointer& section = mSections[i];
        if (ijob > 2) {
            if (ijob > 7) {
                //Calculate lamina stresses
                CalculateLaminaStrains(data);
                CalculateLaminaStresses(data);
            } else {
                // calculate force resultants
                CalculateSectionResponse(data);

                if (ijob > 4) {
                    // Compute stresses
                    CalculateStressesFromForceResultants(data.generalizedStresses,
                                                         section->GetThickness(GetProperties()));
                }
            }
            DecimalCorrection(data.generalizedStresses);
        }

        // save the results
        DecimalCorrection(data.generalizedStrains);

        // now the results are in the element coordinate system.
        // if necessary, rotate the results in the section (local)
        // coordinate system
        if (section->GetOrientationAngle() != 0.0 && !bGlobal) {
            if (ijob > 7) {
                section->GetRotationMatrixForGeneralizedStresses(-(section->GetOrientationAngle()), R);
                for (unsigned int i = 0; i < data.rlaminateStresses.size(); i++) {
                    data.rlaminateStresses[i] = prod(R, data.rlaminateStresses[i]);
                }

                section->GetRotationMatrixForGeneralizedStrains(-(section->GetOrientationAngle()), R);
                for (unsigned int i = 0; i < data.rlaminateStrains.size(); i++) {
                    data.rlaminateStrains[i] = prod(R, data.rlaminateStrains[i]);
                }
            } else if (ijob > 2) {
                section->GetRotationMatrixForGeneralizedStresses
                (-(section->GetOrientationAngle()), R);
                data.generalizedStresses = prod(R, data.generalizedStresses);
            } else {
                section->GetRotationMatrixForGeneralizedStrains
                (-(section->GetOrientationAngle()), R);
                data.generalizedStrains = prod(R, data.generalizedStrains);
            }
        }

        Matrix& iValue = rValues[i];
        if (iValue.size1() != 3 || iValue.size2() != 3) {
            iValue.resize(3, 3, false);
        }

        if (ijob == 1) { // strains
            iValue(0, 0) = data.generalizedStrains(0);
            iValue(1, 1) = data.generalizedStrains(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(2);
            iValue(0, 2) = iValue(2, 0) = 0;
            iValue(1, 2) = iValue(2, 1) = 0;
        } else if (ijob == 2) { // curvatures
            iValue(0, 0) = data.generalizedStrains(3);
            iValue(1, 1) = data.generalizedStrains(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(5);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        } else if (ijob == 3) { // forces
            iValue(0, 0) = data.generalizedStresses(0);
            iValue(1, 1) = data.generalizedStresses(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(2);
            iValue(0, 2) = iValue(2, 0) = 0;
            iValue(1, 2) = iValue(2, 1) = 0;
        } else if (ijob == 4) { // moments
            iValue(0, 0) = data.generalizedStresses(3);
            iValue(1, 1) = data.generalizedStresses(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(5);
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        } else if (ijob == 5) { // SHELL_STRESS_TOP_SURFACE
            iValue(0, 0) = data.generalizedStresses(0) +
                           data.generalizedStresses(3);
            iValue(1, 1) = data.generalizedStresses(1) +
                           data.generalizedStresses(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.generalizedStresses[2] + data.generalizedStresses[5];
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        } else if (ijob == 6) { // SHELL_STRESS_MIDDLE_SURFACE
            iValue(0, 0) = data.generalizedStresses(0);
            iValue(1, 1) = data.generalizedStresses(1);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.generalizedStresses[2];
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        } else if (ijob == 7) { // SHELL_STRESS_BOTTOM_SURFACE
            iValue(0, 0) = data.generalizedStresses(0) -
                           data.generalizedStresses(3);
            iValue(1, 1) = data.generalizedStresses(1) -
                           data.generalizedStresses(4);
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.generalizedStresses[2] - data.generalizedStresses[5];
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        } else if (ijob == 8) { // SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE
            iValue(0, 0) =
                data.rlaminateStresses[data.rlaminateStresses.size() - 1][0];
            iValue(1, 1) =
                data.rlaminateStresses[data.rlaminateStresses.size() - 1][1];
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) =
                               data.rlaminateStresses[data.rlaminateStresses.size() - 1][2];
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        } else if (ijob == 9) { // SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE
            iValue(0, 0) = data.rlaminateStresses[0][0];
            iValue(1, 1) = data.rlaminateStresses[0][1];
            iValue(2, 2) = 0.0;
            iValue(0, 1) = iValue(1, 0) = data.rlaminateStresses[0][2];
            iValue(0, 2) = iValue(2, 0) = 0.0;
            iValue(1, 2) = iValue(2, 1) = 0.0;
        } else if (ijob == 99) { // SHELL_ORTHOTROPIC_4PLY_THROUGH_THICKNESS
            // Testing variable to get lamina stress/strain values
            // on the 8 surfaces of a 4 ply laminate

            int surface = 0; // start from top ply top surface
            // Output global results sequentially
            for (SizeType row = 0; row < 3; row++) {
                for (SizeType col = 0; col < 3; col++) {
                    if (surface > 7) {
                        iValue(row, col) = 0.0;
                    } else {
                        iValue(row, col) = data.rlaminateStresses[surface][2];
                    }
                    surface++;
                }
            }

            bool tsai_wu_thru_output = true;
            if (tsai_wu_thru_output) {
                // Must use non-global stresses for this!!!
                std::vector<Matrix> Laminae_Strengths =
                    std::vector<Matrix>(section->NumberOfPlies());
                for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++) {
                    Laminae_Strengths[ply].resize(3, 3, 0.0);
                    Laminae_Strengths[ply].clear();
                }
                const PropertiesType& props = GetProperties();
                section->GetLaminaeStrengths(Laminae_Strengths, props);
                Vector ply_orientation(section->NumberOfPlies());
                section->GetLaminaeOrientation(GetProperties(), ply_orientation);

                // Rotate lamina stress from section CS
                // to lamina angle to lamina material principal directions
                for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++) {
                    double total_rotation = -ply_orientation[ply]; // already rotated to section CS
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
                for (unsigned int ply = 0; ply < section->NumberOfPlies(); ply++) {
                    Vector lamina_stress_top = Vector(data.rlaminateStresses[2 * ply]);
                    Vector lamina_stress_bottom = Vector(data.rlaminateStresses[2 * ply + 1]);
                    // for testing, we want both top and bottom results, so trick the function

                    // top surface
                    data.rlaminateStresses[2 * ply + 1] = Vector(lamina_stress_top);
                    temp_tsai_wu = CalculateTsaiWuPlaneStress(data, Laminae_Strengths[ply], ply);
                    tsai_output[2 * ply] = temp_tsai_wu;

                    // bottom surface
                    data.rlaminateStresses[2 * ply] = Vector(lamina_stress_bottom);
                    data.rlaminateStresses[2 * ply + 1] = Vector(lamina_stress_bottom);
                    temp_tsai_wu = CalculateTsaiWuPlaneStress(data, Laminae_Strengths[ply], ply);
                    tsai_output[2 * ply+1] = temp_tsai_wu;
                }

                //dump into results

                int surface = 0; // start from top ply top surface
                // Output global results sequentially
                for (SizeType row = 0; row < 3; row++) {
                    for (SizeType col = 0; col < 3; col++) {
                        if (surface > 7) {
                            iValue(row, col) = 0.0;
                        } else {
                            iValue(row, col) = tsai_output[surface];
                        }
                        surface++;
                    }
                }
            }

        }

        // if requested, rotate the results in the global coordinate system
        if (bGlobal) {
            const Matrix& RG = referenceCoordinateSystem.Orientation();
            noalias(aux33) = prod(trans(RG), iValue);
            noalias(iValue) = prod(aux33, RG);
        }
    } // Gauss Loop

    return true;
}

// =========================================================================
//
// CalculationData
//
// =========================================================================

ShellThinElement3D4N::CalculationData::CalculationData
(const ShellQ4_LocalCoordinateSystem& localcoordsys,
 const ShellQ4_LocalCoordinateSystem& refcoordsys,
 const ProcessInfo& rCurrentProcessInfo)
    : LCS(localcoordsys)
    , LCS0(refcoordsys)
    , CurrentProcessInfo(rCurrentProcessInfo)
{
}

ShellCrossSection::SectionBehaviorType ShellThinElement3D4N::GetSectionBehavior() const
{
    return ShellCrossSection::Thin;
}

// =========================================================================
//
// Class ShellThinElement3D4N - Serialization
//
// =========================================================================

void ShellThinElement3D4N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseShellElement);
    bool is_corotational = (nullptr != dynamic_cast<ShellQ4_CorotationalCoordinateTransformation*>(mpCoordinateTransformation.get()));
    rSerializer.save("is_corotational", is_corotational);
    rSerializer.save("CTr", *mpCoordinateTransformation);
}

void ShellThinElement3D4N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseShellElement);
    bool is_corotational;
    rSerializer.load("is_corotational", is_corotational);
    if (is_corotational) {
        mpCoordinateTransformation = Kratos::make_shared<ShellQ4_CorotationalCoordinateTransformation>(pGetGeometry());
    } else {
        mpCoordinateTransformation = Kratos::make_shared<ShellQ4_CoordinateTransformation>(pGetGeometry());
    }
    rSerializer.load("CTr", *mpCoordinateTransformation);
}
}
