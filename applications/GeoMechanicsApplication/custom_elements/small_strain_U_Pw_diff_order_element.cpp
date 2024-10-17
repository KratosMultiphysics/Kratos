// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// Project includes
#include "containers/array_1d.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "includes/cfd_variables.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_retention/retention_law_factory.h"
#include "custom_utilities/constitutive_law_utilities.hpp"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "stress_state_policy.h"

namespace Kratos
{
Element::Pointer SmallStrainUPwDiffOrderElement::Create(IndexType               NewId,
                                                        NodesArrayType const&   ThisNodes,
                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallStrainUPwDiffOrderElement(
        NewId, GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

Element::Pointer SmallStrainUPwDiffOrderElement::Create(IndexType               NewId,
                                                        GeometryType::Pointer   pGeom,
                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallStrainUPwDiffOrderElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

int SmallStrainUPwDiffOrderElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();

    if (rGeom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    // check pressure geometry pointer
    KRATOS_DEBUG_ERROR_IF_NOT(mpPressureGeometry) << "Pressure Geometry is not defined\n";

    // verify that the variables are correctly initialized
    // Verify specific properties
    const PropertiesType& rProp = this->GetProperties();

    if (!rProp.Has(IGNORE_UNDRAINED))
        KRATOS_ERROR << "IGNORE_UNDRAINED does not exist in the parameter list" << this->Id() << std::endl;

    if (!rProp[IGNORE_UNDRAINED]) {
        if (!rProp.Has(PERMEABILITY_XX) || rProp[PERMEABILITY_XX] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_XX has Key zero, is not defined or "
                            "has an invalid value at element"
                         << this->Id() << std::endl;

        if (!rProp.Has(PERMEABILITY_YY) || rProp[PERMEABILITY_YY] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_YY has Key zero, is not defined or "
                            "has an invalid value at element"
                         << this->Id() << std::endl;

        if (!rProp.Has(PERMEABILITY_XY) || rProp[PERMEABILITY_XY] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_XY has Key zero, is not defined or "
                            "has an invalid value at element"
                         << this->Id() << std::endl;

        if (rGeom.WorkingSpaceDimension() > 2) {
            if (!rProp.Has(PERMEABILITY_ZZ) || rProp[PERMEABILITY_ZZ] < 0.0)
                KRATOS_ERROR << "PERMEABILITY_ZZ has Key zero, is not defined "
                                "or has an invalid value at element"
                             << this->Id() << std::endl;

            if (!rProp.Has(PERMEABILITY_YZ) || rProp[PERMEABILITY_YZ] < 0.0)
                KRATOS_ERROR << "PERMEABILITY_YZ has Key zero, is not defined "
                                "or has an invalid value at element"
                             << this->Id() << std::endl;

            if (!rProp.Has(PERMEABILITY_ZX) || rProp[PERMEABILITY_ZX] < 0.0)
                KRATOS_ERROR << "PERMEABILITY_ZX has Key zero, is not defined "
                                "or has an invalid value at element"
                             << this->Id() << std::endl;
        }
    }

    // verify that the dofs exist
    for (unsigned int i = 0; i < rGeom.size(); ++i) {
        if (!rGeom[i].SolutionStepsDataHas(DISPLACEMENT))
            KRATOS_ERROR << "missing variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;

        if (!rGeom[i].HasDofFor(DISPLACEMENT_X) || !rGeom[i].HasDofFor(DISPLACEMENT_Y) ||
            !rGeom[i].HasDofFor(DISPLACEMENT_Z))
            KRATOS_ERROR << "missing one of the dofs for the variable "
                            "DISPLACEMENT on node "
                         << rGeom[i].Id() << std::endl;

        if (!rGeom[i].SolutionStepsDataHas(WATER_PRESSURE))
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << rGeom[i].Id() << std::endl;

        if (!rGeom[i].HasDofFor(WATER_PRESSURE))
            KRATOS_ERROR << "missing the dof for the variable WATER_PRESSURE "
                            "on node "
                         << rGeom[i].Id() << std::endl;
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(rProp.Has(CONSTITUTIVE_LAW))
        << "Constitutive law not provided for property " << rProp.Id() << std::endl;

    // verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    rProp.GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for (unsigned int i = 0; i < LawFeatures.mStrainMeasures.size(); ++i) {
        if (LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
            correct_strain_measure = true;
    }

    if (!correct_strain_measure)
        KRATOS_ERROR << "constitutive law is not compatible with the element "
                        "type StrainMeasure_Infinitesimal "
                     << this->Id() << std::endl;

    rProp.GetValue(CONSTITUTIVE_LAW)->Check(rProp, rGeom, rCurrentProcessInfo);

    // Verify that the constitutive law has the correct dimension
    const SizeType strainSize = this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    if (rGeom.WorkingSpaceDimension() > 2) {
        KRATOS_ERROR_IF_NOT(strainSize == VOIGT_SIZE_3D)
            << "Wrong constitutive law used. This is a 3D element! expected "
               "strain size is "
            << VOIGT_SIZE_3D << " But received: " << strainSize << " in element id: " << this->Id()
            << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strainSize == VOIGT_SIZE_2D_PLANE_STRAIN)
            << "Wrong constitutive law used. This is a 2D element! expected "
               "strain size is "
            << VOIGT_SIZE_2D_PLANE_STRAIN << " But received: " << strainSize
            << " in element id: " << this->Id() << std::endl;
    }

    return 0;

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);

    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); // Note: this is for nonlocal damage

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    RetentionLaw::Parameters RetentionParameters(GetProperties());

    const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(deformation_gradients);
    const auto strain_vectors = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        GetStressStatePolicy().GetVoigtSize());

    const auto number_of_integration_points = GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
    for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);
        Variables.B            = b_matrices[GPoint];
        Variables.F            = deformation_gradients[GPoint];
        Variables.StrainVector = strain_vectors[GPoint];

        ConstitutiveLawUtilities::SetConstitutiveParameters(
            ConstitutiveParameters, Variables.StrainVector, Variables.ConstitutiveMatrix, Variables.Nu,
            Variables.DNu_DX, Variables.F, determinants_of_deformation_gradients[GPoint]);

        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->InitializeMaterialResponseCauchy(ConstitutiveParameters);

        mRetentionLawVector[GPoint]->InitializeSolutionStep(RetentionParameters);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geom             = GetGeometry();
    const auto          integration_method = this->GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geom.IntegrationPoints(integration_method);
    const auto Np_container = mpPressureGeometry->ShapeFunctionsValues(integration_method);

    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Np_container, this->GetPressureSolutionVector());
    const auto degrees_saturation = this->CalculateDegreesOfSaturation(fluid_pressures);

    const auto solid_densities =
        GeoTransportEquationUtilities::CalculateSoilDensities(degrees_saturation, GetProperties());

    const auto det_Js_initial_configuration =
        GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(r_geom, integration_method);

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(integration_points, det_Js_initial_configuration);

    const auto mass_matrix_u = GeoEquationOfMotionUtilities::CalculateMassMatrix(
        r_geom.WorkingSpaceDimension(), r_geom.PointsNumber(), integration_points.size(),
        r_geom.ShapeFunctionsValues(integration_method), solid_densities, integration_coefficients);

    rMassMatrix = ZeroMatrix(GetNumberOfDOF(), GetNumberOfDOF());
    GeoElementUtilities::AssembleUUBlockMatrix(rMassMatrix, mass_matrix_u);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    RetentionLaw::Parameters RetentionParameters(GetProperties());

    const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(deformation_gradients);
    const auto strain_vectors = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        GetStressStatePolicy().GetVoigtSize());

    const auto number_of_integration_points = GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
    for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);
        Variables.B            = b_matrices[GPoint];
        Variables.F            = deformation_gradients[GPoint];
        Variables.StrainVector = strain_vectors[GPoint];

        ConstitutiveLawUtilities::SetConstitutiveParameters(
            ConstitutiveParameters, Variables.StrainVector, Variables.ConstitutiveMatrix, Variables.Nu,
            Variables.DNu_DX, Variables.F, determinants_of_deformation_gradients[GPoint]);

        // Compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
        mStateVariablesFinalized[GPoint] =
            mConstitutiveLawVector[GPoint]->GetValue(STATE_VARIABLES, mStateVariablesFinalized[GPoint]);

        mRetentionLawVector[GPoint]->FinalizeSolutionStep(RetentionParameters);
    }

    // Assign pressure values to the intermediate nodes for post-processing
    if (!GetProperties()[IGNORE_UNDRAINED]) AssignPressureToIntermediateNodes();

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType&  rGeom     = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumDim    = rGeom.WorkingSpaceDimension();

    switch (NumUNodes) {
    case 6: // 2D T6P3
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
        ThreadSafeNodeWrite(rGeom[3], WATER_PRESSURE, 0.5 * (p0 + p1));
        ThreadSafeNodeWrite(rGeom[4], WATER_PRESSURE, 0.5 * (p1 + p2));
        ThreadSafeNodeWrite(rGeom[5], WATER_PRESSURE, 0.5 * (p2 + p0));
        break;
    }
    case 8: // 2D Q8P4
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
        ThreadSafeNodeWrite(rGeom[4], WATER_PRESSURE, 0.5 * (p0 + p1));
        ThreadSafeNodeWrite(rGeom[5], WATER_PRESSURE, 0.5 * (p1 + p2));
        ThreadSafeNodeWrite(rGeom[6], WATER_PRESSURE, 0.5 * (p2 + p3));
        ThreadSafeNodeWrite(rGeom[7], WATER_PRESSURE, 0.5 * (p3 + p0));
        break;
    }
    case 9: // 2D Q9P4
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
        ThreadSafeNodeWrite(rGeom[4], WATER_PRESSURE, 0.5 * (p0 + p1));
        ThreadSafeNodeWrite(rGeom[5], WATER_PRESSURE, 0.5 * (p1 + p2));
        ThreadSafeNodeWrite(rGeom[6], WATER_PRESSURE, 0.5 * (p2 + p3));
        ThreadSafeNodeWrite(rGeom[7], WATER_PRESSURE, 0.5 * (p3 + p0));
        ThreadSafeNodeWrite(rGeom[8], WATER_PRESSURE, 0.25 * (p0 + p1 + p2 + p3));
        break;
    }
    case 10: // 3D T10P4  //2D T10P6
    {
        if (NumDim == 3) {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4], WATER_PRESSURE, 0.5 * (p0 + p1));
            ThreadSafeNodeWrite(rGeom[5], WATER_PRESSURE, 0.5 * (p1 + p2));
            ThreadSafeNodeWrite(rGeom[6], WATER_PRESSURE, 0.5 * (p2 + p0));
            ThreadSafeNodeWrite(rGeom[7], WATER_PRESSURE, 0.5 * (p0 + p3));
            ThreadSafeNodeWrite(rGeom[8], WATER_PRESSURE, 0.5 * (p1 + p3));
            ThreadSafeNodeWrite(rGeom[9], WATER_PRESSURE, 0.5 * (p2 + p3));
        } else if (NumDim == 2) {
            constexpr double c1 = 1.0 / 9.0;
            const double     p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double     p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double     p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double     p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            const double     p4 = rGeom[4].FastGetSolutionStepValue(WATER_PRESSURE);
            const double     p5 = rGeom[5].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[0], WATER_PRESSURE, p0);
            ThreadSafeNodeWrite(rGeom[1], WATER_PRESSURE, p1);
            ThreadSafeNodeWrite(rGeom[2], WATER_PRESSURE, p2);
            ThreadSafeNodeWrite(rGeom[3], WATER_PRESSURE, (2.0 * p0 - p1 + 8.0 * p3) * c1);
            ThreadSafeNodeWrite(rGeom[4], WATER_PRESSURE, (2.0 * p1 - p0 + 8.0 * p3) * c1);
            ThreadSafeNodeWrite(rGeom[5], WATER_PRESSURE, (2.0 * p1 - p2 + 8.0 * p4) * c1);
            ThreadSafeNodeWrite(rGeom[6], WATER_PRESSURE, (2.0 * p2 - p1 + 8.0 * p4) * c1);
            ThreadSafeNodeWrite(rGeom[7], WATER_PRESSURE, (2.0 * p2 - p0 + 8.0 * p5) * c1);
            ThreadSafeNodeWrite(rGeom[8], WATER_PRESSURE, (2.0 * p0 - p2 + 8.0 * p5) * c1);
            ThreadSafeNodeWrite(rGeom[9], WATER_PRESSURE, (4.0 * (p3 + p4 + p5) - (p0 + p1 + p2)) * c1);
        }
        break;
    }
    case 15: // 2D T15P10
    {
        constexpr double c1 = 0.0390625;
        const double     p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p4 = rGeom[4].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p5 = rGeom[5].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p6 = rGeom[6].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p7 = rGeom[7].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p8 = rGeom[8].FastGetSolutionStepValue(WATER_PRESSURE);
        const double     p9 = rGeom[9].FastGetSolutionStepValue(WATER_PRESSURE);
        ThreadSafeNodeWrite(rGeom[0], WATER_PRESSURE, p0);
        ThreadSafeNodeWrite(rGeom[1], WATER_PRESSURE, p1);
        ThreadSafeNodeWrite(rGeom[2], WATER_PRESSURE, p2);
        ThreadSafeNodeWrite(rGeom[3], WATER_PRESSURE, (3.0 * p0 + p1 + 27.0 * p3 - 5.4 * p4) * c1);
        ThreadSafeNodeWrite(rGeom[4], WATER_PRESSURE, (14.4 * (p3 + p4) - 1.6 * (p0 + p1)) * c1);
        ThreadSafeNodeWrite(rGeom[5], WATER_PRESSURE, (3.0 * p1 + p0 + 27.0 * p4 - 5.4 * p3) * c1);
        ThreadSafeNodeWrite(rGeom[6], WATER_PRESSURE, (3.0 * p1 + p2 + 27.0 * p5 - 5.4 * p6) * c1);
        ThreadSafeNodeWrite(rGeom[7], WATER_PRESSURE, (14.4 * (p5 + p6) - 1.6 * (p1 + p2)) * c1);
        ThreadSafeNodeWrite(rGeom[8], WATER_PRESSURE, (3.0 * p2 + p1 + 27.0 * p6 - 5.4 * p5) * c1);
        ThreadSafeNodeWrite(rGeom[9], WATER_PRESSURE, (3.0 * p2 + p0 + 27.0 * p7 - 5.4 * p8) * c1);
        ThreadSafeNodeWrite(rGeom[10], WATER_PRESSURE, (14.4 * (p7 + p8) - 1.6 * (p0 + p2)) * c1);
        ThreadSafeNodeWrite(rGeom[11], WATER_PRESSURE, (3.0 * p0 + p2 + 27.0 * p8 - 5.4 * p7) * c1);
        ThreadSafeNodeWrite(
            rGeom[12], WATER_PRESSURE,
            (p1 + p2 + 7.2 * (p3 + p8) - 3.6 * (p4 + p7) - 1.8 * (p5 + p6) + 21.6 * p9 - 1.6 * p0) * c1);
        ThreadSafeNodeWrite(
            rGeom[13], WATER_PRESSURE,
            (p0 + p2 + 7.2 * (p4 + p5) - 3.6 * (p3 + p6) - 1.8 * (p7 + p8) + 21.6 * p9 - 1.6 * p1) * c1);
        ThreadSafeNodeWrite(
            rGeom[14], WATER_PRESSURE,
            (p0 + p1 + 7.2 * (p6 + p7) - 3.6 * (p5 + p8) - 1.8 * (p3 + p4) + 21.6 * p9 - 1.6 * p2) * c1);
        break;
    }
    case 20: // 3D H20P8
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p4 = rGeom[4].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p5 = rGeom[5].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p6 = rGeom[6].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p7 = rGeom[7].FastGetSolutionStepValue(WATER_PRESSURE);
        // edges -- bottom
        ThreadSafeNodeWrite(rGeom[8], WATER_PRESSURE, 0.5 * (p0 + p1));
        ThreadSafeNodeWrite(rGeom[9], WATER_PRESSURE, 0.5 * (p1 + p2));
        ThreadSafeNodeWrite(rGeom[10], WATER_PRESSURE, 0.5 * (p2 + p3));
        ThreadSafeNodeWrite(rGeom[11], WATER_PRESSURE, 0.5 * (p3 + p0));
        // edges -- middle
        ThreadSafeNodeWrite(rGeom[12], WATER_PRESSURE, 0.5 * (p4 + p0));
        ThreadSafeNodeWrite(rGeom[13], WATER_PRESSURE, 0.5 * (p5 + p1));
        ThreadSafeNodeWrite(rGeom[14], WATER_PRESSURE, 0.5 * (p6 + p2));
        ThreadSafeNodeWrite(rGeom[15], WATER_PRESSURE, 0.5 * (p7 + p3));
        // edges -- top
        ThreadSafeNodeWrite(rGeom[16], WATER_PRESSURE, 0.5 * (p4 + p5));
        ThreadSafeNodeWrite(rGeom[17], WATER_PRESSURE, 0.5 * (p5 + p6));
        ThreadSafeNodeWrite(rGeom[18], WATER_PRESSURE, 0.5 * (p6 + p7));
        ThreadSafeNodeWrite(rGeom[19], WATER_PRESSURE, 0.5 * (p7 + p4));
        break;
    }
    case 27: // 3D H27P8
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p4 = rGeom[4].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p5 = rGeom[5].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p6 = rGeom[6].FastGetSolutionStepValue(WATER_PRESSURE);
        const double p7 = rGeom[7].FastGetSolutionStepValue(WATER_PRESSURE);
        // edges -- bottom
        ThreadSafeNodeWrite(rGeom[8], WATER_PRESSURE, 0.5 * (p0 + p1));
        ThreadSafeNodeWrite(rGeom[9], WATER_PRESSURE, 0.5 * (p1 + p2));
        ThreadSafeNodeWrite(rGeom[10], WATER_PRESSURE, 0.5 * (p2 + p3));
        ThreadSafeNodeWrite(rGeom[11], WATER_PRESSURE, 0.5 * (p3 + p0));
        // edges -- middle
        ThreadSafeNodeWrite(rGeom[12], WATER_PRESSURE, 0.5 * (p4 + p0));
        ThreadSafeNodeWrite(rGeom[13], WATER_PRESSURE, 0.5 * (p5 + p1));
        ThreadSafeNodeWrite(rGeom[14], WATER_PRESSURE, 0.5 * (p6 + p2));
        ThreadSafeNodeWrite(rGeom[15], WATER_PRESSURE, 0.5 * (p7 + p3));
        // edges -- top
        ThreadSafeNodeWrite(rGeom[16], WATER_PRESSURE, 0.5 * (p4 + p5));
        ThreadSafeNodeWrite(rGeom[17], WATER_PRESSURE, 0.5 * (p5 + p6));
        ThreadSafeNodeWrite(rGeom[18], WATER_PRESSURE, 0.5 * (p6 + p7));
        ThreadSafeNodeWrite(rGeom[19], WATER_PRESSURE, 0.5 * (p7 + p0));
        // face centers
        ThreadSafeNodeWrite(rGeom[20], WATER_PRESSURE, 0.25 * (p0 + p1 + p2 + p3));
        ThreadSafeNodeWrite(rGeom[21], WATER_PRESSURE, 0.25 * (p0 + p1 + p4 + p5));
        ThreadSafeNodeWrite(rGeom[22], WATER_PRESSURE, 0.25 * (p1 + p2 + p5 + p6));
        ThreadSafeNodeWrite(rGeom[23], WATER_PRESSURE, 0.25 * (p2 + p3 + p6 + p7));
        ThreadSafeNodeWrite(rGeom[24], WATER_PRESSURE, 0.25 * (p3 + p0 + p7 + p4));
        ThreadSafeNodeWrite(rGeom[25], WATER_PRESSURE, 0.25 * (p4 + p5 + p6 + p7));
        // element center
        ThreadSafeNodeWrite(rGeom[26], WATER_PRESSURE, 0.125 * (p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7));
        break;
    }
    default:
        KRATOS_ERROR << "Unexpected geometry type for different order "
                        "interpolation element"
                     << this->Id() << std::endl;
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                  const std::vector<Vector>& rValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == CAUCHY_STRESS_VECTOR) {
        KRATOS_ERROR_IF(rValues.size() != mStressVector.size())
            << "Unexpected number of values for "
               "SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints"
            << std::endl;
        std::copy(rValues.begin(), rValues.end(), mStressVector.begin());
    } else {
        KRATOS_ERROR_IF(rValues.size() < mConstitutiveLawVector.size())
            << "Insufficient number of values for "
               "SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints"
            << std::endl;
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            mConstitutiveLawVector[GPoint]->SetValue(rVariable, rValues[GPoint], rCurrentProcessInfo);
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
                                                                  std::vector<int>&    rValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto number_of_integration_points =
        GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

    rValues.resize(number_of_integration_points);
    for (auto i = SizeType{0}; i < number_of_integration_points; ++i) {
        rValues[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                  std::vector<double>&    rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const auto number_of_integration_points = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

    rOutput.resize(number_of_integration_points);

    if (rVariable == VON_MISES_STRESS) {
        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStress(mStressVector[GPoint]);
        }
    } else if (rVariable == MEAN_EFFECTIVE_STRESS) {
        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateMeanStress(mStressVector[GPoint]);
        }
    } else if (rVariable == MEAN_STRESS) {
        std::vector<Vector> StressVector;
        CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateMeanStress(StressVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_VON_MISES_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStrain(StrainVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_VOLUMETRIC_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateTrace(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VON_MISES_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStrain(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VOLUMETRIC_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateTrace(StrainVector[GPoint]);
        }
    } else if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION ||
               rVariable == BISHOP_COEFFICIENT || rVariable == DERIVATIVE_OF_SATURATION ||
               rVariable == RELATIVE_PERMEABILITY) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        RetentionLaw::Parameters RetentionParameters(GetProperties());

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            // Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            RetentionParameters.SetFluidPressure(GeoTransportEquationUtilities::CalculateFluidPressure(
                Variables.Np, Variables.PressureVector));

            if (rVariable == DEGREE_OF_SATURATION)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateSaturation(RetentionParameters);
            else if (rVariable == EFFECTIVE_SATURATION)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateEffectiveSaturation(RetentionParameters);
            else if (rVariable == BISHOP_COEFFICIENT)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(RetentionParameters);
            else if (rVariable == DERIVATIVE_OF_SATURATION)
                rOutput[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(RetentionParameters);
            else if (rVariable == RELATIVE_PERMEABILITY)
                rOutput[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);
        }
    } else if (rVariable == HYDRAULIC_HEAD) {
        const double          NumericalLimit = std::numeric_limits<double>::epsilon();
        const PropertiesType& rProp          = this->GetProperties();

        // Defining the shape functions, the Jacobian and the shape functions local gradients Containers
        const Matrix&  NContainer = rGeom.ShapeFunctionsValues(this->GetIntegrationMethod());
        const SizeType NumUNodes  = rGeom.PointsNumber();

        // Defining necessary variables
        Vector NodalHydraulicHead = ZeroVector(NumUNodes);
        for (unsigned int node = 0; node < NumUNodes; ++node) {
            Vector NodeVolumeAcceleration(3);
            noalias(NodeVolumeAcceleration) = rGeom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
            const double g = norm_2(NodeVolumeAcceleration);
            if (g > NumericalLimit) {
                const double FluidWeight = g * rProp[DENSITY_WATER];

                Vector NodeCoordinates(3);
                noalias(NodeCoordinates) = rGeom[node].Coordinates();
                Vector NodeVolumeAccelerationUnitVector(3);
                noalias(NodeVolumeAccelerationUnitVector) = NodeVolumeAcceleration / g;

                const double WaterPressure = rGeom[node].FastGetSolutionStepValue(WATER_PRESSURE);
                NodalHydraulicHead[node] = -inner_prod(NodeCoordinates, NodeVolumeAccelerationUnitVector) -
                                           PORE_PRESSURE_SIGN_FACTOR * WaterPressure / FluidWeight;
            } else {
                NodalHydraulicHead[node] = 0.0;
            }
        }

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; GPoint++) {
            double HydraulicHead = 0.0;
            for (unsigned int node = 0; node < NumUNodes; ++node)
                HydraulicHead += NContainer(GPoint, node) * NodalHydraulicHead[node];

            rOutput[GPoint] = HydraulicHead;
        }
    } else {
        for (unsigned int i = 0; i < number_of_integration_points; ++i) {
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                  std::vector<array_1d<double, 3>>& rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const auto number_of_integration_points = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    rOutput.resize(number_of_integration_points);

    if (rVariable == FLUID_FLUX_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        const auto strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
        auto relative_permeability_values =
            CalculateRelativePermeabilityValues(GeoTransportEquationUtilities::CalculateFluidPressures(
                Variables.NpContainer, Variables.PressureVector));
        const auto permeability_update_factors =
            GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(strain_vectors, GetProperties());
        std::transform(relative_permeability_values.cbegin(), relative_permeability_values.cend(),
                       permeability_update_factors.cbegin(), relative_permeability_values.begin(),
                       std::multiplies<>{});

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            // compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables, GPoint);
            Variables.B = b_matrices[GPoint];

            // Compute FluidFlux vector q [m/s]
            const SizeType Dim       = r_geometry.WorkingSpaceDimension();
            const SizeType NumUNodes = r_geometry.PointsNumber();

            Vector   BodyAcceleration = ZeroVector(Dim);
            SizeType Index            = 0;
            for (SizeType i = 0; i < NumUNodes; ++i) {
                for (unsigned int idim = 0; idim < Dim; ++idim)
                    BodyAcceleration[idim] += Variables.Nu[i] * Variables.BodyAcceleration[Index++];
            }

            const auto relative_permeability = relative_permeability_values[GPoint];

            // Compute strain, need to update porosity
            Variables.F            = deformation_gradients[GPoint];
            Variables.StrainVector = strain_vectors[GPoint];

            Vector GradPressureTerm(Dim);
            noalias(GradPressureTerm) = prod(trans(Variables.DNp_DX), Variables.PressureVector);
            noalias(GradPressureTerm) +=
                PORE_PRESSURE_SIGN_FACTOR * GetProperties()[DENSITY_WATER] * BodyAcceleration;

            Vector AuxFluidFlux = ZeroVector(Dim);
            AuxFluidFlux        = PORE_PRESSURE_SIGN_FACTOR * Variables.DynamicViscosityInverse *
                           relative_permeability * prod(Variables.IntrinsicPermeability, GradPressureTerm);

            Vector FluidFlux = ZeroVector(3);
            for (unsigned int idim = 0; idim < Dim; ++idim)
                FluidFlux[idim] = AuxFluidFlux[idim];

            if (rOutput[GPoint].size() != 3) rOutput[GPoint].resize(3, false);

            rOutput[GPoint] = FluidFlux;
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                  std::vector<Vector>&    rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();

    rOutput.resize(rGeom.IntegrationPointsNumber(this->GetIntegrationMethod()));

    if (rVariable == CAUCHY_STRESS_VECTOR) {
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size() != mStressVector[GPoint].size())
                rOutput[GPoint].resize(mStressVector[GPoint].size(), false);

            rOutput[GPoint] = mStressVector[GPoint];
        }
    } else if (rVariable == TOTAL_STRESS_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, GetProperties(), rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
        std::vector<Matrix> constitutive_matrices;
        this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                             Variables.NuContainer, Variables.DNu_DXContainer,
                                             strain_vectors, mStressVector, constitutive_matrices);
        const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
            constitutive_matrices, this->GetProperties());
        const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
            Variables.NpContainer, Variables.PressureVector);
        const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] =
                mStressVector[GPoint] + PORE_PRESSURE_SIGN_FACTOR * biot_coefficients[GPoint] *
                                            bishop_coefficients[GPoint] * fluid_pressures[GPoint] *
                                            GetStressStatePolicy().GetVoigtVector();
        }
    } else if (rVariable == ENGINEERING_STRAIN_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            noalias(Variables.Nu) = row(Variables.NuContainer, GPoint);

            Matrix J0;
            Matrix InvJ0;
            this->CalculateDerivativesOnInitialConfiguration(
                Variables.detJInitialConfiguration, J0, InvJ0, Variables.DNu_DXInitialConfiguration, GPoint);

            // Calculating operator B
            Variables.B = this->CalculateBMatrix(Variables.DNu_DXInitialConfiguration, Variables.Nu);

            // Compute infinitesimal strain
            Variables.StrainVector =
                StressStrainUtilities::CalculateCauchyStrain(Variables.B, Variables.DisplacementVector);

            if (rOutput[GPoint].size() != Variables.StrainVector.size())
                rOutput[GPoint].resize(Variables.StrainVector.size(), false);

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        rOutput                          = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                  std::vector<Matrix>&    rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rOutput.resize(GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()));

    if (rVariable == CAUCHY_STRESS_TENSOR) {
        std::vector<Vector> StressVector;
        this->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
        }
    } else if (rVariable == TOTAL_STRESS_TENSOR) {
        std::vector<Vector> StressVector;
        this->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                  VectorType&        rRightHandSideVector,
                                                  const ProcessInfo& rCurrentProcessInfo,
                                                  bool               CalculateStiffnessMatrixFlag,
                                                  bool               CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const PropertiesType&                           rProp = this->GetProperties();
    const GeometryType&                             rGeom = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);

    // Stiffness matrix is needed to calculate Biot coefficient
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    RetentionLaw::Parameters RetentionParameters(rProp);

    const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJuContainer);

    const auto det_Js_initial_configuration =
        GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(rGeom, this->GetIntegrationMethod());

    const auto integration_coefficients_on_initial_configuration =
        this->CalculateIntegrationCoefficients(IntegrationPoints, det_Js_initial_configuration);

    const auto deformation_gradients = CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NuContainer, Variables.DNu_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);
    const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
        constitutive_matrices, this->GetProperties());
    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Variables.NpContainer, Variables.PressureVector);
    const auto degrees_of_saturation     = CalculateDegreesOfSaturation(fluid_pressures);
    const auto derivatives_of_saturation = CalculateDerivativesOfSaturation(fluid_pressures);
    const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
        biot_coefficients, degrees_of_saturation, derivatives_of_saturation, rProp);
    auto       relative_permeability_values = CalculateRelativePermeabilityValues(fluid_pressures);
    const auto permeability_update_factors  = GetOptionalPermeabilityUpdateFactors(strain_vectors);
    std::transform(permeability_update_factors.cbegin(), permeability_update_factors.cend(),
                   relative_permeability_values.cbegin(), relative_permeability_values.begin(),
                   std::multiplies<>{});

    const auto bishop_coefficients = CalculateBishopCoefficients(fluid_pressures);

    for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);
        Variables.B                  = b_matrices[GPoint];
        Variables.F                  = deformation_gradients[GPoint];
        Variables.StrainVector       = strain_vectors[GPoint];
        Variables.ConstitutiveMatrix = constitutive_matrices[GPoint];

        Variables.RelativePermeability = relative_permeability_values[GPoint];
        Variables.BishopCoefficient    = bishop_coefficients[GPoint];

        Variables.BiotCoefficient        = biot_coefficients[GPoint];
        Variables.BiotModulusInverse     = biot_moduli_inverse[GPoint];
        Variables.DegreeOfSaturation     = degrees_of_saturation[GPoint];
        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        Variables.IntegrationCoefficientInitialConfiguration =
            integration_coefficients_on_initial_configuration[GPoint];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

std::vector<double> SmallStrainUPwDiffOrderElement::GetOptionalPermeabilityUpdateFactors(const std::vector<Vector>& rStrainVectors) const
{
    return GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(rStrainVectors, GetProperties());
}

std::vector<double> SmallStrainUPwDiffOrderElement::CalculateDerivativesOfSaturation(const std::vector<double>& rFluidPressures)
{
    KRATOS_ERROR_IF(rFluidPressures.size() != mRetentionLawVector.size());
    std::vector<double> result;

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};
    std::transform(rFluidPressures.begin(), rFluidPressures.end(), mRetentionLawVector.begin(),
                   std::back_inserter(result), [&retention_law_params](auto fluid_pressure, auto pRetentionLaw) {
        retention_law_params.SetFluidPressure(fluid_pressure);
        return pRetentionLaw->CalculateDerivativeOfSaturation(retention_law_params);
    });

    return result;
}

std::vector<double> SmallStrainUPwDiffOrderElement::CalculateDegreesOfSaturation(const std::vector<double>& rFluidPressures)
{
    KRATOS_ERROR_IF(rFluidPressures.size() != mRetentionLawVector.size());
    std::vector<double> result;

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};
    std::transform(rFluidPressures.begin(), rFluidPressures.end(), mRetentionLawVector.begin(),
                   std::back_inserter(result), [&retention_law_params](auto fluid_pressure, auto pRetentionLaw) {
        retention_law_params.SetFluidPressure(fluid_pressure);
        return pRetentionLaw->CalculateSaturation(retention_law_params);
    });

    return result;
}

void SmallStrainUPwDiffOrderElement::CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix,
                                                                      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();

    // Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NuContainer, Variables.DNu_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJuContainer);

    const auto stiffness_matrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        b_matrices, constitutive_matrices, integration_coefficients);

    GeoElementUtilities::AssembleUUBlockMatrix(rStiffnessMatrix, stiffness_matrix);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::InitializeElementVariables(ElementVariables& rVariables,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom      = GetGeometry();
    const SizeType      NumUNodes  = rGeom.PointsNumber();
    const SizeType      NumPNodes  = mpPressureGeometry->PointsNumber();
    const SizeType      NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());
    const SizeType      Dim        = rGeom.WorkingSpaceDimension();

    // Variables at all integration points
    rVariables.NuContainer.resize(NumGPoints, NumUNodes, false);
    rVariables.NuContainer = rGeom.ShapeFunctionsValues(this->GetIntegrationMethod());

    rVariables.NpContainer.resize(NumGPoints, NumPNodes, false);
    rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues(this->GetIntegrationMethod());

    rVariables.Nu.resize(NumUNodes, false);
    rVariables.Np.resize(NumPNodes, false);

    rVariables.DNu_DXContainer.resize(NumGPoints, false);
    for (SizeType i = 0; i < NumGPoints; ++i)
        ((rVariables.DNu_DXContainer)[i]).resize(NumUNodes, Dim, false);
    rVariables.DNu_DX.resize(NumUNodes, Dim, false);
    rVariables.DNu_DXInitialConfiguration.resize(NumUNodes, Dim, false);
    rVariables.detJuContainer.resize(NumGPoints, false);

    try {
        rGeom.ShapeFunctionsIntegrationPointsGradients(
            rVariables.DNu_DXContainer, rVariables.detJuContainer, this->GetIntegrationMethod());
    } catch (Kratos::Exception& e) {
        KRATOS_INFO("Original error message") << e.what() << std::endl;
#ifdef KRATOS_COMPILED_IN_WINDOWS
        KRATOS_INFO("Error in calculation of dNu/dx. Most probably the element is "
                    "distorted. Element ID: ")
            << this->Id() << std::endl;
#endif
        KRATOS_ERROR << "In calculation of dNu/dx. Most probably the element "
                        "is distorted. Element ID: "
                     << this->Id() << std::endl;
    }

    (rVariables.DNp_DXContainer).resize(NumGPoints, false);
    for (SizeType i = 0; i < NumGPoints; ++i)
        ((rVariables.DNp_DXContainer)[i]).resize(NumPNodes, Dim, false);
    (rVariables.DNp_DX).resize(NumPNodes, Dim, false);
    Vector detJpContainer = ZeroVector(NumGPoints);

    try {
        mpPressureGeometry->ShapeFunctionsIntegrationPointsGradients(
            rVariables.DNp_DXContainer, detJpContainer, this->GetIntegrationMethod());
    } catch (Kratos::Exception& e) {
        KRATOS_INFO("Original error message") << e.what() << std::endl;
#ifdef KRATOS_COMPILED_IN_WINDOWS
        KRATOS_INFO("Error in calculation of dNp/dx. Most probably the element is "
                    "distorted. Element ID: ")
            << this->Id() << std::endl;
#endif
        KRATOS_ERROR << "In calculation of dNp/dx. Most probably the element "
                        "is distorted. Element ID: "
                     << this->Id() << std::endl;
    }

    // Variables computed at each integration point
    const SizeType VoigtSize = this->GetStressStatePolicy().GetVoigtSize();

    rVariables.B.resize(VoigtSize, NumUNodes * Dim, false);
    noalias(rVariables.B) = ZeroMatrix(VoigtSize, NumUNodes * Dim);

    rVariables.StrainVector.resize(VoigtSize, false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);

    rVariables.StressVector.resize(VoigtSize, false);

    // Needed parameters for consistency with the general constitutive law
    rVariables.F.resize(Dim, Dim, false);
    noalias(rVariables.F) = identity_matrix<double>(Dim);

    // Nodal variables
    this->InitializeNodalVariables(rVariables);

    // Properties variables
    this->InitializeProperties(rVariables);

    // ProcessInfo variables
    rVariables.VelocityCoefficient   = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    // Retention law
    rVariables.DegreeOfSaturation   = 1.0;
    rVariables.RelativePermeability = 1.0;
    rVariables.BishopCoefficient    = 1.0;

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::InitializeNodalVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    const SizeType      NumPNodes = mpPressureGeometry->PointsNumber();

    Vector BodyAccelerationAux = ZeroVector(3);
    rVariables.BodyAcceleration.resize(NumUNodes * Dim, false);
    rVariables.DisplacementVector.resize(NumUNodes * Dim, false);
    rVariables.VelocityVector.resize(NumUNodes * Dim, false);

    for (SizeType i = 0; i < NumUNodes; ++i) {
        SizeType Local_i    = i * Dim;
        BodyAccelerationAux = rGeom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        rVariables.BodyAcceleration[Local_i]   = BodyAccelerationAux[0];
        rVariables.DisplacementVector[Local_i] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
        rVariables.VelocityVector[Local_i]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);

        rVariables.BodyAcceleration[Local_i + 1] = BodyAccelerationAux[1];
        rVariables.DisplacementVector[Local_i + 1] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        rVariables.VelocityVector[Local_i + 1] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);

        if (Dim > 2) {
            rVariables.BodyAcceleration[Local_i + 2] = BodyAccelerationAux[2];
            rVariables.DisplacementVector[Local_i + 2] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
            rVariables.VelocityVector[Local_i + 2] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z);
        }
    }

    rVariables.PressureVector.resize(NumPNodes, false);
    rVariables.PressureDtVector.resize(NumPNodes, false);
    rVariables.DeltaPressureVector.resize(NumPNodes, false);
    for (SizeType i = 0; i < NumPNodes; ++i) {
        rVariables.PressureVector[i]      = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.PressureDtVector[i]    = rGeom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
        rVariables.DeltaPressureVector[i] = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE) -
                                            rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE, 1);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::InitializeProperties(ElementVariables& rVariables)
{
    KRATOS_TRY

    const SizeType        dimension = GetGeometry().WorkingSpaceDimension();
    const PropertiesType& rProp     = this->GetProperties();

    rVariables.IgnoreUndrained = rProp[IGNORE_UNDRAINED];
    rVariables.UseHenckyStrain = false;
    if (rProp.Has(USE_HENCKY_STRAIN)) rVariables.UseHenckyStrain = rProp[USE_HENCKY_STRAIN];

    rVariables.ConsiderGeometricStiffness = false;
    if (rProp.Has(CONSIDER_GEOMETRIC_STIFFNESS))
        rVariables.ConsiderGeometricStiffness = rProp[CONSIDER_GEOMETRIC_STIFFNESS];

    rVariables.DynamicViscosityInverse = 1.0 / rProp[DYNAMIC_VISCOSITY];
    // Setting the intrinsic permeability matrix
    (rVariables.IntrinsicPermeability).resize(dimension, dimension, false);
    rVariables.IntrinsicPermeability(0, 0) = rProp[PERMEABILITY_XX];
    rVariables.IntrinsicPermeability(1, 1) = rProp[PERMEABILITY_YY];
    rVariables.IntrinsicPermeability(0, 1) = rProp[PERMEABILITY_XY];
    rVariables.IntrinsicPermeability(1, 0) = rVariables.IntrinsicPermeability(0, 1);

    if (dimension == 3) {
        rVariables.IntrinsicPermeability(2, 2) = rProp[PERMEABILITY_ZZ];
        rVariables.IntrinsicPermeability(2, 0) = rProp[PERMEABILITY_ZX];
        rVariables.IntrinsicPermeability(1, 2) = rProp[PERMEABILITY_YZ];
        rVariables.IntrinsicPermeability(0, 2) = rVariables.IntrinsicPermeability(2, 0);
        rVariables.IntrinsicPermeability(2, 1) = rVariables.IntrinsicPermeability(1, 2);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateKinematics(ElementVariables& rVariables, unsigned int GPoint)

{
    KRATOS_TRY

    // Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Nu) = row(rVariables.NuContainer, GPoint);
    noalias(rVariables.Np) = row(rVariables.NpContainer, GPoint);

    noalias(rVariables.DNu_DX) = rVariables.DNu_DXContainer[GPoint];
    noalias(rVariables.DNp_DX) = rVariables.DNp_DXContainer[GPoint];

    rVariables.detJ = rVariables.detJuContainer[GPoint];

    Matrix J0;
    Matrix InvJ0;
    this->CalculateDerivativesOnInitialConfiguration(rVariables.detJInitialConfiguration, J0, InvJ0,
                                                     rVariables.DNu_DXInitialConfiguration, GPoint);

    KRATOS_CATCH("")
}

Matrix SmallStrainUPwDiffOrderElement::CalculateBMatrix(const Matrix& rDN_DX, const Vector& rN) const
{
    return this->GetStressStatePolicy().CalculateBMatrix(rDN_DX, rN, this->GetGeometry());
}

std::vector<Matrix> SmallStrainUPwDiffOrderElement::CalculateBMatrices(
    const GeometryType::ShapeFunctionsGradientsType& rDN_DXContainer, const Matrix& rNContainer) const
{
    std::vector<Matrix> result;
    for (unsigned int GPoint = 0; GPoint < rDN_DXContainer.size(); ++GPoint) {
        result.push_back(this->CalculateBMatrix(rDN_DXContainer[GPoint], row(rNContainer, GPoint)));
    }

    return result;
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddLHS(MatrixType&       rLeftHandSideMatrix,
                                                        ElementVariables& rVariables) const
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);

        const auto permeability_matrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix(
            rVariables.DNp_DX, rVariables.DynamicViscosityInverse, rVariables.IntrinsicPermeability,
            rVariables.RelativePermeability, rVariables.IntegrationCoefficient);
        GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, permeability_matrix);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix,
                                                                    const ElementVariables& rVariables) const
{
    KRATOS_TRY

    const auto stiffness_matrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrixGPoint(
        rVariables.B, rVariables.ConstitutiveMatrix, rVariables.IntegrationCoefficient);

    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, stiffness_matrix);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix,
                                                                   const ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix CouplingMatrix = GeoTransportEquationUtilities::CalculateCouplingMatrix(
        rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
        rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);
    GeoElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix, CouplingMatrix);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        Matrix CouplingMatrixT = PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient *
                                 rVariables.VelocityCoefficient * trans(CouplingMatrix);
        GeoElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix, CouplingMatrixT);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                          const ElementVariables& rVariables) const
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    GeoElementUtilities::AssemblePPBlockMatrix(
        rLeftHandSideMatrix, compressibility_matrix * rVariables.DtPressureCoefficient);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddRHS(VectorType&       rRightHandSideVector,
                                                        ElementVariables& rVariables,
                                                        unsigned int      GPoint)
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables, GPoint);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector,
                                                                   const ElementVariables& rVariables,
                                                                   unsigned int GPoint)
{
    KRATOS_TRY

    Vector stiffness_force =
        -1.0 * prod(trans(rVariables.B), mStressVector[GPoint]) * rVariables.IntegrationCoefficient;
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, stiffness_force);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector,
                                                                 ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();

    const auto soil_density = GeoTransportEquationUtilities::CalculateSoilDensity(
        rVariables.DegreeOfSaturation, this->GetProperties());

    Vector   BodyAcceleration = ZeroVector(Dim);
    SizeType Index            = 0;
    for (SizeType i = 0; i < NumUNodes; ++i) {
        for (SizeType idim = 0; idim < Dim; ++idim) {
            BodyAcceleration[idim] += rVariables.Nu[i] * rVariables.BodyAcceleration[Index++];
        }
    }

    for (SizeType i = 0; i < NumUNodes; ++i) {
        Index = i * Dim;
        for (SizeType idim = 0; idim < Dim; ++idim) {
            rRightHandSideVector[Index + idim] += rVariables.Nu[i] * soil_density * BodyAcceleration[idim] *
                                                  rVariables.IntegrationCoefficientInitialConfiguration;
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector,
                                                                  const ElementVariables& rVariables) const
{
    KRATOS_TRY

    const Matrix coupling_matrix =
        (-1.0) * GeoTransportEquationUtilities::CalculateCouplingMatrix(
                     rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
                     rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);
    const Vector coupling_force = prod(coupling_matrix, rVariables.PressureVector);
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, coupling_force);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        const Vector coupling_flow = PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient *
                                     prod(trans(coupling_matrix), rVariables.VelocityVector);
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, coupling_flow);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector,
                                                                        const ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix CompressibilityMatrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);
    Vector CompressibilityFlow = -prod(CompressibilityMatrix, rVariables.PressureDtVector);
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, CompressibilityFlow);

    KRATOS_CATCH("")
}

std::vector<double> SmallStrainUPwDiffOrderElement::CalculateRelativePermeabilityValues(const std::vector<double>& rFluidPressures) const
{
    KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

    auto result = std::vector<double>{};
    std::transform(mRetentionLawVector.begin(), mRetentionLawVector.end(), rFluidPressures.begin(),
                   std::back_inserter(result), [&retention_law_params](auto pRetentionLaw, auto FluidPressure) {
        retention_law_params.SetFluidPressure(FluidPressure);
        return pRetentionLaw->CalculateRelativePermeability(retention_law_params);
    });
    return result;
}

std::vector<double> SmallStrainUPwDiffOrderElement::CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const
{
    KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

    auto result = std::vector<double>{};
    std::transform(mRetentionLawVector.begin(), mRetentionLawVector.end(), rFluidPressures.begin(),
                   std::back_inserter(result), [&retention_law_params](auto pRetentionLaw, auto FluidPressure) {
        retention_law_params.SetFluidPressure(FluidPressure);
        return pRetentionLaw->CalculateBishopCoefficient(retention_law_params);
    });
    return result;
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                                                     const ElementVariables& rVariables) const
{
    KRATOS_TRY

    const Matrix permeability_matrix =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse * rVariables.RelativePermeability *
        prod(rVariables.DNp_DX, Matrix(prod(rVariables.IntrinsicPermeability, trans(rVariables.DNp_DX)))) *
        rVariables.IntegrationCoefficient;

    const Vector permeability_flow = -prod(permeability_matrix, rVariables.PressureVector);

    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, permeability_flow);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                  const ElementVariables& rVariables)
{
    KRATOS_TRY

    const Matrix grad_Np_T_perm =
        rVariables.DynamicViscosityInverse * GetProperties()[DENSITY_WATER] * rVariables.RelativePermeability *
        prod(rVariables.DNp_DX, rVariables.IntrinsicPermeability) * rVariables.IntegrationCoefficient;

    const GeometryType& r_geom      = GetGeometry();
    const SizeType      dimension   = r_geom.WorkingSpaceDimension();
    const SizeType      num_U_nodes = r_geom.PointsNumber();
    const SizeType      num_P_nodes = mpPressureGeometry->PointsNumber();

    Vector body_acceleration = ZeroVector(dimension);

    SizeType index = 0;
    for (SizeType i = 0; i < num_U_nodes; ++i) {
        for (SizeType idim = 0; idim < dimension; ++idim) {
            body_acceleration[idim] += rVariables.Nu[i] * rVariables.BodyAcceleration[index];
            index++;
        }
    }

    for (SizeType i = 0; i < num_P_nodes; ++i) {
        rRightHandSideVector[num_U_nodes * dimension + i] +=
            inner_prod(row(grad_Np_T_perm, i), body_acceleration);
    }

    KRATOS_CATCH("")
}

Vector SmallStrainUPwDiffOrderElement::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
{
    return this->GetStressStatePolicy().CalculateGreenLagrangeStrain(rDeformationGradient);
}

Matrix SmallStrainUPwDiffOrderElement::CalculateDeformationGradient(unsigned int GPoint) const
{
    KRATOS_TRY

    // Calculation of derivative of shape function with respect to reference
    // configuration derivative of shape function (displacement)
    Matrix J0;
    Matrix InvJ0;
    Matrix DNu_DX0;
    double detJ0;
    this->CalculateDerivativesOnInitialConfiguration(detJ0, J0, InvJ0, DNu_DX0, GPoint);

    // Calculating current Jacobian in order to find deformation gradient
    Matrix J;
    Matrix InvJ;
    double detJ;
    this->CalculateJacobianOnCurrentConfiguration(detJ, J, InvJ, GPoint);

    KRATOS_ERROR_IF(detJ < 0.0)
        << "ERROR:: Element " << this->Id() << " is inverted. DetJ: " << detJ << std::endl
        << "This usually indicates that the deformations are too large for the mesh size." << std::endl;

    return prod(J, InvJ0);

    KRATOS_CATCH("")
}

std::vector<Matrix> SmallStrainUPwDiffOrderElement::CalculateDeformationGradients() const
{
    std::vector<Matrix> result;
    for (unsigned int GPoint = 0;
         GPoint < this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()); ++GPoint) {
        result.push_back(CalculateDeformationGradient(GPoint));
    }

    return result;
}

SizeType SmallStrainUPwDiffOrderElement::GetNumberOfDOF() const
{
    return GetGeometry().PointsNumber() * GetGeometry().WorkingSpaceDimension() +
           mpPressureGeometry->PointsNumber();
}

Element::DofsVectorType SmallStrainUPwDiffOrderElement::GetDofs() const
{
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), *mpPressureGeometry,
                                                      GetGeometry().WorkingSpaceDimension());
}

void SmallStrainUPwDiffOrderElement::SetUpPressureGeometryPointer()
{
    const auto& r_geometry        = GetGeometry();
    const auto  number_of_U_nodes = r_geometry.PointsNumber();
    const auto  dimension         = r_geometry.WorkingSpaceDimension();
    switch (number_of_U_nodes) {
    case 6: // 2D T6P3
        mpPressureGeometry = make_shared<Triangle2D3<Node>>(r_geometry(0), r_geometry(1), r_geometry(2));
        break;
    case 8: // 2D Q8P4
        mpPressureGeometry = make_shared<Quadrilateral2D4<Node>>(r_geometry(0), r_geometry(1),
                                                                 r_geometry(2), r_geometry(3));
        break;
    case 9: // 2D Q9P4
        mpPressureGeometry = make_shared<Quadrilateral2D4<Node>>(r_geometry(0), r_geometry(1),
                                                                 r_geometry(2), r_geometry(3));
        break;
    case 10:
        if (dimension == 3) // 3D T10P4
            mpPressureGeometry = make_shared<Tetrahedra3D4<Node>>(r_geometry(0), r_geometry(1),
                                                                  r_geometry(2), r_geometry(3));
        else if (dimension == 2) // 2D T10P6
            mpPressureGeometry = make_shared<Triangle2D6<Node>>(
                r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3), r_geometry(4), r_geometry(5));
        break;
    case 15: // 2D T15P10
        mpPressureGeometry = make_shared<Triangle2D10<Node>>(
            r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3), r_geometry(4),
            r_geometry(5), r_geometry(6), r_geometry(7), r_geometry(8), r_geometry(9));
        break;
    case 20: // 3D H20P8
        mpPressureGeometry =
            make_shared<Hexahedra3D8<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3),
                                            r_geometry(4), r_geometry(5), r_geometry(6), r_geometry(7));
        break;
    case 27: // 3D H27P8
        mpPressureGeometry =
            make_shared<Hexahedra3D8<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3),
                                            r_geometry(4), r_geometry(5), r_geometry(6), r_geometry(7));
        break;
    default:
        KRATOS_ERROR << "Unexpected geometry type for different order interpolation element "
                     << this->Id() << std::endl;
    }
}

Vector SmallStrainUPwDiffOrderElement::GetPressureSolutionVector()
{
    Vector result(mpPressureGeometry->PointsNumber());
    std::transform(this->GetGeometry().begin(),
                   this->GetGeometry().begin() + mpPressureGeometry->PointsNumber(), result.begin(),
                   [](const auto& node) { return node.FastGetSolutionStepValue(WATER_PRESSURE); });
    return result;
}

void SmallStrainUPwDiffOrderElement::CalculateAnyOfMaterialResponse(
    const std::vector<Matrix>&                       rDeformationGradients,
    ConstitutiveLaw::Parameters&                     rConstitutiveParameters,
    const Matrix&                                    rNuContainer,
    const GeometryType::ShapeFunctionsGradientsType& rDNu_DXContainer,
    std::vector<Vector>&                             rStrainVectors,
    std::vector<Vector>&                             rStressVectors,
    std::vector<Matrix>&                             rConstitutiveMatrices)
{
    const SizeType voigt_size =
        GetGeometry().WorkingSpaceDimension() == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN;

    if (rStrainVectors.size() != rDeformationGradients.size()) {
        rStrainVectors.resize(rDeformationGradients.size());
        std::fill(rStrainVectors.begin(), rStrainVectors.end(), ZeroVector(voigt_size));
    }
    if (rStressVectors.size() != rDeformationGradients.size()) {
        rStressVectors.resize(rDeformationGradients.size());
        std::fill(rStressVectors.begin(), rStressVectors.end(), ZeroVector(voigt_size));
    }
    if (rConstitutiveMatrices.size() != rDeformationGradients.size()) {
        rConstitutiveMatrices.resize(rDeformationGradients.size());
        std::fill(rConstitutiveMatrices.begin(), rConstitutiveMatrices.end(), ZeroMatrix(voigt_size, voigt_size));
    }

    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(rDeformationGradients);

    for (unsigned int GPoint = 0; GPoint < rDeformationGradients.size(); ++GPoint) {
        ConstitutiveLawUtilities::SetConstitutiveParameters(
            rConstitutiveParameters, rStrainVectors[GPoint], rConstitutiveMatrices[GPoint],
            row(rNuContainer, GPoint), rDNu_DXContainer[GPoint], rDeformationGradients[GPoint],
            determinants_of_deformation_gradients[GPoint]);
        rConstitutiveParameters.SetStressVector(rStressVectors[GPoint]);

        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(rConstitutiveParameters);
    }
}

} // Namespace Kratos
