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
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_3.h"

// Application includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

// Default Constructor
SmallStrainUPwDiffOrderElement::SmallStrainUPwDiffOrderElement() : Element() {}

// Constructor 1
SmallStrainUPwDiffOrderElement::SmallStrainUPwDiffOrderElement(IndexType             NewId,
                                                               GeometryType::Pointer pGeometry,
                                                               std::unique_ptr<StressStatePolicy> pStressStatePolicy)
    : Element(NewId, pGeometry), mpStressStatePolicy(std::move(pStressStatePolicy))
{
}

// Constructor 2
SmallStrainUPwDiffOrderElement::SmallStrainUPwDiffOrderElement(IndexType               NewId,
                                                               GeometryType::Pointer   pGeometry,
                                                               PropertiesType::Pointer pProperties,
                                                               std::unique_ptr<StressStatePolicy> pStressStatePolicy)
    : Element(NewId, pGeometry, pProperties), mpStressStatePolicy(std::move(pStressStatePolicy))
{
}

SmallStrainUPwDiffOrderElement::~SmallStrainUPwDiffOrderElement() = default;

Element::Pointer SmallStrainUPwDiffOrderElement::Create(IndexType               NewId,
                                                        NodesArrayType const&   ThisNodes,
                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallStrainUPwDiffOrderElement(
        NewId, GetGeometry().Create(ThisNodes), pProperties, mpStressStatePolicy->Clone()));
}

Element::Pointer SmallStrainUPwDiffOrderElement::Create(IndexType               NewId,
                                                        GeometryType::Pointer   pGeom,
                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
        new SmallStrainUPwDiffOrderElement(NewId, pGeom, pProperties, mpStressStatePolicy->Clone()));
}

int SmallStrainUPwDiffOrderElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();

    if (rGeom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

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

void SmallStrainUPwDiffOrderElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType&                             rGeom = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    if (mConstitutiveLawVector.size() != IntegrationPoints.size())
        mConstitutiveLawVector.resize(IntegrationPoints.size());

    if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial(
                GetProperties(), rGeom, row(rGeom.ShapeFunctionsValues(this->GetIntegrationMethod()), i));
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the "
                        "element with ID "
                     << this->Id() << std::endl;

    // Retention law
    if (mRetentionLawVector.size() != IntegrationPoints.size())
        mRetentionLawVector.resize(IntegrationPoints.size());

    for (unsigned int i = 0; i < mRetentionLawVector.size(); ++i) {
        mRetentionLawVector[i] = RetentionLawFactory::Clone(GetProperties());
        mRetentionLawVector[i]->InitializeMaterial(
            GetProperties(), rGeom, row(rGeom.ShapeFunctionsValues(this->GetIntegrationMethod()), i));
    }

    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumDim    = rGeom.WorkingSpaceDimension();

    switch (NumUNodes) {
    case 6: // 2D T6P3
        mpPressureGeometry = make_shared<Triangle2D3<Node>>(rGeom(0), rGeom(1), rGeom(2));
        break;
    case 8: // 2D Q8P4
        mpPressureGeometry = make_shared<Quadrilateral2D4<Node>>(rGeom(0), rGeom(1), rGeom(2), rGeom(3));
        break;
    case 9: // 2D Q9P4
        mpPressureGeometry = make_shared<Quadrilateral2D4<Node>>(rGeom(0), rGeom(1), rGeom(2), rGeom(3));
        break;
    case 10: // 3D T10P4  //2D T10P6
        if (NumDim == 3)
            mpPressureGeometry = make_shared<Tetrahedra3D4<Node>>(rGeom(0), rGeom(1), rGeom(2), rGeom(3));
        else if (NumDim == 2)
            mpPressureGeometry = make_shared<Triangle2D6<Node>>(rGeom(0), rGeom(1), rGeom(2),
                                                                rGeom(3), rGeom(4), rGeom(5));
        break;
    case 15: // 2D T15P10
        mpPressureGeometry =
            make_shared<Triangle2D10<Node>>(rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4),
                                            rGeom(5), rGeom(6), rGeom(7), rGeom(8), rGeom(9));
        break;
    case 20: // 3D H20P8
        mpPressureGeometry = make_shared<Hexahedra3D8<Node>>(
            rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7));
        break;
    case 27: // 3D H27P8
        mpPressureGeometry = make_shared<Hexahedra3D8<Node>>(
            rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7));
        break;
    default:
        KRATOS_ERROR << "Unexpected geometry type for different order "
                        "interpolation element"
                     << this->Id() << std::endl;
    }

    // resize mStressVector:
    const SizeType Dim       = rGeom.WorkingSpaceDimension();
    const SizeType VoigtSize = (Dim == N_DIM_3D ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);
    if (mStressVector.size() != IntegrationPoints.size()) {
        mStressVector.resize(IntegrationPoints.size());
        for (unsigned int i = 0; i < mStressVector.size(); ++i) {
            mStressVector[i].resize(VoigtSize);
            std::fill(mStressVector[i].begin(), mStressVector[i].end(), 0.0);
        }
    }

    if (mStateVariablesFinalized.size() != IntegrationPoints.size())
        mStateVariablesFinalized.resize(IntegrationPoints.size());

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        int nStateVariables = 0;
        nStateVariables = mConstitutiveLawVector[i]->GetValue(NUMBER_OF_UMAT_STATE_VARIABLES, nStateVariables);
        if (nStateVariables > 0) {
            mConstitutiveLawVector[i]->SetValue(STATE_VARIABLES, mStateVariablesFinalized[i], rCurrentProcessInfo);
        }
    }

    mIsInitialised = true;

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    // erasing stress vectors
    for (unsigned int i = 0; i < mStressVector.size(); ++i) {
        mStressVector[i].clear();
    }
    mStressVector.clear();

    for (unsigned int i = 0; i < mStateVariablesFinalized.size(); ++i) {
        mStateVariablesFinalized[i].clear();
    }
    mStateVariablesFinalized.clear();

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);

    // Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // Set constitutive law flags:
    // ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); // Note: this is for nonlocal damage

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(GetProperties(), rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
        // compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // set gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->InitializeMaterialResponseCauchy(ConstitutiveParameters);

        // retention law
        mRetentionLawVector[GPoint]->InitializeSolutionStep(RetentionParameters);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

void SmallStrainUPwDiffOrderElement::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                                                          VectorType&        rRightHandSideVector,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    const SizeType      NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    // Resetting the LHS
    if (rLeftHandSideMatrix.size1() != ElementSize)
        rLeftHandSideMatrix.resize(ElementSize, ElementSize, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(ElementSize, ElementSize);

    // Resetting the RHS
    if (rRightHandSideVector.size() != ElementSize) rRightHandSideVector.resize(ElementSize, false);
    noalias(rRightHandSideVector) = ZeroVector(ElementSize);

    // calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag  = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateLeftHandSide(MatrixType&        rLeftHandSideMatrix,
                                                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    const SizeType      NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    // Resetting the LHS
    if (rLeftHandSideMatrix.size1() != ElementSize)
        rLeftHandSideMatrix.resize(ElementSize, ElementSize, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(ElementSize, ElementSize);

    // calculation flags
    bool       CalculateStiffnessMatrixFlag = true;
    bool       CalculateResidualVectorFlag  = false;
    VectorType tempRightHandSideVector;

    CalculateAll(rLeftHandSideMatrix, tempRightHandSideVector, rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateRightHandSide(VectorType&        rRightHandSideVector,
                                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    const SizeType      NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    // Resetting the RHS
    if (rRightHandSideVector.size() != ElementSize) rRightHandSideVector.resize(ElementSize, false);
    noalias(rRightHandSideVector) = ZeroVector(ElementSize);

    // calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag  = true;
    auto temp                         = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType&                             rGeom     = GetGeometry();
    const SizeType                                  Dim       = rGeom.WorkingSpaceDimension();
    const SizeType                                  NumUNodes = rGeom.PointsNumber();
    const SizeType                                  BlockElementSize = NumUNodes * Dim;
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    Matrix MassMatrixContribution = ZeroMatrix(BlockElementSize, BlockElementSize);

    // Defining shape functions and the determinant of the jacobian at all integration points

    // Loop over integration points
    Matrix Nu               = ZeroMatrix(Dim, NumUNodes * Dim);
    Matrix AuxDensityMatrix = ZeroMatrix(Dim, NumUNodes * Dim);
    Matrix DensityMatrix    = ZeroMatrix(Dim, Dim);

    for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        // compute element kinematics (Np, gradNpT, |J|, B)
        this->CalculateKinematics(Variables, GPoint);

        // calculating weighting coefficient for integration
        Variables.IntegrationCoefficientInitialConfiguration = this->CalculateIntegrationCoefficient(
            IntegrationPoints, GPoint, Variables.detJInitialConfiguration);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->CalculateSoilDensity(Variables);

        // Setting the shape function matrix
        SizeType Index = 0;
        for (SizeType i = 0; i < NumUNodes; ++i) {
            for (SizeType iDim = 0; iDim < Dim; ++iDim) {
                Nu(iDim, Index++) = Variables.Nu(i);
            }
        }

        GeoElementUtilities::AssembleDensityMatrix(DensityMatrix, Variables.Density);

        noalias(AuxDensityMatrix) = prod(DensityMatrix, Nu);

        // Adding contribution to Mass matrix
        noalias(MassMatrixContribution) +=
            prod(trans(Nu), AuxDensityMatrix) * Variables.IntegrationCoefficientInitialConfiguration;
    }

    // Distribute mass block matrix into the elemental matrix
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType ElementSize = BlockElementSize + NumPNodes;

    if (rMassMatrix.size1() != ElementSize || rMassMatrix.size2() != ElementSize)
        rMassMatrix.resize(ElementSize, ElementSize, false);
    noalias(rMassMatrix) = ZeroMatrix(ElementSize, ElementSize);

    for (SizeType i = 0; i < NumUNodes; ++i) {
        SizeType Index_i = i * Dim;

        for (SizeType j = 0; j < NumUNodes; ++j) {
            SizeType Index_j = j * Dim;
            for (SizeType idim = 0; idim < Dim; ++idim) {
                for (SizeType jdim = 0; jdim < Dim; ++jdim) {
                    rMassMatrix(Index_i + idim, Index_j + jdim) +=
                        MassMatrixContribution(Index_i + idim, Index_j + jdim);
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateDampingMatrix(MatrixType&        rDampingMatrix,
                                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Rayleigh Method: Damping Matrix = alpha*M + beta*K
    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    const SizeType      NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    // Compute Mass Matrix
    MatrixType MassMatrix(ElementSize, ElementSize);

    this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

    // Compute Stiffness matrix
    MatrixType StiffnessMatrix(ElementSize, ElementSize);

    this->CalculateMaterialStiffnessMatrix(StiffnessMatrix, rCurrentProcessInfo);

    // Compute Damping Matrix
    if (rDampingMatrix.size1() != ElementSize)
        rDampingMatrix.resize(ElementSize, ElementSize, false);
    noalias(rDampingMatrix) = ZeroMatrix(ElementSize, ElementSize);

    const PropertiesType& rProp = this->GetProperties();

    if (rProp.Has(RAYLEIGH_ALPHA)) noalias(rDampingMatrix) += rProp[RAYLEIGH_ALPHA] * MassMatrix;
    else noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_ALPHA] * MassMatrix;

    if (rProp.Has(RAYLEIGH_BETA)) noalias(rDampingMatrix) += rProp[RAYLEIGH_BETA] * StiffnessMatrix;
    else noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_BETA] * StiffnessMatrix;

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void SmallStrainUPwDiffOrderElement::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractFirstTimeDerivativesOfUPwDofs(GetDofs(), Step);
}

void SmallStrainUPwDiffOrderElement::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSecondTimeDerivativesOfUPwDofs(GetDofs(), Step);
}

void SmallStrainUPwDiffOrderElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(GetProperties(), rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
        // compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // set gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
        mStateVariablesFinalized[GPoint] =
            mConstitutiveLawVector[GPoint]->GetValue(STATE_VARIABLES, mStateVariablesFinalized[GPoint]);

        // retention law
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
        ThreadSafeNodeWrite(rGeom[19], WATER_PRESSURE, 0.5 * (p7 + p0));
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

void SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                                                  const std::vector<double>& rValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
        mConstitutiveLawVector[GPoint]->SetValue(rVariable, rValues[GPoint], rCurrentProcessInfo);
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

void SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                  const std::vector<Matrix>& rValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
        mConstitutiveLawVector[GPoint]->SetValue(rVariable, rValues[GPoint], rCurrentProcessInfo);
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

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                                                  std::vector<ConstitutiveLaw::Pointer>& rValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == CONSTITUTIVE_LAW) {
        if (rValues.size() != mConstitutiveLawVector.size()) {
            rValues.resize(mConstitutiveLawVector.size());
        }
        for (unsigned int i = 0; i < rValues.size(); ++i) {
            rValues[i] = mConstitutiveLawVector[i];
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                  std::vector<double>&    rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType& IntegrationPointsNumber = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

    if (rOutput.size() != IntegrationPointsNumber) rOutput.resize(IntegrationPointsNumber, false);

    if (rVariable == VON_MISES_STRESS) {
        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStress(mStressVector[GPoint]);
        }
    } else if (rVariable == MEAN_EFFECTIVE_STRESS) {
        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateMeanStress(mStressVector[GPoint]);
        }
    } else if (rVariable == MEAN_STRESS) {
        std::vector<Vector> StressVector;
        CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateMeanStress(StressVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_VON_MISES_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStrain(StrainVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_VOLUMETRIC_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateTrace(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VON_MISES_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStrain(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VOLUMETRIC_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateTrace(StrainVector[GPoint]);
        }
    } else if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION ||
               rVariable == BISHOP_COEFFICIENT || rVariable == DERIVATIVE_OF_SATURATION ||
               rVariable == RELATIVE_PERMEABILITY) {
        // Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(GetProperties(), rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mRetentionLawVector.size(); ++GPoint) {
            // Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            Variables.FluidPressure = CalculateFluidPressure(Variables);
            SetRetentionParameters(Variables, RetentionParameters);

            if (rVariable == DEGREE_OF_SATURATION)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateSaturation(RetentionParameters);
            if (rVariable == EFFECTIVE_SATURATION)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateEffectiveSaturation(RetentionParameters);
            if (rVariable == BISHOP_COEFFICIENT)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(RetentionParameters);
            if (rVariable == DERIVATIVE_OF_SATURATION)
                rOutput[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(RetentionParameters);
            if (rVariable == RELATIVE_PERMEABILITY)
                rOutput[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);
        }
    } else if (rVariable == HYDRAULIC_HEAD) {
        const double          NumericalLimit = std::numeric_limits<double>::epsilon();
        const PropertiesType& rProp          = this->GetProperties();
        const SizeType NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

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

        if (rOutput.size() != NumGPoints) rOutput.resize(NumGPoints);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++) {
            double HydraulicHead = 0.0;
            for (unsigned int node = 0; node < NumUNodes; ++node)
                HydraulicHead += NContainer(GPoint, node) * NodalHydraulicHead[node];

            rOutput[GPoint] = HydraulicHead;
        }
    } else {
        for (unsigned int i = 0; i < IntegrationPointsNumber; ++i)
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                  std::vector<array_1d<double, 3>>& rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType& IntegrationPointsNumber = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

    if (rOutput.size() != IntegrationPointsNumber) rOutput.resize(IntegrationPointsNumber);

    if (rVariable == FLUID_FLUX_VECTOR) {
        // Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(GetProperties(), rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            // compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables, GPoint);

            // Compute FluidFlux vector q [m/s]
            const SizeType Dim       = rGeom.WorkingSpaceDimension();
            const SizeType NumUNodes = rGeom.PointsNumber();

            Vector   BodyAcceleration = ZeroVector(Dim);
            SizeType Index            = 0;
            for (SizeType i = 0; i < NumUNodes; ++i) {
                for (unsigned int idim = 0; idim < Dim; ++idim)
                    BodyAcceleration[idim] += Variables.Nu[i] * Variables.BodyAcceleration[Index++];
            }

            CalculateFluidPressure(Variables);
            SetRetentionParameters(Variables, RetentionParameters);

            const double RelativePermeability =
                mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);

            // Compute strain, need to update porosity
            this->CalculateStrain(Variables, GPoint);
            this->CalculatePermeabilityUpdateFactor(Variables);

            Vector GradPressureTerm(Dim);
            noalias(GradPressureTerm) = prod(trans(Variables.DNp_DX), Variables.PressureVector);
            noalias(GradPressureTerm) +=
                PORE_PRESSURE_SIGN_FACTOR * GetProperties()[DENSITY_WATER] * BodyAcceleration;

            Vector AuxFluidFlux = ZeroVector(Dim);
            AuxFluidFlux        = PORE_PRESSURE_SIGN_FACTOR * Variables.DynamicViscosityInverse *
                           RelativePermeability * Variables.PermeabilityUpdateFactor *
                           prod(Variables.IntrinsicPermeability, GradPressureTerm);

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
    const SizeType& IntegrationPointsNumber = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());
    const SizeType Dim = rGeom.WorkingSpaceDimension();

    if (rOutput.size() != IntegrationPointsNumber) rOutput.resize(IntegrationPointsNumber);

    if (rVariable == CAUCHY_STRESS_VECTOR) {
        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size() != mStressVector[GPoint].size())
                rOutput[GPoint].resize(mStressVector[GPoint].size(), false);

            rOutput[GPoint] = mStressVector[GPoint];
        }
    } else if (rVariable == ENGINEERING_STRAIN_VECTOR) {
        // Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            noalias(Variables.Nu) = row(Variables.NuContainer, GPoint);

            Matrix J0;
            Matrix InvJ0;
            this->CalculateDerivativesOnInitialConfiguration(
                Variables.detJInitialConfiguration, J0, InvJ0, Variables.DNu_DXInitialConfiguration, GPoint);

            // Calculating operator B
            this->CalculateBMatrix(Variables.B, Variables.DNu_DXInitialConfiguration, Variables.Nu);

            // Compute infinitesimal strain
            this->CalculateCauchyStrain(Variables);

            if (rOutput[GPoint].size() != Variables.StrainVector.size())
                rOutput[GPoint].resize(Variables.StrainVector.size(), false);

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        // Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            // compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables, GPoint);

            // Compute infinitesimal strain
            this->CalculateStrain(Variables, GPoint);

            if (rOutput[GPoint].size() != Variables.StrainVector.size())
                rOutput[GPoint].resize(Variables.StrainVector.size(), false);

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else if (rVariable == TOTAL_STRESS_VECTOR) {
        // Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, GetProperties(), rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        const PropertiesType& rProp = this->GetProperties();

        const SizeType VoigtSize   = mStressVector[0].size();
        Vector         VoigtVector = ZeroVector(VoigtSize);
        const SizeType StressTensorSize = (Dim == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D);
        for (unsigned int i = 0; i < StressTensorSize; ++i)
            VoigtVector[i] = 1.0;

        // create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(GetProperties(), rCurrentProcessInfo);

        const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

        Vector TotalStressVector(mStressVector[0].size());

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            // compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables, GPoint);

            // Compute infinitesimal strain
            this->CalculateStrain(Variables, GPoint);

            // set gauss points variables to constitutive law parameters
            this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

            // compute constitutive tensor and/or stresses
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            Variables.BiotCoefficient = CalculateBiotCoefficient(Variables, hasBiotCoefficient);

            this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

            noalias(TotalStressVector) = mStressVector[GPoint];
            noalias(TotalStressVector) += PORE_PRESSURE_SIGN_FACTOR * Variables.BiotCoefficient *
                                          Variables.BishopCoefficient * Variables.FluidPressure * VoigtVector;

            if (rOutput[GPoint].size() != TotalStressVector.size()) {
                rOutput[GPoint].resize(TotalStressVector.size(), false);
            }
            rOutput[GPoint] = TotalStressVector;
        }
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

    const GeometryType& rGeom = GetGeometry();
    const auto number_of_integration_points = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());
    const auto dimension = rGeom.WorkingSpaceDimension();

    rOutput.resize(number_of_integration_points);

    if (rVariable == CAUCHY_STRESS_TENSOR) {
        std::vector<Vector> StressVector;

        this->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size2() != dimension)
                rOutput[GPoint].resize(dimension, dimension, false);

            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size2() != dimension)
                rOutput[GPoint].resize(dimension, dimension, false);

            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size2() != dimension)
                rOutput[GPoint].resize(dimension, dimension, false);

            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else if (rVariable == TOTAL_STRESS_TENSOR) {
        std::vector<Vector> StressVector;

        this->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size2() != dimension)
                rOutput[GPoint].resize(dimension, dimension, false);

            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
        }
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
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

    const GeometryType&   rGeom = GetGeometry();
    const PropertiesType& rProp = this->GetProperties();

    // Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Stiffness matrix is always needed to calculate Biot coefficient
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(rProp, rCurrentProcessInfo);

    // Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

    for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        // compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // set gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // compute constitutive tensor and/or stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);
        this->CalculatePermeabilityUpdateFactor(Variables);

        // calculating weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

        Variables.IntegrationCoefficientInitialConfiguration = this->CalculateIntegrationCoefficient(
            IntegrationPoints, GPoint, Variables.detJInitialConfiguration);

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
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

    for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        // compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // set gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        // calculating weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

        // Contributions of material stiffness to the left hand side
        this->CalculateAndAddStiffnessMatrix(rStiffnessMatrix, Variables);
    }

    KRATOS_CATCH("")
}

double SmallStrainUPwDiffOrderElement::CalculateBulkModulus(const Matrix& ConstitutiveMatrix) const
{
    KRATOS_TRY

    const SizeType IndexG = ConstitutiveMatrix.size1() - 1;
    return ConstitutiveMatrix(0, 0) - (4.0 / 3.0) * ConstitutiveMatrix(IndexG, IndexG);

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
    const SizeType VoigtSize = (Dim == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);

    rVariables.B.resize(VoigtSize, NumUNodes * Dim, false);
    noalias(rVariables.B) = ZeroMatrix(VoigtSize, NumUNodes * Dim);

    rVariables.StrainVector.resize(VoigtSize, false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);

    rVariables.StressVector.resize(VoigtSize, false);

    // Needed parameters for consistency with the general constitutive law
    rVariables.detF = 1.0;
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
    rVariables.FluidPressure          = 0.0;
    rVariables.DegreeOfSaturation     = 1.0;
    rVariables.DerivativeOfSaturation = 0.0;
    rVariables.RelativePermeability   = 1.0;
    rVariables.BishopCoefficient      = 1.0;

    // permeability change
    rVariables.PermeabilityUpdateFactor = 1.0;

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

double SmallStrainUPwDiffOrderElement::CalculateBiotCoefficient(const ElementVariables& rVariables,
                                                                const bool& hasBiotCoefficient) const
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    // Properties variables
    if (hasBiotCoefficient) {
        return rProp[BIOT_COEFFICIENT];
    } else {
        // calculate Bulk modulus from stiffness matrix
        return 1.0 - CalculateBulkModulus(rVariables.ConstitutiveMatrix) / rProp[BULK_MODULUS_SOLID];
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::InitializeBiotCoefficients(ElementVariables& rVariables,
                                                                const bool& hasBiotCoefficient)
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    rVariables.BiotCoefficient = CalculateBiotCoefficient(rVariables, hasBiotCoefficient);

    if (rProp[IGNORE_UNDRAINED]) {
        rVariables.BiotModulusInverse =
            (rVariables.BiotCoefficient - rProp[POROSITY]) / rProp[BULK_MODULUS_SOLID] + rProp[POROSITY] / TINY;
    } else {
        rVariables.BiotModulusInverse =
            (rVariables.BiotCoefficient - rProp[POROSITY]) / rProp[BULK_MODULUS_SOLID] +
            rProp[POROSITY] / rProp[BULK_MODULUS_FLUID];
    }

    rVariables.BiotModulusInverse *= rVariables.DegreeOfSaturation;
    rVariables.BiotModulusInverse -= rVariables.DerivativeOfSaturation * rProp[POROSITY];

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculatePermeabilityUpdateFactor(ElementVariables& rVariables)
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    if (rProp[PERMEABILITY_CHANGE_INVERSE_FACTOR] > 0.0) {
        const double InverseCK = rProp[PERMEABILITY_CHANGE_INVERSE_FACTOR];
        const double epsV      = StressStrainUtilities::CalculateTrace(rVariables.StrainVector);
        const double ePrevious = rProp[POROSITY] / (1.0 - rProp[POROSITY]);
        const double eCurrent  = (1.0 + ePrevious) * std::exp(epsV) - 1.0;
        const double permLog10 = (eCurrent - ePrevious) * InverseCK;
        rVariables.PermeabilityUpdateFactor = pow(10.0, permLog10);
    } else {
        rVariables.PermeabilityUpdateFactor = 1.0;
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

    // Compute the deformation matrix B
    this->CalculateBMatrix(rVariables.B, rVariables.DNu_DX, rVariables.Nu);

    rVariables.detJ = rVariables.detJuContainer[GPoint];

    Matrix J0;
    Matrix InvJ0;
    this->CalculateDerivativesOnInitialConfiguration(rVariables.detJInitialConfiguration, J0, InvJ0,
                                                     rVariables.DNu_DXInitialConfiguration, GPoint);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateDerivativesOnInitialConfiguration(
    double& detJ, Matrix& J0, Matrix& InvJ0, Matrix& DNu_DX0, unsigned int GPoint) const
{
    KRATOS_TRY

    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    GeometryUtils::JacobianOnInitialConfiguration(rGeom, IntegrationPoints[GPoint], J0);
    const Matrix& DN_De = rGeom.ShapeFunctionsLocalGradients(this->GetIntegrationMethod())[GPoint];
    MathUtils<double>::InvertMatrix(J0, InvJ0, detJ);
    GeometryUtils::ShapeFunctionsGradients(DN_De, InvJ0, DNu_DX0);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateBMatrix(Matrix& rB, const Matrix& DNp_DX, const Vector& Np)
{
    rB = mpStressStatePolicy->CalculateBMatrix(DNp_DX, Np, this->GetGeometry());
}

void SmallStrainUPwDiffOrderElement::SetConstitutiveParameters(ElementVariables& rVariables,
                                                               ConstitutiveLaw::Parameters& rConstitutiveParameters) const
{
    KRATOS_TRY

    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);

    // Needed parameters for consistency with the general constitutive law
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.DNu_DX);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Nu);

    rConstitutiveParameters.SetDeterminantF(rVariables.detF);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);

    KRATOS_CATCH("")
}

double SmallStrainUPwDiffOrderElement::CalculateIntegrationCoefficient(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints, unsigned int GPoint, double detJ)
{
    return mpStressStatePolicy->CalculateIntegrationCoefficient(IntegrationPoints[GPoint], detJ, GetGeometry());
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);
        this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix,
                                                                    ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix StiffnessMatrix =
        prod(trans(rVariables.B), Matrix(prod(rVariables.ConstitutiveMatrix, rVariables.B))) *
        rVariables.IntegrationCoefficient;

    // Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, StiffnessMatrix);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix,
                                                                   ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      VoigtSize = (Dim == N_DIM_3D ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);
    const SizeType StressTensorSize = (Dim == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D);

    Vector VoigtVector = ZeroVector(VoigtSize);

    for (unsigned int i = 0; i < StressTensorSize; ++i)
        VoigtVector[i] = 1.0;

    Matrix CouplingMatrix = PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient * rVariables.BishopCoefficient *
                            prod(trans(rVariables.B), Matrix(outer_prod(VoigtVector, rVariables.Np))) *
                            rVariables.IntegrationCoefficient;

    // Distribute coupling block matrix into the elemental matrix
    GeoElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix, CouplingMatrix);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        Matrix CouplingMatrixT = PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient *
                                 rVariables.VelocityCoefficient * trans(CouplingMatrix);

        // Distribute transposed coupling block matrix into the elemental matrix
        GeoElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix, CouplingMatrixT);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                          ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix CompressibilityMatrix =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.DtPressureCoefficient * rVariables.BiotModulusInverse *
        outer_prod(rVariables.Np, rVariables.Np) * rVariables.IntegrationCoefficient;

    // Distribute compressibility block matrix into the elemental matrix
    GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, CompressibilityMatrix);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                       ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix PermeabilityMatrix =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse *
        rVariables.RelativePermeability * rVariables.PermeabilityUpdateFactor *
        prod(rVariables.DNp_DX, Matrix(prod(rVariables.IntrinsicPermeability, trans(rVariables.DNp_DX)))) *
        rVariables.IntegrationCoefficient;

    // Distribute permeability block matrix into the elemental matrix
    GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, PermeabilityMatrix);

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
                                                                   ElementVariables& rVariables,
                                                                   unsigned int      GPoint)
{
    KRATOS_TRY

    Vector StiffnessForce =
        -1.0 * prod(trans(rVariables.B), mStressVector[GPoint]) * rVariables.IntegrationCoefficient;

    // Distribute stiffness block vector into the elemental vector
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, StiffnessForce);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector,
                                                                 ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();

    this->CalculateSoilDensity(rVariables);

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
            rRightHandSideVector[Index + idim] +=
                rVariables.Nu[i] * rVariables.Density * BodyAcceleration[idim] *
                rVariables.IntegrationCoefficientInitialConfiguration;
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateSoilDensity(ElementVariables& rVariables)
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    rVariables.Density = (rVariables.DegreeOfSaturation * rProp[POROSITY] * rProp[DENSITY_WATER]) +
                         (1.0 - rProp[POROSITY]) * rProp[DENSITY_SOLID];

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector,
                                                                  ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      VoigtSize = (Dim == N_DIM_3D ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);
    const SizeType StressTensorSize = (Dim == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D);

    Vector VoigtVector = ZeroVector(VoigtSize);
    for (SizeType idim = 0; idim < StressTensorSize; ++idim)
        VoigtVector[idim] = 1.0;

    Matrix CouplingMatrix = -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient * rVariables.BishopCoefficient *
                            prod(trans(rVariables.B), Matrix(outer_prod(VoigtVector, rVariables.Np))) *
                            rVariables.IntegrationCoefficient;

    Vector CouplingForce = prod(CouplingMatrix, rVariables.PressureVector);

    // Distribute coupling block vector 1 into the elemental vector
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, CouplingForce);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        Vector CouplingFlow = PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient *
                              prod(trans(CouplingMatrix), rVariables.VelocityVector);

        // Distribute coupling block vector 2 into the elemental vector
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, CouplingFlow);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector,
                                                                        ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix CompressibilityMatrix = -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotModulusInverse *
                                   outer_prod(rVariables.Np, rVariables.Np) * rVariables.IntegrationCoefficient;

    Vector CompressibilityFlow = -prod(CompressibilityMatrix, rVariables.PressureDtVector);

    // Distribute compressibility block vector into the elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, CompressibilityFlow);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                                                     ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix PermeabilityMatrix =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse *
        rVariables.RelativePermeability * rVariables.PermeabilityUpdateFactor *
        prod(rVariables.DNp_DX, Matrix(prod(rVariables.IntrinsicPermeability, trans(rVariables.DNp_DX)))) *
        rVariables.IntegrationCoefficient;

    Vector PermeabilityFlow = -prod(PermeabilityMatrix, rVariables.PressureVector);

    // Distribute permeability block vector into the elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, PermeabilityFlow);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                  ElementVariables& rVariables)
{
    KRATOS_TRY

    Matrix GradNpTPerm = rVariables.DynamicViscosityInverse * GetProperties()[DENSITY_WATER] *
                         rVariables.RelativePermeability * rVariables.PermeabilityUpdateFactor *
                         prod(rVariables.DNp_DX, rVariables.IntrinsicPermeability) *
                         rVariables.IntegrationCoefficient;

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    const SizeType      NumPNodes = mpPressureGeometry->PointsNumber();

    Vector BodyAcceleration = ZeroVector(Dim);

    SizeType Index = 0;
    for (SizeType i = 0; i < NumUNodes; ++i) {
        for (SizeType idim = 0; idim < Dim; ++idim) {
            BodyAcceleration[idim] += rVariables.Nu[i] * rVariables.BodyAcceleration[Index++];
        }
    }

    for (SizeType i = 0; i < NumPNodes; ++i) {
        rRightHandSideVector[NumUNodes * Dim + i] += inner_prod(row(GradNpTPerm, i), BodyAcceleration);
    }

    KRATOS_CATCH("")
}

GeometryData::IntegrationMethod SmallStrainUPwDiffOrderElement::GetIntegrationMethod() const
{
    GeometryData::IntegrationMethod GI_GAUSS;
    const GeometryType&             rGeom     = GetGeometry();
    const SizeType                  TNumNodes = rGeom.PointsNumber();
    //
    switch (TNumNodes) {
    case 3:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
        break;
    case 6:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
        break;
    case 10:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_4;
        break;
    case 15:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_5;
        break;
    default:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
        break;
    }

    return GI_GAUSS;
}

void SmallStrainUPwDiffOrderElement::CalculateStrain(ElementVariables& rVariables, unsigned int GPoint)
{
    if (rVariables.UseHenckyStrain) {
        this->CalculateDeformationGradient(rVariables, GPoint);
        const SizeType Dim       = GetGeometry().WorkingSpaceDimension();
        const SizeType VoigtSize = (Dim == N_DIM_3D ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);
        noalias(rVariables.StrainVector) = StressStrainUtilities::CalculateHenckyStrain(rVariables.F, VoigtSize);
    } else {
        this->CalculateCauchyStrain(rVariables);
    }
}

void SmallStrainUPwDiffOrderElement::CalculateCauchyStrain(ElementVariables& rVariables)
{
    noalias(rVariables.StrainVector) = prod(rVariables.B, rVariables.DisplacementVector);
}

Vector SmallStrainUPwDiffOrderElement::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient)
{
    return mpStressStatePolicy->CalculateGreenLagrangeStrain(rDeformationGradient);
}

void SmallStrainUPwDiffOrderElement::CalculateDeformationGradient(ElementVariables& rVariables, unsigned int GPoint)
{
    KRATOS_TRY

    // calculation of derivative of shape function with respect to reference
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
        << "This usually indicates that the deformations are too large for the "
           "mesh size."
        << std::endl;

    // Deformation gradient
    noalias(rVariables.F) = prod(J, InvJ0);
    rVariables.detF       = MathUtils<double>::Det(rVariables.F);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateJacobianOnCurrentConfiguration(double& detJ,
                                                                             Matrix& rJ,
                                                                             Matrix& rInvJ,
                                                                             unsigned int GPoint) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    rJ = rGeom.Jacobian(rJ, GPoint, this->GetIntegrationMethod());
    MathUtils<double>::InvertMatrix(rJ, rInvJ, detJ);

    KRATOS_CATCH("")
}

double SmallStrainUPwDiffOrderElement::CalculateFluidPressure(const ElementVariables& rVariables) const
{
    KRATOS_TRY

    return inner_prod(rVariables.Np, rVariables.PressureVector);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::SetRetentionParameters(const ElementVariables& rVariables,
                                                            RetentionLaw::Parameters& rRetentionParameters) const
{
    KRATOS_TRY

    rRetentionParameters.SetFluidPressure(rVariables.FluidPressure);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateRetentionResponse(ElementVariables& rVariables,
                                                                RetentionLaw::Parameters& rRetentionParameters,
                                                                unsigned int GPoint)
{
    KRATOS_TRY

    rVariables.FluidPressure = CalculateFluidPressure(rVariables);
    SetRetentionParameters(rVariables, rRetentionParameters);

    rVariables.DegreeOfSaturation = mRetentionLawVector[GPoint]->CalculateSaturation(rRetentionParameters);
    rVariables.DerivativeOfSaturation =
        mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(rRetentionParameters);
    rVariables.RelativePermeability =
        mRetentionLawVector[GPoint]->CalculateRelativePermeability(rRetentionParameters);
    rVariables.BishopCoefficient = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(rRetentionParameters);

    KRATOS_CATCH("")
}

Element::DofsVectorType SmallStrainUPwDiffOrderElement::GetDofs() const
{
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), *mpPressureGeometry,
                                                      GetGeometry().WorkingSpaceDimension());
}

const StressStatePolicy& SmallStrainUPwDiffOrderElement::GetStressStatePolicy() const
{
    return *mpStressStatePolicy;
}

} // Namespace Kratos