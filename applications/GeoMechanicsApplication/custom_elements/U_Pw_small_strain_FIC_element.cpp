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

// Application includes
#include "custom_elements/U_Pw_small_strain_FIC_element.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainFICElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                   NodesArrayType const& ThisNodes,
                                                                   PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainFICElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainFICElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                   GeometryType::Pointer pGeom,
                                                                   PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainFICElement(NewId, pGeom, pProperties,
                                                         this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    UPwBaseElement<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);

    for (unsigned int i = 0; i < TDim; ++i) {
        mNodalConstitutiveTensor[i].resize(VoigtSize);

        for (unsigned int j = 0; j < VoigtSize; j++) {
            for (unsigned int k = 0; k < TNumNodes; k++)
                mNodalConstitutiveTensor[i][j][k] = 0.0;
        }
    }

    for (unsigned int i = 0; i < TDim; ++i) {
        for (unsigned int j = 0; j < TNumNodes; j++)
            mNodalDtStress[i][j] = 0.0;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
int UPwSmallStrainFICElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Verify generic variables
    int ierr = UPwSmallStrainElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& Prop = this->GetProperties();

    // Verify specific properties
    if (Prop[IGNORE_UNDRAINED])
        KRATOS_ERROR << "IGNORE_UNDRAINED cannot be used in FIC elements. Use "
                        "Non FIC elements instead"
                     << this->Id() << std::endl;

    return ierr;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    UPwSmallStrainElement<TDim, TNumNodes>::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    UPwSmallStrainElement<TDim, TNumNodes>::FinalizeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Defining necessary variables
    const GeometryType& Geom       = this->GetGeometry();
    const SizeType      NumGPoints = Geom.IntegrationPointsNumber(mThisIntegrationMethod);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom, this->GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); // Note: this is for nonlocal damage

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Extrapolation variables
    array_1d<Matrix, TDim> ConstitutiveTensorContainer;
    for (unsigned int i = 0; i < TDim; ++i) {
        ConstitutiveTensorContainer[i].resize(NumGPoints, Variables.ConstitutiveMatrix.size1(), false);
    }
    Matrix DtStressContainer(NumGPoints, TDim);

    Vector StressVector(VoigtSize);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        // set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        this->SaveGPConstitutiveTensor(ConstitutiveTensorContainer, Variables.ConstitutiveMatrix, GPoint);

        // Compute DtStress
        noalias(Variables.StrainVector) = prod(Variables.B, Variables.VelocityVector);
        noalias(StressVector) = prod(Variables.ConstitutiveMatrix, Variables.StrainVector);

        this->SaveGPDtStress(DtStressContainer, StressVector, GPoint);
    }

    this->ExtrapolateGPConstitutiveTensor(ConstitutiveTensorContainer);
    this->ExtrapolateGPDtStress(DtStressContainer);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Defining necessary variables
    const GeometryType& Geom       = this->GetGeometry();
    const SizeType      NumGPoints = Geom.IntegrationPointsNumber(mThisIntegrationMethod);

    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom, this->GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); // Note: this is for nonlocal damage

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Containers for extrapolation variables
    Matrix DtStressContainer(NumGPoints, TDim);

    Vector StressVector(VoigtSize);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        // set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute ConstitutiveTensor
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        // Compute DtStress
        noalias(Variables.StrainVector) = prod(Variables.B, Variables.VelocityVector);
        noalias(StressVector) = prod(Variables.ConstitutiveMatrix, Variables.StrainVector);
        this->SaveGPDtStress(DtStressContainer, StressVector, GPoint);
    }

    this->ExtrapolateGPDtStress(DtStressContainer);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::SaveGPConstitutiveTensor(array_1d<Matrix, TDim>& rConstitutiveTensorContainer,
                                                                         const Matrix& ConstitutiveMatrix,
                                                                         const unsigned int& GPoint)
{
    KRATOS_TRY;

    for (unsigned int i = 0; i < TDim; ++i) {
        for (unsigned int j = 0; j < ConstitutiveMatrix.size1(); j++) {
            rConstitutiveTensorContainer[i](GPoint, j) = ConstitutiveMatrix(i, j);
        }
    }

    /* INFO: (Quadrilateral_2D_4 with GI_GAUSS_2)
     *
     *                                ( |D00-0 D01-0 D02-0|   |D10-0 D11-0 D12-0| )
     * rConstitutiveTensorContainer = ( |D00-1 D01-1 D02-1|   |D10-1 D11-1 D12-1| )
     *                                ( |D00-2 D01-2 D02-2|   |D10-2 D11-2 D12-2| )
     *                                ( |D00-3 D01-3 D02-3| , |D10-3 D11-3 D12-3| )
     *
     * D00-0 = D(0,0) at GP 0
     */
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::SaveGPDtStress(Matrix&       rDtStressContainer,
                                                               const Vector& StressVector,
                                                               const unsigned int& GPoint)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < TDim; ++i) {
        rDtStressContainer(GPoint, i) = StressVector[i];
    }

    /* INFO: (Quadrilateral_2D_4 with GI_GAUSS_2)
     *
     *                      |S0-0 S1-0|
     * rDtStressContainer = |S0-1 S1-1|
     *                      |S0-2 S1-2|
     *                      |S0-3 S1-3|
     *
     * S0-0 = S[0] at GP 0
     */

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 3>::ExtrapolateGPConstitutiveTensor(const array_1d<Matrix, 2>& ConstitutiveTensorContainer)
{
    // Triangle_2d_3 with GI_GAUSS_2
    KRATOS_TRY

    const SizeType Dim      = 2;
    const SizeType NumNodes = 3;

    BoundedMatrix<double, NumNodes, NumNodes> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double, NumNodes, VoigtSize> AuxNodalConstitutiveTensor;

    for (unsigned int i = 0; i < Dim; ++i) {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix, ConstitutiveTensorContainer[i]);

        for (unsigned int j = 0; j < VoigtSize; j++)
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor, j);
    }

    // clang-format off
    /* INFO:
     *
     *                            [ ( |D00-0|   |D01-0|   |D02-0| )   ( |D10-0|   |D11-0|   |D12-0| ) ]
     * mNodalConstitutiveTensor = [ ( |D00-1|   |D01-1|   |D02-1| )   ( |D10-1|   |D11-1|   |D12-1| ) ]
     *                            [ ( |D00-2| , |D01-2| , |D02-2| ) , ( |D10-2| , |D11-2| , |D12-2| ) ]
     *
     * D00-0 = D(0,0) at node 0
     */
    // clang-format on

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template <>
void UPwSmallStrainFICElement<2, 4>::ExtrapolateGPConstitutiveTensor(const array_1d<Matrix, 2>& ConstitutiveTensorContainer)
{
    KRATOS_TRY

    // Quadrilateral_2d_4 with GI_GAUSS_2
    const SizeType Dim      = 2;
    const SizeType NumNodes = 4;

    BoundedMatrix<double, NumNodes, NumNodes> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double, NumNodes, VoigtSize> AuxNodalConstitutiveTensor;

    for (unsigned int i = 0; i < Dim; ++i) {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix, ConstitutiveTensorContainer[i]);

        for (unsigned int j = 0; j < VoigtSize; j++)
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor, j);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template <>
void UPwSmallStrainFICElement<3, 4>::ExtrapolateGPConstitutiveTensor(const array_1d<Matrix, 3>& ConstitutiveTensorContainer)
{
    // Tetrahedra_3d_4 with GI_GAUSS_2
    KRATOS_TRY

    const SizeType Dim      = 3;
    const SizeType NumNodes = 4;

    BoundedMatrix<double, NumNodes, NumNodes> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double, NumNodes, VoigtSize> AuxNodalConstitutiveTensor;

    for (unsigned int i = 0; i < Dim; ++i) {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix, ConstitutiveTensorContainer[i]);

        for (unsigned int j = 0; j < VoigtSize; j++)
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor, j);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 8>::ExtrapolateGPConstitutiveTensor(const array_1d<Matrix, 3>& ConstitutiveTensorContainer)
{
    // Hexahedra_3d_8 with GI_GAUSS_2
    KRATOS_TRY

    const SizeType Dim      = 3;
    const SizeType NumNodes = 8;

    BoundedMatrix<double, NumNodes, NumNodes> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double, NumNodes, VoigtSize> AuxNodalConstitutiveTensor;

    for (unsigned int i = 0; i < Dim; ++i) {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix, ConstitutiveTensorContainer[i]);

        for (unsigned int j = 0; j < VoigtSize; j++)
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor, j);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::ExtrapolateGPDtStress(const Matrix& DtStressContainer)
{
    KRATOS_TRY

    BoundedMatrix<double, TNumNodes, TNumNodes> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double, TNumNodes, TDim> AuxNodalDtStress;
    noalias(AuxNodalDtStress) = prod(ExtrapolationMatrix, DtStressContainer);

    for (unsigned int i = 0; i < TDim; ++i)
        noalias(mNodalDtStress[i]) = column(AuxNodalDtStress, i);

    /* INFO:
     *
     *                  ( |S0-0|   |S1-0| )
     * mNodalDtStress = ( |S0-1|   |S1-1| )
     *                  ( |S0-2| , |S1-2| )
     *
     * S0-0 = S[0] at node 0
     */

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                             VectorType& rRightHandSideVector,
                                                             const ProcessInfo& CurrentProcessInfo,
                                                             const bool CalculateStiffnessMatrixFlag,
                                                             const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           Prop = this->GetProperties();
    const GeometryType&                             Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        Geom.IntegrationPoints(mThisIntegrationMethod);
    const SizeType NumGPoints = IntegrationPoints.size();

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom, Prop, CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, CurrentProcessInfo);

    FICElementVariables FICVariables;
    this->InitializeFICElementVariables(FICVariables, Variables.DN_DXContainer, Geom, Prop, CurrentProcessInfo);

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), CurrentProcessInfo);

    const bool hasBiotCoefficient = Prop.Has(BIOT_COEFFICIENT);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute Np, GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        // set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        // Compute ShapeFunctionsSecondOrderGradients
        this->CalculateShapeFunctionsSecondOrderGradients(FICVariables, Variables);

        // Compute constitutive tensor and stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // set shear modulus from stiffness matrix
        FICVariables.ShearModulus = CalculateShearModulus(Variables.ConstitutiveMatrix);

        // calculate Bulk modulus from stiffness matrix
        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

        Variables.IntegrationCoefficientInitialConfiguration = this->CalculateIntegrationCoefficient(
            IntegrationPoints, GPoint, Variables.detJInitialConfiguration);

        if (CalculateStiffnessMatrixFlag) {
            // Contributions to the left hand side
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
            this->CalculateAndAddLHSStabilization(rLeftHandSideMatrix, Variables, FICVariables);
        }

        if (CalculateResidualVectorFlag) {
            // Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
            this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateShearModulus(const Matrix& ConstitutiveMatrix) const
{
    const int IndexG = ConstitutiveMatrix.size1() - 1;
    return ConstitutiveMatrix(IndexG, IndexG);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::InitializeFICElementVariables(
    FICElementVariables&                             rFICVariables,
    const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer,
    const GeometryType&                              Geom,
    const PropertiesType&                            Prop,
    const ProcessInfo&                               CurrentProcessInfo)
{
    KRATOS_TRY

    // Nodal Variables
    this->ExtrapolateShapeFunctionsGradients(rFICVariables.NodalShapeFunctionsGradients, DN_DXContainer);

    // General Variables
    this->CalculateElementLength(rFICVariables.ElementLength, Geom);

    // Variables computed at each GP
    this->InitializeSecondOrderTerms(rFICVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 3>::ExtrapolateShapeFunctionsGradients(
    array_1d<array_1d<double, 6>, 3>&                rNodalShapeFunctionsGradients,
    const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer)
{
    // Triangle_2d_3 with GI_GAUSS_2
    // No necessary
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 4>::ExtrapolateShapeFunctionsGradients(
    array_1d<array_1d<double, 8>, 4>&                rNodalShapeFunctionsGradients,
    const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer)
{
    // Quadrilateral_2d_4 with GI_GAUSS_2
    KRATOS_TRY

    BoundedMatrix<double, 4, 8> ShapeFunctionsGradientsContainer; // NumGPoints X TDim*TNumNodes
    unsigned int                index;

    for (unsigned int i = 0; i < 4; ++i) // NumGPoints
    {
        for (unsigned int j = 0; j < 4; j++) // TNumNodes
        {
            index = j * 2;

            ShapeFunctionsGradientsContainer(i, index)     = DN_DXContainer[i](j, 0);
            ShapeFunctionsGradientsContainer(i, index + 1) = DN_DXContainer[i](j, 1);
        }
    }

    BoundedMatrix<double, 4, 4> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double, 4, 8> AuxNodalShapeFunctionsGradients;
    noalias(AuxNodalShapeFunctionsGradients) = prod(ExtrapolationMatrix, ShapeFunctionsGradientsContainer);

    for (unsigned int i = 0; i < 4; ++i) // TNumNodes
    {
        index = i * 2;

        rNodalShapeFunctionsGradients[i][0] = AuxNodalShapeFunctionsGradients(0, index);
        rNodalShapeFunctionsGradients[i][1] = AuxNodalShapeFunctionsGradients(0, index + 1);
        rNodalShapeFunctionsGradients[i][2] = AuxNodalShapeFunctionsGradients(1, index);
        rNodalShapeFunctionsGradients[i][3] = AuxNodalShapeFunctionsGradients(1, index + 1);
        rNodalShapeFunctionsGradients[i][4] = AuxNodalShapeFunctionsGradients(2, index);
        rNodalShapeFunctionsGradients[i][5] = AuxNodalShapeFunctionsGradients(2, index + 1);
        rNodalShapeFunctionsGradients[i][6] = AuxNodalShapeFunctionsGradients(3, index);
        rNodalShapeFunctionsGradients[i][7] = AuxNodalShapeFunctionsGradients(3, index + 1);
    }

    /* INFO:
     *
     *                                 ( |N0x-0|   |N1x-0|   |N2x-0|   |N3x-0| )
     *                                 ( |N0y-0|   |N1y-0|   |N2y-0|   |N3y-0| )
     *                                 ( |N0x-1|   |N1x-1|   |N2x-1|   |N3x-1| )
     * rNodalShapeFunctionsGradients = ( |N0y-1|   |N1y-1|   |N2y-1|   |N3y-1| )
     *                                 ( |N0x-2|   |N1x-2|   |N2x-2|   |N3x-2| )
     *                                 ( |N0y-2|   |N1y-2|   |N2y-2|   |N3y-2| )
     *                                 ( |N0x-3|   |N1x-3|   |N2x-3|   |N3x-3| )
     *                                 ( |N0y-3| , |N1y-3| , |N2y-3| , |N3y-3| )
     *
     * N0x-0 = aN0/ax at node 0
     */
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 4>::ExtrapolateShapeFunctionsGradients(
    array_1d<array_1d<double, 12>, 4>&               rNodalShapeFunctionsGradients,
    const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer)
{
    // Tetrahedra_3d_4 with GI_GAUSS_2
    // No necessary
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 8>::ExtrapolateShapeFunctionsGradients(
    array_1d<array_1d<double, 24>, 8>&               rNodalShapeFunctionsGradients,
    const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer)
{
    // Hexahedra_3d_8 with GI_GAUSS_2
    KRATOS_TRY

    BoundedMatrix<double, 8, 24> ShapeFunctionsGradientsContainer; // NumGPoints X TDim*TNumNodes
    unsigned int                 index;

    for (unsigned int i = 0; i < 8; ++i) // NumGPoints
    {
        for (unsigned int j = 0; j < 8; j++) // TNumNodes
        {
            index = j * 3;

            ShapeFunctionsGradientsContainer(i, index)     = DN_DXContainer[i](j, 0);
            ShapeFunctionsGradientsContainer(i, index + 1) = DN_DXContainer[i](j, 1);
            ShapeFunctionsGradientsContainer(i, index + 2) = DN_DXContainer[i](j, 2);
        }
    }

    BoundedMatrix<double, 8, 8> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double, 8, 24> AuxNodalShapeFunctionsGradients;
    noalias(AuxNodalShapeFunctionsGradients) = prod(ExtrapolationMatrix, ShapeFunctionsGradientsContainer);

    for (unsigned int i = 0; i < 8; ++i) // TNumNodes
    {
        index = i * 3;

        rNodalShapeFunctionsGradients[i][0]  = AuxNodalShapeFunctionsGradients(0, index);
        rNodalShapeFunctionsGradients[i][1]  = AuxNodalShapeFunctionsGradients(0, index + 1);
        rNodalShapeFunctionsGradients[i][2]  = AuxNodalShapeFunctionsGradients(0, index + 2);
        rNodalShapeFunctionsGradients[i][3]  = AuxNodalShapeFunctionsGradients(1, index);
        rNodalShapeFunctionsGradients[i][4]  = AuxNodalShapeFunctionsGradients(1, index + 1);
        rNodalShapeFunctionsGradients[i][5]  = AuxNodalShapeFunctionsGradients(1, index + 2);
        rNodalShapeFunctionsGradients[i][6]  = AuxNodalShapeFunctionsGradients(2, index);
        rNodalShapeFunctionsGradients[i][7]  = AuxNodalShapeFunctionsGradients(2, index + 1);
        rNodalShapeFunctionsGradients[i][8]  = AuxNodalShapeFunctionsGradients(2, index + 2);
        rNodalShapeFunctionsGradients[i][9]  = AuxNodalShapeFunctionsGradients(3, index);
        rNodalShapeFunctionsGradients[i][10] = AuxNodalShapeFunctionsGradients(3, index + 1);
        rNodalShapeFunctionsGradients[i][11] = AuxNodalShapeFunctionsGradients(3, index + 2);
        rNodalShapeFunctionsGradients[i][12] = AuxNodalShapeFunctionsGradients(4, index);
        rNodalShapeFunctionsGradients[i][13] = AuxNodalShapeFunctionsGradients(4, index + 1);
        rNodalShapeFunctionsGradients[i][14] = AuxNodalShapeFunctionsGradients(4, index + 2);
        rNodalShapeFunctionsGradients[i][15] = AuxNodalShapeFunctionsGradients(5, index);
        rNodalShapeFunctionsGradients[i][16] = AuxNodalShapeFunctionsGradients(5, index + 1);
        rNodalShapeFunctionsGradients[i][17] = AuxNodalShapeFunctionsGradients(5, index + 2);
        rNodalShapeFunctionsGradients[i][18] = AuxNodalShapeFunctionsGradients(6, index);
        rNodalShapeFunctionsGradients[i][19] = AuxNodalShapeFunctionsGradients(6, index + 1);
        rNodalShapeFunctionsGradients[i][20] = AuxNodalShapeFunctionsGradients(6, index + 2);
        rNodalShapeFunctionsGradients[i][21] = AuxNodalShapeFunctionsGradients(7, index);
        rNodalShapeFunctionsGradients[i][22] = AuxNodalShapeFunctionsGradients(7, index + 1);
        rNodalShapeFunctionsGradients[i][23] = AuxNodalShapeFunctionsGradients(7, index + 2);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 3>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = sqrt(4.0 * Geom.Area() / Globals::Pi);
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 4>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = sqrt(4.0 * Geom.Area() / Globals::Pi);
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 4>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = pow((6.0 * Geom.Volume() / Globals::Pi), (1.0 / 3.0));
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 8>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = pow((6.0 * Geom.Volume() / Globals::Pi), (1.0 / 3.0));
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 3>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    KRATOS_TRY

    const SizeType Dim = 2;

    for (unsigned int i = 0; i < Dim; ++i)
        rFICVariables.ConstitutiveTensorGradients[i].resize(VoigtSize);

    rFICVariables.DimVoigtMatrix.resize(Dim, VoigtSize, false);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 4>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    KRATOS_TRY
    const SizeType Dim      = 2;
    const SizeType NumNodes = 4;

    // Voigt identity matrix
    rFICVariables.VoigtMatrix.resize(VoigtSize, VoigtSize, false);
    noalias(rFICVariables.VoigtMatrix) = ZeroMatrix(VoigtSize, VoigtSize);
    rFICVariables.VoigtMatrix(0, 0)    = 1.0;
    rFICVariables.VoigtMatrix(1, 1)    = 1.0;
    rFICVariables.VoigtMatrix(2, 2)    = 0.5;

    for (unsigned int i = 0; i < NumNodes; ++i)
        rFICVariables.ShapeFunctionsSecondOrderGradients[i].resize(VoigtSize, false);

    for (unsigned int i = 0; i < Dim; ++i)
        rFICVariables.ConstitutiveTensorGradients[i].resize(VoigtSize);

    rFICVariables.DimVoigtMatrix.resize(Dim, VoigtSize, false);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 4>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    KRATOS_TRY

    const SizeType Dim = 3;

    for (unsigned int i = 0; i < Dim; ++i)
        rFICVariables.ConstitutiveTensorGradients[i].resize(VoigtSize);

    rFICVariables.DimVoigtMatrix.resize(Dim, VoigtSize, false);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 8>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    KRATOS_TRY

    const SizeType Dim      = 3;
    const SizeType NumNodes = 8;

    // Voigt identity matrix
    rFICVariables.VoigtMatrix.resize(VoigtSize, VoigtSize, false);
    noalias(rFICVariables.VoigtMatrix) = ZeroMatrix(VoigtSize, VoigtSize);
    rFICVariables.VoigtMatrix(0, 0)    = 1.0;
    rFICVariables.VoigtMatrix(1, 1)    = 1.0;
    rFICVariables.VoigtMatrix(2, 2)    = 1.0;
    rFICVariables.VoigtMatrix(3, 3)    = 0.5;
    rFICVariables.VoigtMatrix(4, 4)    = 0.5;
    rFICVariables.VoigtMatrix(5, 5)    = 0.5;

    for (unsigned int i = 0; i < NumNodes; ++i)
        rFICVariables.ShapeFunctionsSecondOrderGradients[i].resize(VoigtSize, false);

    for (unsigned int i = 0; i < Dim; ++i)
        rFICVariables.ConstitutiveTensorGradients[i].resize(VoigtSize);

    rFICVariables.DimVoigtMatrix.resize(Dim, VoigtSize, false);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 3>::CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables,
                                                                                 ElementVariables& rVariables)
{
    // Not necessary
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 4>::CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables,
                                                                                 ElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B), rFICVariables.VoigtMatrix);
    unsigned int index;
    for (unsigned int i = 0; i < 4; ++i) // TNumNodes
    {
        index = 2 * i;

        noalias(rFICVariables.ShapeFunctionsSecondOrderGradients[i]) =
            prod(trans(rVariables.UVoigtMatrix), rFICVariables.NodalShapeFunctionsGradients[i]);

        rFICVariables.StrainGradients(0, index) =
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] +
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][1];
        rFICVariables.StrainGradients(1, index + 1) =
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][1];
        rFICVariables.StrainGradients(0, index + 1) =
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][2];
        rFICVariables.StrainGradients(1, index) = rFICVariables.StrainGradients(0, index + 1);
    }

    /* INFO:
     *
     *                                                    ( |N0xx|   |N1xx|   |N2xx|   |N3xx| )
     * rFICVariables.ShapeFunctionsSecondOrderGradients = ( |N0yy|   |N1yy|   |N2yy|   |N3yy| )
     *                                                    ( |N0xy| , |N1xy| , |N2xy| , |N3xy| )
     *
     * N0xx = a2N0/ax2 at current GP
     */

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 4>::CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables,
                                                                                 ElementVariables& rVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 8>::CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables,
                                                                                 ElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B), rFICVariables.VoigtMatrix);
    unsigned int index;
    for (unsigned int i = 0; i < 8; ++i) // TNumNodes
    {
        index = 3 * i;

        noalias(rFICVariables.ShapeFunctionsSecondOrderGradients[i]) =
            prod(trans(rVariables.UVoigtMatrix), rFICVariables.NodalShapeFunctionsGradients[i]);

        rFICVariables.StrainGradients(0, index) =
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] +
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][1] +
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][2];
        rFICVariables.StrainGradients(1, index + 1) =
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][1] +
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][2];
        rFICVariables.StrainGradients(2, index + 2) =
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] +
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][1] +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][2];
        rFICVariables.StrainGradients(0, index + 1) =
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][3];
        rFICVariables.StrainGradients(1, index) = rFICVariables.StrainGradients(0, index + 1);
        rFICVariables.StrainGradients(1, index + 2) =
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][4];
        rFICVariables.StrainGradients(2, index + 1) = rFICVariables.StrainGradients(1, index + 2);
        rFICVariables.StrainGradients(0, index + 2) =
            0.5 * rFICVariables.ShapeFunctionsSecondOrderGradients[i][5];
        rFICVariables.StrainGradients(2, index) = rFICVariables.StrainGradients(0, index + 2);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateAndAddLHSStabilization(
    MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    KRATOS_TRY;

    this->CalculateAndAddStrainGradientMatrix(rLeftHandSideMatrix, rVariables, rFICVariables);

    this->CalculateAndAddDtStressGradientMatrix(rLeftHandSideMatrix, rVariables, rFICVariables);

    this->CalculateAndAddPressureGradientMatrix(rLeftHandSideMatrix, rVariables, rFICVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 3>::CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix,
                                                                         ElementVariables& rVariables,
                                                                         FICElementVariables& rFICVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 4>::CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix,
                                                                         ElementVariables& rVariables,
                                                                         FICElementVariables& rFICVariables)
{
    KRATOS_TRY

    noalias(rVariables.PUMatrix) =
        PORE_PRESSURE_SIGN_FACTOR * rVariables.VelocityCoefficient * 0.25 *
        rFICVariables.ElementLength * rFICVariables.ElementLength * rVariables.BiotCoefficient *
        prod(rVariables.GradNpT, rFICVariables.StrainGradients) * rVariables.IntegrationCoefficient;

    // Distribute strain gradient matrix into the elemental matrix
    GeoElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix, rVariables.PUMatrix);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template <>
void UPwSmallStrainFICElement<3, 4>::CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix,
                                                                         ElementVariables& rVariables,
                                                                         FICElementVariables& rFICVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 8>::CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix,
                                                                         ElementVariables& rVariables,
                                                                         FICElementVariables& rFICVariables)
{
    KRATOS_TRY

    noalias(rVariables.PUMatrix) =
        PORE_PRESSURE_SIGN_FACTOR * rVariables.VelocityCoefficient * 0.25 *
        rFICVariables.ElementLength * rFICVariables.ElementLength * rVariables.BiotCoefficient *
        prod(rVariables.GradNpT, rFICVariables.StrainGradients) * rVariables.IntegrationCoefficient;

    // Distribute strain gradient matrix into the elemental matrix
    GeoElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix, rVariables.PUMatrix);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateAndAddDtStressGradientMatrix(
    MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    KRATOS_TRY;

    this->CalculateConstitutiveTensorGradients(rFICVariables, rVariables);

    double StabilizationParameter = -PORE_PRESSURE_SIGN_FACTOR * rFICVariables.ElementLength *
                                    rFICVariables.ElementLength * rVariables.BiotCoefficient /
                                    (8.0 * rFICVariables.ShearModulus);

    noalias(rVariables.PUMatrix) = -rVariables.VelocityCoefficient * StabilizationParameter / 3.0 *
                                   prod(rVariables.GradNpT, rFICVariables.DimUMatrix) *
                                   rVariables.IntegrationCoefficient;

    // Distribute DtStressGradient Matrix into the elemental matrix
    GeoElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix, rVariables.PUMatrix);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 3>::CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables,
                                                                          const ElementVariables& Variables)
{
    KRATOS_TRY

    const SizeType Dim      = 2;
    const SizeType NumNodes = 3;

    for (unsigned int i = 0; i < Dim; ++i) {
        for (unsigned int j = 0; j < VoigtSize; j++) {
            for (unsigned int k = 0; k < Dim; k++) {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;

                for (unsigned int l = 0; l < NumNodes; l++)
                    (rFICVariables.ConstitutiveTensorGradients[i][j][k]) +=
                        Variables.GradNpT(l, k) * (mNodalConstitutiveTensor[i][j][l]);
            }
        }
    }

    for (unsigned int i = 0; i < Dim; ++i) {
        for (unsigned int j = 0; j < VoigtSize; j++) {
            rFICVariables.DimVoigtMatrix(i, j) = 0.0;

            for (unsigned int k = 0; k < Dim; k++)
                rFICVariables.DimVoigtMatrix(i, j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix, Variables.B);

    // clang-format off
    /* INFO:
     *
     * rFICVariables.ConstitutiveTensorGradients = [ ( |D00x|   |D01x|   |D02x| )   ( |D10x|   |D11x|   |D12x| ) ]
     *                                             [ ( |D00y| , |D01y| , |D02y| ) , ( |D10y| , |D11y| , |D12y| ) ]
     *
     * rFICVariables.DimVoigtMatrix = |D00x+D10x D01x+D11x D02x+D12x|
     *                                |D00y+D10y D01y+D11y D02y+D12y|
     *
     * D00x = aD(0,0)/ax at current GP
     */
    // clang-format on

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 4>::CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables,
                                                                          const ElementVariables& Variables)
{
    KRATOS_TRY

    const SizeType Dim      = 2;
    const SizeType NumNodes = 4;

    for (unsigned int i = 0; i < Dim; ++i) {
        for (unsigned int j = 0; j < VoigtSize; j++) {
            for (unsigned int k = 0; k < Dim; k++) {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;

                for (unsigned int l = 0; l < NumNodes; l++)
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] +=
                        Variables.GradNpT(l, k) * mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for (unsigned int i = 0; i < Dim; ++i) {
        for (unsigned int j = 0; j < VoigtSize; j++) {
            rFICVariables.DimVoigtMatrix(i, j) = 0.0;

            for (unsigned int k = 0; k < Dim; k++)
                rFICVariables.DimVoigtMatrix(i, j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix, Variables.B);

    // Adding ShapeFunctionsSecondOrderGradients terms
    unsigned int index;

    for (unsigned int i = 0; i < NumNodes; ++i) {
        index = Dim * i;

        rFICVariables.DimUMatrix(0, index) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] *
                (Variables.ConstitutiveMatrix(0, 0) + Variables.ConstitutiveMatrix(1, 0)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][2] *
                (Variables.ConstitutiveMatrix(0, 2) + Variables.ConstitutiveMatrix(1, 2));
        rFICVariables.DimUMatrix(0, index + 1) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] *
                (Variables.ConstitutiveMatrix(0, 2) + Variables.ConstitutiveMatrix(1, 2)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][2] *
                (Variables.ConstitutiveMatrix(0, 1) + Variables.ConstitutiveMatrix(1, 1));
        rFICVariables.DimUMatrix(1, index) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][1] *
                (Variables.ConstitutiveMatrix(0, 2) + Variables.ConstitutiveMatrix(1, 2)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][2] *
                (Variables.ConstitutiveMatrix(0, 0) + Variables.ConstitutiveMatrix(1, 0));
        rFICVariables.DimUMatrix(1, index + 1) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][1] *
                (Variables.ConstitutiveMatrix(0, 1) + Variables.ConstitutiveMatrix(1, 1)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][2] *
                (Variables.ConstitutiveMatrix(0, 2) + Variables.ConstitutiveMatrix(1, 2));
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 4>::CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables,
                                                                          const ElementVariables& Variables)
{
    KRATOS_TRY

    const SizeType Dim      = 3;
    const SizeType NumNodes = 4;

    for (unsigned int i = 0; i < Dim; ++i) {
        for (unsigned int j = 0; j < VoigtSize; j++) {
            for (unsigned int k = 0; k < Dim; k++) {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;

                for (unsigned int l = 0; l < NumNodes; l++)
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] +=
                        Variables.GradNpT(l, k) * mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for (unsigned int i = 0; i < Dim; ++i) {
        for (unsigned int j = 0; j < VoigtSize; j++) {
            rFICVariables.DimVoigtMatrix(i, j) = 0.0;

            for (unsigned int k = 0; k < Dim; k++)
                rFICVariables.DimVoigtMatrix(i, j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix, Variables.B);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 8>::CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables,
                                                                          const ElementVariables& Variables)
{
    KRATOS_TRY

    const SizeType Dim      = 3;
    const SizeType NumNodes = 8;

    for (unsigned int i = 0; i < Dim; ++i) {
        for (unsigned int j = 0; j < VoigtSize; j++) {
            for (unsigned int k = 0; k < Dim; k++) {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;

                for (unsigned int l = 0; l < NumNodes; l++)
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] +=
                        Variables.GradNpT(l, k) * mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for (unsigned int i = 0; i < Dim; ++i) {
        for (unsigned int j = 0; j < VoigtSize; j++) {
            rFICVariables.DimVoigtMatrix(i, j) = 0.0;

            for (unsigned int k = 0; k < Dim; k++)
                rFICVariables.DimVoigtMatrix(i, j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix, Variables.B);

    // Adding ShapeFunctionsSecondOrderGradients terms
    unsigned int index;

    for (unsigned int i = 0; i < NumNodes; ++i) {
        index = Dim * i;

        rFICVariables.DimUMatrix(0, index) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] *
                (Variables.ConstitutiveMatrix(0, 0) + Variables.ConstitutiveMatrix(1, 0) +
                 Variables.ConstitutiveMatrix(2, 0)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][3] *
                (Variables.ConstitutiveMatrix(0, 3) + Variables.ConstitutiveMatrix(1, 3) +
                 Variables.ConstitutiveMatrix(2, 3)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][5] *
                (Variables.ConstitutiveMatrix(0, 5) + Variables.ConstitutiveMatrix(1, 5) +
                 Variables.ConstitutiveMatrix(2, 5));
        rFICVariables.DimUMatrix(0, index + 1) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] *
                (Variables.ConstitutiveMatrix(0, 3) + Variables.ConstitutiveMatrix(1, 3) +
                 Variables.ConstitutiveMatrix(2, 3)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][3] *
                (Variables.ConstitutiveMatrix(0, 1) + Variables.ConstitutiveMatrix(1, 1) +
                 Variables.ConstitutiveMatrix(2, 1)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][5] *
                (Variables.ConstitutiveMatrix(0, 4) + Variables.ConstitutiveMatrix(1, 4) +
                 Variables.ConstitutiveMatrix(2, 4));
        rFICVariables.DimUMatrix(0, index + 2) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][0] *
                (Variables.ConstitutiveMatrix(0, 5) + Variables.ConstitutiveMatrix(1, 5) +
                 Variables.ConstitutiveMatrix(2, 5)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][3] *
                (Variables.ConstitutiveMatrix(0, 4) + Variables.ConstitutiveMatrix(1, 4) +
                 Variables.ConstitutiveMatrix(2, 4)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][5] *
                (Variables.ConstitutiveMatrix(0, 2) + Variables.ConstitutiveMatrix(1, 2) +
                 Variables.ConstitutiveMatrix(2, 2));
        rFICVariables.DimUMatrix(1, index) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][1] *
                (Variables.ConstitutiveMatrix(0, 3) + Variables.ConstitutiveMatrix(1, 3) +
                 Variables.ConstitutiveMatrix(2, 3)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][3] *
                (Variables.ConstitutiveMatrix(0, 0) + Variables.ConstitutiveMatrix(1, 0) +
                 Variables.ConstitutiveMatrix(2, 0)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][4] *
                (Variables.ConstitutiveMatrix(0, 5) + Variables.ConstitutiveMatrix(1, 5) +
                 Variables.ConstitutiveMatrix(2, 5));
        rFICVariables.DimUMatrix(1, index + 1) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][1] *
                (Variables.ConstitutiveMatrix(0, 1) + Variables.ConstitutiveMatrix(1, 1) +
                 Variables.ConstitutiveMatrix(2, 1)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][3] *
                (Variables.ConstitutiveMatrix(0, 3) + Variables.ConstitutiveMatrix(1, 3) +
                 Variables.ConstitutiveMatrix(2, 3)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][4] *
                (Variables.ConstitutiveMatrix(0, 4) + Variables.ConstitutiveMatrix(1, 4) +
                 Variables.ConstitutiveMatrix(2, 4));
        rFICVariables.DimUMatrix(1, index + 2) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][1] *
                (Variables.ConstitutiveMatrix(0, 4) + Variables.ConstitutiveMatrix(1, 4) +
                 Variables.ConstitutiveMatrix(2, 4)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][3] *
                (Variables.ConstitutiveMatrix(0, 5) + Variables.ConstitutiveMatrix(1, 5) +
                 Variables.ConstitutiveMatrix(2, 5)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][4] *
                (Variables.ConstitutiveMatrix(0, 2) + Variables.ConstitutiveMatrix(1, 2) +
                 Variables.ConstitutiveMatrix(2, 2));
        rFICVariables.DimUMatrix(2, index) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][2] *
                (Variables.ConstitutiveMatrix(0, 5) + Variables.ConstitutiveMatrix(1, 5) +
                 Variables.ConstitutiveMatrix(2, 5)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][4] *
                (Variables.ConstitutiveMatrix(0, 3) + Variables.ConstitutiveMatrix(1, 3) +
                 Variables.ConstitutiveMatrix(2, 3)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][5] *
                (Variables.ConstitutiveMatrix(0, 0) + Variables.ConstitutiveMatrix(1, 0) +
                 Variables.ConstitutiveMatrix(2, 0));
        rFICVariables.DimUMatrix(2, index + 1) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][2] *
                (Variables.ConstitutiveMatrix(0, 4) + Variables.ConstitutiveMatrix(1, 4) +
                 Variables.ConstitutiveMatrix(2, 4)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][4] *
                (Variables.ConstitutiveMatrix(0, 1) + Variables.ConstitutiveMatrix(1, 1) +
                 Variables.ConstitutiveMatrix(2, 1)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][5] *
                (Variables.ConstitutiveMatrix(0, 3) + Variables.ConstitutiveMatrix(1, 3) +
                 Variables.ConstitutiveMatrix(2, 3));
        rFICVariables.DimUMatrix(2, index + 2) +=
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][2] *
                (Variables.ConstitutiveMatrix(0, 2) + Variables.ConstitutiveMatrix(1, 2) +
                 Variables.ConstitutiveMatrix(2, 2)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][4] *
                (Variables.ConstitutiveMatrix(0, 4) + Variables.ConstitutiveMatrix(1, 4) +
                 Variables.ConstitutiveMatrix(2, 4)) +
            rFICVariables.ShapeFunctionsSecondOrderGradients[i][5] *
                (Variables.ConstitutiveMatrix(0, 5) + Variables.ConstitutiveMatrix(1, 5) +
                 Variables.ConstitutiveMatrix(2, 5));
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateAndAddPressureGradientMatrix(
    MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    KRATOS_TRY;

    const double SignBiotCoefficient = -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient;

    const double StabilizationParameter = rFICVariables.ElementLength * rFICVariables.ElementLength *
                                          SignBiotCoefficient / (8.0 * rFICVariables.ShearModulus);

    noalias(rVariables.PPMatrix) =
        rVariables.DtPressureCoefficient * StabilizationParameter *
        (SignBiotCoefficient - 2.0 * rFICVariables.ShearModulus * rVariables.BiotModulusInverse /
                                   (3.0 * SignBiotCoefficient)) *
        prod(rVariables.GradNpT, trans(rVariables.GradNpT)) * rVariables.IntegrationCoefficient;

    // Distribute pressure gradient block matrix into the elemental matrix
    GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, rVariables.PPMatrix);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateAndAddRHSStabilization(
    VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    KRATOS_TRY;

    this->CalculateAndAddStrainGradientFlow(rRightHandSideVector, rVariables, rFICVariables);

    this->CalculateAndAddDtStressGradientFlow(rRightHandSideVector, rVariables, rFICVariables);

    this->CalculateAndAddPressureGradientFlow(rRightHandSideVector, rVariables, rFICVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 3>::CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector,
                                                                       ElementVariables& rVariables,
                                                                       FICElementVariables& rFICVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<2, 4>::CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector,
                                                                       ElementVariables& rVariables,
                                                                       FICElementVariables& rFICVariables)
{
    KRATOS_TRY

    noalias(rVariables.PUMatrix) = 0.25 * rFICVariables.ElementLength * rFICVariables.ElementLength *
                                   rVariables.BiotCoefficient * (-PORE_PRESSURE_SIGN_FACTOR) *
                                   prod(rVariables.GradNpT, rFICVariables.StrainGradients) *
                                   rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = prod(rVariables.PUMatrix, rVariables.VelocityVector);

    // Distribute Strain Gradient vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, rVariables.PVector);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 4>::CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector,
                                                                       ElementVariables& rVariables,
                                                                       FICElementVariables& rFICVariables)
{
    // Not necessary
}

//----------------------------------------------------------------------------------------
template <>
void UPwSmallStrainFICElement<3, 8>::CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector,
                                                                       ElementVariables& rVariables,
                                                                       FICElementVariables& rFICVariables)
{
    KRATOS_TRY

    noalias(rVariables.PUMatrix) = 0.25 * rFICVariables.ElementLength * rFICVariables.ElementLength *
                                   rVariables.BiotCoefficient * (-PORE_PRESSURE_SIGN_FACTOR) *
                                   prod(rVariables.GradNpT, rFICVariables.StrainGradients) *
                                   rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = prod(rVariables.PUMatrix, rVariables.VelocityVector);

    // Distribute Strain Gradient vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, rVariables.PVector);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateAndAddDtStressGradientFlow(
    VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    KRATOS_TRY;

    this->CalculateDtStressGradients(rFICVariables, rVariables);

    double StabilizationParameter = rFICVariables.ElementLength * rFICVariables.ElementLength *
                                    (rVariables.BiotCoefficient * (-PORE_PRESSURE_SIGN_FACTOR)) /
                                    (8.0 * rFICVariables.ShearModulus);

    noalias(rVariables.PVector) = StabilizationParameter / 3.0 *
                                  prod(rVariables.GradNpT, rFICVariables.DimVector) *
                                  rVariables.IntegrationCoefficient;

    // Distribute DtStressGradient block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, rVariables.PVector);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateDtStressGradients(FICElementVariables& rFICVariables,
                                                                           const ElementVariables& Variables)
{
    KRATOS_TRY;

    for (unsigned int i = 0; i < TDim; ++i) {
        for (unsigned int j = 0; j < TDim; j++) {
            rFICVariables.DtStressGradients[i][j] = 0.0;

            for (unsigned int k = 0; k < TNumNodes; k++)
                rFICVariables.DtStressGradients[i][j] += Variables.GradNpT(k, j) * mNodalDtStress[i][k];
        }
    }

    for (unsigned int i = 0; i < TDim; ++i) {
        rFICVariables.DimVector[i] = 0.0;

        for (unsigned int j = 0; j < TDim; j++)
            rFICVariables.DimVector[i] += rFICVariables.DtStressGradients[j][i];
    }

    /* INFO: (2D)
     *
     * rFICVariables.DtStressGradients = ( |S0x|   |S1x| )
     *                                   ( |S0y| , |S1y| )
     *
     * rFICVariables.DimVector = |S0x+S1x|
     *                           |S0y+S1y|
     *
     * S0x = aS[0]/ax at current GP
     */

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateAndAddPressureGradientFlow(
    VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    KRATOS_TRY;

    double SignBiotCoefficient    = -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient;
    double StabilizationParameter = rFICVariables.ElementLength * rFICVariables.ElementLength *
                                    SignBiotCoefficient / (8.0 * rFICVariables.ShearModulus);

    noalias(rVariables.PPMatrix) =
        StabilizationParameter *
        (SignBiotCoefficient - 2.0 * rFICVariables.ShearModulus * rVariables.BiotModulusInverse /
                                   (3.0 * SignBiotCoefficient)) *
        prod(rVariables.GradNpT, trans(rVariables.GradNpT)) * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0 * prod(rVariables.PPMatrix, rVariables.DtPressureVector);

    // Distribute PressureGradient block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, rVariables.PVector);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------

template class UPwSmallStrainFICElement<2, 3>;
template class UPwSmallStrainFICElement<2, 4>;
template class UPwSmallStrainFICElement<3, 4>;
template class UPwSmallStrainFICElement<3, 8>;

} // Namespace Kratos
