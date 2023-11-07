// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Gennady Markelov
//

#include "custom_elements/transient_thermal_element.h"
#include "custom_constitutive/thermal_dispersion_law.h"
#include "includes/condition.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(IndexType NewId)
    : Element(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(IndexType NewId,
                                                                  const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(IndexType NewId,
                                                                  GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::~TransientThermalElement() = default;

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientThermalElement(
        NewId, GetGeometry().Create(rThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientThermalElement(NewId, pGeom, pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                          const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const unsigned int N_DOF = GetNumberOfDOF();
    if (rElementalDofList.size() != N_DOF) {
        rElementalDofList.resize(N_DOF);
    }

    const GeometryType& rGeom = GetGeometry();
    for (unsigned int i = 0; i < N_DOF; ++i) {
        rElementalDofList[i] = rGeom[i].pGetDof(TEMPERATURE);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const unsigned int N_DOF = GetNumberOfDOF();
    if (rResult.size() != N_DOF) {
        rResult.resize(N_DOF, false);
    }

    const GeometryType& rGeom = GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = rGeom[i].GetDof(TEMPERATURE).EquationId();
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mIsInitialised = true;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
int TransientThermalElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    CheckDomainSize();
    const GeometryType& rGeom = GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        CheckSolutionStepsData(i, TEMPERATURE);
        CheckSolutionStepsData(i, DT_TEMPERATURE);
        if (!rGeom[i].HasDofFor(TEMPERATURE)) {
            KRATOS_ERROR << "missing degree of freedom for TEMPERATURE on node "
                         << rGeom[i].Id() << std::endl;
        }
    }

    VerifyProperty(DENSITY_WATER);
    VerifyProperty(POROSITY);
    VerifyProperty(SATURATION);
    VerifyProperty(DENSITY_SOLID);
    VerifyProperty(SPECIFIC_HEAT_CAPACITY_WATER);
    VerifyProperty(SPECIFIC_HEAT_CAPACITY_SOLID);
    VerifyProperty(THERMAL_CONDUCTIVITY_WATER);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XX);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_YY);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XY);
    VerifyProperty(LONGITUDINAL_DISPERSIVITY);
    VerifyProperty(TRANSVERSE_DISPERSIVITY);
    VerifyProperty(SOLID_COMPRESSIBILITY);

    if (TDim == 2) {
        auto pos = std::find_if(rGeom.begin(), rGeom.end(),
                                [](const auto& node) { return node.Z() != 0.0; });
        if (pos != rGeom.end()) {
            KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << pos->Id()
                         << std::endl;
        }
    }

    if (TDim > 2) {
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_ZZ);
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_YZ);
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XZ);
    }

    KRATOS_CATCH("");

    return 0;
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();

    ElementVariables Variables;
    InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute GradNpT, B and StrainVector
        CalculateKinematics(Variables, GPoint);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

        if (CalculateStiffnessMatrixFlag) {
            CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
        }

        if (CalculateResidualVectorFlag) {
            CalculateAndAddRHS(rRightHandSideVector, Variables);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::InitializeElementVariables(
    ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    InitializeProperties(rVariables);

    rVariables.DtTemperatureCoefficient = rCurrentProcessInfo[DT_TEMPERATURE_COEFFICIENT];

    InitializeNodalTemperatureVariables(rVariables);

    // Variables computed at each GP
    rVariables.N.resize(TNumNodes, false);
    rVariables.GradNT.resize(TNumNodes, TDim, false);

    const GeometryType& rGeom = GetGeometry();
    const unsigned int NumGPoints =
        rGeom.IntegrationPointsNumber(GetIntegrationMethod());

    // shape functions
    rVariables.NContainer.resize(NumGPoints, TNumNodes, false);
    rVariables.NContainer = rGeom.ShapeFunctionsValues(GetIntegrationMethod());

    // gradient of shape functions and determinant of Jacobian
    rVariables.detJContainer.resize(NumGPoints, false);

    rGeom.ShapeFunctionsIntegrationPointsGradients(
        rVariables.DN_DXContainer, rVariables.detJContainer, GetIntegrationMethod());

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                  ElementVariables& rVariables)
{
    KRATOS_TRY

    CalculateAndAddConductivityMatrix(rLeftHandSideMatrix, rVariables);
    CalculateAndAddCapacityMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddConductivityMatrix(
    MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY

    CalculateConductivityMatrix(rVariables);
    GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(
        rLeftHandSideMatrix, rVariables.ConductivityMatrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddCapacityMatrix(
    MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY

    CalculateCapacityMatrix(rVariables);
    GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(
        rLeftHandSideMatrix, rVariables.CapacityMatrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                  ElementVariables& rVariables)
{
    KRATOS_TRY

    CalculateAndAddCapacityVector(rRightHandSideVector, rVariables);
    CalculateAndAddConductivityVector(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateKinematics(ElementVariables& rVariables,
                                                                   unsigned int PointNumber)
{
    KRATOS_TRY

    // Setting the vector of shape functions and the matrix of the shape functions global gradients
    rVariables.N = row(rVariables.NContainer, PointNumber);
    rVariables.GradNT = rVariables.DN_DXContainer[PointNumber];
    rVariables.detJ = rVariables.detJContainer[PointNumber];

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
unsigned int TransientThermalElement<TDim, TNumNodes>::GetNumberOfDOF() const
{
    return TNumNodes;
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::InitializeNodalTemperatureVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rVariables.TemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
        rVariables.DtTemperatureVector[i] =
            rGeom[i].FastGetSolutionStepValue(DT_TEMPERATURE);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
double TransientThermalElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    unsigned int PointNumber,
    double detJ)
{
    return rIntegrationPoints[PointNumber].Weight() * detJ;
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateCapacityMatrix(ElementVariables& rVariables) const
{
    KRATOS_TRY

    const double cWater = rVariables.Porosity * rVariables.Saturation *
                          rVariables.WaterDensity * rVariables.WaterHeatCapacity;
    const double cSolid = (1.0 - rVariables.Porosity) *
                          rVariables.SolidDensity * rVariables.SolidHeatCapacity;
    noalias(rVariables.CapacityMatrix) =
        (cWater + cSolid) * outer_prod(rVariables.N, rVariables.N) *
        rVariables.IntegrationCoefficient * rVariables.DtTemperatureCoefficient;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateConductivityMatrix(ElementVariables& rVariables)
{
    KRATOS_TRY

    rVariables.ConstitutiveMatrix = ZeroMatrix(TDim, TDim);
    GeoThermalDispersion2DLaw::CalculateThermalDispersionMatrix(
        rVariables.ConstitutiveMatrix, GetProperties());

    BoundedMatrix<double, TDim, TNumNodes> Temp =
        prod(rVariables.ConstitutiveMatrix, trans(rVariables.GradNT));
    noalias(rVariables.ConductivityMatrix) =
        prod(rVariables.GradNT, Temp) * rVariables.IntegrationCoefficient;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddCapacityVector(
    VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    KRATOS_TRY

    CalculateCapacityVector(rVariables);
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(
        rRightHandSideVector, rVariables.CapacityVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateCapacityVector(ElementVariables& rVariables) const
{
    KRATOS_TRY

    rVariables.CapacityMatrix /= rVariables.DtTemperatureCoefficient;
    noalias(rVariables.CapacityVector) =
        -prod(rVariables.CapacityMatrix, rVariables.DtTemperatureVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddConductivityVector(
    VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    KRATOS_TRY

    CalculateConductivityVector(rVariables);
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(
        rRightHandSideVector, rVariables.ConductivityVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateConductivityVector(ElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rVariables.ConductivityVector) =
        -prod(rVariables.ConductivityMatrix, rVariables.TemperatureVector);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::InitializeProperties(ElementVariables& rVariables)
{
    KRATOS_TRY

    const PropertiesType& rProp = GetProperties();

    rVariables.WaterDensity = rProp[DENSITY_WATER];
    rVariables.SolidDensity = rProp[DENSITY_SOLID];
    rVariables.Porosity = rProp[POROSITY];
    rVariables.WaterHeatCapacity = rProp[SPECIFIC_HEAT_CAPACITY_WATER];
    rVariables.SolidHeatCapacity = rProp[SPECIFIC_HEAT_CAPACITY_SOLID];
    rVariables.Saturation = rProp[SATURATION];

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int N_DOF = GetNumberOfDOF();

    if (rLeftHandSideMatrix.size1() != N_DOF) {
        rLeftHandSideMatrix.resize(N_DOF, N_DOF, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(N_DOF, N_DOF);

    if (rRightHandSideVector.size() != N_DOF) {
        rRightHandSideVector.resize(N_DOF, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(N_DOF);

    constexpr bool CalculateStiffnessMatrixFlag = true;
    constexpr bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod TransientThermalElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    using Data = GeometryData::IntegrationMethod;

    switch (TNumNodes) {
    case 3:
        return Data::GI_GAUSS_2;
    case 6:
        return Data::GI_GAUSS_2;
    case 10:
        return Data::GI_GAUSS_4;
    case 15:
        return Data::GI_GAUSS_5;
    default:
        return Data::GI_GAUSS_2;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::VerifyProperty(Kratos::Variable<double>& rVariable) const
{
    const PropertiesType& rProp = GetProperties();
    if (!rProp.Has(rVariable)) {
        KRATOS_ERROR << rVariable.Name()
                     << " does not exist in the material properties." << std::endl;
    }
    else if (rProp[rVariable] < 0.0) {
        KRATOS_ERROR << rVariable.Name() << " has an invalid value at element"
                     << Id() << "." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CheckDomainSize() const
{
    const GeometryType& rGeom = GetGeometry();
    if (rGeom.DomainSize() < 1.0e-15) {
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << Id() << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CheckSolutionStepsData(
    int rId, Kratos::Variable<double>& rVariable) const
{
    const GeometryType& rGeom = GetGeometry();
    if (rGeom[rId].SolutionStepsDataHas(rVariable)) {
        KRATOS_ERROR << "missing variable " << rVariable.Name() << " on node "
                     << rGeom[rId].Id() << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
}

template class TransientThermalElement<2, 3>;
template class TransientThermalElement<2, 4>;
template class TransientThermalElement<2, 6>;
template class TransientThermalElement<2, 8>;
template class TransientThermalElement<2, 9>;
template class TransientThermalElement<2, 10>;
template class TransientThermalElement<2, 15>;
template class TransientThermalElement<3, 4>;
template class TransientThermalElement<3, 8>;
template class TransientThermalElement<3, 10>;
template class TransientThermalElement<3, 20>;
template class TransientThermalElement<3, 27>;

} // Namespace Kratos
