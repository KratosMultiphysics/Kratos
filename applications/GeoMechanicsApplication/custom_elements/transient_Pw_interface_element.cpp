// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// Application includes
#include "custom_elements/transient_Pw_interface_element.hpp"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/interface_element_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientPwInterfaceElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                      NodesArrayType const& ThisNodes,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientPwInterfaceElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties,
        this->GetStressStatePolicy().Clone(), this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientPwInterfaceElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                      GeometryType::Pointer pGeom,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientPwInterfaceElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int TransientPwInterfaceElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& r_properties = this->GetProperties();
    const GeometryType&   r_geometry   = this->GetGeometry();

    KRATOS_ERROR_IF(this->Id() < 1)
        << "Element found with Id 0 or negative, element: " << this->Id() << std::endl;

    CheckUtilities::CheckHasNodalSolutionStepData(
        r_geometry, {std::cref(WATER_PRESSURE), std::cref(DT_WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)});
    CheckUtilities::CheckHasDofs(r_geometry, {std::cref(WATER_PRESSURE)});

    const CheckProperties check_properties("material properties at element", r_properties,
                                           this->Id(), CheckProperties::Bounds::AllInclusive);
    check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive)->Check(MINIMUM_JOINT_WIDTH);
    check_properties.Check(TRANSVERSAL_PERMEABILITY);
    check_properties.Check(BULK_MODULUS_FLUID);
    check_properties.Check(DYNAMIC_VISCOSITY);
    check_properties.CheckAvailabilityOnly(BIOT_COEFFICIENT);
    check_properties.Check(DENSITY_WATER);
    check_properties.Check(BULK_MODULUS_SOLID);
    constexpr auto max_value_porosity = 1.0;
    check_properties.Check(POROSITY, max_value_porosity);

    return ierr;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    UPwBaseElement::Initialize(rCurrentProcessInfo);

    // Compute initial gap of the joint
    this->CalculateInitialGap(this->GetGeometry());

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    // Resizing mass matrix
    if (rMassMatrix.size1() != N_DOF) rMassMatrix.resize(N_DOF, N_DOF, false);
    noalias(rMassMatrix) = ZeroMatrix(N_DOF, N_DOF);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo&)
{
    // Intentionally empty
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo&)
{
    // Intentionally empty
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION || rVariable == BISHOP_COEFFICIENT ||
        rVariable == DERIVATIVE_OF_SATURATION || rVariable == RELATIVE_PERMEABILITY) {
        UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    } else {
        // Variables computed on Lobatto points
        const GeometryType& r_geometry = this->GetGeometry();
        const unsigned int  NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);
        std::vector<double> GPValues(NumGPoints);

        for (unsigned int i = 0; i < NumGPoints; ++i) {
            GPValues[i] = 0.0;
        }

        // Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        this->InterpolateOutputDoubles(rValues, GPValues);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == FLUID_FLUX_VECTOR || rVariable == LOCAL_FLUID_FLUX_VECTOR) {
        UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    } else {
        // Variables computed on Lobatto points
        const GeometryType& r_geometry = this->GetGeometry();
        const unsigned int  NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);
        std::vector<array_1d<double, 3>> GPValues(NumGPoints);

        for (unsigned int i = 0; i < NumGPoints; ++i) {
            GPValues[i] = ZeroVector(3);
        }

        // Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        this->template InterpolateOutputValues<array_1d<double, 3>>(rValues, GPValues);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == PERMEABILITY_MATRIX || rVariable == LOCAL_PERMEABILITY_MATRIX) {
        UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    } else {
        // Printed on standard GiD Gauss points
        const unsigned int OutputGPoints =
            this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        for (unsigned int i = 0; i < OutputGPoints; ++i) {
            rValues[i].resize(TDim, TDim, false);
            noalias(rValues[i]) = ZeroMatrix(TDim, TDim);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnLobattoIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == FLUID_FLUX_VECTOR) {
        const PropertiesType& r_properties = this->GetProperties();
        const GeometryType&   r_geometry   = this->GetGeometry();
        const unsigned int NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

        // Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
            r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        GeometryType::JacobiansType JContainer(NumGPoints);
        r_geometry.Jacobian(JContainer, mThisIntegrationMethod);

        // Defining necessary variables
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, r_geometry);

        BoundedMatrix<double, TNumNodes, TDim> GradNpT;

        array_1d<double, TDim> LocalFluidFlux;
        array_1d<double, TDim> GradPressureTerm;
        array_1d<double, TDim> FluidFlux;
        SFGradAuxVariables     SFGradAuxVars;

        // Element variables
        InterfaceElementVariables Variables;
        this->InitializeElementVariables(Variables, r_geometry, r_properties, rCurrentProcessInfo);

        // VG: Perhaps a new parameter to get join width and not minimum joint width
        const double& JointWidth = r_properties[MINIMUM_JOINT_WIDTH];

        RetentionLaw::Parameters RetentionParameters(this->GetProperties());

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            this->template CalculateShapeFunctionsGradients<BoundedMatrix<double, TNumNodes, TDim>>(
                GradNpT, SFGradAuxVars, JContainer[GPoint], RotationMatrix, DN_DeContainer[GPoint],
                NContainer, JointWidth, GPoint);

            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                Variables.BodyAcceleration, NContainer, Variables.VolumeAcceleration, GPoint);

            InterfaceElementUtilities::FillPermeabilityMatrix(
                Variables.LocalPermeabilityMatrix, JointWidth, r_properties[TRANSVERSAL_PERMEABILITY]);

            noalias(GradPressureTerm) = prod(trans(GradNpT), Variables.PressureVector);
            noalias(GradPressureTerm) +=
                PORE_PRESSURE_SIGN_FACTOR * r_properties[DENSITY_WATER] * Variables.BodyAcceleration;

            noalias(LocalFluidFlux) = PORE_PRESSURE_SIGN_FACTOR * Variables.DynamicViscosityInverse *
                                      Variables.RelativePermeability *
                                      prod(Variables.LocalPermeabilityMatrix, GradPressureTerm);

            noalias(FluidFlux) = prod(trans(RotationMatrix), LocalFluidFlux);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], FluidFlux);
        }
    } else if (rVariable == LOCAL_FLUID_FLUX_VECTOR) {
        const PropertiesType& r_properties = this->GetProperties();
        const GeometryType&   r_geometry   = this->GetGeometry();
        const unsigned int NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

        // Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
            r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        GeometryType::JacobiansType JContainer(NumGPoints);
        r_geometry.Jacobian(JContainer, mThisIntegrationMethod);

        // Defining necessary variables
        array_1d<double, TNumNodes> PressureVector;
        for (unsigned int i = 0; i < TNumNodes; ++i)
            PressureVector[i] = r_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE);

        array_1d<double, TNumNodes * TDim> VolumeAcceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(VolumeAcceleration, r_geometry,
                                                                     VOLUME_ACCELERATION);
        array_1d<double, TDim> BodyAcceleration;

        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, r_geometry);

        // VG: Perhaps a new parameter to get join width and not minimum joint width
        const double& JointWidth = r_properties[MINIMUM_JOINT_WIDTH];

        BoundedMatrix<double, TNumNodes, TDim> GradNpT;
        BoundedMatrix<double, TDim, TDim>      LocalPermeabilityMatrix = ZeroMatrix(TDim, TDim);
        const double           DynamicViscosityInverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
        const double&          FluidDensity            = r_properties[DENSITY_WATER];
        array_1d<double, TDim> LocalFluidFlux;
        array_1d<double, TDim> GradPressureTerm;
        SFGradAuxVariables     SFGradAuxVars;

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            this->template CalculateShapeFunctionsGradients<BoundedMatrix<double, TNumNodes, TDim>>(
                GradNpT, SFGradAuxVars, JContainer[GPoint], RotationMatrix, DN_DeContainer[GPoint],
                NContainer, JointWidth, GPoint);

            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                BodyAcceleration, NContainer, VolumeAcceleration, GPoint);

            InterfaceElementUtilities::FillPermeabilityMatrix(
                LocalPermeabilityMatrix, JointWidth, r_properties[TRANSVERSAL_PERMEABILITY]);

            noalias(GradPressureTerm) = prod(trans(GradNpT), PressureVector);
            noalias(GradPressureTerm) += PORE_PRESSURE_SIGN_FACTOR * FluidDensity * BodyAcceleration;

            noalias(LocalFluidFlux) = -DynamicViscosityInverse * prod(LocalPermeabilityMatrix, GradPressureTerm);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], LocalFluidFlux);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnLobattoIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == PERMEABILITY_MATRIX) {
        const GeometryType&   r_geometry   = this->GetGeometry();
        const PropertiesType& r_properties = this->GetProperties();
        const unsigned int NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

        // Defining necessary variables
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, r_geometry);

        // VG: Perhaps a new parameter to get join width and not minimum joint width
        const double&                     JointWidth = r_properties[MINIMUM_JOINT_WIDTH];
        BoundedMatrix<double, TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim, TDim);
        BoundedMatrix<double, TDim, TDim> PermeabilityMatrix;

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            InterfaceElementUtilities::FillPermeabilityMatrix(
                LocalPermeabilityMatrix, JointWidth, r_properties[TRANSVERSAL_PERMEABILITY]);

            noalias(PermeabilityMatrix) =
                prod(trans(RotationMatrix),
                     BoundedMatrix<double, TDim, TDim>(prod(LocalPermeabilityMatrix, RotationMatrix)));

            rOutput[GPoint].resize(TDim, TDim, false);
            noalias(rOutput[GPoint]) = PermeabilityMatrix;
        }
    } else if (rVariable == LOCAL_PERMEABILITY_MATRIX) {
        const GeometryType&   r_geometry   = this->GetGeometry();
        const PropertiesType& r_properties = this->GetProperties();

        // Defining the shape functions container
        const unsigned int NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

        // VG: Perhaps a new parameter to get join width and not minimum joint width
        const double& JointWidth              = r_properties[MINIMUM_JOINT_WIDTH];
        const double& TransversalPermeability = r_properties[TRANSVERSAL_PERMEABILITY];
        BoundedMatrix<double, TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim, TDim);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            InterfaceElementUtilities::FillPermeabilityMatrix(LocalPermeabilityMatrix, JointWidth,
                                                              TransversalPermeability);

            rOutput[GPoint].resize(TDim, TDim, false);
            noalias(rOutput[GPoint]) = LocalPermeabilityMatrix;
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                                VectorType& rRightHandSideVector,
                                                                const ProcessInfo& CurrentProcessInfo,
                                                                bool CalculateStiffnessMatrixFlag,
                                                                bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           r_properties = this->GetProperties();
    const GeometryType&                             r_geometry   = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(mThisIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    r_geometry.Jacobian(JContainer, mThisIntegrationMethod);
    Vector detJContainer(NumGPoints);
    r_geometry.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    // Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables, r_geometry, r_properties, CurrentProcessInfo);

    // VG: Perhaps a new parameter to get join width and not minimum joint width
    Variables.JointWidth = r_properties[MINIMUM_JOINT_WIDTH];

    // Auxiliary variables
    array_1d<double, TDim> RelDispVector;
    SFGradAuxVariables     SFGradAuxVars;

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    const bool hasBiotCoefficient = r_properties.Has(BIOT_COEFFICIENT);

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, detJContainer);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer, GPoint);

        this->template CalculateShapeFunctionsGradients<Matrix>(
            Variables.GradNpT, SFGradAuxVars, JContainer[GPoint], Variables.RotationMatrix,
            DN_DeContainer[GPoint], NContainer, Variables.JointWidth, GPoint);

        // Compute BodyAcceleration and Permeability Matrix
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, NContainer, Variables.VolumeAcceleration, GPoint);

        InterfaceElementUtilities::FillPermeabilityMatrix(
            Variables.LocalPermeabilityMatrix, Variables.JointWidth, r_properties[TRANSVERSAL_PERMEABILITY]);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::InitializeElementVariables(InterfaceElementVariables& rVariables,
                                                                              const GeometryType& rGeometry,
                                                                              const PropertiesType& rProperties,
                                                                              const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    // Properties variables
    rVariables.IgnoreUndrained = false; // by inheritance? does not have a meaning for a Pw element

    rVariables.DynamicViscosityInverse = 1.0 / rProperties[DYNAMIC_VISCOSITY];

    // ProcessInfo variables
    rVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    // Nodal Variables
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rVariables.PressureVector[i]   = rGeometry[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.DtPressureVector[i] = rGeometry[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }

    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration,
                                                                 rGeometry, VOLUME_ACCELERATION);

    // General Variables
    this->CalculateRotationMatrix(rVariables.RotationMatrix, rGeometry);

    // Variables computed at each GP
    rVariables.Np.resize(TNumNodes, false);
    rVariables.GradNpT.resize(TNumNodes, TDim, false);

    // Auxiliary variables
    noalias(rVariables.LocalPermeabilityMatrix) = ZeroMatrix(TDim, TDim);

    // Retention law
    rVariables.FluidPressure          = 0.0;
    rVariables.DegreeOfSaturation     = 1.0;
    rVariables.DerivativeOfSaturation = 0.0;
    rVariables.RelativePermeability   = 1.0;
    rVariables.BishopCoefficient      = 1.0;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                      InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(
    MatrixType& rLeftHandSideMatrix, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    rLeftHandSideMatrix += compressibility_matrix * rVariables.DtPressureCoefficient * rVariables.JointWidth;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddPermeabilityMatrix(
    MatrixType& rLeftHandSideMatrix, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const auto permeability_matrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
        rVariables.GradNpT, rVariables.DynamicViscosityInverse, rVariables.LocalPermeabilityMatrix,
        rVariables.RelativePermeability * rVariables.JointWidth, rVariables.IntegrationCoefficient);

    rLeftHandSideMatrix += permeability_matrix;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                      InterfaceElementVariables& rVariables,
                                                                      unsigned int GPoint)
{
    KRATOS_TRY

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddCompressibilityFlow(
    VectorType& rRightHandSideVector, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    const array_1d<double, TNumNodes> compressibility_flow =
        -rVariables.JointWidth * prod(compressibility_matrix, rVariables.DtPressureVector);

    rRightHandSideVector += compressibility_flow;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddPermeabilityFlow(
    VectorType& rRightHandSideVector, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const BoundedMatrix<double, TNumNodes, TDim> temp_matrix =
        prod(rVariables.GradNpT, rVariables.LocalPermeabilityMatrix);

    const BoundedMatrix<double, TNumNodes, TNumNodes> permeability_matrix =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse *
        rVariables.RelativePermeability * prod(temp_matrix, trans(rVariables.GradNpT)) *
        rVariables.JointWidth * rVariables.IntegrationCoefficient;

    const array_1d<double, TNumNodes> permeability_flow =
        -1.0 * prod(permeability_matrix, rVariables.PressureVector);

    rRightHandSideVector += permeability_flow;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                                const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const auto fluid_density = this->GetProperties()[DENSITY_WATER];
    const BoundedMatrix<double, TNumNodes, TDim> temp_matrix =
        -PORE_PRESSURE_SIGN_FACTOR * prod(rVariables.GradNpT, rVariables.LocalPermeabilityMatrix) *
        rVariables.JointWidth * rVariables.IntegrationCoefficient;

    const array_1d<double, TNumNodes> fluid_body_flow =
        rVariables.DynamicViscosityInverse * fluid_density * rVariables.RelativePermeability *
        prod(temp_matrix, rVariables.BodyAcceleration);

    rRightHandSideVector += fluid_body_flow;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                              const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                                    const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if (rValues.size() != N_DOF) rValues.resize(N_DOF, false);

    // Why are we constructing a zero vector here?
    unsigned int index = 0;
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[index++] = 0.0;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if (rValues.size() != N_DOF) rValues.resize(N_DOF, false);

    // Why are we constructing a zero vector here?
    unsigned int index = 0;
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[index++] = 0.0;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if (rValues.size() != N_DOF) rValues.resize(N_DOF, false);

    // Why are we constructing a zero vector here?
    unsigned int index = 0;
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[index++] = 0.0;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
std::size_t TransientPwInterfaceElement<TDim, TNumNodes>::GetNumberOfDOF() const
{
    return TNumNodes;
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::DofsVectorType TransientPwInterfaceElement<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(this->GetGeometry(), WATER_PRESSURE);
}

template class TransientPwInterfaceElement<2, 4>;
template class TransientPwInterfaceElement<3, 6>;
template class TransientPwInterfaceElement<3, 8>;

} // Namespace Kratos
