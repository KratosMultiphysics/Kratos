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
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/interface_element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientPwInterfaceElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                      NodesArrayType const& ThisNodes,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientPwInterfaceElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientPwInterfaceElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                      GeometryType::Pointer pGeom,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientPwInterfaceElement(NewId, pGeom, pProperties,
                                                            this->GetStressStatePolicy().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int TransientPwInterfaceElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType&   Geom = this->GetGeometry();

    KRATOS_ERROR_IF(this->Id() < 1)
        << "Element found with Id 0 or negative, element: " << this->Id() << std::endl;

    // Verify dof variables
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        KRATOS_ERROR_IF_NOT(Geom[i].SolutionStepsDataHas(WATER_PRESSURE))
            << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        KRATOS_ERROR_IF_NOT(Geom[i].SolutionStepsDataHas(DT_WATER_PRESSURE))
            << "missing variable DT_WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        KRATOS_ERROR_IF_NOT(Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
            << "missing variable VOLUME_ACCELERATION on node " << Geom[i].Id() << std::endl;

        KRATOS_ERROR_IF_NOT(Geom[i].HasDofFor(WATER_PRESSURE))
            << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;
    }

    // Verify specific properties
    if (Prop.Has(MINIMUM_JOINT_WIDTH) == false || Prop[MINIMUM_JOINT_WIDTH] <= 0.0)
        KRATOS_ERROR << "MINIMUM_JOINT_WIDTH has Key zero, is not defined or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(TRANSVERSAL_PERMEABILITY) == false || Prop[TRANSVERSAL_PERMEABILITY] < 0.0)
        KRATOS_ERROR << "TRANSVERSAL_PERMEABILITY has Key zero, is not defined "
                        "or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(BULK_MODULUS_FLUID) == false || Prop[BULK_MODULUS_FLUID] <= 0.0)
        KRATOS_ERROR << "BULK_MODULUS_FLUID has Key zero, is not defined or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(DYNAMIC_VISCOSITY) == false || Prop[DYNAMIC_VISCOSITY] <= 0.0)
        KRATOS_ERROR << "DYNAMIC_VISCOSITY has Key zero, is not defined or has "
                        "an invalid value at element"
                     << this->Id() << std::endl;

    KRATOS_ERROR_IF_NOT(Prop.Has(BIOT_COEFFICIENT))
        << "BIOT_COEFFICIENT does not exist in the material properties in "
           "element"
        << this->Id() << std::endl;

    // Verify properties
    if (Prop.Has(DENSITY_WATER) == false || Prop[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(BULK_MODULUS_SOLID) == false || Prop[BULK_MODULUS_SOLID] < 0.0)
        KRATOS_ERROR << "BULK_MODULUS_SOLID does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(POROSITY) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY does not exist in the material properties or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    return ierr;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    UPwBaseElement<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);

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
void TransientPwInterfaceElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < mRetentionLawVector.size(); ++GPoint) {
        // Initialize retention law
        mRetentionLawVector[GPoint]->InitializeSolutionStep(RetentionParameters);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < mRetentionLawVector.size(); ++GPoint) {
        // retention law
        mRetentionLawVector[GPoint]->FinalizeSolutionStep(RetentionParameters);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION || rVariable == BISHOP_COEFFICIENT ||
        rVariable == DERIVATIVE_OF_SATURATION || rVariable == RELATIVE_PERMEABILITY) {
        UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    } else {
        // Variables computed on Lobatto points
        const GeometryType& Geom       = this->GetGeometry();
        const unsigned int  NumGPoints = Geom.IntegrationPointsNumber(mThisIntegrationMethod);
        std::vector<double> GPValues(NumGPoints);

        for (unsigned int i = 0; i < NumGPoints; ++i) {
            GPValues[i] = 0.0;
        }

        // Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        this->InterpolateOutputDoubles(rValues, GPValues);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == FLUID_FLUX_VECTOR || rVariable == LOCAL_FLUID_FLUX_VECTOR) {
        UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    } else {
        // Variables computed on Lobatto points
        const GeometryType& Geom       = this->GetGeometry();
        const unsigned int  NumGPoints = Geom.IntegrationPointsNumber(mThisIntegrationMethod);
        std::vector<array_1d<double, 3>> GPValues(NumGPoints);

        for (unsigned int i = 0; i < NumGPoints; ++i) {
            GPValues[i] = ZeroVector(3);
        }

        // Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        this->template InterpolateOutputValues<array_1d<double, 3>>(rValues, GPValues);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

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
        const PropertiesType& Prop       = this->GetProperties();
        const GeometryType&   Geom       = this->GetGeometry();
        const unsigned int    NumGPoints = Geom.IntegrationPointsNumber(mThisIntegrationMethod);

        // Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = Geom.ShapeFunctionsValues(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
            Geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        GeometryType::JacobiansType JContainer(NumGPoints);
        Geom.Jacobian(JContainer, mThisIntegrationMethod);

        // Defining necessary variables
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, Geom);

        BoundedMatrix<double, TNumNodes, TDim> GradNpT;

        array_1d<double, TDim> LocalFluidFlux;
        array_1d<double, TDim> GradPressureTerm;
        array_1d<double, TDim> FluidFlux;
        SFGradAuxVariables     SFGradAuxVars;

        // Element variables
        InterfaceElementVariables Variables;
        this->InitializeElementVariables(Variables, Geom, Prop, rCurrentProcessInfo);

        // VG: Perhaps a new parameter to get join width and not minimum joint width
        const double& JointWidth = Prop[MINIMUM_JOINT_WIDTH];

        // create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            this->template CalculateShapeFunctionsGradients<BoundedMatrix<double, TNumNodes, TDim>>(
                GradNpT, SFGradAuxVars, JContainer[GPoint], RotationMatrix, DN_DeContainer[GPoint],
                NContainer, JointWidth, GPoint);

            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                Variables.BodyAcceleration, NContainer, Variables.VolumeAcceleration, GPoint);

            InterfaceElementUtilities::FillPermeabilityMatrix(
                Variables.LocalPermeabilityMatrix, JointWidth, Prop[TRANSVERSAL_PERMEABILITY]);

            noalias(GradPressureTerm) = prod(trans(GradNpT), Variables.PressureVector);
            noalias(GradPressureTerm) +=
                PORE_PRESSURE_SIGN_FACTOR * Variables.FluidDensity * Variables.BodyAcceleration;

            noalias(LocalFluidFlux) = PORE_PRESSURE_SIGN_FACTOR * Variables.DynamicViscosityInverse *
                                      Variables.RelativePermeability *
                                      prod(Variables.LocalPermeabilityMatrix, GradPressureTerm);

            noalias(FluidFlux) = prod(trans(RotationMatrix), LocalFluidFlux);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], FluidFlux);
        }
    } else if (rVariable == LOCAL_FLUID_FLUX_VECTOR) {
        const PropertiesType& Prop       = this->GetProperties();
        const GeometryType&   Geom       = this->GetGeometry();
        const unsigned int    NumGPoints = Geom.IntegrationPointsNumber(mThisIntegrationMethod);

        // Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = Geom.ShapeFunctionsValues(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
            Geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        GeometryType::JacobiansType JContainer(NumGPoints);
        Geom.Jacobian(JContainer, mThisIntegrationMethod);

        // Defining necessary variables
        array_1d<double, TNumNodes> PressureVector;
        for (unsigned int i = 0; i < TNumNodes; ++i)
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);

        array_1d<double, TNumNodes * TDim> VolumeAcceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(VolumeAcceleration, Geom, VOLUME_ACCELERATION);
        array_1d<double, TDim> BodyAcceleration;

        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, Geom);

        // VG: Perhaps a new parameter to get join width and not minimum joint width
        const double& JointWidth = Prop[MINIMUM_JOINT_WIDTH];

        BoundedMatrix<double, TNumNodes, TDim> GradNpT;
        BoundedMatrix<double, TDim, TDim>      LocalPermeabilityMatrix = ZeroMatrix(TDim, TDim);
        const double           DynamicViscosityInverse = 1.0 / Prop[DYNAMIC_VISCOSITY];
        const double&          FluidDensity            = Prop[DENSITY_WATER];
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

            InterfaceElementUtilities::FillPermeabilityMatrix(LocalPermeabilityMatrix, JointWidth,
                                                              Prop[TRANSVERSAL_PERMEABILITY]);

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
        const GeometryType&   Geom       = this->GetGeometry();
        const PropertiesType& Prop       = this->GetProperties();
        const unsigned int    NumGPoints = Geom.IntegrationPointsNumber(mThisIntegrationMethod);

        // Defining necessary variables
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, Geom);

        // VG: Perhaps a new parameter to get join width and not minimum joint width
        const double&                     JointWidth              = Prop[MINIMUM_JOINT_WIDTH];
        BoundedMatrix<double, TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim, TDim);
        BoundedMatrix<double, TDim, TDim> PermeabilityMatrix;

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            InterfaceElementUtilities::FillPermeabilityMatrix(LocalPermeabilityMatrix, JointWidth,
                                                              Prop[TRANSVERSAL_PERMEABILITY]);

            noalias(PermeabilityMatrix) =
                prod(trans(RotationMatrix),
                     BoundedMatrix<double, TDim, TDim>(prod(LocalPermeabilityMatrix, RotationMatrix)));

            rOutput[GPoint].resize(TDim, TDim, false);
            noalias(rOutput[GPoint]) = PermeabilityMatrix;
        }
    } else if (rVariable == LOCAL_PERMEABILITY_MATRIX) {
        const GeometryType&   Geom = this->GetGeometry();
        const PropertiesType& Prop = this->GetProperties();

        // Defining the shape functions container
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber(mThisIntegrationMethod);

        // VG: Perhaps a new parameter to get join width and not minimum joint width
        const double&                     JointWidth              = Prop[MINIMUM_JOINT_WIDTH];
        const double&                     TransversalPermeability = Prop[TRANSVERSAL_PERMEABILITY];
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
                                                                const bool CalculateStiffnessMatrixFlag,
                                                                const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           Prop = this->GetProperties();
    const GeometryType&                             Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        Geom.IntegrationPoints(mThisIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(mThisIntegrationMethod);
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        Geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    Geom.Jacobian(JContainer, mThisIntegrationMethod);
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    // Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables, Geom, Prop, CurrentProcessInfo);

    // VG: Perhaps a new parameter to get join width and not minimum joint width
    Variables.JointWidth = Prop[MINIMUM_JOINT_WIDTH];

    // Auxiliary variables
    array_1d<double, TDim> RelDispVector;
    SFGradAuxVariables     SFGradAuxVars;

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), CurrentProcessInfo);

    const bool hasBiotCoefficient = Prop.Has(BIOT_COEFFICIENT);

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
            Variables.LocalPermeabilityMatrix, Variables.JointWidth, Prop[TRANSVERSAL_PERMEABILITY]);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, detJContainer[GPoint]);

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
                                                                              const GeometryType& Geom,
                                                                              const PropertiesType& Prop,
                                                                              const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    // Properties variables
    rVariables.IgnoreUndrained = false; // by inheritance? does not have a meaning for a Pw element

    rVariables.DynamicViscosityInverse = 1.0 / Prop[DYNAMIC_VISCOSITY];
    rVariables.FluidDensity            = Prop[DENSITY_WATER];
    rVariables.SolidDensity            = Prop[DENSITY_SOLID];
    rVariables.Porosity                = Prop[POROSITY];

    // ProcessInfo variables
    rVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    // Nodal Variables
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rVariables.PressureVector[i]   = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.DtPressureVector[i] = Geom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }

    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration,
                                                                 Geom, VOLUME_ACCELERATION);

    // General Variables
    this->CalculateRotationMatrix(rVariables.RotationMatrix, Geom);

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
    KRATOS_TRY;

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(
    MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables)
{
    KRATOS_TRY;

    noalias(rVariables.PPMatrix) =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.DtPressureCoefficient * rVariables.BiotModulusInverse *
        outer_prod(rVariables.Np, rVariables.Np) * rVariables.JointWidth * rVariables.IntegrationCoefficient;

    // Distribute compressibility block matrix into the elemental matrix
    rLeftHandSideMatrix += rVariables.PPMatrix;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                                     InterfaceElementVariables& rVariables)
{
    KRATOS_TRY;

    noalias(rVariables.PDimMatrix) =
        -PORE_PRESSURE_SIGN_FACTOR * prod(rVariables.GradNpT, rVariables.LocalPermeabilityMatrix);

    noalias(rVariables.PPMatrix) = rVariables.DynamicViscosityInverse * rVariables.RelativePermeability *
                                   prod(rVariables.PDimMatrix, trans(rVariables.GradNpT)) *
                                   rVariables.JointWidth * rVariables.IntegrationCoefficient;

    // Distribute permeability block matrix into the elemental matrix
    rLeftHandSideMatrix += rVariables.PPMatrix;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                      InterfaceElementVariables& rVariables,
                                                                      unsigned int GPoint)
{
    KRATOS_TRY;

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddCompressibilityFlow(
    VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables)
{
    KRATOS_TRY;

    noalias(rVariables.PPMatrix) = -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotModulusInverse *
                                   outer_prod(rVariables.Np, rVariables.Np) *
                                   rVariables.JointWidth * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0 * prod(rVariables.PPMatrix, rVariables.DtPressureVector);

    // Distribute compressibility block vector into elemental vector
    rRightHandSideVector += rVariables.PVector;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                                                                   InterfaceElementVariables& rVariables)
{
    KRATOS_TRY;

    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT, rVariables.LocalPermeabilityMatrix);

    noalias(rVariables.PPMatrix) = -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse *
                                   rVariables.RelativePermeability *
                                   prod(rVariables.PDimMatrix, trans(rVariables.GradNpT)) *
                                   rVariables.JointWidth * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0 * prod(rVariables.PPMatrix, rVariables.PressureVector);

    // Distribute permeability block vector into elemental vector
    rRightHandSideVector += rVariables.PVector;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwInterfaceElement<TDim, TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                                InterfaceElementVariables& rVariables)
{
    KRATOS_TRY;

    noalias(rVariables.PDimMatrix) = -PORE_PRESSURE_SIGN_FACTOR *
                                     prod(rVariables.GradNpT, rVariables.LocalPermeabilityMatrix) *
                                     rVariables.JointWidth * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = rVariables.DynamicViscosityInverse * rVariables.FluidDensity *
                                  rVariables.RelativePermeability *
                                  prod(rVariables.PDimMatrix, rVariables.BodyAcceleration);

    // Distribute fluid body flow block vector into elemental vector
    rRightHandSideVector += rVariables.PVector;

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
unsigned int TransientPwInterfaceElement<TDim, TNumNodes>::GetNumberOfDOF() const
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
