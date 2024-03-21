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
//

// Application includes
#include "custom_elements/transient_Pw_well_element.hpp"

namespace Kratos
{
    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer TransientPwWellElement<TDim, TNumNodes>::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TransientPwWellElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer TransientPwWellElement<TDim, TNumNodes>::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TransientPwWellElement(NewId, pGeom, pProperties));
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const unsigned int N_DOF = this->GetNumberOfDOF();
        if (rElementalDofList.size() != N_DOF) {
            rElementalDofList.resize(N_DOF);
        }

        const GeometryType& rGeom = this->GetGeometry();
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rElementalDofList[i] = rGeom[i].pGetDof(WATER_PRESSURE);
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const unsigned int N_DOF = this->GetNumberOfDOF();
        if (rResult.size() != N_DOF) {
            rResult.resize(N_DOF, false);
        }

        const GeometryType& rGeom = this->GetGeometry();
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rResult[i] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const PropertiesType& Prop    = this->GetProperties();
        const GeometryType& Geom      = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

        mIsThermalCoupled       = Geom[0].SolutionStepsDataHas(TEMPERATURE);
        mUpdateDensityViscosity = rCurrentProcessInfo[UPDATE_DENSITY_VISCOSITY];

        mIsInitialised = true;

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    int TransientPwWellElement<TDim, TNumNodes>::Check(
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom   = this->GetGeometry();

        if (Geom.DomainSize() < 1.0e-15)
            KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (Geom[i].SolutionStepsDataHas(WATER_PRESSURE) == false) {
                KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;
            }
            if (Geom[i].SolutionStepsDataHas(DT_WATER_PRESSURE) == false) {
                KRATOS_ERROR << "missing variable DT_WATER_PRESSURE on node " << Geom[i].Id() << std::endl;
            }
            if (Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false) {
                KRATOS_ERROR << "missing variable VOLUME_ACCELERATION on node " << Geom[i].Id() << std::endl;
            }
            if (Geom[i].HasDofFor(WATER_PRESSURE) == false) {
                KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;
            }
        }

        // Verify ProcessInfo variables

        // Verify properties
        if (Prop.Has(DENSITY_WATER) == false || Prop[DENSITY_WATER] < 0.0) {
            KRATOS_ERROR << "DENSITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(VISCOSITY_WATER) == false || Prop[VISCOSITY_WATER] < 0.0) {
            KRATOS_ERROR << "VISCOSITY_WATER does not exist in the material properties or has an "
                            "invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(DYNAMIC_VISCOSITY) == false || Prop[DYNAMIC_VISCOSITY] < 0.0) {
            KRATOS_ERROR << "DYNAMIC_VISCOSITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(TRANSVERSAL_PERMEABILITY) == false || Prop[PERMEABILITY_XX] < 0.0) {
            KRATOS_ERROR << "PERMEABILITY_XX does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(WELL_LENGTH) == false || Prop[WELL_LENGTH] < 0.0) {
            KRATOS_ERROR << "WELL_LENGTH does not exist in the material properties or has an "
                            "invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(WELL_DIAMETER) == false || Prop[WELL_DIAMETER] < 0.0) {
            KRATOS_ERROR << "WELL_DIAMETER does not exist in the material properties or has an "
                            "invalid value at element" << this->Id() << std::endl;
        }
        if (TDim == 2) {
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                if (Geom[i].Z() != 0.0)
                    KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
            }
        }

        KRATOS_CATCH("");

        return 0;
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        // Previous definitions
        const GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->GetIntegrationMethod());
        const unsigned int NumGPoints = IntegrationPoints.size();

        // Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            // Compute GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            // Compute weighting coefficient for integration
            Variables.IntegrationCoefficient =
                this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

            // Contributions to the left hand side
            if (CalculateStiffnessMatrixFlag) {
                this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
            }
        }

        GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, Variables.compressibilityMatrix);
        GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, Variables.DtPressureCoefficient * Variables.permeabilityMatrix);

        if (CalculateResidualVectorFlag) {
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::InitializeElementVariables(
        ElementVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Properties variables
        this->InitializeProperties(rVariables);

        // ProcessInfo variables
        rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

        // Nodal Variables
        this->InitializeNodalPorePressureVariables(rVariables);
        this->InitializeNodalVolumeAccelerationVariables(rVariables);

        // Variables computed at each GP
        rVariables.N.resize(TNumNodes, false);
        rVariables.GradNT.resize(TNumNodes, TDim, false);

        // General Variables
        const GeometryType& Geom      = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

        // shape functions
        (rVariables.NContainer).resize(NumGPoints, TNumNodes, false);
        rVariables.NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());

        // gradient of shape functions and determinant of Jacobian
        rVariables.detJContainer.resize(NumGPoints, false);

        Geom.ShapeFunctionsIntegrationPointsGradients(rVariables.DN_DXContainer,
            rVariables.detJContainer, this->GetIntegrationMethod());

        rVariables.compressibilityMatrix = ZeroMatrix(TNumNodes, TNumNodes);
        rVariables.permeabilityMatrix    = ZeroMatrix(TNumNodes, TNumNodes);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateAndAddLHS(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;

        this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);
        this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;

        this->CalculateCompressibilityMatrix(rVariables);

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateAndAddPermeabilityMatrix(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;

        this->CalculatePermeabilityMatrix(rVariables);

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                     ElementVariables& rVariables)
    {
        KRATOS_TRY;

        this->CalculateAndAddCompressibilityVector(rRightHandSideVector, rVariables);
        this->CalculateAndAddPermeabilityVector(rRightHandSideVector, rVariables);
        //this->CalculateAndAddFluidBodyVector(rRightHandSideVector, rVariables);   //TODO: is it rho g ?

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateKinematics(
        ElementVariables& rVariables,
        unsigned int PointNumber)
    {
        KRATOS_TRY

        // Setting the vector of shape functions and the matrix of the shape functions global gradients
        rVariables.N = row(rVariables.NContainer, PointNumber);
        rVariables.GradNT = rVariables.DN_DXContainer[PointNumber];
        rVariables.detJ = rVariables.detJContainer[PointNumber];

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    unsigned int TransientPwWellElement<TDim, TNumNodes>::GetNumberOfDOF() const
    {
        return TNumNodes;
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::InitializeNodalPorePressureVariables(
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        const GeometryType& rGeom = this->GetGeometry();

        // Nodal Variables
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rVariables.pressureVector[i]   = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
            rVariables.DtPressureVector[i] = rGeom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::InitializeNodalVolumeAccelerationVariables(
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        const GeometryType& rGeom = this->GetGeometry();

        // Nodal Variables
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration,
                                                                     rGeom, VOLUME_ACCELERATION);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    double TransientPwWellElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        unsigned int PointNumber,
        double detJ)
    {
        return IntegrationPoints[PointNumber].Weight() * detJ;
    }


    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateCompressibilityMatrix(
        ElementVariables& rVariables)
    {
        KRATOS_TRY;

        rVariables.compressibilityMatrix = -PORE_PRESSURE_SIGN_FACTOR * rVariables.DtPressureCoefficient *
            rVariables.BiotModulusInverse * outer_prod(rVariables.N, rVariables.N) *
            rVariables.IntegrationCoefficient;

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculatePermeabilityMatrix(
        ElementVariables& rVariables)
    {
        KRATOS_TRY;

        Matrix PDimMatrix = -PORE_PRESSURE_SIGN_FACTOR * prod(rVariables.GradNT, rVariables.permeabilityMatrix);

        rVariables.permeabilityMatrix =
            (1.0 / rVariables.dynamicViscosity) *
            prod(PDimMatrix, trans(rVariables.GradNT)) * rVariables.IntegrationCoefficient;

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateAndAddCompressibilityVector(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        this->CalculateCompressibilityVector(rVariables);
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector,
            rVariables.compressibilityVector);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateAndAddPermeabilityVector(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        this->CalculatePermeabilityVector(rVariables);
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector,
            rVariables.permeabilityVector);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateCompressibilityVector(
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        rVariables.compressibilityVector = -prod(rVariables.compressibilityMatrix, rVariables.DtPressureVector);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculatePermeabilityVector(ElementVariables& rVariables)
    {
        KRATOS_TRY

        rVariables.permeabilityVector = -prod(rVariables.permeabilityMatrix, rVariables.pressureVector);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::InitializeProperties(ElementVariables& rVariables)
    {
        KRATOS_TRY

        const PropertiesType& rProp = this->GetProperties();

        rVariables.waterDensity      = rProp[DENSITY_WATER];
        rVariables.dynamicViscosity  = rProp[VISCOSITY_WATER];
        rVariables.wellLength        = rProp[WELL_LENGTH];
        rVariables.wellDiameter      = rProp[WELL_DIAMETER];

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------

template class TransientPwWellElement<2, 2>;
template class TransientPwWellElement<2, 3>;
template class TransientPwWellElement<2, 4>;
template class TransientPwWellElement<2, 5>;

} // Namespace Kratos