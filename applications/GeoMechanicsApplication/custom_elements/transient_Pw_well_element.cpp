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
        if (Prop.Has(DYNAMIC_VISCOSITY) == false || Prop[DYNAMIC_VISCOSITY] < 0.0) {
            KRATOS_ERROR << "DYNAMIC_VISCOSITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(TRANSVERSAL_PERMEABILITY) == false || Prop[TRANSVERSAL_PERMEABILITY] < 0.0) {
            KRATOS_ERROR << "TRANSVERSAL_PERMEABILITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(WELL_LENGTH) == false || Prop[WELL_LENGTH] < 0.0) {
            KRATOS_ERROR << "WELL_LENGTH does not exist in the material properties or has an "
                            "invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(WELL_DIAMETER) == false || Prop[WELL_DIAMETER] < 0.0) {
            KRATOS_ERROR << "WELL_DIAMETER does not exist in the material properties or has an "
                            "invalid value at element" << this->Id() << std::endl;
        }
        if (Prop.Has(WATER_COMPRESSIBILITY) == false || Prop[WATER_COMPRESSIBILITY] < 0.0) {
            KRATOS_ERROR << "WATER_COMPRESSIBILITY does not exist in the material properties or has an "
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





        GeometryType::JacobiansType JContainer(NumGPoints);
        for (unsigned int i = 0; i < NumGPoints; ++i)
            JContainer[i].resize(TDim, 1, false);
        rGeom.Jacobian(JContainer, this->GetIntegrationMethod());







        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            // Compute GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            // Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = 
                this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

            //Variables.IntegrationCoefficient = 
            //    this->CalculateIntegrationCoefficient(JContainer[GPoint], IntegrationPoints[GPoint].Weight());

            // Contributions to the left hand side
            if (CalculateStiffnessMatrixFlag) {
                this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
            }
        }

        GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, Variables.DtPressureCoefficient * Variables.compressibilityMatrix);
        GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix,  Variables.permeabilityMatrix);

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
        rVariables.GradNT.resize(TNumNodes, 1, false);

        // General Variables
        const GeometryType& Geom      = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

        // shape functions
        rVariables.NContainer.resize(NumGPoints, TNumNodes, false);
        rVariables.NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());

        // gradient of shape functions and determinant of Jacobian
        rVariables.detJContainer.resize(NumGPoints, false);

        //Geom.ShapeFunctionsIntegrationPointsGradients(rVariables.DN_DXContainer,
        //    rVariables.detJContainer, this->GetIntegrationMethod());
        //

        Geom.DeterminantOfJacobian(rVariables.detJContainer, this->GetIntegrationMethod());
        rVariables.DN_DXContainer = Geom.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

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
        this->CalculateAndAddBodyForceVector(rRightHandSideVector, rVariables);

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
        GeoElementUtilities::GetNodalVariableVector<1, TNumNodes>(rVariables.VolumeAcceleration,
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

//----------------------------------------------------------------------------------------
    template <unsigned int TDim, unsigned int TNumNodes>
    double TransientPwWellElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
        const Matrix& Jacobian,
        const double& Weight)
    {
        double dx_dxi = Jacobian(0, 0);
        double dy_dxi = Jacobian(1, 0);

        double ds = sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);

        return ds * Weight;
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateCompressibilityMatrix(
        ElementVariables& rVariables)
    {
        KRATOS_TRY;

        this->CalculateCompressibilityFactor(rVariables);

        rVariables.compressibilityMatrix += rVariables.compressibilityFactor * outer_prod(rVariables.N, rVariables.N) *
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

        this->CalculatePermeabilityTensor(rVariables);

        Matrix Temp = -prod(rVariables.GradNT, rVariables.permeabilityTensor);

        rVariables.permeabilityMatrix += rVariables.waterDensity / rVariables.dynamicViscosity *
                                        prod(Temp, trans(rVariables.GradNT)) *
                                        rVariables.IntegrationCoefficient ;

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
    void TransientPwWellElement<TDim, TNumNodes>::CalculateCompressibilityVector(ElementVariables& rVariables)
    {
        KRATOS_TRY

        rVariables.compressibilityVector = -prod(rVariables.compressibilityMatrix, rVariables.DtPressureVector);

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
        rVariables.dynamicViscosity  = rProp[DYNAMIC_VISCOSITY];
        rVariables.wellLength        = rProp[WELL_LENGTH];
        rVariables.wellDiameter      = rProp[WELL_DIAMETER];
        rVariables.waterCompressibility = rProp[WATER_COMPRESSIBILITY];

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateCompressibilityFactor(ElementVariables& rVariables)
    {
        const PropertiesType& rProp = this->GetProperties();
        rVariables.compressibilityFactor =
            1.0 / (rProp[DENSITY_WATER] * rProp[WELL_LENGTH] * 9.81) + rProp[WATER_COMPRESSIBILITY];
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculatePermeabilityTensor(ElementVariables& rVariables)
    {
        KRATOS_TRY;

        const PropertiesType& rProp         = this->GetProperties();

        rVariables.permeabilityTensor(0, 0) = rProp[WELL_DIAMETER] * rProp[WELL_DIAMETER] / 32.0;
        //rVariables.permeabilityTensor(1, 1) = rProp[TRANSVERSAL_PERMEABILITY];
        //rVariables.permeabilityTensor(0, 1) = 0.0;
        //rVariables.permeabilityTensor(1, 0) = 0.0;

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateAndAddBodyForceVector(VectorType& rRightHandSideVector,
                                                                                 ElementVariables& rVariables)
    {
        KRATOS_TRY;

        this->CalculateBodyForceVector(rVariables);
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.bodyForceVector);

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateBodyForceVector(ElementVariables& rVariables)
    {
        KRATOS_TRY;

        const GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->GetIntegrationMethod());
        const unsigned int NumGPoints = IntegrationPoints.size();

        const PropertiesType& rProp = this->GetProperties();

        rVariables.bodyForceVector = ZeroVector(TNumNodes);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            this->CalculateKinematics(rVariables, GPoint);
            rVariables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, rVariables.detJ);
            Matrix Temp = -prod(rVariables.GradNT, rVariables.permeabilityTensor);
            rVariables.bodyForceVector += rVariables.waterDensity / rVariables.dynamicViscosity *
                                          prod(Temp, rVariables.waterDensity * rVariables.VolumeAcceleration) *
                                          rVariables.IntegrationCoefficient;
        }
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    void TransientPwWellElement<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                                        VectorType& rRightHandSideVector,
                                                                        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int N_DOF = this->GetNumberOfDOF();

        // Resetting the LHS
        if (rLeftHandSideMatrix.size1() != N_DOF) {
            rLeftHandSideMatrix.resize(N_DOF, N_DOF, false);
        }
        noalias(rLeftHandSideMatrix) = ZeroMatrix(N_DOF, N_DOF);

        // Resetting the RHS
        if (rRightHandSideVector.size() != N_DOF) {
            rRightHandSideVector.resize(N_DOF, false);
        }
        noalias(rRightHandSideVector) = ZeroVector(N_DOF);

        // calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag  = true;

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                     CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template <unsigned int TDim, unsigned int TNumNodes>
    GeometryData::IntegrationMethod TransientPwWellElement<TDim, TNumNodes>::GetIntegrationMethod() const
    {
        GeometryData::IntegrationMethod GI_GAUSS;
        //
        switch (TNumNodes) {
        case 2:
            GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
            break;
        case 3:
            GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
            break;
        case 4:
            GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_3;
            break;
        case 5:
            GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_5;
            break;
        default:
            GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
            break;
        }
        return GI_GAUSS;
    }

    //----------------------------------------------------------------------------------------------------

template class TransientPwWellElement<2, 2>;
template class TransientPwWellElement<2, 3>;
template class TransientPwWellElement<2, 4>;
template class TransientPwWellElement<2, 5>;

} // Namespace Kratos