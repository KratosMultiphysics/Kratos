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
#include "custom_elements/transient_thermal_element.hpp"
#include "custom_constitutive/thermal_dispersion_2D_law.hpp"

namespace Kratos
{

    template<unsigned int TDim, unsigned int TNumNodes>
    TransientThermalElement<TDim, TNumNodes> ::
    TransientThermalElement(IndexType NewId) : Element(NewId) {}

    /// Constructor using an array of nodes
    template<unsigned int TDim, unsigned int TNumNodes>
    TransientThermalElement<TDim, TNumNodes> ::
    TransientThermalElement(IndexType NewId,
        const NodesArrayType& ThisNodes) : Element(NewId, ThisNodes) {}

    /// Constructor using Geometry
    template<unsigned int TDim, unsigned int TNumNodes>
    TransientThermalElement<TDim, TNumNodes> ::
    TransientThermalElement(IndexType NewId,
        GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {}

    /// Constructor using Properties
    template<unsigned int TDim, unsigned int TNumNodes>
    TransientThermalElement<TDim, TNumNodes> ::
    TransientThermalElement(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties) {}

    /// Destructor
    template<unsigned int TDim, unsigned int TNumNodes>
    TransientThermalElement<TDim, TNumNodes> ::
    ~TransientThermalElement() = default;

    // ============================================================================================
	// ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TransientThermalElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TransientThermalElement(NewId, pGeom, pProperties));
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const unsigned int N_DOF = this->GetNumberOfDOF();
        if (rElementalDofList.size() != N_DOF) {
            rElementalDofList.resize(N_DOF);
        }

        const GeometryType& rGeom = this->GetGeometry();
        for (unsigned int i = 0; i < N_DOF; ++i) {
            rElementalDofList[i] = rGeom[i].pGetDof(TEMPERATURE);
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const unsigned int N_DOF = this->GetNumberOfDOF();
        if (rResult.size() != N_DOF) {
            rResult.resize(N_DOF, false);
        }

        const GeometryType& rGeom = this->GetGeometry();
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rResult[index++] = rGeom[i].GetDof(TEMPERATURE).EquationId();
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::Initialize(
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

        mIsInitialised = true;

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    int TransientThermalElement<TDim, TNumNodes>::Check(
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

    	const PropertiesType& rProp = this->GetProperties();
        const GeometryType& rGeom = this->GetGeometry();

        if (rGeom.DomainSize() < 1.0e-15) {
            KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;
        }

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (!rGeom[i].SolutionStepsDataHas(TEMPERATURE)) {
                KRATOS_ERROR << "missing variable TEMPERATURE on node " << rGeom[i].Id() << std::endl;
            }
            if (!rGeom[i].SolutionStepsDataHas(DT_TEMPERATURE)) {
                KRATOS_ERROR << "missing variable DT_TEMPERATURE on node " << rGeom[i].Id() << std::endl;
            }
            if (!rGeom[i].HasDofFor(TEMPERATURE)) {
                KRATOS_ERROR << "missing degree of freedom for TEMPERATURE on node " << rGeom[i].Id() << std::endl;
            }
        }

        // Verify properties
        if (!rProp.Has(DENSITY_WATER) || rProp[DENSITY_WATER] < 0.0) {
            KRATOS_ERROR << "DENSITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(POROSITY) || rProp[POROSITY] < 0.0 || rProp[POROSITY] > 1.0) {
            KRATOS_ERROR << "POROSITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(SATURATION) || rProp[SATURATION] < 0.0 || rProp[SATURATION] > 1.0) {
            KRATOS_ERROR << "SATURATION does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(DENSITY_SOLID) || rProp[DENSITY_SOLID] < 0.0) {
            KRATOS_ERROR << "DENSITY_SOLID does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(SPECIFIC_HEAT_CAPACITY_WATER) || rProp[SPECIFIC_HEAT_CAPACITY_WATER] < 0.0) {
            KRATOS_ERROR << "SPECIFIC_HEAT_CAPACITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(SPECIFIC_HEAT_CAPACITY_SOLID) || rProp[SPECIFIC_HEAT_CAPACITY_SOLID] < 0.0) {
            KRATOS_ERROR << "SPECIFIC_HEAT_CAPACITY_SOLID does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(THERMAL_CONDUCTIVITY_WATER) || rProp[THERMAL_CONDUCTIVITY_WATER] < 0.0) {
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(THERMAL_CONDUCTIVITY_SOLID_XX) || rProp[THERMAL_CONDUCTIVITY_SOLID_XX] < 0.0) {
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_XX does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(THERMAL_CONDUCTIVITY_SOLID_YY) || rProp[THERMAL_CONDUCTIVITY_SOLID_YY] < 0.0) {
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_YY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(THERMAL_CONDUCTIVITY_SOLID_XY) || rProp[THERMAL_CONDUCTIVITY_SOLID_XY] < 0.0) {
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_XY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(LONGITUDINAL_DISPERSIVITY) || rProp[LONGITUDINAL_DISPERSIVITY] < 0.0) {
            KRATOS_ERROR << "LONGITUDINAL_DISPERSIVITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(TRANSVERSE_DISPERSIVITY) || rProp[TRANSVERSE_DISPERSIVITY] < 0.0) {
            KRATOS_ERROR << "TRANSVERSE_DISPERSIVITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }
        if (!rProp.Has(SOLID_COMPRESSIBILITY) || rProp[SOLID_COMPRESSIBILITY] < 0.0) {
            KRATOS_ERROR << "SOLID_COMPRESSIBILITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        }

        if (TDim == 2) {
            auto pos = std::find_if(rGeom.begin(), rGeom.end(), [](const auto& node) {return node.Z() != 0.0;});
            if (pos != rGeom.end()) {
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << pos->Id() << std::endl;
            }
        }

        KRATOS_CATCH("");

        return 0;
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        //Previous definitions
        const GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->GetIntegrationMethod());
        const unsigned int NumGPoints = IntegrationPoints.size();

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        //Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            //Compute GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient =
                this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

            //Contributions to the left hand side
            if (CalculateStiffnessMatrixFlag) {
                this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
            }
            //Contributions to the right hand side
            if (CalculateResidualVectorFlag) {
                this->CalculateAndAddRHS(rRightHandSideVector, Variables);
            }
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::InitializeElementVariables(
        ElementVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

    	//Properties variables
    	this->InitializeProperties(rVariables);

    	//ProcessInfo variables
        rVariables.dtTemperatureCoefficient = rCurrentProcessInfo[DT_TEMPERATURE_COEFFICIENT];

        //Nodal Variables
        this->InitializeNodalTemperatureVariables(rVariables);

        //Variables computed at each GP
        rVariables.N.resize(TNumNodes, false);
        rVariables.GradNT.resize(TNumNodes, TDim, false);

        //General Variables
        const GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

        // shape functions
        rVariables.NContainer.resize(NumGPoints, TNumNodes, false);
        rVariables.NContainer = rGeom.ShapeFunctionsValues(this->GetIntegrationMethod());

        // gradient of shape functions and determinant of Jacobian
        rVariables.detJContainer.resize(NumGPoints, false);

        rGeom.ShapeFunctionsIntegrationPointsGradients(rVariables.DN_DXContainer,
            rVariables.detJContainer,
            this->GetIntegrationMethod());

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddLHS(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        this->CalculateAndAddConductivityMatrix(rLeftHandSideMatrix, rVariables);
        this->CalculateAndAddCapacityMatrix(rLeftHandSideMatrix, rVariables);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddConductivityMatrix(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        this->CalculateConductivityMatrix(rVariables);
        GeoElementUtilities::
            AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, rVariables.conductivityMatrix);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddCapacityMatrix(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        this->CalculateCapacityMatrix(rVariables);
        GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, rVariables.capacityMatrix);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddRHS(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        this->CalculateAndAddCapacityVector(rRightHandSideVector, rVariables);
        this->CalculateAndAddConductivityVector(rRightHandSideVector, rVariables);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateKinematics(
        ElementVariables& rVariables,
        unsigned int PointNumber)
    {
        KRATOS_TRY

    	//Setting the vector of shape functions and the matrix of the shape functions global gradients
    	rVariables.N = row(rVariables.NContainer, PointNumber);
        rVariables.GradNT = rVariables.DN_DXContainer[PointNumber];
        rVariables.detJ = rVariables.detJContainer[PointNumber];

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    unsigned int TransientThermalElement<TDim, TNumNodes>::GetNumberOfDOF() const
    {
        return TNumNodes;
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::InitializeNodalTemperatureVariables(
        ElementVariables& rVariables)
    {
        KRATOS_TRY

    	const GeometryType& rGeom = this->GetGeometry();

        //Nodal Variables
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rVariables.temperatureVector[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
            rVariables.dtTemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(DT_TEMPERATURE);
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    double TransientThermalElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        unsigned int PointNumber,
        double detJ)
    {
        return IntegrationPoints[PointNumber].Weight() * detJ;
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateCapacityMatrix(
        ElementVariables& rVariables) const
    {
        KRATOS_TRY

        const double cWater = rVariables.porosity * rVariables.saturation
            * rVariables.waterDensity * rVariables.waterHeatCapacity;
        const double cSolid = (1.0 - rVariables.porosity) * rVariables.solidDensity
            * rVariables.solidHeatCapacity;
        noalias(rVariables.capacityMatrix) = (cWater + cSolid) * outer_prod(rVariables.N, rVariables.N)
            * rVariables.IntegrationCoefficient
            * rVariables.dtTemperatureCoefficient;

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateConductivityMatrix(
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        rVariables.constitutiveMatrix = ZeroMatrix(TDim, TDim);
        const Properties& rProp = this->GetProperties();
        GeoThermalDispersion2DLaw::CalculateThermalDispersionMatrix(rVariables.constitutiveMatrix, rProp);

        BoundedMatrix<double, TDim, TNumNodes> Temp = prod(rVariables.constitutiveMatrix, trans(rVariables.GradNT));
        noalias(rVariables.conductivityMatrix) = prod(rVariables.GradNT, Temp) * rVariables.IntegrationCoefficient;

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddCapacityVector(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        this->CalculateCapacityVector(rVariables);
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.capacityVector);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateCapacityVector(
        ElementVariables& rVariables) const
    {
        KRATOS_TRY

        rVariables.capacityMatrix /= rVariables.dtTemperatureCoefficient;
        noalias(rVariables.capacityVector) = - prod(rVariables.capacityMatrix, rVariables.dtTemperatureVector);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddConductivityVector(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        this->CalculateConductivityVector(rVariables);
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.conductivityVector);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    void TransientThermalElement<TDim, TNumNodes>::CalculateConductivityVector(
        ElementVariables& rVariables)
    {
        KRATOS_TRY

        noalias(rVariables.conductivityVector) = - prod(rVariables.conductivityMatrix, rVariables.temperatureVector);

        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::InitializeProperties(
        ElementVariables& rVariables)
    {
        KRATOS_TRY

    	const PropertiesType& rProp = this->GetProperties();

        rVariables.waterDensity = rProp[DENSITY_WATER];
        rVariables.solidDensity = rProp[DENSITY_SOLID];
        rVariables.porosity = rProp[POROSITY];
        rVariables.waterHeatCapacity = rProp[SPECIFIC_HEAT_CAPACITY_WATER];
        rVariables.solidHeatCapacity = rProp[SPECIFIC_HEAT_CAPACITY_SOLID];
        rVariables.saturation = rProp[SATURATION];

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

    	const unsigned int N_DOF = this->GetNumberOfDOF();

        //Resetting the LHS
        if (rLeftHandSideMatrix.size1() != N_DOF) {
            rLeftHandSideMatrix.resize(N_DOF, N_DOF, false);
        }
        noalias(rLeftHandSideMatrix) = ZeroMatrix(N_DOF, N_DOF);

        //Resetting the RHS
        if (rRightHandSideVector.size() != N_DOF) {
            rRightHandSideVector.resize(N_DOF, false);
        }
        noalias(rRightHandSideVector) = ZeroVector(N_DOF);

        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll(rLeftHandSideMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            CalculateStiffnessMatrixFlag,
            CalculateResidualVectorFlag);

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    GeometryData::IntegrationMethod TransientThermalElement<TDim, TNumNodes>::GetIntegrationMethod() const
    {
        GeometryData::IntegrationMethod GI_GAUSS;
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

    // ============================================================================================
    // ============================================================================================
    template class TransientThermalElement<2, 3>;
    template class TransientThermalElement<2, 4>;
    template class TransientThermalElement<2, 6>;
    template class TransientThermalElement<2, 8>;
    template class TransientThermalElement<2, 9>;
    template class TransientThermalElement<2,10>;
    template class TransientThermalElement<2,15>;
    template class TransientThermalElement<3, 4>;
    template class TransientThermalElement<3, 8>;
    template class TransientThermalElement<3,10>;
    template class TransientThermalElement<3,20>;
    template class TransientThermalElement<3,27>;

} // Namespace Kratos
