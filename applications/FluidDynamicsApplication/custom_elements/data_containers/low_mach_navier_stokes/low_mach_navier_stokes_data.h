//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes


// External includes


// Project includes

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/data_containers/fluid_element_data.h"
#include "utilities/element_size_calculator.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< size_t TDim, size_t TNumNodes >
class LowMachNavierStokesData : public FluidElementData<TDim,TNumNodes, true>
{
public:

    ///@name Type Definitions
    ///@{

    using NodalScalarData = typename FluidElementData<TDim,TNumNodes, true>::NodalScalarData;
    using NodalVectorData = typename FluidElementData<TDim,TNumNodes, true>::NodalVectorData;
    using ShapeFunctionsType = typename FluidElementData<TDim,TNumNodes, true>::ShapeFunctionsType;
    using MatrixRowType = typename FluidElementData<TDim,TNumNodes, true>::MatrixRowType;

    static constexpr std::size_t BlockSize = TDim + 2;

    ///@}
    ///@name Public Members
    ///@{

    NodalVectorData Velocity;
    NodalVectorData VelocityOldStep1;
    NodalVectorData VelocityOldStep2;

    NodalVectorData MeshVelocity;

    NodalVectorData BodyForce;

    NodalScalarData HeatFlux;

    NodalScalarData Pressure;
    NodalScalarData PressureOldStep1;
    NodalScalarData PressureOldStep2;

    NodalScalarData Temperature;
    NodalScalarData TemperatureOldStep1;
    NodalScalarData TemperatureOldStep2;

    NodalScalarData Density;

    double ThermodynamicPressure;

    double ThermodynamicPressureDerivative;

    double Conductivity;

    double SpecificHeat;

    double HeatCapacityRatio;

    double ThermalExpansionCoefficient;

    double DeltaTime; // Time increment

    double DynamicTau; // Dynamic tau considered in ASGS stabilization coefficients

    double bdf0;
    double bdf1;
    double bdf2;

    double ElementSize;

    ///@}
    ///@name Public Operations
    ///@{

    void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) override
    {
        // Base class Initialize manages constitutive law parameters
        FluidElementData<TDim,TNumNodes, true>::Initialize(rElement,rProcessInfo);

        // Fill historical nodal data
        const auto& r_geometry = rElement.GetGeometry();
        this->FillFromHistoricalNodalData(Velocity, VELOCITY, r_geometry);
        this->FillFromHistoricalNodalData(VelocityOldStep1, VELOCITY, r_geometry, 1);
        this->FillFromHistoricalNodalData(VelocityOldStep2, VELOCITY, r_geometry, 2);
        this->FillFromHistoricalNodalData(MeshVelocity, MESH_VELOCITY, r_geometry);
        this->FillFromHistoricalNodalData(HeatFlux, HEAT_FLUX, r_geometry);
        this->FillFromHistoricalNodalData(BodyForce, BODY_FORCE, r_geometry);
        this->FillFromHistoricalNodalData(Pressure, PRESSURE, r_geometry);
        this->FillFromHistoricalNodalData(PressureOldStep1, PRESSURE, r_geometry, 1);
        this->FillFromHistoricalNodalData(PressureOldStep2, PRESSURE, r_geometry, 2);
        this->FillFromHistoricalNodalData(Temperature, TEMPERATURE, r_geometry);
        this->FillFromHistoricalNodalData(TemperatureOldStep1, TEMPERATURE, r_geometry, 1);
        this->FillFromHistoricalNodalData(TemperatureOldStep2, TEMPERATURE, r_geometry, 2);
        this->FillFromHistoricalNodalData(Density, DENSITY, r_geometry);

        // Fill data from properties
        const auto& r_properties = rElement.GetProperties();
        this->FillFromProperties(Conductivity, CONDUCTIVITY, r_properties);
        this->FillFromProperties(SpecificHeat, SPECIFIC_HEAT, r_properties);
        this->FillFromProperties(HeatCapacityRatio, HEAT_CAPACITY_RATIO, r_properties);

        // Fill data from ProcessInfo container
        this->FillFromProcessInfo(DeltaTime, DELTA_TIME, rProcessInfo);
        this->FillFromProcessInfo(DynamicTau, DYNAMIC_TAU, rProcessInfo);
        if (rProcessInfo.Has(PRESSURE)) {
            this->FillFromProcessInfo(ThermodynamicPressure, PRESSURE, rProcessInfo);
            ThermodynamicPressureDerivative = 0.0;
        } else {
            ThermodynamicPressure = 0.0;
            ThermodynamicPressureDerivative = 0.0;
        }

        const auto& r_BDF_vector = rProcessInfo[BDF_COEFFICIENTS];
        bdf0 = r_BDF_vector[0];
        bdf1 = r_BDF_vector[1];
        bdf2 = r_BDF_vector[2];

        // Calculate element characteristic size
        ElementSize = ElementSizeCalculator<TDim,TNumNodes>::MinimumElementSize(r_geometry);
    }

    void UpdateGeometryValues(
        const unsigned int IntegrationPointIndex,
        double NewWeight,
        const MatrixRowType& rN,
        const BoundedMatrix<double, TNumNodes, TDim>& rDN_DX) override
    {
        FluidElementData<TDim,TNumNodes,true>::UpdateGeometryValues(IntegrationPointIndex,NewWeight,rN,rDN_DX);

    }

    static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
    {
        const auto& r_geometry = rElement.GetGeometry();
        for (unsigned int i = 0; i < TNumNodes; i++) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_geometry[i]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_geometry[i]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE, r_geometry[i]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_geometry[i]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE, r_geometry[i]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_geometry[i]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_geometry[i]);
        }

        return 0;
    }

    ///@}

};

///@}

///@}

}