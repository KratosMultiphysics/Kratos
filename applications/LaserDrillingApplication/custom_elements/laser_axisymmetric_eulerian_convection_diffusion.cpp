//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/eulerian_conv_diff.h"
#include "laserdrilling_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "laser_axisymmetric_eulerian_convection_diffusion.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
void LaserAxisymmetricEulerianConvectionDiffusionElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = this->GetGeometry();
    const auto integration_points = r_geom.IntegrationPoints(mMyIntegrationMethod);
    const SizeType n_gauss = integration_points.size();

    if (rOutput.size() != n_gauss) {

        rOutput.resize(n_gauss, false);
    }

    // Initialize output value
    for (IndexType g = 0; g < n_gauss; ++g) {

        rOutput[g] = 0.0;
    }

    if (rVariable == THERMAL_ENERGY) {

        // Initialize thermal_energy value
        double thermal_energy = 0.0;

        // Initialize element data container
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

        // Fill element data container with nodal data
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Calculate kinematics
        Vector det_J_vect;
        ShapeFunctionsGradientsType DN_DX;
        const auto N = r_geom.ShapeFunctionsValues(mMyIntegrationMethod);
        r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mMyIntegrationMethod);

        // Gauss points loop
        array_1d<double,TNumNodes> N_g;
        for (IndexType g = 0; g < n_gauss; ++g) {

            // Get Gauss point data
            noalias(N_g) = row(N, g);

            // Calculate radius and temperature at each Gauss point
            double Radius = 0.0;
            double T = 0.0;
            for (IndexType i = 0; i < TNumNodes; ++i) {
                // Gauss point radius
                Radius += N_g[i] * r_geom[i].Y();
                T += N_g[i] * Variables.phi[i]; // This is the unknown variable which, by default, is the Temperature variable
            }

            // Calculate axisymmetric integration weight
            const double w_g = 2.0 * Globals::Pi * Radius * integration_points[g].Weight() * det_J_vect[g];

            thermal_energy += T * Variables.specific_heat * Variables.density * w_g;
        }

        this->SetValue(THERMAL_ENERGY, thermal_energy);

        // Set elemental thermal_energy 
        for (IndexType g = 0; g < n_gauss; ++g) {

            rOutput[g] = thermal_energy; // This is an elemental variable, and thus it is constant for all GPs
        }

    } else if (rVariable == THERMAL_ENERGY_PER_VOLUME) {
        
        // Initialize thermal_energy per unit volume value
        double thermal_energy = 0.0;
        double decomposed_elemental_volume = 0.0;

        // Initialize element data container
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

        // Fill element data container with nodal data
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Calculate kinematics
        Vector det_J_vect;
        ShapeFunctionsGradientsType DN_DX;
        const auto N = r_geom.ShapeFunctionsValues(mMyIntegrationMethod);
        r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mMyIntegrationMethod);

        // Gauss points loop
        array_1d<double,TNumNodes> N_g;
        for (IndexType g = 0; g < n_gauss; ++g) {

            // Get Gauss point data
            noalias(N_g) = row(N, g);

            // Calculate radius and temperature at each Gauss point
            double Radius = 0.0;
            double T = 0.0;
            for (IndexType i = 0; i < TNumNodes; ++i) {
                // Gauss point radius
                Radius += N_g[i] * r_geom[i].Y();
                T += N_g[i] * Variables.phi[i]; // This is the unknown variable which, by default, is the Temperature variable
            }

            // Calculate axisymmetric integration weight
            const double w_g = 2.0 * Globals::Pi * Radius * integration_points[g].Weight() * det_J_vect[g];

            thermal_energy += T * Variables.specific_heat * Variables.density * w_g;
            decomposed_elemental_volume += w_g;
        }

        // Set elemental thermal_energy 
        for (IndexType g = 0; g < n_gauss; ++g) {

            rOutput[g] = thermal_energy/decomposed_elemental_volume; // This is an elemental variable, and thus it is constant for all GPs
        }
    } else if (rVariable == TEMPERATURE) {

        double temperature = 0.0;

        // Initialize element data container
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables, rCurrentProcessInfo);

        // Fill element data container with nodal data
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Calculate kinematics
        Vector det_J_vect;
        ShapeFunctionsGradientsType DN_DX;
        const auto N = r_geom.ShapeFunctionsValues(mMyIntegrationMethod);
        r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mMyIntegrationMethod);

        // Gauss points loop
        array_1d<double,TNumNodes> N_g;
        for (IndexType g = 0; g < n_gauss; ++g) {

            // Get Gauss point data
            noalias(N_g) = row(N, g);

            // Calculate radius and temperature at each Gauss point
            double Radius = 0.0;
            double T = 0.0;
            for (IndexType i = 0; i < TNumNodes; ++i) {
                // Gauss point radius
                Radius += N_g[i] * r_geom[i].Y();
                T += N_g[i] * Variables.phi[i]; // This is the unknown variable which, by default, is the Temperature variable
            }

            temperature = T;
        }

        this->SetValue(TEMPERATURE, temperature);

        // Set elemental temperature 
        for (IndexType g = 0; g < n_gauss; ++g) {

            rOutput[g] = temperature; // This is an elemental variable, and thus it is constant for all GPs
        }
    } else if (rVariable == THERMAL_DECOMPOSITION) {

        double temperature = 0.0;
        double thermal_energy = 0.0;
        double decomposed_elemental_volume = 0.0;

        // Initialize element data container
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables, rCurrentProcessInfo);

        // Fill element data container with nodal data
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Calculate kinematics
        Vector det_J_vect;
        ShapeFunctionsGradientsType DN_DX;
        const auto N = r_geom.ShapeFunctionsValues(mMyIntegrationMethod);
        r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mMyIntegrationMethod);

        // Gauss points loop
        array_1d<double,TNumNodes> N_g;
        for (IndexType g = 0; g < n_gauss; ++g) {

            // Get Gauss point data
            noalias(N_g) = row(N, g);
            double Radius = 0.0;
            double T = 0.0;
            for (IndexType i = 0; i < TNumNodes; ++i) {
                Radius += N_g[i] * r_geom[i].Y();
                T += N_g[i] * Variables.phi[i]; // This is the unknown variable which, by default, is the Temperature variable
            }
            temperature = T;
            const double w_g = 2.0 * Globals::Pi * Radius * integration_points[g].Weight() * det_J_vect[g];

            thermal_energy += T * Variables.specific_heat * Variables.density * w_g;
            decomposed_elemental_volume += w_g;
        }

        double thermal_energy_per_volume_threshold = rCurrentProcessInfo[THERMAL_ENERGY_PER_VOLUME_THRESHOLD];
        double alpha_threshold = rCurrentProcessInfo[ALPHA_THRESHOLD];
        std::string law = rCurrentProcessInfo[DECOMPOSITION_LAW];

        double old_alpha = this->GetValue(THERMAL_DECOMPOSITION);
        double alpha = ComputeAlpha(law, temperature, old_alpha, rCurrentProcessInfo);

        // Set elemental alpha 
        for (IndexType g = 0; g < n_gauss; ++g) {
            rOutput[g] = alpha; // This is an elemental variable, and thus it is constant for all GPs
        }

        // Here we assume a single Gauss point
        this->SetValue(THERMAL_DECOMPOSITION, alpha);
        double thermal_energy_per_volume = thermal_energy / decomposed_elemental_volume;

        if (thermal_energy_per_volume > thermal_energy_per_volume_threshold) {
            if (alpha > alpha_threshold) {
                this->Set(ACTIVE, false);
            }
        }
    } else if (rVariable == DECOMPOSED_ELEMENTAL_VOLUME) {

        if (this->IsActive()) return;

        double decomposed_elemental_volume = 0.0;

        // Initialize element data container
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables, rCurrentProcessInfo);

        // Fill element data container with nodal data
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Calculate kinematics
        Vector det_J_vect;
        ShapeFunctionsGradientsType DN_DX;
        const auto N = r_geom.ShapeFunctionsValues(mMyIntegrationMethod);
        r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mMyIntegrationMethod);

        // Gauss points loop
        array_1d<double,TNumNodes> N_g;
        for (IndexType g = 0; g < n_gauss; ++g) {
            // Get Gauss point data
            noalias(N_g) = row(N, g);
            double Radius = 0.0;
            for (IndexType i = 0; i < TNumNodes; ++i) {
                Radius += N_g[i] * r_geom[i].Y();
            }
            const double w_g = 2.0 * Globals::Pi * Radius * integration_points[g].Weight() * det_J_vect[g];
            decomposed_elemental_volume += w_g;
        }

        // Set elemental volume 
        for (IndexType g = 0; g < n_gauss; ++g) {
            rOutput[g] = decomposed_elemental_volume; // This is an elemental variable, and thus it is constant for all GPs
        }
    } else if (rVariable == ELEMENTAL_VOLUME) {

        double elemental_volume = 0.0;

        // Initialize element data container
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables, rCurrentProcessInfo);

        // Fill element data container with nodal data
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Calculate kinematics
        Vector det_J_vect;
        ShapeFunctionsGradientsType DN_DX;
        const auto N = r_geom.ShapeFunctionsValues(mMyIntegrationMethod);
        r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mMyIntegrationMethod);

        // Gauss points loop
        array_1d<double,TNumNodes> N_g;
        for (IndexType g = 0; g < n_gauss; ++g) {
            // Get Gauss point data
            noalias(N_g) = row(N, g);
            double Radius = 0.0;
            for (IndexType i = 0; i < TNumNodes; ++i) {
                Radius += N_g[i] * r_geom[i].Y();
            }
            const double w_g = 2.0 * Globals::Pi * Radius * integration_points[g].Weight() * det_J_vect[g];
            elemental_volume += w_g;
        }

        // Set elemental volume 
        for (IndexType g = 0; g < n_gauss; ++g) {
            rOutput[g] = elemental_volume; // This is an elemental variable, and thus it is constant for all GPs
        }
    } else if (rVariable == ENERGY_PER_VOLUME) {

        double elem_energy_per_volume = 0.0;

        // Initialize element data container
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables, rCurrentProcessInfo);

        // Fill element data container with nodal data
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Calculate kinematics
        Vector det_J_vect;
        ShapeFunctionsGradientsType DN_DX;
        const auto N = r_geom.ShapeFunctionsValues(mMyIntegrationMethod);
        r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mMyIntegrationMethod);

        // Gauss points loop
        array_1d<double,TNumNodes> N_g;
        for (IndexType g = 0; g < n_gauss; ++g) {
            // Get Gauss point data
            noalias(N_g) = row(N, g);
            double energy_per_volume = 0.0;
            for (IndexType i = 0; i < TNumNodes; ++i) {
                energy_per_volume += N_g[i] * r_geom[i].GetValue(ENERGY_PER_VOLUME);
            }
            elem_energy_per_volume += energy_per_volume;
        }

        // Set elemental volume 
        for (IndexType g = 0; g < n_gauss; ++g) {
            rOutput[g] = elem_energy_per_volume; // This is an elemental variable, and thus it is constant for all GPs
        }
    } else if (rVariable == ENTHALPY_ENERGY_PER_VOLUME) {

        double elem_energy_per_volume = 0.0;

        // Initialize element data container
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables, rCurrentProcessInfo);

        // Fill element data container with nodal data
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Calculate kinematics
        Vector det_J_vect;
        ShapeFunctionsGradientsType DN_DX;
        const auto N = r_geom.ShapeFunctionsValues(mMyIntegrationMethod);
        r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mMyIntegrationMethod);

        // Gauss points loop
        array_1d<double,TNumNodes> N_g;
        for (IndexType g = 0; g < n_gauss; ++g) {
            // Get Gauss point data
            noalias(N_g) = row(N, g);
            double energy_per_volume = 0.0;
            for (IndexType i = 0; i < TNumNodes; ++i) {
                energy_per_volume += N_g[i] * r_geom[i].GetValue(ENTHALPY_ENERGY_PER_VOLUME);
            }
            elem_energy_per_volume += energy_per_volume;
        }

        // Set elemental volume 
        for (IndexType g = 0; g < n_gauss; ++g) {
            rOutput[g] = elem_energy_per_volume; // This is an elemental variable, and thus it is constant for all GPs
        }
    }
}

double ComputeAlpha(std::string law, double temperature, double old_alpha, const ProcessInfo& rCurrentProcessInfo)
{
    if (law == "Prout-Tompkins") {
        double delta_time = rCurrentProcessInfo[DELTA_TIME];
        double A  = rCurrentProcessInfo[A_COEFFICIENT];
        double Ea = rCurrentProcessInfo[EA_COEFFICIENT];
        double R  = rCurrentProcessInfo[R_GAS_CONSTANT];
        double m  = rCurrentProcessInfo[M_COEFFICIENT];
        double n  = rCurrentProcessInfo[N_COEFFICIENT];
        double new_alpha = old_alpha + A * std::exp(-Ea / (R * temperature)) * std::pow(1.0 - old_alpha, n) * std::pow(old_alpha, m) * delta_time;
        return new_alpha;
    } else if (law == "custom") {
        double time = rCurrentProcessInfo[TIME];
        double K1 = rCurrentProcessInfo[DECOMPOSITION_LAW_CONSTANT_1];
        double K2 = rCurrentProcessInfo[DECOMPOSITION_LAW_CONSTANT_2];
        double T0 = rCurrentProcessInfo[DECOMPOSITION_LAW_REFERENCE_TEMPERATURE];
        return std::exp(-K1 / time) / (1.0 + std::exp(-K2 * (temperature - T0)));
    } else {
        return 0.0;
    }
}

template class LaserAxisymmetricEulerianConvectionDiffusionElement<2,3>;
template class LaserAxisymmetricEulerianConvectionDiffusionElement<2,4>;

}
