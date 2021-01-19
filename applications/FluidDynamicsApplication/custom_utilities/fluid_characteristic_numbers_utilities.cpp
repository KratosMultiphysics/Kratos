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
//


// System includes


// External includes


// Project includes
#include "includes/cfd_variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "fluid_characteristic_numbers_utilities.h"


namespace Kratos
{

    void FluidCharacteristicNumbersUtilities::CalculateLocalCFL(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        // Get the projected element size function according to the corresponding geometry
        // Note that in here it is assumed that all the elements in the model part feature the same geometry
        const auto& r_geom = rModelPart.ElementsBegin()->GetGeometry();
        ElementSizeFunctionType minimum_h_func = FluidCharacteristicNumbersUtilities::GetMinimumElementSizeFunction(r_geom);

        // Calculate the CFL number in each element
        const double current_dt = rModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        block_for_each(rModelPart.Elements(), [&](Element& rElement){
            const double element_cfl = FluidCharacteristicNumbersUtilities::CalculateElementCFL(rElement, minimum_h_func, current_dt);
            rElement.SetValue(CFL_NUMBER, element_cfl);
        });

        KRATOS_CATCH("")
    }

    double FluidCharacteristicNumbersUtilities::CalculateElementCFL(
        const Element &rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt)
    {
        // Calculate the midpoint velocity
        const auto& r_geometry = rElement.GetGeometry();
        const unsigned int n_nodes = r_geometry.PointsNumber();
        array_1d<double,3> element_vel = r_geometry[0].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 1; i < n_nodes; ++i) {
            element_vel += r_geometry[i].FastGetSolutionStepValue(VELOCITY);
        }
        element_vel /= static_cast<double>(n_nodes);

        // Calculate element CFL
        const double h_min = rElementSizeCalculator(r_geometry);
        const double elem_cfl = norm_2(element_vel) * Dt / h_min;

        return elem_cfl;
    }

    template<bool ConsiderArtificialMagnitudes>
    double FluidCharacteristicNumbersUtilities::CalculateElementPrandtlNumber(const Element &rElement)
    {
        // Get diffusivevalues
        const double cp = rElement.GetProperties().GetValue(SPECIFIC_HEAT);
        const auto diff_values = GetDiffusivityValues<ConsiderArtificialMagnitudes>(rElement);

        // Calculate Prandtl number
        const double elem_prandtl = cp * std::get<0>(diff_values) / std::get<1>(diff_values);

        return elem_prandtl;
    }

    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementPecletNumbers(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator)
    {
        // Calculate the midpoint velocity
        const auto& r_geometry = rElement.GetGeometry();
        const unsigned int n_nodes = r_geometry.PointsNumber();
        array_1d<double,3> element_vel = r_geometry[0].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 1; i < n_nodes; ++i) {
            element_vel += r_geometry[i].FastGetSolutionStepValue(VELOCITY);
        }
        element_vel /= static_cast<double>(n_nodes);

        // Calculate midpoint density
        const double rho = AuxiliaryGetDensity<DensityIsNodal>(rElement);

        // Calculate diffusive magnitudes
        const auto diff_values = GetDiffusivityValues<ConsiderArtificialMagnitudes>(rElement);

        // Calculate element average size
        const double h_avg = rElementSizeCalculator(r_geometry);

        // Calculate Peclet numbers
        const double aux = h_avg * norm_2(element_vel) * rho;
        const double cp = rElement.GetProperties().GetValue(SPECIFIC_HEAT);
        const double mu_Pe = aux / std::get<0>(diff_values);
        const double k_Pe =  aux * cp / std::get<1>(diff_values);
        return std::make_pair(mu_Pe, k_Pe);
    }

    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    double FluidCharacteristicNumbersUtilities::CalculateElementViscousPecletNumber(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator)
    {
        // Calculate the midpoint velocity
        const auto& r_geometry = rElement.GetGeometry();
        const unsigned int n_nodes = r_geometry.PointsNumber();
        array_1d<double,3> element_vel = r_geometry[0].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 1; i < n_nodes; ++i) {
            element_vel += r_geometry[i].FastGetSolutionStepValue(VELOCITY);
        }
        element_vel /= static_cast<double>(n_nodes);

        // Calculate midpoint density
        const double rho = AuxiliaryGetDensity<DensityIsNodal>(rElement);

        // Calculate dynamic viscosity
        const double mu = GetDynamicViscosityValue<ConsiderArtificialMagnitudes>(rElement);

        // Calculate element average size
        const double h_avg = rElementSizeCalculator(r_geometry);

        // Calculate Peclet numbers
        const double mu_Pe = h_avg * norm_2(element_vel) * rho / mu;
        return mu_Pe;
    }


    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    double FluidCharacteristicNumbersUtilities::CalculateElementThermalPecletNumber(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator)
    {
        // Calculate the midpoint velocity
        const auto& r_geometry = rElement.GetGeometry();
        const unsigned int n_nodes = r_geometry.PointsNumber();
        array_1d<double,3> element_vel = r_geometry[0].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 1; i < n_nodes; ++i) {
            element_vel += r_geometry[i].FastGetSolutionStepValue(VELOCITY);
        }
        element_vel /= static_cast<double>(n_nodes);

        // Calculate midpoint density
        const double rho = AuxiliaryGetDensity<DensityIsNodal>(rElement);

        // Calculate diffusive magnitudes
        const double k = GetConductivityValue<ConsiderArtificialMagnitudes>(rElement);

        // Calculate element average size
        const double h_avg = rElementSizeCalculator(r_geometry);

        // Calculate Peclet numbers
        const double aux = h_avg * norm_2(element_vel) * rho;
        const double cp = rElement.GetProperties().GetValue(SPECIFIC_HEAT);
        const double k_Pe =  aux * cp / k;
        return k_Pe;
    }

    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt)
    {
        // Calculate midpoint density
        const double rho = AuxiliaryGetDensity<DensityIsNodal>(rElement);

        // Calculate diffusive magnitudes
        const auto diff_values = GetDiffusivityValues<ConsiderArtificialMagnitudes>(rElement);

        // Calculate element average size
        const double h = rElementSizeCalculator(rElement.GetGeometry());

        // Calculate Fourier numbers
        const double aux = Dt / (rho * std::pow(h,2));
        const double cp = rElement.GetProperties().GetValue(SPECIFIC_HEAT);
        const double visc_Fo = aux * std::get<0>(diff_values);
        const double thermal_Fo =  aux * std::get<1>(diff_values) / cp;
        return std::make_pair(visc_Fo, thermal_Fo);
    }

    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    double FluidCharacteristicNumbersUtilities::CalculateElementViscousFourierNumber(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt)
    {
        // Calculate midpoint density
        const double rho = AuxiliaryGetDensity<DensityIsNodal>(rElement);

        // Calculate dynamic viscosity
        const double mu = GetDynamicViscosityValue<ConsiderArtificialMagnitudes>(rElement);

        // Calculate element average size
        const double h = rElementSizeCalculator(rElement.GetGeometry());

        // Calculate viscous Fourier number
        const double visc_Fo = Dt * mu / (rho * std::pow(h,2));
        return visc_Fo;
    }


    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    double FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt)
    {
        // Calculate midpoint density
        const double rho = AuxiliaryGetDensity<DensityIsNodal>(rElement);

        // Calculate diffusive magnitudes
        const double k = GetConductivityValue<ConsiderArtificialMagnitudes>(rElement);

        // Calculate element average size
        const double h = rElementSizeCalculator(rElement.GetGeometry());

        // Calculate thermal Fourier number
        const double cp = rElement.GetProperties().GetValue(SPECIFIC_HEAT);
        const double thermal_Fo =  Dt * k / (rho * cp * std::pow(h,2));
        return thermal_Fo;
    }

    double FluidCharacteristicNumbersUtilities::CalculateElementMachNumber(const Element &rElement)
    {
        // Get midpoint values
        const auto& r_geom = rElement.GetGeometry();
        const unsigned int n_nodes = r_geom.PointsNumber();
        double midpoint_c = r_geom[0].GetValue(SOUND_VELOCITY);
        array_1d<double,3> midpoint_v = r_geom[0].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i_node = 1; i_node < n_nodes; ++i_node) {
            midpoint_c += r_geom[i_node].GetValue(SOUND_VELOCITY);
            midpoint_v += r_geom[i_node].FastGetSolutionStepValue(VELOCITY);
        }
        midpoint_c /= static_cast<double>(n_nodes);
        midpoint_v /= static_cast<double>(n_nodes);

        // Calculate Prandtl number
        const double elem_mach = norm_2(midpoint_v) / midpoint_c;

        return elem_mach;
    }

    typename FluidCharacteristicNumbersUtilities::ElementSizeFunctionType FluidCharacteristicNumbersUtilities::GetAverageElementSizeFunction(const Geometry<Node<3>>& rGeometry)
    {
        ElementSizeFunctionType average_h_func;
        const auto geometry_type = rGeometry.GetGeometryType();
        switch (geometry_type) {
            case GeometryData::Kratos_Triangle2D3:
                average_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<2,3>::AverageElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Quadrilateral2D4:
                average_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<2,4>::AverageElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Tetrahedra3D4:
                average_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<3,4>::AverageElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Quadrilateral3D8:
                average_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<3,8>::AverageElementSize(rGeometry);};
                break;
            default:
                KRATOS_ERROR << "Non supported geometry type." << std::endl;
        }

        return average_h_func;
    }

    typename FluidCharacteristicNumbersUtilities::ElementSizeFunctionType FluidCharacteristicNumbersUtilities::GetMinimumElementSizeFunction(const Geometry<Node<3>>& rGeometry)
    {
        ElementSizeFunctionType min_h_func;
        const auto geometry_type = rGeometry.GetGeometryType();
        switch (geometry_type) {
            case GeometryData::Kratos_Triangle2D3:
                min_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<2,3>::MinimumElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Quadrilateral2D4:
                min_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<2,4>::MinimumElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Tetrahedra3D4:
                min_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<3,4>::MinimumElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Quadrilateral3D8:
                min_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<3,8>::MinimumElementSize(rGeometry);};
                break;
            default:
                KRATOS_ERROR << "Non supported geometry type." << std::endl;
        }

        return min_h_func;
    }

    template<>
    double FluidCharacteristicNumbersUtilities::AuxiliaryGetDensity<false>(const Element& rElement)
    {
        // Retrieve density from properties
        return rElement.GetProperties().GetValue(DENSITY);
    }

    template<>
    double FluidCharacteristicNumbersUtilities::AuxiliaryGetDensity<true>(const Element& rElement)
    {
        // Calculate midpoint nodal density
        const auto& r_geom = rElement.GetGeometry();
        const unsigned int n_nodes = r_geom.PointsNumber();
        double rho = r_geom[0].FastGetSolutionStepValue(DENSITY);
        for (unsigned int i_node = 1; i_node < n_nodes; ++i_node) {
            rho += r_geom[i_node].FastGetSolutionStepValue(DENSITY);
        }
        rho /= static_cast<double>(n_nodes);
        return rho;
    }

    template<>
    std::pair<double,double> FluidCharacteristicNumbersUtilities::GetDiffusivityValues<false>(const Element& rElement)
    {
        // Get fluid properties
        const auto& r_properties = rElement.GetProperties();
        const double k = r_properties.GetValue(CONDUCTIVITY);
        const double mu = r_properties.GetValue(DYNAMIC_VISCOSITY);

        // Return values
        return std::make_pair(mu, k);
    }

    template<>
    std::pair<double,double> FluidCharacteristicNumbersUtilities::GetDiffusivityValues<true>(const Element& rElement)
    {
        // Get fluid properties
        const auto& r_properties = rElement.GetProperties();
        const double k = r_properties.GetValue(CONDUCTIVITY);
        const double mu = r_properties.GetValue(DYNAMIC_VISCOSITY);

        // Get midpoint artificial magnitudes
        const auto& r_geom = rElement.GetGeometry();
        const unsigned int n_nodes = r_geom.PointsNumber();
        double art_k = r_geom[0].GetValue(ARTIFICIAL_CONDUCTIVITY);
        double art_mu = r_geom[0].GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
        for (unsigned int i_node = 1; i_node < n_nodes; ++i_node) {
            art_k += r_geom[i_node].GetValue(ARTIFICIAL_CONDUCTIVITY);
            art_mu += r_geom[i_node].GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
        }
        art_k /= static_cast<double>(n_nodes);
        art_mu /= static_cast<double>(n_nodes);

        // Return values
        return std::make_pair(mu + art_mu, k + art_k);
    }

    template<>
    double FluidCharacteristicNumbersUtilities::GetDynamicViscosityValue<false>(const Element& rElement)
    {
        // Get fluid properties
        return rElement.GetProperties().GetValue(DYNAMIC_VISCOSITY);
    }

    template<>
    double FluidCharacteristicNumbersUtilities::GetDynamicViscosityValue<true>(const Element& rElement)
    {
        // Get fluid properties
        const auto& r_properties = rElement.GetProperties();
        const double mu = r_properties.GetValue(DYNAMIC_VISCOSITY);

        // Get midpoint artificial magnitudes
        const auto& r_geom = rElement.GetGeometry();
        const unsigned int n_nodes = r_geom.PointsNumber();
        double art_mu = r_geom[0].GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
        for (unsigned int i_node = 1; i_node < n_nodes; ++i_node) {
            art_mu += r_geom[i_node].GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
        }
        art_mu /= static_cast<double>(n_nodes);

        // Return values
        return mu + art_mu;
    }

    template<>
    double FluidCharacteristicNumbersUtilities::GetConductivityValue<false>(const Element& rElement)
    {
        // Get fluid properties
        return rElement.GetProperties().GetValue(CONDUCTIVITY);
    }

    template<>
    double FluidCharacteristicNumbersUtilities::GetConductivityValue<true>(const Element& rElement)
    {
        // Get fluid properties
        const auto& r_properties = rElement.GetProperties();
        const double k = r_properties.GetValue(CONDUCTIVITY);

        // Get midpoint artificial magnitudes
        const auto& r_geom = rElement.GetGeometry();
        const unsigned int n_nodes = r_geom.PointsNumber();
        double art_k = r_geom[0].GetValue(ARTIFICIAL_CONDUCTIVITY);
        for (unsigned int i_node = 1; i_node < n_nodes; ++i_node) {
            art_k += r_geom[i_node].GetValue(ARTIFICIAL_CONDUCTIVITY);
        }
        art_k /= static_cast<double>(n_nodes);

        // Return values
        return k + art_k;
    }

// Template instantiation
template double FluidCharacteristicNumbersUtilities::CalculateElementPrandtlNumber<true>(const Element&);
template double FluidCharacteristicNumbersUtilities::CalculateElementPrandtlNumber<false>(const Element&);

template double FluidCharacteristicNumbersUtilities::CalculateElementViscousPecletNumber<true, true>(const Element&, const ElementSizeFunctionType&);
template double FluidCharacteristicNumbersUtilities::CalculateElementViscousPecletNumber<true, false>(const Element&, const ElementSizeFunctionType&);
template double FluidCharacteristicNumbersUtilities::CalculateElementViscousPecletNumber<false, true>(const Element&, const ElementSizeFunctionType&);
template double FluidCharacteristicNumbersUtilities::CalculateElementViscousPecletNumber<false, false>(const Element&, const ElementSizeFunctionType&);

template double FluidCharacteristicNumbersUtilities::CalculateElementThermalPecletNumber<true, true>(const Element&, const ElementSizeFunctionType&);
template double FluidCharacteristicNumbersUtilities::CalculateElementThermalPecletNumber<true, false>(const Element&, const ElementSizeFunctionType&);
template double FluidCharacteristicNumbersUtilities::CalculateElementThermalPecletNumber<false, true>(const Element&, const ElementSizeFunctionType&);
template double FluidCharacteristicNumbersUtilities::CalculateElementThermalPecletNumber<false, false>(const Element&, const ElementSizeFunctionType&);

template std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementPecletNumbers<true, true>(const Element&, const ElementSizeFunctionType&);
template std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementPecletNumbers<true, false>(const Element&, const ElementSizeFunctionType&);
template std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementPecletNumbers<false, true>(const Element&, const ElementSizeFunctionType&);
template std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementPecletNumbers<false, false>(const Element&, const ElementSizeFunctionType&);

template double FluidCharacteristicNumbersUtilities::CalculateElementViscousFourierNumber<true, true>(const Element&, const ElementSizeFunctionType&, const double);
template double FluidCharacteristicNumbersUtilities::CalculateElementViscousFourierNumber<true, false>(const Element&, const ElementSizeFunctionType&, const double);
template double FluidCharacteristicNumbersUtilities::CalculateElementViscousFourierNumber<false, true>(const Element&, const ElementSizeFunctionType&, const double);
template double FluidCharacteristicNumbersUtilities::CalculateElementViscousFourierNumber<false, false>(const Element&, const ElementSizeFunctionType&, const double);

template double FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<true, true>(const Element&, const ElementSizeFunctionType&, const double);
template double FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<true, false>(const Element&, const ElementSizeFunctionType&, const double);
template double FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<false, true>(const Element&, const ElementSizeFunctionType&, const double);
template double FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<false, false>(const Element&, const ElementSizeFunctionType&, const double);

template std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<true, true>(const Element&, const ElementSizeFunctionType&, const double);
template std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<true, false>(const Element&, const ElementSizeFunctionType&, const double);
template std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<false, true>(const Element&, const ElementSizeFunctionType&, const double);
template std::pair<double,double> FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<false, false>(const Element&, const ElementSizeFunctionType&, const double);

} // namespace Kratos.