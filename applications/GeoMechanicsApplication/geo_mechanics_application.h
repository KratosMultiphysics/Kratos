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
//                   Mohamed Nabi
//

#pragma once

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/smart_pointers.h"

// Application includes
#include "geo_mechanics_application_variables.h"

// conditions
#include "custom_conditions/Pw_normal_flux_condition.hpp"
#include "custom_conditions/T_microclimate_flux_condition.h"
#include "custom_conditions/T_normal_flux_condition.h"
#include "custom_conditions/U_Pw_face_load_condition.hpp"
#include "custom_conditions/U_Pw_face_load_interface_condition.hpp"
#include "custom_conditions/U_Pw_force_condition.hpp"
#include "custom_conditions/U_Pw_normal_face_load_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_FIC_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_interface_condition.hpp"
#include "custom_conditions/U_Pw_normal_lysmer_absorbing_condition.hpp"
#include "custom_conditions/axisymmetric_U_Pw_normal_face_load_condition.hpp"
#include "custom_conditions/axisymmetric_line_normal_fluid_flux_2D_diff_order_condition.hpp"
#include "custom_conditions/axisymmetric_line_normal_load_2D_diff_order_condition.hpp"
#include "custom_conditions/line_load_2D_diff_order_condition.hpp"
#include "custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.hpp"
#include "custom_conditions/line_normal_load_2D_diff_order_condition.hpp"
#include "custom_conditions/surface_load_3D_diff_order_condition.hpp"
#include "custom_conditions/surface_normal_fluid_flux_3D_diff_order_condition.hpp"
#include "custom_conditions/surface_normal_load_3D_diff_order_condition.hpp"

// Geometries
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_interface_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_2d_4.h"
#include "geometries/line_2d_5.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/prism_interface_3d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/quadrilateral_interface_2d_4.h"
#include "geometries/quadrilateral_interface_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_15.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"

// elements
#include "custom_elements/U_Pw_small_strain_FIC_element.hpp"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/U_Pw_small_strain_interface_element.hpp"
#include "custom_elements/U_Pw_small_strain_link_interface_element.hpp"
#include "custom_elements/U_Pw_updated_lagrangian_FIC_element.hpp"
#include "custom_elements/U_Pw_updated_lagrangian_element.hpp"
#include "custom_elements/drained_U_Pw_small_strain_element.hpp"
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_elements/steady_state_Pw_element.hpp"
#include "custom_elements/steady_state_Pw_interface_element.hpp"
#include "custom_elements/steady_state_Pw_piping_element.hpp"
#include "custom_elements/transient_Pw_element.hpp"
#include "custom_elements/transient_Pw_interface_element.hpp"
#include "custom_elements/undrained_U_Pw_small_strain_element.hpp"
#include "custom_elements/updated_lagrangian_U_Pw_diff_order_element.hpp"

// Element policies
#include "custom_elements/axisymmetric_stress_state.h"
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/three_dimensional_stress_state.h"

// geo structural element
#include "custom_elements/geo_cable_element.hpp"
#include "custom_elements/geo_cr_beam_element_2D2N.hpp"
#include "custom_elements/geo_cr_beam_element_3D2N.hpp"
#include "custom_elements/geo_cr_beam_element_linear_2D2N.hpp"
#include "custom_elements/geo_cr_beam_element_linear_3D2N.hpp"
#include "custom_elements/geo_curved_beam_element.hpp"
#include "custom_elements/geo_linear_truss_element.hpp"
#include "custom_elements/geo_truss_element.hpp"
#include "custom_elements/transient_thermal_element.h"

// constitutive models
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "custom_constitutive/elastic_isotropic_K0_3d_law.h"
#include "custom_constitutive/linear_elastic_2D_beam_law.h"
#include "custom_constitutive/linear_elastic_2D_interface_law.h"
#include "custom_constitutive/linear_elastic_3D_interface_law.h"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"
#include "custom_constitutive/linear_elastic_plane_strain_K0_law.h"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.h"
#include "custom_constitutive/small_strain_udsm_2D_interface_law.hpp"
#include "custom_constitutive/small_strain_udsm_2D_plane_strain_law.hpp"
#include "custom_constitutive/small_strain_udsm_3D_interface_law.hpp"
#include "custom_constitutive/small_strain_udsm_3D_law.hpp"
#include "custom_constitutive/small_strain_umat_2D_interface_law.hpp"
#include "custom_constitutive/small_strain_umat_2D_plane_strain_law.hpp"
#include "custom_constitutive/small_strain_umat_3D_interface_law.hpp"
#include "custom_constitutive/small_strain_umat_3D_law.hpp"
#include "custom_constitutive/thermal_dispersion_law.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeoMechanicsApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosGeoMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosGeoMechanicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    KratosGeoMechanicsApplication();
    ~KratosGeoMechanicsApplication() override = default;
    KratosGeoMechanicsApplication(const KratosGeoMechanicsApplication&) = delete;
    KratosGeoMechanicsApplication& operator=(const KratosGeoMechanicsApplication&) = delete;
    KratosGeoMechanicsApplication(KratosGeoMechanicsApplication&&) = delete;
    KratosGeoMechanicsApplication& operator=(KratosGeoMechanicsApplication&&) = delete;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosGeoMechanicsApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in KratosGeoMechanicsApplication")
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() )

        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    // elements
    // transient one-phase flow elements:
    const TransientPwElement<2, 3> mTransientPwElement2D3N { 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)) , std::make_unique<PlaneStrainStressState>()};
    const TransientPwElement<2, 4> mTransientPwElement2D4N { 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<PlaneStrainStressState>()};
    const TransientPwElement<2, 6> mTransientPwElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)) , std::make_unique<PlaneStrainStressState>()};
    const TransientPwElement<2, 8> mTransientPwElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<PlaneStrainStressState>()};
    const TransientPwElement<2, 9> mTransientPwElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)) , std::make_unique<PlaneStrainStressState>()};
    const TransientPwElement<2,10> mTransientPwElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<PlaneStrainStressState>() };
    const TransientPwElement<2,15> mTransientPwElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<PlaneStrainStressState>() };
    const TransientPwElement<3, 4> mTransientPwElement3D4N { 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<ThreeDimensionalStressState>()};
    const TransientPwElement<3, 8> mTransientPwElement3D8N { 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<ThreeDimensionalStressState>()};
    const TransientPwElement<3,10> mTransientPwElement3D10N{ 0, Kratos::make_shared< Tetrahedra3D10   <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<ThreeDimensionalStressState>() };
    const TransientPwElement<3,20> mTransientPwElement3D20N{ 0, Kratos::make_shared< Hexahedra3D20    <NodeType> >(Element::GeometryType::PointsArrayType(20)), std::make_unique<ThreeDimensionalStressState>() };
    const TransientPwElement<3,27> mTransientPwElement3D27N{ 0, Kratos::make_shared< Hexahedra3D27    <NodeType> >(Element::GeometryType::PointsArrayType(27)), std::make_unique<ThreeDimensionalStressState>() };

    const TransientPwInterfaceElement<2,4> mTransientPwInterfaceElement2D4N{ 0, Kratos::make_shared< QuadrilateralInterface2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const TransientPwInterfaceElement<3,6> mTransientPwInterfaceElement3D6N{ 0, Kratos::make_shared< PrismInterface3D6         <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<ThreeDimensionalStressState>() };
    const TransientPwInterfaceElement<3,8> mTransientPwInterfaceElement3D8N{ 0, Kratos::make_shared< HexahedraInterface3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    // Steady-State one-phase flow elements:
    const SteadyStatePwElement<2, 3> mSteadyStatePwElement2D3N { 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)) , std::make_unique<PlaneStrainStressState>()};
    const SteadyStatePwElement<2, 4> mSteadyStatePwElement2D4N { 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<PlaneStrainStressState>()};
    const SteadyStatePwElement<2, 6> mSteadyStatePwElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)) , std::make_unique<PlaneStrainStressState>()};
    const SteadyStatePwElement<2, 8> mSteadyStatePwElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<PlaneStrainStressState>()};
    const SteadyStatePwElement<2, 9> mSteadyStatePwElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)) , std::make_unique<PlaneStrainStressState>()};
    const SteadyStatePwElement<2,10> mSteadyStatePwElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<PlaneStrainStressState>() };
    const SteadyStatePwElement<2,15> mSteadyStatePwElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<PlaneStrainStressState>() };
    const SteadyStatePwElement<3, 4> mSteadyStatePwElement3D4N { 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<ThreeDimensionalStressState>()};
    const SteadyStatePwElement<3, 8> mSteadyStatePwElement3D8N { 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<ThreeDimensionalStressState>()};
    const SteadyStatePwElement<3,10> mSteadyStatePwElement3D10N{ 0, Kratos::make_shared< Tetrahedra3D10   <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<ThreeDimensionalStressState>() };
    const SteadyStatePwElement<3,20> mSteadyStatePwElement3D20N{ 0, Kratos::make_shared< Hexahedra3D20    <NodeType> >(Element::GeometryType::PointsArrayType(20)), std::make_unique<ThreeDimensionalStressState>() };
    const SteadyStatePwElement<3,27> mSteadyStatePwElement3D27N{ 0, Kratos::make_shared< Hexahedra3D27    <NodeType> >(Element::GeometryType::PointsArrayType(27)), std::make_unique<ThreeDimensionalStressState>() };

    const SteadyStatePwInterfaceElement<2,4> mSteadyStatePwInterfaceElement2D4N{ 0, Kratos::make_shared< QuadrilateralInterface2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const SteadyStatePwInterfaceElement<3,6> mSteadyStatePwInterfaceElement3D6N{ 0, Kratos::make_shared< PrismInterface3D6         <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<ThreeDimensionalStressState>() };
    const SteadyStatePwInterfaceElement<3,8> mSteadyStatePwInterfaceElement3D8N{ 0, Kratos::make_shared< HexahedraInterface3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    const SteadyStatePwPipingElement<2,4> mSteadyStatePwPipingElement2D4N{ 0, Kratos::make_shared< QuadrilateralInterface2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const SteadyStatePwPipingElement<3,6> mSteadyStatePwPipingElement3D6N{ 0, Kratos::make_shared< PrismInterface3D6         <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<ThreeDimensionalStressState>() };
    const SteadyStatePwPipingElement<3,8> mSteadyStatePwPipingElement3D8N{ 0, Kratos::make_shared< HexahedraInterface3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    // small strain elements:
    const UPwSmallStrainElement<2, 3> mUPwSmallStrainElement2D3N { 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)) , std::make_unique<PlaneStrainStressState>()};
    const UPwSmallStrainElement<2, 4> mUPwSmallStrainElement2D4N { 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<PlaneStrainStressState>()};
    const UPwSmallStrainElement<2, 6> mUPwSmallStrainElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)) , std::make_unique<PlaneStrainStressState>()};
    const UPwSmallStrainElement<2, 8> mUPwSmallStrainElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<PlaneStrainStressState>()};
    const UPwSmallStrainElement<2, 9> mUPwSmallStrainElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)) , std::make_unique<PlaneStrainStressState>()};
    const UPwSmallStrainElement<2,10> mUPwSmallStrainElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<PlaneStrainStressState>() };
    const UPwSmallStrainElement<2,15> mUPwSmallStrainElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<PlaneStrainStressState>() };
    const UPwSmallStrainElement<3, 4> mUPwSmallStrainElement3D4N { 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<ThreeDimensionalStressState>()};
    const UPwSmallStrainElement<3, 8> mUPwSmallStrainElement3D8N { 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<ThreeDimensionalStressState>()};
    const UPwSmallStrainElement<3,10> mUPwSmallStrainElement3D10N{ 0, Kratos::make_shared< Tetrahedra3D10   <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<ThreeDimensionalStressState>() };
    const UPwSmallStrainElement<3,20> mUPwSmallStrainElement3D20N{ 0, Kratos::make_shared< Hexahedra3D20    <NodeType> >(Element::GeometryType::PointsArrayType(20)), std::make_unique<ThreeDimensionalStressState>() };
    const UPwSmallStrainElement<3,27> mUPwSmallStrainElement3D27N{ 0, Kratos::make_shared< Hexahedra3D27    <NodeType> >(Element::GeometryType::PointsArrayType(27)), std::make_unique<ThreeDimensionalStressState>() };

    // small strain drained elements:
    const DrainedUPwSmallStrainElement<2,3> mDrainedUPwSmallStrainElement2D3N{ 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)), std::make_unique<PlaneStrainStressState>() };
    const DrainedUPwSmallStrainElement<2,4> mDrainedUPwSmallStrainElement2D4N{ 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const DrainedUPwSmallStrainElement<3,4> mDrainedUPwSmallStrainElement3D4N{ 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<ThreeDimensionalStressState>() };
    const DrainedUPwSmallStrainElement<3,8> mDrainedUPwSmallStrainElement3D8N{ 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    // small strain undrained elements:
    const UndrainedUPwSmallStrainElement<2,3> mUndrainedUPwSmallStrainElement2D3N{ 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)), std::make_unique<PlaneStrainStressState>() };
    const UndrainedUPwSmallStrainElement<2,4> mUndrainedUPwSmallStrainElement2D4N{ 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const UndrainedUPwSmallStrainElement<3,4> mUndrainedUPwSmallStrainElement3D4N{ 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<ThreeDimensionalStressState>() };
    const UndrainedUPwSmallStrainElement<3,8> mUndrainedUPwSmallStrainElement3D8N{ 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    // FIC elements
    const UPwSmallStrainFICElement<2,3> mUPwSmallStrainFICElement2D3N{ 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)), std::make_unique<PlaneStrainStressState>() };
    const UPwSmallStrainFICElement<2,4> mUPwSmallStrainFICElement2D4N{ 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const UPwSmallStrainFICElement<3,4> mUPwSmallStrainFICElement3D4N{ 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<ThreeDimensionalStressState>() };
    const UPwSmallStrainFICElement<3,8> mUPwSmallStrainFICElement3D8N{ 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    // Small strain different order elements
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<PlaneStrainStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<PlaneStrainStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)), std::make_unique<PlaneStrainStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<PlaneStrainStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<PlaneStrainStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D10N{ 0, Kratos::make_shared< Tetrahedra3D10   <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<ThreeDimensionalStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D20N{ 0, Kratos::make_shared< Hexahedra3D20    <NodeType> >(Element::GeometryType::PointsArrayType(20)), std::make_unique<ThreeDimensionalStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D27N{ 0, Kratos::make_shared< Hexahedra3D27    <NodeType> >(Element::GeometryType::PointsArrayType(27)), std::make_unique<ThreeDimensionalStressState>() };

    // small strain axisymmtric elements:
    const UPwSmallStrainElement<2, 3> mUPwSmallStrainAxisymmetricElement2D3N { 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)) , std::make_unique<AxisymmetricStressState>()};
    const UPwSmallStrainElement<2, 4> mUPwSmallStrainAxisymmetricElement2D4N { 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<AxisymmetricStressState>()};
    const UPwSmallStrainElement<2, 6> mUPwSmallStrainAxisymmetricElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)) , std::make_unique<AxisymmetricStressState>()};
    const UPwSmallStrainElement<2, 8> mUPwSmallStrainAxisymmetricElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<AxisymmetricStressState>()};
    const UPwSmallStrainElement<2, 9> mUPwSmallStrainAxisymmetricElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)) , std::make_unique<AxisymmetricStressState>()};
    const UPwSmallStrainElement<2,10> mUPwSmallStrainAxisymmetricElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<AxisymmetricStressState>() };
    const UPwSmallStrainElement<2,15> mUPwSmallStrainAxisymmetricElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<AxisymmetricStressState>() };

    const UPwSmallStrainFICElement<2,3> mUPwSmallStrainAxisymmetricFICElement2D3N{ 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)), std::make_unique<AxisymmetricStressState>() };
    const UPwSmallStrainFICElement<2,4> mUPwSmallStrainAxisymmetricFICElement2D4N{ 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<AxisymmetricStressState>() };

    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderAxisymmetricElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<AxisymmetricStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderAxisymmetricElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<AxisymmetricStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderAxisymmetricElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)), std::make_unique<AxisymmetricStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderAxisymmetricElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<AxisymmetricStressState>() };
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderAxisymmetricElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<AxisymmetricStressState>() };

    // interface elements
    const UPwSmallStrainInterfaceElement<2,4> mUPwSmallStrainInterfaceElement2D4N{ 0, Kratos::make_shared< QuadrilateralInterface2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const UPwSmallStrainInterfaceElement<3,6> mUPwSmallStrainInterfaceElement3D6N{ 0, Kratos::make_shared< PrismInterface3D6         <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<ThreeDimensionalStressState>() };
    const UPwSmallStrainInterfaceElement<3,8> mUPwSmallStrainInterfaceElement3D8N{ 0, Kratos::make_shared< HexahedraInterface3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    const UPwSmallStrainLinkInterfaceElement<2,4> mUPwSmallStrainLinkInterfaceElement2D4N{ 0, Kratos::make_shared< QuadrilateralInterface2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const UPwSmallStrainLinkInterfaceElement<3,6> mUPwSmallStrainLinkInterfaceElement3D6N{ 0, Kratos::make_shared< PrismInterface3D6         <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<ThreeDimensionalStressState>() };
    const UPwSmallStrainLinkInterfaceElement<3,8> mUPwSmallStrainLinkInterfaceElement3D8N{ 0, Kratos::make_shared< HexahedraInterface3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    // Updated-Lagrangian elements:
    const UPwUpdatedLagrangianElement<2, 3> mUPwUpdatedLagrangianElement2D3N { 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)) , std::make_unique<PlaneStrainStressState>()};
    const UPwUpdatedLagrangianElement<2, 4> mUPwUpdatedLagrangianElement2D4N { 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<PlaneStrainStressState>()};
    const UPwUpdatedLagrangianElement<2, 6> mUPwUpdatedLagrangianElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)) , std::make_unique<PlaneStrainStressState>()};
    const UPwUpdatedLagrangianElement<2, 8> mUPwUpdatedLagrangianElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<PlaneStrainStressState>()};
    const UPwUpdatedLagrangianElement<2, 9> mUPwUpdatedLagrangianElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)) , std::make_unique<PlaneStrainStressState>()};
    const UPwUpdatedLagrangianElement<2,10> mUPwUpdatedLagrangianElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<PlaneStrainStressState>() };
    const UPwUpdatedLagrangianElement<2,15> mUPwUpdatedLagrangianElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<PlaneStrainStressState>() };
    const UPwUpdatedLagrangianElement<3, 4> mUPwUpdatedLagrangianElement3D4N { 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<ThreeDimensionalStressState>()};
    const UPwUpdatedLagrangianElement<3, 8> mUPwUpdatedLagrangianElement3D8N { 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<ThreeDimensionalStressState>()};
    const UPwUpdatedLagrangianElement<3,10> mUPwUpdatedLagrangianElement3D10N{ 0, Kratos::make_shared< Tetrahedra3D10   <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<ThreeDimensionalStressState>() };
    const UPwUpdatedLagrangianElement<3,20> mUPwUpdatedLagrangianElement3D20N{ 0, Kratos::make_shared< Hexahedra3D20    <NodeType> >(Element::GeometryType::PointsArrayType(20)), std::make_unique<ThreeDimensionalStressState>() };
    const UPwUpdatedLagrangianElement<3,27> mUPwUpdatedLagrangianElement3D27N{ 0, Kratos::make_shared< Hexahedra3D27    <NodeType> >(Element::GeometryType::PointsArrayType(27)), std::make_unique<ThreeDimensionalStressState>() };

    const UPwUpdatedLagrangianFICElement<2,3> mUPwUpdatedLagrangianFICElement2D3N{ 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)), std::make_unique<PlaneStrainStressState>() };
    const UPwUpdatedLagrangianFICElement<2,4> mUPwUpdatedLagrangianFICElement2D4N{ 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<PlaneStrainStressState>() };
    const UPwUpdatedLagrangianFICElement<3,4> mUPwUpdatedLagrangianFICElement3D4N{ 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<ThreeDimensionalStressState>() };
    const UPwUpdatedLagrangianFICElement<3,8> mUPwUpdatedLagrangianFICElement3D8N{ 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<ThreeDimensionalStressState>() };

    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<PlaneStrainStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<PlaneStrainStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)), std::make_unique<PlaneStrainStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<PlaneStrainStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<PlaneStrainStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement3D10N{ 0, Kratos::make_shared< Tetrahedra3D10   <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<ThreeDimensionalStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement3D20N{ 0, Kratos::make_shared< Hexahedra3D20    <NodeType> >(Element::GeometryType::PointsArrayType(20)), std::make_unique<ThreeDimensionalStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement3D27N{ 0, Kratos::make_shared< Hexahedra3D27    <NodeType> >(Element::GeometryType::PointsArrayType(27)), std::make_unique<ThreeDimensionalStressState>() };

    // Updated-Lagrangian axisymmetric elements
    const UPwUpdatedLagrangianElement<2, 3> mUPwUpdatedLagrangianAxisymmetricElement2D3N { 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)) , std::make_unique<AxisymmetricStressState>()};
    const UPwUpdatedLagrangianElement<2, 4> mUPwUpdatedLagrangianAxisymmetricElement2D4N { 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)) , std::make_unique<AxisymmetricStressState>()};
    const UPwUpdatedLagrangianElement<2, 6> mUPwUpdatedLagrangianAxisymmetricElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)) , std::make_unique<AxisymmetricStressState>()};
    const UPwUpdatedLagrangianElement<2, 8> mUPwUpdatedLagrangianAxisymmetricElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)) , std::make_unique<AxisymmetricStressState>()};
    const UPwUpdatedLagrangianElement<2, 9> mUPwUpdatedLagrangianAxisymmetricElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)) , std::make_unique<AxisymmetricStressState>()};
    const UPwUpdatedLagrangianElement<2,10> mUPwUpdatedLagrangianAxisymmetricElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<AxisymmetricStressState>() };
    const UPwUpdatedLagrangianElement<2,15> mUPwUpdatedLagrangianAxisymmetricElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<AxisymmetricStressState>() };

    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType(6)), std::make_unique<AxisymmetricStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType(8)), std::make_unique<AxisymmetricStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType(9)), std::make_unique<AxisymmetricStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)), std::make_unique<AxisymmetricStressState>() };
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)), std::make_unique<AxisymmetricStressState>() };

    const UPwUpdatedLagrangianFICElement<2,3> mUPwUpdatedLagrangianAxisymmetricFICElement2D3N{ 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType(3)), std::make_unique<AxisymmetricStressState>() };
    const UPwUpdatedLagrangianFICElement<2,4> mUPwUpdatedLagrangianAxisymmetricFICElement2D4N{ 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType(4)), std::make_unique<AxisymmetricStressState>() };

    // geo structural element
    const GeoCrBeamElement2D2N       mGeoCrBeamElement2D2N      { 0, Kratos::make_shared< Line2D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoCrBeamElement3D2N       mGeoCrBeamElement3D2N      { 0, Kratos::make_shared< Line3D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoCrBeamElementLinear2D2N mGeoCrBeamElementLinear2D2N{ 0, Kratos::make_shared< Line2D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoCrBeamElementLinear3D2N mGeoCrBeamElementLinear3D2N{ 0, Kratos::make_shared< Line3D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoTrussElement<2,2>       mGeoTrussElement2D2N       { 0, Kratos::make_shared< Line2D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoTrussElement<3,2>       mGeoTrussElement3D2N       { 0, Kratos::make_shared< Line3D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoLinearTrussElement<2,2> mGeoLinearTrussElement2D2N { 0, Kratos::make_shared< Line2D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoLinearTrussElement<3,2> mGeoLinearTrussElement3D2N { 0, Kratos::make_shared< Line3D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoCableElement<2,2>       mGeoCableElement2D2N       { 0, Kratos::make_shared< Line2D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoCableElement<3,2>       mGeoCableElement3D2N       { 0, Kratos::make_shared< Line3D2 <NodeType> >(Element::GeometryType::PointsArrayType(2)) };
    const GeoCurvedBeamElement<2,3>  mGeoCurvedBeamElement2D3N  { 0, Kratos::make_shared< Line2D3 <NodeType> >(Element::GeometryType::PointsArrayType(3)) };

    // transient one-phase temperature elements:
    const TransientThermalElement<2, 3> mTransientThermalElement2D3N { 0, Kratos::make_shared< Triangle2D3      <NodeType> >(Element::GeometryType::PointsArrayType( 3)) };
    const TransientThermalElement<2, 6> mTransientThermalElement2D6N { 0, Kratos::make_shared< Triangle2D6      <NodeType> >(Element::GeometryType::PointsArrayType( 6)) };
    const TransientThermalElement<2,10> mTransientThermalElement2D10N{ 0, Kratos::make_shared< Triangle2D10     <NodeType> >(Element::GeometryType::PointsArrayType(10)) };
    const TransientThermalElement<2,15> mTransientThermalElement2D15N{ 0, Kratos::make_shared< Triangle2D15     <NodeType> >(Element::GeometryType::PointsArrayType(15)) };
    const TransientThermalElement<2, 4> mTransientThermalElement2D4N { 0, Kratos::make_shared< Quadrilateral2D4 <NodeType> >(Element::GeometryType::PointsArrayType( 4)) };
    const TransientThermalElement<2, 8> mTransientThermalElement2D8N { 0, Kratos::make_shared< Quadrilateral2D8 <NodeType> >(Element::GeometryType::PointsArrayType( 8)) };
    const TransientThermalElement<2, 9> mTransientThermalElement2D9N { 0, Kratos::make_shared< Quadrilateral2D9 <NodeType> >(Element::GeometryType::PointsArrayType( 9)) };
    const TransientThermalElement<3, 4> mTransientThermalElement3D4N { 0, Kratos::make_shared< Tetrahedra3D4    <NodeType> >(Element::GeometryType::PointsArrayType( 4)) };
    const TransientThermalElement<3,10> mTransientThermalElement3D10N{ 0, Kratos::make_shared< Tetrahedra3D10   <NodeType> >(Element::GeometryType::PointsArrayType(10)) };
    const TransientThermalElement<3, 8> mTransientThermalElement3D8N { 0, Kratos::make_shared< Hexahedra3D8     <NodeType> >(Element::GeometryType::PointsArrayType( 8)) };
    const TransientThermalElement<3,20> mTransientThermalElement3D20N{ 0, Kratos::make_shared< Hexahedra3D20    <NodeType> >(Element::GeometryType::PointsArrayType(20)) };
    const TransientThermalElement<3,27> mTransientThermalElement3D27N{ 0, Kratos::make_shared< Hexahedra3D27    <NodeType> >(Element::GeometryType::PointsArrayType(27)) };

    // conditions
    const UPwForceCondition<2,1> mUPwForceCondition2D1N{ 0, Kratos::make_shared< Point2D <NodeType> >(Condition::GeometryType::PointsArrayType(1)) };
    const UPwForceCondition<3,1> mUPwForceCondition3D1N{ 0, Kratos::make_shared< Point3D <NodeType> >(Condition::GeometryType::PointsArrayType(1)) };

    const UPwFaceLoadCondition<2,2> mUPwFaceLoadCondition2D2N{ 0, Kratos::make_shared< Line2D2          <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const UPwFaceLoadCondition<2,3> mUPwFaceLoadCondition2D3N{ 0, Kratos::make_shared< Line2D3          <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwFaceLoadCondition<2,4> mUPwFaceLoadCondition2D4N{ 0, Kratos::make_shared< Line2D4          <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const UPwFaceLoadCondition<2,5> mUPwFaceLoadCondition2D5N{ 0, Kratos::make_shared< Line2D5          <NodeType> >(Condition::GeometryType::PointsArrayType(5)) };
    const UPwFaceLoadCondition<3,3> mUPwFaceLoadCondition3D3N{ 0, Kratos::make_shared< Triangle3D3      <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwFaceLoadCondition<3,4> mUPwFaceLoadCondition3D4N{ 0, Kratos::make_shared< Quadrilateral3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };

    const UPwNormalFaceLoadCondition<2,2> mUPwNormalFaceLoadCondition2D2N{ 0, Kratos::make_shared< Line2D2          <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const UPwNormalFaceLoadCondition<2,3> mUPwNormalFaceLoadCondition2D3N{ 0, Kratos::make_shared< Line2D3          <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwNormalFaceLoadCondition<2,4> mUPwNormalFaceLoadCondition2D4N{ 0, Kratos::make_shared< Line2D4          <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const UPwNormalFaceLoadCondition<2,5> mUPwNormalFaceLoadCondition2D5N{ 0, Kratos::make_shared< Line2D5          <NodeType> >(Condition::GeometryType::PointsArrayType(5)) };
    const UPwNormalFaceLoadCondition<3,3> mUPwNormalFaceLoadCondition3D3N{ 0, Kratos::make_shared< Triangle3D3      <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwNormalFaceLoadCondition<3,4> mUPwNormalFaceLoadCondition3D4N{ 0, Kratos::make_shared< Quadrilateral3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };

    const UPwNormalFluxCondition<2,2> mUPwNormalFluxCondition2D2N{ 0, Kratos::make_shared< Line2D2          <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const UPwNormalFluxCondition<2,3> mUPwNormalFluxCondition2D3N{ 0, Kratos::make_shared< Line2D3          <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwNormalFluxCondition<2,4> mUPwNormalFluxCondition2D4N{ 0, Kratos::make_shared< Line2D4          <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const UPwNormalFluxCondition<2,5> mUPwNormalFluxCondition2D5N{ 0, Kratos::make_shared< Line2D5          <NodeType> >(Condition::GeometryType::PointsArrayType(5)) };
    const UPwNormalFluxCondition<3,3> mUPwNormalFluxCondition3D3N{ 0, Kratos::make_shared< Triangle3D3      <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwNormalFluxCondition<3,4> mUPwNormalFluxCondition3D4N{ 0, Kratos::make_shared< Quadrilateral3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };

    const PwNormalFluxCondition<2,2> mPwNormalFluxCondition2D2N{ 0, Kratos::make_shared< Line2D2          <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const PwNormalFluxCondition<2,3> mPwNormalFluxCondition2D3N{ 0, Kratos::make_shared< Line2D3          <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const PwNormalFluxCondition<2,4> mPwNormalFluxCondition2D4N{ 0, Kratos::make_shared< Line2D4          <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const PwNormalFluxCondition<2,5> mPwNormalFluxCondition2D5N{ 0, Kratos::make_shared< Line2D5          <NodeType> >(Condition::GeometryType::PointsArrayType(5)) };
    const PwNormalFluxCondition<3,3> mPwNormalFluxCondition3D3N{ 0, Kratos::make_shared< Triangle3D3      <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const PwNormalFluxCondition<3,4> mPwNormalFluxCondition3D4N{ 0, Kratos::make_shared< Quadrilateral3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };

    const UPwFaceLoadInterfaceCondition<2,2> mUPwFaceLoadInterfaceCondition2D2N{ 0, Kratos::make_shared< Line2D2                   <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const UPwFaceLoadInterfaceCondition<3,4> mUPwFaceLoadInterfaceCondition3D4N{ 0, Kratos::make_shared< QuadrilateralInterface3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };

    const UPwNormalFluxInterfaceCondition<2,2> mUPwNormalFluxInterfaceCondition2D2N{ 0, Kratos::make_shared< Line2D2                   <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const UPwNormalFluxInterfaceCondition<3,4> mUPwNormalFluxInterfaceCondition3D4N{ 0, Kratos::make_shared< QuadrilateralInterface3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };

    const UPwNormalFluxFICCondition<2,2> mUPwNormalFluxFICCondition2D2N{ 0, Kratos::make_shared< Line2D2          <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const UPwNormalFluxFICCondition<3,3> mUPwNormalFluxFICCondition3D3N{ 0, Kratos::make_shared< Triangle3D3      <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwNormalFluxFICCondition<3,4> mUPwNormalFluxFICCondition3D4N{ 0, Kratos::make_shared< Quadrilateral3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };

    const LineLoad2DDiffOrderCondition mLineLoadDiffOrderCondition2D3N{ 0, Kratos::make_shared< Line2D3 <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const LineLoad2DDiffOrderCondition mLineLoadDiffOrderCondition2D4N{ 0, Kratos::make_shared< Line2D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const LineLoad2DDiffOrderCondition mLineLoadDiffOrderCondition2D5N{ 0, Kratos::make_shared< Line2D5 <NodeType> >(Condition::GeometryType::PointsArrayType(5)) };

    const LineNormalLoad2DDiffOrderCondition mLineNormalLoadDiffOrderCondition2D3N{ 0, Kratos::make_shared< Line2D3 <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const LineNormalLoad2DDiffOrderCondition mLineNormalLoadDiffOrderCondition2D4N{ 0, Kratos::make_shared< Line2D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const LineNormalLoad2DDiffOrderCondition mLineNormalLoadDiffOrderCondition2D5N{ 0, Kratos::make_shared< Line2D5 <NodeType> >(Condition::GeometryType::PointsArrayType(5)) };

    const LineNormalFluidFlux2DDiffOrderCondition mLineNormalFluidFluxDiffOrderCondition2D3N{ 0, Kratos::make_shared< Line2D3 <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };

    const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D6N{ 0, Kratos::make_shared< Triangle3D6      <NodeType> >(Condition::GeometryType::PointsArrayType(6)) };
    const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D8N{ 0, Kratos::make_shared< Quadrilateral3D8 <NodeType> >(Condition::GeometryType::PointsArrayType(8)) };
    const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D9N{ 0, Kratos::make_shared< Quadrilateral3D9 <NodeType> >(Condition::GeometryType::PointsArrayType(9)) };

    const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D6N{ 0, Kratos::make_shared< Triangle3D6      <NodeType> >(Condition::GeometryType::PointsArrayType(6)) };
    const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D8N{ 0, Kratos::make_shared< Quadrilateral3D8 <NodeType> >(Condition::GeometryType::PointsArrayType(8)) };
    const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D9N{ 0, Kratos::make_shared< Quadrilateral3D9 <NodeType> >(Condition::GeometryType::PointsArrayType(9)) };

    const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D6N{ 0, Kratos::make_shared< Triangle3D6      <NodeType> >(Condition::GeometryType::PointsArrayType(6)) };
    const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D8N{ 0, Kratos::make_shared< Quadrilateral3D8 <NodeType> >(Condition::GeometryType::PointsArrayType(8)) };
    const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D9N{ 0, Kratos::make_shared< Quadrilateral3D9 <NodeType> >(Condition::GeometryType::PointsArrayType(9)) };

    const AxisymmetricUPwNormalFaceLoadCondition<2,2> mAxisymmetricUPwNormalFaceLoadCondition2D2N{ 0, Kratos::make_shared< Line2D2 <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const AxisymmetricUPwNormalFaceLoadCondition<2,3> mAxisymmetricUPwNormalFaceLoadCondition2D3N{ 0, Kratos::make_shared< Line2D3 <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const AxisymmetricUPwNormalFaceLoadCondition<2,4> mAxisymmetricUPwNormalFaceLoadCondition2D4N{ 0, Kratos::make_shared< Line2D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const AxisymmetricUPwNormalFaceLoadCondition<2,5> mAxisymmetricUPwNormalFaceLoadCondition2D5N{ 0, Kratos::make_shared< Line2D5 <NodeType> >(Condition::GeometryType::PointsArrayType(5)) };

    const AxisymmetricLineNormalLoad2DDiffOrderCondition mAxisymmetricLineNormalLoadDiffOrderCondition2D3N{ 0, Kratos::make_shared< Line2D3<NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const AxisymmetricLineNormalLoad2DDiffOrderCondition mAxisymmetricLineNormalLoadDiffOrderCondition2D4N{ 0, Kratos::make_shared< Line2D4<NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const AxisymmetricLineNormalLoad2DDiffOrderCondition mAxisymmetricLineNormalLoadDiffOrderCondition2D5N{ 0, Kratos::make_shared< Line2D5<NodeType> >(Condition::GeometryType::PointsArrayType(5)) };

    const AxisymmetricLineNormalFluidFlux2DDiffOrderCondition mAxisymmetricLineNormalFluidFluxDiffOrderCondition2D3N{ 0, Kratos::make_shared< Line2D3<NodeType> >(Condition::GeometryType::PointsArrayType(3)) };

    const UPwLysmerAbsorbingCondition<2,2> mUPwLysmerAbsorbingCondition2D2N{ 0, Kratos::make_shared< Line2D2          <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const UPwLysmerAbsorbingCondition<2,3> mUPwLysmerAbsorbingCondition2D3N{ 0, Kratos::make_shared< Line2D3          <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwLysmerAbsorbingCondition<3,3> mUPwLysmerAbsorbingCondition3D3N{ 0, Kratos::make_shared< Triangle3D3      <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const UPwLysmerAbsorbingCondition<3,4> mUPwLysmerAbsorbingCondition3D4N{ 0, Kratos::make_shared< Quadrilateral3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };

    const GeoTNormalFluxCondition<2, 2> mGeoTNormalFluxCondition2D2N{ 0, Kratos::make_shared< Line2D2          <NodeType> >(Condition::GeometryType::PointsArrayType(2)) };
    const GeoTNormalFluxCondition<2, 3> mGeoTNormalFluxCondition2D3N{ 0, Kratos::make_shared< Line2D3          <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const GeoTNormalFluxCondition<2, 4> mGeoTNormalFluxCondition2D4N{ 0, Kratos::make_shared< Line2D4          <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const GeoTNormalFluxCondition<2, 5> mGeoTNormalFluxCondition2D5N{ 0, Kratos::make_shared< Line2D5          <NodeType> >(Condition::GeometryType::PointsArrayType(5)) };
    const GeoTNormalFluxCondition<3, 3> mGeoTNormalFluxCondition3D3N{ 0, Kratos::make_shared< Triangle3D3      <NodeType> >(Condition::GeometryType::PointsArrayType(3)) };
    const GeoTNormalFluxCondition<3, 6> mGeoTNormalFluxCondition3D6N{ 0, Kratos::make_shared< Triangle3D6      <NodeType> >(Condition::GeometryType::PointsArrayType(6)) };
    const GeoTNormalFluxCondition<3, 4> mGeoTNormalFluxCondition3D4N{ 0, Kratos::make_shared< Quadrilateral3D4 <NodeType> >(Condition::GeometryType::PointsArrayType(4)) };
    const GeoTNormalFluxCondition<3, 8> mGeoTNormalFluxCondition3D8N{ 0, Kratos::make_shared< Quadrilateral3D8 <NodeType> >(Condition::GeometryType::PointsArrayType(8)) };
    const GeoTNormalFluxCondition<3, 9> mGeoTNormalFluxCondition3D9N{ 0, Kratos::make_shared< Quadrilateral3D9 <NodeType> >(Condition::GeometryType::PointsArrayType(9)) };

    const GeoTMicroClimateFluxCondition<2, 2> mGeoTMicroClimateFluxCondition2D2N{0, Kratos::make_shared<Line2D2         <NodeType>>(Condition::GeometryType::PointsArrayType(2))};
    const GeoTMicroClimateFluxCondition<2, 3> mGeoTMicroClimateFluxCondition2D3N{0, Kratos::make_shared<Line2D3         <NodeType>>(Condition::GeometryType::PointsArrayType(3))};
    const GeoTMicroClimateFluxCondition<2, 4> mGeoTMicroClimateFluxCondition2D4N{0, Kratos::make_shared<Line2D4         <NodeType>>(Condition::GeometryType::PointsArrayType(4))};
    const GeoTMicroClimateFluxCondition<2, 5> mGeoTMicroClimateFluxCondition2D5N{0, Kratos::make_shared<Line2D5         <NodeType>>(Condition::GeometryType::PointsArrayType(5))};
    const GeoTMicroClimateFluxCondition<3, 3> mGeoTMicroClimateFluxCondition3D3N{0, Kratos::make_shared<Triangle3D3     <NodeType>>(Condition::GeometryType::PointsArrayType(3))};
    const GeoTMicroClimateFluxCondition<3, 6> mGeoTMicroClimateFluxCondition3D6N{0, Kratos::make_shared<Triangle3D6     <NodeType>>(Condition::GeometryType::PointsArrayType(6))};
    const GeoTMicroClimateFluxCondition<3, 4> mGeoTMicroClimateFluxCondition3D4N{0, Kratos::make_shared<Quadrilateral3D4<NodeType>>(Condition::GeometryType::PointsArrayType(4))};
    const GeoTMicroClimateFluxCondition<3, 8> mGeoTMicroClimateFluxCondition3D8N{0, Kratos::make_shared<Quadrilateral3D8<NodeType>>(Condition::GeometryType::PointsArrayType(8))};
    const GeoTMicroClimateFluxCondition<3, 9> mGeoTMicroClimateFluxCondition3D9N{0, Kratos::make_shared<Quadrilateral3D9<NodeType>>(Condition::GeometryType::PointsArrayType(9))};

    // constitutive models
    const BilinearCohesive3DLaw             mBilinearCohesive3DLaw;
    const BilinearCohesive2DLaw             mBilinearCohesive2DLaw;
    const LinearPlaneStrainK0Law            mLinearPlaneStrainK0Law;
    const GeoLinearElasticPlaneStrain2DLaw  mLinearElasticPlaneStrain2DLaw;
    const ElasticIsotropicK03DLaw           mElasticIsotropicK03DLaw;
    const GeoLinearElasticPlaneStress2DLaw  mLinearElasticPlaneStress2DLaw;

    const SmallStrainUDSM3DLaw            mSmallStrainUDSM3DLaw{};
    const SmallStrainUDSM2DPlaneStrainLaw mSmallStrainUDSM2DPlaneStrainLaw{};
    const SmallStrainUDSM2DInterfaceLaw   mSmallStrainUDSM2DInterfaceLaw{};
    const SmallStrainUDSM3DInterfaceLaw   mSmallStrainUDSM3DInterfaceLaw{};

    const SmallStrainUMAT3DLaw            mSmallStrainUMAT3DLaw{};
    const SmallStrainUMAT2DPlaneStrainLaw mSmallStrainUMAT2DPlaneStrainLaw{};
    const SmallStrainUMAT2DInterfaceLaw   mSmallStrainUMAT2DInterfaceLaw{};
    const SmallStrainUMAT3DInterfaceLaw   mSmallStrainUMAT3DInterfaceLaw{};

    const LinearElastic2DInterfaceLaw     mLinearElastic2DInterfaceLaw;
    const LinearElastic3DInterfaceLaw     mLinearElastic3DInterfaceLaw;

    const LinearElastic2DBeamLaw          mLinearElastic2DBeamLaw;

    const GeoThermalDispersionLaw mGeoThermalDispersion2DLaw{ConstitutiveLaw::SizeType(2)};
    const GeoThermalDispersionLaw mGeoThermalDispersion3DLaw{ConstitutiveLaw::SizeType(3)};

    ///@}

}; // Class KratosGeoMechanicsApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos
