//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_ENERGY_SPLIT_ELEMENTS_PROCESS_H_INCLUDED)
#define KRATOS_ENERGY_SPLIT_ELEMENTS_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/gid_io.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/divide_tetrahedra_3d_4.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"


namespace Kratos
{

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
class EnergyCheckProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EnergyCheckProcess
    KRATOS_CLASS_POINTER_DEFINITION(EnergyCheckProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    // GaussPointErrorProcess(
    //     ModelPart& rModelPart,
    //     unsigned int DomainSize)
    //     : mrModelPart(rModelPart),
    //       mDomainSize(DomainSize)
    // {
    // }

    EnergyCheckProcess(
        ModelPart &rModelPart,
        unsigned int DomainSize,
        const std::string &FileName)
        : mrModelPart(rModelPart),
          mDomainSize(DomainSize),
          mFileName(FileName)
    {
    }

    /// Destructor.
    ~EnergyCheckProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{





    void WritingFile(const double Potential_energy_water_test,const double Total_cinematic_energy_water,const double Total_energy_water_system,const double Potential_energy_air_test,const double Total_cinematic_air_energy,const double Total_energy_air_system,const double total_energy_system ,const double air_zg,const double water_zg,const double interface_area,const double total_area_negative,const double total_area_positive)
    {
        KRATOS_WATCH("I am writting into file energy check")
        const double current_time = mrModelPart.GetProcessInfo()[TIME];
        std::ofstream MyFile(mFileName, std::ios::app);
        if (MyFile.is_open())
        {
            std::cout<< "success" << std::endl;
        }
        else{
            std::cout << "no funciona" << std::endl;
        }
        //std::string output_line = "total_energy_system\t\tSTEP \n";

        std::string output_line = std::to_string(current_time) + "\t";
        output_line += std::to_string(Potential_energy_water_test) + "\t";
        output_line += std::to_string(Total_cinematic_energy_water) + "\t";
        output_line += std::to_string(Total_energy_water_system) + "\t";
        output_line += std::to_string(Potential_energy_air_test) + "\t";
        output_line += std::to_string(Total_cinematic_air_energy) + "\t";
        output_line += std::to_string(Total_energy_air_system) + "\t";
        output_line += std::to_string(total_energy_system) + "\t";
        output_line += std::to_string(air_zg) + "\t";
        output_line += std::to_string(water_zg) + "\t";
        output_line += std::to_string(interface_area)+"\t";
        output_line += std::to_string(total_area_negative)+"\t";
        output_line += std::to_string(total_area_positive)+"\n";

        MyFile << output_line ;

    }


    double PotentialEnergyTest(const double weighted_elevation,
        const double volume,
        const double Density )

    {

        const double gravity =9.81;
        const double elevation_cg=weighted_elevation/volume;
        KRATOS_WATCH(elevation_cg)
        KRATOS_WATCH(volume)
        return (Density * gravity * elevation_cg*volume);

    }

    double PotentialEnergyRuben(
        const double zg,
        const double volume,
        const double density,
        const double gravity)

    {
        return density * volume * gravity * zg;
    }

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
     return "GaussPointErrorProcessAnalytical";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "GaussPointErrorProcessAnalytical";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    void Execute()
    {

        // Variables initialization
        double cut_area_negative = 0.0;
        double cut_area_positive = 0.0;
        double total_area_negative = 0.0;
        double total_area_positive = 0.0;
        double Total_energy_air_splitted = 0.0;
        double Total_energy_water_splitted = 0.0;
        double Total_energy_air = 0.0;
        double Total_energy_water = 0.0;
        double Total_energy_air_entire_element = 0;
        double Total_energy_water_entire_element = 0;
        double Total_cinematic_energy_water = 0;
        double Total_cinematic_energy_air = 0;
        double Total_hydrostatic_energy_water = 0;
        double Total_hydrostatic_energy_air = 0;
        double Total_cinematic_energy_water_splitted = 0;
        double Total_hydrostatic_energy_water_splitted = 0;
        double Total_cinematic_energy_water_not_splitted = 0;
        double Total_hydrostatic_energy_water_not_splitted = 0;
        double Total_cinematic_energy_air_splitted = 0;
        double Total_hydrostatic_energy_air_splitted = 0;
        double Total_cinematic_energy_air_not_splitted = 0;
        double Total_hydrostatic_energy_air_not_splitted = 0;
        double gravity = 9.81;
        double Total_hydrostatic_energy_water_elemental = 0;
        double Total_potential_energy_water_fluid = 0;
        double Total_hydrostatic_energy_air_elemental = 0;
        double Total_potential_energy_air_fluid = 0;
        const double Water_Density = 1e+3;
        const double Air_Density = 1.22;
        double interface_area = 0.0;
        double air_zg = 0.0;
        double water_zg = 0.0;

        if (mDomainSize == 2)
        // Total Energy check for 2D cases[Fluid1_energy and FLuid2_energy]: A loop in each element of the mesh considering splitted elements (air/water elements)
        {
            double total_energy = 0.0;
            double total_potential_energy = 0.0;
            double total_kinematic_energy = 0.0;

            for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem)
            {
                auto &r_geom = it_elem->GetGeometry();
                const auto elem_dist = SetDistanceVector(*it_elem);

                // Split element case
                if (this->IsCut(elem_dist))
                {
                    // Generate the splitting pattern

                    DivideTriangle2D3<Node>::Pointer p_divide_util = Kratos::make_shared<DivideTriangle2D3<Node>>(r_geom, elem_dist);
                    p_divide_util->GenerateDivision();

                    // Generate the modified shape functions util
                    auto p_geom = it_elem->pGetGeometry();
                    ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, elem_dist);

                    Matrix N_neg_side;
                    Vector w_gauss_neg_side;
                    Matrix N_neg_side_interface;
                    Vector w_gauss_interface;

                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side_interface;
                    p_ausas_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(N_neg_side, DN_DX_neg_side, w_gauss_neg_side, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    p_ausas_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(N_neg_side_interface, DN_DX_neg_side_interface, w_gauss_interface, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    // Save the coordinates of all the subdivision Gauss pts.nl
                    const unsigned int n_neg_gauss_pt = w_gauss_neg_side.size();
                    const unsigned int n_interface_side_gauss_pt = w_gauss_interface.size();
                    for (unsigned int i_gauss_interface = 0; i_gauss_interface < n_interface_side_gauss_pt; ++i_gauss_interface)
                    {
                        interface_area += w_gauss_interface[i_gauss_interface];
                        ;
                    }

                    for (unsigned int i_gauss = 0; i_gauss < n_neg_gauss_pt; ++i_gauss)
                    {
                        const double Water_Density = 1e+3;
                        Vector N_gauss = row(N_neg_side, i_gauss);
                        const double w_gauss = w_gauss_neg_side[i_gauss];
                        const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Water_Density);

                        // const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                        water_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                        // const double total_energy_gauss = kinetic_energy + potential_energy;

                        total_area_negative += w_gauss;

                        // Total_energy_water += total_energy_gauss;
                        // Total_energy_water_splitted+= total_energy_gauss;
                        Total_cinematic_energy_water += kinetic_energy;
                    }
                    Matrix N_pos_side;
                    Vector w_gauss_pos_side;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_pos_side;
                    p_ausas_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(N_pos_side, DN_DX_pos_side, w_gauss_pos_side, GeometryData::IntegrationMethod::GI_GAUSS_2);

                    // Save the coordinates of all the subdivision Gauss pts.
                    const unsigned int n_pos_gauss_pt = w_gauss_pos_side.size();

                    for (unsigned int i_gauss = 0; i_gauss < n_pos_gauss_pt; ++i_gauss)
                    {
                        const double Air_Density = 1.22;
                        Vector N_gauss = row(N_pos_side, i_gauss);
                        const double w_gauss = w_gauss_pos_side[i_gauss];
                        const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                        const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                        const double total_energy_gauss = kinetic_energy + potential_energy;
                        air_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);

                        total_area_positive += w_gauss;

                        Total_energy_air += total_energy_gauss;
                        Total_energy_air_splitted += total_energy_gauss;
                        Total_cinematic_energy_air += kinetic_energy;

                        // Total_hydrostatic_energy_air+=potential_energy;
                    }
                }
                // No split element case
                else
                {
                    // Get geometry data
                    Vector jac_vect;
                    r_geom.DeterminantOfJacobian(jac_vect, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    auto N = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);
                    const auto &r_integrations_points = r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
                    const unsigned int n_gauss = r_integrations_points.size();
                    // Check if the element is in the positive (fluid) region
                    if (this->IsPositive(elem_dist))
                    {
                        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
                        {
                            const double Air_Density = 1.22;
                            Vector N_gauss = row(N, i_gauss);
                            const double w_gauss = r_integrations_points[i_gauss].Weight() * jac_vect[i_gauss];
                            const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                            const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                            const double total_energy_gauss = kinetic_energy + potential_energy;
                            air_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                            total_area_positive += w_gauss;

                            Total_energy_air += total_energy_gauss;
                            Total_energy_air_entire_element += total_energy_gauss;
                            Total_cinematic_energy_air += kinetic_energy;
                            Total_hydrostatic_energy_air += potential_energy;
                        }
                    }
                    else
                    {
                        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
                        {
                            const double Water_Density = 1e+3;
                            Vector N_gauss = row(N, i_gauss);
                            const double w_gauss = r_integrations_points[i_gauss].Weight() * jac_vect[i_gauss];
                            const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                            const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                            const double total_energy_gauss = kinetic_energy + potential_energy;

                            water_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                            total_area_negative += w_gauss;

                            Total_energy_water += total_energy_gauss;
                            Total_energy_water_entire_element += total_energy_gauss;
                            Total_cinematic_energy_water += kinetic_energy;
                            Total_hydrostatic_energy_water += potential_energy;
                        }
                    }
                }
            }
        }

        else if (mDomainSize == 3)
        {


            double total_energy = 0.0;
            double total_potential_energy = 0.0;
            double total_kinematic_energy = 0.0;

            for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem)
            {
                auto &r_geom = it_elem->GetGeometry();
                const auto elem_dist = SetDistanceVector(*it_elem);

                // Split element case
                if (this->IsCut(elem_dist))
                {
                    // Generate the splitting pattern
                    DivideTetrahedra3D4<Node>::Pointer p_divide_util = Kratos::make_shared<DivideTetrahedra3D4<Node>>(r_geom, elem_dist);
                    p_divide_util->GenerateDivision();

                    // Generate the modified shape functions util
                    auto p_geom = it_elem->pGetGeometry();
                    ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, elem_dist);

                    Matrix N_neg_side;
                    Vector w_gauss_neg_side;
                    Matrix N_neg_side_interface;
                    Vector w_gauss_interface;

                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side_interface;
                    p_ausas_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(N_neg_side, DN_DX_neg_side, w_gauss_neg_side, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    p_ausas_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(N_neg_side_interface, DN_DX_neg_side_interface, w_gauss_interface, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    // Save the coordinates of all the subdivision Gauss pts.
                    const unsigned int n_neg_gauss_pt = w_gauss_neg_side.size();
                    const unsigned int n_interface_side_gauss_pt = w_gauss_interface.size();
                    for (unsigned int i_gauss_interface = 0; i_gauss_interface < n_interface_side_gauss_pt; ++i_gauss_interface)
                    {
                        interface_area += w_gauss_interface[i_gauss_interface];
                        ;
                    }

                    for (unsigned int i_gauss = 0; i_gauss < n_neg_gauss_pt; ++i_gauss)
                    {
                        const double Water_Density = 1e+3;
                        Vector N_gauss = row(N_neg_side, i_gauss);
                        const double w_gauss = w_gauss_neg_side[i_gauss];
                        const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Water_Density);

                        // const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                        water_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                        // const double total_energy_gauss = kinetic_energy + potential_energy;

                        total_area_negative += w_gauss;

                        // Total_energy_water += total_energy_gauss;
                        // Total_energy_water_splitted+= total_energy_gauss;
                        Total_cinematic_energy_water += kinetic_energy;
                    }
                    Matrix N_pos_side;
                    Vector w_gauss_pos_side;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_pos_side;
                    p_ausas_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(N_pos_side, DN_DX_pos_side, w_gauss_pos_side, GeometryData::IntegrationMethod::GI_GAUSS_2);

                    // Save the coordinates of all the subdivision Gauss pts.
                    const unsigned int n_pos_gauss_pt = w_gauss_pos_side.size();

                    for (unsigned int i_gauss = 0; i_gauss < n_pos_gauss_pt; ++i_gauss)
                    {
                        const double Air_Density = 1.22;
                        Vector N_gauss = row(N_pos_side, i_gauss);
                        const double w_gauss = w_gauss_pos_side[i_gauss];
                        const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                        const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                        const double total_energy_gauss = kinetic_energy + potential_energy;
                        air_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);

                        total_area_positive += w_gauss;

                        Total_energy_air += total_energy_gauss;
                        Total_energy_air_splitted += total_energy_gauss;
                        Total_cinematic_energy_air += kinetic_energy;

                        // Total_hydrostatic_energy_air+=potential_energy;
                    }
                }
                // No split element case
                else
                {
                    // Get geometry data
                    Vector jac_vect;
                    r_geom.DeterminantOfJacobian(jac_vect, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    auto N = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);
                    const auto &r_integrations_points = r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
                    const unsigned int n_gauss = r_integrations_points.size();
                    // Check if the element is in the positive (fluid) region
                    if (this->IsPositive(elem_dist))
                    {
                        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
                        {
                            const double Air_Density = 1.22;
                            Vector N_gauss = row(N, i_gauss);
                            const double w_gauss = r_integrations_points[i_gauss].Weight() * jac_vect[i_gauss];
                            const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                            const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                            const double total_energy_gauss = kinetic_energy + potential_energy;
                            air_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                            total_area_positive += w_gauss;

                            Total_energy_air += total_energy_gauss;
                            Total_energy_air_entire_element += total_energy_gauss;
                            Total_cinematic_energy_air += kinetic_energy;
                            Total_hydrostatic_energy_air += potential_energy;
                        }
                    }
                    else
                    {
                        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
                        {
                            const double Water_Density = 1e+3;
                            Vector N_gauss = row(N, i_gauss);
                            const double w_gauss = r_integrations_points[i_gauss].Weight() * jac_vect[i_gauss];
                            const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                            const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                            const double total_energy_gauss = kinetic_energy + potential_energy;

                            water_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                            total_area_negative += w_gauss;

                            Total_energy_water += total_energy_gauss;
                            Total_energy_water_entire_element += total_energy_gauss;
                            Total_cinematic_energy_water += kinetic_energy;
                            Total_hydrostatic_energy_water += potential_energy;
                        }
                    }
                }
            }
        }

        air_zg /= total_area_positive;
        water_zg /= total_area_negative;

        KRATOS_WATCH(air_zg)
        KRATOS_WATCH(water_zg)

        // const double Potential_energy_water_test=PotentialEnergyTest(total_area_negative,Total_hydrostatic_energy_water_elemental,Water_Density);
        // const double Potential_energy_air_test=PotentialEnergyTest(total_area_positive,Total_hydrostatic_energy_air_elemental,Air_Density);
        const double Potential_energy_water_test = PotentialEnergyRuben(total_area_negative, water_zg, Water_Density, gravity);
        const double Potential_energy_air_test = PotentialEnergyRuben(total_area_positive, air_zg, Air_Density, gravity);

        const double Total_energy_air_system = Potential_energy_air_test + Total_cinematic_energy_air;
        const double Total_energy_water_system = Potential_energy_water_test + Total_cinematic_energy_water;
        const double total_energy_system = Total_energy_air_system + Total_energy_water_system;
        // Print the obtained values
        std::cout << std::endl;
        std::cout << "Total_potential_water_energy" << Potential_energy_water_test << std::endl;
        std::cout << "Total_cinematic_water_energy" << Total_cinematic_energy_water << std::endl;
        std::cout << "Total_potential_air_energy" << Potential_energy_air_test << std::endl;
        std::cout << "Total_cinematic_air_energy" << Total_cinematic_energy_air << std::endl;
        std::cout << "Total_energy_air:" << Total_energy_air_system << std::endl;
        std::cout << "Total_energy_water:" << Total_energy_water_system << std::endl;
        std::cout << "total_energy_system:" << total_energy_system << std::endl;

        std::cout << std::endl;

        WritingFile(Potential_energy_water_test, Total_cinematic_energy_water, Total_energy_water_system, Potential_energy_air_test, Total_cinematic_energy_air, Total_energy_air_system, total_energy_system, air_zg, water_zg, interface_area, total_area_negative, total_area_positive);


    }

    ///@}
    ///@name Friends
    ///@{

    Vector EnergyCalculator()
    {


        // Variables initialization
        double cut_area_negative = 0.0;
        double cut_area_positive = 0.0;
        double total_area_negative = 0.0;
        double total_area_positive = 0.0;
        double Total_energy_air_splitted = 0.0;
        double Total_energy_water_splitted = 0.0;
        double Total_energy_air = 0.0;
        double Total_energy_water = 0.0;
        double Total_energy_air_entire_element = 0;
        double Total_energy_water_entire_element = 0;
        double Total_cinematic_energy_water = 0;
        double Total_cinematic_energy_air = 0;
        double Total_hydrostatic_energy_water = 0;
        double Total_hydrostatic_energy_air = 0;
        double Total_cinematic_energy_water_splitted = 0;
        double Total_hydrostatic_energy_water_splitted = 0;
        double Total_cinematic_energy_water_not_splitted = 0;
        double Total_hydrostatic_energy_water_not_splitted = 0;
        double Total_cinematic_energy_air_splitted = 0;
        double Total_hydrostatic_energy_air_splitted = 0;
        double Total_cinematic_energy_air_not_splitted = 0;
        double Total_hydrostatic_energy_air_not_splitted = 0;
        double gravity = 9.81;
        double Total_hydrostatic_energy_water_elemental = 0;
        double Total_potential_energy_water_fluid = 0;
        double Total_hydrostatic_energy_air_elemental = 0;
        double Total_potential_energy_air_fluid = 0;
        const double Water_Density = 1e+3;
        const double Air_Density = 1.22;
        double interface_area = 0.0;
        double air_zg = 0.0;
        double water_zg = 0.0;

        if (mDomainSize == 2)
        // Total Energy check for 2D cases[Fluid1_energy and FLuid2_energy]: A loop in each element of the mesh considering splitted elements (air/water elements)
        {
            double total_energy = 0.0;
            double total_potential_energy = 0.0;
            double total_kinematic_energy = 0.0;

            for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem)
            {
                auto &r_geom = it_elem->GetGeometry();
                const auto elem_dist = SetDistanceVector(*it_elem);

                // Split element case
                if (this->IsCut(elem_dist))
                {
                    // Generate the splitting pattern

                    DivideTriangle2D3<Node>::Pointer p_divide_util = Kratos::make_shared<DivideTriangle2D3<Node>>(r_geom, elem_dist);
                    p_divide_util->GenerateDivision();

                    // Generate the modified shape functions util
                    auto p_geom = it_elem->pGetGeometry();
                    ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, elem_dist);

                    Matrix N_neg_side;
                    Vector w_gauss_neg_side;
                    Matrix N_neg_side_interface;
                    Vector w_gauss_interface;

                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side_interface;
                    p_ausas_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(N_neg_side, DN_DX_neg_side, w_gauss_neg_side, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    p_ausas_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(N_neg_side_interface, DN_DX_neg_side_interface, w_gauss_interface, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    // Save the coordinates of all the subdivision Gauss pts.nl
                    const unsigned int n_neg_gauss_pt = w_gauss_neg_side.size();
                    const unsigned int n_interface_side_gauss_pt = w_gauss_interface.size();
                    for (unsigned int i_gauss_interface = 0; i_gauss_interface < n_interface_side_gauss_pt; ++i_gauss_interface)
                    {
                        interface_area += w_gauss_interface[i_gauss_interface];
                        ;
                    }

                    for (unsigned int i_gauss = 0; i_gauss < n_neg_gauss_pt; ++i_gauss)
                    {
                        const double Water_Density = 1e+3;
                        Vector N_gauss = row(N_neg_side, i_gauss);
                        const double w_gauss = w_gauss_neg_side[i_gauss];
                        const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Water_Density);

                        // const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                        water_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                        // const double total_energy_gauss = kinetic_energy + potential_energy;

                        total_area_negative += w_gauss;

                        // Total_energy_water += total_energy_gauss;
                        // Total_energy_water_splitted+= total_energy_gauss;
                        Total_cinematic_energy_water += kinetic_energy;
                    }
                    Matrix N_pos_side;
                    Vector w_gauss_pos_side;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_pos_side;
                    p_ausas_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(N_pos_side, DN_DX_pos_side, w_gauss_pos_side, GeometryData::IntegrationMethod::GI_GAUSS_2);

                    // Save the coordinates of all the subdivision Gauss pts.
                    const unsigned int n_pos_gauss_pt = w_gauss_pos_side.size();

                    for (unsigned int i_gauss = 0; i_gauss < n_pos_gauss_pt; ++i_gauss)
                    {
                        const double Air_Density = 1.22;
                        Vector N_gauss = row(N_pos_side, i_gauss);
                        const double w_gauss = w_gauss_pos_side[i_gauss];
                        const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                        const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                        const double total_energy_gauss = kinetic_energy + potential_energy;
                        air_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);

                        total_area_positive += w_gauss;

                        Total_energy_air += total_energy_gauss;
                        Total_energy_air_splitted += total_energy_gauss;
                        Total_cinematic_energy_air += kinetic_energy;

                        // Total_hydrostatic_energy_air+=potential_energy;
                    }
                }
                // No split element case
                else
                {
                    // Get geometry data
                    Vector jac_vect;
                    r_geom.DeterminantOfJacobian(jac_vect, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    auto N = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);
                    const auto &r_integrations_points = r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
                    const unsigned int n_gauss = r_integrations_points.size();
                    // Check if the element is in the positive (fluid) region
                    if (this->IsPositive(elem_dist))
                    {
                        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
                        {
                            const double Air_Density = 1.22;
                            Vector N_gauss = row(N, i_gauss);
                            const double w_gauss = r_integrations_points[i_gauss].Weight() * jac_vect[i_gauss];
                            const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                            const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                            const double total_energy_gauss = kinetic_energy + potential_energy;
                            air_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                            total_area_positive += w_gauss;

                            Total_energy_air += total_energy_gauss;
                            Total_energy_air_entire_element += total_energy_gauss;
                            Total_cinematic_energy_air += kinetic_energy;
                            Total_hydrostatic_energy_air += potential_energy;
                        }
                    }
                    else
                    {
                        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
                        {
                            const double Water_Density = 1e+3;
                            Vector N_gauss = row(N, i_gauss);
                            const double w_gauss = r_integrations_points[i_gauss].Weight() * jac_vect[i_gauss];
                            const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                            const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                            const double total_energy_gauss = kinetic_energy + potential_energy;

                            water_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                            total_area_negative += w_gauss;

                            Total_energy_water += total_energy_gauss;
                            Total_energy_water_entire_element += total_energy_gauss;
                            Total_cinematic_energy_water += kinetic_energy;
                            Total_hydrostatic_energy_water += potential_energy;
                        }
                    }
                }
            }
        }

        else if (mDomainSize == 3)
        {

            double total_energy = 0.0;
            double total_potential_energy = 0.0;
            double total_kinematic_energy = 0.0;

            for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem)
            {
                auto &r_geom = it_elem->GetGeometry();
                const auto elem_dist = SetDistanceVector(*it_elem);

                // Split element case
                if (this->IsCut(elem_dist))
                {
                    // Generate the splitting pattern
                    DivideTetrahedra3D4<Node>::Pointer p_divide_util = Kratos::make_shared<DivideTetrahedra3D4<Node>>(r_geom, elem_dist);
                    p_divide_util->GenerateDivision();

                    // Generate the modified shape functions util
                    auto p_geom = it_elem->pGetGeometry();
                    ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, elem_dist);

                    Matrix N_neg_side;
                    Vector w_gauss_neg_side;
                    Matrix N_neg_side_interface;
                    Vector w_gauss_interface;

                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side_interface;
                    p_ausas_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(N_neg_side, DN_DX_neg_side, w_gauss_neg_side, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    p_ausas_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(N_neg_side_interface, DN_DX_neg_side_interface, w_gauss_interface, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    // Save the coordinates of all the subdivision Gauss pts.
                    const unsigned int n_neg_gauss_pt = w_gauss_neg_side.size();
                    const unsigned int n_interface_side_gauss_pt = w_gauss_interface.size();
                    for (unsigned int i_gauss_interface = 0; i_gauss_interface < n_interface_side_gauss_pt; ++i_gauss_interface)
                    {
                        interface_area += w_gauss_interface[i_gauss_interface];
                        ;
                    }

                    for (unsigned int i_gauss = 0; i_gauss < n_neg_gauss_pt; ++i_gauss)
                    {
                        const double Water_Density = 1e+3;
                        Vector N_gauss = row(N_neg_side, i_gauss);
                        const double w_gauss = w_gauss_neg_side[i_gauss];
                        const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Water_Density);

                        // const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                        water_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                        // const double total_energy_gauss = kinetic_energy + potential_energy;

                        total_area_negative += w_gauss;

                        // Total_energy_water += total_energy_gauss;
                        // Total_energy_water_splitted+= total_energy_gauss;
                        Total_cinematic_energy_water += kinetic_energy;
                    }
                    Matrix N_pos_side;
                    Vector w_gauss_pos_side;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_pos_side;
                    p_ausas_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(N_pos_side, DN_DX_pos_side, w_gauss_pos_side, GeometryData::IntegrationMethod::GI_GAUSS_2);

                    // Save the coordinates of all the subdivision Gauss pts.
                    const unsigned int n_pos_gauss_pt = w_gauss_pos_side.size();

                    for (unsigned int i_gauss = 0; i_gauss < n_pos_gauss_pt; ++i_gauss)
                    {
                        const double Air_Density = 1.22;
                        Vector N_gauss = row(N_pos_side, i_gauss);
                        const double w_gauss = w_gauss_pos_side[i_gauss];
                        const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                        const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                        const double total_energy_gauss = kinetic_energy + potential_energy;
                        air_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);

                        total_area_positive += w_gauss;

                        Total_energy_air += total_energy_gauss;
                        Total_energy_air_splitted += total_energy_gauss;
                        Total_cinematic_energy_air += kinetic_energy;

                        // Total_hydrostatic_energy_air+=potential_energy;
                    }
                }
                // No split element case
                else
                {
                    // Get geometry data
                    Vector jac_vect;
                    r_geom.DeterminantOfJacobian(jac_vect, GeometryData::IntegrationMethod::GI_GAUSS_2);
                    auto N = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);
                    const auto &r_integrations_points = r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
                    const unsigned int n_gauss = r_integrations_points.size();
                    // Check if the element is in the positive (fluid) region
                    if (this->IsPositive(elem_dist))
                    {
                        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
                        {
                            const double Air_Density = 1.22;
                            Vector N_gauss = row(N, i_gauss);
                            const double w_gauss = r_integrations_points[i_gauss].Weight() * jac_vect[i_gauss];
                            const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                            const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Air_Density);
                            const double total_energy_gauss = kinetic_energy + potential_energy;
                            air_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                            total_area_positive += w_gauss;

                            Total_energy_air += total_energy_gauss;
                            Total_energy_air_entire_element += total_energy_gauss;
                            Total_cinematic_energy_air += kinetic_energy;
                            Total_hydrostatic_energy_air += potential_energy;
                        }
                    }
                    else
                    {
                        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
                        {
                            const double Water_Density = 1e+3;
                            Vector N_gauss = row(N, i_gauss);
                            const double w_gauss = r_integrations_points[i_gauss].Weight() * jac_vect[i_gauss];
                            const double kinetic_energy = CalculateGaussPointKinematicEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                            const double potential_energy = CalculateGaussPointPotentialEnergy(r_geom, N_gauss, w_gauss, Water_Density);
                            const double total_energy_gauss = kinetic_energy + potential_energy;

                            water_zg += CalculateGaussPointPotentialEnergyGravityCenter(r_geom, N_gauss, w_gauss);
                            total_area_negative += w_gauss;

                            Total_energy_water += total_energy_gauss;
                            Total_energy_water_entire_element += total_energy_gauss;
                            Total_cinematic_energy_water += kinetic_energy;
                            Total_hydrostatic_energy_water += potential_energy;
                        }
                    }
                }
            }
        }

        air_zg /= total_area_positive;
        water_zg /= total_area_negative;

        KRATOS_WATCH(air_zg)
        KRATOS_WATCH(water_zg)

        // const double Potential_energy_water_test=PotentialEnergyTest(total_area_negative,Total_hydrostatic_energy_water_elemental,Water_Density);
        // const double Potential_energy_air_test=PotentialEnergyTest(total_area_positive,Total_hydrostatic_energy_air_elemental,Air_Density);
        const double Potential_energy_water_test = PotentialEnergyRuben(total_area_negative, water_zg, Water_Density, gravity);
        const double Potential_energy_air_test = PotentialEnergyRuben(total_area_positive, air_zg, Air_Density, gravity);

        const double Total_energy_air_system = Potential_energy_air_test + Total_cinematic_energy_air;
        const double Total_energy_water_system = Potential_energy_water_test + Total_cinematic_energy_water;
        const double total_energy_system = Total_energy_air_system + Total_energy_water_system;
        // Print the obtained values
        std::cout << std::endl;
        std::cout << "Total_potential_water_energy" << Potential_energy_water_test << std::endl;
        std::cout << "Total_cinematic_water_energy" << Total_cinematic_energy_water << std::endl;
        std::cout << "Total_potential_air_energy" << Potential_energy_air_test << std::endl;
        std::cout << "Total_cinematic_air_energy" << Total_cinematic_energy_air << std::endl;
        std::cout << "Total_energy_air:" << Total_energy_air_system << std::endl;
        std::cout << "Total_energy_water:" << Total_energy_water_system << std::endl;
        std::cout << "total_energy_system:" << total_energy_system << std::endl;

        std::cout << std::endl;

        // WritingFile(Potential_energy_water_test, Total_cinematic_energy_water, Total_energy_water_system, Potential_energy_air_test, Total_cinematic_energy_air, Total_energy_air_system, total_energy_system, air_zg, water_zg, interface_area, total_area_negative, total_area_positive);
        Vector WaterEnergy = ZeroVector(5);

        WaterEnergy[0] = Total_energy_water_system;
        WaterEnergy[1] = Potential_energy_water_test;
        WaterEnergy[2] = Total_cinematic_energy_water;
        WaterEnergy[3] = water_zg;
        WaterEnergy[4] = interface_area;
        return WaterEnergy;
    }
    double CalculateInterfaceArea()
    {
        // Inicializar el área de la interfaz
        double interface_area = 0.0;

        // Recorrer todos los elementos del ModelPart
        for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem)
        {
            auto &r_geom = it_elem->GetGeometry();
            const auto elem_dist = SetDistanceVector(*it_elem);

            // Comprobar si el elemento está cortado por la interfaz
            if (this->IsCut(elem_dist))
            {
                // Utilidades específicas para 2D o 3D
                if (mDomainSize == 2)
                {
                    // Crear la utilidad para dividir triángulos
                    DivideTriangle2D3<Node>::Pointer p_divide_util = Kratos::make_shared<DivideTriangle2D3<Node>>(r_geom, elem_dist);
                    p_divide_util->GenerateDivision();

                    // Crear las funciones de forma modificadas
                    auto p_geom = it_elem->pGetGeometry();
                    ModifiedShapeFunctions::Pointer p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, elem_dist);

                    // Calcular el área de la interfaz
                    Matrix N_neg_side_interface;
                    Vector w_gauss_interface;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side_interface;

                    p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                        N_neg_side_interface, DN_DX_neg_side_interface, w_gauss_interface, GeometryData::IntegrationMethod::GI_GAUSS_2);

                    // Sumar el área de la interfaz en este elemento
                    const unsigned int n_interface_side_gauss_pt = w_gauss_interface.size();
                    for (unsigned int i_gauss_interface = 0; i_gauss_interface < n_interface_side_gauss_pt; ++i_gauss_interface)
                    {
                        interface_area += w_gauss_interface[i_gauss_interface];
                    }
                }
                else if (mDomainSize == 3)
                {
                    // Crear la utilidad para dividir tetraedros
                    DivideTetrahedra3D4<Node>::Pointer p_divide_util = Kratos::make_shared<DivideTetrahedra3D4<Node>>(r_geom, elem_dist);
                    p_divide_util->GenerateDivision();

                    // Crear las funciones de forma modificadas
                    auto p_geom = it_elem->pGetGeometry();
                    ModifiedShapeFunctions::Pointer p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, elem_dist);

                    // Calcular el área de la interfaz
                    Matrix N_neg_side_interface;
                    Vector w_gauss_interface;
                    Geometry<Node>::ShapeFunctionsGradientsType DN_DX_neg_side_interface;

                    p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                        N_neg_side_interface, DN_DX_neg_side_interface, w_gauss_interface, GeometryData::IntegrationMethod::GI_GAUSS_2);

                    // Sumar el área de la interfaz en este elemento
                    const unsigned int n_interface_side_gauss_pt = w_gauss_interface.size();
                    for (unsigned int i_gauss_interface = 0; i_gauss_interface < n_interface_side_gauss_pt; ++i_gauss_interface)
                    {
                        interface_area += w_gauss_interface[i_gauss_interface];
                    }
                }
            }
        }

        // Devolver el área de la interfaz total
        return interface_area;
    }

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    bool IsCut(const Vector& rElementalDistances)
    {
        unsigned int n_pos = 0;
        unsigned int n_neg = 0;

        for (unsigned int i_node = 0; i_node < rElementalDistances.size(); ++i_node){
            if(rElementalDistances[i_node] > 0.0) {
                n_pos++;
            } else {
                n_neg++;
            }
        }

        if (n_pos != 0 && n_neg != 0){
            return true;
        } else {
            return false;
        }
    }

    bool IsPositive(const Vector& rElementalDistances)
    {
        unsigned int n_pos = 0;
        const unsigned int n_nodes = rElementalDistances.size();
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node){
            if(rElementalDistances[i_node] > 0.0) {
                n_pos++;
            }
        }

        if (n_pos == n_nodes){
            return true;
        } else {
            return false;
        }
    }

    Vector SetDistanceVector(const Element& rElement)
    {
        const auto& r_geom = rElement.GetGeometry();
        const unsigned int n_nodes = r_geom.PointsNumber();
        Vector distances_vector(n_nodes);
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            distances_vector[i_node] = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
        }
        return distances_vector;
    }

    double CalculateGaussPointKinematicEnergy(
        const Geometry<Node>& rGeometry,
        const Vector& rN,
        const double Weight,
        const double Density)
    {
        // Gauss pt. interpolation of velocity
        array_1d<double,3> v_gauss = ZeroVector(3);
        const unsigned int n_nodes = rGeometry.PointsNumber();
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            v_gauss += rN[i_node] * rGeometry[i_node].FastGetSolutionStepValue(VELOCITY);
        }

        // Calculate and return Gauss pt. kinetic energy
        const double squared_v_norm = inner_prod(v_gauss, v_gauss);
        return Weight * Density * squared_v_norm / 2.0;
    }

    double CalculateGaussPointPotentialEnergy(
        const Geometry<Node>& rGeometry,
        const Vector& rN,
        const double Weight,
        const double Density)
    {   // Gauss pt. interpolation of elevation
        const double gravity = 9.81;
        double z_interpolated = 0.0;
        const unsigned int n_nodes = rGeometry.PointsNumber();
// It should be improve in order to have z_interpolated depending on the directions of the problem gravity
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node){
            z_interpolated += rGeometry[i_node].Z() * rN[i_node];
        }
        // Calculate and return Gauss pt. potential energy
        return (Density * gravity * z_interpolated) * Weight;

    }

        double CalculateGaussPointPotentialEnergyGravityCenter(
        const Geometry<Node>& rGeometry,
        const Vector& rN,
        const double Weight)
    {
        double z_interpolated = 0.0;
        const unsigned int n_nodes = rGeometry.PointsNumber();
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node){
            z_interpolated += rGeometry[i_node].Z() * rN[i_node];
        }

        return (z_interpolated) * Weight;

    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    unsigned int mDomainSize;

    std::string mFileName;
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

    /// Assignment operator.
    EnergyCheckProcess&operator=(EnergyCheckProcess const &rOther);

    /// Copy constructor.
    //GaussPointErrorProcessAnalytical(GaussPointErrorProcessAnalytical const& rOther);


    ///@}

}; // Class GaussPointErrorProcessAnalytical

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
// inline std::istream &operator>>(std::istream &rIStream,
//                                 EnergyCheckProcess&rThis);

// /// output stream function
// inline std::ostream &operator<<(std::ostream &rOStream,
//                                 const EnergyCheckProcess&rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}


}  // namespace Kratos.

#endif // KRATOS_GAUSS_POINT_ERROR_PROCESS_H_INCLUDED  defined