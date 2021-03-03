
Skip to content
Pull requests
Issues
Marketplace
Explore
@uxuech
Learn Git and GitHub without any code!

Using the Hello World guide, you’ll start a branch, write comments, and open a pull request.
KratosMultiphysics /
Kratos

51
551

    132

Code
Issues 206
Pull requests 170
Discussions
Actions
Projects 28
Wiki
Security

    Insights

Kratos/applications/FluidDynamicsApplication/custom_processes/energy_splitelements_process.h
@uxuech
uxuech energy_check_process
Latest commit bacdd0e 2 days ago
History
1 contributor
764 lines (614 sloc) 31.3 KB
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

#if !defined(KRATOS_ENERGY_SPLIT_LEMENTS_PROCESS_H_INCLUDED)
#define KRATOS_ENERGY_SPLIT_ELEMENTS_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

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
        unsigned int DomainSize)
        : mrModelPart(rModelPart),
          mDomainSize(DomainSize)
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

    void Execute() override
    {
        KRATOS_TRY;

        // Variables initialization
        double cut_area_negative= 0.0;
        double cut_area_positive= 0.0;
        double total_area_negative= 0.0;
        double total_area_positive= 0.0;
        double Total_energy_air_splitted=0.0;
        double Total_energy_water_splitted = 0.0;
        double Total_energy_air = 0.0;
        double Total_energy_water= 0.0;

        if(mDomainSize == 2)
        // Total Energy check for 2D cases[Fluid1_energy and FLuid2_energy]: A loop in each element of the mesh considering splitted elements (air/water elements)
        { 
            for(auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem) 
            {
                auto &r_geom = it_elem->GetGeometry();
                const auto elem_dist = SetDistanceVector(*it_elem);

                const auto &r_p_0 = r_geom[0].FastGetSolutionStepValue(PRESSURE);
                const auto &r_p_1 = r_geom[1].FastGetSolutionStepValue(PRESSURE);
                const auto &r_p_2 = r_geom[2].FastGetSolutionStepValue(PRESSURE);
                const auto &r_vel_0 = r_geom[0].FastGetSolutionStepValue(VELOCITY);
                const auto &r_vel_1 = r_geom[1].FastGetSolutionStepValue(VELOCITY);
                const auto &r_vel_2 = r_geom[2].FastGetSolutionStepValue(VELOCITY);

                // const bool is_valid = aux_n == r_geom.PointsNumber() ? true : false;
                const bool is_valid = true;
                if (is_valid) {
                    if (this->IsCut(elem_dist)){
                        // Generate the splitting pattern
                        DivideTriangle2D3::Pointer p_divide_util = Kratos::make_shared<DivideTriangle2D3>(r_geom, elem_dist);
                        p_divide_util->GenerateDivision();

                        // Generate the modified shape functions util
                        auto p_geom = it_elem->pGetGeometry();
                        ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, elem_dist);

                        Matrix N_neg_side;
                        Vector w_gauss_neg_side;
                        Geometry<Node<3>>::ShapeFunctionsGradientsType DN_DX_neg_side;
                        p_ausas_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(N_neg_side, DN_DX_neg_side, w_gauss_neg_side, GeometryData::GI_GAUSS_2);

                        // Save the coordinates of all the subdivision Gauss pts.
                        const unsigned int n_neg_gauss_pt = w_gauss_neg_side.size();
                        const unsigned int n_neg_sbdiv = (p_divide_util->mNegativeSubdivisions).size();

                        for (unsigned int i_sbdiv = 0; i_sbdiv < n_neg_sbdiv; ++i_sbdiv){
                            auto p_sbdiv = (p_divide_util->mNegativeSubdivisions)[i_sbdiv];
                            auto N = p_sbdiv->ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                            for (unsigned int i_gauss = 0;  i_gauss < 3; ++i_gauss){
                                array_1d<double, 3> v_solu;
                                array_1d<double, 3> energy_solu;
                                const double p_solu = r_p_0*N_neg_side(i_gauss,0) + r_p_1*N_neg_side(i_gauss,1) + r_p_2*N_neg_side(i_gauss,2);
                                v_solu[0] = r_vel_0[0]*N_neg_side(i_gauss,0) + r_vel_1[0]*N_neg_side(i_gauss,1) + r_vel_2[0]*N_neg_side(i_gauss,2);
                                v_solu[1] = r_vel_0[1]*N_neg_side(i_gauss,0) + r_vel_1[1]*N_neg_side(i_gauss,1) + r_vel_2[1]*N_neg_side(i_gauss,2);
                                v_solu[2] = r_vel_0[2]*N_neg_side(i_gauss,0) + r_vel_1[2]*N_neg_side(i_gauss,1) + r_vel_2[2]*N_neg_side(i_gauss,2);

                                const double Water_Density=1e+3;
                                energy_solu[0] =( (std::pow(v_solu[0],2) * Water_Density) / 2) + p_solu;
                                energy_solu[1] = ((std::pow(v_solu[1], 2) * Water_Density) / 2) + p_solu;
                                energy_solu[2] = ((std::pow(v_solu[2] , 2) * Water_Density) / 2) + p_solu;
                                const double Total_energy_water = std::pow(std::pow(energy_solu[0], 2) + std::pow(energy_solu[1], 2) + std::pow(energy_solu[2] , 2 ) , 0.5);

                                Total_energy_water_splitted += w_gauss_neg_side[i_gauss] * Total_energy_water;
                                cut_area_negative += w_gauss_neg_side[i_gauss];
                                total_area_negative += w_gauss_neg_side[i_gauss];
                            }
                        }

                        
                        Matrix N_pos_side;
                        Vector w_gauss_pos_side;
                        Geometry<Node<3>>::ShapeFunctionsGradientsType DN_DX_pos_side;
                        p_ausas_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(N_pos_side, DN_DX_pos_side, w_gauss_pos_side, GeometryData::GI_GAUSS_2);

                        // Save the coordinates of all the subdivision Gauss pts.
                        const unsigned int n_pos_gauss_pt = w_gauss_pos_side.size();
                        const unsigned int n_pos_sbdiv = (p_divide_util->mPositiveSubdivisions).size();

                        for (unsigned int i_sbdiv = 0; i_sbdiv < n_pos_sbdiv; ++i_sbdiv)
                        {
                                auto p_sbdiv = (p_divide_util->mPositiveSubdivisions)[i_sbdiv];
                                auto N = p_sbdiv->ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                                for (unsigned int i_gauss = 0; i_gauss < 3; ++i_gauss)
                                {
                                    array_1d<double, 3> v_solu;
                                    array_1d<double, 3> energy_solu;
                                    const double p_solu = r_p_0 * N_pos_side(i_gauss, 0) + r_p_1 * N_pos_side(i_gauss, 1) + r_p_2 * N_pos_side(i_gauss, 2);
                                    v_solu[0] = r_vel_0[0] * N_pos_side(i_gauss, 0) + r_vel_1[0] * N_pos_side(i_gauss, 1) + r_vel_2[0] * N_pos_side(i_gauss, 2);
                                    v_solu[1] = r_vel_0[1] * N_pos_side(i_gauss, 0) + r_vel_1[1] * N_pos_side(i_gauss, 1) + r_vel_2[1] * N_pos_side(i_gauss, 2);
                                    v_solu[2] = r_vel_0[2] * N_pos_side(i_gauss, 0) + r_vel_1[2] * N_pos_side(i_gauss, 1) + r_vel_2[2] * N_pos_side(i_gauss, 2);
                                    const double Air_Density = 1.22;
                                    energy_solu[0] = ((std::pow(v_solu[0], 2) * Air_Density) / 2) + p_solu;
                                    energy_solu[1] = ((std::pow(v_solu[1], 2) * Air_Density) / 2) + p_solu;
                                    energy_solu[2] = ((std::pow(v_solu[2], 2) * Air_Density) / 2) + p_solu;
                                    const double Air_Energy_total = std::pow(std::pow(energy_solu[0], 2) + std::pow(energy_solu[1], 2) + std::pow(energy_solu[2], 2), 0.5);

                                    Total_energy_air_splitted += w_gauss_pos_side[i_gauss] * Air_Energy_total;
                                    Total_energy_air += w_gauss_pos_side[i_gauss] * Air_Energy_total;

                                    // Add the local Gauss contribution to the areas
                                    cut_area_positive += w_gauss_pos_side[i_gauss];
                                    total_area_positive += w_gauss_pos_side[i_gauss];
                                }




                            }

                            
                        }

                    // No split element case
                    else {
                        // Check if the element is in the positive (fluid) region
                        if (this->IsPositive(elem_dist)){
                            // Get geometry data
                            Vector jac_vect;
                            jac_vect = r_geom.DeterminantOfJacobian(jac_vect, GeometryData::GI_GAUSS_2);
                            auto N = r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                            for (unsigned int i_gauss = 0; i_gauss < 3; ++i_gauss){
                                // Obtained solution
                                array_1d<double, 3> v_solu;
                                array_1d<double, 3> energy_solu;
                                const double p_solu = r_p_0*N(i_gauss,0) + r_p_1*N(i_gauss,1) + r_p_2*N(i_gauss,2);
                                v_solu[0] = r_vel_0[0]*N(i_gauss,0) + r_vel_1[0]*N(i_gauss,1) + r_vel_2[0]*N(i_gauss,2);
                                v_solu[1] = r_vel_0[1]*N(i_gauss,0) + r_vel_1[1]*N(i_gauss,1) + r_vel_2[1]*N(i_gauss,2);
                                v_solu[2] = r_vel_0[2]*N(i_gauss,0) + r_vel_1[2]*N(i_gauss,1) + r_vel_2[2]*N(i_gauss,2);

                                const double Air_Density=1.22;
                                energy_solu[0] = ((std::pow(v_solu[0], 2) * Air_Density) / 2) + p_solu;
                                energy_solu[1] = ((std::pow(v_solu[1], 2) * Air_Density) / 2) + p_solu;
                                energy_solu[2] = ((std::pow(v_solu[2], 2) * Air_Density) / 2) + p_solu;
                                const double Air_energy = std::pow(std::pow(energy_solu[0], 2) + std::pow(energy_solu[1], 2) + std::pow(energy_solu[2], 2), 0.5);
            
                                Total_energy_air += Air_energy * (jac_vect[i_gauss] / 6.0);
                                // Add the local Gauss contribution to the areas
                                total_area_positive += (jac_vect[i_gauss] / 6.0);
                            }
                        }
                        
                        else if (this->IsNegative(elem_dist)){
                            // Get geometry data
                            Vector jac_vect;
                            jac_vect = r_geom.DeterminantOfJacobian(jac_vect, GeometryData::GI_GAUSS_2);
                            auto N = r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

                            // Obtained solution
                            for (unsigned int i_gauss = 0; i_gauss < 3; ++i_gauss)
                            {
                                array_1d<double, 3> v_solu;
                                array_1d<double, 3> energy_solu;

                                const double p_solu = r_p_0 * N(i_gauss, 0) + r_p_1 * N(i_gauss, 1) + r_p_2 * N(i_gauss, 2);
                                v_solu[0] = r_vel_0[0] * N(i_gauss, 0) + r_vel_1[0] * N(i_gauss, 1) + r_vel_2[0] * N(i_gauss, 2);
                                v_solu[1] = r_vel_0[1] * N(i_gauss, 0) + r_vel_1[1] * N(i_gauss, 1) + r_vel_2[1] * N(i_gauss, 2);
                                v_solu[2] = r_vel_0[2] * N(i_gauss, 0) + r_vel_1[2] * N(i_gauss, 1) + r_vel_2[2] * N(i_gauss, 2);

                                const double Water_Density = 1e+3;
                                energy_solu[0] = ((std::pow(v_solu[0], 2) * Water_Density) / 2) + p_solu;
                                energy_solu[1] = ((std::pow(v_solu[1], 2) * Water_Density) / 2) + p_solu;
                                energy_solu[2] = ((std::pow(v_solu[2], 2) * Water_Density) / 2) + p_solu;
                                const double energy_water = std::pow(std::pow(energy_solu[0], 2) + std::pow(energy_solu[1], 2) + std::pow(energy_solu[2], 2), 0.5);
                                Total_energy_water += energy_water * (jac_vect[i_gauss] / 6.0);
                                // Add the local Gauss contribution to the areas
                                total_area_negative += (jac_vect[i_gauss] / 6.0);
                            }
                        }
                    }
                }
        }
        }
        else if(mDomainSize == 3)
        {
            for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem)
            {
                auto &r_geom = it_elem->GetGeometry();
                const auto elem_dist = SetDistanceVector(*it_elem);

                const auto &r_p_0 = r_geom[0].FastGetSolutionStepValue(PRESSURE);
                const auto &r_p_1 = r_geom[1].FastGetSolutionStepValue(PRESSURE);
                const auto &r_p_2 = r_geom[2].FastGetSolutionStepValue(PRESSURE);
                const auto &r_p_3 = r_geom[3].FastGetSolutionStepValue(PRESSURE);

                const auto &r_vel_0 = r_geom[0].FastGetSolutionStepValue(VELOCITY);
                const auto &r_vel_1 = r_geom[1].FastGetSolutionStepValue(VELOCITY);
                const auto &r_vel_2 = r_geom[2].FastGetSolutionStepValue(VELOCITY);
                const auto &r_vel_3 = r_geom[3].FastGetSolutionStepValue(VELOCITY);

                // const bool is_valid = aux_n == r_geom.PointsNumber() ? true : false;
                const bool is_valid = true;

                // Split element case
                if (is_valid)
                {
                    if (this->IsCut3D(elem_dist))
                    {
                        // Generate the splitting pattern
                        DivideTetrahedra3D4::Pointer p_divide_util = Kratos::make_shared<DivideTetrahedra3D4>(r_geom, elem_dist);
                        p_divide_util->GenerateDivision();

                        // Generate the modified shape functions util
                        auto p_geom = it_elem->pGetGeometry();
                        ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, elem_dist);

                        Matrix N_neg_side;
                        Vector w_gauss_neg_side;
                        Geometry<Node<3>>::ShapeFunctionsGradientsType DN_DX_neg_side;
                        p_ausas_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(N_neg_side, DN_DX_neg_side, w_gauss_neg_side, GeometryData::GI_GAUSS_2);

                        // Save the coordinates of all the subdivision Gauss pts.
                        const unsigned int n_neg_gauss_pt = w_gauss_neg_side.size();
                        const unsigned int n_neg_sbdiv = (p_divide_util->mNegativeSubdivisions).size();

                        for (unsigned int i_sbdiv = 0; i_sbdiv < n_neg_sbdiv; ++i_sbdiv)
                        {
                            auto p_sbdiv = (p_divide_util->mNegativeSubdivisions)[i_sbdiv];
                            auto N = p_sbdiv->ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                            for (unsigned int i_gauss = 0; i_gauss < 4; ++i_gauss)
                            {
                                array_1d<double, 3> v_solu;
                                array_1d<double, 3> energy_solu;
                                const double p_solu = r_p_0 * N_neg_side(i_gauss, 0) + r_p_1 * N_neg_side(i_gauss, 1) + r_p_2 * N_neg_side(i_gauss, 2) + r_p_3 * N_neg_side(i_gauss, 3);
                                v_solu[0] = r_vel_0[0] * N_neg_side(i_gauss, 0) + r_vel_1[0] * N_neg_side(i_gauss, 1) + r_vel_2[0] * N_neg_side(i_gauss, 2) + r_vel_3[0] * N_neg_side(i_gauss, 3);
                                v_solu[1] = r_vel_0[1] * N_neg_side(i_gauss, 0) + r_vel_1[1] * N_neg_side(i_gauss, 1) + r_vel_2[1] * N_neg_side(i_gauss, 2) + r_vel_3[1] * N_neg_side(i_gauss, 3);
                                v_solu[2] = r_vel_0[2] * N_neg_side(i_gauss, 0) + r_vel_1[2] * N_neg_side(i_gauss, 1) + r_vel_2[2] * N_neg_side(i_gauss, 2) + r_vel_3[2] * N_neg_side(i_gauss, 3);

                                const double Water_Density = 1e+3;

                                energy_solu[0] = ((std::pow(v_solu[0], 2) * Water_Density) / 2) + p_solu;
                                energy_solu[1] = ((std::pow(v_solu[1], 2) * Water_Density) / 2) + p_solu;
                                energy_solu[2] = ((std::pow(v_solu[2], 2) * Water_Density) / 2) + p_solu;
                                const double Total_energy_water = std::pow(std::pow(energy_solu[0], 2) + std::pow(energy_solu[1], 2) + std::pow(energy_solu[2], 2), 0.5);

                                Total_energy_water_splitted += w_gauss_neg_side[i_gauss] * Total_energy_water;
                                cut_area_negative += w_gauss_neg_side[i_gauss];
                                total_area_negative += w_gauss_neg_side[i_gauss];
                            }

                        }
                        Matrix N_pos_side;
                        Vector w_gauss_pos_side;
                        Geometry<Node<3>>::ShapeFunctionsGradientsType DN_DX_pos_side;
                        p_ausas_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(N_pos_side, DN_DX_pos_side, w_gauss_pos_side, GeometryData::GI_GAUSS_2);

                        // Save the coordinates of all the subdivision Gauss pts.
                        const unsigned int n_pos_gauss_pt = w_gauss_pos_side.size();
                        const unsigned int n_pos_sbdiv = (p_divide_util->mPositiveSubdivisions).size();

                        for (unsigned int i_sbdiv = 0; i_sbdiv < n_pos_sbdiv; ++i_sbdiv)
                        {
                            auto p_sbdiv = (p_divide_util->mPositiveSubdivisions)[i_sbdiv];
                            auto N = p_sbdiv->ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                            for (unsigned int i_gauss = 0; i_gauss < 4; ++i_gauss)
                            {
                                array_1d<double, 3> v_solu;
                                array_1d<double, 3> energy_solu;
                                const double p_solu = r_p_0 * N_pos_side(i_gauss, 0) + r_p_1 * N_pos_side(i_gauss, 1) + r_p_2 * N_pos_side(i_gauss, 2)+ r_p_3 * N_pos_side(i_gauss, 3);
                                v_solu[0] = r_vel_0[0] * N_pos_side(i_gauss, 0) + r_vel_1[0] * N_pos_side(i_gauss, 1) + r_vel_2[0] * N_pos_side(i_gauss, 2) + r_vel_3[0] * N_pos_side(i_gauss, 3);
                                v_solu[1] = r_vel_0[1] * N_pos_side(i_gauss, 0) + r_vel_1[1] * N_pos_side(i_gauss, 1) + r_vel_2[1] * N_pos_side(i_gauss, 2)+ r_vel_3[1] * N_pos_side(i_gauss, 3);
                                v_solu[2] = r_vel_0[2] * N_pos_side(i_gauss, 0) + r_vel_1[2] * N_pos_side(i_gauss, 1) + r_vel_2[2] * N_pos_side(i_gauss, 2)+ r_vel_3[2] * N_pos_side(i_gauss, 3);
                                const double Air_Density = 1.22;
                                energy_solu[0] = ((std::pow(v_solu[0], 2) * Air_Density) / 2) + p_solu;
                                energy_solu[1] = ((std::pow(v_solu[1], 2) * Air_Density) / 2) + p_solu;
                                energy_solu[2] = ((std::pow(v_solu[2], 2) * Air_Density) / 2) + p_solu;
                                const double energy_air = std::pow(std::pow(energy_solu[0], 2) + std::pow(energy_solu[1], 2) + std::pow(energy_solu[2], 2), 0.5);
                                

                                Total_energy_air_splitted += w_gauss_pos_side[i_gauss] * energy_air;
                                Total_energy_air += w_gauss_pos_side[i_gauss] * energy_air;

                                // Add the local Gauss contribution to the areas
                                cut_area_positive += w_gauss_pos_side[i_gauss];
                                total_area_positive += w_gauss_pos_side[i_gauss];
                    }
                    }
                    }

                    // No split element case
                else
                {

                    // Check if the element is in the positive (fluid) region
                    if (this->IsPositive3D(elem_dist))
                    {
                        
                        // Get geometry data
                        Vector jac_vect;
                        jac_vect = r_geom.DeterminantOfJacobian(jac_vect, GeometryData::GI_GAUSS_2);
                        auto N = r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                        for (unsigned int i_gauss = 0; i_gauss < 4; ++i_gauss){
                            // Obtained solution
                            array_1d<double, 3> v_solu;
                            array_1d<double, 3> energy_solu;
                            const double p_solu = r_p_0 * N(i_gauss, 0) + r_p_1 * N(i_gauss, 1) + r_p_2 * N(i_gauss, 2) + r_p_3 * N(i_gauss, 3);
                            v_solu[0] = r_vel_0[0] * N(i_gauss, 0) + r_vel_1[0] * N(i_gauss, 1) + r_vel_2[0] * N(i_gauss, 2)+ r_vel_3[0] * N(i_gauss, 3);
                            v_solu[1] = r_vel_0[1] * N(i_gauss, 0) + r_vel_1[1] * N(i_gauss, 1) + r_vel_2[1] * N(i_gauss, 2)+r_vel_3[1] * N(i_gauss, 3);
                            v_solu[2] = r_vel_0[2] * N(i_gauss, 0) + r_vel_1[2] * N(i_gauss, 1) + r_vel_2[2] * N(i_gauss, 2) + r_vel_3[2] * N(i_gauss, 3);

                            const double Air_Density = 1.22;
                            energy_solu[0] = ((std::pow(v_solu[0], 2) * Air_Density) / 2) + p_solu;
                            energy_solu[1] = ((std::pow(v_solu[1], 2) * Air_Density) / 2) + p_solu;
                            energy_solu[2] = ((std::pow(v_solu[2], 2) * Air_Density) / 2) + p_solu;
                            const double energy_air = std::pow(std::pow(energy_solu[0], 2) + std::pow(energy_solu[1], 2) + std::pow(energy_solu[2], 2), 0.5);
                            Total_energy_air += energy_air * (jac_vect[i_gauss] / 6.0);
                            // Add the local Gauss contribution to the areas
                            total_area_positive += (jac_vect[i_gauss] / 6.0);
                            }
                    }
                

                    else if (this->IsNegative3D(elem_dist))
                    {
                        // Get geometry data
                        Vector jac_vect;
                        jac_vect = r_geom.DeterminantOfJacobian(jac_vect, GeometryData::GI_GAUSS_2);
                        auto N = r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

                        // Obtained solution
                        for (unsigned int i_gauss = 0; i_gauss < 4; ++i_gauss){
                            array_1d<double, 3> v_solu;
                            array_1d<double, 3> energy_solu;
                            const double p_solu = r_p_0 * N(i_gauss, 0) + r_p_1 * N(i_gauss, 1) + r_p_2 * N(i_gauss, 2) + r_p_3 * N(i_gauss, 3);
                            v_solu[0] = r_vel_0[0] * N(i_gauss, 0) + r_vel_1[0] * N(i_gauss, 1) + r_vel_2[0] * N(i_gauss, 2) + r_vel_3[0] * N(i_gauss, 3);
                            v_solu[1] = r_vel_0[1] * N(i_gauss, 0) + r_vel_1[1] * N(i_gauss, 1) + r_vel_2[1] * N(i_gauss, 2) + r_vel_3[1] * N(i_gauss, 3);
                            v_solu[2] = r_vel_0[2] * N(i_gauss, 0) + r_vel_1[2] * N(i_gauss, 1) + r_vel_2[2] * N(i_gauss, 2) + r_vel_3[2] * N(i_gauss, 3);

                            const double Water_Density = 1e03;
                            energy_solu[0] = ((std::pow(v_solu[0], 2) * Water_Density) / 2) + p_solu;
                            energy_solu[1] = ((std::pow(v_solu[1], 2) * Water_Density) / 2) + p_solu;
                            energy_solu[2] = ((std::pow(v_solu[2], 2) * Water_Density) / 2) + p_solu;
                            const double energy_water = std::pow(std::pow(energy_solu[0], 2) + std::pow(energy_solu[1], 2) + std::pow(energy_solu[2], 2), 0.5);
                            Total_energy_water += energy_water * (jac_vect[i_gauss] / 6.0);
                            // Add the local Gauss contribution to the areas
                            total_area_negative += (jac_vect[i_gauss] / 6.0);
                        }
                    }
                }   
        }
        }
        }
     
        // Print the obtained values
        std::cout << std::endl;
        std::cout << "Water_area_splited:" << cut_area_negative << std::endl;
        std::cout << "Air_area_splited" << cut_area_positive << std::endl;
        std::cout << "Water_area:" << total_area_negative << std::endl;
        std::cout << "Air_area:" << total_area_positive << std::endl;
        std::cout << "Total_energy_air_splitted:" << Total_energy_air_splitted << std::endl;
        std::cout << "Total_energy_water_splitted:" << Total_energy_water_splitted << std::endl;
        std::cout << "Total_energy_air:" << Total_energy_air << std::endl;
        std::cout << "Total_energy_water:" << Total_energy_water << std::endl;
        std::cout << std::endl;

        KRATOS_CATCH("");
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

    ///@}
    ///@name Friends
    ///@{


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

        for (unsigned int i_node = 0; i_node < 3; ++i_node){
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

        for (unsigned int i_node = 0; i_node < 3; ++i_node){
            if(rElementalDistances[i_node] > 0.0) {
                n_pos++;
            }
        }

        if (n_pos == 3){
            return true;
        } else {
            return false;
        }
    }

    bool IsNegative(const Vector &rElementalDistances)
    {
        unsigned int n_neg = 0;

        for (unsigned int i_node = 0; i_node < 3; ++i_node)
        {
            if (rElementalDistances[i_node] > 0.0)
            {
                n_neg++;
            }
        }

        if (n_neg == 3)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    bool IsCut3D(const Vector &rElementalDistances)
    {
        unsigned int n_pos = 0;
        unsigned int n_neg = 0;

        for (unsigned int i_node = 0; i_node < 4; ++i_node)
        {
            if (rElementalDistances[i_node] > 0.0)
            {
                n_pos++;
            }
            else
            {
                n_neg++;
            }
        }

        if (n_pos != 0 && n_neg != 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    bool IsPositive3D(const Vector &rElementalDistances)
    {
        unsigned int n_pos = 0;

        for (unsigned int i_node = 0; i_node < 4; ++i_node)
        {
            if (rElementalDistances[i_node] > 0.0)
            {
                n_pos++;
            }
        }

        if (n_pos == 4)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    bool IsNegative3D(const Vector &rElementalDistances)
    {
        unsigned int n_neg = 0;

        for (unsigned int i_node = 0; i_node < 4; ++i_node)
        {
            if (rElementalDistances[i_node] > 0.0)
            {
                n_neg++;
            }
        }

        if (n_neg == 4)
        {
            return true;
        }
        else
        {
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
inline std::istream &operator>>(std::istream &rIStream,
                                EnergyCheckProcess&rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const EnergyCheckProcess&rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GAUSS_POINT_ERROR_PROCESS_H_INCLUDED  defined

    © 2021 GitHub, Inc.
    Terms
    Privacy
    Security
    Status
    Docs

    Contact GitHub
    Pricing
    API
    Training
    Blog
    About

