//    |  /           |
//    ' /   __| _` | __|  _ \   __|w
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Montanino
//

#pragma once

// System includes
#include <fstream>
#include <iostream>
#include <cmath>

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"

// Application includes
#include "shallow_water_application_variables.h"

namespace Kratos
{

class InterpolateSwToPfemUtility
{

protected:

/// Basic Structs for the utility ---------------------------------------------------------------------------------------------------------------------------------------------

    struct UtilityVariables
    {
        double X_max, X_min, Y_max, Y_min, Z_max, Z_min;
        int NRows, NColumns, NSections;
        double RowSize, ColumnSize, SectionSize;
    };

public:

    KRATOS_CLASS_POINTER_DEFINITION( InterpolateSwToPfemUtility );

    /// Constructor
    InterpolateSwToPfemUtility() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~InterpolateSwToPfemUtility() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InterpolateVariables (ModelPart& rPfemModelPart, ModelPart& rSWModelPart)
    {
        
        const std::size_t NNodesSW = static_cast<int>(rSWModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_beginSW = rSWModelPart.NodesBegin();
        
        const std::size_t NNodesPfem = static_cast<int>(rPfemModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_beginPfem = rPfemModelPart.NodesBegin();

        double distance, distance_to_first_neigh_node, distance_to_second_neigh_node;
        double ratio, ratio1,ratio2;
        double tol = 1e-4;
        double Tol = 5.0;
        double Tol2 = Tol*Tol;
                
        int     first_neigh_node, second_neigh_node;
        
        double             height, vertical_velocity, topography;
        array_1d<double,3> velocity;

        for(std::size_t i = 0; i < NNodesPfem; i++) {
            ModelPart::NodesContainerType::iterator it_node = node_beginPfem + i;

            distance_to_first_neigh_node = 1e32;
            distance_to_second_neigh_node = 1e32;

            for (std::size_t j = 0; j < NNodesSW; j++){
            ModelPart::NodesContainerType::iterator jt_node = node_beginSW + j;

                distance = (it_node->X() - jt_node->X())*(it_node->X() - jt_node->X()) + (it_node->Y() - jt_node->Y())*(it_node->Y() - jt_node->Y());
                if (distance < distance_to_first_neigh_node){
                       distance_to_first_neigh_node = distance;
                       first_neigh_node =  j;
                }
            }

            for (std::size_t j = 0; j < NNodesSW; j++){
            ModelPart::NodesContainerType::iterator jt_node = node_beginSW + j;

                distance = (it_node->X() - jt_node->X())*(it_node->X() - jt_node->X()) + (it_node->Y() - jt_node->Y())*(it_node->Y() - jt_node->Y());
                if ((distance < distance_to_second_neigh_node) & (distance > distance_to_first_neigh_node)){
                       distance_to_second_neigh_node = distance;
                       second_neigh_node =  j;
                }
            }
            
            ModelPart::NodesContainerType::iterator jt_node1 = node_beginSW + first_neigh_node;
            ModelPart::NodesContainerType::iterator jt_node2 = node_beginSW + second_neigh_node;

            double dx = std::abs(jt_node1->X() - jt_node2->X());
            double dy = std::abs(jt_node1->Y() - jt_node2->Y());
            
            if (distance_to_first_neigh_node*distance_to_first_neigh_node > Tol2*(dx*dx + dy*dy)){
                KRATOS_WARNING("Interpolate_sw_to_pfem") << "PFEM node quite far from SW output interface " << std::endl; 
            }

            if (((dx > tol) & (dy < tol)) | ((dx < tol) & (dy > tol))){
                if ((dx > tol) & (dy < tol)){
                    ratio = (it_node->X() - jt_node2->X())/(jt_node1->X() - jt_node2->X());
                } else if ((dx < tol) & (dy > tol)){
                    ratio = (it_node->Y() - jt_node2->Y())/(jt_node1->Y() - jt_node2->Y());
                }

                velocity = ratio*(jt_node1->GetValue(VELOCITY) - jt_node2->GetValue(VELOCITY)) + jt_node2->GetValue(VELOCITY);
                height   = ratio*(jt_node1->GetValue(HEIGHT) - jt_node2->GetValue(HEIGHT)) + jt_node2->GetValue(HEIGHT); 
                vertical_velocity = ratio*(jt_node1->GetValue(VERTICAL_VELOCITY) - jt_node2->GetValue(VERTICAL_VELOCITY)) + jt_node2->GetValue(VERTICAL_VELOCITY); 
                topography = ratio*(jt_node1->GetValue(TOPOGRAPHY) - jt_node2->GetValue(TOPOGRAPHY)) + jt_node2->GetValue(TOPOGRAPHY); 
            } else {
                ratio1 = (it_node->X() - jt_node2->X())/(jt_node1->X() - jt_node2->X());
                ratio2 = (it_node->Y() - jt_node2->Y())/(jt_node1->Y() - jt_node2->Y());

                velocity = ratio1*(jt_node1->GetValue(VELOCITY) - jt_node2->GetValue(VELOCITY)) + jt_node2->GetValue(VELOCITY);
                height   = ratio1*(jt_node1->GetValue(HEIGHT) - jt_node2->GetValue(HEIGHT)) + jt_node2->GetValue(HEIGHT); 
                vertical_velocity = ratio1*(jt_node1->GetValue(VERTICAL_VELOCITY) - jt_node2->GetValue(VERTICAL_VELOCITY)) + jt_node2->GetValue(VERTICAL_VELOCITY); 
                topography = ratio1*(jt_node1->GetValue(TOPOGRAPHY) - jt_node2->GetValue(TOPOGRAPHY)) + jt_node2->GetValue(TOPOGRAPHY); 

                velocity        += ratio2*(jt_node1->GetValue(VELOCITY) - jt_node2->GetValue(VELOCITY)) + jt_node2->GetValue(VELOCITY);
                height          += ratio2*(jt_node1->GetValue(HEIGHT) - jt_node2->GetValue(HEIGHT)) + jt_node2->GetValue(HEIGHT); 
                vertical_velocity += ratio2*(jt_node1->GetValue(VERTICAL_VELOCITY) - jt_node2->GetValue(VERTICAL_VELOCITY)) + jt_node2->GetValue(VERTICAL_VELOCITY); 
                topography      += ratio2*(jt_node1->GetValue(TOPOGRAPHY) - jt_node2->GetValue(TOPOGRAPHY)) + jt_node2->GetValue(TOPOGRAPHY);

                velocity            = 0.5*velocity;
                height              = 0.5*height;
                vertical_velocity   = 0.5*vertical_velocity;
                topography          = 0.5*topography;
            }

            it_node->GetValue(HEIGHT) = height;
            it_node->GetValue(VELOCITY) = velocity;
            it_node->GetValue(VERTICAL_VELOCITY) = vertical_velocity;
            it_node->GetValue(TOPOGRAPHY) = topography;

        }
    };

protected:

    /// Member Variables



///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
/// Common --------------------------------------------------------------------------------------------------------------------------------------------------------------------------



}; // Class InterpolateSwToPfemUtility

} // namespace Kratos.


