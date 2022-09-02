//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/geometrical_object.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/hexahedra_3d_8.h"
#include "utilities/qef_utility.h"
#include "utilities/parallel_utilities.h"
#include "custom_modelers/voxel_mesh_generator_modeler.h"


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

/**
 * @class Dual countouring mesher
 * @ingroup KratosCore
 * @brief Utilities to re-mesh a shape using dual countouring
 */
class DualCountouringMesher
{
public:

    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;

    /// Pointer definition of VoxelInsideVolume
    KRATOS_CLASS_POINTER_DEFINITION( DualCountouringMesher );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    DualCountouringMesher(){}

    /// Destructor
    virtual ~DualCountouringMesher(){}

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Creates a mesh adapting to the shape of the rSkinPart passed
    * @param rSkinModelPart Geometrical Object containing rSkinPart and its information
    * @param rModel the model that will be containing the resulting mesh (the resulting mesh 
    * will be a modelPart named "fited_mesh")
    */
    void DualCountourAdaptativeRemesh(
        ModelPart& rSkinModelPart,
        ModelPart& rVoxelMesh,
        ModelPart& rFitedMesh) 
    {     
        array_1d<double,3> cell_size(3); 
        array_1d<double,3> min_bounding_box(3);
        array_1d<double,3> max_bounding_box(3);

        for(int i = 0; i < 3; i++) {
            std::vector<double> planes = GetKeyPlanes(i);
            cell_size[i] = planes[1] - planes[0];
            min_bounding_box[i] = planes[0];
            max_bounding_box[i] = planes[planes.size() -1];
        }
        
        //cell size tiene que venir de rVoxelMesh
        //modificar constructor per a desierd boundingbox
        GeometricalObjectsBins voxel_bin(rSkinModelPart.ElementsBegin(), rSkinModelPart.ElementsEnd(),cell_size,min_bounding_box,max_bounding_box);
        const auto& number_of_cells = voxel_bin.GetNumberOfCells();

        rFitedMesh.AddNodalSolutionStepVariable(DISTANCE); 
        //const auto& cell_size = voxel_bin.GetCellSizes();

        KRATOS_WATCH(voxel_bin.GetCellSizes());
        KRATOS_WATCH(voxel_bin.GetBoundingBox().GetMinPoint());
        KRATOS_WATCH(voxel_bin.GetBoundingBox().GetMaxPoint());

        //KRATOS_CHECK_EQUAL(voxel_mesh.Elements().size(),30*30*30);
        

        for (std::size_t i = 0; i < number_of_cells[0]; i++) {
            for (std::size_t j = 0; j < number_of_cells[1]; j++) {
                for (std::size_t k = 0; k < number_of_cells[2]; k++) {
                    BoundingBox<Point> box = voxel_bin.GetCellBoundingBox(i,j,k);
                    std::vector<GeometricalObject*> triangles =  voxel_bin.GetCell(i,j,k);
                    
                    int new_id = i + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1; 

                    array_1d<double,3> qef = QuadraticErrorFunction::QuadraticErrorFunctionPoint(box,triangles);
                    rFitedMesh.CreateNewNode(new_id, qef[0], qef[1], qef[2]);
                    rFitedMesh.pGetNode(new_id)->FastGetSolutionStepValue(DISTANCE) = 1;  //?                  
                }
            }
        }
        
        Properties::Pointer p_properties(new Properties(0)); 

        for (std::size_t i = 1; i < number_of_cells[0]; i++) {
            for (std::size_t j = 1; j < number_of_cells[1]; j++) {
                for (std::size_t k = 1; k < number_of_cells[2]; k++) {
                    rFitedMesh.CreateNewElement( "Element3D8N", 
                        i + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1,  //id
                        {(i-1) + (j-1) * number_of_cells[0] + (k-1) * number_of_cells[1] * number_of_cells[0] + 1,  //nodes
                        i + (j-1) * number_of_cells[0] + (k-1) * number_of_cells[1] * number_of_cells[0] + 1,
                        i + j * number_of_cells[0] + (k-1) * number_of_cells[1] * number_of_cells[0] + 1,
                        (i-1) + j * number_of_cells[0] + (k-1) * number_of_cells[1] * number_of_cells[0] + 1,
                        (i-1) + (j-1) * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1,
                        i + (j-1) * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1,
                        i + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1,
                        (i-1) + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1},
                        p_properties);      
                }
            }
        }

        KRATOS_WATCH('fin');

        /********************************************************************************************
        Properties::Pointer p_properties_1(new Properties(0)); 
        rVoxelMesh.CreateNewElement("Element3D8N", 8, 
            {id+1, id+2, id+3, id+4, id+5, id+6, id+7, id+8}, p_properties_1);    

        *******************************************************************************************/

    }

     void ThisIsTheRemix(
        ModelPart& rSkinModelPart,
        ModelPart& rVoxelMesh,
        ModelPart& rFitedMesh) 
    {   }

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    //This functions implementation should be moveed to VoxelUtilities
    double MeanDistance(GeometryType& rVoxel) {
        double mean = 0;
        int nodes_inside = 0;
        for (NodeType& node : rVoxel.Points()) {
            double dist = node.FastGetSolutionStepValue(DISTANCE);
            if (dist > 0) {
                nodes_inside++;
                mean += dist;
            }
        };

        if(nodes_inside > 0) return mean/nodes_inside;
        else return -1;
    }

}; /* Class DualCountouringMesher */

}  /* namespace Kratos.*/