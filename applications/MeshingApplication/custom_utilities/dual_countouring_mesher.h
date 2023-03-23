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

/** This class does not work in a suitable way. TO-DO:
 * The point resulting from calculating the qef is not constrained to the cell (it can be anywhere in space, leading to too 
 * deformated hexs). This should be solved in the kratos/ulitilies/qef_utility.cpp, adding some kind of restriction.
 * It should be checked if hexs have a real hex shape, so that the are valid. If they are not, qef point of that voxel should be 
 * moved a minimum distance or recalculated with penalizations so that the hex is valid. 
 *  
 * */


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
 * @brief Utilities to re-mesh a shape (skin) using dual countouring
 */
class DualCountouringMesher : public VoxelMeshGeneratorModeler
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

    /// Pointer definition of DualCountouringMesher
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

    //Constructor
    DualCountouringMesher(
        Model& rModel, Parameters rParameters = Parameters()) 
        : VoxelMeshGeneratorModeler(rModel, rParameters ) {

            bool cells_in_touch = false;
            bool cells_inside = false;
            for(auto parameters : rParameters["coloring_settings_list"]){
                cells_in_touch = cells_in_touch || (parameters["type"].GetString() == "cells_in_touch");
                cells_inside = cells_inside || (parameters["type"].GetString() == "cells_in_touch");
            }

            KRATOS_ERROR_IF(!cells_in_touch) 
                << "Dual Mesher Function only works including cells_in_touch in the colouring settings" << std::endl;

            KRATOS_ERROR_IF(!cells_inside) 
                << "Dual Mesher Function only works including cells_inside in the colouring settings" << std::endl;
        }

    /// Destructor
    virtual ~DualCountouringMesher(){}

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Creates a mesh adapting to the shape of the rSkinPart passed
    * @param rFitedMesh The resulting modelPart containing the mesh adapted to the surface
    * @note How the function works
    * First of all we get the cell_size and BoundingBox from the voxel mesh that th modeler created. From these, we create
    * a GeometricalObjectsBins (with found cellsize and BoundingBox), so that we have acces to the SkinModelPart elements
    * that intersect each of the voxels of the voxel mesh. From this intersections, the QEF point is calculated and a new 
    * node is created with ID corresponding to the cell position. Finally, these new nodes are connected into new elements 
    * in such a way that the new element is only created if all the found qef points come from cells that were either in 
    * touch or inside the surface of the volume we wanted to fit.
    */
    void DualCountourAdaptativeRemesh(
        ModelPart& rFitedMesh) 
    {     
        KRATOS_ERROR_IF(!rFitedMesh.HasNodalSolutionStepVariable(DISTANCE)) 
            << "Input mesh is expected to have DISTANCE as nodal variable" << std::endl;
        
        array_1d<double,3> cell_size(3); 
        array_1d<double,3> min_bounding_box(3);
        array_1d<double,3> max_bounding_box(3);

        for(int i = 0; i < 3; i++) {
            std::vector<double> planes = GetKeyPlanes(i);
            cell_size[i] = planes[1] - planes[0];
            min_bounding_box[i] = planes[0];
            max_bounding_box[i] = planes[planes.size() -1];
        }
        
        GeometricalObjectsBins voxel_bin(mpInputModelPart->ElementsBegin(), mpInputModelPart->ElementsEnd(),cell_size,min_bounding_box,max_bounding_box);

        const auto& number_of_cells = voxel_bin.GetNumberOfCells();
        
        for (std::size_t i = 0; i < number_of_cells[0]; i++) {
            for (std::size_t j = 0; j < number_of_cells[1]; j++) {
                for (std::size_t k = 0; k < number_of_cells[2]; k++) {
                    std::vector<GeometricalObject*> triangles =  voxel_bin.GetCell(i,j,k);
                    int new_id = i + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1; 

                    array_1d<double,3> qef;
                    if (triangles.size() == 0) {
                        //no triangles --> qef = center 
                        qef(0) = min_bounding_box[0] + (i + 0.5)* cell_size(0);
                        qef(1) = min_bounding_box[1] + (j + 0.5)* cell_size(1);
                        qef(2) = min_bounding_box[2] + (k + 0.5)* cell_size(2);
                    } else {
                        BoundingBox<Point> box = voxel_bin.GetCellBoundingBox(i,j,k);
                        qef = QuadraticErrorFunction::QuadraticErrorFunctionPoint(box,triangles); 
                    }

                    rFitedMesh.CreateNewNode(new_id, qef[0], qef[1], qef[2]);
                    if (mIsInside[new_id]) {
                        rFitedMesh.pGetNode(new_id)->FastGetSolutionStepValue(DISTANCE) = 1;                   
                    } else {
                        rFitedMesh.pGetNode(new_id)->FastGetSolutionStepValue(DISTANCE) = 0;
                    }
                }
            }
        }
        
        Properties::Pointer p_properties(new Properties(0)); 

        for (std::size_t i = 1; i < number_of_cells[0]; i++) {
            for (std::size_t j = 1; j < number_of_cells[1]; j++) {
                for (std::size_t k = 1; k < number_of_cells[2]; k++) {
                    if (CheckAllNeighbourQEFAreInside(i,j,k)) {
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
                             
                            //Here it should be checked that the created element is valid 
                    }     
                }
            }
        }
    }

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

        std::vector<int> mIsInside; //bools
        double mInsideColor;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    bool CheckAllNeighbourQEFAreInside(int i, int j, int k) {
        bool res = (mColors.GetElementalColor(i-1,j-1,k-1) == mInsideColor) &&
                    (mColors.GetElementalColor(i,j-1,k-1) == mInsideColor) &&
                    (mColors.GetElementalColor(i,j,k-1) == mInsideColor) &&
                    (mColors.GetElementalColor(i-1,j,k-1) == mInsideColor) &&
                    (mColors.GetElementalColor(i-1,j-1,k) == mInsideColor) &&
                    (mColors.GetElementalColor(i,j-1,k) == mInsideColor) &&
                    (mColors.GetElementalColor(i,j,k) == mInsideColor) &&
                    (mColors.GetElementalColor(i-1,j,k) == mInsideColor);
        return res;
    }

    void GenerateElementsWithCellColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters) override {
        double inside_color = EntityGeneratorParameters["color"].GetDouble();
        mInsideColor = inside_color;
        std::size_t properties_id = EntityGeneratorParameters["properties_id"].GetInt();

        if(!rTheVolumeModelPart.HasProperties(properties_id)){
            rTheVolumeModelPart.CreateNewProperties(properties_id);
        }    
        Properties::Pointer p_properties = rTheVolumeModelPart.pGetProperties(properties_id);

        std::size_t cell_index = 0;
        array_1d<std::size_t, 3> number_of_cells;
        for(int i = 0 ; i < 3 ; i++){
            number_of_cells[i]=mKeyPlanes[i].size() - 1;
        }

        ModelPart::NodesContainerType new_nodes;
        ModelPart::ElementsContainerType new_elements;
        auto& r_prototype_element = KratosComponents<Element>::Get("Element3D8N");
            
        mIsInside.resize(number_of_cells[0]*number_of_cells[1]*number_of_cells[2]);

        Element::NodesArrayType cell_nodes(8);
        for (std::size_t k = 0; k < number_of_cells[2]; k++) {
            for (std::size_t j = 0; j < number_of_cells[1]; j++) {
                for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                    if(mColors.GetElementalColor(i,j,k) == inside_color){
                        cell_nodes(0) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j  , k);
                        cell_nodes(1) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j  , k);
                        cell_nodes(2) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j+1, k);
                        cell_nodes(3) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j+1, k);
                        cell_nodes(4) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j  , k+1);
                        cell_nodes(5) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j  , k+1);
                        cell_nodes(6) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j+1, k+1);
                        cell_nodes(7) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j+1, k+1);

                        //create the new element
                        Element::Pointer p_element = r_prototype_element.Create(mStartElementId + cell_index, cell_nodes, p_properties);
                        new_elements.push_back(p_element);
                        cell_index++;

                        mIsInside[i + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1] = 1;
                    } else {
                        mIsInside[i + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1] = 0;
                    }
                }
            }
        }
    
        rTheVolumeModelPart.AddNodes(new_nodes.begin(), new_nodes.end());
        rTheVolumeModelPart.AddElements(new_elements.begin(), new_elements.end());
   }


}; /* Class DualCountouringMesher */

}  /* namespace Kratos.*/
