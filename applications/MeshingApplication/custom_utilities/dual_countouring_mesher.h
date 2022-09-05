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

    //Constructor
    DualCountouringMesher(
        Model& rModel, Parameters rParameters = Parameters()) 
        : VoxelMeshGeneratorModeler(rModel, rParameters ) {}

    /// Destructor
    virtual ~DualCountouringMesher(){}

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Creates a mesh adapting to the shape of the rSkinPart passed
    * @param rFitedMesh The resulting modelPart containing the mesh adapted to the surface
    */
    void DualCountourAdaptativeRemesh(
        ModelPart& rFitedMesh) 
    {     
        rFitedMesh.AddNodalSolutionStepVariable(DISTANCE); //This should be smth like Kratos_error_if_not(rFitedMesh.Has(DISTANCE))
        
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
        
        //KRATOS_WATCH(number_of_cells);

        for (std::size_t i = 0; i < number_of_cells[0]; i++) {
            for (std::size_t j = 0; j < number_of_cells[1]; j++) {
                for (std::size_t k = 0; k < number_of_cells[2]; k++) {
                    BoundingBox<Point> box = voxel_bin.GetCellBoundingBox(i,j,k);
                    std::vector<GeometricalObject*> triangles =  voxel_bin.GetCell(i,j,k);
                    int new_id = i + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1; 

                    array_1d<double,3> qef = QuadraticErrorFunction::QuadraticErrorFunctionPoint(box,triangles);
                    KRATOS_WATCH(qef);
                    rFitedMesh.CreateNewNode(new_id, qef[0], qef[1], qef[2]);
                    if (mIsInside[i + j * number_of_cells[0] + k * number_of_cells[1] * number_of_cells[0] + 1]) {
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

        std::vector<int> mIsInside; //bools

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

    void GenerateElementsWithCellColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters) override {
         double inside_color = EntityGeneratorParameters["color"].GetDouble();
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

        //KRATOS_WATCH(mIsInside);
    
        rTheVolumeModelPart.AddNodes(new_nodes.begin(), new_nodes.end());
        rTheVolumeModelPart.AddElements(new_elements.begin(), new_elements.end());
   }


}; /* Class DualCountouringMesher */

}  /* namespace Kratos.*/