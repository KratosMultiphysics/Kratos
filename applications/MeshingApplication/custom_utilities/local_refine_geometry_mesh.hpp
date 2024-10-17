// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:             BSD License
//                       Kratos default license: kratos/license.txt
//
//  Main authors:        Nelson Lafontaine
//                       Jordi Cotela Dalmau
//                       Riccardo Rossi
//  Co-authors:          Vicente Mataix Ferrandiz
//

#pragma once


// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

// Extrenal includes

// Project includes
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
#include "geometries/geometry.h"

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
class KRATOS_API(MESHING_APPLICATION) LocalRefineGeometryMesh
{
public:

    ///@name Type Definitions
    ///@{
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef std::vector<Matrix> Matrix_Order_Tensor;
    typedef std::vector<Vector> Vector_Order_Tensor;
    typedef std::vector<Vector_Order_Tensor> Node_Vector_Order_Tensor;
    typedef Node PointType;
    typedef Node ::Pointer PointPointerType;
    typedef std::vector<PointType::Pointer> PointVector;
    typedef PointVector::iterator PointIterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LocalRefineGeometryMesh(ModelPart& model_part) : mModelPart(model_part)
    {

    }

    /// Destructor
    virtual ~LocalRefineGeometryMesh()
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * Refine the mesh locally, call all the commands necessaries to compute the remeshing
    * @param refine_on_reference: Boolean that defines if refine or not the mesh according to the reference
    * @param interpolate_internal_variables: Boolean that defines if to interpolate or not the internal variables
    */

    virtual void LocalRefineMesh(
        bool refine_on_reference,
        bool interpolate_internal_variables
        );

    /**
    * This function initialises the matrix Cord
    * @return Coord: The matrix that stores all the index of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */

    virtual void CSRRowMatrix(
        ModelPart& this_model_part,
        compressed_matrix<int>& Coord
        );

    /**
    * This functions looks for potential edges that could be refined
    * @return Coord: The matrix that stores all the index of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */

    virtual void SearchEdgeToBeRefined(
        ModelPart& this_model_part,
        compressed_matrix<int>& Coord
        );

    /**
    * It creates the list of new nodes
    * @return this_model_part: The model part of the model (it is the input too)
    * @return Coord: The matrix that stores all the index of the geometry
    * @return List_New_Nodes: List that contents the index of the new nodes to be created
    * @return Position_Node: The vector that contents the position in the edge of the new nodes
    */

    virtual void CreateListOfNewNodes(
        ModelPart& this_model_part,
        compressed_matrix<int>& Coord,
        std::vector<int> &List_New_Nodes,
        std::vector<array_1d<int, 2 > >& Position_Node
        );

    /**
    * Computes the coordinates of the new nodes in the center of the edges
    * Insert the news nodes in the model part and interpolate the variables
    * @param List_New_Nodes: List that contents the index of the new nodes to be created
    * @param Position_Node: The vector that contents the position in the edge of the new nodes
    * @return this_model_part: The model part of the model (it is the input too)
    */

    virtual void CalculateCoordinateAndInsertNewNodes(
        ModelPart& this_model_part,
        const std::vector<array_1d<int, 2 > >& Position_Node,
        const std::vector<int> &List_New_Nodes
        );

    /**
    * It erases the old elements and it creates the new ones
    * @param Coord: The coordinates of the element
    * @param New_Elements: The new elements created
    * @param interpolate_internal_variables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return this_model_part: The model part of the model (it is the input too)
    */

    virtual void EraseOldElementAndCreateNewElement(
        ModelPart& this_model_part,
        const compressed_matrix<int>& Coord,
        PointerVector< Element >& New_Elements,
        bool interpolate_internal_variables
        );

    /**
    * Remove the old conditions and creates new ones
    * @param Coord: The coordinates of the nodes of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */

    virtual void EraseOldConditionsAndCreateNew(
        ModelPart& this_model_part,
        const compressed_matrix<int>& Coord
        );

    /**
    * It calculates the new edges of the new elements
    * first it calculates the new edges correspondign to the lower face (as a triangle),
    * later it added to the upper face
    * @param geom: The prism element geometry
    * @param edge_ids: The ids of the edges
    * @return aux: The vector that includes the index of the new edges
    */

    virtual void CalculateEdges(
        Element::GeometryType& geom,
        const compressed_matrix<int>& Coord,
        int* edge_ids,
        std::vector<int> & aux
        );

    /**
    * It creates a partition of the process between the different threads
    * @param number_of_threads: Number the threads considered in the computation
    * @param number_of_rows:
    * @return partitions: The vector that contents the partitions corresponding to each thread
    */

    inline void CreatePartition(
      unsigned int number_of_threads,
      const int number_of_rows,
      vector<unsigned int>& partitions
      );

    /**
    * Interpolates the internal variables
    * @param number_elem: Number of elements
    * @param father_elem: Father element (the original one)
    * @param child_elem: Child element (the new ones created)
    * @param rCurrentProcessInfo: The model part process info
    */

    template<typename TGeometricalObjectPointerType>
    void InterpolateInteralVariables(
            const int& number_elem,
            const TGeometricalObjectPointerType father_elem,
            TGeometricalObjectPointerType child_elem,
            const ProcessInfo& rCurrentProcessInfo
            )
    {
        // NOTE: Right now there is not an interpolation at all, it just copying the values
        std::vector<Vector> values;
        father_elem->CalculateOnIntegrationPoints(INTERNAL_VARIABLES, values, rCurrentProcessInfo);
        child_elem->SetValuesOnIntegrationPoints(INTERNAL_VARIABLES, values, rCurrentProcessInfo);
    }

    virtual void UpdateSubModelPartNodes(ModelPart &rModelPart);

    virtual void ResetFatherNodes(ModelPart &rModelPart);

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

        ModelPart& mModelPart;       /// The model part to be refined
        int mCurrentRefinementLevel; /// The current refinement level
        std::unordered_map<std::size_t, unsigned int> mMapNodeIdToPos;
        std::vector<std::size_t> mMapPosToNodeId;
    ///@}
    ///@name Protected Operators
    ///@{

    template<typename TIteratorType>
    void SearchEdgeToBeRefinedGeneric(
            TIteratorType GeometricalObjectsBegin,
            TIteratorType GeometricalObjectsEnd,
            compressed_matrix<int>& rCoord
					  )
    {
        KRATOS_TRY;

        for (TIteratorType it = GeometricalObjectsBegin; it != GeometricalObjectsEnd; ++it) {
            if (it->GetValue(SPLIT_ELEMENT)) {
                auto& r_geom = it->GetGeometry(); // Nodes of the element
                for (unsigned int i = 0; i < r_geom.size(); i++) {
                    int index_i = mMapNodeIdToPos[r_geom[i].Id()];
                    for (unsigned int j = 0; j < r_geom.size(); j++) {
                        int index_j = mMapNodeIdToPos[r_geom[j].Id()];
                        if (index_j > index_i)
                        {
                            rCoord(index_i, index_j) = -2;
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Protected Operations
    ///@{

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

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{
        
    ///@}
};

} // namespace Kratos.
