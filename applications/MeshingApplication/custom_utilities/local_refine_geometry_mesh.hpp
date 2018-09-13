// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//    Co-authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_LOCAL_REFINE_GEOMETRY_MESH)
#define  KRATOS_LOCAL_REFINE_GEOMETRY_MESH
#ifdef _OPENMP
#include <omp.h>
#endif

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

// Extrenal includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/dof.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/constitutive_law.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "processes/node_erase_process.h"

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
class LocalRefineGeometryMesh
{
public:

    ///@name Type Definitions
    ///@{
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef boost::numeric::ublas::vector<Matrix> Matrix_Order_Tensor;
    typedef boost::numeric::ublas::vector<Vector> Vector_Order_Tensor;
    typedef boost::numeric::ublas::vector<Vector_Order_Tensor> Node_Vector_Order_Tensor;
    typedef Node < 3 > PointType;
    typedef Node < 3 > ::Pointer PointPointerType;
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
    ~LocalRefineGeometryMesh()
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

    void LocalRefineMesh(
            bool refine_on_reference,
            bool interpolate_internal_variables);

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
            boost::numeric::ublas::vector<int> &List_New_Nodes,
            boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node
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
            const boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node,
            const boost::numeric::ublas::vector<int> &List_New_Nodes
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
    * This process renumerates the elements and nodes
    * @param New_Elements: Pointers to the new elements created
    * @return this_model_part: The model part of the model (it is the input too)
    */

    virtual void RenumeringElementsAndNodes(
            ModelPart& this_model_part,
            PointerVector< Element >& New_Elements
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

    virtual void InterpolateInteralVariables(
            const int& number_elem,
            const Element::Pointer father_elem,
            Element::Pointer child_elem,
            ProcessInfo& rCurrentProcessInfo
            );

    virtual void UpdateSubModelPartNodes(ModelPart &rModelPart);

    virtual void ResetFatherNodes(ModelPart &rModelPart);

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

    ModelPart& mModelPart;

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

#endif // KRATOS_LOCAL_REFINE_GEOMETRY_MESH  defined
