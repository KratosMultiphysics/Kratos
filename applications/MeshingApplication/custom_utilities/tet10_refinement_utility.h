// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//

#if !defined(KRATOS_TET10_REFINEMENT_UTILITY)
#define  KRATOS_TET10_REFINEMENT_UTILITY

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

/* Project includes */
#include "includes/node.h"
#include "custom_utilities/local_refine_tetrahedra_mesh.hpp"
#include "containers/model.h"
#include "includes/element.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_10.h"

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
class Tet10RefinementUtility : public LocalRefineTetrahedraMesh
{
public:

    ///@name Type Definitions
    ///@{
    ///@}

    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;

    ///@name Life Cycle
    ///@{

    /// Default constructors
    Tet10RefinementUtility(ModelPart& model_part) : LocalRefineTetrahedraMesh(model_part)
    {

    }

    /// Destructor
    ~Tet10RefinementUtility() //TODO maybe {}
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
        /*
    void LocalRefineMesh(
        bool refine_on_reference,
        bool interpolate_internal_variables) {
            //for (auto element : mModelPart.Elements()) element->SetValue(SPLIT_ELEMENT,true);
            //super.LocalRefineMesh();
        } */
  
    
protected:
    ///@name Protected static Member Variables
    ///@{
    //int mPreviousRefinementLevel;
    //int mCurrentRefinementLevel;
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

     void EraseOldElementAndCreateNewElement(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& NewElements,
            bool interpolate_internal_variables
    ) override
    {
	ElementsArrayType& rElements = this_model_part.Elements();
    unsigned int CurrentId = (rElements.end() -1)->Id();
    PointsArrayType newNodes;

    for (ElementsArrayType::iterator it = rElements.begin(); it != rElements.end(); it++) {
        Properties::Pointer p_properties_1 = it->pGetProperties();
        auto nodes = it->pGetGeometry()->Points();
        newNodes.push_back(&nodes[0]);
        newNodes.push_back(&nodes[1]);
        newNodes.push_back(&nodes[2]);
        newNodes.push_back(&nodes[3]);

        getNewNodes(Coord, newNodes);

        Tetrahedra3D10<NodeType> geometry(newNodes);
        
    }
    }

    void getNewNodes(const compressed_matrix<int>& Coord, PointsArrayType& rNewNodes) {

    }

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

#endif // KRATOS_TET10_REFINEMENT_UTILITY  defined
