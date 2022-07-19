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
        
    void LocalRefineTet10Mesh(bool interpolate_internal_variables) {
            for (auto element : mModelPart.Elements()) element.SetValue(SPLIT_ELEMENT,true);
            LocalRefineMesh(false, interpolate_internal_variables);
        } 
  
    
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
        /**
    * It erases the old elements and it creates the new ones
    * @param Coord: The compressed matrix containing at (i,j) the id of the node created between nodes i,j
    * @param New_Elements: The new elements created
    * @param interpolate_internal_variables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void EraseOldElementAndCreateNewElement(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& NewElements,
            bool interpolate_internal_variables
    ) override
    {
	ElementsArrayType& rElements = this_model_part.Elements();
	ElementsArrayType::iterator it_begin = rElements.ptr_begin();
	ElementsArrayType::iterator it_end = rElements.ptr_end();
	unsigned int to_be_deleted = 0;
	unsigned int large_id = (rElements.end() - 1)->Id() * 15;
	unsigned int current_id = (rElements.end() - 1)->Id() + 1;
	bool create_element = false;
	int number_elem = 0;
	int splitted_edges = 0;
	int internal_node = 0;

	const ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();
	int edge_ids[6];
	int t[56];
	std::vector<int> aux;

        std::cout << "****************** REFINING MESH ******************" << std::endl;
        std::cout << "OLD NUMBER ELEMENTS: " << rElements.size() << std::endl;

	for (int & i : t)
	{
	    i = -1;
	}

	for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
	{
	  Element::GeometryType& geom = it->GetGeometry();
          CalculateEdges(geom, Coord, edge_ids, aux);

	  create_element = TetrahedraSplit::Split_Tetrahedra(edge_ids, t, &number_elem, &splitted_edges, &internal_node);

	  if (create_element == true)
	  {
	      to_be_deleted++;

          GlobalPointersVector< Element >& rChildElements = it->GetValue(NEIGHBOUR_ELEMENTS);
          // We will use this flag to identify the element later, when operating on
          // SubModelParts. Note that fully refined elements already have this flag set
          // to true, but this is not the case for partially refined element, so we set it here.
          it->SetValue(SPLIT_ELEMENT,true);
          rChildElements.resize(0);

	      // Create the new connectivity
	      
		 
		  unsigned int i0 = aux[t[0]];
		  unsigned int i1 = aux[t[1]];
		  unsigned int i2 = aux[t[2]];
		  unsigned int i3 = aux[t[3]];
          unsigned int i4 = aux[t[4]];
		  unsigned int i5 = aux[t[5]];
		  unsigned int i6 = aux[t[6]];
		  unsigned int i7 = aux[t[7]];
          unsigned int i8 = aux[t[8]];
		  unsigned int i9 = aux[t[9]];
          

		  Tetrahedra3D10<Node < 3 > > geom(
		      this_model_part.Nodes()(i0),
		      this_model_part.Nodes()(i1),
		      this_model_part.Nodes()(i2),
		      this_model_part.Nodes()(i3),
              this_model_part.Nodes()(i4),
		      this_model_part.Nodes()(i5),
		      this_model_part.Nodes()(i6),
		      this_model_part.Nodes()(i7),
              this_model_part.Nodes()(i8),
		      this_model_part.Nodes()(i9)
		  );

		  // Generate new element by cloning the base one
		  Element::Pointer p_element;
          const Element& rElem = KratosComponents<Element>::Get("Element3D10N");
		  p_element = rElem.Create(current_id, geom, it->pGetProperties());
		  p_element->Initialize(rCurrentProcessInfo);
		  p_element->InitializeSolutionStep(rCurrentProcessInfo);
		  p_element->FinalizeSolutionStep(rCurrentProcessInfo);

		  // Setting the internal variables in the child elem
		  if (interpolate_internal_variables == true)
		  {
		      InterpolateInteralVariables(number_elem, *it.base(), p_element, rCurrentProcessInfo);
		  }

		  // Transfer elemental variables
		  p_element->Data() = it->Data();
		  p_element->GetValue(SPLIT_ELEMENT) = false;
		  NewElements.push_back(p_element);

          rChildElements.push_back( Element::WeakPointer(p_element) );

		  current_id++;
	      it->SetId(large_id);
	      large_id++;
	  }

	  for (unsigned int i = 0; i < 32; i++)
	  {
	      t[i] = -1;
	  }
      }

      /* Adding news elements to the model part */
      for (PointerVector< Element >::iterator it_new = NewElements.begin(); it_new != NewElements.end(); it_new++)
      {
	  rElements.push_back(*(it_new.base()));
      }

      /* All of the elements to be erased are at the end */
      rElements.Sort();

      /* Now remove all of the "old" elements */
      rElements.erase(this_model_part.Elements().end() - to_be_deleted, this_model_part.Elements().end());

      std::cout << "NEW NUMBER ELEMENTS: " << rElements.size() << std::endl;


      // Now update the elements in SubModelParts
      if (NewElements.size() > 0)
      {
          for (ModelPart::SubModelPartIterator iSubModelPart = this_model_part.SubModelPartsBegin();
                  iSubModelPart != this_model_part.SubModelPartsEnd(); iSubModelPart++)
          {
              to_be_deleted = 0;
              NewElements.clear();

              // Create list of new elements in SubModelPart
              // Count how many elements will be removed
              for (ModelPart::ElementIterator iElem = iSubModelPart->ElementsBegin();
                      iElem != iSubModelPart->ElementsEnd(); iElem++)
              {
                  if( iElem->GetValue(SPLIT_ELEMENT) )
                  {
                      to_be_deleted++;
                      GlobalPointersVector< Element >& rChildElements = iElem->GetValue(NEIGHBOUR_ELEMENTS);

                      for ( auto iChild = rChildElements.ptr_begin();
                              iChild != rChildElements.ptr_end(); iChild++ )
                      {
                          NewElements.push_back((*iChild)->shared_from_this());
                      }
                  }
              }

              // Add new elements to SubModelPart
              iSubModelPart->Elements().reserve( iSubModelPart->Elements().size() + NewElements.size() );
              for (PointerVector< Element >::iterator it_new = NewElements.begin();
                      it_new != NewElements.end(); it_new++)
              {
                  iSubModelPart->Elements().push_back(*(it_new.base()));
              }

              // Delete old elements
              iSubModelPart->Elements().Sort();
              iSubModelPart->Elements().erase(iSubModelPart->Elements().end() - to_be_deleted, iSubModelPart->Elements().end());

              KRATOS_WATCH(iSubModelPart->Info());
              KRATOS_WATCH(to_be_deleted);
              KRATOS_WATCH(iSubModelPart->Elements().size());
              KRATOS_WATCH(this_model_part.Elements().size());
          }
      }
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
