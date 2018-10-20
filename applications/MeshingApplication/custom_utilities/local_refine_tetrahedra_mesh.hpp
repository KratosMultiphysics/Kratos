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

#if !defined(KRATOS_LOCAL_REFINE_TETRAHEDRA_MESH)
#define  KRATOS_LOCAL_REFINE_TETRAHEDRA_MESH

#ifdef _OPENMP
#include <omp.h>
#endif

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

/* Project includes */
#include "geometries/tetrahedra_3d_4.h"
#include "custom_utilities/local_refine_geometry_mesh.hpp"
#include "utilities/split_tetrahedra.h"
#include "utilities/split_triangle.c"

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
class LocalRefineTetrahedraMesh : public LocalRefineGeometryMesh
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LocalRefineTetrahedraMesh(ModelPart& model_part) : LocalRefineGeometryMesh(model_part)
    {

    }

    /// Destructor
    ~LocalRefineTetrahedraMesh()
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * Insert the news nodes in the center of elements and interopolate the variables.
    * @param geom: The curret geometry
    * @return model_part: The model part of the model (it is the input too)
    */

    unsigned int CreateCenterNode(Geometry<Node < 3 > >& geom, ModelPart& model_part)
    {
        // Determine a new unique id
        unsigned int new_id = (model_part.NodesEnd() - 1)->Id() + 1;

        if( model_part.Nodes().find(new_id) != model_part.NodesEnd() )
	{
            KRATOS_THROW_ERROR(std::logic_error, "adding a center node with an already existing id","");
	}

	// Determine the coordinates of the new node
	double X = (geom[0].X() + geom[1].X() + geom[2].X() + geom[3].X()) / 4.0;
        double Y = (geom[0].Y() + geom[1].Y() + geom[2].Y() + geom[3].Y()) / 4.0;
        double Z = (geom[0].Z() + geom[1].Z() + geom[2].Z() + geom[3].Z()) / 4.0;

        double X0 = (geom[0].X0() + geom[1].X0() + geom[2].X0() + geom[3].X0()) / 4.0;
        double Y0 = (geom[0].Y0() + geom[1].Y0() + geom[2].Y0() + geom[3].Y0()) / 4.0;
        double Z0 = (geom[0].Z0() + geom[1].Z0() + geom[2].Z0() + geom[3].Z0()) / 4.0;

        // Generate the new node
        Node < 3 > ::Pointer pnode = model_part.CreateNewNode(new_id, X, Y, Z);

        unsigned int buffer_size = model_part.NodesBegin()->GetBufferSize();
        pnode->SetBufferSize(buffer_size);

        pnode->X0() = X0;
        pnode->Y0() = Y0;
        pnode->Z0() = Z0;

        // Add the dofs
        Node < 3 > ::DofsContainerType& reference_dofs = (model_part.NodesBegin())->GetDofs();
        unsigned int step_data_size = model_part.GetNodalSolutionStepDataSize();

        for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
        {
            Node < 3 > ::DofType& rDof = *iii;
            Node < 3 > ::DofType::Pointer p_new_dof = pnode->pAddDof(rDof);

            // The variables are left as free for the internal node
            (p_new_dof)->FreeDof();
        }

        // Intepolating the data
        for (unsigned int step = 0; step < buffer_size; step++)
        {
            double* new_step_data = pnode->SolutionStepData().Data(step);
            double* step_data1 = geom[0].SolutionStepData().Data(step);
            double* step_data2 = geom[1].SolutionStepData().Data(step);
            double* step_data3 = geom[2].SolutionStepData().Data(step);
            double* step_data4 = geom[3].SolutionStepData().Data(step);

            // Copying this data in the position of the vector we are interested in
            for (unsigned int j = 0; j < step_data_size; j++)
            {
                new_step_data[j] = 0.25 * (step_data1[j] + step_data2[j] + step_data3[j] + step_data4[j]);
            }
        }

        pnode->GetValue(FATHER_NODES).resize(0);
        pnode->GetValue(FATHER_NODES).push_back( Node< 3 >::WeakPointer( geom(0) ) );
        pnode->GetValue(FATHER_NODES).push_back( Node< 3 >::WeakPointer( geom(1) ) );
        pnode->GetValue(FATHER_NODES).push_back( Node< 3 >::WeakPointer( geom(2) ) );
        pnode->GetValue(FATHER_NODES).push_back( Node< 3 >::WeakPointer( geom(3) ) );

        return new_id;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It erases the old elements and it creates the new ones
    * @param Coord: The coordinates of the element
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

	ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();
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

	  if (internal_node == 1)
	  {
	      // Generate new internal node
	      aux[10] = CreateCenterNode(geom, this_model_part);

	      bool verified = false;
	      for(int iii = 0; iii < number_elem*4; iii++)
	      {
		  if(t[iii] == 10)
		  {
		      verified = true;
		  }
	      }

	      if(verified == false)
	      {
		  KRATOS_WATCH(number_elem);
		  for(int iii = 0; iii < number_elem * 4; iii++)
		  {
		      std::cout << t[iii] << std::endl;
		  }

		  KRATOS_THROW_ERROR(std::logic_error,"internal node is created but not used","");
	      }

	  }

	  if (create_element == true)
	  {
	      to_be_deleted++;

          WeakPointerVector< Element >& rChildElements = it->GetValue(NEIGHBOUR_ELEMENTS);
          // We will use this flag to identify the element later, when operating on
          // SubModelParts. Note that fully refined elements already have this flag set
          // to true, but this is not the case for partially refined element, so we set it here.
          it->SetValue(SPLIT_ELEMENT,true);
          rChildElements.resize(0);

	      // Create the new connectivity
	      for (int i = 0; i < number_elem; i++)
	      {
		  unsigned int base = i * 4;
		  unsigned int i0 = aux[t[base]];
		  unsigned int i1 = aux[t[base + 1]];
		  unsigned int i2 = aux[t[base + 2]];
		  unsigned int i3 = aux[t[base + 3]];

		  Tetrahedra3D4<Node < 3 > > geom(
		      this_model_part.Nodes()(i0),
		      this_model_part.Nodes()(i1),
		      this_model_part.Nodes()(i2),
		      this_model_part.Nodes()(i3)
		  );

		  // Generate new element by cloning the base one
		  Element::Pointer p_element;
		  p_element = it->Create(current_id, geom, it->pGetProperties());
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
	      }
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
                      WeakPointerVector< Element >& rChildElements = iElem->GetValue(NEIGHBOUR_ELEMENTS);

                      for ( auto iChild = rChildElements.ptr_begin();
                              iChild != rChildElements.ptr_end(); iChild++ )
                      {
                          NewElements.push_back( (*iChild).lock() );
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

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Remove the old conditions and creates new ones
    * @param Coord: The coordinates of the nodes of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void EraseOldConditionsAndCreateNew(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord
            ) override
    {
        KRATOS_TRY;

        PointerVector< Condition > NewConditions;

        ConditionsArrayType& rConditions = this_model_part.Conditions();

        if(rConditions.size() > 0)
        {
            ConditionsArrayType::iterator it_begin = rConditions.ptr_begin();
            ConditionsArrayType::iterator it_end = rConditions.ptr_end();
            unsigned int to_be_deleted = 0;
            unsigned int large_id = (rConditions.end() - 1)->Id() * 7;
            int  edge_ids[3];
            int  t[12];
            int  number_elem             = 0;
            int  splitted_edges  = 0;
            int  nint            = 0;
            array_1d<int, 6> aux;

            ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();

            unsigned int current_id = (rConditions.end() - 1)->Id() + 1;
            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                Condition::GeometryType& geom = it->GetGeometry();

                if (geom.size() == 3)
                {
                    CalculateEdgesFaces(geom, Coord, edge_ids, aux);

                    // Create the new conditions
                    bool create_condition =  Split_Triangle(edge_ids, t, &number_elem, &splitted_edges, &nint);

                    if(create_condition==true)
                    {
                        WeakPointerVector< Condition >& rChildConditions = it->GetValue(NEIGHBOUR_CONDITIONS);
                        // We will use this flag to identify the condition later, when operating on
                        // SubModelParts.
                        it->SetValue(SPLIT_ELEMENT,true);
                        rChildConditions.resize(0);
                        to_be_deleted++;

                        for(int i = 0; i < number_elem; i++)
                        {
                            unsigned int base = i * 3;
                            unsigned int i0   = aux[t[base]];
                            unsigned int i1   = aux[t[base+1]];
                            unsigned int i2   = aux[t[base+2]];

                            Triangle3D3<Node<3> > newgeom(
                                    this_model_part.Nodes()(i0),
                                    this_model_part.Nodes()(i1),
                                    this_model_part.Nodes()(i2)
                                    );

                            // Generate new condition by cloning the base one
                            Condition::Pointer pcond;
                            pcond = it->Create(current_id, newgeom, it->pGetProperties());
                            pcond ->Initialize(rCurrentProcessInfo);
                            pcond ->InitializeSolutionStep(rCurrentProcessInfo);
                            pcond ->FinalizeSolutionStep(rCurrentProcessInfo);

                            // Transfer condition variables
                            pcond->Data() = it->Data();
                            pcond->GetValue(SPLIT_ELEMENT) = false;
                            NewConditions.push_back(pcond);

                            rChildConditions.push_back( Condition::WeakPointer( pcond ) );

                            current_id++;
                        }
                        it->SetId(large_id);
                        large_id++;
                    }
                }
            }

            /* All of the conditions to be erased are at the end */
            this_model_part.Conditions().Sort();

            /* Now remove all of the "old" conditions*/
            this_model_part.Conditions().erase(this_model_part.Conditions().end() - to_be_deleted, this_model_part.Conditions().end());

            unsigned int total_size = this_model_part.Conditions().size()+ NewConditions.size();
            this_model_part.Conditions().reserve(total_size);

            /// Add the new Conditions to the ModelPart
            for (auto iCond = NewConditions.ptr_begin();
                    iCond != NewConditions.ptr_end(); iCond++)
            {
                this_model_part.Conditions().push_back( *iCond );
            }

            /* Renumber id */
            unsigned int my_index = 1;
            for(ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); it++)
            {
                it->SetId(my_index++);
            }


            // Now update the conditions in SubModelParts
            if (NewConditions.size() > 0)
            {
                for (ModelPart::SubModelPartIterator iSubModelPart = this_model_part.SubModelPartsBegin();
                        iSubModelPart != this_model_part.SubModelPartsEnd(); iSubModelPart++)
                {
                    to_be_deleted = 0;
                    NewConditions.clear();

                    // Create list of new conditions in SubModelPart
                    // Count how many conditions will be removed
                    for (ModelPart::ConditionIterator iCond = iSubModelPart->ConditionsBegin();
                            iCond != iSubModelPart->ConditionsEnd(); iCond++)
                    {
                        if( iCond->GetValue(SPLIT_ELEMENT) )
                        {
                            to_be_deleted++;
                            WeakPointerVector< Condition >& rChildConditions = iCond->GetValue(NEIGHBOUR_CONDITIONS);

                            for ( auto iChild = rChildConditions.ptr_begin();
                                    iChild != rChildConditions.ptr_end(); iChild++ )
                            {
                                NewConditions.push_back( (*iChild).lock() );
                            }
                        }
                    }

                    // Add new conditions to SubModelPart
                    iSubModelPart->Conditions().reserve( iSubModelPart->Conditions().size() + NewConditions.size() );
                    for (PointerVector< Condition >::iterator it_new = NewConditions.begin();
                            it_new != NewConditions.end(); it_new++)
                    {
                        iSubModelPart->Conditions().push_back(*(it_new.base()));
                    }

                    // Delete old conditions
                    iSubModelPart->Conditions().Sort();
                    iSubModelPart->Conditions().erase(iSubModelPart->Conditions().end() - to_be_deleted, iSubModelPart->Conditions().end());

                    KRATOS_WATCH(iSubModelPart->Info());
                    KRATOS_WATCH(to_be_deleted);
                    KRATOS_WATCH(iSubModelPart->Conditions().size());
                    KRATOS_WATCH(this_model_part.Conditions().size());
                }
            }
        }
        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It calculates the new edges of the new tetrahedras,
    * first it calculates the new edges correspondign to the lower face (as a triangle),
    * later it added to the upper face
    * @param geom: The tetrahedra element geometry
    * @param edge_ids: The ids of the edges
    * @return aux: The vector that includes the index of the new edges
    */

    void CalculateEdges(
            Element::GeometryType& geom,
            const compressed_matrix<int>& Coord,
            int* edge_ids,
            std::vector<int> & aux
            ) override
    {
      aux.resize(11, false);

      int index_0 = geom[0].Id() - 1;
      int index_1 = geom[1].Id() - 1;
      int index_2 = geom[2].Id() - 1;
      int index_3 = geom[3].Id() - 1;

      // Put the global ids in aux
      aux[0] = geom[0].Id();
      aux[1] = geom[1].Id();
      aux[2] = geom[2].Id();
      aux[3] = geom[3].Id();

      if (index_0 > index_1)
      {
	  aux[4] = Coord(index_1, index_0);
      }
      else
      {
	  aux[4] = Coord(index_0, index_1);
      }

      if (index_0 > index_2)
      {
	  aux[5] = Coord(index_2, index_0);
      }
      else
      {
	  aux[5] = Coord(index_0, index_2);
      }

      if (index_0 > index_3)
      {
	  aux[6] = Coord(index_3, index_0);
      }
      else
      {
	  aux[6] = Coord(index_0, index_3);
      }

      if (index_1 > index_2)
      {
	  aux[7] = Coord(index_2, index_1);
      }
      else
      {
	  aux[7] = Coord(index_1, index_2);
      }

      if (index_1 > index_3)
      {
	  aux[8] = Coord(index_3, index_1);
      }
      else
      {
	  aux[8] = Coord(index_1, index_3);
      }

      if (index_2 > index_3)
      {
	  aux[9] = Coord(index_3, index_2);
      }
      else
      {
	  aux[9] = Coord(index_2, index_3);
      }

      // Edge 01
      if (aux[4] < 0)
      {
	  if (index_0 > index_1)
	  {
	    edge_ids[0] = 0;
	  }
	  else
	  {
	    edge_ids[0] = 1;
	  }
      }
      else
      {
	  edge_ids[0] = 4;
      }

      // Edge 02
      if (aux[5] < 0)
	  if (index_0 > index_2)
	  {
	    edge_ids[1] = 0;
	  }
	  else
	  {
	    edge_ids[1] = 2;
	  }
      else
      {
	  edge_ids[1] = 5;
      }

      // Edge 03
      if (aux[6] < 0)
      {
	  if (index_0 > index_3)
	  {
	    edge_ids[2] = 0;
	  }
	  else
	  {
	    edge_ids[2] = 3;
	  }
      }
      else
      {
	  edge_ids[2] = 6;
      }

      // Edge 12
      if (aux[7] < 0)
      {
	  if (index_1 > index_2)
	  {
	    edge_ids[3] = 1;
	  }
	  else
	  {
	    edge_ids[3] = 2;
	  }
      }
      else
      {
	  edge_ids[3] = 7;
      }

      // Edge 13
      if (aux[8] < 0)
	  if (index_1 > index_3)
	  {
	    edge_ids[4] = 1;
	  }
	  else
	  {
	    edge_ids[4] = 3;
	  }
      else
      {
	  edge_ids[4] = 8;
      }

      // Edge 23
      if (aux[9] < 0)
      {
	  if (index_2 > index_3)
	  {
	    edge_ids[5] = 2;
	  }
	  else
	  {
	    edge_ids[5] = 3;
	  }
      }
      else
      {
	  edge_ids[5] = 9;
      }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It calculates the new edges of the new triangles,
    * first it calculates the new edges correspondign to the lower face (as a triangle),
    * later it added to the upper face
    * @param geom: The prism element geometry
    * @param edge_ids: The ids of the edges
    * @return aux: The vector that includes the index of the new edges
    */

    void  CalculateEdgesFaces(Element::GeometryType& geom,
		      const compressed_matrix<int>& Coord,
		      int*  edge_ids,
		      array_1d<int, 6>& aux
		      )
    {
        int index_0 = geom[0].Id() - 1;
        int index_1 = geom[1].Id() - 1;
        int index_2 = geom[2].Id() - 1;

        aux[0] = geom[0].Id();
        aux[1] = geom[1].Id();
        aux[2] = geom[2].Id();

        if (index_0 > index_1)
	{
            aux[3] = Coord(index_1, index_0);
	}
        else
	{
            aux[3] = Coord(index_0, index_1);
	}


        if (index_1 > index_2)
	{
            aux[4] = Coord(index_2, index_1);
	}
        else
	{
            aux[4] = Coord(index_1, index_2);
	}


        if (index_2 > index_0)
	{
            aux[5] = Coord(index_0, index_2);
	}
        else
	{
            aux[5] = Coord(index_2, index_0);
	}

        // Edge 01
        if (aux[3] < 0)
	{
            if (index_0 > index_1)
	    {
	      edge_ids[0] = 0;
	    }
            else
	    {
	      edge_ids[0] = 1;
	    }
	}
        else
	{
            edge_ids[0] = 3;
	}

        // Edge 12
        if (aux[4] < 0)
	{
            if (index_1 > index_2)
	    {
	      edge_ids[1] = 1;
	    }
            else
	    {
	      edge_ids[1] = 2;
	    }
	}
        else
	{
            edge_ids[1] = 4;
	}

        // Edge 20
        if (aux[5] < 0)
	{
            if (index_2 > index_0)
	    {
	      edge_ids[2] = 2;
	    }
            else
	    {
	      edge_ids[2] = 0;
	    }
	}
        else
	{
            edge_ids[2] = 5;
	}
    }

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

#endif // KRATOS_LOCAL_REFINE_TETRAHEDRA_MESH  defined
