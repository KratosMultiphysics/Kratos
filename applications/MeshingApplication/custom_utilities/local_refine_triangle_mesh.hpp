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

#if !defined(KRATOS_LOCAL_REFINE_TRIANGLE_MESH)
#define  KRATOS_LOCAL_REFINE_TRIANGLE_MESH

#ifdef _OPENMP
#include <omp.h>
#endif

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

/* Project includes */
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "custom_utilities/local_refine_geometry_mesh.hpp"
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
class LocalRefineTriangleMesh : public LocalRefineGeometryMesh
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LocalRefineTriangleMesh(ModelPart& model_part) : LocalRefineGeometryMesh(model_part)
    {

    }

    /// Destructor
    ~LocalRefineTriangleMesh()
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * Computes the coordinate of the baricenter node of the element (mean of the faces's baricenter)
    * Insert the news nodes in the center of elements and interopolate the variables.
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void CalculateCoordinateCenterNodeAndInsertNewNodes(ModelPart& this_model_part)
    {
        array_1d<double, 3 > Coord_Node_1;
        array_1d<double, 3 > Coord_Node_2;
        array_1d<double, 3 > Coord_Node_3;
        array_1d<double, 3 > Coordinate_center_node;
        std::vector<int> node_center;

        NodesArrayType& pNodes = this_model_part.Nodes();
        int Id_Center = pNodes.size() + 1;
        ElementsArrayType& rElements = this_model_part.Elements();
        ElementsArrayType::iterator it_begin = rElements.ptr_begin();
        ElementsArrayType::iterator it_end = rElements.ptr_end();

        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
        {
            Element::GeometryType& geom = it->GetGeometry();
            noalias(Coord_Node_1) = geom[0].Coordinates();
            noalias(Coord_Node_2) = geom[1].Coordinates();
            noalias(Coord_Node_3) = geom[2].Coordinates();

            unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
            Node < 3 > ::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();
            noalias(Coordinate_center_node) = 0.33333333333333333 * (Coord_Node_1 + Coord_Node_2 + Coord_Node_3);

            // Inserting the news node in the model part
            Node < 3 > ::Pointer pnode = this_model_part.CreateNewNode(Id_Center, Coordinate_center_node[0], Coordinate_center_node[1], Coordinate_center_node[2]);
            pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());

            pnode->X0() = 0.33333333333333333 * (geom[0].X0() + geom[1].X0() + geom[2].X0());
            pnode->Y0() = 0.33333333333333333 * (geom[0].Y0() + geom[1].Y0() + geom[2].Y0());
            pnode->Z0() = 0.33333333333333333 * (geom[0].Z0() + geom[1].Z0() + geom[2].Z0());

            for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
            {
                Node < 3 > ::DofType& rDof = *iii;
                Node < 3 > ::DofType::Pointer p_new_dof = pnode->pAddDof(rDof);
                if (geom[0].IsFixed(iii->GetVariable()) == true && geom[1].IsFixed(iii->GetVariable()) == true && geom[2].IsFixed(iii->GetVariable()))
                    (p_new_dof)->FixDof();
                else
                {
                    (p_new_dof)->FreeDof();
                }

            }

            // Interpolating the data
            unsigned int buffer_size = pnode->GetBufferSize();
            for (unsigned int step = 0; step < buffer_size; step++)
            {
                double* new_step_data = pnode->SolutionStepData().Data(step);
                double* step_data1 = geom[0].SolutionStepData().Data(step);
                double* step_data2 = geom[1].SolutionStepData().Data(step);
                double* step_data3 = geom[2].SolutionStepData().Data(step);

                // Copying this data in the position of the vector we are interested in
                for (unsigned int j = 0; j < step_data_size; j++)
                {
                    new_step_data[j] = 0.333333333333333333 * (step_data1[j] + step_data2[j] + step_data3[j]);
                }
            }

            node_center.push_back(Id_Center);
            Id_Center++;
        }
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
            PointerVector< Element >& New_Elements,
            bool interpolate_internal_variables
    ) override
    {
        ElementsArrayType& rElements = this_model_part.Elements();
        ElementsArrayType::iterator it_begin = rElements.ptr_begin();
        ElementsArrayType::iterator it_end = rElements.ptr_end();

        unsigned int to_be_deleted = 0;
        unsigned int large_id = (rElements.end() - 1)->Id() * 10;
        bool create_element = false;
        int edge_ids[3];
        int t[12];
        int number_elem = 0;
        int splitted_edges = 0;
        int nint = 0;
        std::vector<int> aux;

        ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();

	std::cout << "****************** REFINING MESH ******************" << std::endl;
        std::cout << "OLD NUMBER ELEMENTS: " << rElements.size() << std::endl;

        PointerVector< Element > Old_Elements;

        unsigned int current_id = (rElements.end() - 1)->Id() + 1;
        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
        {
            for (int & i : t)
            {
                i = -1;
            }
            Element::GeometryType& geom = it->GetGeometry();
            CalculateEdges(geom, Coord, edge_ids, aux);

            const unsigned int dimension = geom.WorkingSpaceDimension();

            // It creates the new conectivities
            create_element = Split_Triangle(edge_ids, t, &number_elem, &splitted_edges, &nint);

            // It creates the new elements
            if (create_element == true)
            {
                to_be_deleted++;
                for (int i = 0; i < number_elem; i++)
                {

                    unsigned int base = i * 3;
                    unsigned int i0 = aux[t[base]];
                    unsigned int i1 = aux[t[base + 1]];
                    unsigned int i2 = aux[t[base + 2]];

                    if (dimension == 2)
                    {
                        Triangle2D3<Node < 3 > > geom(
                            this_model_part.Nodes()(i0),
                            this_model_part.Nodes()(i1),
                            this_model_part.Nodes()(i2)
                        );

                        Element::Pointer p_element;
                        p_element = it->Create(current_id, geom, it->pGetProperties());
                        p_element->Initialize();
                        p_element->InitializeSolutionStep(rCurrentProcessInfo);
                        p_element->FinalizeSolutionStep(rCurrentProcessInfo);

                        // Setting the internal variables in the child elem
                        if (interpolate_internal_variables == true)
			{
                            InterpolateInteralVariables(number_elem, *it.base(), p_element, rCurrentProcessInfo);
			}

                        // Transfer elemental variables
                        p_element->Data() = it->Data();
                        //const unsigned int& level = it->GetValue(REFINEMENT_LEVEL);
                        p_element->GetValue(SPLIT_ELEMENT) = false;
                        //p_element->SetValue(REFINEMENT_LEVEL, 1);
                        New_Elements.push_back(p_element);
                    }
                    else
                    {
                        Triangle3D3<Node < 3 > > geom(
                            this_model_part.Nodes()(i0),
                            this_model_part.Nodes()(i1),
                            this_model_part.Nodes()(i2)
                        );

                        Element::Pointer p_element;
                        p_element = it->Create(current_id, geom, it->pGetProperties());
                        p_element->Initialize();
                        p_element->InitializeSolutionStep(rCurrentProcessInfo);
                        p_element->FinalizeSolutionStep(rCurrentProcessInfo);

                        // Setting the internal variables in the child elem
                        if (interpolate_internal_variables == true)
                            InterpolateInteralVariables(number_elem, *it.base(), p_element, rCurrentProcessInfo);

                        // Transfer elemental variables
                        p_element->Data() = it->Data();
                        p_element->GetValue(SPLIT_ELEMENT) = false;
                        New_Elements.push_back(p_element);
                    }

                    current_id++;
                }
                it->SetId(large_id);
                large_id++;
            }

        }

        /* Adding news elements to the model part */
        for (PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new != New_Elements.end(); it_new++)
        {
            rElements.push_back(*(it_new.base()));
        }

        /* All of the elements to be erased are at the end */
        rElements.Sort();

        /* Now remove all of the "old" elements */
        rElements.erase(this_model_part.Elements().end() - to_be_deleted, this_model_part.Elements().end());

        std::cout << "NEW NUMBER ELEMENTS: " << rElements.size() << std::endl;
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

        PointerVector< Condition > New_Conditions;

        ConditionsArrayType& rConditions = this_model_part.Conditions();

        if (rConditions.size() > 0)
        {
            ConditionsArrayType::iterator it_begin = rConditions.ptr_begin();
            ConditionsArrayType::iterator it_end   = rConditions.ptr_end();
            unsigned int to_be_deleted = 0;
            unsigned int large_id = (rConditions.end() - 1)->Id() * 7;

            ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();

	    const unsigned int dimension = it_begin->GetGeometry().WorkingSpaceDimension();

            unsigned int current_id = (rConditions.end() - 1)->Id() + 1;

	    for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it)
	    {
                Condition::GeometryType& geom = it->GetGeometry();

                if (geom.size() == 2)
                {
                    int index_0 = geom[0].Id() - 1;
                    int index_1 = geom[1].Id() - 1;
                    int new_id;

                    if (index_0 > index_1)
		    {
                        new_id = Coord(index_1, index_0);
		    }
                    else
		    {
                        new_id = Coord(index_0, index_1);
		    }

                    if (new_id > 0) // We need to create a new condition
                    {
                        to_be_deleted++;
                        if (dimension == 2)
                        {
                            Line2D2<Node < 3 > > newgeom1(
                                this_model_part.Nodes()(geom[0].Id()),
                                this_model_part.Nodes()(new_id)
                            );

                            Line2D2<Node < 3 > > newgeom2(
                                this_model_part.Nodes()(new_id),
                                this_model_part.Nodes()(geom[1].Id())
                            );

                            Condition::Pointer pcond1 = it->Create(current_id++, newgeom1, it->pGetProperties());
                            Condition::Pointer pcond2 = it->Create(current_id++, newgeom2, it->pGetProperties());

                            pcond1->Data() = it->Data();
                            pcond2->Data() = it->Data();

                            New_Conditions.push_back(pcond1);
                            New_Conditions.push_back(pcond2);
                        }
                        else
                        {
                            Line3D2<Node < 3 > > newgeom1(
                                this_model_part.Nodes()(geom[0].Id()),
                                this_model_part.Nodes()(new_id)
                            );

                            Line3D2<Node < 3 > > newgeom2(
                                this_model_part.Nodes()(new_id),
                                this_model_part.Nodes()(geom[1].Id())
                            );

                            Condition::Pointer pcond1 = it->Create(current_id++, newgeom1, it->pGetProperties());
                            Condition::Pointer pcond2 = it->Create(current_id++, newgeom2, it->pGetProperties());

                            pcond1->Data() = it->Data();
                            pcond2->Data() = it->Data();

                            New_Conditions.push_back(pcond1);
                            New_Conditions.push_back(pcond2);
                        }

                        it->SetId(large_id);
                        large_id++;
                    }
                }
            }

            /* All of the elements to be erased are at the end */
            this_model_part.Conditions().Sort();

            /* Remove all of the "old" elements */
            this_model_part.Conditions().erase(this_model_part.Conditions().end() - to_be_deleted, this_model_part.Conditions().end());

            unsigned int total_size = this_model_part.Conditions().size() + New_Conditions.size();
            this_model_part.Conditions().reserve(total_size);

            /* Adding news elements to the model part */
            for (PointerVector< Condition >::iterator it_new = New_Conditions.begin(); it_new != New_Conditions.end(); it_new++)
            {
                it_new->Initialize();
                it_new->InitializeSolutionStep(rCurrentProcessInfo);
                it_new->FinalizeSolutionStep(rCurrentProcessInfo);
                this_model_part.Conditions().push_back(*(it_new.base()));
            }

            /* Renumber */
            unsigned int my_index = 1;
            for (ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); it++)
            {
                it->SetId(my_index++);
            }

        }

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It calculates the new edges of the new triangles,
    * first it calculates the new edges correspondign to the lower face (as a triangle),
    * later it added to the upper face
    * @param geom: The triangle element geometry
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
        aux.resize(6, false);

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

#endif // KRATOS_LOCAL_REFINE_TRIANGLE_MESH  defined
