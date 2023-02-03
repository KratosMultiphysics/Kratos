#if !defined(KRATOS_LOCAL_REFINE_TRIANGLE_MESH_GENERIC)
#define  KRATOS_LOCAL_REFINE_TRIANGLE_MESH_GENERIC

#ifdef _OPENMP
#include <omp.h>
#endif

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

/* Project includes */
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "custom_utilities/local_refine_geometry_mesh.hpp"
#include "utilities/split_triangle.h"

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

template<class TGeometricalObjectType, typename TArrayType>
class LocalRefineTriangleMeshGeneric : public LocalRefineGeometryMesh
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LocalRefineTriangleMeshGeneric(ModelPart& model_part) : LocalRefineGeometryMesh(model_part)
    {

    }

    /// Destructor
    ~LocalRefineTriangleMeshGeneric() override
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It erases the old elements and it creates the new ones
    * @param rCoord: The coordinates of the element
    * @param New_Elements: The new elements created
    * @param interpolate_internal_variables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void EraseOldObjetcsAndCreateNewObjects(
            ModelPart& rThisModelPart,
            TArrayType& rObjects,
            const compressed_matrix<int>& rCoord,
            PointerVector< TGeometricalObjectType >& rNew_Objects,
            bool interpolate_internal_variables
    )
    {

        typename TArrayType::iterator GeometricalObjectsBegin = rObjects.ptr_begin();
        typename TArrayType::iterator GeometricalObjectsEnd  = rObjects.ptr_end();

        
        unsigned int to_be_deleted = 0;
        unsigned int large_id = (GeometricalObjectsEnd - 1)->Id() * 10;
        bool create_object = false;
        int edge_ids[3];
        int t[12];
        int number_object = 0;
        int splitted_edges = 0;
        int nint = 0;
        std::vector<int> aux;

        const ProcessInfo& rCurrentProcessInfo = rThisModelPart.GetProcessInfo();

	    std::cout << "****************** REFINING MESH ******************" << std::endl;
        std::cout << "OLD NUMBER OBJECTS: " << GeometricalObjectsEnd - GeometricalObjectsBegin  << std::endl;

        PointerVector< TGeometricalObjectType > Old_Objects;

        unsigned int current_id = (GeometricalObjectsEnd - 1)->Id() + 1;
        for (TArrayType::iterator it = GeometricalObjectsBegin; it != GeometricalObjectsEnd; ++it)
        {
            for (int & i : t)
            {
                i = -1;
            }
            TGeometricalObjectType::GeometryType& geom = it->GetGeometry();
            CalculateEdges(geom, rCoord, edge_ids, aux);

            const unsigned int dimension = geom.WorkingSpaceDimension();

            // It creates the new conectivities
            create_object = TriangleSplit::Split_Triangle(edge_ids, t, &number_object, &splitted_edges, &nint);

            // It creates the new objects
            if (create_object == true)
            {
                to_be_deleted++;
                GlobalPointersVector<TGeometricalObjectType> &rChildObjects = GetNeighbour(it); 
                rChildObjects.resize(0);
                for (int i = 0; i < number_object; i++)
                {

                    unsigned int base = i * 3;
                    unsigned int i0 = aux[t[base]];
                    unsigned int i1 = aux[t[base + 1]];
                    unsigned int i2 = aux[t[base + 2]];

                    if (dimension == 2)
                    {
                        Triangle2D3<Node < 3 > > geom(
                            rThisModelPart.Nodes()(i0),
                            rThisModelPart.Nodes()(i1),
                            rThisModelPart.Nodes()(i2)
                        );

                        TGeometricalObjectType::Pointer p_object;
                        p_object = it->Create(current_id, geom, it->pGetProperties());
                        p_object->Initialize(rCurrentProcessInfo);
                        p_object->InitializeSolutionStep(rCurrentProcessInfo);
                        p_object->FinalizeSolutionStep(rCurrentProcessInfo);

                        // Setting the internal variables in the child elem
                        if (interpolate_internal_variables == true)
			            {
                            InterpolateInteralVariables(number_object, *it.base(), p_object, rCurrentProcessInfo);
			            }

                        // Transfer elemental variables
                        p_object->GetData() = it->GetData();
                        //const unsigned int& level = it->GetValue(REFINEMENT_LEVEL);
                        p_object->GetValue(SPLIT_ELEMENT) = false;
                        //p_element->SetValue(REFINEMENT_LEVEL, 1);
                        rNew_Objects.push_back(p_object);
                        rChildObjects.push_back( TGeometricalObjectType::WeakPointer(p_object) );
                    }
                    else
                    {
                        Triangle3D3<Node < 3 > > geom(
                            rThisModelPart.Nodes()(i0),
                            rThisModelPart.Nodes()(i1),
                            rThisModelPart.Nodes()(i2)
                        );

                        TGeometricalObjectType::Pointer p_object;
                        p_object = it->Create(current_id, geom, it->pGetProperties());
                        p_object->Initialize(rCurrentProcessInfo);
                        p_object->InitializeSolutionStep(rCurrentProcessInfo);
                        p_object->FinalizeSolutionStep(rCurrentProcessInfo);

                        // Setting the internal variables in the child elem
                        if (interpolate_internal_variables == true)
                            InterpolateInteralVariables(number_object, *it.base(), p_object, rCurrentProcessInfo);

                        // Transfer elemental variables
                        p_object->GetData() = it->GetData();
                        p_object->GetValue(SPLIT_ELEMENT) = false;
                        rNew_Objects.push_back(p_object);
                        rChildObjects.push_back( TGeometricalObjectType::WeakPointer(p_object) );
                    }

                    current_id++;
                }
                it->SetId(large_id);
                large_id++;
            }

        }

        /* Adding news elements to the model part */
        for (PointerVector< TGeometricalObjectType >::iterator it_new = rNew_Objects.begin(); it_new != rNew_Objects.end(); it_new++)
        {
            rObjects.push_back(*(it_new.base()));
        }

        /* All of the elements to be erased are at the end */
        rObjects.Sort();

        /* Now remove all of the "old" elements */
        rObjects.erase(rObjects.end() - to_be_deleted, rObjects.end());

        std::cout << "NEW NUMBER OBJECTS: " << rObjects.size() << std::endl;
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
            Geometry<Node<3>>& geom,
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

        if (index_0 > index_1) {
            aux[3] = Coord(index_1, index_0);
	    }else{
            aux[3] = Coord(index_0, index_1);
	    }

        if (index_1 > index_2) {
            aux[4] = Coord(index_2, index_1);
	    }else{
            aux[4] = Coord(index_1, index_2);
	    }

        if (index_2 > index_0){
            aux[5] = Coord(index_0, index_2);
	    }else{
            aux[5] = Coord(index_2, index_0);
	    }

        // Edge 01
        if (aux[3] < 0)	{
            if (index_0 > index_1){
	            edge_ids[0] = 0;
	        }else{
	            edge_ids[0] = 1;
	        }
	    }else{
            edge_ids[0] = 3;
	    }

        // Edge 12
        if (aux[4] < 0){
            if (index_1 > index_2) {
	            edge_ids[1] = 1;
	        }else{
	            edge_ids[1] = 2;
	        }
	    }else{
            edge_ids[1] = 4;
	    }

        // Edge 20
        if (aux[5] < 0) {
            if (index_2 > index_0) {
                edge_ids[2] = 2;
            }else{
	            edge_ids[2] = 0;
	        }
	    }else{
            edge_ids[2] = 5;
	    }
    }

    GlobalPointersVector<Element>& GetNeighbour(ElementsArrayType::iterator it_elem)
    {
        return it_elem->GetValue(NEIGHBOUR_ELEMENTS);
    }

    GlobalPointersVector<Condition>& GetNeighbour(ConditionsArrayType::iterator it_cond)
    {
        return it_cond->GetValue(NEIGHBOUR_CONDITIONS);
    }
    
};

} // namespace Kratos.

#endif // KRATOS_LOCAL_REFINE_TRIANGLE_MESH  defined
