// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pablo Agustin Becker
//

#pragma once

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

// External includes

// Project includes
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
class LocalRefineTriangleMeshGeneric 
    : public LocalRefineGeometryMesh
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

    /**
    * It erases the old elements and it creates the new ones
    * @param rModelPart: The model part of the model (it is the input too)
    * @param rCoord: The coordinates of the element
    * @param New_Elements: The new elements created
    * @param InterpolateInternalVariables: A boolean that defines if it is necessary to interpolate the internal variables
    */
    void EraseOldObjetcsAndCreateNewObjects(
        ModelPart& rThisModelPart,
        TArrayType& rObjects,
        const compressed_matrix<int>& rCoord,
        PointerVector< TGeometricalObjectType >& rNewObjects,
        bool InterpolateInternalVariables
        )
    {
        typename TArrayType::iterator GeometricalObjectsBegin = rObjects.ptr_begin();
        typename TArrayType::iterator GeometricalObjectsEnd  = rObjects.ptr_end();

        unsigned int to_be_deleted = 0;
        unsigned int large_id = (GeometricalObjectsEnd - 1)->Id() * 10;
        bool create_object = false;
        int EdgeIds[3];
        int t[12];
        int number_object = 0;
        int splitted_edges = 0;
        int nint = 0;
        std::vector<int> rAux;

        const ProcessInfo& rCurrentProcessInfo = rThisModelPart.GetProcessInfo();

	    KRATOS_INFO("LocalRefineTriangleMeshGeneric") << "****************** REFINING MESH ******************" << std::endl;
        KRATOS_INFO("LocalRefineTriangleMeshGeneric") << "OLD NUMBER OBJECTS: " << GeometricalObjectsEnd - GeometricalObjectsBegin  << std::endl;

        PointerVector< TGeometricalObjectType > Old_Objects;

        unsigned int current_id = (GeometricalObjectsEnd - 1)->Id() + 1;
        for (typename TArrayType::iterator it = GeometricalObjectsBegin; it != GeometricalObjectsEnd; ++it) {
            for (int & i : t) {
                i = -1;
            }
            typename TGeometricalObjectType::GeometryType& rGeom = it->GetGeometry();
            CalculateEdges(rGeom, rCoord, EdgeIds, rAux);

            const unsigned int dimension = rGeom.WorkingSpaceDimension();

            // It creates the new connectivities
            create_object = TriangleSplit::Split_Triangle(EdgeIds, t, &number_object, &splitted_edges, &nint);

            // It creates the new objects
            if (create_object) {
                to_be_deleted++;
                auto& r_child_objects = GetNeighbour(it); 
                r_child_objects.resize(0);
                for (int i = 0; i < number_object; i++) {
                    const unsigned int base = i * 3;
                    const unsigned int i0 = rAux[t[base]];
                    const unsigned int i1 = rAux[t[base + 1]];
                    const unsigned int i2 = rAux[t[base + 2]];

                    if (dimension == 2) {
                        Triangle2D3<Node > rGeom(
                            rThisModelPart.Nodes()(i0),
                            rThisModelPart.Nodes()(i1),
                            rThisModelPart.Nodes()(i2)
                        );

                        typename TGeometricalObjectType::Pointer p_object;
                        p_object = it->Create(current_id, rGeom, it->pGetProperties());
                        p_object->Initialize(rCurrentProcessInfo);
                        p_object->InitializeSolutionStep(rCurrentProcessInfo);
                        p_object->FinalizeSolutionStep(rCurrentProcessInfo);

                        // Setting the internal variables in the child elem
                        if (InterpolateInternalVariables) {
                            InterpolateInteralVariables(number_object, *it.base(), p_object, rCurrentProcessInfo);
			            }

                        // Transfer elemental variables
                        p_object->GetData() = it->GetData();
                        //const unsigned int& level = it->GetValue(REFINEMENT_LEVEL);
                        p_object->GetValue(SPLIT_ELEMENT) = false;
                        //p_element->SetValue(REFINEMENT_LEVEL, 1);
                        rNewObjects.push_back(p_object);
                        r_child_objects.push_back(typename TGeometricalObjectType::WeakPointer(p_object) );
                    } else {
                        Triangle3D3<Node > rGeom(
                            rThisModelPart.Nodes()(i0),
                            rThisModelPart.Nodes()(i1),
                            rThisModelPart.Nodes()(i2)
                        );

                        typename TGeometricalObjectType::Pointer p_object;
                        p_object = it->Create(current_id, rGeom, it->pGetProperties());
                        p_object->Initialize(rCurrentProcessInfo);
                        p_object->InitializeSolutionStep(rCurrentProcessInfo);
                        p_object->FinalizeSolutionStep(rCurrentProcessInfo);

                        // Setting the internal variables in the child elem
                        if (InterpolateInternalVariables == true)
                            InterpolateInteralVariables(number_object, *it.base(), p_object, rCurrentProcessInfo);

                        // Transfer elemental variables
                        p_object->GetData() = it->GetData();
                        p_object->GetValue(SPLIT_ELEMENT) = false;
                        rNewObjects.push_back(p_object);
                        r_child_objects.push_back(typename TGeometricalObjectType::WeakPointer(p_object) );
                    }

                    current_id++;
                }
                it->SetId(large_id);
                large_id++;
            }

        }

        /* Adding news elements to the model part */
        for (auto it_new = rNewObjects.begin(); it_new != rNewObjects.end(); it_new++) {
            rObjects.push_back(*(it_new.base()));
        }

        /* All of the elements to be erased are at the end */
        rObjects.Sort();

        /* Now remove all of the "old" elements */
        rObjects.erase(rObjects.end() - to_be_deleted, rObjects.end());

        KRATOS_INFO("LocalRefineTriangleMeshGeneric") << "NEW NUMBER OBJECTS: " << rObjects.size() << std::endl;
    }

    /**
    * It calculates the new edges of the new triangles,
    * first it calculates the new edges correspondingn to the lower face (as a triangle),
    * later it added to the upper face
    * @param rGeom: The triangle element geometry
    * @param EdgeIds: The ids of the edges
    * @return rAux: The vector that includes the index of the new edges
    */
    void CalculateEdges(
        Geometry<Node>& rGeom,
        const compressed_matrix<int>& rCoord,
        int* EdgeIds,
        std::vector<int>& rAux
        ) override
    {
        rAux.resize(6, false);

        const std::size_t index_0 = mMapNodeIdToPos[rGeom[0].Id()];
        const std::size_t index_1 = mMapNodeIdToPos[rGeom[1].Id()];
        const std::size_t index_2 = mMapNodeIdToPos[rGeom[2].Id()];

        rAux[0] = rGeom[0].Id();
        rAux[1] = rGeom[1].Id();
        rAux[2] = rGeom[2].Id();

        if (index_0 > index_1) {
            rAux[3] = rCoord(index_1, index_0);
	    } else {
            rAux[3] = rCoord(index_0, index_1);
	    }

        if (index_1 > index_2) {
            rAux[4] = rCoord(index_2, index_1);
	    } else {
            rAux[4] = rCoord(index_1, index_2);
	    }

        if (index_2 > index_0){
            rAux[5] = rCoord(index_0, index_2);
	    } else {
            rAux[5] = rCoord(index_2, index_0);
	    }

        // Edge 01
        if (rAux[3] < 0)	{
            if (index_0 > index_1){
	            EdgeIds[0] = 0;
	        } else {
	            EdgeIds[0] = 1;
	        }
	    } else {
            EdgeIds[0] = 3;
	    }

        // Edge 12
        if (rAux[4] < 0){
            if (index_1 > index_2) {
	            EdgeIds[1] = 1;
	        } else {
	            EdgeIds[1] = 2;
	        }
	    } else {
            EdgeIds[1] = 4;
	    }

        // Edge 20
        if (rAux[5] < 0) {
            if (index_2 > index_0) {
                EdgeIds[2] = 2;
            } else {
	            EdgeIds[2] = 0;
	        }
	    } else {
            EdgeIds[2] = 5;
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
