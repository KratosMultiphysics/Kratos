// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <unordered_set>

// External includes
/* The includes related with the PMMG library */
#if !defined(PMMG_INCLUDES)
#define PMMG_INCLUDES
#include "parmmg/libparmmg.h"
#endif /* PMMG_INCLUDES defined */

// Project includes
#include "containers/model.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "custom_utilities/pmmg_utilities.h"
#include "meshing_application_variables.h"

// NOTE: The following contains the license of the PMMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

namespace Kratos
{

template<>
const SizeType PMMGMeshInfo<PMMGLibrary::PMMG3D>::NumberFirstTypeConditions() const
{
    return NumberOfTriangles;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
const SizeType PMMGMeshInfo<PMMGLibrary::PMMG3D>::NumberSecondTypeConditions() const
{
    return NumberOfQuadrilaterals;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
const SizeType PMMGMeshInfo<PMMGLibrary::PMMG3D>::NumberFirstTypeElements() const
{
    return NumberOfTetrahedra;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
const SizeType PMMGMeshInfo<PMMGLibrary::PMMG3D>::NumberSecondTypeElements() const
{
    return NumberOfPrism;
}

template struct PMMGMeshInfo<PMMGLibrary::PMMG3D>;

/***********************************************************************************/
/***********************************************************************************/

// The member variables related with the PMMG library
PMMG_pParMesh mParMmgMesh;      /// The mesh data from PMMG

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::SetEchoLevel(const SizeType EchoLevel)
{
    mEchoLevel = EchoLevel;
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
SizeType ParMmgUtilities<TPMMGLibrary>::GetEchoLevel()
{
    return mEchoLevel;
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::SetDiscretization(const DiscretizationOption Discretization)
{
    mDiscretization = Discretization;
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
DiscretizationOption ParMmgUtilities<TPMMGLibrary>::GetDiscretization()
{
    return mDiscretization;
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::SetRemoveRegions(const bool RemoveRegions)
{
    mRemoveRegions = RemoveRegions;
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
bool ParMmgUtilities<TPMMGLibrary>::GetRemoveRegions()
{
    return mRemoveRegions;
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::PrintAndGetParMmgMeshInfo(PMMGMeshInfo<TPMMGLibrary>& rPMMGMeshInfo)
{
    int np,ne,nprism,nt,nquad,na;
    KRATOS_ERROR_IF(PMMG_Get_meshSize(mParMmgMesh,&np,&ne,&nprism,&nt,&nquad,&na) != 1) << "Unable to get mesh size" << std::endl;

    /* Warning: mesh groups must be merged on each local partition */
    rPMMGMeshInfo.NumberOfNodes = np;
    if (TPMMGLibrary == PMMGLibrary::PMMG3D) { // 3D
        rPMMGMeshInfo.NumberOfTriangles = nt;
        rPMMGMeshInfo.NumberOfQuadrilaterals = nquad;
        rPMMGMeshInfo.NumberOfTetrahedra = ne;
        rPMMGMeshInfo.NumberOfPrism = nprism;
    }

    KRATOS_INFO_IF("ParMmgUtilities", mEchoLevel > 0) << "\tNodes created: " << rPMMGMeshInfo.NumberOfNodes << std::endl;
    if (TPMMGLibrary == PMMGLibrary::PMMG3D) { // 3D
        KRATOS_INFO_IF("ParMmgUtilities", mEchoLevel > 0) <<
        "Conditions created: " << rPMMGMeshInfo.NumberOfTriangles + rPMMGMeshInfo.NumberOfQuadrilaterals << "\n\tTriangles: " << rPMMGMeshInfo.NumberOfTriangles << "\tQuadrilaterals: " << rPMMGMeshInfo.NumberOfQuadrilaterals << "\n" <<
        "Elements created: " << rPMMGMeshInfo.NumberOfTetrahedra + rPMMGMeshInfo.NumberOfPrism << "\n\tTetrahedron: " << rPMMGMeshInfo.NumberOfTetrahedra << "\tPrisms: " << rPMMGMeshInfo.NumberOfPrism << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
IndexVectorType ParMmgUtilities<TPMMGLibrary>::FindDuplicateNodeIds(const ModelPart& rModelPart)
{
    DoubleVectorMapType node_map;

    IndexVectorType nodes_to_remove_ids;

    DoubleVectorType coords(Dimension);

    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    for(SizeType i = 0; i < r_nodes_array.size(); ++i) {
        auto it_node = it_node_begin + i;

        const array_1d<double, 3>& r_coordinates = it_node->Coordinates();

        for(IndexType i_coord = 0; i_coord < Dimension; i_coord++)
            coords[i_coord] = r_coordinates[i_coord];

        auto& r_count = node_map[coords];
        r_count += 1;

        if (r_count > 1) {
            nodes_to_remove_ids.push_back(it_node->Id());
            KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 0) << "The mode " << it_node->Id() <<  " is repeated"<< std::endl;
        }
    }

    return nodes_to_remove_ids;
}
/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType ParMmgUtilities<PMMGLibrary::PMMG3D>::CheckFirstTypeConditions()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_triangles(3);

    IndexVectorType conditions_to_remove;

    int np,ne,nprism,nt,nquad,na;
    KRATOS_ERROR_IF(PMMG_Get_meshSize(mParMmgMesh,&np,&ne,&nprism,&nt,&nquad,&na) !=1 ) << "Unable to get mesh size" << std::endl;

    for(int i = 0; i < nt; ++i) {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;

        KRATOS_ERROR_IF(PMMG_Get_triangle(mParMmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 ) << "Unable to get triangle" << std::endl;

        ids_triangles[0] = mLocalToGlobal[vertex_0];
        ids_triangles[1] = mLocalToGlobal[vertex_1];
        ids_triangles[2] = mLocalToGlobal[vertex_2];

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_triangles.begin(), ids_triangles.end());

        auto& r_count = triangle_map[ids_triangles];
        r_count += 1;

        if (r_count > 1)
            conditions_to_remove.push_back(i + 1);
    }

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType ParMmgUtilities<PMMGLibrary::PMMG3D>::CheckSecondTypeConditions()
{
    IndexVectorMapType quadrilateral_map;

    IndexVectorType ids_quadrialteral(4);

    IndexVectorType conditions_to_remove;

    int np,ne,nprism,nt,nquad,na;
    KRATOS_ERROR_IF(PMMG_Get_meshSize(mParMmgMesh,&np,&ne,&nprism,&nt,&nquad,&na) !=1 ) << "Unable to get mesh size" << std::endl;

    for(int i = 0; i < nquad; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, prop_id, is_required;

        KRATOS_ERROR_IF(PMMG_Get_quadrilateral(mParMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &prop_id, &is_required) != 1 ) << "Unable to get quadrilateral" << std::endl;

        ids_quadrialteral[0] = mLocalToGlobal[vertex_0];
        ids_quadrialteral[1] = mLocalToGlobal[vertex_1];
        ids_quadrialteral[2] = mLocalToGlobal[vertex_2];
        ids_quadrialteral[3] = mLocalToGlobal[vertex_3];

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_quadrialteral.begin(), ids_quadrialteral.end());

        auto& r_count = quadrilateral_map[ids_quadrialteral];
        r_count += 1;

        if (r_count > 1)
            conditions_to_remove.push_back(i + 1);
    }

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType ParMmgUtilities<PMMGLibrary::PMMG3D>::CheckFirstTypeElements()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_tetrahedron(4);

    IndexVectorType elements_to_remove;

    int np,ne,nprism,nt,nquad,na;
    KRATOS_ERROR_IF(PMMG_Get_meshSize(mParMmgMesh,&np,&ne,&nprism,&nt,&nquad,&na) !=1 ) << "Unable to get mesh size" << std::endl;

    for(int i = 0; i < ne; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, prop_id, is_required;

        KRATOS_ERROR_IF(PMMG_Get_tetrahedron(mParMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &prop_id, &is_required) != 1 ) << "Unable to get tetrahedron" << std::endl;

        ids_tetrahedron[0] = mLocalToGlobal[vertex_0];
        ids_tetrahedron[1] = mLocalToGlobal[vertex_1];
        ids_tetrahedron[2] = mLocalToGlobal[vertex_2];
        ids_tetrahedron[3] = mLocalToGlobal[vertex_3];

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_tetrahedron.begin(), ids_tetrahedron.end());

        auto& r_count = triangle_map[ids_tetrahedron];
        r_count += 1;

        if (r_count > 1)
            elements_to_remove.push_back(i + 1);
    }

    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType ParMmgUtilities<PMMGLibrary::PMMG3D>::CheckSecondTypeElements()
{
    IndexVectorMapType prism_map;

    IndexVectorType ids_prisms(6);

    IndexVectorType elements_to_remove;

    int np,ne,nprism,nt,nquad,na;
    KRATOS_ERROR_IF(PMMG_Get_meshSize(mParMmgMesh,&np,&ne,&nprism,&nt,&nquad,&na) !=1 ) << "Unable to get mesh size" << std::endl;

    for(int i = 0; i < nprism; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, prop_id, is_required;

        KRATOS_ERROR_IF(PMMG_Get_prism(mParMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &vertex_4, &vertex_5, &prop_id, &is_required) != 1 ) << "Unable to get prism" << std::endl;

        ids_prisms[0] = mLocalToGlobal[vertex_0];
        ids_prisms[1] = mLocalToGlobal[vertex_1];
        ids_prisms[2] = mLocalToGlobal[vertex_2];
        ids_prisms[3] = mLocalToGlobal[vertex_3];
        ids_prisms[4] = mLocalToGlobal[vertex_4];
        ids_prisms[5] = mLocalToGlobal[vertex_5];

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_prisms.begin(), ids_prisms.end());

        auto& r_count = prism_map[ids_prisms];
        r_count += 1;

        if (r_count > 1)
            elements_to_remove.push_back(i + 1);
    }

    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/


template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::BlockNode(const IndexType iNode)
{
    KRATOS_ERROR_IF(PMMG_Set_requiredVertex(mParMmgMesh, iNode) != 1 ) << "Unable to block vertex" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::BlockCondition(const IndexType iCondition)
{
    KRATOS_ERROR_IF(PMMG_Set_requiredTriangle(mParMmgMesh, iCondition) != 1 ) << "Unable to block edge" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/


template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::BlockElement(const IndexType iElement)
{
    KRATOS_ERROR_IF(PMMG_Set_requiredTetrahedron(mParMmgMesh, iElement) != 1 ) << "Unable to block tetrahedron" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Node<3>::Pointer ParMmgUtilities<PMMGLibrary::PMMG3D>::CreateNode(
    ModelPart& rModelPart,
    const IndexType iNode,
    int& Ref,
    int& IsRequired
    )
{
    double coord_0, coord_1, coord_2;
    int is_corner;

    KRATOS_ERROR_IF(PMMG_Get_vertex(mParMmgMesh, &coord_0, &coord_1, &coord_2, &Ref, &is_corner, &IsRequired) != 1 ) << "Unable to get vertex" << std::endl;

    // std::cout<< iNode << " " << coord_0 <<  " " << coord_1 << " " << coord_2 << " " << Ref << std::endl;

    NodeType::Pointer p_node = rModelPart.CreateNewNode(iNode, coord_0, coord_1, coord_2);

    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Condition::Pointer ParMmgUtilities<PMMGLibrary::PMMG3D>::CreateFirstTypeCondition(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
    const IndexType CondId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    // We create the default one
    Condition::Pointer p_condition = nullptr;

    int vertex_0, vertex_1, vertex_2;

    KRATOS_ERROR_IF(PMMG_Get_triangle(mParMmgMesh, &vertex_0, &vertex_1, &vertex_2, &Ref, &IsRequired) != 1 ) << "Unable to get triangle" << std::endl;

    // Sometimes PMMG creates conditions where there are not, then we skip
    Properties::Pointer p_prop = nullptr;
    Condition::Pointer p_base_condition = nullptr;

    if (rMapPointersRefCondition[Ref].get() == nullptr) {
        if (mDiscretization != DiscretizationOption::ISOSURFACE) { // The ISOSURFACE method creates new conditions from scratch, so we allow no previous Properties
            // KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "Condition. Null pointer returned" << std::endl;
            return p_condition;
        } else {
            p_prop = rModelPart.pGetProperties(0);
            PointerVector<NodeType> dummy_nodes (3);
            p_base_condition = KratosComponents<Condition>::Get("SurfaceCondition3D3N").Create(0, dummy_nodes, p_prop);
            p_base_condition->Set(MARKER);
        }
    } else {
        p_base_condition = rMapPointersRefCondition[Ref];
        p_prop = p_base_condition->pGetProperties();
    }

    // FIXME: This is not the correct solution to the problem, I asked in the PMMG Forum
    if (mLocalToGlobal[vertex_0] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_1] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_2] == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> condition_nodes (3);
        condition_nodes[0] = rModelPart.pGetNode(mLocalToGlobal[vertex_0]);
        condition_nodes[1] = rModelPart.pGetNode(mLocalToGlobal[vertex_1]);
        condition_nodes[2] = rModelPart.pGetNode(mLocalToGlobal[vertex_2]);

        p_condition = p_base_condition->Create(CondId, PointerVector<NodeType>{condition_nodes}, p_prop);
        if (p_base_condition->Is(MARKER)) p_condition->Set(MARKER);
    } else if (mEchoLevel > 2)
        KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "Condition creation avoided" << std::endl;

    if (p_condition != nullptr) KRATOS_ERROR_IF(p_condition->GetGeometry().Area() < ZeroTolerance) << "Creating a almost zero or negative area condition" << std::endl;
    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Condition::Pointer ParMmgUtilities<PMMGLibrary::PMMG3D>::CreateSecondTypeCondition(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
    const IndexType CondId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    Condition::Pointer p_condition = nullptr;

    int vertex_0, vertex_1, vertex_2, vertex_3;

    KRATOS_ERROR_IF(PMMG_Get_quadrilateral(mParMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &Ref, &IsRequired) != 1 ) << "Unable to get quadrilateral" << std::endl;

    // Sometimes PMMG creates conditions where there are not, then we skip
    if (rMapPointersRefCondition[Ref].get() == nullptr) {
        KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "Condition. Null pointer returned" << std::endl;
        return p_condition;
    }

    // FIXME: This is not the correct solution to the problem, I asked in the PMMG Forum
    if (mLocalToGlobal[vertex_0] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_1] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_2] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_3] == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> condition_nodes (4);
        condition_nodes[0] = rModelPart.pGetNode(mLocalToGlobal[vertex_0]);
        condition_nodes[1] = rModelPart.pGetNode(mLocalToGlobal[vertex_1]);
        condition_nodes[2] = rModelPart.pGetNode(mLocalToGlobal[vertex_2]);
        condition_nodes[3] = rModelPart.pGetNode(mLocalToGlobal[vertex_3]);

        p_condition = rMapPointersRefCondition[Ref]->Create(CondId, PointerVector<NodeType>{condition_nodes}, rMapPointersRefCondition[Ref]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "Condition creation avoided" << std::endl;

    if (p_condition != nullptr) KRATOS_ERROR_IF(p_condition->GetGeometry().Area() < ZeroTolerance) << "Creating a almost zero or negative area condition" << std::endl;
    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Element::Pointer ParMmgUtilities<PMMGLibrary::PMMG3D>::CreateFirstTypeElement(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
    const IndexType ElemId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    Element::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2, vertex_3;

    KRATOS_ERROR_IF(PMMG_Get_tetrahedron(mParMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &Ref, &IsRequired) != 1 ) << "Unable to get tetrahedron" << std::endl;

    // std::cout << "elem id id id id " << ElemId << " "<< vertex_0 << " "<< vertex_1 << " "<< vertex_2 << " "<< vertex_3 << std::endl;
    // std::cout << "elem id id id id " << ElemId << " "<< mLocalToGlobal[vertex_0] << " "<< mLocalToGlobal[vertex_1] << " "<< mLocalToGlobal[vertex_2] << " "<< mLocalToGlobal[vertex_3] << std::endl;

    if (mDiscretization == DiscretizationOption::ISOSURFACE) {
        // The existence of a _nullptr_ indicates an element that was removed. This is not an alarming indicator.
        if (rMapPointersRefElement[Ref].get() == nullptr) {
            // KRATOS_INFO("ParMmgUtilities") << "Element has been removed from domain. Ok." << std::endl;
            return p_element;
        } else {
            if (mLocalToGlobal[vertex_0] == 0) SkipCreation = true;
            if (mLocalToGlobal[vertex_1] == 0) SkipCreation = true;
            if (mLocalToGlobal[vertex_2] == 0) SkipCreation = true;
            if (mLocalToGlobal[vertex_3] == 0) SkipCreation = true;
            if (!SkipCreation) {
                std::vector<NodeType::Pointer> element_nodes (4);
                element_nodes[0] = rModelPart.pGetNode(mLocalToGlobal[vertex_0]);
                element_nodes[1] = rModelPart.pGetNode(mLocalToGlobal[vertex_1]);
                element_nodes[2] = rModelPart.pGetNode(mLocalToGlobal[vertex_2]);
                element_nodes[3] = rModelPart.pGetNode(mLocalToGlobal[vertex_3]);
                p_element = rMapPointersRefElement[Ref]->Create(ElemId, PointerVector<NodeType>{element_nodes}, rMapPointersRefElement[Ref]->pGetProperties());

                // Setting inside flag
                if (Ref == 2) {
                    p_element->Set(INSIDE, true);
                } else if (Ref == 3) {
                    p_element->Set(INSIDE, false);
                    if (mRemoveRegions) p_element->Set(TO_ERASE, true);
                }
            }
        }
    } else {
        // The existence of a _nullptr_ indicates a missing element. Two options are possible: error or replacement
        Properties::Pointer p_prop = nullptr;
        Element::Pointer p_base_element = nullptr;

        // Sometimes PMMG creates elements where there are not, then we skip
        if (rMapPointersRefElement[Ref].get() == nullptr) {
            KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "Element. Null pointer returned" << std::endl;
            return p_element;
        } else {
            p_base_element = rMapPointersRefElement[Ref];
            p_prop = p_base_element->pGetProperties();
        }

        // FIXME: This is not the correct solution to the problem, I asked in the PMMG Forum
        if (mLocalToGlobal[vertex_0] == 0) SkipCreation = true;
        if (mLocalToGlobal[vertex_1] == 0) SkipCreation = true;
        if (mLocalToGlobal[vertex_2] == 0) SkipCreation = true;
        if (mLocalToGlobal[vertex_3] == 0) SkipCreation = true;

        if (!SkipCreation) {
            std::vector<NodeType::Pointer> element_nodes (4);
            element_nodes[0] = rModelPart.pGetNode(mLocalToGlobal[vertex_0]);
            element_nodes[1] = rModelPart.pGetNode(mLocalToGlobal[vertex_1]);
            element_nodes[2] = rModelPart.pGetNode(mLocalToGlobal[vertex_2]);
            element_nodes[3] = rModelPart.pGetNode(mLocalToGlobal[vertex_3]);

            p_element = p_base_element->Create(ElemId, PointerVector<NodeType>{element_nodes}, p_prop);
        } else if (mEchoLevel > 2)
            KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "Element creation avoided" << std::endl;
    }

    if (p_element!= nullptr) KRATOS_ERROR_IF(p_element->GetGeometry().Volume() < ZeroTolerance) << "Creating a almost zero or negative volume element" << std::endl;
    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Element::Pointer ParMmgUtilities<PMMGLibrary::PMMG3D>::CreateSecondTypeElement(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
    const IndexType ElemId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    Element::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5;

    KRATOS_ERROR_IF(PMMG_Get_prism(mParMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &vertex_4, &vertex_5, &Ref, &IsRequired) != 1 ) << "Unable to get prism" << std::endl;

    // Sometimes PMMG creates elements where there are not, then we skip
    if (rMapPointersRefElement[Ref].get() == nullptr) {
        KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "Element. Null pointer returned" << std::endl;
        return p_element;
    }

    // FIXME: This is not the correct solution to the problem, I asked in the PMMG Forum
    if (mLocalToGlobal[vertex_0] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_1] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_2] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_3] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_4] == 0) SkipCreation = true;
    if (mLocalToGlobal[vertex_5] == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> element_nodes (6);
        element_nodes[0] = rModelPart.pGetNode(mLocalToGlobal[vertex_0]);
        element_nodes[1] = rModelPart.pGetNode(mLocalToGlobal[vertex_1]);
        element_nodes[2] = rModelPart.pGetNode(mLocalToGlobal[vertex_2]);
        element_nodes[3] = rModelPart.pGetNode(mLocalToGlobal[vertex_3]);
        element_nodes[4] = rModelPart.pGetNode(mLocalToGlobal[vertex_4]);
        element_nodes[5] = rModelPart.pGetNode(mLocalToGlobal[vertex_5]);

        p_element = rMapPointersRefElement[Ref]->Create(ElemId, PointerVector<NodeType>{element_nodes}, rMapPointersRefElement[Ref]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "Element creation avoided" << std::endl;

    if (p_element!= nullptr) KRATOS_ERROR_IF(p_element->GetGeometry().Volume() < ZeroTolerance) << "Creating a almost zero or negative volume element" << std::endl;
    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::InitMesh()
{
    mParMmgMesh = nullptr;

    // We init the PMMG mesh and sol
    if (mDiscretization == DiscretizationOption::STANDARD) {
        PMMG_Init_parMesh( PMMG_ARG_start, PMMG_ARG_ppParMesh, &mParMmgMesh, PMMG_ARG_pMesh, PMMG_ARG_pMet, PMMG_ARG_dim, 3, PMMG_ARG_MPIComm, MPI_COMM_WORLD, PMMG_ARG_end);
    } else {
        KRATOS_ERROR << "Discretization type: " << static_cast<int>(mDiscretization) << " not fully implemented" << std::endl;
    }

    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::InitVerbosity()
{
    /* We set the PMMG verbosity */
    int verbosity_mmg;
    if (mEchoLevel == 0)
        verbosity_mmg = -1;
    else if (mEchoLevel == 1)
        verbosity_mmg = 0; // NOTE: This way just the essential info from PMMG will be printed, but the custom message will appear
    else if (mEchoLevel == 2)
        verbosity_mmg = 3;
    else if (mEchoLevel == 3)
        verbosity_mmg = 5;
    else
        verbosity_mmg = 10;

    InitVerbosityParameter(verbosity_mmg);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::InitAPIModeParameter(const IndexType APImode)
{
    KRATOS_ERROR_IF( !PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_APImode, APImode) ) << "Unable to set API mode" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::InitNodeGloNumParameter(const IndexType nodeGloNum)
{
    KRATOS_ERROR_IF( !PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_nodeGloNum, nodeGloNum) ) << "Unable to set node global numbering" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::InitVerbosityParameter(const IndexType VerbosityPMMG)
{
    KRATOS_ERROR_IF( !PMMG_Set_iparameter(mParMmgMesh, PMMG_IPARAM_verbose, VerbosityPMMG) ) << "Unable to set verbosity" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetMeshSize(PMMGMeshInfo<PMMGLibrary::PMMG3D>& rPMMGMeshInfo)
{
    //Give the size of the mesh: NumNodes Vertex, num_elements tetra and prism, NumArrayConditions triangles and quadrilaterals, 0 edges (3D)
    KRATOS_ERROR_IF( PMMG_Set_meshSize(mParMmgMesh, rPMMGMeshInfo.NumberOfNodes, rPMMGMeshInfo.NumberOfTetrahedra, rPMMGMeshInfo.NumberOfPrism, rPMMGMeshInfo.NumberOfTriangles, rPMMGMeshInfo.NumberOfQuadrilaterals, 0) != 1 ) << "Unable to set mesh size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetSolSizeScalar(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( PMMG_Set_metSize(mParMmgMesh, MMG5_Vertex, NumNodes, MMG5_Scalar) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetSolSizeVector(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( PMMG_Set_metSize(mParMmgMesh, MMG5_Vertex, NumNodes, MMG5_Vector) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetSolSizeTensor(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( PMMG_Set_metSize(mParMmgMesh,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetDispSizeVector(const SizeType NumNodes)
{
    // TODO: Reactivate when dependency problem is solved
//     KRATOS_ERROR_IF( PMMG_Set_iparameter(mParMmgMesh,mParMmgDisp,PMMG_IPARAM_lag, 1) != 1 ) << "Unable to set lagrangian movement" << std::endl;
    // KRATOS_ERROR_IF( PMMG_Set_solSize(mParMmgMesh,mParMmgDisp,MMG5_Vertex,NumNodes,MMG5_Vector) != 1 ) << "Unable to set displacement size" << std::endl;
    KRATOS_ERROR << "SetDispSizeVector Not Yet Implemented" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::CheckMeshData()
{
    KRATOS_ERROR << "CheckMeshData Not Yet Implemented" << std::endl;
//     if (mDiscretization == DiscretizationOption::LAGRANGIAN) {
//         KRATOS_ERROR_IF( PMMG_Check(mParMmgMesh, mParMmgSol) != 1 ) << "Wrong solution data" << std::endl;
//         KRATOS_ERROR_IF( PMMG_Chk_meshData(mParMmgMesh, mParMmgDisp) != 1 ) << "Wrong displacement data" << std::endl;
//     } else {
//         KRATOS_ERROR_IF( PMMG_Chk_meshData(mParMmgMesh, mParMmgSol) != 1 ) << "Wrong mesh data" << std::endl;
//     }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::InputMesh(const std::string& rOutputName)
{
     const std::string mesh_name = rOutputName + ".mesh";
     const char* mesh_file = mesh_name.c_str();

    // // a)  Give the ouptut mesh name using PMMG_Set_inputMeshName
     PMMG_Set_inputMeshName(mParMmgMesh, mesh_file);

    // // b) function calling
     KRATOS_INFO_IF("ParMmgUtilities", PMMG_loadMesh_distributed(mParMmgMesh, mesh_file) != 1) << "UNABLE TO LOAD MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::InputSol(const std::string& rInputName)
{
     const std::string sol_name = rInputName + ".sol";
     const char* sol_file = sol_name.c_str();

    // // a)  Give the ouptut mesh name using PMMG_Set_inputSolName (by default, the mesh is saved in the "mesh.o.mesh") file
     PMMG_Set_inputSolsName(mParMmgMesh,  sol_file);

    // // b) function calling
     KRATOS_INFO_IF("ParMmgUtilities", PMMG_loadSol_centralized(mParMmgMesh,  sol_file) != 1) << "UNABLE TO READ METRIC" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::OutputMesh(const std::string& rOutputName)
{
     const std::string mesh_name = rOutputName + ".mesh";
     const char* mesh_file = mesh_name.c_str();

    // // a)  Give the ouptut mesh name using PMMG_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh") file
     PMMG_Set_outputMeshName(mParMmgMesh,mesh_file);

    // // b) function calling
     KRATOS_INFO_IF("ParMmgUtilities", PMMG_saveMesh_centralized(mParMmgMesh,mesh_file) != 1) << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::OutputSol(const std::string& rOutputName)
{
     const std::string sol_name = rOutputName + ".sol";
     const char* sol_file = sol_name.c_str();

    // // a)  Give the ouptut sol name using PMMG_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file
     PMMG_Set_outputSolsName(mParMmgMesh,  sol_file);

    // // b) Function calling
     KRATOS_INFO_IF("ParMmgUtilities", PMMG_saveMet_centralized(mParMmgMesh, sol_file) != 1)<< "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::OutputDisplacement(const std::string& rOutputName)
{
    KRATOS_ERROR << "OutputDisplacement Not Yet Implemented" << std::endl;
//     const std::string sol_name = rOutputName + ".disp.sol";
//     const char* sol_file = sol_name.c_str();

    // // a)  Give the ouptut sol name using PMMG_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file
//     PMMG_Set_outputSolName(mParMmgMesh, mParMmgDisp, sol_file);

    // // b) Function calling
//     KRATOS_INFO_IF("ParMmgUtilities", PMMG_saveSol(mParMmgMesh,mParMmgDisp, sol_file) != 1)<< "UNABLE TO SAVE DISPLACEMENT" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::OutputReferenceEntitities(
    const std::string& rOutputName,
    const std::unordered_map<IndexType,Condition::Pointer>& rRefCondition,
    const std::unordered_map<IndexType,Element::Pointer>& rRefElement
    )
{
    /* ELEMENTS */
    std::string element_name;
    Parameters elem_reference_json;
    for (auto& r_elem : rRefElement) {
        CompareElementsAndConditionsUtility::GetRegisteredName(*(r_elem.second), element_name);
        const std::string name = std::to_string(r_elem.first);
        elem_reference_json.AddEmptyValue(name);
        elem_reference_json[name].SetString(element_name);
    }

    const std::string& r_elem_json_text = elem_reference_json.PrettyPrintJsonString();

    std::filebuf elem_buffer;
    elem_buffer.open(rOutputName + ".elem.ref.json",std::ios::out);
    std::ostream elem_os(&elem_buffer);
    elem_os << r_elem_json_text;
    elem_buffer.close();

    /* CONDITIONS */
    std::string condition_name;
    Parameters cond_reference_json;
    for (auto& r_cond : rRefCondition) {
        CompareElementsAndConditionsUtility::GetRegisteredName(*(r_cond.second), condition_name);
        const std::string name = std::to_string(r_cond.first);
        cond_reference_json.AddEmptyValue(name);
        cond_reference_json[name].SetString(condition_name);
    }

    const std::string& r_cond_json_text = cond_reference_json.PrettyPrintJsonString();

    std::filebuf cond_buffer;
    cond_buffer.open(rOutputName + ".cond.ref.json",std::ios::out);
    std::ostream cond_os(&cond_buffer);
    cond_os << r_cond_json_text;
    cond_buffer.close();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::FreeAll()
{
    PMMG_Free_all(PMMG_ARG_start, PMMG_ARG_ppParMesh, &mParMmgMesh, PMMG_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::PMMGLibCallMetric(Parameters ConfigurationParameters)
{
    KRATOS_TRY;

    /* Advanced configurations */
    // Global hausdorff value (default value = 0.01) applied on the whole boundary
    if (ConfigurationParameters["advanced_parameters"]["force_hausdorff_value"].GetBool()) {
        if ( PMMG_Set_dparameter(mParMmgMesh,PMMG_DPARAM_hausd, ConfigurationParameters["advanced_parameters"]["hausdorff_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the Hausdorff parameter" << std::endl;
    }

    // Avoid/allow point relocation
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_nomove, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_move_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to fix the nodes" << std::endl;

    // Avoid/allow surface modifications
    if (static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_surf_mesh"].GetBool()) == 1) KRATOS_ERROR << "Trying to do surface" << std::endl;
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_nosurf, 1) != 1 )
        KRATOS_ERROR << "Unable to set no surfacic modifications" << std::endl;

    // Don't insert nodes on mesh
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_noinsert, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_insert_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no insertion/suppression point" << std::endl;

    // Don't swap mesh
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_noswap, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_swap_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // Number Of iterations
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_niter, static_cast<int>(ConfigurationParameters["advanced_parameters"]["niter"].GetInt())) != 1 )
        KRATOS_ERROR << "Unable to set number of remeshing iterations" << std::endl;

    // Mesh Size
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_meshSize, static_cast<int>(ConfigurationParameters["advanced_parameters"]["meshSize"].GetInt())) != 1 )
        KRATOS_ERROR << "Unable to set target mesh size of Mmg" << std::endl;

    // Mesh Ratio
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_metisRatio, static_cast<int>(ConfigurationParameters["advanced_parameters"]["metisRatio"].GetInt())) != 1 )
        KRATOS_ERROR << "Unable to set wanted ratio # mesh / # metis super nodes" << std::endl;

    // hgradreq ( To Be added)
    // if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_hgradreq, static_cast<int>(ConfigurationParameters["advanced_parameters"]["hgradreq"].GetDouble())) != 1 )
    //     KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // API Mode: Forced to "1" all the time
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_APImode, static_cast<int>(PMMG_APIDISTRIB_nodes)) != 1 )
        KRATOS_ERROR << "Unable to set initialize parallel library through interface faces or nodes" << std::endl;

    // Node global numbering: Forced to "1" all the time
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_nodeGloNum,1) != 1 )
        KRATOS_ERROR << "Unable to set initialize parallel library with node global numbering" << std::endl;

    // Set the angle detection
    const bool deactivate_detect_angle = ConfigurationParameters["advanced_parameters"]["deactivate_detect_angle"].GetBool();
    if ( deactivate_detect_angle) {
        if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_angle, static_cast<int>(!deactivate_detect_angle)) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection on" << std::endl;
    }

    // Set the gradation
    if (ConfigurationParameters["advanced_parameters"]["force_gradation_value"].GetBool()) {
        if ( PMMG_Set_dparameter(mParMmgMesh,PMMG_DPARAM_hgrad, ConfigurationParameters["advanced_parameters"]["gradation_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set gradation" << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_min"].GetBool()) {
        if ( PMMG_Set_dparameter(mParMmgMesh,PMMG_DPARAM_hmin, ConfigurationParameters["force_sizes"]["minimal_size"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the minimal edge size " << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_max"].GetBool()) {
        if ( PMMG_Set_dparameter(mParMmgMesh,PMMG_DPARAM_hmax, ConfigurationParameters["force_sizes"]["maximal_size"].GetDouble()) != 1 ) {
            KRATOS_ERROR << "Unable to set the maximal edge size " << std::endl;
        }
    }

    // Actually computing remesh
    int ier;
    KRATOS_INFO_IF("", mEchoLevel > 0) << "HEY 1.1" << std::endl;

    ier = PMMG_parmmglib_distributed(mParMmgMesh);

    KRATOS_INFO_IF("", mEchoLevel > 0) << "HEY 1.2" << std::endl;

    if ( ier == PMMG_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF PMMG3DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == PMMG_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF PMMG3DLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const IndexType Color,
    const IndexType Index
    )
{
    KRATOS_ERROR_IF( PMMG_Set_vertex(mParMmgMesh, X, Y, Z, Color, Index) != 1 ) << "Unable to set vertex" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetConditions(
    GeometryType& rGeometry,
    const IndexType Color,
    const IndexType Index
    )
{
    if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point3D) // Point
        KRATOS_ERROR << "ERROR:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2) { // Line
        KRATOS_ERROR << "Kratos_Line3D2 remeshing pending to be implemented" << std::endl;
    } else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {// Triangle
        const IndexType id_1 = rGeometry[0].Id(); // First node Id
        const IndexType id_2 = rGeometry[1].Id(); // Second node Id
        const IndexType id_3 = rGeometry[2].Id(); // Third node Id

        KRATOS_ERROR_IF( PMMG_Set_triangle(mParMmgMesh, id_1, id_2, id_3, Color, Index) != 1 ) << "Unable to set triangle" << std::endl;

        // Set fixed boundary
        bool blocked_1 = false;
        if (rGeometry[0].IsDefined(BLOCKED))
            blocked_1 = rGeometry[0].Is(BLOCKED);
        bool blocked_2 = false;
        if (rGeometry[1].IsDefined(BLOCKED))
            blocked_2 = rGeometry[1].Is(BLOCKED);
        bool blocked_3 = false;
        if (rGeometry[2].IsDefined(BLOCKED))
            blocked_3 = rGeometry[2].Is(BLOCKED);

        if (blocked_1 && blocked_2 && blocked_3) BlockCondition(Index);
    } else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) { // Quadrilaterals
        KRATOS_ERROR << "Not yet implemented" << std::endl;
        //KRATOS_ERROR_IF( PMMG_Set_quadrilateral(mParMmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 ) << "Unable to set quadrilateral" << std::endl;
    } else {
        const SizeType size_geometry = rGeometry.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << rGeometry.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetElements(
    GeometryType& rGeometry,
    const IndexType Color,
    const IndexType Index
    )
{
    const IndexType id_1 = local_node_id[rGeometry[0].Id()]; // First node Id
    const IndexType id_2 = local_node_id[rGeometry[1].Id()]; // Second node Id
    const IndexType id_3 = local_node_id[rGeometry[2].Id()]; // Third node Id
    const IndexType id_4 = local_node_id[rGeometry[3].Id()]; // Fourth node Id

    if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) { // Tetrahedron
        KRATOS_ERROR_IF( PMMG_Set_tetrahedron(mParMmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 ) << "Unable to set tetrahedron" << std::endl;
    } else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) { // Prisms
        const SizeType size_geometry = rGeometry.size();
        KRATOS_ERROR << "ERROR: PRISM NON IMPLEMENTED IN THE LIBRARY " << size_geometry << std::endl;
    } else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) { // Hexaedron
        const SizeType size_geometry = rGeometry.size();
        KRATOS_ERROR << "ERROR: HEXAEDRON NON IMPLEMENTED IN THE LIBRARY " << size_geometry << std::endl;
    } else {
        const SizeType size_geometry = rGeometry.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetMetricScalar(
    const double Metric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( PMMG_Set_scalarMet( mParMmgMesh, Metric, local_node_id[NodeId]) != 1 ) << "Unable to set scalar metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetMetricVector(
    const array_1d<double, 3>& rMetric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( PMMG_Set_vectorMet( mParMmgMesh, rMetric[0], rMetric[1], rMetric[2], local_node_id[NodeId]) != 1 ) << "Unable to set vector metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetMetricTensor(
    const array_1d<double, 6>& rMetric,
    const IndexType NodeId
    )
{
    // KRATOS_ERROR << "Not yet implemented" << std::endl;
    // The order is XX, XY, XZ, YY, YZ, ZZ
    KRATOS_ERROR_IF( PMMG_Set_tensorMet( mParMmgMesh, rMetric[0], rMetric[3], rMetric[5], rMetric[1], rMetric[4], rMetric[2], local_node_id[NodeId]) != 1 ) << "Unable to set tensor metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetDisplacementVector(
    const array_1d<double, 3>& rDisplacement,
    const IndexType NodeId
    )
{
    KRATOS_ERROR << "Not yet implemented" << std::endl;
    // KRATOS_ERROR_IF( PMMG_Set_vectorSol(mParMmgDisp, rDisplacement[0], rDisplacement[1], rDisplacement[2], NodeId) != 1 ) << "Unable to set vector displacement" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::GetMetricScalar(double& rMetric)
{

    KRATOS_ERROR_IF( PMMG_Get_scalarMet( mParMmgMesh, &rMetric) != 1 ) << "Unable to get scalar metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::GetMetricVector(array_1d<double, 3>& rMetric)
{
    KRATOS_ERROR_IF( PMMG_Get_vectorMet( mParMmgMesh, &rMetric[0], &rMetric[1], &rMetric[2]) != 1 ) << "Unable to get vector metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::GetMetricTensor(array_1d<double, 6>& rMetric)
{
    // The order is XX, XY, XZ, YY, YZ, ZZ
    KRATOS_ERROR_IF( PMMG_Get_tensorMet( mParMmgMesh, &rMetric[0], &rMetric[3], &rMetric[5], &rMetric[1], &rMetric[4], &rMetric[2]) != 1 ) << "Unable to get tensor metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::GetDisplacementVector(array_1d<double, 3>& rDisplacement)
{
    KRATOS_ERROR << "Not yet implemented" << std::endl;
    // KRATOS_ERROR_IF( PMMG_Get_vectorSol(mParMmgDisp, &rDisplacement[0], &rDisplacement[1], &rDisplacement[2]) != 1 ) << "Unable to get vector displacement" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::ReorderAllIds(ModelPart& rModelPart)
{
    // Iterate over nodes
    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    for(IndexType i = 0; i < r_nodes_array.size(); ++i)
        (it_node_begin + i)->SetId(i + 1);

    // Iterate over conditions
    auto& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    for(IndexType i = 0; i < r_conditions_array.size(); ++i)
        (it_cond_begin + i)->SetId(i + 1);

    // Iterate over elements
    auto& r_elements_array = rModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();
    for(IndexType i = 0; i < r_elements_array.size(); ++i)
        (it_elem_begin + i)->SetId(i + 1);
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::GenerateMeshDataFromModelPart(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,std::vector<std::string>>& rColors,
    ColorsMapType& rColorMapCondition,
    ColorsMapType& rColorMapElement,
    const FrameworkEulerLagrange Framework,
    const bool CollapsePrismElements
    )
{
    // Before computing colors we do some check and throw a warning to get the user informed
    const std::vector<std::string> sub_model_part_names = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(rModelPart);

    for (auto sub_model_part_name : sub_model_part_names) {
        ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(rModelPart, sub_model_part_name);

        KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 0 && (r_sub_model_part.NumberOfNodes() > 0 && (r_sub_model_part.NumberOfConditions() == 0 && r_sub_model_part.NumberOfElements() == 0))) <<
        "The submodelpart: " << sub_model_part_name << " contains only nodes and no geometries (conditions/elements)." << std::endl <<
        "It is not guaranteed that the submodelpart will be preserved." << std::endl <<
        "PLEASE: Add some \"dummy\" conditions to the submodelpart to preserve it" << std::endl;
    }

    /////////* MESH FILE */////////
    // Build mesh in MMG5 format //

    // Iterate over components
    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    auto& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    auto& r_elements_array = rModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();

    // The following nodes will be remeshed
    std::unordered_set<IndexType> remeshed_nodes;

    /* Manually set of the mesh */
    PMMGMeshInfo<TPMMGLibrary> mmg_mesh_info;
    if (TPMMGLibrary == PMMGLibrary::PMMG3D) { // 3D
        /* Conditions */
        std::size_t num_tri = 0, num_quad = 0;
        for(IndexType i = 0; i < r_conditions_array.size(); ++i) {
            auto it_cond = it_cond_begin + i;

            if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) { // Triangles
                for (auto& r_node : it_cond->GetGeometry())
                    remeshed_nodes.insert(r_node.Id());
                num_tri += 1;
            } else if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) { // Quadrilaterals
                for (auto& r_node : it_cond->GetGeometry())
                    remeshed_nodes.insert(r_node.Id());
                num_quad += 1;
            } else {
                it_cond->Set(OLD_ENTITY, true);
                KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "WARNING:: YOUR GEOMETRY CONTAINS A CONDITION WITH " << it_cond->GetGeometry().PointsNumber() <<" NODES THAT CAN NOT BE REMESHED" << std::endl;
            }
        }

        mmg_mesh_info.NumberOfTriangles = num_tri;
        mmg_mesh_info.NumberOfQuadrilaterals = num_quad;

        KRATOS_INFO_IF("ParMmgUtilities", ((num_tri + num_quad) < r_conditions_array.size()) && mEchoLevel > 0) <<
        "Number of Conditions: " << r_conditions_array.size() << " Number of Triangles: " << num_tri << " Number of Quadrilaterals: " << num_quad << std::endl;

        /* Elements */
        std::size_t num_tetra = 0, num_prisms = 0;
        for(IndexType i = 0; i < r_elements_array.size(); ++i) {
            auto it_elem = it_elem_begin + i;

            if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) { // Tetrahedron
                for (auto& r_node : it_elem->GetGeometry())
                    remeshed_nodes.insert(r_node.Id());
                num_tetra += 1;
            } else if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) { // Prisms
                if (CollapsePrismElements) {
                    KRATOS_INFO_IF("ParMmgUtilities", mEchoLevel > 1) << "Prismatic element " << it_elem->Id() << " will be collapsed to a triangle" << std::endl;
                } else {
                    for (auto& r_node : it_elem->GetGeometry())
                        remeshed_nodes.insert(r_node.Id());
                    num_prisms += 1;
                }
            } else {
                it_elem->Set(OLD_ENTITY, true);
                KRATOS_WARNING_IF("ParMmgUtilities", mEchoLevel > 1) << "WARNING:: YOUR GEOMETRY CONTAINS AN ELEMENT WITH " << it_elem->GetGeometry().PointsNumber() <<" NODES CAN NOT BE REMESHED" << std::endl;
            }
        }

        mmg_mesh_info.NumberOfTetrahedra = num_tetra;
        mmg_mesh_info.NumberOfPrism = num_prisms;

        KRATOS_INFO_IF("ParMmgUtilities", ((num_tetra + num_prisms) < r_elements_array.size()) && mEchoLevel > 0) <<
        "Number of Elements: " << r_elements_array.size() << " Number of Tetrahedron: " << num_tetra << " Number of Prisms: " << num_prisms << std::endl;
    }

    // Set flag on nodes
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;
        if (remeshed_nodes.find(it_node->Id()) == remeshed_nodes.end()) {
            it_node->Set(OLD_ENTITY, true);
        }
    }

    mmg_mesh_info.NumberOfNodes = remeshed_nodes.size();
    SetMeshSize(mmg_mesh_info);

    // We reorder the ids to avoid conflicts with the rest (using as reference the OLD_ENTITY)
    IndexType counter_to_remesh = 0;
    #pragma omp parallel for reduction(+: counter_to_remesh)
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            ++counter_to_remesh;
        }
    }

    local_node_id = std::map<int, int>();
    local_elem_id = std::map<int, int>();
    local_cond_id = std::map<int, int>();

    // RESETING THE ID OF THE NODES (important for non consecutive meshes)
    IndexType counter_remesh = 0;
    IndexType counter_not_remesh = counter_to_remesh + 1;
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        local_node_id[it_node->Id()] = ++counter_remesh;
        // const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
        // if (!old_entity) {
        //     it_node->SetId(counter_remesh);
        //     ++counter_remesh;
        // } else {
        //     it_node->SetId(counter_not_remesh);
        //     ++counter_not_remesh;
        // }
    }
    counter_to_remesh = 0;
    #pragma omp parallel for reduction(+: counter_to_remesh)
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = it_cond_begin + i;

        const bool old_entity = it_cond->IsDefined(OLD_ENTITY) ? it_cond->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            ++counter_to_remesh;
        }
    }
    // RESETING THE ID OF THE CONDITIONS (important for non consecutive meshes)
    counter_remesh = 0;
    counter_not_remesh = counter_to_remesh + 1;
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = it_cond_begin + i;

        local_cond_id[it_cond->Id()] = ++counter_remesh;
        // const bool old_entity = it_cond->IsDefined(OLD_ENTITY) ? it_cond->Is(OLD_ENTITY) : false;
        // if (!old_entity) {
        //     it_cond->SetId(counter_remesh);
        //     ++counter_remesh;
        // } else {
        //     it_cond->SetId(counter_not_remesh);
        //     ++counter_not_remesh;
        // }
    }
    counter_to_remesh = 0;
    #pragma omp parallel for reduction(+: counter_to_remesh)
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;

        const bool old_entity = it_elem->IsDefined(OLD_ENTITY) ? it_elem->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            ++counter_to_remesh;
        }
    }
    // RESETING THE ID OF THE ELEMENTS (important for non consecutive meshes)
    counter_remesh = 0;
    counter_not_remesh = counter_to_remesh + 1;
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;

        local_elem_id[it_elem->Id()] = ++counter_remesh;
        // const bool old_entity = it_elem->IsDefined(OLD_ENTITY) ? it_elem->Is(OLD_ENTITY) : false;
        // if (!old_entity) {
        //     it_elem->SetId(counter_remesh);
        //     ++counter_remesh;
        // } else {
        //     it_elem->SetId(counter_not_remesh);
        //     ++counter_not_remesh;
        // }
    }

    // Now we compute the colors
    rColors.clear();
    ColorsMapType nodes_colors, cond_colors, elem_colors;
    AssignUniqueModelPartCollectionTagUtility model_part_collections(rModelPart);
    model_part_collections.ComputeTags(nodes_colors, cond_colors, elem_colors, rColors);

    // The ISOSURFACE has some reserved Ids. We reassign
    if (mDiscretization == DiscretizationOption::ISOSURFACE) {
        // Create auxiliar model part
        if (!rModelPart.HasSubModelPart("SKIN_ISOSURFACE")) {
            rModelPart.CreateSubModelPart("SKIN_ISOSURFACE");
        }

        // Do some checks
        bool id_2_exists = false;
        bool id_3_exists = false;
        bool id_10_exists = false;
        IndexType max_index = 0;
        for (auto& r_color : rColors) {
            const IndexType index = r_color.first;
            if (index == 2) {
                id_2_exists = true;
            } else if (index == 3) {
                id_3_exists = true;
            } else if (index == 10) {
                id_10_exists = true;
            }
            if (index > max_index)
                max_index = index;
        }

        // Identify the submodelparts with elements
        std::unordered_set<std::string> auxiliar_set_elements;
        // We build a set for all the model parts containing elements
        for (auto& r_elem_color : elem_colors) {
            const IndexType color = r_elem_color.second;
            const auto& r_sub_model_parts_names = rColors[color];
            auxiliar_set_elements.insert(r_sub_model_parts_names.begin(), r_sub_model_parts_names.end());
        }
        std::vector<std::string> auxiliar_vector_elements(auxiliar_set_elements.size());
        std::copy(auxiliar_set_elements.begin(), auxiliar_set_elements.end(), std::back_inserter(auxiliar_vector_elements));

        // Move the map to the end
        if (id_2_exists) {
            // New group
            rColors.insert(IndexStringVectorPairType(max_index + 1, rColors[2]));
            rColors.erase(2);

            // Reassign
            for (auto& r_nodes_color : nodes_colors) {
                IndexType& r_color = r_nodes_color.second;
                if (r_color == 2) {
                    r_color = max_index + 1;
                }
            }
            for (auto& r_cond_color : cond_colors) {
                IndexType& r_color = r_cond_color.second;
                if (r_color == 2) {
                    r_color = max_index + 1;
                }
            }
            for (auto& r_elem_color : elem_colors) {
                IndexType& r_color = r_elem_color.second;
                if (r_color == 2) {
                    r_color = max_index + 1;
                }
            }

        }
        if (id_3_exists) {
            // New group
            rColors.insert(IndexStringVectorPairType(max_index + 2, rColors[3]));
            rColors.erase(3);

            // Reassign
            for (auto& r_nodes_color : nodes_colors) {
                IndexType& r_color = r_nodes_color.second;
                if (r_color == 3) {
                    r_color = max_index + 2;
                }
            }
            for (auto& r_cond_color : cond_colors) {
                IndexType& r_color = r_cond_color.second;
                if (r_color == 3) {
                    r_color = max_index + 2;
                }
            }
            for (auto& r_elem_color : elem_colors) {
                IndexType& r_color = r_elem_color.second;
                if (r_color == 3) {
                    r_color = max_index + 2;
                }
            }

        }
        if (id_10_exists) {
            // New group
            rColors.insert(IndexStringVectorPairType(max_index + 3, rColors[10]));
            rColors.erase(10);

            // Reassign
            for (auto& r_nodes_color : nodes_colors) {
                IndexType& r_color = r_nodes_color.second;
                if (r_color == 10) {
                    r_color = max_index + 3;
                }
            }
            for (auto& r_cond_color : cond_colors) {
                IndexType& r_color = r_cond_color.second;
                if (r_color == 10) {
                    r_color = max_index + 3;
                }
            }
            for (auto& r_elem_color : elem_colors) {
                IndexType& r_color = r_elem_color.second;
                if (r_color == 10) {
                    r_color = max_index + 3;
                }
            }
        }

        // Fill the model parts
        rColors.insert(IndexStringVectorPairType(2, auxiliar_vector_elements));
        rColors.insert(IndexStringVectorPairType(3, auxiliar_vector_elements));
        std::vector<std::string> auxiliar_name_vector (1, "SKIN_ISOSURFACE");
        rColors.insert(IndexStringVectorPairType(10, auxiliar_name_vector));
    }

    /* Nodes */
    #pragma omp parallel for firstprivate(nodes_colors)
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            const array_1d<double, 3>& r_coordinates = Framework == FrameworkEulerLagrange::LAGRANGIAN ? it_node->GetInitialPosition() : it_node->Coordinates();
            //SetNodes(r_coordinates[0], r_coordinates[1], r_coordinates[2], nodes_colors[it_node->Id()], it_node->Id());
            SetNodes(r_coordinates[0], r_coordinates[1], r_coordinates[2], nodes_colors[it_node->Id()], local_node_id[it_node->Id()]);

            bool blocked = false;
            if (it_node->IsDefined(BLOCKED))
                blocked = it_node->Is(BLOCKED);
            if (blocked)
                BlockNode(it_node->Id());
        }
    }

    /* Conditions */
    #pragma omp parallel for firstprivate(cond_colors)
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)  {
        auto it_cond = it_cond_begin + i;

        const bool old_entity = it_cond->IsDefined(OLD_ENTITY) ? it_cond->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            SetConditions(it_cond->GetGeometry(), cond_colors[it_cond->Id()], local_cond_id[it_cond->Id()]);

            bool blocked = false;
            if (it_cond->IsDefined(BLOCKED))
                blocked = it_cond->Is(BLOCKED);
            if (blocked)
                BlockCondition(it_cond->Id());
        }
    }

    /* Elements */
    #pragma omp parallel for firstprivate(elem_colors)
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;

        const bool old_entity = it_elem->IsDefined(OLD_ENTITY) ? it_elem->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            SetElements(it_elem->GetGeometry(), elem_colors[it_elem->Id()], local_elem_id[it_elem->Id()]);

            bool blocked = false;
            if (it_elem->IsDefined(BLOCKED))
                blocked = it_elem->Is(BLOCKED);
            if (blocked)
                BlockElement(it_elem->Id());
        }
    }

    // Create auxiliar colors maps
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)  {
        auto it_cond = it_cond_begin + i;
        const IndexType cond_id = it_cond->Id();
        const IndexType color = cond_colors[cond_id];
        if ((rColorMapCondition.find(color) == rColorMapCondition.end()))
            rColorMapCondition.insert (IndexPairType(color,cond_id));
    }
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;
        const IndexType elem_id = it_elem->Id();
        const IndexType color = elem_colors[elem_id];
        if ((rColorMapElement.find(color) == rColorMapElement.end()))
            rColorMapElement.insert (IndexPairType(color,elem_id));
    }

    // Add missing entities
    Model& r_model = rModelPart.GetModel();
    for (auto& r_color : rColors) {
        const IndexType color = r_color.first;
        if (color != 0 && r_color.second.size() == 1) { // Not including main model part, and adding only simple model parts
            for (auto& r_sub_model_part_name : r_color.second) {
                ModelPart& r_sub_model_part = r_model.GetModelPart(rModelPart.Name() + "." + r_sub_model_part_name);
                if ((rColorMapCondition.find(color) == rColorMapCondition.end())) {
                    if (r_sub_model_part.NumberOfConditions() > 0) {
                        const IndexType cond_id = r_sub_model_part.Conditions().begin()->Id();
                        rColorMapCondition.insert (IndexPairType(color,cond_id));
                    }
                }
                if ((rColorMapElement.find(color) == rColorMapElement.end())) {
                    if (r_sub_model_part.NumberOfElements() > 0) {
                        const IndexType elem_id = r_sub_model_part.Elements().begin()->Id();
                        rColorMapElement.insert (IndexPairType(color,elem_id));
                    }
                }
            }
        }
    }
}

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::GenerateParallelInterfaces(
    ModelPart& rModelPart
    )
{
    auto& neighbour_indices = rModelPart.GetCommunicator().NeighbourIndices();
    std::size_t icomm,ncomm;

    ncomm = 0;
    for(std::size_t i = 0; i < neighbour_indices.size(); i++)
        if ((neighbour_indices[i]) >= 0) ncomm++;
    PMMG_Set_numberOfNodeCommunicators(mParMmgMesh, ncomm);

    icomm = 0;
    for(std::size_t i = 0; i < neighbour_indices.size(); i++) {
        std::vector<int> globalId(0);
        std::vector<int> localId(0);
        globalId.reserve(rModelPart.GetCommunicator().LocalMesh(i).NumberOfNodes() + rModelPart.GetCommunicator().GhostMesh(i).NumberOfNodes());
        localId.reserve(rModelPart.GetCommunicator().LocalMesh(i).NumberOfNodes() + rModelPart.GetCommunicator().GhostMesh(i).NumberOfNodes());

        if ((neighbour_indices[i]) >= 0) {
            for(auto& node: rModelPart.GetCommunicator().LocalMesh(i).Nodes()) {
                globalId.push_back(node.Id());
                localId.push_back(local_node_id[node.Id()]);
            }

            for(auto& node: rModelPart.GetCommunicator().GhostMesh(i).Nodes()) {
                globalId.push_back(node.Id());
                localId.push_back(local_node_id[node.Id()]);
            }
            PMMG_Set_ithNodeCommunicatorSize(mParMmgMesh, icomm, rModelPart.GetCommunicator().NeighbourIndices()[i], rModelPart.GetCommunicator().LocalMesh(i).NumberOfNodes() + rModelPart.GetCommunicator().GhostMesh(i).NumberOfNodes());
            PMMG_Set_ithNodeCommunicator_nodes(mParMmgMesh, icomm++, localId.data(), globalId.data(), 1); // Last parameters shoould be 1
        }
    }
}

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::PrintParallelInterfaces(
    ModelPart& rModelPart
    )
{
    auto& neighbour_indices = rModelPart.GetCommunicator().NeighbourIndices();

    for(std::size_t i = 0; i < neighbour_indices.size(); i++) {
        std::vector<int> globalId(0);
        std::vector<int> localId(0);
        globalId.reserve(rModelPart.GetCommunicator().LocalMesh(i).NumberOfNodes() + rModelPart.GetCommunicator().GhostMesh(i).NumberOfNodes());
        localId.reserve(rModelPart.GetCommunicator().LocalMesh(i).NumberOfNodes() + rModelPart.GetCommunicator().GhostMesh(i).NumberOfNodes());

        if ((neighbour_indices[i]) >= 0) {
            for(auto& node: rModelPart.GetCommunicator().LocalMesh(i).Nodes()) {
                globalId.push_back(node.Id());
                localId.push_back(local_node_id[node.Id()]);
            }

            for(auto& node: rModelPart.GetCommunicator().GhostMesh(i).Nodes()) {
                globalId.push_back(node.Id());
                localId.push_back(local_node_id[node.Id()]);
            }
        }

        std::cout << "List of nodes from process " << rModelPart.GetCommunicator().GetDataCommunicator().Rank() << " to process " << neighbour_indices[i] << std::endl;
        KRATOS_WATCH(globalId)
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::GenerateReferenceMaps(
    ModelPart& rModelPart,
    const ColorsMapType& rColorMapCondition,
    const ColorsMapType& rColorMapElement,
    std::unordered_map<IndexType,Condition::Pointer>& rRefCondition,
    std::unordered_map<IndexType,Element::Pointer>& rRefElement
    )
{
    /* We clone the first condition and element of each type (we will assume that each sub model part has just one kind of condition, in my opinion it is quite reccomended to create more than one sub model part if you have more than one element or condition) */

    auto& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    auto& r_elements_array = rModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();

    if (r_conditions_array.size() > 0) {
        const std::string type_name = (Dimension == 2) ? "Condition2D2N" : (TPMMGLibrary == PMMGLibrary::PMMG3D) ? "SurfaceCondition3D3N" : "Condition3D2N";
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(type_name);
        rRefCondition[0] = r_clone_condition.Create(0, r_clone_condition.pGetGeometry(), it_cond_begin->pGetProperties());
    }
    if (r_elements_array.size() > 0) {
        rRefElement[0] = it_elem_begin->Create(0, it_elem_begin->GetGeometry(), it_elem_begin->pGetProperties());
    }

    // Now we add the reference elements and conditions
    for (auto& ref_cond : rColorMapCondition) {
        Condition::Pointer p_cond = rModelPart.pGetCondition(ref_cond.second);
        rRefCondition[ref_cond.first] = p_cond->Create(0, p_cond->GetGeometry(), p_cond->pGetProperties());
    }
    for (auto& ref_elem : rColorMapElement) {
        Element::Pointer p_elem = rModelPart.pGetElement(ref_elem.second);
        rRefElement[ref_elem.first] = p_elem->Create(0, p_elem->GetGeometry(), p_elem->pGetProperties());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::GenerateSolDataFromModelPart(ModelPart& rModelPart)
{
    // Iterate in the nodes
    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    // Set size of the solution
    const Variable<TensorArrayType>& r_tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_" + std::to_string(Dimension)+"D");

    if (it_node_begin->Has(r_tensor_variable)) {
        SetSolSizeTensor(r_nodes_array.size());
    } else {
        SetSolSizeScalar(r_nodes_array.size());
    }

    // In case of considering metric tensor
    if (it_node_begin->Has(r_tensor_variable)) {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;

            const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
            if (!old_entity) {
                KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(r_tensor_variable)) << "METRIC_TENSOR_" + std::to_string(Dimension) + "D  not defined for node " << it_node->Id() << std::endl;

                // We get the metric
                 const TensorArrayType& r_metric = it_node->GetValue(r_tensor_variable);

                // We set the metric
                SetMetricTensor(r_metric, it_node->Id());
            }
        }
    } else {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;

            const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
            if (!old_entity) {
                KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(METRIC_SCALAR)) << "METRIC_SCALAR not defined for node " << it_node->Id() << std::endl;

                // We get the metric
                const double metric = it_node->GetValue(METRIC_SCALAR);

                // We set the metric
                SetMetricScalar(metric, it_node->Id());
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::GenerateDisplacementDataFromModelPart(ModelPart& rModelPart)
{
    // Iterate in the nodes
    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    // Set size of the solution
    SetDispSizeVector(r_nodes_array.size());

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            // We get the displacement
            const array_1d<double, 3>& r_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            // We set the displacement
            SetDisplacementVector(r_displacement, it_node->Id());
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::WriteMeshDataToModelPart(
    ModelPart& rModelPart,
    const std::unordered_map<IndexType,std::vector<std::string>>& rColors,
    const NodeType::DofsContainerType& rDofs,
    const PMMGMeshInfo<TPMMGLibrary>& rPMMGMeshInfo,
    std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
    std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement
    )
{
    int **idx_node_loc, **owner, **idx_node_glob;
    int nunique, ntot;
    int n_node_comm,n_face_comm,*nitem_node_comm,*nitem_face_comm;
    int *color_node, *color_face,**face_owner,nunique_face,ntot_face;
    int ier = PMMG_Get_numberOfNodeCommunicators(mParMmgMesh, &n_node_comm);

    color_node      = (int *) malloc(n_node_comm*sizeof(int));
    nitem_node_comm = (int *) malloc(n_node_comm*sizeof(int));
    for(IndexType icomm = 0; icomm < n_node_comm; icomm++ )
        ier = PMMG_Get_ithNodeCommunicatorSize(mParMmgMesh, icomm,
                                            &color_node[icomm],
                                            &nitem_node_comm[icomm]);

    idx_node_loc  = (int **) malloc(n_node_comm*sizeof(int *));
    owner  = (int **) malloc(n_node_comm*sizeof(int *));
    idx_node_glob  = (int **) malloc(n_node_comm*sizeof(int *));
    for(IndexType icomm = 0; icomm < n_node_comm; icomm++ ) {
        idx_node_loc[icomm]  = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
        owner[icomm]  = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
        idx_node_glob[icomm]  = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));

    }

    int err_node_comm = PMMG_Get_NodeCommunicator_nodes(mParMmgMesh, idx_node_loc);
    int err_owners = PMMG_Get_NodeCommunicator_owners(mParMmgMesh, owner, idx_node_glob, &nunique, &ntot);

    const int rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
    const int size = rModelPart.GetCommunicator().GetDataCommunicator().Size();
    std::vector<int> array_of_local_nodes(size,0);
    std::vector<int> reduced_array_of_local_nodes(size,0);
    if (rank!=size) {
        array_of_local_nodes[rank+1]=rPMMGMeshInfo.NumberOfNodes-nunique;
    }

    rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(array_of_local_nodes, reduced_array_of_local_nodes);

    KRATOS_WATCH(array_of_local_nodes)
    KRATOS_WATCH(reduced_array_of_local_nodes)
    for (int i =1; i<size;i++){
        reduced_array_of_local_nodes[i] += reduced_array_of_local_nodes[i-1];
    }
    KRATOS_WATCH(reduced_array_of_local_nodes)


    std::map<int,int> local_to_partition_index;

    int errglonum;
    for (IndexType i_node = 1; i_node <= rPMMGMeshInfo.NumberOfNodes; ++i_node) {
        errglonum = PMMG_Get_vertexGloNum(mParMmgMesh,&mLocalToGlobal[i_node],&local_to_partition_index[i_node]);
    }
//    for(IndexType icomm = 0; icomm < n_node_comm; icomm++ ) {
//        for(IndexType inode = 0; inode < nitem_node_comm[icomm]; inode++ ) {
//            idx_node_glob[icomm][inode] = mLocalToGlobal[idx_node_loc[icomm][inode]];
//        }
//    }
//    printf("MYRANK %d\n",mParMmgMesh->myrank);
//    PMMG_printCommunicator(mParMmgMesh, PMMG_APIDISTRIB_nodes, idx_node_loc, idx_node_glob,NULL);



    // Create a new model part // TODO: Use a different kind of element for each submodelpart (in order to be able of remeshing more than one kind o element or condition)
    std::unordered_map<IndexType, IndexVectorType> color_nodes, first_color_cond, second_color_cond, first_color_elem, second_color_elem;

    // The tempotal store of
    ConditionsArrayType created_conditions_vector;
    ElementsArrayType created_elements_vector;

    // Auxiliar values
    int ref, is_required;

    /* NODES */ // TODO: ADD OMP
    for (IndexType i_node = 1; i_node <= rPMMGMeshInfo.NumberOfNodes; ++i_node) {
        int id_to_write = mLocalToGlobal[i_node];
        int partition_index = local_to_partition_index[i_node];
        // int id_to_write = i_node;

        NodeType::Pointer p_node = CreateNode(rModelPart, id_to_write, ref, is_required);
        if (id_to_write==2972) {
            std::cout << "rank id partition: "<<rank << " " <<p_node->Id() << " " <<  partition_index  <<std::endl;

            KRATOS_WATCH(partition_index)
        }

        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = partition_index;

        KRATOS_ERROR_IF(partition_index<0) << "PARTITION INDEX NEGATIVE: " << partition_index << "FOR NODE: " << p_node->Id() << std::endl;
        KRATOS_ERROR_IF(partition_index>size) << "PARTITION GREATER THAN SIZE: " << partition_index << "FOR NODE: " << p_node->Id() << std::endl;

        // Set the DOFs in the nodes
        for (auto it_dof = rDofs.begin(); it_dof != rDofs.end(); ++it_dof)
            p_node->pAddDof(**it_dof);

        if (ref != 0) color_nodes[static_cast<IndexType>(ref)].push_back(id_to_write);// NOTE: ref == 0 is the MainModelPart
    }

    /* CONDITIONS */ // TODO: ADD OMP
    if (rMapPointersRefCondition.size() > 0) {
        IndexType cond_id = 1;

        IndexType counter_first_cond = 0;
        const IndexVectorType first_condition_to_remove = CheckFirstTypeConditions();
        for (IndexType i_cond = 1; i_cond <= rPMMGMeshInfo.NumberFirstTypeConditions(); ++i_cond) {
            bool skip_creation = false;
            if (counter_first_cond < first_condition_to_remove.size()) {
                if (first_condition_to_remove[counter_first_cond] == i_cond) {
                    skip_creation = true;
                    counter_first_cond += 1;
                }
            }

            Condition::Pointer p_condition = CreateFirstTypeCondition(rModelPart, rMapPointersRefCondition, cond_id, ref, is_required, skip_creation);

            if (p_condition.get() != nullptr) {
                created_conditions_vector.push_back(p_condition);
//                 rModelPart.AddCondition(p_condition);
                if (ref != 0) first_color_cond[static_cast<IndexType>(ref)].push_back(cond_id);// NOTE: ref == 0 is the MainModelPart
                cond_id += 1;
            }
        }

        IndexType counter_second_cond = 0;
        const IndexVectorType second_condition_to_remove = CheckSecondTypeConditions();
        for (IndexType i_cond = 1; i_cond <= rPMMGMeshInfo.NumberSecondTypeConditions(); ++i_cond) {
            bool skip_creation = false;
            if (counter_second_cond < second_condition_to_remove.size()) {
                if (second_condition_to_remove[counter_second_cond] == i_cond) {
                    skip_creation = true;
                    counter_second_cond += 1;
                }
            }
            Condition::Pointer p_condition = CreateSecondTypeCondition(rModelPart, rMapPointersRefCondition, cond_id, ref, is_required, skip_creation);

            if (p_condition.get() != nullptr) {
                created_conditions_vector.push_back(p_condition);
//                 rModelPart.AddCondition(p_condition);
                if (ref != 0) second_color_cond[static_cast<IndexType>(ref)].push_back(cond_id);// NOTE: ref == 0 is the MainModelPart
                cond_id += 1;
            }
        }
    }

    /* ELEMENTS */ // TODO: ADD OMP
    if (rMapPointersRefElement.size() > 0) {
        IndexType elem_id = 1;

        IndexType counter_first_elem = 0;
        const IndexVectorType first_elements_to_remove = CheckFirstTypeElements();
        for (IndexType i_elem = 1; i_elem <= rPMMGMeshInfo.NumberFirstTypeElements(); ++i_elem) {
            bool skip_creation = false;
            if (counter_first_elem < first_elements_to_remove.size()) {
                if (first_elements_to_remove[counter_first_elem] == i_elem) {
                    skip_creation = true;
                    counter_first_elem += 1;
                }
            }

            Element::Pointer p_element = CreateFirstTypeElement(rModelPart, rMapPointersRefElement, elem_id, ref, is_required, skip_creation);

            if (p_element.get() != nullptr) {
                created_elements_vector.push_back(p_element);
//                 rModelPart.AddElement(p_element);
                if (ref != 0) first_color_elem[static_cast<IndexType>(ref)].push_back(elem_id);// NOTE: ref == 0 is the MainModelPart
                elem_id += 1;
            }
        }

        IndexType counter_second_elem = 0;
        const IndexVectorType second_elements_to_remove = CheckSecondTypeElements();
        for (IndexType i_elem = 1; i_elem <= rPMMGMeshInfo.NumberSecondTypeElements(); ++i_elem) {
            bool skip_creation = false;
            if (counter_second_elem < second_elements_to_remove.size()) {
                if (second_elements_to_remove[counter_second_elem] == i_elem) {
                    skip_creation = true;
                    counter_second_elem += 1;
                }
            }

            Element::Pointer p_element = CreateSecondTypeElement(rModelPart, rMapPointersRefElement, elem_id, ref, is_required,skip_creation);

            if (p_element.get() != nullptr) {
                created_elements_vector.push_back(p_element);
//                 rModelPart.AddElement(p_element);
                if (ref != 0) second_color_elem[static_cast<IndexType>(ref)].push_back(elem_id);// NOTE: ref == 0 is the MainModelPart
                elem_id += 1;
            }
        }
    }

    // Finally we add the conditions and elements to the main model part
    rModelPart.AddConditions(created_conditions_vector.begin(), created_conditions_vector.end());
    rModelPart.AddElements(created_elements_vector.begin(), created_elements_vector.end());

    // We add nodes, conditions and elements to the sub model parts
    for (auto& r_color_list : rColors) {
        const IndexType key = r_color_list.first;

        if (key != 0) {// NOTE: key == 0 is the MainModelPart
            for (auto sub_model_part_name : r_color_list.second) {
                ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(rModelPart, sub_model_part_name);

                if (color_nodes.find(key) != color_nodes.end()) r_sub_model_part.AddNodes(color_nodes[key]);
                if (first_color_cond.find(key) != first_color_cond.end()) r_sub_model_part.AddConditions(first_color_cond[key]);
                if (second_color_cond.find(key) != second_color_cond.end()) r_sub_model_part.AddConditions(second_color_cond[key]);
                if (first_color_elem.find(key) != first_color_elem.end()) r_sub_model_part.AddElements(first_color_elem[key]);
                if (second_color_elem.find(key) != second_color_elem.end()) r_sub_model_part.AddElements(second_color_elem[key]);
            }
        } else if (mDiscretization == DiscretizationOption::ISOSURFACE) {
            if (rModelPart.HasSubModelPart("AUXILIAR_ISOSURFACE_MODEL_PART")) {
                auto& r_sub_model_part = rModelPart.GetSubModelPart("AUXILIAR_ISOSURFACE_MODEL_PART");
                r_sub_model_part.AddConditions(first_color_cond[0]);
                r_sub_model_part.AddConditions(second_color_cond[0]);
            } else {
                if (first_color_cond[0].size() + second_color_cond[0].size() > 0) {
                auto& r_sub_model_part = rModelPart.CreateSubModelPart("AUXILIAR_ISOSURFACE_MODEL_PART");
                    r_sub_model_part.AddConditions(first_color_cond[0]);
                    r_sub_model_part.AddConditions(second_color_cond[0]);
                }
            }
        }
    }

    // In case of need to remove regions we remove the unused elements
    if (mRemoveRegions) {
        rModelPart.RemoveElementsFromAllLevels(TO_ERASE);
    }

    // TODO: Add OMP
    // NOTE: We add the nodes from the elements and conditions to the respective submodelparts
    const std::vector<std::string> sub_model_part_names = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(rModelPart);

    for (auto sub_model_part_name : sub_model_part_names) {
        ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(rModelPart, sub_model_part_name);

        std::unordered_set<IndexType> node_ids;

        auto& r_sub_conditions_array = r_sub_model_part.Conditions();
        const SizeType sub_num_conditions = r_sub_conditions_array.end() - r_sub_conditions_array.begin();

        for(IndexType i = 0; i < sub_num_conditions; ++i)  {
            auto it_cond = r_sub_conditions_array.begin() + i;
            auto& r_cond_geom = it_cond->GetGeometry();

            for (SizeType i_node = 0; i_node < r_cond_geom.size(); ++i_node)
                node_ids.insert(r_cond_geom[i_node].Id());
        }

        auto& r_sub_elements_array = r_sub_model_part.Elements();
        const SizeType sub_num_elements = r_sub_elements_array.end() - r_sub_elements_array.begin();

        for(IndexType i = 0; i < sub_num_elements; ++i) {
            auto it_elem = r_sub_elements_array.begin() + i;
            auto& r_elem_geom = it_elem->GetGeometry();

            for (SizeType i_node = 0; i_node < r_elem_geom.size(); ++i_node)
                node_ids.insert(r_elem_geom[i_node].Id());
        }

        IndexVectorType vector_ids;
        std::copy(node_ids.begin(), node_ids.end(), std::back_inserter(vector_ids));
        r_sub_model_part.AddNodes(vector_ids);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::WriteSolDataToModelPart(ModelPart& rModelPart)
{
    // Iterate in the nodes
    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    const Variable<TensorArrayType>& r_tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_" + std::to_string(Dimension)+"D");

    // Auxilia metric
    TensorArrayType metric;

    // #pragma omp parallel for firstprivate(metric)
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        // We get the metric
        GetMetricTensor(metric);

        it_node->SetValue(r_tensor_variable, metric);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::WriteReferenceEntitities(
    ModelPart& rModelPart,
    const std::string& rFilename,
    std::unordered_map<IndexType,Condition::Pointer>& rRefCondition,
    std::unordered_map<IndexType,Element::Pointer>& rRefElement
    )
{
    // Getting auxiliar properties
    auto p_auxiliar_prop = rModelPart.CreateNewProperties(0);

    /* Elements */
    std::ifstream elem_infile(rFilename + ".elem.ref.json");
    KRATOS_ERROR_IF_NOT(elem_infile.good()) << "References elements file: " << rFilename  + ".json" << " cannot be found" << std::endl;
    std::stringstream elem_buffer;
    elem_buffer << elem_infile.rdbuf();
    Parameters elem_ref_json(elem_buffer.str());
    for (auto it_param = elem_ref_json.begin(); it_param != elem_ref_json.end(); ++it_param) {
        const std::size_t key = std::stoi(it_param.name());;
        Element const& r_clone_element = KratosComponents<Element>::Get(it_param->GetString());
        rRefElement[key] = r_clone_element.Create(0, r_clone_element.pGetGeometry(), p_auxiliar_prop);
    }

    /* Conditions */
    std::ifstream cond_infile(rFilename + ".cond.ref.json");
    KRATOS_ERROR_IF_NOT(cond_infile.good()) << "References conditions file: " << rFilename  + ".json" << " cannot be found" << std::endl;
    std::stringstream cond_buffer;
    cond_buffer << cond_infile.rdbuf();
    Parameters cond_ref_json(cond_buffer.str());
    for (auto it_param = cond_ref_json.begin(); it_param != cond_ref_json.end(); ++it_param) {
        const std::size_t key = std::stoi(it_param.name());;
        Condition const& r_clone_element = KratosComponents<Condition>::Get(it_param->GetString());
        rRefCondition[key] = r_clone_element.Create(0, r_clone_element.pGetGeometry(), p_auxiliar_prop);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::CreateAuxiliarSubModelPartForFlags(ModelPart& rModelPart)
{
    ModelPart& r_auxiliar_model_part = rModelPart.CreateSubModelPart("AUXILIAR_MODEL_PART_TO_LATER_REMOVE");

    const auto& r_flags = KratosComponents<Flags>::GetComponents();

    for (auto& r_flag : r_flags) {
        const std::string name_sub_model = "FLAG_" + r_flag.first;
        if (name_sub_model.find("NOT") == std::string::npos && name_sub_model.find("ALL") == std::string::npos) { // Avoiding inactive flags
            r_auxiliar_model_part.CreateSubModelPart(name_sub_model);
            ModelPart& r_auxiliar_sub_model_part = r_auxiliar_model_part.GetSubModelPart(name_sub_model);
            FastTransferBetweenModelPartsProcess(r_auxiliar_sub_model_part, rModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, *(r_flag.second)).Execute();
            // If the number of elements transfered is 0 we remove the model part
            if (r_auxiliar_sub_model_part.NumberOfNodes() == 0
            && r_auxiliar_sub_model_part.NumberOfElements() == 0
            && r_auxiliar_sub_model_part.NumberOfConditions() == 0) {
                r_auxiliar_model_part.RemoveSubModelPart(name_sub_model);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::AssignAndClearAuxiliarSubModelPartForFlags(ModelPart& rModelPart)
{
    const auto& r_flags = KratosComponents<Flags>::GetComponents();

    ModelPart& r_auxiliar_model_part = rModelPart.GetSubModelPart("AUXILIAR_MODEL_PART_TO_LATER_REMOVE");
    for (auto& r_flag : r_flags) {
        const std::string name_sub_model = "FLAG_" + r_flag.first;
        if (r_auxiliar_model_part.HasSubModelPart(name_sub_model)) {
            ModelPart& r_auxiliar_sub_model_part = r_auxiliar_model_part.GetSubModelPart(name_sub_model);
            VariableUtils().SetFlag(*(r_flag.second), true, r_auxiliar_sub_model_part.Nodes());
            VariableUtils().SetFlag(*(r_flag.second), true, r_auxiliar_sub_model_part.Conditions());
            VariableUtils().SetFlag(*(r_flag.second), true, r_auxiliar_sub_model_part.Elements());
        }
    }

    rModelPart.RemoveSubModelPart("AUXILIAR_MODEL_PART_TO_LATER_REMOVE");
}

/***********************************************************************************/
/***********************************************************************************/

template class ParMmgUtilities<PMMGLibrary::PMMG3D>;

}// namespace Kratos.
