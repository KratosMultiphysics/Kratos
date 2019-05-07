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

// External includes
/* The includes related with the MMG library */
#if !defined(MMG_INCLUDES)
#define MMG_INCLUDES
#include "mmg/libmmg.h"
#include "mmg/mmg2d/libmmg2d.h"
#include "mmg/mmg3d/libmmg3d.h"
#include "mmg/mmgs/libmmgs.h"
#endif /* MMG_INCLUDES defined */

// Project includes
#include "custom_utilities/mmg_utilities.h"

// NOTE: The following contains the license of the MMG library
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
SizeType MMGMeshInfo<MMGLibrary::MMG2D>::NumberFirstTypeConditions()
{
    return NumberOfLines;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMG3D>::NumberFirstTypeConditions()
{
    return NumberOfTriangles;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMGS>::NumberFirstTypeConditions()
{
    return NumberOfLines;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMG2D>::NumberSecondTypeConditions()
{
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMG3D>::NumberSecondTypeConditions()
{
    return NumberOfQuadrilaterals;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMGS>::NumberSecondTypeConditions()
{
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMG2D>::NumberFirstTypeElements()
{
    return NumberOfTriangles;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMG3D>::NumberFirstTypeElements()
{
    return NumberOfTetrahedra;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMGS>::NumberFirstTypeElements()
{
    return NumberOfTriangles;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMG2D>::NumberSecondTypeElements()
{
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMG3D>::NumberSecondTypeElements()
{
    return NumberOfPrism;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType MMGMeshInfo<MMGLibrary::MMGS>::NumberSecondTypeElements()
{
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template struct MMGMeshInfo<MMGLibrary::MMG2D>;
template struct MMGMeshInfo<MMGLibrary::MMG3D>;
template struct MMGMeshInfo<MMGLibrary::MMGS>;

/***********************************************************************************/
/***********************************************************************************/

// The member variables related with the MMG library
MMG5_pMesh mMmgMesh; /// The mesh data from MMG
MMG5_pSol  mMmgSol;  /// The metric variable for MMG
// MMG5_pSol  mmgDisp; /// The displacement variable for MMG

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgUtilities<TMMGLibrary>::SetEchoLevel(const SizeType EchoLevel)
{
    mEchoLevel = EchoLevel;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
SizeType MmgUtilities<TMMGLibrary>::GetEchoLevel()
{
    return mEchoLevel;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgUtilities<TMMGLibrary>::SetDiscretization(const DiscretizationOption Discretization)
{
    mDiscretization = Discretization;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
DiscretizationOption MmgUtilities<TMMGLibrary>::GetDiscretization()
{
    return mDiscretization;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgUtilities<TMMGLibrary>::SetRemoveRegions(const SizeType RemoveRegions)
{
    mRemoveRegions = RemoveRegions;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
SizeType MmgUtilities<TMMGLibrary>::GetRemoveRegions()
{
    return mRemoveRegions;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgUtilities<TMMGLibrary>::PrintAndGetMmgMeshInfo(MMGMeshInfo<TMMGLibrary>& rMMGMeshInfo)
{
    rMMGMeshInfo.NumberOfNodes = mMmgMesh->np;
    if (TMMGLibrary == MMGLibrary::MMG2D) { // 2D
        rMMGMeshInfo.NumberOfLines = mMmgMesh->na;
    } else if (TMMGLibrary == MMGLibrary::MMG3D) { // 3D
        rMMGMeshInfo.NumberOfTriangles = mMmgMesh->nt;
        rMMGMeshInfo.NumberOfQuadrilaterals = mMmgMesh->nquad;
    } else { // Surfaces
        rMMGMeshInfo.NumberOfLines = mMmgMesh->na;
    }
    if (TMMGLibrary == MMGLibrary::MMG2D) { // 2D
        rMMGMeshInfo.NumberOfTriangles = mMmgMesh->nt;
    } else if (TMMGLibrary == MMGLibrary::MMG3D) { // 3D
        rMMGMeshInfo.NumberOfTetrahedra = mMmgMesh->ne;
        rMMGMeshInfo.NumberOfPrism = mMmgMesh->nprism;
    } else { // Surfaces
        rMMGMeshInfo.NumberOfTriangles = mMmgMesh->nt;
    }

    KRATOS_INFO_IF("MmgUtilities", mEchoLevel > 0) << "\tNodes created: " << rMMGMeshInfo.NumberOfNodes << std::endl;
    if (TMMGLibrary == MMGLibrary::MMG2D) { // 2D
        KRATOS_INFO_IF("MmgUtilities", mEchoLevel > 0) <<
        "Conditions created: " << rMMGMeshInfo.NumberOfLines << "\n" <<
        "Elements created: " << rMMGMeshInfo.NumberOfTriangles << std::endl;
    } else if (TMMGLibrary == MMGLibrary::MMG3D) { // 3D
        KRATOS_INFO_IF("MmgUtilities", mEchoLevel > 0) <<
        "Conditions created: " << rMMGMeshInfo.NumberOfTriangles + rMMGMeshInfo.NumberOfQuadrilaterals << "\n\tTriangles: " << rMMGMeshInfo.NumberOfTriangles << "\tQuadrilaterals: " << rMMGMeshInfo.NumberOfQuadrilaterals << "\n" <<
        "Elements created: " << rMMGMeshInfo.NumberOfTetrahedra + rMMGMeshInfo.NumberOfPrism << "\n\tTetrahedron: " << rMMGMeshInfo.NumberOfTetrahedra << "\tPrisms: " << rMMGMeshInfo.NumberOfPrism << std::endl;
    } else { // Surfaces
        KRATOS_INFO_IF("MmgUtilities", mEchoLevel > 0) <<
        "Conditions created: " << rMMGMeshInfo.NumberOfLines << "\n" <<
        "Elements created: " << rMMGMeshInfo.NumberOfTriangles << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
IndexVectorType MmgUtilities<TMMGLibrary>::FindDuplicateNodeIds(ModelPart& rModelPart)
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
            KRATOS_WARNING_IF("MmgUtilities", mEchoLevel > 0) << "The mode " << it_node->Id() <<  " is repeated"<< std::endl;
        }
    }

    return nodes_to_remove_ids;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgUtilities<MMGLibrary::MMG2D>::CheckFirstTypeConditions()
{
    IndexVectorMapType edge_map;

    IndexVectorType ids(2);

    IndexVectorType conditions_to_remove;

    // Iterate in the conditions
    for(int i = 0; i < mMmgMesh->na; ++i) {
        int edge_0, edge_1, prop_id, is_ridge, is_required;

        KRATOS_ERROR_IF(MMG2D_Get_edge(mMmgMesh, &edge_0, &edge_1, &prop_id, &is_ridge, &is_required) != 1 ) << "Unable to get edge" << std::endl;

        ids[0] = edge_0;
        ids[1] = edge_1;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        auto& r_count = edge_map[ids];
        r_count += 1;

        if (r_count > 1)
            conditions_to_remove.push_back(i + 1);
    }

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgUtilities<MMGLibrary::MMG3D>::CheckFirstTypeConditions()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_triangles(3);

    IndexVectorType conditions_to_remove;

    for(int i = 0; i < mMmgMesh->nt; ++i) {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;

        KRATOS_ERROR_IF(MMG3D_Get_triangle(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 ) << "Unable to get triangle" << std::endl;

        ids_triangles[0] = vertex_0;
        ids_triangles[1] = vertex_1;
        ids_triangles[2] = vertex_2;

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
IndexVectorType MmgUtilities<MMGLibrary::MMGS>::CheckFirstTypeConditions()
{
    IndexVectorMapType edge_map;

    IndexVectorType ids(2);

    IndexVectorType conditions_to_remove;

    // Iterate in the conditions
    for(int i = 0; i < mMmgMesh->na; ++i) {
        int edge_0, edge_1, prop_id, is_ridge, is_required;

        KRATOS_ERROR_IF(MMGS_Get_edge(mMmgMesh, &edge_0, &edge_1, &prop_id, &is_ridge, &is_required) != 1 ) << "Unable to get edge" << std::endl;

        ids[0] = edge_0;
        ids[1] = edge_1;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        auto& r_count = edge_map[ids];
        r_count += 1;

        if (r_count > 1)
            conditions_to_remove.push_back(i + 1);
    }

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgUtilities<MMGLibrary::MMG2D>::CheckSecondTypeConditions()
{
    IndexVectorType conditions_to_remove(0);

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgUtilities<MMGLibrary::MMG3D>::CheckSecondTypeConditions()
{
    IndexVectorMapType quadrilateral_map;

    IndexVectorType ids_quadrialteral(4);

    IndexVectorType conditions_to_remove;

    for(int i = 0; i < mMmgMesh->nquad; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, prop_id, is_required;

        KRATOS_ERROR_IF(MMG3D_Get_quadrilateral(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &prop_id, &is_required) != 1 ) << "Unable to get quadrilateral" << std::endl;

        ids_quadrialteral[0] = vertex_0;
        ids_quadrialteral[1] = vertex_1;
        ids_quadrialteral[2] = vertex_2;
        ids_quadrialteral[3] = vertex_3;

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
IndexVectorType MmgUtilities<MMGLibrary::MMGS>::CheckSecondTypeConditions()
{
    IndexVectorType conditions_to_remove(0);

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgUtilities<MMGLibrary::MMG2D>::CheckFirstTypeElements()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_triangles(3);

    IndexVectorType elements_to_remove;

    // Iterate in the elements
    for(int i = 0; i < mMmgMesh->nt; ++i) {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;

        KRATOS_ERROR_IF(MMG2D_Get_triangle(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 ) << "Unable to get triangle" << std::endl;

        ids_triangles[0] = vertex_0;
        ids_triangles[1] = vertex_1;
        ids_triangles[2] = vertex_2;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_triangles.begin(), ids_triangles.end());

        auto& r_count = triangle_map[ids_triangles];
        r_count += 1;

        if (r_count > 1)
            elements_to_remove.push_back(i + 1);
    }

    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgUtilities<MMGLibrary::MMG3D>::CheckFirstTypeElements()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_tetrahedron(4);

    IndexVectorType elements_to_remove;

    for(int i = 0; i < mMmgMesh->ne; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, prop_id, is_required;

        KRATOS_ERROR_IF(MMG3D_Get_tetrahedron(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &prop_id, &is_required) != 1 ) << "Unable to get tetrahedron" << std::endl;

        ids_tetrahedron[0] = vertex_0;
        ids_tetrahedron[1] = vertex_1;
        ids_tetrahedron[2] = vertex_2;
        ids_tetrahedron[3] = vertex_3;

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
IndexVectorType MmgUtilities<MMGLibrary::MMGS>::CheckFirstTypeElements()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_triangles(3);

    IndexVectorType elements_to_remove;

    // Iterate in the elements
    for(int i = 0; i < mMmgMesh->nt; ++i) {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;

        KRATOS_ERROR_IF(MMGS_Get_triangle(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 ) << "Unable to get triangle" << std::endl;

        ids_triangles[0] = vertex_0;
        ids_triangles[1] = vertex_1;
        ids_triangles[2] = vertex_2;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_triangles.begin(), ids_triangles.end());

        auto& r_count = triangle_map[ids_triangles];
        r_count += 1;

        if (r_count > 1)
            elements_to_remove.push_back(i + 1);
    }

    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgUtilities<MMGLibrary::MMG2D>::CheckSecondTypeElements()
{
    IndexVectorType elements_to_remove(0);
    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgUtilities<MMGLibrary::MMG3D>::CheckSecondTypeElements()
{
    IndexVectorMapType prism_map;

    IndexVectorType ids_prisms(6);

    IndexVectorType elements_to_remove;

    for(int i = 0; i < mMmgMesh->nprism; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, prop_id, is_required;

        KRATOS_ERROR_IF(MMG3D_Get_prism(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &vertex_4, &vertex_5, &prop_id, &is_required) != 1 ) << "Unable to get prism" << std::endl;

        ids_prisms[0] = vertex_0;
        ids_prisms[1] = vertex_1;
        ids_prisms[2] = vertex_2;
        ids_prisms[3] = vertex_3;
        ids_prisms[4] = vertex_4;
        ids_prisms[5] = vertex_5;

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
IndexVectorType MmgUtilities<MMGLibrary::MMGS>::CheckSecondTypeElements()
{
    IndexVectorType elements_to_remove(0);
    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::BlockNode(const IndexType iNode)
{
    KRATOS_ERROR_IF(MMG2D_Set_requiredVertex(mMmgMesh, iNode) != 1 ) << "Unable to block vertex" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/


template<>
void MmgUtilities<MMGLibrary::MMG3D>::BlockNode(const IndexType iNode)
{
    KRATOS_ERROR_IF(MMG3D_Set_requiredVertex(mMmgMesh, iNode) != 1 ) << "Unable to block vertex" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::BlockNode(const IndexType iNode)
{
    KRATOS_ERROR_IF(MMGS_Set_requiredVertex(mMmgMesh, iNode) != 1 ) << "Unable to block vertex" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::BlockCondition(const IndexType iCondition)
{
    KRATOS_ERROR_IF(MMG2D_Set_requiredEdge(mMmgMesh, iCondition) != 1 ) << "Unable to block edge" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::BlockCondition(const IndexType iCondition)
{
    KRATOS_ERROR_IF(MMG3D_Set_requiredTriangle(mMmgMesh, iCondition) != 1 ) << "Unable to block edge" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::BlockCondition(const IndexType iCondition)
{
    KRATOS_ERROR_IF(MMGS_Set_requiredEdge(mMmgMesh, iCondition) != 1 ) << "Unable to block edge" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::BlockElement(const IndexType iElement)
{
    KRATOS_ERROR_IF(MMG2D_Set_requiredTriangle(mMmgMesh, iElement) != 1 ) << "Unable to block triangle" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/


template<>
void MmgUtilities<MMGLibrary::MMG3D>::BlockElement(const IndexType iElement)
{
    KRATOS_ERROR_IF(MMG3D_Set_requiredTetrahedron(mMmgMesh, iElement) != 1 ) << "Unable to block tetrahedron" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::BlockElement(const IndexType iElement)
{
    KRATOS_ERROR_IF(MMGS_Set_requiredTriangle(mMmgMesh, iElement) != 1 ) << "Unable to block triangle" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
NodeType::Pointer MmgUtilities<MMGLibrary::MMG2D>::CreateNode(
    ModelPart& rModelPart,
    const IndexType iNode,
    int& Ref,
    int& IsRequired
    )
{
    double coord_0, coord_1;
    int is_corner;

    KRATOS_ERROR_IF(MMG2D_Get_vertex(mMmgMesh, &coord_0, &coord_1, &Ref, &is_corner, &IsRequired) != 1 ) << "Unable to get vertex" << std::endl;

    NodeType::Pointer p_node = rModelPart.CreateNewNode(iNode, coord_0, coord_1, 0.0);

    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
NodeType::Pointer MmgUtilities<MMGLibrary::MMG3D>::CreateNode(
    ModelPart& rModelPart,
    const IndexType iNode,
    int& Ref,
    int& IsRequired
    )
{
    double coord_0, coord_1, coord_2;
    int is_corner;

    KRATOS_ERROR_IF(MMG3D_Get_vertex(mMmgMesh, &coord_0, &coord_1, &coord_2, &Ref, &is_corner, &IsRequired) != 1 ) << "Unable to get vertex" << std::endl;

    NodeType::Pointer p_node = rModelPart.CreateNewNode(iNode, coord_0, coord_1, coord_2);

    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
NodeType::Pointer MmgUtilities<MMGLibrary::MMGS>::CreateNode(
    ModelPart& rModelPart,
    const IndexType iNode,
    int& Ref,
    int& IsRequired
    )
{
    double coord_0, coord_1, coord_2;
    int is_corner;

    KRATOS_ERROR_IF(MMGS_Get_vertex(mMmgMesh, &coord_0, &coord_1, &coord_2, &Ref, &is_corner, &IsRequired) != 1 ) << "Unable to get vertex" << std::endl;

    NodeType::Pointer p_node = rModelPart.CreateNewNode(iNode, coord_0, coord_1, coord_2);

    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Condition::Pointer MmgUtilities<MMGLibrary::MMG2D>::CreateFirstTypeCondition(
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

    int edge_0, edge_1, is_ridge;

    KRATOS_ERROR_IF(MMG2D_Get_edge(mMmgMesh, &edge_0, &edge_1, &Ref, &is_ridge, &IsRequired) != 1 ) << "Unable to get edge" << std::endl;

    // Sometimes MMG creates conditions where there are not, then we skip
    Properties::Pointer p_prop = nullptr;
    Condition::Pointer p_base_condition = nullptr;
    if (rMapPointersRefCondition[Ref].get() == nullptr) {
        if (mDiscretization != DiscretizationOption::ISOSURFACE) { // The ISOSURFACE method creates new conditions from scratch, so we allow no previous Properties
            KRATOS_WARNING("MmgUtilities") << "Condition. Null pointer returned" << std::endl;
            return p_condition;
        } else {
            p_prop = rModelPart.pGetProperties(0);
            PointerVector<NodeType> dummy_nodes (2);
            p_base_condition = KratosComponents<Condition>::Get("LineCondition2D2N").Create(0, dummy_nodes, p_prop);
            p_base_condition->Set(MARKER);
        }
    } else {
        p_base_condition = rMapPointersRefCondition[Ref];
        p_prop = p_base_condition->pGetProperties();
    }

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (edge_0 == 0) SkipCreation = true;
    if (edge_1 == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> condition_nodes (2);
        condition_nodes[0] = rModelPart.pGetNode(edge_0);
        condition_nodes[1] = rModelPart.pGetNode(edge_1);

        p_condition = p_base_condition->Create(CondId, PointerVector<NodeType>{condition_nodes}, p_prop);
        if (p_base_condition->Is(MARKER)) p_condition->Set(MARKER);
    } else if (mEchoLevel > 2)
        KRATOS_INFO("MmgUtilities") << "Condition creation avoided" << std::endl;

    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Condition::Pointer MmgUtilities<MMGLibrary::MMG3D>::CreateFirstTypeCondition(
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

    KRATOS_ERROR_IF(MMG3D_Get_triangle(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &Ref, &IsRequired) != 1 ) << "Unable to get triangle" << std::endl;

    // Sometimes MMG creates conditions where there are not, then we skip
    Properties::Pointer p_prop = nullptr;
    Condition::Pointer p_base_condition = nullptr;

    if (rMapPointersRefCondition[Ref].get() == nullptr) {
        if (mDiscretization != DiscretizationOption::ISOSURFACE) { // The ISOSURFACE method creates new conditions from scratch, so we allow no previous Properties
            KRATOS_WARNING("MmgUtilities") << "Condition. Null pointer returned" << std::endl;
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

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> condition_nodes (3);
        condition_nodes[0] = rModelPart.pGetNode(vertex_0);
        condition_nodes[1] = rModelPart.pGetNode(vertex_1);
        condition_nodes[2] = rModelPart.pGetNode(vertex_2);

        p_condition = p_base_condition->Create(CondId, PointerVector<NodeType>{condition_nodes}, p_prop);
        if (p_base_condition->Is(MARKER)) p_condition->Set(MARKER);
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgUtilities") << "Condition creation avoided" << std::endl;

    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Condition::Pointer MmgUtilities<MMGLibrary::MMGS>::CreateFirstTypeCondition(
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

    int edge_0, edge_1, is_ridge;

    KRATOS_ERROR_IF(MMGS_Get_edge(mMmgMesh, &edge_0, &edge_1, &Ref, &is_ridge, &IsRequired) != 1 ) << "Unable to get edge" << std::endl;

    // Sometimes MMG creates conditions where there are not, then we skip
    if (rMapPointersRefCondition[Ref].get() == nullptr) {
        KRATOS_WARNING("MmgUtilities") << "Condition. Null pointer returned" << std::endl;
        return p_condition;
    }

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (edge_0 == 0) SkipCreation = true;
    if (edge_1 == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> condition_nodes (2);
        condition_nodes[0] = rModelPart.pGetNode(edge_0);
        condition_nodes[1] = rModelPart.pGetNode(edge_1);

        p_condition = rMapPointersRefCondition[Ref]->Create(CondId, PointerVector<NodeType>{condition_nodes}, rMapPointersRefCondition[Ref]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_INFO("MmgUtilities") << "Condition creation avoided" << std::endl;

    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Condition::Pointer MmgUtilities<MMGLibrary::MMG2D>::CreateSecondTypeCondition(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
    const IndexType CondId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Condition::Pointer MmgUtilities<MMGLibrary::MMG3D>::CreateSecondTypeCondition(
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

    KRATOS_ERROR_IF(MMG3D_Get_quadrilateral(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &Ref, &IsRequired) != 1 ) << "Unable to get quadrilateral" << std::endl;

    // Sometimes MMG creates conditions where there are not, then we skip
    if (rMapPointersRefCondition[Ref].get() == nullptr) {
        KRATOS_WARNING("MmgUtilities") << "Condition. Null pointer returned" << std::endl;
        return p_condition;
    }

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    if (vertex_3 == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> condition_nodes (4);
        condition_nodes[0] = rModelPart.pGetNode(vertex_0);
        condition_nodes[1] = rModelPart.pGetNode(vertex_1);
        condition_nodes[2] = rModelPart.pGetNode(vertex_2);
        condition_nodes[3] = rModelPart.pGetNode(vertex_3);

        p_condition = rMapPointersRefCondition[Ref]->Create(CondId, PointerVector<NodeType>{condition_nodes}, rMapPointersRefCondition[Ref]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgUtilities") << "Condition creation avoided" << std::endl;

    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Condition::Pointer MmgUtilities<MMGLibrary::MMGS>::CreateSecondTypeCondition(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
    const IndexType CondId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Element::Pointer MmgUtilities<MMGLibrary::MMG2D>::CreateFirstTypeElement(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
    const IndexType ElemId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    Element::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2;

    KRATOS_ERROR_IF(MMG2D_Get_triangle(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &Ref, &IsRequired) != 1 ) << "Unable to get triangle" << std::endl;

    if (mRemoveRegions && mDiscretization == DiscretizationOption::ISOSURFACE) {
        // The existence of a _nullptr_ indicates an element that was removed. This is not an alarming indicator.
        if (rMapPointersRefElement[Ref].get() == nullptr) {
            // KRATOS_INFO("MmgUtilities") << "Element has been removed from domain. Ok." << std::endl;
            return p_element;
        } else {
            // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
            if (vertex_0 == 0) SkipCreation = true;
            if (vertex_1 == 0) SkipCreation = true;
            if (vertex_2 == 0) SkipCreation = true;
            if (!SkipCreation) {
                std::vector<NodeType::Pointer> element_nodes (3);
                element_nodes[0] = rModelPart.pGetNode(vertex_0);
                element_nodes[1] = rModelPart.pGetNode(vertex_1);
                element_nodes[2] = rModelPart.pGetNode(vertex_2);
                p_element = rMapPointersRefElement[Ref]->Create(ElemId, PointerVector<NodeType>{element_nodes}, rMapPointersRefElement[Ref]->pGetProperties());
            }
        }
    } else {
        // The existence of a _nullptr_ indicates a missing element. Two options are possible: error or replacement
        Properties::Pointer p_prop = nullptr;
        Element::Pointer p_base_element = nullptr;

        // Sometimes MMG creates elements where there are not, then we skip
        if (rMapPointersRefElement[Ref].get() == nullptr) {
            if (mDiscretization != DiscretizationOption::ISOSURFACE) { // The ISOSURFACE method creates new conditions from scratch, so we allow no previous Properties
                KRATOS_WARNING("MmgUtilities") << "Element. Null pointer returned" << std::endl;
                return p_element;
            } else {
                p_prop = rModelPart.pGetProperties(0);
                PointerVector<NodeType> dummy_nodes (3);
                p_base_element = KratosComponents<Element>::Get("Element2D3N").Create(0, dummy_nodes, p_prop);
            }
        } else {
            p_base_element = rMapPointersRefElement[Ref];
            p_prop = p_base_element->pGetProperties();
        }

        // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
        if (vertex_0 == 0) SkipCreation = true;
        if (vertex_1 == 0) SkipCreation = true;
        if (vertex_2 == 0) SkipCreation = true;

        if (!SkipCreation) {
            std::vector<NodeType::Pointer> element_nodes (3);
            element_nodes[0] = rModelPart.pGetNode(vertex_0);
            element_nodes[1] = rModelPart.pGetNode(vertex_1);
            element_nodes[2] = rModelPart.pGetNode(vertex_2);

            p_element = p_base_element->Create(ElemId, PointerVector<NodeType>{element_nodes}, p_prop);
        } else if (mEchoLevel > 2)
            KRATOS_WARNING("MmgUtilities") << "Element creation avoided" << std::endl;
    }

    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Element::Pointer MmgUtilities<MMGLibrary::MMG3D>::CreateFirstTypeElement(
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

    KRATOS_ERROR_IF(MMG3D_Get_tetrahedron(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &Ref, &IsRequired) != 1 ) << "Unable to get tetrahedron" << std::endl;

    if (mRemoveRegions && mDiscretization == DiscretizationOption::ISOSURFACE) {
        // The existence of a _nullptr_ indicates an element that was removed. This is not an alarming indicator.
        if (rMapPointersRefElement[Ref].get() == nullptr) {
            // KRATOS_INFO("MmgUtilities") << "Element has been removed from domain. Ok." << std::endl;
            return p_element;
        } else {
            if (vertex_0 == 0) SkipCreation = true;
            if (vertex_1 == 0) SkipCreation = true;
            if (vertex_2 == 0) SkipCreation = true;
            if (vertex_3 == 0) SkipCreation = true;
            if (!SkipCreation) {
                std::vector<NodeType::Pointer> element_nodes (4);
                element_nodes[0] = rModelPart.pGetNode(vertex_0);
                element_nodes[1] = rModelPart.pGetNode(vertex_1);
                element_nodes[2] = rModelPart.pGetNode(vertex_2);
                element_nodes[3] = rModelPart.pGetNode(vertex_3);
                p_element = rMapPointersRefElement[Ref]->Create(ElemId, PointerVector<NodeType>{element_nodes}, rMapPointersRefElement[Ref]->pGetProperties());
            }
        }
    } else {
        // The existence of a _nullptr_ indicates a missing element. Two options are possible: error or replacement
        Properties::Pointer p_prop = nullptr;
        Element::Pointer p_base_element = nullptr;

        // Sometimes MMG creates elements where there are not, then we skip
        if (rMapPointersRefElement[Ref].get() == nullptr) {
            if (mDiscretization != DiscretizationOption::ISOSURFACE) { // The ISOSURFACE method creates new conditions from scratch, so we allow no previous Properties
                KRATOS_WARNING("MmgUtilities") << "Element. Null pointer returned" << std::endl;
                return p_element;
            } else {
                p_prop = rModelPart.pGetProperties(0);
                PointerVector<NodeType> dummy_nodes (4);
                p_base_element = KratosComponents<Element>::Get("Element3D4N").Create(0, dummy_nodes, p_prop);
            }
        } else {
            p_base_element = rMapPointersRefElement[Ref];
            p_prop = p_base_element->pGetProperties();
        }

        // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
        if (vertex_0 == 0) SkipCreation = true;
        if (vertex_1 == 0) SkipCreation = true;
        if (vertex_2 == 0) SkipCreation = true;
        if (vertex_3 == 0) SkipCreation = true;

        if (!SkipCreation) {
            std::vector<NodeType::Pointer> element_nodes (4);
            element_nodes[0] = rModelPart.pGetNode(vertex_0);
            element_nodes[1] = rModelPart.pGetNode(vertex_1);
            element_nodes[2] = rModelPart.pGetNode(vertex_2);
            element_nodes[3] = rModelPart.pGetNode(vertex_3);

            p_element = p_base_element->Create(ElemId, PointerVector<NodeType>{element_nodes}, p_prop);
        } else if (mEchoLevel > 2)
            KRATOS_WARNING("MmgUtilities") << "Element creation avoided" << std::endl;
    }

    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Element::Pointer MmgUtilities<MMGLibrary::MMGS>::CreateFirstTypeElement(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
    const IndexType ElemId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    Element::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2;

    KRATOS_ERROR_IF(MMGS_Get_triangle(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &Ref, &IsRequired) != 1 ) << "Unable to get triangle" << std::endl;

    // Sometimes MMG creates elements where there are not, then we skip
    if (rMapPointersRefElement[Ref].get() == nullptr) {
        KRATOS_WARNING("MmgUtilities") << "Element. Null pointer returned" << std::endl;
        return p_element;
    }

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> element_nodes (3);
        element_nodes[0] = rModelPart.pGetNode(vertex_0);
        element_nodes[1] = rModelPart.pGetNode(vertex_1);
        element_nodes[2] = rModelPart.pGetNode(vertex_2);

        p_element = rMapPointersRefElement[Ref]->Create(ElemId, PointerVector<NodeType>{element_nodes}, rMapPointersRefElement[Ref]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgUtilities") << "Element creation avoided" << std::endl;

    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Element::Pointer MmgUtilities<MMGLibrary::MMG2D>::CreateSecondTypeElement(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
    const IndexType ElemId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Element::Pointer MmgUtilities<MMGLibrary::MMG3D>::CreateSecondTypeElement(
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

    KRATOS_ERROR_IF(MMG3D_Get_prism(mMmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &vertex_4, &vertex_5, &Ref, &IsRequired) != 1 ) << "Unable to get prism" << std::endl;

    // Sometimes MMG creates elements where there are not, then we skip
    if (rMapPointersRefElement[Ref].get() == nullptr) {
        KRATOS_WARNING("MmgUtilities") << "Element. Null pointer returned" << std::endl;
        return p_element;
    }

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    if (vertex_3 == 0) SkipCreation = true;
    if (vertex_4 == 0) SkipCreation = true;
    if (vertex_5 == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> element_nodes (6);
        element_nodes[0] = rModelPart.pGetNode(vertex_0);
        element_nodes[1] = rModelPart.pGetNode(vertex_1);
        element_nodes[2] = rModelPart.pGetNode(vertex_2);
        element_nodes[3] = rModelPart.pGetNode(vertex_3);
        element_nodes[4] = rModelPart.pGetNode(vertex_4);
        element_nodes[5] = rModelPart.pGetNode(vertex_5);

        p_element = rMapPointersRefElement[Ref]->Create(ElemId, PointerVector<NodeType>{element_nodes}, rMapPointersRefElement[Ref]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgUtilities") << "Element creation avoided" << std::endl;

    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Element::Pointer MmgUtilities<MMGLibrary::MMGS>::CreateSecondTypeElement(
    ModelPart& rModelPart,
    std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
    const IndexType ElemId,
    int& Ref,
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::InitMesh()
{
    mMmgMesh = nullptr;
    mMmgSol = nullptr;
//     MmgDisp = nullptr;

    // We init the MMG mesh and sol
    if (mDiscretization == DiscretizationOption::STANDARD) {
        MMG2D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mMmgMesh, MMG5_ARG_ppMet, &mMmgSol, MMG5_ARG_end);
    } else if (mDiscretization == DiscretizationOption::ISOSURFACE) {
        MMG2D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mMmgMesh, MMG5_ARG_ppLs, &mMmgSol, MMG5_ARG_end);
    } else {
        KRATOS_ERROR << "Discretization type: " << static_cast<int>(mDiscretization) << " not fully implemented" << std::endl;
    }

    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::InitMesh()
{
    mMmgMesh = nullptr;
    mMmgSol = nullptr;
//     MmgDisp = nullptr;

    // We init the MMG mesh and sol
    if (mDiscretization == DiscretizationOption::STANDARD) {
        MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mMmgMesh, MMG5_ARG_ppMet, &mMmgSol, MMG5_ARG_end);
//     } else if (Discretization == DiscretizationOption::LAGRANGIAN) {
//         MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mMmgMesh, MMG5_ARG_ppMet, &mMmgSol, MMG5_ARG_ppDisp, &MmgDisp, MMG5_ARG_end);
    } else if (mDiscretization == DiscretizationOption::ISOSURFACE) {
        MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mMmgMesh, MMG5_ARG_ppLs, &mMmgSol, MMG5_ARG_end);
    } else {
        KRATOS_ERROR << "Discretization type: " << static_cast<int>(mDiscretization) << " not fully implemented" << std::endl;
    }

    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::InitMesh()
{
    mMmgMesh = nullptr;
    mMmgSol = nullptr;
//     MmgDisp = nullptr;

    // We init the MMG mesh and sol
    if (mDiscretization == DiscretizationOption::STANDARD) {
        MMGS_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mMmgMesh, MMG5_ARG_ppMet, &mMmgSol, MMG5_ARG_end);
    } else if (mDiscretization == DiscretizationOption::ISOSURFACE) {
        MMGS_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mMmgMesh, MMG5_ARG_ppLs, &mMmgSol, MMG5_ARG_end);
    } else {
        KRATOS_ERROR << "Discretization type: " << static_cast<int>(mDiscretization) << " not fully implemented" << std::endl;
    }

    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgUtilities<TMMGLibrary>::InitVerbosity()
{
    /* We set the MMG verbosity */
    int verbosity_mmg;
    if (mEchoLevel == 0)
        verbosity_mmg = -1;
    else if (mEchoLevel == 1)
        verbosity_mmg = 0; // NOTE: This way just the essential info from MMG will be printed, but the custom message will appear
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
void MmgUtilities<MMGLibrary::MMG2D>::InitVerbosityParameter(const IndexType VerbosityMMG)
{
    KRATOS_ERROR_IF( !MMG2D_Set_iparameter(mMmgMesh,mMmgSol,MMG2D_IPARAM_verbose, VerbosityMMG) ) << "Unable to set verbosity" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::InitVerbosityParameter(const IndexType VerbosityMMG)
{
    KRATOS_ERROR_IF( !MMG3D_Set_iparameter(mMmgMesh,mMmgSol,MMG3D_IPARAM_verbose, VerbosityMMG) ) << "Unable to set verbosity" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::InitVerbosityParameter(const IndexType VerbosityMMG)
{
    KRATOS_ERROR_IF( !MMGS_Set_iparameter(mMmgMesh,mMmgSol,MMGS_IPARAM_verbose, VerbosityMMG) ) << "Unable to set verbosity" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetMeshSize(MMGMeshInfo<MMGLibrary::MMG2D>& rMMGMeshInfo)
{
    //Give the size of the mesh: NumNodes vertices, num_elements triangles, num_conditions edges (2D)
    KRATOS_ERROR_IF( MMG2D_Set_meshSize(mMmgMesh, rMMGMeshInfo.NumberOfNodes, rMMGMeshInfo.NumberOfTriangles, rMMGMeshInfo.NumberOfLines) != 1 ) << "Unable to set mesh size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetMeshSize(MMGMeshInfo<MMGLibrary::MMG3D>& rMMGMeshInfo)
{
    //Give the size of the mesh: NumNodes Vertex, num_elements tetra and prism, NumArrayConditions triangles and quadrilaterals, 0 edges (3D)
    KRATOS_ERROR_IF( MMG3D_Set_meshSize(mMmgMesh, rMMGMeshInfo.NumberOfNodes, rMMGMeshInfo.NumberOfTetrahedra, rMMGMeshInfo.NumberOfPrism, rMMGMeshInfo.NumberOfTriangles, rMMGMeshInfo.NumberOfQuadrilaterals, 0) != 1 ) << "Unable to set mesh size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetMeshSize(MMGMeshInfo<MMGLibrary::MMGS>& rMMGMeshInfo)
{
    //Give the size of the mesh: NumNodes vertices, num_elements triangles, num_conditions edges (3D)
    KRATOS_ERROR_IF( MMGS_Set_meshSize(mMmgMesh, rMMGMeshInfo.NumberOfNodes, rMMGMeshInfo.NumberOfTriangles, rMMGMeshInfo.NumberOfLines) != 1 ) << "Unable to set mesh size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetSolSizeScalar(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMG2D_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Scalar) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetSolSizeScalar(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMG3D_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Scalar) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetSolSizeScalar(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMGS_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Scalar) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetSolSizeVector(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMG2D_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Vector) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetSolSizeVector(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMG3D_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Vector) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetSolSizeVector(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMGS_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Vector) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetSolSizeTensor(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMG2D_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetSolSizeTensor(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMG3D_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetSolSizeTensor(const SizeType NumNodes)
{
    KRATOS_ERROR_IF( MMGS_Set_solSize(mMmgMesh,mMmgSol,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 ) << "Unable to set metric size" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::CheckMeshData()
{
    KRATOS_ERROR_IF( MMG2D_Chk_meshData(mMmgMesh, mMmgSol) != 1 ) << "Wrong mesh data" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::CheckMeshData()
{
    KRATOS_ERROR_IF( MMG3D_Chk_meshData(mMmgMesh, mMmgSol) != 1 ) << "Wrong mesh data" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::CheckMeshData()
{
    KRATOS_ERROR_IF( MMGS_Chk_meshData(mMmgMesh, mMmgSol) != 1 ) << "Wrong mesh data" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::InputMesh(const std::string& rOutputName)
{
    const std::string mesh_name = rOutputName + ".mesh";
    const char* mesh_file = mesh_name.c_str();

    // a)  Give the ouptut mesh name using MMG2D_Set_inputMeshName
    MMG2D_Set_inputMeshName(mMmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMG2D_loadMesh(mMmgMesh,mesh_file) != 1) << "UNABLE TO LOAD MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::InputMesh(const std::string& rOutputName)
{
    const std::string mesh_name = rOutputName + ".mesh";
    const char* mesh_file = mesh_name.c_str();

    // a)  Give the ouptut mesh name using MMG3D_Set_inputMeshName
    MMG3D_Set_inputMeshName(mMmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMG3D_loadMesh(mMmgMesh,mesh_file) != 1) << "UNABLE TO LOAD MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::InputMesh(const std::string& rOutputName)
{
    const std::string mesh_name = rOutputName + ".mesh";
    const char* mesh_file = mesh_name.c_str();

    // a)  Give the ouptut mesh name using MMGS_Set_inputMeshName
    MMGS_Set_inputMeshName(mMmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMGS_loadMesh(mMmgMesh,mesh_file) != 1) << "UNABLE TO LOAD MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::InputSol(const std::string& rInputName)
{
    const std::string sol_name = rInputName + ".sol";
    const char* sol_file = sol_name.c_str();

    // a)  Give the ouptut mesh name using MMG2D_Set_inputSolName (by default, the mesh is saved in the "mesh.o.mesh") file
    MMG2D_Set_inputSolName(mMmgMesh, mMmgSol, sol_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMG2D_loadSol(mMmgMesh, mMmgSol, sol_file) != 1) << "UNABLE TO READ METRIC" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::InputSol(const std::string& rInputName)
{
    const std::string sol_name = rInputName + ".sol";
    const char* sol_file = sol_name.c_str();

    // a)  Give the ouptut mesh name using MMG3D_Set_inputSolName (by default, the mesh is saved in the "mesh.o.mesh") file
    MMG3D_Set_inputSolName(mMmgMesh, mMmgSol, sol_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMG3D_loadSol(mMmgMesh, mMmgSol, sol_file) != 1) << "UNABLE TO READ METRIC" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::InputSol(const std::string& rInputName)
{
    const std::string sol_name = rInputName + ".sol";
    const char* sol_file = sol_name.c_str();

    // a)  Give the ouptut mesh name using MMGS_Set_inputSolName (by default, the mesh is saved in the "mesh.o.mesh") file
    MMGS_Set_inputSolName(mMmgMesh, mMmgSol, sol_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMGS_loadSol(mMmgMesh, mMmgSol, sol_file) != 1) << "UNABLE TO READ METRIC" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::OutputMesh(
    const std::string& rOutputName,
    const bool PostOutput,
    const IndexType Step
    )
{
    const std::string mesh_name = rOutputName + "_step=" + std::to_string(Step) + (PostOutput ? ".o" : "") + ".mesh";
    const char* mesh_file = mesh_name.c_str();

    // a)  Give the ouptut mesh name using MMG2D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh") file
    MMG2D_Set_outputMeshName(mMmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMG2D_saveMesh(mMmgMesh,mesh_file) != 1) << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::OutputMesh(
    const std::string& rOutputName,
    const bool PostOutput,
    const IndexType Step
    )
{
    const std::string mesh_name = rOutputName + "_step=" + std::to_string(Step) + (PostOutput ? ".o" : "") + ".mesh";
    const char* mesh_file = mesh_name.c_str();

    // a)  Give the ouptut mesh name using MMG3D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh") file
    MMG3D_Set_outputMeshName(mMmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMG3D_saveMesh(mMmgMesh,mesh_file) != 1) << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::OutputMesh(
    const std::string& rOutputName,
    const bool PostOutput,
    const IndexType Step
    )
{
    const std::string mesh_name = rOutputName + "_step=" + std::to_string(Step) + (PostOutput ? ".o" : "") + ".mesh";
    const char* mesh_file = mesh_name.c_str();

    // a)  Give the ouptut mesh name using MMGS_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh") file
    MMGS_Set_outputMeshName(mMmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgUtilities", MMGS_saveMesh(mMmgMesh,mesh_file) != 1) << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::OutputSol(
    const std::string& rOutputName,
    const bool PostOutput,
    const IndexType Step
    )
{
    const std::string sol_name = rOutputName + "_step=" + std::to_string(Step) + (PostOutput ? ".o" : "") + ".sol";
    const char* sol_file = sol_name.c_str();

    // a)  Give the ouptut sol name using MMG2D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file
    MMG2D_Set_outputSolName(mMmgMesh, mMmgSol, sol_file);

    // b) Function calling
    KRATOS_INFO_IF("MmgUtilities", MMG2D_saveSol(mMmgMesh, mMmgSol, sol_file) != 1) << "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::OutputSol(
    const std::string& rOutputName,
    const bool PostOutput,
    const IndexType Step
    )
{
    const std::string sol_name = rOutputName + "_step=" + std::to_string(Step) + (PostOutput ? ".o" : "") + ".sol";
    const char* sol_file = sol_name.c_str();

    // a)  Give the ouptut sol name using MMG3D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file
    MMG3D_Set_outputSolName(mMmgMesh, mMmgSol, sol_file);

    // b) Function calling
    KRATOS_INFO_IF("MmgUtilities", MMG3D_saveSol(mMmgMesh,mMmgSol, sol_file) != 1)<< "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::OutputSol(
    const std::string& rOutputName,
    const bool PostOutput,
    const IndexType Step
    )
{
    const std::string sol_name = rOutputName + "_step=" + std::to_string(Step) + (PostOutput ? ".o" : "") + ".sol";
    const char* sol_file = sol_name.c_str();

    // a)  Give the ouptut sol name using MMGS_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file
    MMGS_Set_outputSolName(mMmgMesh, mMmgSol, sol_file);

    // b) Function calling
    KRATOS_INFO_IF("MmgUtilities", MMGS_saveSol(mMmgMesh,mMmgSol, sol_file) != 1)<< "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::FreeAll()
{
    MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mMmgMesh,MMG5_ARG_ppMet,&mMmgSol,MMG5_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::FreeAll()
{
    MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mMmgMesh,MMG5_ARG_ppMet,&mMmgSol,MMG5_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::FreeAll()
{
    MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mMmgMesh,MMG5_ARG_ppMet,&mMmgSol,MMG5_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::MMGLibCallMetric(Parameters ConfigurationParameters)
{
    KRATOS_TRY;

    /* Advanced configurations */
    // Global hausdorff value (default value = 0.01) applied on the whole boundary
    if (ConfigurationParameters["advanced_parameters"]["force_hausdorff_value"].GetBool()) {
        if ( MMG2D_Set_dparameter(mMmgMesh,mMmgSol,MMG2D_DPARAM_hausd, ConfigurationParameters["advanced_parameters"]["hausdorff_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the Hausdorff parameter" << std::endl;
    }

    // Avoid/allow point relocation
    if ( MMG2D_Set_iparameter(mMmgMesh,mMmgSol,MMG2D_IPARAM_nomove, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_move_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to fix the nodes" << std::endl;

    // Avoid/allow surface modifications
    if ( MMG2D_Set_iparameter(mMmgMesh,mMmgSol,MMG2D_IPARAM_nosurf, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_surf_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no surfacic modifications" << std::endl;

    // Don't insert nodes on mesh
    if ( MMG2D_Set_iparameter(mMmgMesh,mMmgSol,MMG2D_IPARAM_noinsert, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_insert_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no insertion/suppression point" << std::endl;

    // Don't swap mesh
    if ( MMG2D_Set_iparameter(mMmgMesh,mMmgSol,MMG2D_IPARAM_noswap, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_swap_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // Set the angle detection
    const bool deactivate_detect_angle = ConfigurationParameters["advanced_parameters"]["deactivate_detect_angle"].GetBool();
    if ( deactivate_detect_angle) {
        if ( MMG2D_Set_iparameter(mMmgMesh,mMmgSol,MMG2D_IPARAM_angle, static_cast<int>(!deactivate_detect_angle)) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection on" << std::endl;
    }

    // Set the gradation
    if (ConfigurationParameters["advanced_parameters"]["force_gradation_value"].GetBool()) {
        if ( MMG2D_Set_dparameter(mMmgMesh,mMmgSol,MMG2D_DPARAM_hgrad, ConfigurationParameters["advanced_parameters"]["gradation_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set gradation" << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_min"].GetBool()) {
        if ( MMG2D_Set_dparameter(mMmgMesh,mMmgSol,MMG2D_DPARAM_hmin, ConfigurationParameters["force_sizes"]["minimal_size"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the minimal edge size " << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_max"].GetBool()) {
        if ( MMG2D_Set_dparameter(mMmgMesh,mMmgSol,MMG2D_DPARAM_hmax, ConfigurationParameters["force_sizes"]["maximal_size"].GetDouble()) != 1 ) {
            KRATOS_ERROR << "Unable to set the maximal edge size " << std::endl;
        }
    }

    const int ier = MMG2D_mmg2dlib(mMmgMesh, mMmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG2DLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::MMGLibCallIsoSurface(Parameters ConfigurationParameters)
{
    KRATOS_TRY;

    /**------------------- Level set discretization option ---------------------*/
    /* Ask for level set discretization */
    KRATOS_ERROR_IF( MMG2D_Set_iparameter(mMmgMesh,mMmgSol,MMG2D_IPARAM_iso, 1) != 1 ) << "Unable to ask for level set discretization" << std::endl;

    /** (Not mandatory): check if the number of given entities match with mesh size */
    KRATOS_ERROR_IF( MMG2D_Chk_meshData(mMmgMesh,mMmgSol) != 1 ) << "Unable to check if the number of given entities match with mesh size" << std::endl;

    /**------------------- Level set discretization ---------------------------*/

//     /* Debug mode ON (default value = OFF) */
//     KRATOS_ERROR_IF( MMG2D_Set_iparameter(mMmgMesh,mMmgSol,MMG2D_IPARAM_debug, 1) != 1 ) << "Unable to set on debug mode" << std::endl;

    const int ier = MMG2D_mmg2dls(mMmgMesh, mMmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG2DLS: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG2DLS. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::MMGLibCallMetric(Parameters ConfigurationParameters)
{
    KRATOS_TRY;

    /* Advanced configurations */
    // Global hausdorff value (default value = 0.01) applied on the whole boundary
    if (ConfigurationParameters["advanced_parameters"]["force_hausdorff_value"].GetBool()) {
        if ( MMG3D_Set_dparameter(mMmgMesh,mMmgSol,MMG3D_DPARAM_hausd, ConfigurationParameters["advanced_parameters"]["hausdorff_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the Hausdorff parameter" << std::endl;
    }

    // Avoid/allow point relocation
    if ( MMG3D_Set_iparameter(mMmgMesh,mMmgSol,MMG3D_IPARAM_nomove, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_move_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to fix the nodes" << std::endl;

    // Avoid/allow surface modifications
    if ( MMG3D_Set_iparameter(mMmgMesh,mMmgSol,MMG3D_IPARAM_nosurf, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_surf_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no surfacic modifications" << std::endl;

    // Don't insert nodes on mesh
    if ( MMG3D_Set_iparameter(mMmgMesh,mMmgSol,MMG3D_IPARAM_noinsert, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_insert_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no insertion/suppression point" << std::endl;

    // Don't swap mesh
    if ( MMG3D_Set_iparameter(mMmgMesh,mMmgSol,MMG3D_IPARAM_noswap, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_swap_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // Set the angle detection
    const bool deactivate_detect_angle = ConfigurationParameters["advanced_parameters"]["deactivate_detect_angle"].GetBool();
    if ( deactivate_detect_angle) {
        if ( MMG3D_Set_iparameter(mMmgMesh,mMmgSol,MMG3D_IPARAM_angle, static_cast<int>(!deactivate_detect_angle)) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection on" << std::endl;
    }

    // Set the gradation
    if (ConfigurationParameters["advanced_parameters"]["force_gradation_value"].GetBool()) {
        if ( MMG3D_Set_dparameter(mMmgMesh,mMmgSol,MMG3D_DPARAM_hgrad, ConfigurationParameters["advanced_parameters"]["gradation_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set gradation" << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_min"].GetBool()) {
        if ( MMG3D_Set_dparameter(mMmgMesh,mMmgSol,MMG3D_DPARAM_hmin, ConfigurationParameters["force_sizes"]["minimal_size"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the minimal edge size " << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_max"].GetBool()) {
        if ( MMG3D_Set_dparameter(mMmgMesh,mMmgSol,MMG3D_DPARAM_hmax, ConfigurationParameters["force_sizes"]["maximal_size"].GetDouble()) != 1 ) {
            KRATOS_ERROR << "Unable to set the maximal edge size " << std::endl;
        }
    }

    const int ier = MMG3D_mmg3dlib(mMmgMesh, mMmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG3DLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::MMGLibCallIsoSurface(Parameters ConfigurationParameters)
{
    KRATOS_TRY;

    /**------------------- Level set discretization option ---------------------*/
    /* Ask for level set discretization */
    KRATOS_ERROR_IF( MMG3D_Set_iparameter(mMmgMesh,mMmgSol,MMG3D_IPARAM_iso, 1) != 1 ) << "Unable to ask for level set discretization" << std::endl;

    /** (Not mandatory): check if the number of given entities match with mesh size */
    KRATOS_ERROR_IF( MMG3D_Chk_meshData(mMmgMesh,mMmgSol) != 1 ) << "Unable to check if the number of given entities match with mesh size" << std::endl;

    /**------------------- Including surface options ---------------------------*/

    // Global hausdorff value (default value = 0.01) applied on the whole boundary
    if (ConfigurationParameters["advanced_parameters"]["force_hausdorff_value"].GetBool()) {
        if ( MMG3D_Set_dparameter(mMmgMesh,mMmgSol,MMG3D_DPARAM_hausd, ConfigurationParameters["advanced_parameters"]["hausdorff_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the Hausdorff parameter" << std::endl;
    }

    // Set the gradation
    if (ConfigurationParameters["advanced_parameters"]["force_gradation_value"].GetBool()) {
        if ( MMG3D_Set_dparameter(mMmgMesh,mMmgSol,MMG3D_DPARAM_hgrad, ConfigurationParameters["advanced_parameters"]["gradation_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set gradation" << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_min"].GetBool()) {
        if ( MMG3D_Set_dparameter(mMmgMesh,mMmgSol,MMG3D_DPARAM_hmin, ConfigurationParameters["force_sizes"]["minimal_size"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the minimal edge size " << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_max"].GetBool()) {
        if ( MMG3D_Set_dparameter(mMmgMesh,mMmgSol,MMG3D_DPARAM_hmax, ConfigurationParameters["force_sizes"]["maximal_size"].GetDouble()) != 1 ) {
            KRATOS_ERROR << "Unable to set the maximal edge size " << std::endl;
        }
    }

    /**------------------- level set discretization ---------------------------*/
//     /* Debug mode ON (default value = OFF) */
//     KRATOS_ERROR_IF( MMG3D_Set_iparameter(mMmgMesh,mMmgSol,MMG3D_IPARAM_debug, 1) != 1 ) << "Unable to set on debug mode" << std::endl;

    const int ier = MMG3D_mmg3dls(mMmgMesh, mMmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG3DLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::MMGLibCallMetric(Parameters ConfigurationParameters)
{
    KRATOS_TRY;

    /* Advanced configurations */
    // Global hausdorff value (default value = 0.01) applied on the whole boundary
    if (ConfigurationParameters["advanced_parameters"]["force_hausdorff_value"].GetBool()) {
        if ( MMGS_Set_dparameter(mMmgMesh,mMmgSol,MMGS_DPARAM_hausd, ConfigurationParameters["advanced_parameters"]["hausdorff_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the Hausdorff parameter" << std::endl;
    }

    // Avoid/allow point relocation
    if ( MMGS_Set_iparameter(mMmgMesh,mMmgSol,MMGS_IPARAM_nomove, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_move_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to fix the nodes" << std::endl;

    // Don't insert nodes on mesh
    if ( MMGS_Set_iparameter(mMmgMesh,mMmgSol,MMGS_IPARAM_noinsert, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_insert_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no insertion/suppression point" << std::endl;

    // Don't swap mesh
    if ( MMGS_Set_iparameter(mMmgMesh,mMmgSol,MMGS_IPARAM_noswap, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_swap_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // Disabled/enabled normal regularization
    if ( MMGS_Set_iparameter(mMmgMesh,mMmgSol,MMGS_IPARAM_nreg, static_cast<int>(ConfigurationParameters["advanced_parameters"]["normal_regularization_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable disabled/enabled normal regularization " << std::endl;

    // Set the angle detection
    const bool deactivate_detect_angle = ConfigurationParameters["advanced_parameters"]["deactivate_detect_angle"].GetBool();
    if ( deactivate_detect_angle) {
        if ( MMGS_Set_iparameter(mMmgMesh,mMmgSol,MMGS_IPARAM_angle, static_cast<int>(!deactivate_detect_angle)) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection on" << std::endl;
    }

    // Set the gradation
    if (ConfigurationParameters["advanced_parameters"]["force_gradation_value"].GetBool()) {
        if ( MMGS_Set_dparameter(mMmgMesh,mMmgSol,MMGS_DPARAM_hgrad, ConfigurationParameters["advanced_parameters"]["gradation_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set gradation" << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_min"].GetBool()) {
        if ( MMGS_Set_dparameter(mMmgMesh,mMmgSol,MMGS_DPARAM_hmin, ConfigurationParameters["force_sizes"]["minimal_size"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the minimal edge size " << std::endl;
    }

    // Minimal edge size
    if (ConfigurationParameters["force_sizes"]["force_max"].GetBool()) {
        if ( MMGS_Set_dparameter(mMmgMesh,mMmgSol,MMGS_DPARAM_hmax, ConfigurationParameters["force_sizes"]["maximal_size"].GetDouble()) != 1 ) {
            KRATOS_ERROR << "Unable to set the maximal edge size " << std::endl;
        }
    }

    // Compute remesh
    const int ier = MMGS_mmgslib(mMmgMesh, mMmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMGSLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMGSLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::MMGLibCallIsoSurface(Parameters ConfigurationParameters)
{
    KRATOS_TRY;

    /**------------------- Level set discretization option ---------------------*/
    /* Ask for level set discretization */
    KRATOS_ERROR_IF( MMGS_Set_iparameter(mMmgMesh,mMmgSol,MMGS_IPARAM_iso, 1) != 1 ) << "Unable to ask for level set discretization" << std::endl;

    /** (Not mandatory): check if the number of given entities match with mesh size */
    KRATOS_ERROR_IF( MMGS_Chk_meshData(mMmgMesh,mMmgSol) != 1 ) << "Unable to check if the number of given entities match with mesh size" << std::endl;

    /**------------------- level set discretization ---------------------------*/
//     /* Debug mode ON (default value = OFF) */
//     KRATOS_ERROR_IF( MMGS_Set_iparameter(mMmgMesh,mMmgSol,MMGS_IPARAM_debug, 1) != 1 ) << "Unable to set on debug mode" << std::endl;

    const int ier = MMGS_mmgsls(mMmgMesh, mMmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMGSLS: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMGSLS. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const IndexType Color,
    const IndexType Index
    )
{
    KRATOS_ERROR_IF( MMG2D_Set_vertex(mMmgMesh, X, Y, Color, Index) != 1 ) << "Unable to set vertex" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const IndexType Color,
    const IndexType Index
    )
{
    KRATOS_ERROR_IF( MMG3D_Set_vertex(mMmgMesh, X, Y, Z, Color, Index) != 1 ) << "Unable to set vertex" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const IndexType Color,
    const IndexType Index
    )
{
    KRATOS_ERROR_IF( MMGS_Set_vertex(mMmgMesh, X, Y, Z, Color, Index) != 1 ) << "Unable to set vertex" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetConditions(
    GeometryType& rGeometry,
    const IndexType Color,
    const IndexType Index
    )
{
    if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point2D) // Point
        KRATOS_ERROR << "ERROR:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2) { // Line
        const IndexType id_1 = rGeometry[0].Id(); // First node id
        const IndexType id_2 = rGeometry[1].Id(); // Second node id

        KRATOS_ERROR_IF( MMG2D_Set_edge(mMmgMesh, id_1, id_2, Color, Index) != 1 ) << "Unable to set edge" << std::endl;

        // Set fixed boundary
        bool blocked_1 = false;
        if (rGeometry[0].IsDefined(BLOCKED))
            blocked_1 = rGeometry[0].Is(BLOCKED);
        bool blocked_2 = false;
        if (rGeometry[1].IsDefined(BLOCKED))
            blocked_2 = rGeometry[1].Is(BLOCKED);

        if (blocked_1 && blocked_2) BlockCondition(Index);
    } else {
        const IndexType size_geometry = rGeometry.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << rGeometry.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetConditions(
    GeometryType& rGeometry,
    const IndexType Color,
    const IndexType Index
    )
{
    if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point3D) // Point
        KRATOS_ERROR << "ERROR:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2) { // Line
        KRATOS_ERROR << "Kratos_Line3D2 remeshing pending to be implemented" << std::endl;
//         const IndexType id1 = rGeometry[0].Id(); // First node id
//         const IndexType id2 = rGeometry[1].Id(); // Second node id
//
//         KRATOS_ERROR_IF( MMG3D_Set_edge(mMmgMesh, id1, id2, Color, Index) != 1 ) << "Unable to set edge" << std::endl;
//
//         // Set fixed boundary
//         bool blocked_1 = false;
//         if (rGeometry[0].IsDefined(BLOCKED))
//             blocked_1 = rGeometry[0].Is(BLOCKED);
//         bool blocked_2 = false;
//         if (rGeometry[1].IsDefined(BLOCKED))
//             blocked_2 = rGeometry[1].Is(BLOCKED);
//
//         if ((blocked_1 && blocked_2)) KRATOS_ERROR_IF( MMG3D_Set_requiredEdge(mMmgMesh, Index) != 1 ) << "Unable to block edge" << std::endl;
    } else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {// Triangle
        const IndexType id_1 = rGeometry[0].Id(); // First node Id
        const IndexType id_2 = rGeometry[1].Id(); // Second node Id
        const IndexType id_3 = rGeometry[2].Id(); // Third node Id

        KRATOS_ERROR_IF( MMG3D_Set_triangle(mMmgMesh, id_1, id_2, id_3, Color, Index) != 1 ) << "Unable to set triangle" << std::endl;

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
        const IndexType id_1 = rGeometry[0].Id(); // First node Id
        const IndexType id_2 = rGeometry[1].Id(); // Second node Id
        const IndexType id_3 = rGeometry[2].Id(); // Third node Id
        const IndexType id_4 = rGeometry[3].Id(); // Fourth node Id

        KRATOS_ERROR_IF( MMG3D_Set_quadrilateral(mMmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 ) << "Unable to set quadrilateral" << std::endl;
    } else {
        const SizeType size_geometry = rGeometry.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << rGeometry.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetConditions(
    GeometryType& rGeometry,
    const IndexType Color,
    const IndexType Index
    )
{
    if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point3D) // Point
        KRATOS_ERROR << "ERROR:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2) { // Line
        const IndexType id_1 = rGeometry[0].Id(); // First node id
        const IndexType id_2 = rGeometry[1].Id(); // Second node id

        KRATOS_ERROR_IF( MMGS_Set_edge(mMmgMesh, id_1, id_2, Color, Index) != 1 ) << "Unable to set edge" << std::endl;

        // Set fixed boundary
        bool blocked_1 = false;
        if (rGeometry[0].IsDefined(BLOCKED))
            blocked_1 = rGeometry[0].Is(BLOCKED);
        bool blocked_2 = false;
        if (rGeometry[1].IsDefined(BLOCKED))
            blocked_2 = rGeometry[1].Is(BLOCKED);

        if (blocked_1 && blocked_2) BlockCondition(Index);
    } else {
        const IndexType size_geometry = rGeometry.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << rGeometry.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetElements(
    GeometryType& rGeometry,
    const IndexType Color,
    const IndexType Index
    )
{
    const IndexType id_1 = rGeometry[0].Id(); // First node Id
    const IndexType id_2 = rGeometry[1].Id(); // Second node Id
    const IndexType id_3 = rGeometry[2].Id(); // Third node Id

    KRATOS_ERROR_IF( MMG2D_Set_triangle(mMmgMesh, id_1, id_2, id_3, Color, Index) != 1 ) << "Unable to set triangle" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetElements(
    GeometryType& rGeometry,
    const IndexType Color,
    const IndexType Index
    )
{
    const IndexType id_1 = rGeometry[0].Id(); // First node Id
    const IndexType id_2 = rGeometry[1].Id(); // Second node Id
    const IndexType id_3 = rGeometry[2].Id(); // Third node Id
    const IndexType id_4 = rGeometry[3].Id(); // Fourth node Id

    if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) { // Tetrahedron
        KRATOS_ERROR_IF( MMG3D_Set_tetrahedron(mMmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 ) << "Unable to set tetrahedron" << std::endl;
    } else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) { // Prisms
        const IndexType id_5 = rGeometry[4].Id(); // 5th node Id
        const IndexType id_6 = rGeometry[5].Id(); // 6th node Id

        KRATOS_ERROR_IF( MMG3D_Set_prism(mMmgMesh, id_1, id_2, id_3, id_4, id_5, id_6, Color, Index) != 1 ) << "Unable to set prism" << std::endl;
    } else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) { // Hexaedron
//         const IndexType id_5 = rGeometry[4].Id(); // 5th node Id
//         const IndexType id_6 = rGeometry[5].Id(); // 6th node Id
//         const IndexType id_6 = rGeometry[7].Id(); // 7th node Id
//         const IndexType id_6 = rGeometry[8].Id(); // 8th node Id

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
void MmgUtilities<MMGLibrary::MMGS>::SetElements(
    GeometryType& rGeometry,
    const IndexType Color,
    const IndexType Index
    )
{
    const IndexType id_1 = rGeometry[0].Id(); // First node Id
    const IndexType id_2 = rGeometry[1].Id(); // Second node Id
    const IndexType id_3 = rGeometry[2].Id(); // Third node Id

    KRATOS_ERROR_IF( MMGS_Set_triangle(mMmgMesh, id_1, id_2, id_3, Color, Index) != 1 ) << "Unable to set triangle" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetMetricScalar(

    const double Metric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( MMG2D_Set_scalarSol(mMmgSol, Metric, NodeId) != 1 ) << "Unable to set scalar metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetMetricScalar(
    const double Metric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( MMG3D_Set_scalarSol(mMmgSol, Metric, NodeId) != 1 ) << "Unable to set scalar metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetMetricScalar(
    const double Metric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( MMGS_Set_scalarSol(mMmgSol, Metric, NodeId) != 1 ) << "Unable to set scalar metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetMetricVector(
    const array_1d<double, 2>& Metric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( MMG2D_Set_vectorSol(mMmgSol, Metric[0], Metric[1], NodeId) != 1 ) << "Unable to set vector metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetMetricVector(
    const array_1d<double, 3>& Metric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( MMG3D_Set_vectorSol(mMmgSol, Metric[0], Metric[1], Metric[2], NodeId) != 1 ) << "Unable to set vector metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetMetricVector(
    const array_1d<double, 3>& Metric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( MMGS_Set_vectorSol(mMmgSol, Metric[0], Metric[1], Metric[2], NodeId) != 1 ) << "Unable to set vector metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG2D>::SetMetricTensor(
    const array_1d<double, 3>& Metric,
    const IndexType NodeId
    )
{
    // The order is XX, XY, YY
    KRATOS_ERROR_IF( MMG2D_Set_tensorSol(mMmgSol, Metric[0], Metric[2], Metric[1], NodeId) != 1 ) << "Unable to set tensor metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMG3D>::SetMetricTensor(
    const array_1d<double, 6>& Metric,
    const IndexType NodeId
    )
{
    // The order is XX, XY, XZ, YY, YZ, ZZ
    KRATOS_ERROR_IF( MMG3D_Set_tensorSol(mMmgSol, Metric[0], Metric[3], Metric[5], Metric[1], Metric[4], Metric[2], NodeId) != 1 ) << "Unable to set tensor metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgUtilities<MMGLibrary::MMGS>::SetMetricTensor(
    const array_1d<double, 6>& Metric,
    const IndexType NodeId
    )
{
    // The order is XX, XY, XZ, YY, YZ, ZZ
    KRATOS_ERROR_IF( MMGS_Set_tensorSol(mMmgSol, Metric[0], Metric[3], Metric[5], Metric[1], Metric[4], Metric[2], NodeId) != 1 ) << "Unable to set tensor metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template class MmgUtilities<MMGLibrary::MMG2D>;
template class MmgUtilities<MMGLibrary::MMG3D>;
template class MmgUtilities<MMGLibrary::MMGS>;

}// namespace Kratos.
