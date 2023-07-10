// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Marc Nunez
//                   Carlos Roig
//                   Vicente Mataix Ferrandiz
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
#include "custom_utilities/parmmg/pmmg_utilities.h"
#include "mpi/includes/mpi_data_communicator.h"


// NOTE: The following contains the license of the PMMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017- .
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
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

template struct PMMGMeshInfo<PMMGLibrary::PMMG3D>;

/***********************************************************************************/
/***********************************************************************************/

// The member variables related with the PMMG library
PMMG_pParMesh mParMmgMesh;      /// The mesh data from PMMG

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::PrintAndGetParMmgMeshInfo(PMMGMeshInfo<TPMMGLibrary>& rPMMGMeshInfo)
{
    int np,ne,nprism,nt,nquad,na;
    KRATOS_ERROR_IF(PMMG_Get_meshSize(mParMmgMesh,&np,&ne,&nprism,&nt,&nquad,&na) != 1) << "Unable to get mesh size" << std::endl;

    /* Warning: mesh groups must be merged on each local partition */
    rPMMGMeshInfo.NumberOfNodes = np;
    if constexpr (TPMMGLibrary == PMMGLibrary::PMMG3D) { // 3D
        rPMMGMeshInfo.NumberOfTriangles = nt;
        rPMMGMeshInfo.NumberOfQuadrilaterals = nquad;
        rPMMGMeshInfo.NumberOfTetrahedra = ne;
        rPMMGMeshInfo.NumberOfPrism = nprism;
    }

    KRATOS_INFO_IF("ParMmgUtilities", GetEchoLevel() > 0) << "\tNodes created: " << rPMMGMeshInfo.NumberOfNodes << std::endl;
    if constexpr (TPMMGLibrary == PMMGLibrary::PMMG3D) { // 3D
        KRATOS_INFO_IF("ParMmgUtilities", GetEchoLevel() > 0) <<
        "Conditions created: " << rPMMGMeshInfo.NumberOfTriangles + rPMMGMeshInfo.NumberOfQuadrilaterals << "\n\tTriangles: " << rPMMGMeshInfo.NumberOfTriangles << "\tQuadrilaterals: " << rPMMGMeshInfo.NumberOfQuadrilaterals << "\n" <<
        "Elements created: " << rPMMGMeshInfo.NumberOfTetrahedra + rPMMGMeshInfo.NumberOfPrism << "\n\tTetrahedron: " << rPMMGMeshInfo.NumberOfTetrahedra << "\tPrisms: " << rPMMGMeshInfo.NumberOfPrism << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
IndexVectorType ParMmgUtilities<TPMMGLibrary>::FindDuplicateNodeIds(const ModelPart& rModelPart)
{
    return BaseType::FindDuplicateNodeIds(rModelPart);
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

        ids_triangles[0] = mLocalToGlobalNodePostMap[vertex_0];
        ids_triangles[1] = mLocalToGlobalNodePostMap[vertex_1];
        ids_triangles[2] = mLocalToGlobalNodePostMap[vertex_2];

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

        ids_tetrahedron[0] = mLocalToGlobalNodePostMap[vertex_0];
        ids_tetrahedron[1] = mLocalToGlobalNodePostMap[vertex_1];
        ids_tetrahedron[2] = mLocalToGlobalNodePostMap[vertex_2];
        ids_tetrahedron[3] = mLocalToGlobalNodePostMap[vertex_3];

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
Node::Pointer ParMmgUtilities<PMMGLibrary::PMMG3D>::CreateNode(
    ModelPart& rModelPart,
    const IndexType iNode,
    int& Ref,
    int& IsRequired
    )
{
    double coord_0, coord_1, coord_2;
    int is_corner;

    KRATOS_ERROR_IF(PMMG_Get_vertex(mParMmgMesh, &coord_0, &coord_1, &coord_2, &Ref, &is_corner, &IsRequired) != 1 ) << "Unable to get vertex" << std::endl;

    NodeType::Pointer p_node = Kratos::make_intrusive<NodeType>(iNode, coord_0, coord_1, coord_2);
    // Giving model part's variables list to the node
    p_node->SetSolutionStepVariablesList(rModelPart.pGetNodalSolutionStepVariablesList());
    //set buffer size
    p_node->SetBufferSize(rModelPart.GetBufferSize());

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

    if (rMapPointersRefCondition.find(Ref)==rMapPointersRefCondition.end()) {
        // Reference received that is not in the initial map of conditions. Returning a nullptr
        return nullptr;
    }

    if (rMapPointersRefCondition[Ref].get() == nullptr) {
        if (GetDiscretization() != DiscretizationOption::ISOSURFACE) { // The ISOSURFACE method creates new conditions from scratch, so we allow no previous Properties
            KRATOS_WARNING_IF("ParMmgUtilities", GetEchoLevel() > 1) << "Condition. Null reference-pointer returned. Rank " <<
                rModelPart.GetCommunicator().MyPID() << " Id: " << CondId << " Ref: " << Ref<< std::endl;

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
    if (mLocalToGlobalNodePostMap[vertex_0] == 0) SkipCreation = true;
    if (mLocalToGlobalNodePostMap[vertex_1] == 0) SkipCreation = true;
    if (mLocalToGlobalNodePostMap[vertex_2] == 0) SkipCreation = true;

    if (!SkipCreation) {
        std::vector<NodeType::Pointer> condition_nodes (3);
        condition_nodes[0] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_0]);
        condition_nodes[1] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_1]);
        condition_nodes[2] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_2]);

        p_condition = p_base_condition->Create(CondId, PointerVector<NodeType>{condition_nodes}, p_prop);
        if (p_base_condition->Is(MARKER)) p_condition->Set(MARKER);
    } else if (GetEchoLevel() > 2)
        KRATOS_WARNING_IF("ParMmgUtilities", GetEchoLevel() > 1) << "Condition creation avoided" << std::endl;

    if (p_condition != nullptr && Ref != 0) KRATOS_ERROR_IF(p_condition->GetGeometry().Area() < ZeroTolerance) << "Creating a almost zero or negative area condition with area " <<  p_condition->GetGeometry().Area()  << " and vertices: " << mLocalToGlobalNodePostMap[vertex_0] << " " << mLocalToGlobalNodePostMap[vertex_1] << " " << mLocalToGlobalNodePostMap[vertex_2]<< std::endl;
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

    if (GetDiscretization() == DiscretizationOption::ISOSURFACE) {
        // The existence of a _nullptr_ indicates an element that was removed. This is not an alarming indicator.
        if (rMapPointersRefElement[Ref].get() == nullptr) {
            // KRATOS_INFO("ParMmgUtilities") << "Element has been removed from domain. Ok." << std::endl;
            return p_element;
        } else {
            if (mLocalToGlobalNodePostMap[vertex_0] == 0) SkipCreation = true;
            if (mLocalToGlobalNodePostMap[vertex_1] == 0) SkipCreation = true;
            if (mLocalToGlobalNodePostMap[vertex_2] == 0) SkipCreation = true;
            if (mLocalToGlobalNodePostMap[vertex_3] == 0) SkipCreation = true;
            if (!SkipCreation) {
                std::vector<NodeType::Pointer> element_nodes (4);
                element_nodes[0] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_0]);
                element_nodes[1] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_1]);
                element_nodes[2] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_2]);
                element_nodes[3] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_3]);
                p_element = rMapPointersRefElement[Ref]->Create(ElemId, PointerVector<NodeType>{element_nodes}, rMapPointersRefElement[Ref]->pGetProperties());

                // Setting inside flag
                if (Ref == 2) {
                    p_element->Set(INSIDE, true);
                } else if (Ref == 3) {
                    p_element->Set(INSIDE, false);
                    if (GetRemoveRegions()) p_element->Set(TO_ERASE, true);
                }
            }
        }
    } else {
        // The existence of a _nullptr_ indicates a missing element. Two options are possible: error or replacement
        Properties::Pointer p_prop = nullptr;
        Element::Pointer p_base_element = nullptr;

        // Sometimes PMMG creates elements where there are not, then we skip
        if (rMapPointersRefElement[Ref].get() == nullptr) {
            return p_element;
        } else {
            p_base_element = rMapPointersRefElement[Ref];
            p_prop = p_base_element->pGetProperties();
        }

        // FIXME: This is not the correct solution to the problem, I asked in the PMMG Forum
        if (mLocalToGlobalNodePostMap[vertex_0] == 0) SkipCreation = true;
        if (mLocalToGlobalNodePostMap[vertex_1] == 0) SkipCreation = true;
        if (mLocalToGlobalNodePostMap[vertex_2] == 0) SkipCreation = true;
        if (mLocalToGlobalNodePostMap[vertex_3] == 0) SkipCreation = true;

        if (!SkipCreation) {
            std::vector<NodeType::Pointer> element_nodes (4);
            element_nodes[0] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_0]);
            element_nodes[1] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_1]);
            element_nodes[2] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_2]);
            element_nodes[3] = rModelPart.pGetNode(mLocalToGlobalNodePostMap[vertex_3]);

            p_element = p_base_element->Create(ElemId, PointerVector<NodeType>{element_nodes}, p_prop);
        } else if (GetEchoLevel() > 2)
            KRATOS_WARNING_IF("ParMmgUtilities", GetEchoLevel() > 1) << "Element creation avoided" << std::endl;
    }

    if (p_element!= nullptr) KRATOS_ERROR_IF(p_element->GetGeometry().Volume() < ZeroTolerance) << "Creating a almost zero or negative volume element" << std::endl;
    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::InitMesh(const DataCommunicator& rDataCommunicator)
{
    mParMmgMesh = nullptr;

    // We init the PMMG mesh and sol
    if (GetDiscretization() == DiscretizationOption::STANDARD) {
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "ParMMG requires a distributed DataCommunicator!" << std::endl;
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDefinedOnThisRank()) << "This rank is not part of this MPI_Comm!" << std::endl;
        MPI_Comm the_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
        PMMG_Init_parMesh( PMMG_ARG_start, PMMG_ARG_ppParMesh, &mParMmgMesh, PMMG_ARG_pMesh, PMMG_ARG_pMet, PMMG_ARG_dim, 3, PMMG_ARG_MPIComm, the_mpi_comm, PMMG_ARG_end);
    } else {
        KRATOS_ERROR << "Discretization type: " << static_cast<int>(GetDiscretization()) << " not fully implemented" << std::endl;
    }

    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::InitVerbosity()
{
    BaseType::InitVerbosity();
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
    KRATOS_ERROR_IF( !PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_globalNum, nodeGloNum) ) << "Unable to set node global numbering" << std::endl;
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
    KRATOS_ERROR << "SetDispSizeVector not yet implemented in ParMmg" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::CheckMeshData()
{
    KRATOS_ERROR << "CheckMeshData not yet implemented in ParMmg" << std::endl;
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

    // // b) function calling
    KRATOS_INFO_IF("ParMmgUtilities", PMMG_saveMesh_distributed(mParMmgMesh,mesh_file) != 1) << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::OutputSol(const std::string& rOutputName)
{
     const std::string sol_name = rOutputName + ".sol";
     const char* sol_file = sol_name.c_str();

    // // b) Function calling
    KRATOS_INFO_IF("ParMmgUtilities", PMMG_saveMet_distributed(mParMmgMesh, sol_file) != 1)<< "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::OutputDisplacement(const std::string& rOutputName)
{
    KRATOS_ERROR << "OutputDisplacement is not yet implemented in ParMmg" << std::endl;
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
    BaseType::OutputReferenceEntitities(rOutputName, rRefCondition, rRefElement);
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
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_nosurf, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_surf_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no surfacic modifications" << std::endl;

    // Don't insert nodes on mesh
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_noinsert, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_insert_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no insertion/suppression point" << std::endl;

    // Don't swap mesh
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_noswap, static_cast<int>(ConfigurationParameters["advanced_parameters"]["no_swap_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // Number Of iterations
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_niter, static_cast<int>(ConfigurationParameters["advanced_parameters"]["number_of_iterations"].GetInt())) != 1 )
        KRATOS_ERROR << "Unable to set number of remeshing iterations" << std::endl;

    // Mesh Size
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_meshSize, static_cast<int>(ConfigurationParameters["advanced_parameters"]["mesh_size"].GetInt())) != 1 )
        KRATOS_ERROR << "Unable to set target mesh size of Mmg" << std::endl;

    // Mesh Ratio
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_metisRatio, static_cast<int>(ConfigurationParameters["advanced_parameters"]["metis_ratio"].GetInt())) != 1 )
        KRATOS_ERROR << "Unable to set wanted ratio # mesh / # metis super nodes" << std::endl;

    // hgradreq ( To Be added)
    // if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_hgradreq, static_cast<int>(ConfigurationParameters["advanced_parameters"]["hgradreq"].GetDouble())) != 1 )
    //     KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // API Mode: Forced to "1" all the time
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_APImode, static_cast<int>(PMMG_APIDISTRIB_nodes)) != 1 )
        KRATOS_ERROR << "Unable to set initialize parallel library through interface faces or nodes" << std::endl;

    // Node global numbering: Forced to "1" all the time
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_globalNum,1) != 1 )
        KRATOS_ERROR << "Unable to set initialize parallel library with node global numbering" << std::endl;

    // Set the angle detection
    const bool deactivate_detect_angle = ConfigurationParameters["advanced_parameters"]["deactivate_detect_angle"].GetBool();
    if ( deactivate_detect_angle) {
        if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_angle, static_cast<int>(!deactivate_detect_angle)) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection on" << std::endl;
    }

    // Set the value for angle detection (default 45°)
    if (ConfigurationParameters["advanced_parameters"]["force_angle_detection_value"].GetBool()) {
        if ( PMMG_Set_dparameter(mParMmgMesh,PMMG_DPARAM_angleDetection, ConfigurationParameters["advanced_parameters"]["angle_detection_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection value" << std::endl;
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
    int ier = PMMG_parmmglib_distributed(mParMmgMesh);

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
        const IndexType id_1 = mGlobalToLocalNodePreMap[rGeometry[0].Id()]; // First node Id
        const IndexType id_2 = mGlobalToLocalNodePreMap[rGeometry[1].Id()]; // Second node Id
        const IndexType id_3 = mGlobalToLocalNodePreMap[rGeometry[2].Id()]; // Third node Id

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
        KRATOS_ERROR << "Not yet implemented in ParMmg" << std::endl;
        //KRATOS_ERROR_IF( PMMG_Set_quadrilateral(mParMmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 ) << "Unable to set quadrilateral" << std::endl;
    } else {
        const SizeType size_geometry = rGeometry.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << static_cast<int>(rGeometry.GetGeometryType()) << std::endl;
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
    const IndexType id_1 = mGlobalToLocalNodePreMap[rGeometry[0].Id()]; // First node Id
    const IndexType id_2 = mGlobalToLocalNodePreMap[rGeometry[1].Id()]; // Second node Id
    const IndexType id_3 = mGlobalToLocalNodePreMap[rGeometry[2].Id()]; // Third node Id
    const IndexType id_4 = mGlobalToLocalNodePreMap[rGeometry[3].Id()]; // Fourth node Id

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
    KRATOS_ERROR_IF( PMMG_Set_scalarMet( mParMmgMesh, Metric, mGlobalToLocalNodePreMap[NodeId]) != 1 ) << "Unable to set scalar metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetMetricVector(
    const array_1d<double, 3>& rMetric,
    const IndexType NodeId
    )
{
    KRATOS_ERROR_IF( PMMG_Set_vectorMet( mParMmgMesh, rMetric[0], rMetric[1], rMetric[2], mGlobalToLocalNodePreMap[NodeId]) != 1 ) << "Unable to set vector metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetMetricTensor(
    const array_1d<double, 6>& rMetric,
    const IndexType NodeId
    )
{
    // The order is XX, XY, XZ, YY, YZ, ZZ
    KRATOS_ERROR_IF( PMMG_Set_tensorMet( mParMmgMesh, rMetric[0], rMetric[3], rMetric[5], rMetric[1], rMetric[4], rMetric[2], mGlobalToLocalNodePreMap[NodeId]) != 1 ) << "Unable to set tensor metric" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ParMmgUtilities<PMMGLibrary::PMMG3D>::SetDisplacementVector(
    const array_1d<double, 3>& rDisplacement,
    const IndexType NodeId
    )
{
    KRATOS_ERROR << "Not yet implemented in ParMmg" << std::endl;
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
    KRATOS_ERROR << "Not yet implemented in ParMmg" << std::endl;
    // KRATOS_ERROR_IF( PMMG_Get_vectorSol(mParMmgDisp, &rDisplacement[0], &rDisplacement[1], &rDisplacement[2]) != 1 ) << "Unable to get vector displacement" << std::endl;
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

        KRATOS_WARNING_IF("ParMmgUtilities", GetEchoLevel() > 0 && (r_sub_model_part.NumberOfNodes() > 0 && (r_sub_model_part.NumberOfConditions() == 0 && r_sub_model_part.NumberOfElements() == 0))) <<
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
    std::unordered_set<IndexType> total_nodes_to_remesh;

    /* Gathering LOCAL mesh info */
    PMMGMeshInfo<TPMMGLibrary> pmmg_mesh_info;
    if constexpr (TPMMGLibrary == PMMGLibrary::PMMG3D) { // 3D
        /* Conditions */
        std::size_t num_tri = 0, num_quad = 0;
        for(IndexType i = 0; i < r_conditions_array.size(); ++i) {
            auto it_cond = it_cond_begin + i;
            if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) { // Triangles
                for (auto& r_node : it_cond->GetGeometry()) {
                    total_nodes_to_remesh.insert(r_node.Id());
                }
                num_tri += 1;
            } else {
                KRATOS_ERROR << "ParMmg currently only supports triangles on conditions. Your geometry type was: " << static_cast<int>((it_cond->GetGeometry()).GetGeometryType()) <<  std::endl;
            }
        }

        pmmg_mesh_info.NumberOfTriangles = num_tri;
        pmmg_mesh_info.NumberOfQuadrilaterals = num_quad;

        KRATOS_INFO_IF("ParMmgUtilities", ((num_tri + num_quad) < r_conditions_array.size()) && GetEchoLevel() > 0) <<
        "Number of Conditions: " << r_conditions_array.size() << " Number of Triangles: " << num_tri << " Number of Quadrilaterals: " << num_quad << std::endl;

        /* Elements */
        std::size_t num_tetra = 0, num_prisms = 0;
        for(IndexType i = 0; i < r_elements_array.size(); ++i) {
            auto it_elem = it_elem_begin + i;

            if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) { // Tetrahedron
                for (auto& r_node : it_elem->GetGeometry()) {
                    total_nodes_to_remesh.insert(r_node.Id());
                }
                num_tetra += 1;
            } else {
                KRATOS_ERROR << "ParMmg currently only supports tetrahedras on elements. Your geometry type was: " << static_cast<int>((it_elem->GetGeometry()).GetGeometryType()) <<  std::endl;
            }
        }

        pmmg_mesh_info.NumberOfTetrahedra = num_tetra;
        pmmg_mesh_info.NumberOfPrism = num_prisms;

        KRATOS_INFO_IF("ParMmgUtilities", ((num_tetra + num_prisms) < r_elements_array.size()) && GetEchoLevel() > 0) <<
        "Number of Elements: " << r_elements_array.size() << " Number of Tetrahedron: " << num_tetra << " Number of Prisms: " << num_prisms << std::endl;
    }

    pmmg_mesh_info.NumberOfNodes = total_nodes_to_remesh.size();
    SetMeshSize(pmmg_mesh_info);

    mGlobalToLocalNodePreMap = std::unordered_map<int, int>();
    mGlobalToLocalElemPreMap = std::unordered_map<int, int>();
    mGlobalToLocalCondPreMap = std::unordered_map<int, int>();

    IndexType counter_remesh = 0;
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;
        mGlobalToLocalNodePreMap[it_node->Id()] = ++counter_remesh;
    }
    counter_remesh = 0;
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = it_cond_begin + i;
        mGlobalToLocalCondPreMap[it_cond->Id()] = ++counter_remesh;
    }
    counter_remesh = 0;
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;
        mGlobalToLocalElemPreMap[it_elem->Id()] = ++counter_remesh;
    }

    // Now we compute the colors
    rColors.clear();
    ColorsMapType nodes_colors, cond_colors, elem_colors;
    AssignUniqueModelPartCollectionTagUtility model_part_collections(rModelPart);
    model_part_collections.ComputeTags(nodes_colors, cond_colors, elem_colors, rColors);

    /* Nodes */
    block_for_each(r_nodes_array, nodes_colors,
        [this,&Framework](NodeType& rNode, ColorsMapType& nodes_colors) {

        const array_1d<double, 3>& r_coordinates = Framework == FrameworkEulerLagrange::LAGRANGIAN ? rNode.GetInitialPosition() : rNode.Coordinates();
        SetNodes(r_coordinates[0], r_coordinates[1], r_coordinates[2], nodes_colors[rNode.Id()], mGlobalToLocalNodePreMap[rNode.Id()]);
    });

    /* Conditions */
    block_for_each(r_conditions_array, cond_colors,
        [this](Condition& rCondition, ColorsMapType& cond_colors) {

        SetConditions(rCondition.GetGeometry(), cond_colors[rCondition.Id()], mGlobalToLocalCondPreMap[rCondition.Id()]);
    });

    /* Elements */
    block_for_each(r_elements_array, elem_colors,
        [this](Element& rElement, ColorsMapType& elem_colors) {

        SetElements(rElement.GetGeometry(), elem_colors[rElement.Id()], mGlobalToLocalElemPreMap[rElement.Id()]);
    });

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
    std::size_t icomm, ncomm;

    ncomm = 0;
    for(std::size_t i = 0; i < neighbour_indices.size(); i++) {
        if ((neighbour_indices[i]) >= 0) {
            ncomm++;
        }
    }
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
                localId.push_back(mGlobalToLocalNodePreMap[node.Id()]);
            }

            for(auto& node: rModelPart.GetCommunicator().GhostMesh(i).Nodes()) {
                globalId.push_back(node.Id());
                localId.push_back(mGlobalToLocalNodePreMap[node.Id()]);
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
                localId.push_back(mGlobalToLocalNodePreMap[node.Id()]);
            }

            for(auto& node: rModelPart.GetCommunicator().GhostMesh(i).Nodes()) {
                globalId.push_back(node.Id());
                localId.push_back(mGlobalToLocalNodePreMap[node.Id()]);
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
        rRefCondition[0] = r_clone_condition.Create(0, it_cond_begin->pGetGeometry(), it_cond_begin->pGetProperties());
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
    BaseType::GenerateSolDataFromModelPart(rModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::GenerateDisplacementDataFromModelPart(ModelPart& rModelPart)
{
    KRATOS_ERROR << "Lagrangian approach not supported by ParMmg" << std::endl;
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
    const IndexType rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
    const IndexType size = rModelPart.GetCommunicator().GetDataCommunicator().Size();

    rModelPart.GetCommunicator().LocalMesh().Nodes().clear();
    rModelPart.GetCommunicator().InterfaceMesh().Nodes().clear();
    rModelPart.GetCommunicator().GhostMesh().Nodes().clear();

    std::vector<int> array_of_local_elements(size,0);
    std::vector<int> reduced_array_of_local_elements(size,0);
    if (rank != (size-1)) {
        array_of_local_elements[rank+1]=rPMMGMeshInfo.NumberFirstTypeElements();
    }

    rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(array_of_local_elements, reduced_array_of_local_elements);

    for (IndexType i = 1; i < size;i++){
        reduced_array_of_local_elements[i] += reduced_array_of_local_elements[i-1];
    }

    std::unordered_map<int,int> local_id_to_partition_index;

    for (IndexType i_node = 1; i_node <= rPMMGMeshInfo.NumberOfNodes; ++i_node) {
        PMMG_Get_vertexGloNum(mParMmgMesh,&mLocalToGlobalNodePostMap[i_node],&local_id_to_partition_index[i_node]);
    }

    // Create a new model part // TODO: Use a different kind of element for each submodelpart (in order to be able of remeshing more than one kind o element or condition)
    std::unordered_map<IndexType, IndexVectorType> color_nodes, first_color_cond, second_color_cond, first_color_elem, second_color_elem;

    // The tempotal store of
    ConditionsArrayType created_conditions_vector;
    ElementsArrayType created_elements_vector;
    NodesArrayType created_nodes_vector;

    // Auxiliar values
    int ref, is_required;

    /* NODES */ // TODO: ADD OMP
    for (IndexType i_node = 1; i_node <= rPMMGMeshInfo.NumberOfNodes; ++i_node) {

        KRATOS_ERROR_IF(mLocalToGlobalNodePostMap[i_node]==0) << "0 global index at local node: " << i_node <<  std::endl;

        NodeType::Pointer p_node = CreateNode(rModelPart, mLocalToGlobalNodePostMap[i_node], ref, is_required);

        IndexType partition_index = local_id_to_partition_index[i_node];
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = partition_index;

        KRATOS_ERROR_IF(partition_index<0) << "Negative partition index:: " << partition_index << " at node: " << p_node->Id() << std::endl;
        KRATOS_ERROR_IF(partition_index>size-1) << "Partition index is gretear than greater than the MPI-Size: " << partition_index << " at node: " << p_node->Id() << std::endl;

        // Set the DOFs in the nodes
        for (auto it_dof = rDofs.begin(); it_dof != rDofs.end(); ++it_dof)
            p_node->pAddDof(**it_dof);

        created_nodes_vector.push_back(p_node);
        if (ref != 0) color_nodes[static_cast<IndexType>(ref)].push_back(mLocalToGlobalNodePostMap[i_node]);// NOTE: ref == 0 is the MainModelPart
    }
    rModelPart.AddNodes(created_nodes_vector.begin(), created_nodes_vector.end());

    /* CONDITIONS */ // TODO: ADD OMP
    if (rMapPointersRefCondition.size() > 0) {
        IndexType counter_first_cond = 0;

        std::unordered_map<IndexType,IndexType> local_to_global_triangles;
        std::unordered_map<IndexType,IndexType> local_to_owner_triangles;
        for (IndexType i_cond = 1; i_cond <= rPMMGMeshInfo.NumberFirstTypeConditions(); ++i_cond) {

            int global_id, owner;
            KRATOS_ERROR_IF(PMMG_Get_triangleGloNum(mParMmgMesh, &global_id, &owner) !=1)<< "Unable to get triangle global id" << std::endl;
            local_to_global_triangles[i_cond] = global_id;
            local_to_owner_triangles[i_cond] = owner;
        }
        const IndexVectorType first_condition_to_remove = CheckFirstTypeConditions();
        for (IndexType i_cond = 1; i_cond <= rPMMGMeshInfo.NumberFirstTypeConditions(); ++i_cond) {
            bool skip_creation = false;
            if (counter_first_cond < first_condition_to_remove.size()) {
                if (first_condition_to_remove[counter_first_cond] == i_cond) {
                    skip_creation = true;
                    counter_first_cond += 1;
                }
            }
            IndexType id_to_use = local_to_global_triangles[i_cond];

            Condition::Pointer p_condition = CreateFirstTypeCondition(rModelPart, rMapPointersRefCondition, id_to_use, ref, is_required, skip_creation);

            if (p_condition.get() != nullptr && ref != 0) {
                created_conditions_vector.push_back(p_condition);
                first_color_cond[static_cast<IndexType>(ref)].push_back(id_to_use);// NOTE: ref == 0 is the MainModelPart
            }
        }
    }
    rModelPart.AddConditions(created_conditions_vector.begin(), created_conditions_vector.end());

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

            Element::Pointer p_element = CreateFirstTypeElement(rModelPart, rMapPointersRefElement, elem_id+reduced_array_of_local_elements[rank], ref, is_required, skip_creation);

            if (p_element.get() != nullptr) {
                created_elements_vector.push_back(p_element);
                if (ref != 0) first_color_elem[static_cast<IndexType>(ref)].push_back(elem_id+reduced_array_of_local_elements[rank]);// NOTE: ref == 0 is the MainModelPart
                elem_id += 1;
            }
        }
    }
    rModelPart.AddElements(created_elements_vector.begin(), created_elements_vector.end());

    // We add nodes, conditions and elements to the sub model parts

    for (auto& r_color_list : rColors) {
        const IndexType key = r_color_list.first;


        if (key != 0) {// NOTE: key == 0 is the MainModelPart
            for (auto sub_model_part_name : r_color_list.second) {
                ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(rModelPart, sub_model_part_name);

                if (color_nodes.find(key) != color_nodes.end())
                    r_sub_model_part.AddNodes(color_nodes[key]);
                if (first_color_cond.find(key) != first_color_cond.end())
                    r_sub_model_part.AddConditions(first_color_cond[key]);
                if (second_color_cond.find(key) != second_color_cond.end())
                    r_sub_model_part.AddConditions(second_color_cond[key]);
                if (first_color_elem.find(key) != first_color_elem.end())
                    r_sub_model_part.AddElements(first_color_elem[key]);
                if (second_color_elem.find(key) != second_color_elem.end())
                    r_sub_model_part.AddElements(second_color_elem[key]);
            }
        }
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

        std::unordered_map<IndexType, std::vector<IndexType>> pi_to_vector_node;
        //get nodes we need from other ranks
        for(int i_node = 0; i_node < static_cast<int>(r_sub_model_part.Nodes().size()); ++i_node ){
            auto it_node = r_sub_model_part.NodesBegin() + i_node;
            IndexType partition_index = it_node->FastGetSolutionStepValue(PARTITION_INDEX);
            if (partition_index != rank) {
                pi_to_vector_node[partition_index].push_back(it_node->Id());
            }
        }

        std::vector<IndexType> receive_nodes_all_ranks;
        for (IndexType i_rank=0; i_rank < size; i_rank++ ) {
            std::vector<IndexType> receive_nodes;
            if (i_rank != rank){
               receive_nodes = rModelPart.GetCommunicator().GetDataCommunicator().SendRecv(pi_to_vector_node[i_rank], i_rank, i_rank);
               receive_nodes_all_ranks.insert(receive_nodes_all_ranks.end(), receive_nodes.begin(), receive_nodes.end());
            }
        }

        for (IndexType i_node = 0; i_node<receive_nodes_all_ranks.size(); i_node++)
        {
            if (r_sub_model_part.Nodes().find(receive_nodes_all_ranks[i_node]) != r_sub_model_part.Nodes().end()) {
                receive_nodes_all_ranks.erase(receive_nodes_all_ranks.begin() + i_node);
            }
        }

        r_sub_model_part.AddNodes(receive_nodes_all_ranks);


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
    BaseType::WriteReferenceEntitities(rModelPart, rFilename, rRefCondition, rRefElement);
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::CreateAuxiliarSubModelPartForFlags(ModelPart& rModelPart)
{
    KRATOS_TRY;

    auto& data_comm = rModelPart.GetCommunicator().GetDataCommunicator();

    ModelPart& r_auxiliar_model_part = rModelPart.CreateSubModelPart("AUXILIAR_MODEL_PART_TO_LATER_REMOVE");

    const auto& r_flags = KratosComponents<Flags>::GetComponents();

    for (auto& r_flag : r_flags) {
        const std::string name_sub_model = "FLAG_" + r_flag.first;
        if (name_sub_model.find("NOT") == std::string::npos && name_sub_model.find("ALL") == std::string::npos) { // Avoiding inactive flags
            r_auxiliar_model_part.CreateSubModelPart(name_sub_model);
            ModelPart& r_auxiliar_sub_model_part = r_auxiliar_model_part.GetSubModelPart(name_sub_model);
            FastTransferBetweenModelPartsProcess(r_auxiliar_sub_model_part, rModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, *(r_flag.second)).Execute();

            const bool is_aux_sub_model_empty = r_auxiliar_sub_model_part.NumberOfNodes() == 0
                && r_auxiliar_sub_model_part.NumberOfElements() == 0
                && r_auxiliar_sub_model_part.NumberOfConditions() == 0;

            // Remove sub model part if it is empty in all ranks
            if (data_comm.MinAll(is_aux_sub_model_empty)) {
                r_auxiliar_model_part.RemoveSubModelPart(name_sub_model);
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgUtilities<TPMMGLibrary>::AssignAndClearAuxiliarSubModelPartForFlags(ModelPart& rModelPart)
{
   BaseType::AssignAndClearAuxiliarSubModelPartForFlags(rModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
std::unordered_map<int, int> ParMmgUtilities<TPMMGLibrary>::GetNodalLocalToGlobalMap() {
    return mLocalToGlobalNodePostMap;
}

template class ParMmgUtilities<PMMGLibrary::PMMG3D>;

}// namespace Kratos.
