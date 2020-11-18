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

#if !defined(PMMG_INCLUDES)
#define PMMG_INCLUDES
#include "parmmg/libparmmg.h"
#endif /* PMMG_INCLUDES defined */

// System includes

// External includes

// Project includes
#include "custom_utilities/parmmg/pmmg_utilities.h"

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

template<>
SizeType PMMGMeshInfo<PMMGLibrary::PMMG3D>::NumberFirstTypeConditions() const
{
    return NumberOfTriangles;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType PMMGMeshInfo<PMMGLibrary::PMMG3D>::NumberSecondTypeConditions() const
{
    return NumberOfQuadrilaterals;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType PMMGMeshInfo<PMMGLibrary::PMMG3D>::NumberFirstTypeElements() const
{
    return NumberOfTetrahedra;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType PMMGMeshInfo<PMMGLibrary::PMMG3D>::NumberSecondTypeElements() const
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
    KRATOS_ERROR << "OutputDisplacement not yet Implemented in ParMmg" << std::endl;
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
    if ( PMMG_Set_iparameter(mParMmgMesh,PMMG_IPARAM_globalNum,1) != 1 )
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

    ier = PMMG_parmmglib_distributed(mParMmgMesh);

    if ( ier == PMMG_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF PMMG3DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == PMMG_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF PMMG3DLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

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
    KRATOS_ERROR << "PMMG_Set_vectorSol not yet implemented in ParMmg" << std::endl;
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
    KRATOS_ERROR << "PMMG_Get_vectorSol not yet implemented in ParMmg" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template class ParMmgUtilities<PMMGLibrary::PMMG3D>;

}// namespace Kratos.
