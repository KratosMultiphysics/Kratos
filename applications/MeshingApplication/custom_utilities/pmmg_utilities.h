// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Marc Nunez, Carlos Roig, Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PMMG_UTILITIES)
#define KRATOS_PMMG_UTILITIES

// System includes

// External includes

// Project includes
#include "meshing_application.h"

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
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// Index definition
    typedef std::size_t                                               IndexType;

    /// Size definition
    typedef std::size_t                                                SizeType;

    /// Index vector
    typedef std::vector<IndexType>                              IndexVectorType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @struct PMMGMeshInfo
 * @ingroup MeshingApplication
 * @brief Stores the ParMmg mesh information
 * @author Vicente Mataix Ferrandiz
 */
template<PMMGLibrary TPMMGLibrary>
struct PMMGMeshInfo
{
    // Info stored
    SizeType NumberOfNodes;
    SizeType NumberOfLines;
    SizeType NumberOfTriangles;
    SizeType NumberOfQuadrilaterals;
    SizeType NumberOfPrism;
    SizeType NumberOfTetrahedra;

    /**
     * @brief It returns the number of the first type of conditions
     */
    const SizeType NumberFirstTypeConditions() const;

    /**
     * @brief It returns the number of the second type of conditions
     */
    const SizeType NumberSecondTypeConditions() const;

    /**
     * @brief It returns the number of the first type of elements
     */
    const SizeType NumberFirstTypeElements() const;

    /**
     * @brief It returns the number of the second type of elements
     */
    const SizeType NumberSecondTypeElements() const;
};

/**
 * @class ParMmgUtilities
 * @ingroup MeshingApplication
 * @brief Provides the Kratos interface to the PMMG library API
 * @author Vicente Mataix Ferrandiz
 */
template<PMMGLibrary TPMMGLibrary>
class KRATOS_API(MESHING_APPLICATION) ParMmgUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParMmgUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ParMmgUtilities);

    /// Node definition
    typedef Node <3>                                                   NodeType;
    // Geometry definition
    typedef Geometry<NodeType>                                     GeometryType;

    /// Spatial dimension
    static constexpr SizeType Dimension = 3;

    /// The type of array considered for the tensor
    typedef typename std::conditional<Dimension == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    /// Double vector
    typedef std::vector<double> DoubleVectorType;

    /// Double vector map
    typedef std::unordered_map<DoubleVectorType, IndexType, KeyHasherRange<DoubleVectorType>, KeyComparorRange<DoubleVectorType> > DoubleVectorMapType;

    /// Index vector map
    typedef std::unordered_map<IndexVectorType, IndexType, KeyHasherRange<IndexVectorType>, KeyComparorRange<IndexVectorType> > IndexVectorMapType;

    /// Colors map
    typedef std::unordered_map<IndexType,IndexType> ColorsMapType;

    /// Index pair
    typedef std::pair<IndexType,IndexType> IndexPairType;

    /// Index and string vector pair
    typedef std::pair<IndexType, std::vector<std::string>> IndexStringVectorPairType;

    /// Definition of the zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief This method sets the echo level
     * @param[in] EchoLevel Sets the echo level
     */
    void SetEchoLevel(const SizeType EchoLevel);

    /**
     * @brief This method gets the echo level
     * @return mEchoLevel Gets the echo level
     */
    SizeType GetEchoLevel();

      /**
     * @brief This method sets the discretization method
     * @param[in] Discretization Sets the discretization method
     */
    void SetDiscretization(const DiscretizationOption Discretization);

    /**
     * @brief This method gets the discretization method
     * @return mDiscretization Gets the discretization method
     */
    DiscretizationOption GetDiscretization();

    /**
     * @brief This method sets if the regions must be removed
     * @param[in] RemoveRegions Sets if the regions must be removed
     */
    void SetRemoveRegions(const bool RemoveRegions);

    /**
     * @brief This method gets if the regions must be removed
     * @return mRemoveRegions Gets if the regions must be removed
     */
    bool GetRemoveRegions();

    /**
     * @brief It prints info about the current mesh
     * @param[in,out] rPMMGMeshInfo The number of nodes, conditions and elements
     */
    void PrintAndGetParMmgMeshInfo(PMMGMeshInfo<TPMMGLibrary>& rPMMGMeshInfo);

    /**
     * @brief Returns a vector of ids of spatially repeated nodes
     * @param[in,out] rModelPart The model part whose nodes are checked
     */
    IndexVectorType FindDuplicateNodeIds(const ModelPart& rModelPart);

    /**
     * @brief Returns a vector of ids of repeated conditions
     * @details For 2D and surface meshes it returns the ids of repeated edges found in the PMMG mesh data structure. In 3D it returns ids of triangles.
     */
    IndexVectorType CheckFirstTypeConditions();

    /**
     * @brief Returns a vector of ids of repeated conditions
     * @details For 3D meshes it returns the ids of repeated quadrilaterals found in the PMMG mesh data structure. Otherwise, it returns an empty vector.
     */
    IndexVectorType CheckSecondTypeConditions();

    /**
     * @brief Returns a vector of ids of repeated elements
     * @details For 2D and surface meshes it returns the ids of repeated triangles found in the PMMG mesh data structure. In 3D it returns ids of tetrahedra.
     */
    IndexVectorType CheckFirstTypeElements();

    /**
     * @brief Returns a vector of ids of repeated elements
     * @details For 3D meshes it returns the ids of repeated prisms found in the PMMG mesh data structure. Otherwise, it returns an empty vector.
     */
    IndexVectorType CheckSecondTypeElements();

    /**
     * @brief It blocks certain nodes before remesh the model
     * @param[in] iNode The index of the node
     */
    void BlockNode(const IndexType iNode);

    /**
     * @brief It blocks certain conditions before remesh the model
     * @param[in] iCondition The index of the condition
     */
    void BlockCondition(const IndexType iCondition);

    /**
     * @brief It blocks certain elements before remesh the model
     * @param[in] iElement The index of the element
     */
    void BlockElement(const IndexType iElement);

     /**
     * @brief Here the API mode is set using the API
     * @param[in] API Mode sets the mode in which the parallel communicator works
     */
    void InitAPIModeParameter(const IndexType VerbosityPMMG);

    /**
     * @brief Ask for output node global numbering
     * @param[in] 1 to ask it, 0 to skip it
     */
    void InitNodeGloNumParameter(const IndexType nodeGloNum);

    /**
     * @brief Here the verbosity is set using the API
     * @param[in] VerbosityPMMG The equivalent verbosity level in the PMMG API
     */
    void InitVerbosityParameter(const IndexType VerbosityPMMG);

    /**
     * @brief This sets the size of the mesh
     * @param[in,out] rPMMGMeshInfo The number of nodes, conditions and elements
     */
    void SetMeshSize(PMMGMeshInfo<TPMMGLibrary>& rPMMGMeshInfo);

    /**
     * @brief This sets the size of the solution for the scalar case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeScalar(const SizeType NumNodes);

    /**
     * @brief This sets the size of the solution for the vector case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeVector(const SizeType NumNodes);

    /**
     * @brief This sets the size of the solution for the tensor case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeTensor(const SizeType NumNodes);

    /**
     * @brief This sets the size of the displacement for lagrangian movement
     * @param[in] NumNodes Number of nodes
     */
    void SetDispSizeVector(const SizeType NumNodes);

    /**
     * @brief This checks the mesh data and prints if it is OK
     */
    void CheckMeshData();

      /**
     * @brief This sets the output mesh
     * @param[in] rInputName The input name
     */
    void InputMesh(const std::string& rInputName);

    /**
     * @brief This sets the output sol
     * @param[in] rInputName The input name
     */
    void InputSol(const std::string& rInputName);

    /**
     * @brief This sets the output mesh
     * @param[in] rOutputName The output name
     */
    void OutputMesh(const std::string& rOutputName);

    /**
     * @brief This sets the output sol
     * @param[in] rOutputName The output name
     */
    void OutputSol(const std::string& rOutputName);

    /**
     * @brief This sets the output displacement
     * @param[in] rOutputName The output name
     */
    void OutputDisplacement(const std::string& rOutputName);

    /**
     * @brief This method generates the maps of reference for conditions and elements from an existing json
     * @param[in] rOutputName The name of the files
     * @param[in] rRefCondition The conditions of reference
     * @param[in] rRefElement The elements of reference
     */
    void OutputReferenceEntitities(
        const std::string& rOutputName,
        const std::unordered_map<IndexType,Condition::Pointer>& rRefCondition,
        const std::unordered_map<IndexType,Element::Pointer>& rRefElement
        );

    /**
     * @brief This frees the PMMG structures
     */
    void FreeAll();

    /**
     * @brief This loads the solution
     */
    void PMMGLibCallMetric(Parameters ConfigurationParameters);

    /**
     * @brief This loads the solution
     */
    void PMMGLibCallIsoSurface(Parameters ConfigurationParameters);

    /**
     * @brief This sets the nodes of the mesh
     * @param[in] X Coordinate X
     * @param[in] Y Coordinate Y
     * @param[in] Z Coordinate Z
     * @param[in] Color Reference of the node(submodelpart)
     * @param[in] Index The index number of the node
     */
    void SetNodes(
        const double X,
        const double Y,
        const double Z,
        const IndexType Color,
        const IndexType Index
        );

    /**
     * @brief This sets the conditions of the mesh
     * @param[in,out] rGeometry The geometry of the condition
     * @param[in] Color Reference of the node(submodelpart)
     * @param[in] Index The index number of the node
     */
    void SetConditions(
        GeometryType& rGeometry,
        const IndexType Color,
        const IndexType Index
        );

    /**
     * @brief This sets elements of the mesh
     * @param[in,out] rGeometry The geometry of the element
     * @param[in] Color Reference of the node(submodelpart)
     * @param[in] Index The index number of the node
     */
    void SetElements(
        GeometryType& rGeometry,
        const IndexType Color,
        const IndexType Index
        );

    /**
     * @brief This function is used to set the metric scalar
     * @param[in] Metric The inverse of the size node
     * @param[in] NodeId The id of the node
     */
    void SetMetricScalar(
        const double Metric,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to set the metric vector (x, y, z)
     * @param[in] rMetric This array contains the components of the metric vector
     * @param[in] NodeId The id of the node
     */
    void SetMetricVector(
        const array_1d<double, Dimension>& rMetric,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to set the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param[in] rMetric This array contains the components of the metric tensor in the PMMG defined order
     * @param[in] NodeId The id of the node
     */
    void SetMetricTensor(
        const TensorArrayType& rMetric,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to set the displacement vector (x, y, z)
     * @param[in] rMetric This array contains the components of the displacement vector
     * @param[in] NodeId The id of the node
     */
    void SetDisplacementVector(
        const array_1d<double, 3>& rDisplacement,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to retrieve the metric scalar
     * @param[in,out] rMetric The inverse of the size node
     */
    void GetMetricScalar(double& rMetric);

    /**
     * @brief This function is used to retrieve the metric vector (x, y, z)
     * @param[in,out] rMetric This array contains the components of the metric vector
     */
    void GetMetricVector(array_1d<double, Dimension>& rMetric);

    /**
     * @brief This function is used to retrieve the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param[in,out] rMetric This array contains the components of the metric tensor in the PMMG defined order
     */
    void GetMetricTensor(TensorArrayType& rMetric);

    /**
     * @brief This function is used to retrieve the displacement vector (x, y, z)
     * @param[in,out] rDisplacement This array contains the components of the displacement vector
     */
    void GetDisplacementVector(array_1d<double, 3>& rDisplacement);

    /**
     * @brief This function reorder the nodes, conditions and elements to avoid problems with non-consecutive ids
     * @param[in,out] rModelPart The model part of interest to study
     */

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const // override
    {
        return "ParMmgUtilities";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const // override
    {
        rOStream << "ParMmgUtilities";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const // override
    {
    }

protected:

    ///@name Protected Member Variables
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
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    SizeType mEchoLevel = 0;                                               /// The echo level of the utilities
    bool mRemoveRegions = false;                                           /// Cutting-out specified regions during surface remeshing
    DiscretizationOption mDiscretization = DiscretizationOption::STANDARD; /// Discretization The discretization type
    std::map<int, int> mGlobalToLocalNodePreMap;                           /// Map of nodal global ids to local ids before remeshing
    std::map<int, int> mGlobalToLocalElemPreMap;                           /// Map of elemental global ids to local ids before remeshing
    std::map<int, int> mGlobalToLocalCondPreMap;                           /// Map of condition global ids to local ids before remeshing
    std::map<int, int> mLocalToGlobalNodePostMap;                          /// Map of nodal local ids to global ids after remeshing

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
    ///@name Un accessible methods
    ///@{

//     /// Assignment operator.
//     ParMmgUtilities& operator=(ParMmgUtilities const& rOther);

//     /// Copy constructor.
//     ParMmgUtilities(ParMmgUtilities const& rOther);

    ///@}

};// class ParMmgUtilities
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<PMMGLibrary TPMMGLibrary>
inline std::istream& operator >> (std::istream& rIStream,
                                  ParMmgUtilities<TPMMGLibrary>& rThis);

/// output stream function
template<PMMGLibrary TPMMGLibrary>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ParMmgUtilities<TPMMGLibrary>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.
#endif /* KRATOS_PMMG_UTILITIES defined */
