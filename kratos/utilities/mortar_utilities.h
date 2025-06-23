//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MORTAR_UTILITIES)
#define KRATOS_MORTAR_UTILITIES

// System includes
#include <numeric>
#include <unordered_map>

// External includes

// Project includes
#include "includes/variables.h"
#include "includes/node.h"
#include "geometries/geometry.h"

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

class ModelPart; // forward-declaring to not having to include it here

/**
 * @brief This struct is used in order to identify when using the historical and non historical variables
 */
struct MortarUtilitiesSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};

/**
 * @namespace MortarUtilities
 * @ingroup KratosCore
 * @brief This is a class that provides auxiliar utilities for the mortar integration
 * @details This is a class that provides auxiliar utilities for the mortar integration. Many methods
 * in the following class are templatizaded and with explicit instantations delclared.
 * @note Check the documentation for more details
 * @author Vicente Mataix Ferrandiz
 * Contact: vmataix@cimne.upc.edu
 */
namespace MortarUtilities
{
    ///@name Type Definitions
    ///@{

    // Some geometrical definitions
    typedef Node                                              NodeType;
    typedef Point                                               PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;

    /// Definition of geometries
    typedef Geometry<NodeType>                               GeometryType;
    typedef Geometry<PointType>                         GeometryPointType;

    /// Index type definition
    typedef std::size_t                                         IndexType;

    /// Size type definition
    typedef std::size_t                                          SizeType;

    /// A map for integers
    typedef std::unordered_map<IndexType, IndexType>               IntMap;

    ///@}
    ///@name  Functions
    ///@{

    /**
     * @brief This functions checks if the length of the line is to short, with the potential of provoque ill condition in the dual LM formulation
     * @param rGeometryLine The line to be checked
     * @param Tolerance The threshold length
     * @return True if the line is too short, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) LengthCheck(
        const GeometryPointType& rGeometryLine,
        const double Tolerance = 1.0e-6
        );

    /**
     * @brief This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param rGeometryTriangle The triangle to be checked
     * @return True if the triangle is in bad shape, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) HeronCheck(const GeometryPointType& rGeometryTriangle);

    /**
     * @brief This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param rPointOrig1 The triangle first point
     * @param rPointOrig2 The triangle second point
     * @param rPointOrig3 The triangle third point
     * @return True if the triangle is in bad shape, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) HeronCheck(
        const PointType& rPointOrig1,
        const PointType& rPointOrig2,
        const PointType& rPointOrig3
        );

    /**
     * @brief This function rotates to align the projected points to a parallel plane to XY
     * @param rPointToRotate The points from the origin geometry and the the point rotated
     * @param rPointReferenceRotation The center point used as reference to rotate
     * @param rSlaveTangentXi The first tangent vector of the slave condition
     * @param rSlaveTangentEta The second tangent vector of the slave condition
     * @param Inversed If we rotate to the XY or we recover from XY
     */
    void KRATOS_API(KRATOS_CORE) RotatePoint(
        PointType& rPointToRotate,
        const PointType& rPointReferenceRotation,
        const array_1d<double, 3>& rSlaveTangentXi,
        const array_1d<double, 3>& rSlaveTangentEta,
        const bool Inversed
        );

    /**
     * @brief This function calculates the r_normal in a specific GP with a given shape function
     * @param rN The shape function considered
     * @param rGeometry The geometry of condition of interest
     * @return The r_normal in the GP
     */
    array_1d<double,3> KRATOS_API(KRATOS_CORE) GaussPointUnitNormal(
        const Vector& rN,
        const GeometryType& rGeometry
        );

    /**
     * @brief This function gives you the indexes needed to order a vector
     * @param rThisVector The vector to order
     * @return idx The vector of indexes
     */
    template <typename TType>
    std::vector<std::size_t> SortIndexes(const std::vector<TType> &rThisVector) {
        // Initialize original index locations
        std::vector<std::size_t> idx(rThisVector.size());
        iota(idx.begin(), idx.end(), 0);

        // Sort indexes based on comparing values in rThisVector
        std::sort(idx.begin(), idx.end(),
            [&rThisVector](std::size_t i1, std::size_t i2) {return rThisVector[i1] < rThisVector[i2];});

        return idx;
    }

    /**
     * @brief It computes the mean of the normal in the condition in all the nodes
     * @param rModelPart The model part to compute
     * @param ComputeConditions If computed over conditions or elements
     */
    void KRATOS_API(KRATOS_CORE) ComputeNodesMeanNormalModelPart(
        ModelPart& rModelPart,
        const bool ComputeConditions = true
        );

    /**
     * @brief It computes the tangent in all the nodes of the model part
     * @param rModelPart The model part to compute
     * @param pSlipVariable The pointer to the slip variable
     * @param SlipCoefficient The slip contribution
     * @param SlipAlways Uses the slip even in case that LM are available
     */
    void KRATOS_API(KRATOS_CORE) ComputeNodesTangentModelPart(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>* pSlipVariable = NULL,
        const double SlipCoefficient = 1.0,
        const bool SlipAlways = false
        );

    /**
     * @brief It computes the tangent in all the nodes of the model part from its normal
     * @param rModelPart The model part to compute
     */
    void KRATOS_API(KRATOS_CORE) ComputeNodesTangentFromNormalModelPart(ModelPart& rModelPart);

    /**
     * @brief It computes the tangent on the given node using the normal provided
     * @param rNode The node where to compute the tangent
     * @param rNormal The normal vector
     * @param Dimension The current working dimension
     */
    void KRATOS_API(KRATOS_CORE) ComputeTangentsFromNormal(
        NodeType& rNode,
        const array_1d<double, 3>& rNormal,
        const std::size_t Dimension = 3
        );

    /**
     * @brief It computes the tangent on the given node using the LM direction and Slip direction
     * @param rNode The node where to compute the tangent
     * @param StepLM The considered step slip
     * @param pSlipVariable The pointer to the slip variable
     * @param SlipCoefficient The slip contribution
     * @param Dimension The current working dimension
     */
    void KRATOS_API(KRATOS_CORE) ComputeTangentNodeWithLMAndSlip(
        NodeType& rNode,
        const std::size_t StepLM = 0,
        const Variable<array_1d<double, 3>>* pSlipVariable = NULL,
        const double SlipCoefficient = 1.0,
        const std::size_t Dimension = 3
        );

    /**
     * @brief It computes the tangent on the given node using the Slip direction
     * @param rNode The node where to compute the tangent
     * @param StepLM The considered step slip
     * @param pSlipVariable The pointer to the slip variable
     * @param SlipCoefficient The slip contribution
     * @param Dimension The current working dimension
     */
    void KRATOS_API(KRATOS_CORE) ComputeTangentNodeWithSlip(
        NodeType& rNode,
        const std::size_t StepLM = 0,
        const Variable<array_1d<double, 3>>* pSlipVariable = NULL,
        const double SlipCoefficient = 1.0,
        const std::size_t Dimension = 3
        );

    /**
     * @brief It inverts the order of the nodes in the conditions of a model part in order to invert the normal when certain flag is active
     * @param rContainer Reference to the objective container
     * @param Flag The flag of the entities inverted
     */
    template<class TContainerType>
    void InvertNormalForFlag(
        TContainerType& rContainer,
        const Flags Flag
        )
    {
        bool to_invert = false;
        const auto it_cont_begin = rContainer.begin();
        #pragma omp parallel for firstprivate(to_invert)
        for(int i = 0; i < static_cast<int>(rContainer.size()); ++i) {
            auto it_cont = it_cont_begin + i;
            to_invert = Flag == Flags() ? true : it_cont->IsDefined(Flag) ? it_cont->Is(Flag) : false;

            if (to_invert) {
                GeometryType& r_geometry = it_cont->GetGeometry();

                auto& data_geom = r_geometry.GetContainer();
                std::reverse(data_geom.begin(), data_geom.end());
            }
        }
    }

    /**
     * @brief It inverts the order of the nodes in the conditions of a model part in order to invert the normal
     * @param rContainer Reference to the objective container
     */
    template<class TContainerType>
    void InvertNormal(TContainerType& rContainer)
    {
        InvertNormalForFlag(rContainer, Flags());
    }

    /**
     * @brief It calculates the matrix of coordinates of a geometry
     * @param rGeometry The geometry to calculate
     * @param Current If we calculate the Current coordinates or the initial ones
     * @param Step The time step where it is computed
     * @return coordinates The matrix containing the coordinates of the geometry
     */

    template< SizeType TDim, SizeType TNumNodes>
    BoundedMatrix<double, TNumNodes, TDim> GetCoordinates(
        const GeometryType& rGeometry,
        const bool Current = true,
        const IndexType Step = 0
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, TDim> coordinates;
        array_1d<double, 3> coord;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            if (Current) {
                coord = rGeometry[i_node].Coordinates();
            } else {
                coord = rGeometry[i_node].GetInitialPosition();

                if (Step > 0)
                    coord += rGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, Step);
            }

            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof)
                coordinates(i_node, i_dof) = coord[i_dof];
        }

        return coordinates;
    }

    /**
     * @brief It calculates the matrix containing the tangent vector TANGENT_XI
     * @param rGeometry The geometry to calculate
     * @return tangent_matrix The matrix containing the tangent vectors of the LM
     */
    template< SizeType TNumNodes, SizeType TDim>
    BoundedMatrix<double, TNumNodes, TDim> ComputeTangentMatrix(const GeometryType& rGeometry)
    {
        // Tangent matrix
        BoundedMatrix<double, TNumNodes, TDim> tangent_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const auto& r_node = rGeometry[i_node];
            const auto& r_tangent = r_node.GetValue(TANGENT_XI);
            for (std::size_t i_dof = 0; i_dof < TDim; ++i_dof) {
                tangent_matrix(i_node, i_dof) = r_tangent[i_dof];
            }
        }

        return tangent_matrix;
    }

    /**
     * @brief It calculates the vector of an historical variable of a geometry
     * @param rGeometry The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param Step The step where it is computed
     * @return var_vector The vector containing the variables of the geometry
     */
    template< SizeType TNumNodes, class TVarType = Variable<double>>
    array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& rGeometry,
        const TVarType& rVariable,
        const IndexType Step
        ) {
        /* DEFINITIONS */
        array_1d<double, TNumNodes> var_vector;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector[i_node] = rGeometry[i_node].FastGetSolutionStepValue(rVariable, Step);

        return var_vector;
    }

    /**
     * @brief It calculates the vector of an historical variable of a geometry
     * @param rGeometry The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param Step The step where it is computed
     * @return var_vector The vector containing the variables of the geometry
     */
    template< SizeType TNumNodes, class TVarType = Variable<double> >
    BoundedMatrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& rGeometry,
        const TVarType& rVariable,
        const unsigned int Step
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, 1> var_vector;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector(i_node, 0) = rGeometry[i_node].FastGetSolutionStepValue(rVariable, Step);

        return var_vector;
    }

    /**
     * @brief It calculates the vector of a non-historical variable of a geometry
     * @param rGeometry The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_vector The vector containing the variables of the geometry
     */
    template< SizeType TNumNodes, class TVarType = Variable<double> >
    array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& rGeometry,
        const TVarType& rVariable
        ) {
        /* DEFINITIONS */
        array_1d<double, TNumNodes> var_vector;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector[i_node] = rGeometry[i_node].GetValue(rVariable);

        return var_vector;
    }

    /**
     * @brief It calculates the vector of a non-historical variable of a geometry
     * @param rGeometry The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_vector The vector containing the variables of the geometry
     */
    template< SizeType TNumNodes, class TVarType = Variable<double> >
    BoundedMatrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& rGeometry,
        const TVarType& rVariable
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, 1> var_vector;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector(i_node, 0) = rGeometry[i_node].GetValue(rVariable);

        return var_vector;
    }

    /**
     * @brief It calculates the matrix of a variable of a geometry
     * @param rGeometry The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param Step The step where it is computed
     * @return var_matrix The matrix containing the variables of the geometry
     */
    template< SizeType TDim, SizeType TNumNodes>
    BoundedMatrix<double, TNumNodes, TDim> GetVariableMatrix(
        const GeometryType& rGeometry,
        const Variable<array_1d<double,3> >& rVariable,
        const unsigned int Step
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, TDim> var_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& value = rGeometry[i_node].FastGetSolutionStepValue(rVariable, Step);
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof)
                var_matrix(i_node, i_dof) = value[i_dof];
        }

        return var_matrix;
    }

    /**
     * @brief It calculates the matrix of a non-historical variable of a geometry
     * @param rGeometry The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_matrix The matrix containing the variables of the geometry
     */
    template< SizeType TDim, SizeType TNumNodes>
    BoundedMatrix<double, TNumNodes, TDim> GetVariableMatrix(
        const GeometryType& rGeometry,
        const Variable<array_1d<double,3> >& rVariable
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, TDim> var_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& value = rGeometry[i_node].GetValue(rVariable);
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof)
                var_matrix(i_node, i_dof) = value[i_dof];
        }

        return var_matrix;
    }

    /**
     * @brief It calculates the matrix containing the absolute value of another matrix
     * @param rInputMatrix The original matrix
     * @return AbsMatrix The matrix containing the absolute value of another matrix
     */
    template< SizeType TDim, SizeType TNumNodes>
    BoundedMatrix<double, TNumNodes, TDim> GetAbsMatrix(const BoundedMatrix<double, TNumNodes, TDim>& rInputMatrix) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, TDim> AbsMatrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof)
                AbsMatrix(i_node, i_dof) = std::abs(rInputMatrix(i_node, i_dof));
        }

        return AbsMatrix;
    }

    /**
     * @brief This method gives the size to be computed
     */
    template< SizeType TDim, class TVarType>
    unsigned int SizeToCompute()
    {
       if (typeid(TVarType) == typeid(Variable<array_1d<double, 3>>))
           return TDim;

       return 1;
    }

    /**
     * @brief This method resets the value
     * @param rThisModelPart The model part to update
     * @param rThisVariable The variable to set
     */
    template< class TVarType, bool THistorical>
    void KRATOS_API(KRATOS_CORE) ResetValue(
        ModelPart& rThisModelPart,
        const TVarType& rThisVariable
        );

    /**
     * @brief This method resets the auxiliar value
     * @param rThisModelPart The model part to update
     */
    template< class TVarType>
    void KRATOS_API(KRATOS_CORE) ResetAuxiliarValue(ModelPart& rThisModelPart);

    /**
     * @brief This method returns the auxiliar variable
     * @return The auxiliar variable
     */
    template< class TVarType>
    const std::string KRATOS_API(KRATOS_CORE) GetAuxiliarVariable();

    /**
     * @brief This method returns the auxiliar variable
     * @param rThisNode Reference to the node of interest
     * @param iSize The Index of the component
     * @return The value of the auxiliar variable
     */
    template< class TVarType>
    double KRATOS_API(KRATOS_CORE) GetAuxiliarValue(
        NodeType& rThisNode,
        const std::size_t iSize
        );

    /**
     * @brief This method adds the value
     * @param rThisGeometry The geometrty to update
     * @param rThisVariable The variable to set
     * @param rThisValue The matrix to be updated
     */
    template< class TVarType, bool THistorical>
    void KRATOS_API(KRATOS_CORE) MatrixValue(
        const GeometryType& rThisGeometry,
        const TVarType& rThisVariable,
        Matrix& rThisValue
        );

    /**
     * @brief This method adds the value
     * @warning This operation is not threadsafe
     * @param rThisGeometry The geometrty to update
     * @param rThisVariable The variable to set
     * @param rThisValue The matrix to be updated
     */
    template< class TVarType, bool THistorical>
    void KRATOS_API(KRATOS_CORE) AddValue(
        GeometryType& rThisGeometry,
        const TVarType& rThisVariable,
        const Matrix& rThisValue
        );

    /**
     * @brief This method adds the value
     * @param rThisNode The node to update
     * @param rThisVariable The variable to set
     */
    template< class TVarType, bool THistorical>
    void KRATOS_API(KRATOS_CORE) AddAreaWeightedNodalValue(
        NodeType& rThisNode,
        const TVarType& rThisVariable,
        const double RefArea = 1.0,
        const double Tolerance = 1.0e-4
        );

    /**
     * @brief This method updates the database in the amster side
     * @param rThisModelPart The model part
     * @param rThisVariable The variable to set
     * @param rDx The vector with the increment of the value
     * @param Index The index used in the  case of a vector variable
     * @param rConectivityDatabase The database that will be used to assemble the system
     */
    template< class TVarType, bool THistorical>
    void KRATOS_API(KRATOS_CORE) UpdateDatabase(
        ModelPart& rThisModelPart,
        const TVarType& rThisVariable,
        Vector& rDx,
        const std::size_t Index,
        IntMap& rConectivityDatabase
        );
};// namespace MortarUtilities
} // namespace Kratos
#endif /* KRATOS_MORTAR_UTILITIES defined */
