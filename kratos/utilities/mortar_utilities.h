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
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "includes/enums.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

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

/**
 * @class MortarUtilities
 * @ingroup KratosCore
 * @brief This is a class that provides auxiliar utilities for the mortar integration
 * @details This is a class that provides auxiliar utilities for the mortar integration. Many methods
 * in the following class are templatizaded and with explicit instantations delclared.
 * @note Check the documentation for more details
 * @author Vicente Mataix Ferrandiz
 * Contact: vmataix@cimne.upc.edu
 */
class MortarUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MortarUtilities
    KRATOS_CLASS_POINTER_DEFINITION( MortarUtilities );

    // Some geometrical definitions
    typedef Node<3>                                              NodeType;
    typedef Point                                               PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;

    /// Definition of geometries
    typedef Geometry<NodeType>                               GeometryType;
    typedef Geometry<PointType>                         GeometryPointType;

    /// The integration method type
    typedef GeometryData::IntegrationMethod             IntegrationMethod;

    /// The containers of the components of the model parts
    typedef ModelPart::NodesContainerType                  NodesArrayType;
    typedef ModelPart::ConditionsContainerType        ConditionsArrayType;

    /// Index type definition
    typedef std::size_t                                         IndexType;

    /// Size type definition
    typedef std::size_t                                          SizeType;

    /// A map for integers
    typedef std::unordered_map<IndexType, IndexType>               IntMap;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MortarUtilities(){}

    /// Destructor.
    virtual ~MortarUtilities(){}

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
    ///@name Operations
    ///@{

    /**
     * @brief This functions checks if the length of the line is to short, with the potential of provoque ill condition in the dual LM formulation
     * @param rGeometryLine The line to be checked
     * @param Tolerance The threshold length
     * @return True if the line is too short, false otherwise
     */
    static inline bool LengthCheck(
        const GeometryPointType& rGeometryLine,
        const double Tolerance = 1.0e-6
        ) {
        const double lx = rGeometryLine[0].X() - rGeometryLine[1].X();
        const double ly = rGeometryLine[0].Y() - rGeometryLine[1].Y();

        const double length = std::sqrt(lx * lx + ly * ly);

        return (length < Tolerance) ? true : false;
    }

    /**
     * @brief This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param rGeometryTriangle The triangle to be checked
     * @return True if the triangle is in bad shape, false otherwise
     */
    static inline bool HeronCheck(const GeometryPointType& rGeometryTriangle) {
        return HeronCheck(rGeometryTriangle[0], rGeometryTriangle[1], rGeometryTriangle[2]);
    }

    /**
     * @brief This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param rPointOrig1 The triangle first point
     * @param rPointOrig2 The triangle second point
     * @param rPointOrig3 The triangle third point
     * @return True if the triangle is in bad shape, false otherwise
     */
    static inline bool HeronCheck(
        const PointType& rPointOrig1,
        const PointType& rPointOrig2,
        const PointType& rPointOrig3
        )
    {
        const double a = MathUtils<double>::Norm3(rPointOrig1.Coordinates()-rPointOrig2.Coordinates());
        const double b = MathUtils<double>::Norm3(rPointOrig2.Coordinates()-rPointOrig3.Coordinates());
        const double c = MathUtils<double>::Norm3(rPointOrig3.Coordinates()-rPointOrig1.Coordinates());

        const double s = 0.5 * (a + b + c);
        const double A2 = s * (s - a) * (s - b) * (s - c);

        const bool Check = A2 <= 0.0 ? true : false;  // We consider as bad shaped the ones with no area or negative A2 (semiperimeter smaller than any side)

//         // Debug
//         KRATOS_INFO("Check") << Check << " A2: " << A2 << std::endl;
//         if (Check == true) {
//             KRATOS_WARNING("Bad shape") << "Warning:: The triangle is in bad shape" << std::endl;
//             KRATOS_INFO("Mathematica triangle") << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << rPointOrig1.X() << "," << rPointOrig1.Y() << "," << rPointOrig1.Z()  << "},{" << rPointOrig2.X() << "," << rPointOrig2.Y() << "," << rPointOrig2.Z()  << "},{" << rPointOrig3.X() << "," << rPointOrig3.Y() << "," << rPointOrig3.Z()  << "}}]}]" << std::endl;
//         }

        return Check;
    }

    /**
     * @brief This function rotates to align the projected points to a parallel plane to XY
     * @param rPointToRotate The points from the origin geometry and the the point rotated
     * @param rPointReferenceRotation The center point used as reference to rotate
     * @param rSlaveTangentXi The first tangent vector of the slave condition
     * @param rSlaveTangentEta The second tangent vector of the slave condition
     * @param Inversed If we rotate to the XY or we recover from XY
     */
    static inline void RotatePoint(
        PointType& rPointToRotate,
        const PointType& rPointReferenceRotation,
        const array_1d<double, 3>& rSlaveTangentXi,
        const array_1d<double, 3>& rSlaveTangentEta,
        const bool Inversed
        )
    {
        // We move to the (0,0,0)
        PointType aux_point_to_rotate;
        aux_point_to_rotate.Coordinates() = rPointToRotate.Coordinates() - rPointReferenceRotation.Coordinates();

        BoundedMatrix<double, 3, 3> rotation_matrix = ZeroMatrix(3, 3);

        if (Inversed == false) {
            for (IndexType i = 0; i < 3; ++i) {
                rotation_matrix(0, i) = rSlaveTangentXi[i];
                rotation_matrix(1, i) = rSlaveTangentEta[i];
            }
        } else {
            for (IndexType i = 0; i < 3; ++i) {
                rotation_matrix(i, 0) = rSlaveTangentXi[i];
                rotation_matrix(i, 1) = rSlaveTangentEta[i];
            }
        }

        rPointToRotate.Coordinates() = prod(rotation_matrix, aux_point_to_rotate) + rPointReferenceRotation.Coordinates();
    }

    /**
     * @brief This function calculates the r_normal in a specific GP with a given shape function
     * @param rN The shape function considered
     * @param rGeometry The geometry of condition of interest
     * @return The r_normal in the GP
     */
    static inline array_1d<double,3> GaussPointUnitNormal(
        const Vector& rN,
        const GeometryType& rGeometry
        ) {
        array_1d<double,3> r_normal = ZeroVector(3);
        for( IndexType i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node )
            r_normal += rN[i_node] * rGeometry[i_node].FastGetSolutionStepValue(NORMAL);

        const double this_norm = norm_2(r_normal);

        KRATOS_DEBUG_ERROR_IF(this_norm < std::numeric_limits<double>::epsilon()) << "Zero norm r_normal vector. Norm:" << this_norm << std::endl;

        r_normal /= this_norm;

        return r_normal;
    }

    /**
     * @brief This function gives you the indexes needed to order a vector
     * @param rThisVector The vector to order
     * @return idx The vector of indexes
     */
    template <typename TType>
    static std::vector<std::size_t> SortIndexes(const std::vector<TType> &rThisVector) {
        // Initialize original index locations
        std::vector<std::size_t> idx(rThisVector.size());
        iota(idx.begin(), idx.end(), 0);

        // Sort indexes based on comparing values in rThisVector
        std::sort(idx.begin(), idx.end(),
            [&rThisVector](std::size_t i1, std::size_t i2) {return rThisVector[i1] < rThisVector[i2];});

        return idx;
    }

    /**
     * @brief It computes the mean of the r_normal in the condition in all the nodes
     * @param rModelPart The model part to compute
     */
    static inline void ComputeNodesMeanNormalModelPart(ModelPart& rModelPart) {
        // Check NORMAL is available
        KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(NORMAL)) << "NORMAL is not available on the solution step data variable database" << std::endl;
        
        // We iterate over nodes
        NodesArrayType& r_nodes_array = rModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();
        const int num_nodes = static_cast<int>(r_nodes_array.size());

        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        // Reset NORMAL
        VariableUtils().SetVectorVar(NORMAL, zero_array, r_nodes_array);

        // Sum all the nodes normals
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
        const auto it_cond_begin = r_conditions_array.begin();

        // Declare auxiliar coordinates
        CoordinatesArrayType aux_coords;

        #pragma omp parallel for firstprivate(aux_coords)
        for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;
            const GeometryType& r_geometry = it_cond->GetGeometry();

            // Set condition normal
            r_geometry.PointLocalCoordinates(aux_coords, r_geometry.Center());
            it_cond->SetValue(NORMAL, r_geometry.UnitNormal(aux_coords));
        }

        // Adding the normal contribution of each node
        for(Condition& r_cond : r_conditions_array) {
            GeometryType& r_geometry = r_cond.GetGeometry();

            // Iterate over nodes
            for (NodeType& r_node : r_geometry) {
                r_geometry.PointLocalCoordinates(aux_coords, r_node.Coordinates());
                noalias(r_node.FastGetSolutionStepValue(NORMAL)) += r_geometry.UnitNormal(aux_coords);
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            array_1d<double, 3>& r_normal = it_node->FastGetSolutionStepValue(NORMAL);
            const double norm_normal = norm_2(r_normal);

            if (norm_normal > std::numeric_limits<double>::epsilon()) r_normal /= norm_normal;
            else KRATOS_ERROR_IF(it_node->Is(INTERFACE)) << "ERROR:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
        }
    }

    /**
     * @brief It inverts the order of the nodes in the conditions of a model part in order to invert the r_normal
     * @param rContainer reference to the objective container
     */
    template<class TContainerType>
    static inline void InvertNormal(TContainerType& rContainer) {
        const auto it_cont_begin = rContainer.begin();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rContainer.size()); ++i) {
            auto it_cont = it_cont_begin + i;
            GeometryType& r_geometry = it_cont->GetGeometry();

            auto& data_geom = r_geometry.GetContainer();
            std::reverse(data_geom.begin(), data_geom.end());
        }
    }

    /**
     * @brief It calculates the matrix of coordinates of a geometry
     * @param rGeometry The geometry to calculate
     * @param Current If we calculate the Current coordinates or the initial ones
     * @param Step The time step where it is computed
     * @return coordinates The matrix containing the coordinates of the geometry
     */

    template< SizeType TDim, SizeType TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TDim> GetCoordinates(
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
     * @brief It calculates the matrix containing the tangent vector of the LM (for frictional contact)
     * @param rGeometry The geometry to calculate
     * @return tangent_matrix The matrix containing the tangent vectors of the LM
     */
    template< SizeType TNumNodes, SizeType TDim>
    static inline BoundedMatrix<double, TNumNodes, TDim> ComputeTangentMatrix(const GeometryType& rGeometry) {
        /* DEFINITIONS */
        // Zero tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();
        // Tangent matrix
        BoundedMatrix<double, TNumNodes, TDim> tangent_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& r_lm = rGeometry[i_node].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
            if (norm_2(r_lm) > zero_tolerance) { // Non zero LM
                const array_1d<double, 3>& r_normal = rGeometry[i_node].FastGetSolutionStepValue(NORMAL);
                const array_1d<double, 3> tangent_lm = r_lm - inner_prod(r_lm, r_normal) * r_normal;
                if (norm_2(tangent_lm) > zero_tolerance) {
                    const array_1d<double, 3> tangent = tangent_lm/norm_2(tangent_lm);
                    for (std::size_t i_dof = 0; i_dof < TDim; ++i_dof)
                        tangent_matrix(i_node, i_dof) = tangent[i_dof];
                } else {
                    for (std::size_t i_dof = 0; i_dof < TDim; ++i_dof)
                        tangent_matrix(i_node, i_dof) = 0.0;
                }
            } else { // In case of zero LM
                for (std::size_t i_dof = 0; i_dof < TDim; ++i_dof)
                    tangent_matrix(i_node, i_dof) = 0.0;
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
    static inline array_1d<double, TNumNodes> GetVariableVector(
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
    static inline BoundedMatrix<double, TNumNodes, 1> GetVariableVectorMatrix(
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
    static inline array_1d<double, TNumNodes> GetVariableVector(
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
    static inline BoundedMatrix<double, TNumNodes, 1> GetVariableVectorMatrix(
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
    static inline BoundedMatrix<double, TNumNodes, TDim> GetVariableMatrix(
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
    static inline BoundedMatrix<double, TNumNodes, TDim> GetVariableMatrix(
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
    static inline BoundedMatrix<double, TNumNodes, TDim> GetAbsMatrix(const BoundedMatrix<double, TNumNodes, TDim>& rInputMatrix) {
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
    static inline unsigned int SizeToCompute()
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
    template< class TVarType, HistoricalValues THist>
    static inline void ResetValue(
        ModelPart& rThisModelPart,
        TVarType& rThisVariable
        );

    /**
     * @brief This method resets the auxiliar value
     * @param rThisModelPart The model part to update
     */
    template< class TVarType>
    static inline void ResetAuxiliarValue(ModelPart& rThisModelPart);

    /**
     * @brief This method returns the auxiliar variable
     */
    template< class TVarType>
    static inline TVarType GetAuxiliarVariable();

    /**
     * @brief This method returns the auxiliar variable
     */
    template< class TVarType>
    static inline double GetAuxiliarValue(
        NodeType::Pointer pThisNode,
        unsigned int iSize
        );

    /**
     * @brief This method adds the value
     * @param rThisGeometry The geometrty to update
     * @param rThisVariable The variable to set
     * @param rThisValue The matrix to be updated
     */
    template< class TVarType, HistoricalValues THist>
    static inline void MatrixValue(
        GeometryType& rThisGeometry,
        TVarType& rThisVariable,
        Matrix& rThisValue
        );

    /**
     * @brief This method adds the value
     * @warning This operation is not threadsafe
     * @param rThisGeometry The geometrty to update
     * @param rThisVariable The variable to set
     * @param rThisValue The matrix to be updated
     */
    template< class TVarType, HistoricalValues THist>
    static inline void AddValue(
        GeometryType& rThisGeometry,
        TVarType& rThisVariable,
        const Matrix& rThisValue
        );

    /**
     * @brief This method adds the value
     * @param pThisNode The node to update
     * @param rThisVariable The variable to set
     */
    template< class TVarType, HistoricalValues THist>
    static inline void AddAreaWeightedNodalValue(
        NodeType::Pointer pThisNode,
        TVarType& rThisVariable,
        const double RefArea = 1.0,
        const double Tolerance = 1.0e-4
        );

    /**
     * @brief This method updates the database in the amster side
     * @param rThisModelPart The model part
     * @param rThisVariable The variable to set
     * @param Dx The vector with the increment of the value
     * @param Index The index used in the  case of a vector variable
     * @param ConectivityDatabase The database that will be used to assemble the system
     */
    template< class TVarType, HistoricalValues THist>
    static inline void UpdateDatabase(
        ModelPart& rThisModelPart,
        TVarType& rThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        );

private:
};// class MortarUtilities

///@name Explicit Specializations
///@{

template<>
inline void MortarUtilities::ResetValue<Variable<double>, Historical>(
        ModelPart& rThisModelPart,
        Variable<double>& rThisVariable
        )
{
    NodesArrayType& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetScalarVar(rThisVariable, 0.0, r_nodes_array);
}

template<>
inline void MortarUtilities::ResetValue<Variable<array_1d<double, 3>>, Historical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& rThisVariable
        ) {
    NodesArrayType& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetVectorVar(rThisVariable, ZeroVector(3), r_nodes_array);
}

template<>
inline void MortarUtilities::ResetValue<Variable<double>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<double>& rThisVariable
        ) {
    NodesArrayType& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(rThisVariable, 0.0, r_nodes_array);
}

template<>
inline void MortarUtilities::ResetValue<Variable<array_1d<double, 3>>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& rThisVariable
        ) {
    const array_1d<double, 3> zero_array = ZeroVector(3);
    NodesArrayType& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(rThisVariable, zero_array, r_nodes_array);
}

template<>
inline void MortarUtilities::ResetAuxiliarValue<Variable<double>>(ModelPart& rThisModelPart) {
    NodesArrayType& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(NODAL_MAUX, 0.0, r_nodes_array);
}

template<>
inline void MortarUtilities::ResetAuxiliarValue<Variable<array_1d<double, 3>>>(ModelPart& rThisModelPart) {
    const array_1d<double, 3> zero_array = ZeroVector(3);
    NodesArrayType& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(NODAL_VAUX, zero_array, r_nodes_array);
}

template< >
inline Variable<double> MortarUtilities::GetAuxiliarVariable<Variable<double>>() {
    return NODAL_MAUX;
}

template< >
inline Variable<array_1d<double, 3>> MortarUtilities::GetAuxiliarVariable<Variable<array_1d<double, 3>>>() {
    return NODAL_VAUX;
}

template< >
inline double MortarUtilities::GetAuxiliarValue<Variable<double>>(
    Node<3>::Pointer pThisNode,
    unsigned int iSize
    ) {
    return pThisNode->GetValue(NODAL_MAUX);
}

template< >
inline double MortarUtilities::GetAuxiliarValue<Variable<array_1d<double, 3>>>(
    Node<3>::Pointer pThisNode,
    unsigned int iSize
    ) {
    switch ( iSize ) {
        case 0:
            return pThisNode->GetValue(NODAL_VAUX_X);
        case 1:
            return pThisNode->GetValue(NODAL_VAUX_Y);
        case 2:
            return pThisNode->GetValue(NODAL_VAUX_Z);
        default:
            return 0.0;
    }

    return 0.0;
}

template<>
inline void MortarUtilities::MatrixValue<Variable<double>, Historical>(
        GeometryType& rThisGeometry,
        Variable<double>& rThisVariable,
        Matrix& rThisValue
        ) {
    if (rThisValue.size1() != rThisGeometry.size() || rThisValue.size2() != 1)
        rThisValue.resize(rThisGeometry.size(), 1, false);

    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node)
        rThisValue(i_node, 0) = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
}

template<>
inline void MortarUtilities::MatrixValue<Variable<array_1d<double, 3>>, Historical>(
        GeometryType& rThisGeometry,
        Variable<array_1d<double, 3>>& rThisVariable,
        Matrix& rThisValue
        ) {
    const std::size_t num_nodes = rThisGeometry.size();
    const std::size_t dimension = rThisGeometry.WorkingSpaceDimension();
    if (rThisValue.size1() != num_nodes || rThisValue.size2() != dimension)
        rThisValue.resize(num_nodes, dimension, false);

    for (IndexType i_node = 0; i_node < num_nodes; ++i_node) {
        const auto& rvalue = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        for (IndexType i_dim = 0; i_dim < dimension; ++i_dim)
            rThisValue(i_node, i_dim) = rvalue[i_dim];
    }
}
template<>
inline void MortarUtilities::MatrixValue<Variable<double>, NonHistorical>(
        GeometryType& rThisGeometry,
        Variable<double>& rThisVariable,
        Matrix& rThisValue
        ) {
    if (rThisValue.size1() != rThisGeometry.size() || rThisValue.size2() != 1)
        rThisValue.resize(rThisGeometry.size(), 1, false);

    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node)
        rThisValue(i_node, 0) = rThisGeometry[i_node].GetValue(rThisVariable);
}

template<>
inline void MortarUtilities::MatrixValue<Variable<array_1d<double, 3>>, NonHistorical>(
        GeometryType& rThisGeometry,
        Variable<array_1d<double, 3>>& rThisVariable,
        Matrix& rThisValue
        ) {
    const std::size_t num_nodes = rThisGeometry.size();
    const std::size_t dimension = rThisGeometry.WorkingSpaceDimension();
    if (rThisValue.size1() != num_nodes || rThisValue.size2() != dimension)
        rThisValue.resize(num_nodes, dimension, false);

    for (IndexType i_node = 0; i_node < num_nodes; ++i_node) {
        const auto& rvalue = rThisGeometry[i_node].GetValue(rThisVariable);
        for (IndexType i_dim = 0; i_dim < dimension; ++i_dim)
            rThisValue(i_node, i_dim) = rvalue[i_dim];
    }
}

template<>
inline void MortarUtilities::AddValue<Variable<double>, Historical>(
        GeometryType& rThisGeometry,
        Variable<double>& rThisVariable,
        const Matrix& rThisValue
        ) {
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        double& aux_value = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        #pragma omp atomic
        aux_value += rThisValue(i_node, 0);
    }
}

template<>
inline void MortarUtilities::AddValue<Variable<array_1d<double, 3>>, Historical>(
        GeometryType& rThisGeometry,
        Variable<array_1d<double, 3>>& rThisVariable,
        const Matrix& rThisValue
        ) {
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        auto& aux_vector = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        for (unsigned int i_dim = 0; i_dim < rThisGeometry.WorkingSpaceDimension(); ++i_dim) {
            double& aux_value = aux_vector[i_dim];
            #pragma omp atomic
            aux_value += rThisValue(i_node, i_dim);
        }
    }
}
template<>
inline void MortarUtilities::AddValue<Variable<double>, NonHistorical>(
        GeometryType& rThisGeometry,
        Variable<double>& rThisVariable,
        const Matrix& rThisValue
        ) {
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        double& aux_value = rThisGeometry[i_node].GetValue(rThisVariable);
        #pragma omp atomic
        aux_value += rThisValue(i_node, 0);
    }
}

template<>
inline void MortarUtilities::AddValue<Variable<array_1d<double, 3>>, NonHistorical>(
        GeometryType& rThisGeometry,
        Variable<array_1d<double, 3>>& rThisVariable,
        const Matrix& rThisValue
        ) {
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        auto& aux_vector = rThisGeometry[i_node].GetValue(rThisVariable);
        for (unsigned int i_dim = 0; i_dim < rThisGeometry.WorkingSpaceDimension(); ++i_dim) {
            double& aux_value = aux_vector[i_dim];
            #pragma omp atomic
            aux_value += rThisValue(i_node, i_dim);
        }
    }
}

template<>
inline void MortarUtilities::AddAreaWeightedNodalValue<Variable<double>, Historical>(
        Node<3>::Pointer pThisNode,
        Variable<double>& rThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG
    if (null_area) KRATOS_WARNING("WARNING:: NODE OF NULL AREA.") << " ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    double& aux_value = pThisNode->FastGetSolutionStepValue(rThisVariable);
    #pragma omp atomic
    aux_value += area_coeff * pThisNode->GetValue(NODAL_MAUX);
}

template<>
inline void MortarUtilities::AddAreaWeightedNodalValue<Variable<array_1d<double, 3>>, Historical>(
        Node<3>::Pointer pThisNode,
        Variable<array_1d<double, 3>>& rThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG
    if (null_area) KRATOS_WARNING("WARNING:: NODE OF NULL AREA.") << " ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    auto& aux_vector = pThisNode->FastGetSolutionStepValue(rThisVariable);
    const auto& nodal_vaux = pThisNode->GetValue(NODAL_VAUX);
    for (IndexType i = 0; i < 3; ++i) {
        double& aux_value = aux_vector[i];
        #pragma omp atomic
        aux_value += area_coeff * nodal_vaux[i];
    }
}

template<>
inline void MortarUtilities::AddAreaWeightedNodalValue<Variable<double>, NonHistorical>(
        Node<3>::Pointer pThisNode,
        Variable<double>& rThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG
    if (null_area) KRATOS_WARNING("WARNING:: NODE OF NULL AREA.") << " ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    double& aux_value = pThisNode->GetValue(rThisVariable);
    #pragma omp atomic
    aux_value += area_coeff * pThisNode->GetValue(NODAL_MAUX);
}

template<>
inline void MortarUtilities::AddAreaWeightedNodalValue<Variable<array_1d<double, 3>>, NonHistorical>(
        Node<3>::Pointer pThisNode,
        Variable<array_1d<double, 3>>& rThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG
    if (null_area) KRATOS_WARNING("WARNING:: NODE OF NULL AREA.") << " ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    auto& aux_vector = pThisNode->GetValue(rThisVariable);
    const auto& nodal_vaux = pThisNode->GetValue(NODAL_VAUX);
    for (IndexType i = 0; i < 3; ++i) {
        double& aux_value = aux_vector[i];
        #pragma omp atomic
        aux_value += area_coeff * nodal_vaux[i];
    }
}

template<>
inline void MortarUtilities::UpdateDatabase<Variable<double>, Historical>(
        ModelPart& rThisModelPart,
        Variable<double>& rThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        p_node->FastGetSolutionStepValue(rThisVariable) += Dx[i];
    }
}

template<>
inline void MortarUtilities::UpdateDatabase<Variable<array_1d<double, 3>>, Historical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& rThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        auto& value = p_node->FastGetSolutionStepValue(rThisVariable);
        value[Index] += Dx[i];
    }
}
template<>
inline void MortarUtilities::UpdateDatabase<Variable<double>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<double>& rThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        p_node->GetValue(rThisVariable) += Dx[i];
    }
}

template<>
inline void MortarUtilities::UpdateDatabase<Variable<array_1d<double, 3>>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& rThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        auto& value = p_node->GetValue(rThisVariable);
        value[Index] += Dx[i];
    }
}

}
#endif /* KRATOS_MORTAR_UTILITIES defined */
