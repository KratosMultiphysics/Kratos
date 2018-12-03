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
     * @brief Project a point over a line/plane following an arbitrary direction
     * @param Geom The geometry where to be projected
     * @param PointDestiny The point to be projected
     * @param PointProjected The point pojected over the plane
     * @param Normal The normal of the geometry
     * @param Vector The direction to project
     * @return Distance The distance between surfaces
     */

    KRATOS_DEPRECATED_MESSAGE("Method moved to geometrical_projection_utilities.h. Please update your declaration") static inline double FastProjectDirection(
        const GeometryType& Geom,
        const PointType& PointDestiny,
        PointType& PointProjected,
        const array_1d<double,3>& Normal,
        const array_1d<double,3>& Vector
        )
    {
        // Zero tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // We define the distance
        double distance = 0.0;

        const array_1d<double,3> vector_points = Geom[0].Coordinates() - PointDestiny.Coordinates();

        if( norm_2( Vector ) < zero_tolerance && norm_2( Normal ) > zero_tolerance ) {
            distance = inner_prod(vector_points, Normal)/norm_2(Normal);

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
            KRATOS_WARNING("Warning: Zero projection vector.") << " Projection using the condition vector instead." << std::endl;
        } else if (std::abs(inner_prod(Vector, Normal) ) > zero_tolerance) {
            distance = inner_prod(vector_points, Normal)/inner_prod(Vector, Normal);

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
        } else {
            PointProjected.Coordinates() = PointDestiny.Coordinates();
            KRATOS_WARNING("Warning: The line and the plane are coplanar.")  << " Something wrong happened " << std::endl;
        }

        return distance;
    }

    /**
     * @brief Project a point over a plane (avoiding some steps)
     * @param PointOrigin A point in the plane
     * @param PointDestiny The point to be projected
     * @param Normal The normal of the plane
     * @param Distance The distance to the projection
     * @return PointProjected The point pojected over the plane
     */

    KRATOS_DEPRECATED_MESSAGE("Method moved to geometrical_projection_utilities.h. Please update your declaration") static inline PointType FastProject(
        const PointType& PointOrigin,
        const PointType& PointDestiny,
        const array_1d<double,3>& Normal,
        double& Distance
        )
    {
        const array_1d<double,3> vector_points = PointDestiny.Coordinates() - PointOrigin.Coordinates();

        Distance = inner_prod(vector_points, Normal);

        PointType point_projected;
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        point_projected.Coordinates() = PointDestiny.Coordinates() - Normal * Distance;
#else
        noalias(point_projected.Coordinates()) = PointDestiny.Coordinates() - Normal * Distance;
#endif // ifdef KRATOS_USE_AMATRIX

        return point_projected;
    }

    /**
     * @brief Projects iteratively to get the coordinate
     * @param GeomOrigin The origin geometry
     * @param PointDestiny The destination point
     * @param ResultingPoint The distance between the point and the plane
     * @return Inside True is inside, false not
     */

    KRATOS_DEPRECATED_MESSAGE("Method moved to geometrical_projection_utilities.h. Please update your declaration") static inline bool ProjectIterativeLine2D(
        GeometryType& GeomOrigin,
        const GeometryType::CoordinatesArrayType& PointDestiny,
        GeometryType::CoordinatesArrayType& ResultingPoint,
        const array_1d<double, 3>& Normal,
        const double Tolerance = 1.0e-8,
        double DeltaXi = 0.5
        )
    {
//         ResultingPoint.clear();

        double old_delta_xi = 0.0;

        array_1d<double, 3> current_global_coords;

        array_1d<array_1d<double, 3>, 2> normals;
        normals[0] = GeomOrigin[0].FastGetSolutionStepValue(NORMAL);
        normals[1] = GeomOrigin[1].FastGetSolutionStepValue(NORMAL);

        BoundedMatrix<double,2,2> X;
        BoundedMatrix<double,2,1> DN;
        for(unsigned int i=0; i<2;++i) {
            X(0,i) = GeomOrigin[i].X();
            X(1,i) = GeomOrigin[i].Y();
        }

        Matrix J = ZeroMatrix( 1, 1 );

        //Newton iteration:

        const unsigned int max_iter = 20;

        for ( unsigned int k = 0; k < max_iter; ++k ) {
            array_1d<double, 2> N_origin;
            N_origin[0] = 0.5 * ( 1.0 - ResultingPoint[0]);
            N_origin[1] = 0.5 * ( 1.0 + ResultingPoint[0]);

            array_1d<double,3> normal_xi = N_origin[0] * normals[0] + N_origin[1] * normals[1];
            normal_xi /= norm_2(normal_xi);

            current_global_coords = ZeroVector(3);
            for( IndexType i_node = 0; i_node < 2; ++i_node )
                current_global_coords += N_origin[i_node] * GeomOrigin[i_node].Coordinates();

            const array_1d<double,3> VectorPoints = GeomOrigin.Center() - PointDestiny;
            const double distance = inner_prod(VectorPoints, Normal)/inner_prod(-normal_xi, Normal);
            const array_1d<double, 3> current_destiny_global_coords = PointDestiny - normal_xi * distance;

            // Derivatives of shape functions
            Matrix ShapeFunctionsGradients;
            ShapeFunctionsGradients = GeomOrigin.ShapeFunctionsLocalGradients(ShapeFunctionsGradients, ResultingPoint );

        #ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
            DN = prod(X,ShapeFunctionsGradients);

            J = prod(trans(DN),DN); // TODO: Add the non linearity concerning the normal
        #else
            noalias(DN) = prod(X,ShapeFunctionsGradients);

            noalias(J) = prod(trans(DN),DN); // TODO: Add the non linearity concerning the normal
        #endif // ifdef KRATOS_USE_AMATRIX

            const array_1d<double, 3>  temp = current_destiny_global_coords - current_global_coords;
            const Vector RHS = prod(trans(DN), subrange(temp,0,2));

            old_delta_xi = DeltaXi;
            DeltaXi = RHS[0]/J(0, 0);

            ResultingPoint[0] += DeltaXi;

            if (ResultingPoint[0] <= -1.0)
                ResultingPoint[0] = -1.0;
            else if (ResultingPoint[0] >= 1.0)
                ResultingPoint[0] = 1.0;

            if ( std::abs(DeltaXi - old_delta_xi) < Tolerance )
                return true;
        }

        return false;
    }

    /**
     * @brief This functions checks if the length of the line is to short, with the potential of provoque ill condition in the dual LM formulation
     * @param GeometryLine The line to be checked
     * @param Tolerance The threshold length
     * @return True if the line is too short, false otherwise
     */

    static inline bool LengthCheck(
        const GeometryPointType& GeometryLine,
        const double Tolerance = 1.0e-6
        ) {
        const double lx = GeometryLine[0].X() - GeometryLine[1].X();
        const double ly = GeometryLine[0].Y() - GeometryLine[1].Y();

        const double length = std::sqrt(lx * lx + ly * ly);

        return (length < Tolerance) ? true : false;
    }

    /**
     * @brief This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param GeometryTriangle The triangle to be checked
     * @return True if the triangle is in bad shape, false otherwise
     */

    static inline bool HeronCheck(const GeometryPointType& GeometryTriangle) {
        return HeronCheck(GeometryTriangle[0], GeometryTriangle[1], GeometryTriangle[2]);
    }

    /**
     * @brief This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param PointOrig1 The triangle first point
     * @param PointOrig2 The triangle second point
     * @param PointOrig3 The triangle third point
     * @return True if the triangle is in bad shape, false otherwise
     */

    static inline bool HeronCheck(
        const PointType& PointOrig1,
        const PointType& PointOrig2,
        const PointType& PointOrig3
        )
    {
        const double a = MathUtils<double>::Norm3(PointOrig1.Coordinates()-PointOrig2.Coordinates());
        const double b = MathUtils<double>::Norm3(PointOrig2.Coordinates()-PointOrig3.Coordinates());
        const double c = MathUtils<double>::Norm3(PointOrig3.Coordinates()-PointOrig1.Coordinates());

        const double s = 0.5 * (a + b + c);
        const double A2 = s * (s - a) * (s - b) * (s - c);

        const bool Check = A2 <= 0.0 ? true : false;  // We consider as bad shaped the ones with no area or negative A2 (semiperimeter smaller than any side)

//         // Debug
//         KRATOS_INFO("Check") << Check << " A2: " << A2 << std::endl;
//         if (Check == true) {
//             KRATOS_WARNING("Bad shape") << "Warning:: The triangle is in bad shape" << std::endl;
//             KRATOS_INFO("Mathematica triangle") << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << PointOrig1.X() << "," << PointOrig1.Y() << "," << PointOrig1.Z()  << "},{" << PointOrig2.X() << "," << PointOrig2.Y() << "," << PointOrig2.Z()  << "},{" << PointOrig3.X() << "," << PointOrig3.Y() << "," << PointOrig3.Z()  << "}}]}]" << std::endl;
//         }

        return Check;
    }

    /**
     * @brief This function rotates to align the projected points to a parallel plane to XY
     * @param PointToRotate The points from the origin geometry and the the point rotated
     * @param PointReferenceRotation The center point used as reference to rotate
     * @param SlaveTangentXi The first tangent vector of the slave condition
     * @param SlaveTangentEta The second tangent vector of the slave condition
     * @param Inversed If we rotate to the XY or we recover from XY
     */

    static inline void RotatePoint(
        PointType& PointToRotate,
        const PointType& PointReferenceRotation,
        const array_1d<double, 3>& SlaveTangentXi,
        const array_1d<double, 3>& SlaveTangentEta,
        const bool Inversed
        )
    {
        // We move to the (0,0,0)
        PointType aux_point_to_rotate;
        aux_point_to_rotate.Coordinates() = PointToRotate.Coordinates() - PointReferenceRotation.Coordinates();

        BoundedMatrix<double, 3, 3> rotation_matrix = ZeroMatrix(3, 3);

        if (Inversed == false) {
            for (IndexType i = 0; i < 3; ++i) {
                rotation_matrix(0, i) = SlaveTangentXi[i];
                rotation_matrix(1, i) = SlaveTangentEta[i];
            }
        } else {
            for (IndexType i = 0; i < 3; ++i) {
                rotation_matrix(i, 0) = SlaveTangentXi[i];
                rotation_matrix(i, 1) = SlaveTangentEta[i];
            }
        }

        PointToRotate.Coordinates() = prod(rotation_matrix, aux_point_to_rotate) + PointReferenceRotation.Coordinates();
    }

    /**
     * @brief This function calculates the normal in a specific GP with a given shape function
     * @param N The shape function considered
     * @param Geom The geometry of condition of interest
     * @return The normal in the GP
     */

    static inline array_1d<double,3> GaussPointUnitNormal(
        const Vector& N,
        const GeometryType& Geom
        ) {
        array_1d<double,3> normal = ZeroVector(3);
        for( IndexType i_node = 0; i_node < Geom.PointsNumber(); ++i_node )
            normal += N[i_node] * Geom[i_node].FastGetSolutionStepValue(NORMAL);

        const double this_norm = norm_2(normal);

        KRATOS_DEBUG_ERROR_IF(this_norm < std::numeric_limits<double>::epsilon()) << "Zero norm normal vector. Norm:" << this_norm << std::endl;

        normal /= this_norm;

        return normal;
    }

    /**
     * @brief This function gives you the indexes needed to order a vector
     * @param ThisVector The vector to order
     * @return idx The vector of indexes
     */

    template <typename TType>
    static std::vector<std::size_t> SortIndexes(const std::vector<TType> &ThisVector) {
        // Initialize original index locations
        std::vector<std::size_t> idx(ThisVector.size());
        iota(idx.begin(), idx.end(), 0);

        // Sort indexes based on comparing values in ThisVector
        std::sort(idx.begin(), idx.end(),
            [&ThisVector](std::size_t i1, std::size_t i2) {return ThisVector[i1] < ThisVector[i2];});

        return idx;
    }

    /**
     * @brief It computes the mean of the normal in the condition in all the nodes
     * @param rModelPart The model part to compute
     */

    static inline void ComputeNodesMeanNormalModelPart(ModelPart& rModelPart) {
        NodesArrayType& nodes_array = rModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());

        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
            (nodes_array.begin() + i)->FastGetSolutionStepValue(NORMAL) = zero_array;
#else
            noalias((nodes_array.begin() + i)->FastGetSolutionStepValue(NORMAL)) = zero_array;
#endif // ifdef KRATOS_USE_AMATRIX

        // Sum all the nodes normals
        ConditionsArrayType& conditions_array = rModelPart.Conditions();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
            auto it_cond = conditions_array.begin() + i;
            GeometryType& this_geometry = it_cond->GetGeometry();

            // Aux coordinates
            CoordinatesArrayType aux_coords;
            aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_geometry.Center());

            it_cond->SetValue(NORMAL, this_geometry.UnitNormal(aux_coords));

            const SizeType number_nodes = this_geometry.PointsNumber();

            for (IndexType i = 0; i < number_nodes; ++i) {
                auto& this_node = this_geometry[i];
                aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_node.Coordinates());
                const array_1d<double, 3> normal = this_geometry.UnitNormal(aux_coords);
                auto& aux_normal = this_node.FastGetSolutionStepValue(NORMAL);
                for (IndexType index = 0; index < 3; ++index) {
                    #pragma omp atomic
                    aux_normal[index] += normal[index];
                }
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = nodes_array.begin() + i;

            array_1d<double, 3>& normal = it_node->FastGetSolutionStepValue(NORMAL);
            const double norm_normal = norm_2(normal);

            if (norm_normal > std::numeric_limits<double>::epsilon()) normal /= norm_normal;
            else KRATOS_ERROR_IF(it_node->Is(INTERFACE)) << "ERROR:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
        }
    }

    /**
     * @brief It inverts the order of the nodes in the conditions of a model part in order to invert the normal
     * @param rContainer reference to the objective container
     */

    template<class TContainerType>
    static inline void InvertNormal(TContainerType& rContainer) {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rContainer.size()); ++i) {
            auto it_cont = rContainer.begin() + i;
            GeometryType& this_geometry = it_cont->GetGeometry();

            auto& data_geom = this_geometry.GetContainer();
            std::reverse(data_geom.begin(), data_geom.end());
        }
    }

    /**
     * @brief It calculates the matrix of coordinates of a geometry
     * @param ThisNodes The geometry to calculate
     * @param Current If we calculate the Current coordinates or the initial ones
     * @param Step The time step where it is computed
     * @return coordinates The matrix containing the coordinates of the geometry
     */

    template< SizeType TDim, SizeType TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TDim> GetCoordinates(
        const GeometryType& ThisNodes,
        const bool Current = true,
        const IndexType Step = 0
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, TDim> coordinates;
        array_1d<double, 3> coord;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            if (Current) {
                coord = ThisNodes[i_node].Coordinates();
            } else {
                coord = ThisNodes[i_node].GetInitialPosition();

                if (Step > 0)
                    coord += ThisNodes[i_node].FastGetSolutionStepValue(DISPLACEMENT, Step);
            }

            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof)
                coordinates(i_node, i_dof) = coord[i_dof];
        }

        return coordinates;
    }

    /**
     * @brief It calculates the matrix containing the tangent vector of the LM (for frictional contact)
     * @param ThisNodes The geometry to calculate
     * @return tangent_matrix The matrix containing the tangent vectors of the LM
     */

    template< SizeType TNumNodes, SizeType TDim>
    static inline BoundedMatrix<double, TNumNodes, TDim> ComputeTangentMatrix(const GeometryType& ThisNodes) {
        /* DEFINITIONS */
        // Zero tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();
        // Tangent matrix
        BoundedMatrix<double, TNumNodes, TDim> tangent_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& lm = ThisNodes[i_node].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
            if (norm_2(lm) > zero_tolerance) { // Non zero LM
                const array_1d<double, 3>& normal = ThisNodes[i_node].FastGetSolutionStepValue(NORMAL);
                const array_1d<double, 3> tangent_lm = lm - inner_prod(lm, normal) * normal;
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
     * @param ThisNodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param Step The step where it is computed
     * @return var_vector The vector containing the variables of the geometry
     */

    template< SizeType TNumNodes, class TVarType = Variable<double>>
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& ThisNodes,
        const TVarType& rVariable,
        const IndexType Step
        ) {
        /* DEFINITIONS */
        array_1d<double, TNumNodes> var_vector;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector[i_node] = ThisNodes[i_node].FastGetSolutionStepValue(rVariable, Step);

        return var_vector;
    }

    /**
     * @brief It calculates the vector of an historical variable of a geometry
     * @param ThisNodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param Step The step where it is computed
     * @return var_vector The vector containing the variables of the geometry
     */

    template< SizeType TNumNodes, class TVarType = Variable<double> >
    static inline BoundedMatrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& ThisNodes,
        const TVarType& rVariable,
        const unsigned int Step
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, 1> var_vector;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector(i_node, 0) = ThisNodes[i_node].FastGetSolutionStepValue(rVariable, Step);

        return var_vector;
    }

    /**
     * @brief It calculates the vector of a non-historical variable of a geometry
     * @param ThisNodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_vector The vector containing the variables of the geometry
     */

    template< SizeType TNumNodes, class TVarType = Variable<double> >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& ThisNodes,
        const TVarType& rVariable
        ) {
        /* DEFINITIONS */
        array_1d<double, TNumNodes> var_vector;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector[i_node] = ThisNodes[i_node].GetValue(rVariable);

        return var_vector;
    }

    /**
     * @brief It calculates the vector of a non-historical variable of a geometry
     * @param ThisNodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_vector The vector containing the variables of the geometry
     */

    template< SizeType TNumNodes, class TVarType = Variable<double> >
    static inline BoundedMatrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& ThisNodes,
        const TVarType& rVariable
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, 1> var_vector;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector(i_node, 0) = ThisNodes[i_node].GetValue(rVariable);

        return var_vector;
    }

    /**
     * @brief It calculates the matrix of a variable of a geometry
     * @param Nodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param Step The step where it is computed
     * @return var_matrix The matrix containing the variables of the geometry
     */

    template< SizeType TDim, SizeType TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TDim> GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVariable,
        const unsigned int Step
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, TDim> var_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& value = Nodes[i_node].FastGetSolutionStepValue(rVariable, Step);
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof)
                var_matrix(i_node, i_dof) = value[i_dof];
        }

        return var_matrix;
    }

    /**
     * @brief It calculates the matrix of a non-historical variable of a geometry
     * @param Nodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_matrix The matrix containing the variables of the geometry
     */

    template< SizeType TDim, SizeType TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TDim> GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVariable
        ) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, TDim> var_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& value = Nodes[i_node].GetValue(rVariable);
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof)
                var_matrix(i_node, i_dof) = value[i_dof];
        }

        return var_matrix;
    }

    /**
     * @brief It calculates the matrix containing the absolute value of another matrix
     * @param InputMatrix The original matrix
     * @return AbsMatrix The matrix containing the absolute value of another matrix
     */

    template< SizeType TDim, SizeType TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TDim> GetAbsMatrix(const BoundedMatrix<double, TNumNodes, TDim>& InputMatrix) {
        /* DEFINITIONS */
        BoundedMatrix<double, TNumNodes, TDim> AbsMatrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof)
                AbsMatrix(i_node, i_dof) = std::abs(InputMatrix(i_node, i_dof));
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
        Node<3>::Pointer pThisNode,
        unsigned int iSize
        );

    /**
     * @brief This method adds the value
     * @param rThisGeometry The geometrty to update
     * @param rThisVariable The variable to set
     */
    template< class TVarType, HistoricalValues THist>
    static inline void MatrixValue(
        GeometryType& rThisGeometry,
        TVarType& rThisVariable,
        Matrix& ThisValue
        );

    /**
     * @brief This method adds the value
     * @warning This operation is not threadsafe
     * @param rThisGeometry The geometrty to update
     * @param rThisVariable The variable to set
     */
    template< class TVarType, HistoricalValues THist>
    static inline void AddValue(
        GeometryType& rThisGeometry,
        TVarType& rThisVariable,
        const Matrix& ThisValue
        );

    /**
     * @brief This method adds the value
     * @param pThisNode The node to update
     * @param rThisVariable The variable to set
     */
    template< class TVarType, HistoricalValues THist>
    static inline void AddAreaWeightedNodalValue(
        Node<3>::Pointer pThisNode,
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
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetScalarVar(rThisVariable, 0.0, nodes_array);
}

template<>
inline void MortarUtilities::ResetValue<Variable<array_1d<double, 3>>, Historical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& rThisVariable
        ) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetVectorVar(rThisVariable, ZeroVector(3), nodes_array);
}

template<>
inline void MortarUtilities::ResetValue<Variable<double>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<double>& rThisVariable
        ) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(rThisVariable, 0.0, nodes_array);
}

template<>
inline void MortarUtilities::ResetValue<Variable<array_1d<double, 3>>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& rThisVariable
        ) {
    const array_1d<double, 3> zero_array = ZeroVector(3);
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(rThisVariable, zero_array, nodes_array);
}

template<>
inline void MortarUtilities::ResetAuxiliarValue<Variable<double>>(ModelPart& rThisModelPart) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(NODAL_MAUX, 0.0, nodes_array);
}

template<>
inline void MortarUtilities::ResetAuxiliarValue<Variable<array_1d<double, 3>>>(ModelPart& rThisModelPart) {
    const array_1d<double, 3> zero_array = ZeroVector(3);
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(NODAL_VAUX, zero_array, nodes_array);
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
        Matrix& ThisValue
        ) {
    if (ThisValue.size1() != rThisGeometry.size() || ThisValue.size2() != 1)
        ThisValue.resize(rThisGeometry.size(), 1, false);

    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node)
        ThisValue(i_node, 0) = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
}

template<>
inline void MortarUtilities::MatrixValue<Variable<array_1d<double, 3>>, Historical>(
        GeometryType& rThisGeometry,
        Variable<array_1d<double, 3>>& rThisVariable,
        Matrix& ThisValue
        ) {
    const std::size_t num_nodes = rThisGeometry.size();
    const std::size_t dimension = rThisGeometry.WorkingSpaceDimension();
    if (ThisValue.size1() != num_nodes || ThisValue.size2() != dimension)
        ThisValue.resize(num_nodes, dimension, false);

    for (IndexType i_node = 0; i_node < num_nodes; ++i_node) {
        const auto& rvalue = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        for (IndexType i_dim = 0; i_dim < dimension; ++i_dim)
            ThisValue(i_node, i_dim) = rvalue[i_dim];
    }
}
template<>
inline void MortarUtilities::MatrixValue<Variable<double>, NonHistorical>(
        GeometryType& rThisGeometry,
        Variable<double>& rThisVariable,
        Matrix& ThisValue
        ) {
    if (ThisValue.size1() != rThisGeometry.size() || ThisValue.size2() != 1)
        ThisValue.resize(rThisGeometry.size(), 1, false);

    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node)
        ThisValue(i_node, 0) = rThisGeometry[i_node].GetValue(rThisVariable);
}

template<>
inline void MortarUtilities::MatrixValue<Variable<array_1d<double, 3>>, NonHistorical>(
        GeometryType& rThisGeometry,
        Variable<array_1d<double, 3>>& rThisVariable,
        Matrix& ThisValue
        ) {
    const std::size_t num_nodes = rThisGeometry.size();
    const std::size_t dimension = rThisGeometry.WorkingSpaceDimension();
    if (ThisValue.size1() != num_nodes || ThisValue.size2() != dimension)
        ThisValue.resize(num_nodes, dimension, false);

    for (IndexType i_node = 0; i_node < num_nodes; ++i_node) {
        const auto& rvalue = rThisGeometry[i_node].GetValue(rThisVariable);
        for (IndexType i_dim = 0; i_dim < dimension; ++i_dim)
            ThisValue(i_node, i_dim) = rvalue[i_dim];
    }
}

template<>
inline void MortarUtilities::AddValue<Variable<double>, Historical>(
        GeometryType& rThisGeometry,
        Variable<double>& rThisVariable,
        const Matrix& ThisValue
        ) {
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        double& aux_value = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        #pragma omp atomic
        aux_value += ThisValue(i_node, 0);
    }
}

template<>
inline void MortarUtilities::AddValue<Variable<array_1d<double, 3>>, Historical>(
        GeometryType& rThisGeometry,
        Variable<array_1d<double, 3>>& rThisVariable,
        const Matrix& ThisValue
        ) {
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        auto& aux_vector = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        for (unsigned int i_dim = 0; i_dim < rThisGeometry.WorkingSpaceDimension(); ++i_dim) {
            double& aux_value = aux_vector[i_dim];
            #pragma omp atomic
            aux_value += ThisValue(i_node, i_dim);
        }
    }
}
template<>
inline void MortarUtilities::AddValue<Variable<double>, NonHistorical>(
        GeometryType& rThisGeometry,
        Variable<double>& rThisVariable,
        const Matrix& ThisValue
        ) {
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        double& aux_value = rThisGeometry[i_node].GetValue(rThisVariable);
        #pragma omp atomic
        aux_value += ThisValue(i_node, 0);
    }
}

template<>
inline void MortarUtilities::AddValue<Variable<array_1d<double, 3>>, NonHistorical>(
        GeometryType& rThisGeometry,
        Variable<array_1d<double, 3>>& rThisVariable,
        const Matrix& ThisValue
        ) {
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        auto& aux_vector = rThisGeometry[i_node].GetValue(rThisVariable);
        for (unsigned int i_dim = 0; i_dim < rThisGeometry.WorkingSpaceDimension(); ++i_dim) {
            double& aux_value = aux_vector[i_dim];
            #pragma omp atomic
            aux_value += ThisValue(i_node, i_dim);
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
