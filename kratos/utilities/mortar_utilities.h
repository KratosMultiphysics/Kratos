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

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "includes/model_part.h"
#include "geometries/point.h"
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
    
class MortarUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef Node<3>                                              NodeType;
    typedef Point<3>                                            PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;
    typedef Geometry<NodeType>                               GeometryType;
    typedef Geometry<PointType>                         GeometryPointType;
    typedef GeometryData::IntegrationMethod             IntegrationMethod;
    typedef ModelPart::NodesContainerType                  NodesArrayType;
    typedef ModelPart::ConditionsContainerType        ConditionsArrayType;
    
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
     * Project a point over a line/plane following an arbitrary direction
     * @param Geom: The geometry where to be projected
     * @param PointDestiny: The point to be projected
     * @param Normal: The normal of the geometry
     * @param Vector: The direction to project
     * @return PointProjected: The point pojected over the plane
     * @return Distance: The distnace between surfaces
     */

    static inline double FastProjectDirection(
        const GeometryType& Geom,
        const PointType& PointDestiny,
        PointType& PointProjected,
        const array_1d<double,3>& Normal,
        const array_1d<double,3>& Vector
        )
    {    
        // We define the tolerance
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        // We define the distance
        double distance = 0.0;
        
        const array_1d<double,3> vector_points = Geom[0].Coordinates() - PointDestiny.Coordinates();

        if( norm_2( Vector ) < tolerance && norm_2( Normal ) > tolerance )
        {
            distance = inner_prod(vector_points, Normal)/norm_2(Normal);

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
            std::cout << " :: Warning: Zero projection vector. Projection using the condition vector instead." << std::endl;
        }
        else if (std::abs(inner_prod(Vector, Normal) ) > tolerance)
        {
            distance = inner_prod(vector_points, Normal)/inner_prod(Vector, Normal); 

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
        }
        else
        {
            PointProjected.Coordinates() = PointDestiny.Coordinates();
            std::cout << " The line and the plane are coplanar, something wrong happened " << std::endl;
        }
        
        return distance;
    }
    
    /**
     * Project a point over a plane (avoiding some steps)
     * @param PointOrigin: A point in the plane
     * @param PointDestiny: The point to be projected
     * @param Normal: The normal of the plane
     * @return PointProjected: The point pojected over the plane
     */
    
    static inline PointType FastProject(
        const PointType& PointOrigin,
        const PointType& PointDestiny,
        const array_1d<double,3>& Normal
        )
    {
        array_1d<double,3> vector_points = PointDestiny.Coordinates() - PointOrigin.Coordinates();

        const double distance = inner_prod(vector_points, Normal); 
        
        PointType point_projected;
        point_projected.Coordinates() = PointDestiny.Coordinates() - Normal * distance;
        
        return point_projected;
    }
    
    /**
     * Projects iteratively to get the coordinate
     * @param GeomOrigin: The origin geometry
     * @param PointDestiny: The destination point
     * @return ResultingPoint: The distance between the point and the plane
     * @return Inside: True is inside, false not
     */
    
    static inline bool ProjectIterativeLine2D(
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
        normals[0] = GeomOrigin[0].GetValue(NORMAL);
        normals[1] = GeomOrigin[1].GetValue(NORMAL);
        
        bounded_matrix<double,2,2> X;
        bounded_matrix<double,2,1> DN;
        for(unsigned int i=0; i<2;i++)
        {
            X(0,i) = GeomOrigin[i].X();
            X(1,i) = GeomOrigin[i].Y();
        }

        Matrix J = ZeroMatrix( 1, 1 );
        
        //Newton iteration:

        const unsigned int max_iter = 20;

        for ( unsigned int k = 0; k < max_iter; k++ )
        {
            array_1d<double, 2> N_origin;
            N_origin[0] = 0.5 * ( 1.0 - ResultingPoint[0]);
            N_origin[1] = 0.5 * ( 1.0 + ResultingPoint[0]);
            
            array_1d<double,3> normal_xi = ZeroVector(3);
            for( unsigned int i_node = 0; i_node < 2; ++i_node )
            {
                normal_xi += N_origin[i_node] * normals[i_node]; 
            }
            
            normal_xi = normal_xi/norm_2(normal_xi); 
            
            current_global_coords = ZeroVector( 3 );
            for( unsigned int i_node = 0; i_node < 2; ++i_node )
            {
                current_global_coords += N_origin[i_node] * GeomOrigin[i_node].Coordinates(); 
            }
            
            const array_1d<double,3> VectorPoints = GeomOrigin.Center() - PointDestiny;
            const double distance = inner_prod(VectorPoints, Normal)/inner_prod(-normal_xi, Normal); 
            const array_1d<double, 3> current_destiny_global_coords = PointDestiny - normal_xi * distance;
            
            // Derivatives of shape functions
            Matrix ShapeFunctionsGradients;
            ShapeFunctionsGradients = GeomOrigin.ShapeFunctionsLocalGradients(ShapeFunctionsGradients, ResultingPoint );
            noalias(DN) = prod(X,ShapeFunctionsGradients);

            noalias(J) = prod(trans(DN),DN); // TODO: Add the non linearity concerning the normal
            Vector RHS = prod(trans(DN),subrange(current_destiny_global_coords - current_global_coords,0,2));
            
            old_delta_xi = DeltaXi;
            DeltaXi = RHS[0]/J(0, 0);
            
            ResultingPoint[0] += DeltaXi;
            
            if (ResultingPoint[0] <= -1.0)
            {
                ResultingPoint[0] = -1.0;
            }
            else if (ResultingPoint[0] >= 1.0)
            {
                ResultingPoint[0] = 1.0;
            }
            
            if ( std::abs(DeltaXi - old_delta_xi) < Tolerance )
            {
                return true;
            }
        }
        
        return false;
    }
    
    /**
     * This functions checks if the length of the line is to short, with the potential of provoque ill condition in the dual LM formulation
     * @param GeometryLine: The line to be checked
     * @param Tolerance: The threshold length
     * @return True if the line is too short, false otherwise
     */
    
    static inline bool LengthCheck(
        const GeometryPointType& GeometryLine,
        const double Tolerance = 1.0e-6
        )
    {
        const double lx = GeometryLine[0].X() - GeometryLine[1].X();
        const double ly = GeometryLine[0].Y() - GeometryLine[1].Y();

        const double length = std::sqrt(lx * lx + ly * ly);
        
        return (length < Tolerance) ? true : false;
    }
    
    /**
     * This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param GeometryTriangle: The triangle to be checked
     * @return True if the triangle is in bad shape, false otherwise
     */
    
    static inline bool HeronCheck(const GeometryPointType GeometryTriangle)
    {
        return HeronCheck(GeometryTriangle[0], GeometryTriangle[1], GeometryTriangle[2]);
    }
    
    /**
     * This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param PointOrig1: The triangle first point
     * @param PointOrig2: The triangle second point
     * @param PointOrig3: The triangle third point
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
//         std::cout << Check << " A2: " << A2 << std::endl;
//         if (Check == true)
//         {
//             std::cout << "Warning:: The triangle is in bad shape" << std::endl;
//             std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << PointOrig1.X() << "," << PointOrig1.Y() << "," << PointOrig1.Z()  << "},{" << PointOrig2.X() << "," << PointOrig2.Y() << "," << PointOrig2.Z()  << "},{" << PointOrig3.X() << "," << PointOrig3.Y() << "," << PointOrig3.Z()  << "}}]}]" << std::endl;
//         }
        
        return Check;
    }

    /**
     * This function rotates to align the projected points to a parallel plane to XY
     * @param PointToRotate: The points from the origin geometry
     * @param PointReferenceRotation: The center point used as reference to rotate
     * @param SlaveNormal: The normal vector of the slave condition
     * @param SlaveTangentXi: The first tangent vector of the slave condition
     * @param SlaveTangentEta: The second tangent vector of the slave condition
     * @param Inversed: If we rotate to the XY or we recover from XY
     * @return PointRotated: The point rotated 
     */
    
    static inline void RotatePoint( 
        PointType& PointToRotate,
        const PointType PointReferenceRotation,
        const array_1d<double, 3> SlaveTangentXi,
        const array_1d<double, 3> SlaveTangentEta,
        const bool Inversed
        )
    {                
        // We move to the (0,0,0)
        PointType aux_point_to_rotate;
        aux_point_to_rotate.Coordinates() = PointToRotate.Coordinates() - PointReferenceRotation.Coordinates();
        
        boost::numeric::ublas::bounded_matrix<double, 3, 3> rotation_matrix = ZeroMatrix(3, 3);
        
        if (Inversed == false)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                rotation_matrix(0, i) = SlaveTangentXi[i];
                rotation_matrix(1, i) = SlaveTangentEta[i];
            }
        }
        else
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                rotation_matrix(i, 0) = SlaveTangentXi[i];
                rotation_matrix(i, 1) = SlaveTangentEta[i];
            }
        }
        
        PointToRotate.Coordinates() = prod(rotation_matrix, aux_point_to_rotate) + PointReferenceRotation.Coordinates();
    }
    
    /**
     * This function calculates the normal in a specific GP with a given shape function
     * @param N: The shape function considered
     * @param Geom: The geometry of condition of interest
     */

    static inline array_1d<double,3> GaussPointNormal(
        const Vector& N,
        const GeometryType& Geom
        )
    {
        array_1d<double,3> normal = ZeroVector(3);
        for( unsigned int i_node = 0; i_node < Geom.PointsNumber(); ++i_node )
        {
            normal += N[i_node] * Geom[i_node].GetValue(NORMAL); 
        }
        
        if (norm_2(normal) > std::numeric_limits<double>::epsilon())
        {
            normal = normal/norm_2(normal); // It is suppossed to be already unitary (just in case)
        }
        
        return normal;
    }
    
    /**
     * This function gives you the indexes needed to order a vector 
     * @param vect: The vector to order
     * @return idx: The vector of indexes
     */
    
    template <typename TType>
    static inline std::vector<std::size_t> SortIndexes(const std::vector<TType> &vect) 
    {
        // Initialize original index locations
        std::vector<std::size_t> idx(vect.size());
        iota(idx.begin(), idx.end(), 0);

        // Sort indexes based on comparing values in vect
        std::sort(idx.begin(), idx.end(),
            [&vect](std::size_t i1, std::size_t i2) {return vect[i1] < vect[i2];});

        return idx;
    }
    
    /**
     * It calculates the matrix of coordinates of a geometry
     * @param nodes: The geometry to calculate
     * @param current: If we calculate the current coordinates or the initial ones
     * @return coordinates: The matrix containing the coordinates of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline bounded_matrix<double, TNumNodes, TDim> GetCoordinates(
        const GeometryType& nodes,
        const bool current = true,
        const unsigned int step = 0
        )
    {
        /* DEFINITIONS */            
        bounded_matrix<double, TNumNodes, TDim> coordinates;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            array_1d<double, 3> coord;
            
            if (current == true)
            {
                coord = nodes[i_node].Coordinates();
            }
            else
            {
                coord = nodes[i_node].GetInitialPosition();
                
                if (step > 0)
                {
                    coord += nodes[i_node].FastGetSolutionStepValue(DISPLACEMENT, step);
                }
            }

            for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
            {
                coordinates(i_node, i_dof) = coord[i_dof];
            }
        }
        
        return coordinates;
    }

    /**
     * It calculates the vector of an historical variable of a geometry
     * @param nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return var_vector: The vector containing the variables of the geometry
     */
    
    template< unsigned int TNumNodes >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& nodes,
        const Variable<double>& rVarName,
        const unsigned int step = 0
        )
    {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            var_vector[i_node] = nodes[i_node].FastGetSolutionStepValue(rVarName, step);
        }
        
        return var_vector;
    }
    
    /**
     * It calculates the vector of an historical variable of a geometry
     * @param nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return var_vector: The vector containing the variables of the geometry
     */
        
    template< unsigned int TNumNodes >
    static inline bounded_matrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& nodes,
        const Variable<double>& rVarName,
        const unsigned int step = 0
        )
    {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, 1> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            var_vector(i_node, 0) = nodes[i_node].FastGetSolutionStepValue(rVarName, step);
        }
        
        return var_vector;
    }

    /**
     * It calculates the vector of a non-historical variable of a geometry
     * @param nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @return var_vector: The vector containing the variables of the geometry
     */
        
    template< unsigned int TNumNodes >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& nodes,
        const Variable<double>& rVarName
        )
    {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            var_vector[i_node] = nodes[i_node].GetValue(rVarName);
        }
        
        return var_vector;
    }
    
    /**
     * It calculates the vector of a non-historical variable of a geometry
     * @param nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @return var_vector: The vector containing the variables of the geometry
     */
    
    template< unsigned int TNumNodes >
    static inline bounded_matrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& nodes,
        const Variable<double>& rVarName
        )
    {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, 1> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            var_vector(i_node, 0) = nodes[i_node].GetValue(rVarName);
        }
        
        return var_vector;
    }
    
    /**
     * It calculates the matrix of a variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return var_matrix: The matrix containing the variables of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVarName,
        const unsigned int step
        )
    {
        /* DEFINITIONS */        
        Matrix var_matrix(TNumNodes, TDim);
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            const array_1d<double, 3> value = Nodes[i_node].FastGetSolutionStepValue(rVarName, step);
            for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
            {
                var_matrix(i_node, i_dof) = value[i_dof];
            }
        }
        
        return var_matrix;
    }

    /**
     * It calculates the matrix of a non-historical variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @return var_matrix: The matrix containing the variables of the geometry
     */
        
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVarName
        )
    {
        /* DEFINITIONS */        
        Matrix var_matrix(TNumNodes, TDim);
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            const array_1d<double, 3>& value = Nodes[i_node].GetValue(rVarName);
            for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
            {
                var_matrix(i_node, i_dof) = value[i_dof];
            }
        }
        
        return var_matrix;
    }
    
    /**
     * It calculates the matrix containing the absolute value of another matrix
     * @param InputMatrix: The original matrix
     * @return AbsMatrix: The matrix containing the absolute value of another matrix
     */
        
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline bounded_matrix<double, TNumNodes, TDim> GetAbsMatrix(const bounded_matrix<double, TNumNodes, TDim> InputMatrix)
    {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, TDim> AbsMatrix;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
            {
                AbsMatrix(i_node, i_dof) = std::abs(InputMatrix(i_node, i_dof));
            }
        }
        
        return AbsMatrix;
    }
    
private:
};// class MortarUtilities

///@name Explicit Specializations
///@{

}
#endif /* KRATOS_MORTAR_UTILITIES defined */
