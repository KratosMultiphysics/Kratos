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
#include <unordered_map>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "includes/enums.h"
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
    
    // Component type
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;  
    
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
    typedef Point                                               PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;
    typedef Geometry<NodeType>                               GeometryType;
    typedef Geometry<PointType>                         GeometryPointType;
    typedef GeometryData::IntegrationMethod             IntegrationMethod;
    typedef ModelPart::NodesContainerType                  NodesArrayType;
    typedef ModelPart::ConditionsContainerType        ConditionsArrayType;
    typedef std::unordered_map<int, int>                           IntMap;
    
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
     * @param Geom The geometry where to be projected
     * @param PointDestiny The point to be projected
     * @param Normal The normal of the geometry
     * @param Vector The direction to project
     * @return PointProjected The point pojected over the plane
     * @return Distance The distnace between surfaces
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

        if( norm_2( Vector ) < tolerance && norm_2( Normal ) > tolerance ) {
            distance = inner_prod(vector_points, Normal)/norm_2(Normal);

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
            std::cout << " :: Warning: Zero projection vector. Projection using the condition vector instead." << std::endl;
        }
        else if (std::abs(inner_prod(Vector, Normal) ) > tolerance) {
            distance = inner_prod(vector_points, Normal)/inner_prod(Vector, Normal); 

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
        }
        else {
            PointProjected.Coordinates() = PointDestiny.Coordinates();
            std::cout << " The line and the plane are coplanar, something wrong happened " << std::endl;
        }
        
        return distance;
    }
    
    /**
     * Project a point over a plane (avoiding some steps)
     * @param PointOrigin A point in the plane
     * @param PointDestiny The point to be projected
     * @param Normal The normal of the plane
     * @return PointProjected The point pojected over the plane
     */
    
    static inline PointType FastProject(
        const PointType& PointOrigin,
        const PointType& PointDestiny,
        const array_1d<double,3>& Normal,
        double& Distance
        )
    {
        array_1d<double,3> vector_points = PointDestiny.Coordinates() - PointOrigin.Coordinates();

        Distance = inner_prod(vector_points, Normal); 
        
        PointType point_projected;
        point_projected.Coordinates() = PointDestiny.Coordinates() - Normal * Distance;
        
        return point_projected;
    }
    
    /**
     * Projects iteratively to get the coordinate
     * @param GeomOrigin The origin geometry
     * @param PointDestiny The destination point
     * @param ResultingPoint The distance between the point and the plane
     * @return Inside True is inside, false not
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
        normals[0] = GeomOrigin[0].FastGetSolutionStepValue(NORMAL);
        normals[1] = GeomOrigin[1].FastGetSolutionStepValue(NORMAL);
        
        bounded_matrix<double,2,2> X;
        bounded_matrix<double,2,1> DN;
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
            
            array_1d<double,3> normal_xi(3, 0.0);
            for( unsigned int i_node = 0; i_node < 2; ++i_node )
                normal_xi += N_origin[i_node] * normals[i_node]; 
            
            normal_xi = normal_xi/norm_2(normal_xi); 
            
            current_global_coords = ZeroVector(3);
            for( unsigned int i_node = 0; i_node < 2; ++i_node )
                current_global_coords += N_origin[i_node] * GeomOrigin[i_node].Coordinates(); 
            
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
                ResultingPoint[0] = -1.0;
            else if (ResultingPoint[0] >= 1.0)
                ResultingPoint[0] = 1.0;
            
            if ( std::abs(DeltaXi - old_delta_xi) < Tolerance )
                return true;
        }
        
        return false;
    }
    
    /**
     * This functions checks if the length of the line is to short, with the potential of provoque ill condition in the dual LM formulation
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
     * This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param GeometryTriangle The triangle to be checked
     * @return True if the triangle is in bad shape, false otherwise
     */
    
    static inline bool HeronCheck(const GeometryPointType& GeometryTriangle) {
        return HeronCheck(GeometryTriangle[0], GeometryTriangle[1], GeometryTriangle[2]);
    }
    
    /**
     * This functions checks if the semiperimeter is smaller than any of the sides of the triangle
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
//         std::cout << Check << " A2: " << A2 << std::endl;
//         if (Check == true) {
//             std::cout << "Warning:: The triangle is in bad shape" << std::endl;
//             std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << PointOrig1.X() << "," << PointOrig1.Y() << "," << PointOrig1.Z()  << "},{" << PointOrig2.X() << "," << PointOrig2.Y() << "," << PointOrig2.Z()  << "},{" << PointOrig3.X() << "," << PointOrig3.Y() << "," << PointOrig3.Z()  << "}}]}]" << std::endl;
//         }
        
        return Check;
    }

    /**
     * This function rotates to align the projected points to a parallel plane to XY
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
        
        bounded_matrix<double, 3, 3> rotation_matrix = ZeroMatrix(3, 3);
        
        if (Inversed == false) {
            for (unsigned int i = 0; i < 3; ++i) {
                rotation_matrix(0, i) = SlaveTangentXi[i];
                rotation_matrix(1, i) = SlaveTangentEta[i];
            }
        }
        else {
            for (unsigned int i = 0; i < 3; ++i) {
                rotation_matrix(i, 0) = SlaveTangentXi[i];
                rotation_matrix(i, 1) = SlaveTangentEta[i];
            }
        }
        
        PointToRotate.Coordinates() = prod(rotation_matrix, aux_point_to_rotate) + PointReferenceRotation.Coordinates();
    }
    
    /**
     * This function calculates the normal in a specific GP with a given shape function
     * @param N The shape function considered
     * @param Geom The geometry of condition of interest
     */

    static inline array_1d<double,3> GaussPointUnitNormal(
        const Vector& N,
        const GeometryType& Geom
        ) {
        array_1d<double,3> normal(3, 0.0);
        for( unsigned int i_node = 0; i_node < Geom.PointsNumber(); ++i_node )
            normal += N[i_node] * Geom[i_node].FastGetSolutionStepValue(NORMAL); 
        
        const double this_norm = norm_2(normal);
        
    #ifdef KRATOS_DEBUG
        const bool not_zero_vector = (this_norm > std::numeric_limits<double>::epsilon());
        if (not_zero_vector == false) KRATOS_ERROR << "Zero norm normal vector. Norm:" << this_norm << std::endl;
    #endif
        
        normal /= this_norm;
        
        return normal;
    }
    
    /**
     * This function gives you the indexes needed to order a vector 
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
     * It computes the mean of the normal in the condition in all the nodes
     * @param rModelPart The model part to compute
     */
    
    static inline void ComputeNodesMeanNormalModelPart(ModelPart& rModelPart) {
        // Tolerance
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        NodesArrayType& nodes_array = rModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) 
            noalias((nodes_array.begin() + i)->FastGetSolutionStepValue(NORMAL)) = ZeroVector(3);
        
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
            
            const unsigned int number_nodes = this_geometry.PointsNumber();
            
            for (unsigned int i = 0; i < number_nodes; ++i) {
                auto& this_node = this_geometry[i];
                aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_node.Coordinates());
                const array_1d<double, 3>& normal = this_geometry.UnitNormal(aux_coords);
                auto& aux_normal = this_node.FastGetSolutionStepValue(NORMAL);
                for (unsigned int index = 0; index < 3; ++index) {
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
            if (norm_normal > tolerance) normal /= norm_normal;
            else KRATOS_ERROR << "WARNING:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
        }
    }
    
    /**
     * It calculates the matrix of coordinates of a geometry
     * @param ThisNodes The geometry to calculate
     * @param Current If we calculate the Current coordinates or the initial ones
     * @return coordinates The matrix containing the coordinates of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline bounded_matrix<double, TNumNodes, TDim> GetCoordinates(
        const GeometryType& ThisNodes,
        const bool Current = true,
        const unsigned int Step = 0
        ) {
        /* DEFINITIONS */            
        bounded_matrix<double, TNumNodes, TDim> coordinates;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
        {
            array_1d<double, 3> coord;
            
            if (Current == true)
                coord = ThisNodes[i_node].Coordinates();
            else {
                coord = ThisNodes[i_node].GetInitialPosition();
                
                if (Step > 0)
                    coord += ThisNodes[i_node].FastGetSolutionStepValue(DISPLACEMENT, Step);
            }

            for (unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
                coordinates(i_node, i_dof) = coord[i_dof];
        }
        
        return coordinates;
    }

    /**
     * It calculates the vector of an historical variable of a geometry
     * @param ThisNodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param Step The step where calculate
     * @return var_vector The vector containing the variables of the geometry
     */
    
    template< unsigned int TNumNodes, class TVarType = Variable<double>>
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& ThisNodes,
        const TVarType& rVariable,
        const unsigned int Step
        ) {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector[i_node] = ThisNodes[i_node].FastGetSolutionStepValue(rVariable, Step);
        
        return var_vector;
    }
    
    /**
     * It calculates the vector of an historical variable of a geometry
     * @param ThisNodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param Step The step where calculate
     * @return var_vector The vector containing the variables of the geometry
     */
        
    template< unsigned int TNumNodes, class TVarType = Variable<double> >
    static inline bounded_matrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& ThisNodes,
        const TVarType& rVariable,
        const unsigned int Step
        ) {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, 1> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector(i_node, 0) = ThisNodes[i_node].FastGetSolutionStepValue(rVariable, Step);

        return var_vector;
    }

    /**
     * It calculates the vector of a non-historical variable of a geometry
     * @param ThisNodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_vector The vector containing the variables of the geometry
     */
        
    template< unsigned int TNumNodes, class TVarType = Variable<double> >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& ThisNodes,
        const TVarType& rVariable
        ) {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector[i_node] = ThisNodes[i_node].GetValue(rVariable);
        
        return var_vector;
    }
    
    /**
     * It calculates the vector of a non-historical variable of a geometry
     * @param ThisNodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_vector The vector containing the variables of the geometry
     */
    
    template< unsigned int TNumNodes, class TVarType = Variable<double> >
    static inline bounded_matrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& ThisNodes,
        const TVarType& rVariable
        ) {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, 1> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            var_vector(i_node, 0) = ThisNodes[i_node].GetValue(rVariable);
        
        return var_vector;
    }
    
    /**
     * It calculates the matrix of a variable of a geometry
     * @param Nodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @param step The step where calculate
     * @return var_matrix The matrix containing the variables of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVariable,
        const unsigned int step
        ) {
        /* DEFINITIONS */        
        Matrix var_matrix(TNumNodes, TDim);
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3> value = Nodes[i_node].FastGetSolutionStepValue(rVariable, step);
            for (unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
                var_matrix(i_node, i_dof) = value[i_dof];
        }
        
        return var_matrix;
    }

    /**
     * It calculates the matrix of a non-historical variable of a geometry
     * @param Nodes The geometry to calculate
     * @param rVariable The name of the variable to calculate
     * @return var_matrix The matrix containing the variables of the geometry
     */
        
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVariable
        ) {
        /* DEFINITIONS */        
        Matrix var_matrix(TNumNodes, TDim);
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& value = Nodes[i_node].GetValue(rVariable);
            for (unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
                var_matrix(i_node, i_dof) = value[i_dof];
        }
        
        return var_matrix;
    }
    
    /**
     * It calculates the matrix containing the absolute value of another matrix
     * @param InputMatrix The original matrix
     * @return AbsMatrix The matrix containing the absolute value of another matrix
     */
        
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline bounded_matrix<double, TNumNodes, TDim> GetAbsMatrix(const bounded_matrix<double, TNumNodes, TDim>& InputMatrix) {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, TDim> AbsMatrix;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
            for (unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
                AbsMatrix(i_node, i_dof) = std::abs(InputMatrix(i_node, i_dof));
        }
        
        return AbsMatrix;
    }
    
    /**
     * This method gives the size to be computed
     */
    template< unsigned int TDim, class TVarType>
    static inline unsigned int SizeToCompute()
    {
       if (typeid(TVarType) == typeid(Variable<array_1d<double, 3>>))
           return TDim;
       
       return 1;
    }
    
    /**
     * This method resets the value
     * @param rThisModelPart The model part to update
     * @param ThisVariable The variable to set
     * @param InvertedPair If the master/slave follows the standard way 
     */
    template< class TVarType, HistoricalValues THist>
    static inline void ResetValue(
        ModelPart& rThisModelPart,
        TVarType& ThisVariable, 
        const bool InvertedPair = false
        );
    
    /**
     * This method resets the auxiliar value
     * @param rThisModelPart The model part to update
     */
    template< class TVarType>
    static inline void ResetAuxiliarValue(ModelPart& rThisModelPart);

    /**
     * This method returns the auxiliar variable
     */
    template< class TVarType>
    static inline TVarType GetAuxiliarVariable();

    /**
     * This method returns the auxiliar variable
     */
    template< class TVarType>
    static inline double GetAuxiliarValue(
        Node<3>::Pointer pThisNode,
        unsigned int iSize
        );
    
    /**
     * This method adds the value
     * @param ThisGeometry The geometrty to update
     * @param ThisVariable The variable to set
     */
    template< class TVarType, HistoricalValues THist>
    static inline void MatrixValue(
        GeometryType& ThisGeometry,
        TVarType& ThisVariable,
        Matrix& ThisValue
        );
    
    /**
     * This method adds the value
     * WARNING This operation is not threadsafe
     * @param ThisGeometry The geometrty to update
     * @param ThisVariable The variable to set
     */
    template< class TVarType, HistoricalValues THist>
    static inline void AddValue(
        GeometryType& ThisGeometry,
        TVarType& ThisVariable,
        const Matrix& ThisValue
        );
    
    /**
     * This method adds the value
     * @param pThisNode The node to update
     * @param ThisVariable The variable to set
     */
    template< class TVarType, HistoricalValues THist>
    static inline void AddAreaWeightedNodalValue(
        Node<3>::Pointer pThisNode,
        TVarType& ThisVariable,
        const double RefArea = 1.0,
        const double Tolerance = 1.0e-4
        );

    /**
     * This method updates the database in the amster side
     * @param rThisModelPart The model part
     * @param ThisVariable The variable to set
     * @param Dx The vector with the increment of the value
     * @param Index The index used in the  case of a vector variable
     * @param ConectivityDatabase The database that will be used to assemble the system
     */
    template< class TVarType, HistoricalValues THist>
    static inline void UpdateDatabase(
        ModelPart& rThisModelPart,
        TVarType& ThisVariable,
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
        Variable<double>& ThisVariable, 
        const bool InvertedPair
        ) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        if (it_node->Is(SLAVE) == !InvertedPair) 
            it_node->FastGetSolutionStepValue(ThisVariable) = 0.0;
    }
}

template<> 
inline void MortarUtilities::ResetValue<ComponentType, Historical>(
        ModelPart& rThisModelPart,
        ComponentType& ThisVariable, 
        const bool InvertedPair
        ) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        if (it_node->Is(SLAVE) == !InvertedPair) 
            it_node->FastGetSolutionStepValue(ThisVariable) = 0.0;
    }
}

template<> 
inline void MortarUtilities::ResetValue<Variable<array_1d<double, 3>>, Historical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& ThisVariable, 
        const bool InvertedPair
        ) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        if (it_node->Is(SLAVE) == !InvertedPair) {
            array_1d<double, 3>& aux_value = it_node->FastGetSolutionStepValue(ThisVariable);
            noalias(aux_value) = ZeroVector(3);
        }
    }
}

template<> 
inline void MortarUtilities::ResetValue<Variable<double>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<double>& ThisVariable, 
        const bool InvertedPair
        ) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        if (it_node->Is(SLAVE) == !InvertedPair) 
            it_node->SetValue(ThisVariable, 0.0);
    }
}

template<> 
inline void MortarUtilities::ResetValue<ComponentType, NonHistorical>(
        ModelPart& rThisModelPart,
        ComponentType& ThisVariable, 
        const bool InvertedPair
        ) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        if (it_node->Is(SLAVE) == !InvertedPair) 
            it_node->SetValue(ThisVariable, 0.0);
    }
}

template<> 
inline void MortarUtilities::ResetValue<Variable<array_1d<double, 3>>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& ThisVariable, 
        const bool InvertedPair
        ) {
    // Zero vector
    const array_1d<double, 3> zero_vector(3, 0.0);
    
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        if (it_node->Is(SLAVE) == !InvertedPair) 
            it_node->SetValue(ThisVariable, zero_vector);
    }
}

template<>
inline void MortarUtilities::ResetAuxiliarValue<Variable<double>>(ModelPart& rThisModelPart) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->SetValue(NODAL_MAUX, 0.0);
    }
}

template<>
inline void MortarUtilities::ResetAuxiliarValue<ComponentType>(ModelPart& rThisModelPart) {
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->SetValue(NODAL_VAUX_X, 0.0);
    }
}

template<>
inline void MortarUtilities::ResetAuxiliarValue<Variable<array_1d<double, 3>>>(ModelPart& rThisModelPart) {
    // Zero vector
    const array_1d<double, 3> zero_vector(3, 0.0);
    
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->SetValue(NODAL_VAUX, zero_vector);
    }
}

template< >
inline Variable<double> MortarUtilities::GetAuxiliarVariable<Variable<double>>() {
    return NODAL_MAUX;
}

template< >
inline ComponentType MortarUtilities::GetAuxiliarVariable<ComponentType>() {
    return NODAL_VAUX_X;
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
inline double MortarUtilities::GetAuxiliarValue<ComponentType>(
    Node<3>::Pointer pThisNode,
    unsigned int iSize
    ) {
    return pThisNode->GetValue(NODAL_VAUX_X);
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
        GeometryType& ThisGeometry,
        Variable<double>& ThisVariable,
        Matrix& ThisValue
        ) {
    if (ThisValue.size1() != ThisGeometry.size() || ThisValue.size2() != 1)
        ThisValue.resize(ThisGeometry.size(), 1, false);
    
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node)
        ThisValue(i_node, 0) = ThisGeometry[i_node].FastGetSolutionStepValue(ThisVariable);
}

template<> 
inline void MortarUtilities::MatrixValue<ComponentType, Historical>(
        GeometryType& ThisGeometry,
        ComponentType& ThisVariable,
        Matrix& ThisValue
        ) {
    if (ThisValue.size1() != ThisGeometry.size() || ThisValue.size2() != 1)
        ThisValue.resize(ThisGeometry.size(), 1, false);
    
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node)
        ThisValue(i_node, 0) = ThisGeometry[i_node].FastGetSolutionStepValue(ThisVariable);
}

template<> 
inline void MortarUtilities::MatrixValue<Variable<array_1d<double, 3>>, Historical>(
        GeometryType& ThisGeometry,
        Variable<array_1d<double, 3>>& ThisVariable,
        Matrix& ThisValue
        ) {    
    const std::size_t num_nodes = ThisGeometry.size();
    const std::size_t dimension = ThisGeometry.WorkingSpaceDimension();
    if (ThisValue.size1() != num_nodes || ThisValue.size2() != dimension)
        ThisValue.resize(num_nodes, dimension, false);
    
    for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
        row(ThisValue, i_node) = subrange(ThisGeometry[i_node].FastGetSolutionStepValue(ThisVariable), 0, dimension);
}
template<> 
inline void MortarUtilities::MatrixValue<Variable<double>, NonHistorical>(
        GeometryType& ThisGeometry,
        Variable<double>& ThisVariable,
        Matrix& ThisValue
        ) {
    if (ThisValue.size1() != ThisGeometry.size() || ThisValue.size2() != 1)
        ThisValue.resize(ThisGeometry.size(), 1, false);
    
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node)
        ThisValue(i_node, 0) = ThisGeometry[i_node].GetValue(ThisVariable);
}

template<> 
inline void MortarUtilities::MatrixValue<ComponentType, NonHistorical>(
        GeometryType& ThisGeometry,
        ComponentType& ThisVariable,
        Matrix& ThisValue
        ) {
    if (ThisValue.size1() != ThisGeometry.size() || ThisValue.size2() != 1)
        ThisValue.resize(ThisGeometry.size(), 1, false);
    
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node)
        ThisValue(i_node, 0) = ThisGeometry[i_node].GetValue(ThisVariable);
}

template<> 
inline void MortarUtilities::MatrixValue<Variable<array_1d<double, 3>>, NonHistorical>(
        GeometryType& ThisGeometry,
        Variable<array_1d<double, 3>>& ThisVariable,
        Matrix& ThisValue
        ) {
    const std::size_t num_nodes = ThisGeometry.size();
    const std::size_t dimension = ThisGeometry.WorkingSpaceDimension();
    if (ThisValue.size1() != num_nodes || ThisValue.size2() != dimension)
        ThisValue.resize(num_nodes, dimension, false);
    
    for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
        row(ThisValue, i_node) = subrange(ThisGeometry[i_node].GetValue(ThisVariable), 0, dimension);
}

template<> 
inline void MortarUtilities::AddValue<Variable<double>, Historical>(
        GeometryType& ThisGeometry,
        Variable<double>& ThisVariable,
        const Matrix& ThisValue
        ) {
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node)
        ThisGeometry[i_node].FastGetSolutionStepValue(ThisVariable) += ThisValue(i_node, 0);
}

template<> 
inline void MortarUtilities::AddValue<ComponentType, Historical>(
        GeometryType& ThisGeometry,
        ComponentType& ThisVariable,
        const Matrix& ThisValue
        ) {
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node)
        ThisGeometry[i_node].FastGetSolutionStepValue(ThisVariable) += ThisValue(i_node, 0);
}

template<> 
inline void MortarUtilities::AddValue<Variable<array_1d<double, 3>>, Historical>(
        GeometryType& ThisGeometry,
        Variable<array_1d<double, 3>>& ThisVariable,
        const Matrix& ThisValue
        ) {
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node) {
        auto& aux_vector = ThisGeometry[i_node].FastGetSolutionStepValue(ThisVariable);
        for (unsigned int i_dim = 0; i_dim < ThisGeometry.WorkingSpaceDimension(); ++i_dim)
            aux_vector[i_dim] += ThisValue(i_node, i_dim);
    }
}
template<> 
inline void MortarUtilities::AddValue<Variable<double>, NonHistorical>(
        GeometryType& ThisGeometry,
        Variable<double>& ThisVariable,
        const Matrix& ThisValue
        ) {
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node)
        ThisGeometry[i_node].GetValue(ThisVariable) += ThisValue(i_node, 0);
}

template<> 
inline void MortarUtilities::AddValue<ComponentType, NonHistorical>(
        GeometryType& ThisGeometry,
        ComponentType& ThisVariable,
        const Matrix& ThisValue
        ) {
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node)
        ThisGeometry[i_node].GetValue(ThisVariable) += ThisValue(i_node, 0);
}

template<> 
inline void MortarUtilities::AddValue<Variable<array_1d<double, 3>>, NonHistorical>(
        GeometryType& ThisGeometry,
        Variable<array_1d<double, 3>>& ThisVariable,
        const Matrix& ThisValue
        ) {
    for (unsigned int i_node = 0; i_node < ThisGeometry.size(); ++i_node) {
        auto& aux_vector = ThisGeometry[i_node].GetValue(ThisVariable);
        for (unsigned int i_dim = 0; i_dim < ThisGeometry.WorkingSpaceDimension(); ++i_dim)
            aux_vector[i_dim] += ThisValue(i_node, i_dim);
    }
}

template<> 
inline void MortarUtilities::AddAreaWeightedNodalValue<Variable<double>, Historical>(
        Node<3>::Pointer pThisNode,
        Variable<double>& ThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG 
    if (null_area) std::cout << "WARNING:: NODE OF NULL AREA. ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    pThisNode->FastGetSolutionStepValue(ThisVariable) += area_coeff * pThisNode->GetValue(NODAL_MAUX);
}

template<> 
inline void MortarUtilities::AddAreaWeightedNodalValue<ComponentType, Historical>(
        Node<3>::Pointer pThisNode,
        ComponentType& ThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG 
    if (null_area) std::cout << "WARNING:: NODE OF NULL AREA. ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    pThisNode->FastGetSolutionStepValue(ThisVariable) += area_coeff * pThisNode->GetValue(NODAL_VAUX_X);
}

template<> 
inline void MortarUtilities::AddAreaWeightedNodalValue<Variable<array_1d<double, 3>>, Historical>(
        Node<3>::Pointer pThisNode,
        Variable<array_1d<double, 3>>& ThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG 
    if (null_area) std::cout << "WARNING:: NODE OF NULL AREA. ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    auto& aux_vector = pThisNode->FastGetSolutionStepValue(ThisVariable);
    aux_vector += area_coeff * pThisNode->GetValue(NODAL_VAUX);
}

template<> 
inline void MortarUtilities::AddAreaWeightedNodalValue<Variable<double>, NonHistorical>(
        Node<3>::Pointer pThisNode,
        Variable<double>& ThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG 
    if (null_area) std::cout << "WARNING:: NODE OF NULL AREA. ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    pThisNode->GetValue(ThisVariable) += area_coeff * pThisNode->GetValue(NODAL_MAUX);
}

template<> 
inline void MortarUtilities::AddAreaWeightedNodalValue<ComponentType, NonHistorical>(
        Node<3>::Pointer pThisNode,
        ComponentType& ThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG 
    if (null_area) std::cout << "WARNING:: NODE OF NULL AREA. ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    pThisNode->GetValue(ThisVariable) += area_coeff * pThisNode->GetValue(NODAL_VAUX_X);
}

template<> 
inline void MortarUtilities::AddAreaWeightedNodalValue<Variable<array_1d<double, 3>>, NonHistorical>(
        Node<3>::Pointer pThisNode,
        Variable<array_1d<double, 3>>& ThisVariable,
        const double RefArea,
        const double Tolerance
        ) {
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG 
    if (null_area) std::cout << "WARNING:: NODE OF NULL AREA. ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    auto& aux_vector = pThisNode->GetValue(ThisVariable);
    aux_vector += area_coeff * pThisNode->GetValue(NODAL_VAUX);
}

template<> 
inline void MortarUtilities::UpdateDatabase<Variable<double>, Historical>(
        ModelPart& rThisModelPart,
        Variable<double>& ThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        p_node->FastGetSolutionStepValue(ThisVariable) += Dx[i];
    }
}

template<> 
inline void MortarUtilities::UpdateDatabase<ComponentType, Historical>(
        ModelPart& rThisModelPart,
        ComponentType& ThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        p_node->FastGetSolutionStepValue(ThisVariable) += Dx[i];
    }
}

template<> 
inline void MortarUtilities::UpdateDatabase<Variable<array_1d<double, 3>>, Historical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& ThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        auto& value = p_node->FastGetSolutionStepValue(ThisVariable); 
        value[Index] += Dx[i];
    }
}
template<> 
inline void MortarUtilities::UpdateDatabase<Variable<double>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<double>& ThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        p_node->GetValue(ThisVariable) += Dx[i];
    }
}

template<> 
inline void MortarUtilities::UpdateDatabase<ComponentType, NonHistorical>(
        ModelPart& rThisModelPart,
        ComponentType& ThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        p_node->GetValue(ThisVariable) += Dx[i];
    }
}

template<> 
inline void MortarUtilities::UpdateDatabase<Variable<array_1d<double, 3>>, NonHistorical>(
        ModelPart& rThisModelPart,
        Variable<array_1d<double, 3>>& ThisVariable,
        Vector& Dx,
        unsigned int Index,
        IntMap& ConectivityDatabase
        ) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        auto& value = p_node->GetValue(ThisVariable); 
        value[Index] += Dx[i];
    }
}

}
#endif /* KRATOS_MORTAR_UTILITIES defined */
