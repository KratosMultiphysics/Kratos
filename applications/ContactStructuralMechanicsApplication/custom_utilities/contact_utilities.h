// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_CONTACT_UTILITIES)
#define KRATOS_CONTACT_UTILITIES

#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application_variables.h"
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
    
class ContactUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Point<3>                                        PointType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef Geometry<PointType>                     GeometryPointType;
    typedef GeometryData::IntegrationMethod         IntegrationMethod;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    
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
        const double Tolerance = std::numeric_limits<double>::epsilon();
        
        // We define the distance
        double Distance = 0.0;
        
        const array_1d<double,3> VectorPoints = Geom[0].Coordinates() - PointDestiny.Coordinates();

        if( norm_2( Vector ) < Tolerance && norm_2( Normal ) > Tolerance )
        {
            Distance = inner_prod(VectorPoints, Normal)/norm_2(Normal);

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * Distance;
            std::cout << " :: Warning: Zero projection vector. Projection using the condition vector instead." << std::endl;
        }
        else if (std::abs(inner_prod(Vector, Normal) ) > Tolerance)
        {
            Distance = inner_prod(VectorPoints, Normal)/inner_prod(Vector, Normal); 

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * Distance;
        }
        else
        {
            PointProjected.Coordinates() = PointDestiny.Coordinates();
            std::cout << " The line and the plane are coplanar, something wrong happened " << std::endl;
        }
        
        return Distance;
    }
    
    /**
     * Project a point over a line/plane following an arbitrary direction
     * @param Geom: The geometry where to be projected
     * @param PointDestiny: The point to be projected
     * @param Vector: The direction to project
     * @return PointProjected: The point pojected over the plane
     * @return Distance: The distance between the point and the plane
     */

    static inline void ProjectDirection(
        const GeometryType& Geom,
        const PointType& PointDestiny,
        PointType& PointProjected,
        double& Distance,
        const array_1d<double,3>& Vector
        )
    {        
        array_1d<double,3> Normal;
        
        GeometryNormal(Normal, Geom);
        
        Distance = FastProjectDirection(Geom, PointDestiny, PointProjected, Normal,Vector);
    }
    
    static inline void ProjectCoordDirection(
        const GeometryType& Geom,
        const GeometryType::CoordinatesArrayType& CoordDestiny,
        GeometryType::CoordinatesArrayType& CoordProjected,
        double& Distance,
        const array_1d<double,3>& Vector
        )
    {        
        array_1d<double,3> Normal;
        
        GeometryNormal(Normal, Geom);
        
        PointType PointDestiny;
        PointDestiny.Coordinates() = CoordDestiny;
        
        PointType PointProjected;
        
        Distance = FastProjectDirection(Geom, PointDestiny, PointProjected, Normal,Vector);
        
        CoordProjected = PointProjected.Coordinates();
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
        const array_1d<double, 3> Normal,
        const double Tolerance = 1.0e-8
        )
    {
        ResultingPoint.clear();
        
        double OldDeltaXi = 0.0;
        double DeltaXi    = 0.5;

        array_1d<double, 3> CurrentGlobalCoords;

        array_1d<array_1d<double, 3>, 2> normals;
        normals[0] = GeomOrigin[0].GetValue(NORMAL);
        normals[1] = GeomOrigin[1].GetValue(NORMAL);
        
        boost::numeric::ublas::bounded_matrix<double,2,2> X;
        boost::numeric::ublas::bounded_matrix<double,2,1> DN;
        for(unsigned int i=0; i<2;i++)
        {
            X(0,i) = GeomOrigin[i].X();
            X(1,i) = GeomOrigin[i].Y();
        }

        Matrix J = ZeroMatrix( 1, 1 );
        
        //Newton iteration:

        unsigned int maxiter = 10;

        for ( unsigned int k = 0; k < maxiter; k++ )
        {
            array_1d<double, 2> NOrigin;
            NOrigin[0] = 0.5 * ( 1.0 - ResultingPoint[0]);
            NOrigin[1] = 0.5 * ( 1.0 + ResultingPoint[0]);
            
            array_1d<double,3> normal_xi = ZeroVector(3);
            for( unsigned int iNode = 0; iNode < 2; ++iNode )
            {
                normal_xi += NOrigin[iNode] * normals[iNode]; 
            }
            
            normal_xi = normal_xi/norm_2(normal_xi); 
            
            CurrentGlobalCoords = ZeroVector( 3 );
            for( unsigned int iNode = 0; iNode < 2; ++iNode )
            {
                CurrentGlobalCoords += NOrigin[iNode] * GeomOrigin[iNode].Coordinates(); 
            }
            
            const array_1d<double,3> VectorPoints = GeomOrigin.Center() - PointDestiny;
            const double dist = inner_prod(VectorPoints, Normal)/inner_prod(-normal_xi, Normal); 
            const array_1d<double, 3> CurrentDestinyGlobalCoords = PointDestiny - normal_xi * dist;
            
//             // Debug
//             KRATOS_WATCH(CurrentGlobalCoords)
//             KRATOS_WATCH(CurrentDestinyGlobalCoords)
            
            // Derivatives of shape functions
            Matrix shape_functions_gradients;
            shape_functions_gradients = GeomOrigin.ShapeFunctionsLocalGradients(shape_functions_gradients, ResultingPoint );
            noalias(DN) = prod(X,shape_functions_gradients);

            noalias(J) = prod(trans(DN),DN); // TODO: Add the non linearity concerning the normal
            Vector RHS = prod(trans(DN),subrange(CurrentDestinyGlobalCoords - CurrentGlobalCoords,0,2));
            
            OldDeltaXi = DeltaXi;
            DeltaXi = RHS[0]/J(0, 0);
            
            ResultingPoint[0] += DeltaXi;
            
            if (ResultingPoint[0] < -1.0)
            {
                ResultingPoint[0] = -1.0;
            }
            else if (ResultingPoint[0] > 1.0)
            {
                ResultingPoint[0] = 1.0;
            }
            
//             // Debug
//             KRATOS_WATCH(ResultingPoint[0])
//             KRATOS_WATCH(DeltaXi)
//             KRATOS_WATCH(OldDeltaXi)
            
            if ( std::abs(DeltaXi - OldDeltaXi) < Tolerance )
            {
                return true;
            }
        }
        
        return false;
    }
    
    /**
     * Project a point over a plane
     * @param PointOrigin: A point in the plane
     * @param PointDestiny: The point to be projected
     * @param Normal: The normal of the plane
     * @return PointProjected: The point pojected over the plane
     * @return Distance: The distance between the point and the plane
     */
    
    static inline void Project(
        const PointType& PointOrigin,
        const PointType& PointDestiny,
        PointType& PointProjected,
        double& Distance,
        const array_1d<double,3>& Normal
        )
    {
        array_1d<double,3> VectorPoints = PointDestiny.Coordinates() - PointOrigin.Coordinates();

        Distance = inner_prod(VectorPoints, Normal); 

        PointProjected.Coordinates() = PointDestiny.Coordinates() - Normal * Distance;
    }
    
    static inline PointType FastProject(
        const PointType& PointOrigin,
        const PointType& PointDestiny,
        const array_1d<double,3>& Normal
        )
    {
        array_1d<double,3> VectorPoints = PointDestiny.Coordinates() - PointOrigin.Coordinates();

        const double Distance = inner_prod(VectorPoints, Normal); 
        
        PointType PointProjected;
        PointProjected.Coordinates() = PointDestiny.Coordinates() - Normal * Distance;
        
        return PointProjected;
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
        const PointType PointOrig1,
        const PointType PointOrig2,
        const PointType PointOrig3
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
     * This function scales the points according to a factor (to increase the bounding box)
     * @param PointToScale: The point to scale
     * @param Center: The reference point
     * @param ScaleFactor: The factor considered to "grow" the node
     */
    
    template<class TPointType>
    static inline void ScaleNode(
        TPointType& PointToScale,
        const PointType& Center,
        const double ScaleFactor
        )
    {
        // We calculate the new distance
        const double Distance = ScaleFactor * DistancePoints(PointToScale.Coordinates(), Center.Coordinates());
        
        // Now the vector between nodes
        array_1d<double, 3> VectorPoints = PointToScale.Coordinates() - Center.Coordinates();
        VectorPoints /= norm_2(VectorPoints);
        
        // Finally we rescale
        PointToScale.Coordinates() = Center.Coordinates() + VectorPoints * Distance;
    }
    
    /**
     * Calculates the distance between nodes
     * @param PointOrigin: The first node
     * @param PointDestiny: The second node
     */
    
    static inline double DistancePoints(
        const GeometryType::CoordinatesArrayType& PointOrigin,
        const GeometryType::CoordinatesArrayType& PointDestiny
        )
    {
        return std::sqrt((PointOrigin[0] - PointDestiny[0]) * (PointOrigin[0] - PointDestiny[0])
                       + (PointOrigin[1] - PointDestiny[1]) * (PointOrigin[1] - PointDestiny[1])
                       + (PointOrigin[2] - PointDestiny[2]) * (PointOrigin[2] - PointDestiny[2]));
    }

    /**
     * This function calculates the center and radius of the geometry of a condition
     * @param Cond: The pointer to the condition of interest
     * @return Center: The center of the condition
     * @return Radius: The radius of the condition
     */

    static inline double CenterAndRadius(
        Condition::Pointer pCond,
        PointType& Center
        )
    {
        double Radius = 0.0;
        Center = pCond->GetGeometry().Center();
        
        // TODO: Add calculation of radius to geometry.h 
        array_1d<double, 3> AuxVector;
        for(unsigned int i = 0; i < pCond->GetGeometry().PointsNumber(); i++)
        {
            noalias(AuxVector) = Center.Coordinates() - pCond->GetGeometry()[i].Coordinates();;
            
            const double AuxRadius = inner_prod(AuxVector, AuxVector);

            if(AuxRadius > Radius)
            {
                Radius = AuxRadius;
            }
        }

        Radius = std::sqrt(Radius);
        
        return Radius;
    }

    /**
     * This function calculates the normal of a condition
     * @param Cond: The pointer to the condition of interest
     */

    static inline void ConditionNormal(Condition::Pointer pCond)
    {
        // TODO: Add calculation of normal to geometry.h 
        array_1d<double,3> & Normal     = pCond->GetValue(NORMAL);
        array_1d<double,3> & TangentXi  = pCond->GetValue(TANGENT_XI);
        array_1d<double,3> & TangentEta = pCond->GetValue(TANGENT_ETA);
        
        GeometryNormal(Normal, TangentXi, TangentEta, pCond->GetGeometry());
    }

    /**
     * This function calculates the normal of a condition
     * @param Cond: The pointer to the condition of interest
     */

    static inline array_1d<double,3> GaussPointNormal(
        const Vector N,
        const GeometryType & Geom
        )
    {
        array_1d<double,3> Normal = ZeroVector(3);
        for( unsigned int iNode = 0; iNode < Geom.PointsNumber(); ++iNode )
        {
            Normal += N[iNode] * Geom[iNode].GetValue(NORMAL); 
        }
        
        const double Tolerance = std::numeric_limits<double>::epsilon();
        
        if (norm_2(Normal) > Tolerance)
        {
            Normal = Normal/norm_2(Normal); // It is suppossed to be already unitary (just in case)
        }
        
        return Normal;
    }
    
    /**
     * This function calculates the normal of a geometry
     * @return Normal: The normal of the geometry
     * @param TangentXi: The first tangent of the geometry
     * @param TangentEta: The second tangent of the geometry
     * @param Geom: The geometry of interest
     */

    static inline void GeometryNormal(
        array_1d<double,3> & Normal,
        array_1d<double,3> & TangentXi,
        array_1d<double,3> & TangentEta,
        const GeometryType & Geom
        )
    {
        noalias(Normal) = ZeroVector(3);
        
        // Geom normal is the sum of all nodal normals
        // nodal normal = tangent_eta (or e3 in 2D) x tangent_xi
        for ( unsigned int i = 0; i < Geom.PointsNumber( ); ++i )
        {
            NodalTangents( TangentXi, TangentEta, Geom, i );
            Normal += MathUtils<double>::CrossProduct( TangentEta, TangentXi );
        }
        
        const double Tolerance = std::numeric_limits<double>::epsilon();
        
        if ( norm_2( Normal ) > Tolerance )
        {
            Normal     /= norm_2( Normal );
            TangentXi  /= norm_2( TangentXi );
            TangentEta /= norm_2( TangentEta );
        }
    }
        
    static inline void GeometryNormal(
        array_1d<double,3> & Normal,
        const GeometryType & Geom
        )
    {
        array_1d<double,3> TangentXi, TangentEta;
        
        GeometryNormal(Normal, TangentXi, TangentEta, Geom);
    }
    
    /**
     * This function calculates the tangents at a given node
     */

    static inline void NodalTangents(
        array_1d<double,3> & t1,    // tangent in xi direction
        array_1d<double,3> & t2,    // tangent in eta direction in 3D or simply e3 cartesian vector in 2D
        const GeometryType & Geom,
        const unsigned int i_node
        )
    {
        const unsigned int dimension = Geom.WorkingSpaceDimension( );
        const unsigned int local_dim = Geom.LocalSpaceDimension( );
        Matrix J_i = ZeroMatrix( dimension, local_dim ); 
        Geom.Jacobian( J_i, Geom[i_node].Coordinates( ) );

        if( dimension == 2 )
        {
            // t1 (axis-direction)
            t1[0] = J_i(0, 0);
            t1[1] = J_i(1, 0);
            t1[2] = 0.0;
            // t2 (z-direction)
            t2[0] = 0.0;
            t2[1] = 0.0;
            t2[2] = 1.0;
        }
        else // We can use just else, if it is not 2D it is 3D 
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                // Using the Jacobian tangent directions 
                t1[i] = J_i(i, 0);
                t2[i] = J_i(i, 1);
            }
        }
        
    }
    
    /**
     * It computes the mean of the normal in the condition in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the normal computed
     */
    
    static inline void ComputeNodesMeanNormalModelPart(ModelPart & rModelPart) 
    {
        // Tolerance
        const double Tolerance = std::numeric_limits<double>::epsilon();

        // Initialize normal vectors
        const array_1d<double,3> ZeroVect = ZeroVector(3);
        
        NodesArrayType& pNode = rModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->GetValue(NODAL_AREA)           = 0.0;
            noalias(itNode->GetValue(NORMAL))      = ZeroVect;
            noalias(itNode->GetValue(TANGENT_XI))  = ZeroVect; 
            noalias(itNode->GetValue(TANGENT_ETA)) = ZeroVect; 
        }
        
        // Sum all the nodes normals
        ConditionsArrayType& pCond = rModelPart.Conditions();
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            
            if (itCond->Is(SLAVE) || itCond->Is(MASTER) || itCond->Is(ACTIVE))
//             if (itCond->Is(ACTIVE) || itCond->Is(MASTER)) // NOTE: This can produce troubles, nodal normals are affected by the neighbour conditions (which can be innactive, think about search for neightbour conditions first)
            {
                ConditionNormal(*(itCond.base()));
                
                const unsigned int NumberNodes = itCond->GetGeometry().PointsNumber();
                const double & rArea = itCond->GetGeometry().Area()/NumberNodes;
                const array_1d<double, 3> & rNormal     = itCond->GetValue(NORMAL);
                const array_1d<double, 3> & rTangentXi  = itCond->GetValue(TANGENT_XI);
                const array_1d<double, 3> & rTangentEta = itCond->GetValue(TANGENT_ETA);
                
                for (unsigned int i = 0; i < NumberNodes; i++)
                {
                    #pragma omp atomic
                    itCond->GetGeometry()[i].GetValue(NODAL_AREA)             += rArea;
                    #pragma omp critical
                    noalias( itCond->GetGeometry()[i].GetValue(NORMAL) )      += rArea * rNormal;
                    #pragma omp critical
                    noalias( itCond->GetGeometry()[i].GetValue(TANGENT_XI) )  += rArea * rTangentXi;
                    #pragma omp critical
                    noalias( itCond->GetGeometry()[i].GetValue(TANGENT_ETA) ) += rArea * rTangentEta;
                }
            }
        }
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            const double TotalArea = itNode->GetValue(NODAL_AREA);
            if (TotalArea > Tolerance)
            {
                itNode->GetValue(NORMAL)      /= TotalArea;
                itNode->GetValue(TANGENT_XI)  /= TotalArea;
                itNode->GetValue(TANGENT_ETA) /= TotalArea;
            }
        }
        
//         // Applied laziness - MUST be calculated BEFORE normalizing the normals
//         ComputeDeltaNodesMeanNormalModelPart( rModelPart ); // NOTE: The area pondering?, after or before

        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            const double NormNormal     = norm_2(itNode->GetValue(NORMAL));
            const double NormTangentXi  = norm_2(itNode->GetValue(TANGENT_XI));
            const double NormTangentEta = norm_2(itNode->GetValue(TANGENT_ETA));
            
            if (NormNormal > Tolerance)
            {
                itNode->GetValue(NORMAL)      /= NormNormal;
                itNode->GetValue(TANGENT_XI)  /= NormTangentXi;
                itNode->GetValue(TANGENT_ETA) /= NormTangentEta;
            }
        }
    }

    /**
     * It computes the directional derivative of the normal in the condition in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the normal computed
     */

    static inline void ComputeDeltaNodesMeanNormalModelPart(ModelPart & rModelPart)
    {
        // TODO: Add parallelization

        // Conditions
        ConditionsArrayType& pCond  = rModelPart.Conditions();
        ConditionsArrayType::iterator it_cond_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_cond_end   = pCond.ptr_end();

        // Tolerance
        const double Tolerance = std::numeric_limits<double>::epsilon();

        // Initialize directional derivative
        const unsigned int dimension = it_cond_begin->WorkingSpaceDimension( ); 
        const Matrix ZeroDeltaNormal = ZeroMatrix( dimension, dimension );

        Matrix Delta_ne_adj  = Matrix(dimension, dimension);
        Matrix Ce = Matrix(dimension, dimension);
        
        const Matrix I = IdentityMatrix(dimension, dimension);

        NodesArrayType& pNode = rModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->GetValue(DELTA_NORMAL) = ZeroDeltaNormal;
        }
        
        // Sum the directional derivatives of the adjacent segments
        for(ConditionsArrayType::iterator itCond = it_cond_begin; itCond!=it_cond_end; itCond++)
        {
            if (itCond->Is(ACTIVE) || itCond->Is(MASTER)) 
            {
                const array_1d<double, 3> ne = itCond->GetValue(NORMAL);   // normalized condition normal
                Matrix ne_o_ne = subrange( outer_prod( ne, ne ), 0, dimension, 0, dimension );

                const unsigned int num_nodes = itCond->GetGeometry( ).PointsNumber();
                if( dimension == 2 )
                {
                    if (num_nodes == 2)
                    {
                        const double ne_norm = itCond->GetGeometry( ).Length( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                        
                        Delta_ne_adj( 0, 0 ) =  0.0;
                        Delta_ne_adj( 0, 1 ) = -1.0;
                        Delta_ne_adj( 1, 0 ) =  1.0;
                        Delta_ne_adj( 1, 1 ) =  0.0;
                        
                        Ce = prod( I - ne_o_ne, Delta_ne_adj ) / ne_norm;     // In 2D, Delta_ne_adj is node-independent => evaluated outside the nodes loop
                        
                        for (unsigned int i = 0; i < num_nodes; i++)
                        {
                            NodeType& node_j = itCond->GetGeometry( )[i];
                            
                            // -/+ 0.5 are the values of DN_Dxi for linear line elements at nodes 1 and 2 - no need to call the function
                            double DN_De_j = 0.0;
                            if( i == 0 )
                            {
                                DN_De_j = -0.5;
                            }
                            else if( i == 1 )
                            {
                                DN_De_j =  0.5;
                            }
                            node_j.GetValue(DELTA_NORMAL) += Ce * DN_De_j;
                        }
                    }
                    else
                    {
                        KRATOS_ERROR << "DELTA_NORMAL is not yet defined for higher order 1D elements. Number of nodes: " << num_nodes << std::endl;
                    }
                }
                else if ( dimension == 3 )
                {
                    const double ne_norm = itCond->GetGeometry( ).Area( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                    
                    for (unsigned int i = 0; i < num_nodes; i++)
                    {
                        Matrix J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
                        array_1d<double, 2> DN_De_j        = ZeroVector( 2 );  // Shape functions derivatives for node j [ DN_Dxi_j, DN_Deta_j ]
                        array_1d<double, 3> local_coords_j = ZeroVector( 3 );
                        
                        NodeType& node_j = itCond->GetGeometry( )[i];
                        
                        if( num_nodes == 3 )    // linear triangle element
                        {
                            if( i == 0 )
                            {
                                local_coords_j[0] = 0.0;
                                local_coords_j[1] = 0.0;
                                DN_De_j( 0 ) = -1.0;
                                DN_De_j( 1 ) = -1.0;
                            }
                            else if( i == 1 )
                            {
                                local_coords_j[0] = 1.0;
                                local_coords_j[1] = 0.0;
                                DN_De_j( 0 ) = 1.0;
                                DN_De_j( 1 ) = 0.0;
                            }
                            else // i == 2
                            {
                                local_coords_j[0] = 0.0;
                                local_coords_j[1] = 1.0;
                                DN_De_j( 0 ) = 0.0;
                                DN_De_j( 1 ) = 1.0;
                            }
                        }
                        else if( num_nodes == 4 )    // linear quad element - FIXME: this will screw things up if the user defines a 4-node tri elem
                        {
                            if( i == 0 )
                            {
                                local_coords_j[0] = -1.0;
                                local_coords_j[1] = -1.0;
                                DN_De_j( 0 ) = -0.5;
                                DN_De_j( 1 ) = -0.5;
                            }
                            else if( i == 1 )
                            {
                                local_coords_j[0] =  1.0;
                                local_coords_j[1] = -1.0;
                                DN_De_j( 0 ) =  0.5;
                                DN_De_j( 1 ) = -0.5;
                            }
                            else if( i == 2 )
                            {
                                local_coords_j[0] =  1.0;
                                local_coords_j[1] =  1.0;
                                DN_De_j( 0 ) =  0.5;
                                DN_De_j( 1 ) =  0.5;
                            }
                            else // i == 3
                            {
                                local_coords_j[0] = -1.0;
                                local_coords_j[1] =  1.0;
                                DN_De_j( 0 ) = -0.5;
                                DN_De_j( 1 ) =  0.5;
                            }
                        }
                        else
                        {
                            KRATOS_ERROR << "DELTA_NORMAL is not yet defined for higher order 2D elements. Number of nodes: " << itCond->GetGeometry( ).PointsNumber() << std::endl;
                        }
                        
                        itCond->GetGeometry( ).Jacobian( J, local_coords_j );
                        
                        /*
                         * NOTES:
                         * Delta_ne_adj here is Delta_n_hat_e in Popp's thesis equation A.3 in appendix A.
                         * In 2D, it is also like this, but I split it into a node-dependent DN_De_j and a node-independent Delta_ne_adj
                         * In 3D, Delta_ne_adj is node-dependent => evaluated completely inside the nodes loop
                         * In both cases, Ce is the elemental contribution in Popp's thesis equation A.2 in appendix A
                         * 
                         * ::::: DERIVATION OF Delta_ne_adj Matrix in 3D ::::: 
                         * 
                         *   [ SUM( N,xi * Delta_x ) x SUM( N,eta * x ) ] + [ SUM( N,xi * x ) x SUM( N,eta * Delta_x ) ]
                         * = [ SUM( N,xi * Delta_x ) x J_eta            ] + [            J_xi x SUM( N,eta * Delta_x ) ]
                         * = [ SUM( N,xi * Delta_x ) x J_eta            ] - [ SUM( N,eta * Delta_x ) x J_xi ]
                         * SUM( N,xi * Delta_x ) is the consistent assembly of N,xi in blocks. Similarily for N,eta
                         * therefore, for node j we only care about N,xi(j) and N,eta(j) nodal blocks
                         * = [ N,xi(j) * Delta_x(j) x J_eta ] - [ N,eta(j) * Delta_x(j) x J_xi ]
                         * = [ N,xi(j) * ( ones(3) x J_eta ) - N,eta(j) *( ones(3) x J_xi ) ] * Delta_x(j)
                         * = Delta_ne_adj * Delta_x
                         */
                        
                        Delta_ne_adj(0,0) = 0.0;
                        Delta_ne_adj(0,1) = +J(2,1) * DN_De_j(0) - J(2,0) * DN_De_j(1); 
                        Delta_ne_adj(0,2) = -J(1,1) * DN_De_j(0) + J(1,0) * DN_De_j(1); 
                        Delta_ne_adj(1,0) = -J(2,1) * DN_De_j(0) + J(2,0) * DN_De_j(1); 
                        Delta_ne_adj(1,1) = 0.0;                   
                        Delta_ne_adj(1,2) = +J(0,1) * DN_De_j(0) - J(0,0) * DN_De_j(1); 
                        Delta_ne_adj(2,0) = +J(1,1) * DN_De_j(0) - J(1,0) * DN_De_j(1); 
                        Delta_ne_adj(2,1) = -J(0,1) * DN_De_j(0) + J(0,0) * DN_De_j(1); 
                        Delta_ne_adj(2,2) = 0.0;
                        
                        Ce = prod( I - ne_o_ne, Delta_ne_adj ) / ne_norm;
                        node_j.GetValue(DELTA_NORMAL) += Ce;
                    }
                }
                else
                {
                    KRATOS_ERROR << "Bad dimension provided to calculate DELTA_NORMAL. Dimension = " << dimension << std::endl;
                }
            }
        }
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            const array_1d<double, 3> & nj = itNode->GetValue(NORMAL); // nodal non-normalized normal (this function is called before normalization)
            
            Matrix nj_o_nj = subrange( outer_prod( nj, nj ), 0, dimension, 0, dimension );
            const double nj_norm = norm_2( nj );
            const double nj_norm_3 = nj_norm * nj_norm * nj_norm;
            
            if ( nj_norm_3 > Tolerance )
            {
                const Matrix Cj = I / nj_norm - nj_o_nj / nj_norm_3;
                itNode->GetValue(DELTA_NORMAL) = prod( Cj, itNode->GetValue(DELTA_NORMAL) );
            }
        }
    }
    
    /**
     * It changes from active to inactive and viceversa the nodes 
     * @param ModelPart: The model part to compute
     * @return The modelparts with the conditions changed
     */
    
    static inline void ReComputeActiveInactive(ModelPart & rModelPart)
    {
        // TODO: If works make it parallell (it is difficult, be careful with the repeated nodes) 
       
        ConditionsArrayType& pConditions = rModelPart.GetSubModelPart("Contact").Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();

//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            bool ConditionIsActive = false; // It is supposed to be always defined, and with this only the slave conditions will be taken into account
            if( (itCond)->IsDefined(ACTIVE) == true)
            {
                ConditionIsActive = (itCond)->Is(ACTIVE);
            }
            
            if ( ConditionIsActive == true )
            {
                // Recompute Active/Inactive nodes
                double cn = 0.0;
                if (itCond->GetProperties().Has(NORMAL_AUGMENTATION_FACTOR) == true)
                {
                    cn = itCond->GetProperties().GetValue(NORMAL_AUGMENTATION_FACTOR); 
                }
                
                double ct = 0.0;
                if (itCond->GetProperties().Has(TANGENT_AUGMENTATION_FACTOR) == true)
                {
                    ct = itCond->GetProperties().GetValue(TANGENT_AUGMENTATION_FACTOR); 
                }

                GeometryType & CondGeometry = itCond->GetGeometry();
                
                for(unsigned int itNode = 0; itNode!=CondGeometry.PointsNumber(); itNode++)
                {
                    if (CondGeometry[itNode].Is(VISITED) == false)
                    {
                        const double mu = itCond->GetProperties().GetValue(FRICTION_COEFFICIENT); 
                        const double gn = CondGeometry[itNode].GetValue(WEIGHTED_GAP);
                        const array_1d<double,3> LagrangeMultiplier = CondGeometry[itNode].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                        const array_1d<double,3> NodalNormal = CondGeometry[itNode].GetValue(NORMAL); 
                        
//                         const double LambdaNormal = inner_prod(LagrangeMultiplier, NodalNormal);
//                         const double AugmentedNormalPressure = LambdaNormal + cn * gn;           
                        
                        double AugmentedNormalPressure = inner_prod(LagrangeMultiplier, NodalNormal);
                        if (gn < 0.0) // NOTE: Penetration
                        {
                            AugmentedNormalPressure += cn * gn;     
                        }
                        
                        if (AugmentedNormalPressure < 0.0) // NOTE: This could be conflictive (< or <=)
                        {
                            CondGeometry[itNode].Set(ACTIVE, true);
                        }
                        else
                        {
                            CondGeometry[itNode].Set(ACTIVE, false);
//                             CondGeometry[itNode].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0) = ZeroVector(3);
                        }
                        
                        const array_1d<double, 3> NodalTangentXi  = CondGeometry[itNode].GetValue(TANGENT_XI); 
                        const array_1d<double, 3> NodalTangentEta = CondGeometry[itNode].GetValue(TANGENT_ETA);
                        const double TangentXiLM  = inner_prod(NodalTangentXi,  LagrangeMultiplier);
                        const double TangentEtaLM = inner_prod(NodalTangentEta, LagrangeMultiplier);
                        const double LambdaTangent = std::sqrt(TangentXiLM * TangentXiLM + TangentEtaLM * TangentEtaLM); 
                        
                        const double AugmentedTangentPressure = std::abs(LambdaTangent + ct * CondGeometry[itNode].GetValue(WEIGHTED_SLIP)) + mu * AugmentedNormalPressure;
                        
                        if (AugmentedTangentPressure < 0.0) // TODO: Check if it is minor equal or just minor
                        {
                            CondGeometry[itNode].Set(SLIP, false);
                        }
                        else
                        {
                            CondGeometry[itNode].Set(SLIP, true);
                        }
                        
                        CondGeometry[itNode].Set(VISITED, true);
                    }
                }
            }
        }
    }
    
    /**
     * It changes from active to inactive and viceversa the nodes 
     * @param ModelPart: The model part to compute
     * @return The modelparts with the conditions changed
     */
    
    static inline void ReComputeActiveInactiveALMFrictionless(ModelPart & rModelPart)
    {
        // TODO: If works make it parallell (it is difficult, be careful with the repeated nodes) 
       
        ConditionsArrayType& pConditions = rModelPart.GetSubModelPart("Contact").Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();

//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            bool ConditionIsActive = false; // It is supposed to be always defined, and with this only the slave conditions will be taken into account
            if( (itCond)->IsDefined(ACTIVE) == true)
            {
                ConditionIsActive = (itCond)->Is(ACTIVE);
            }
            
            if ( ConditionIsActive == true )
            {
                // Recompute Active/Inactive nodes
                double epsilon = 0.0;
                if (itCond->GetProperties().Has(PENALTY_FACTOR) == true)
                {
                    epsilon = itCond->GetProperties().GetValue(PENALTY_FACTOR); 
                }
                
                double k = 1.0;
                if (itCond->GetProperties().Has(SCALE_FACTOR) == true)
                {
                    k = itCond->GetProperties().GetValue(SCALE_FACTOR); 
                }

                bool DeactivateCondition = true; // If all nodes are inactive we deactivate the condition
                GeometryType & CondGeometry = itCond->GetGeometry();
                
                for(unsigned int itNode = 0; itNode!=CondGeometry.PointsNumber(); itNode++)
                {
                    if (CondGeometry[itNode].Is(VISITED) == false)
                    {
                        const double gn = CondGeometry[itNode].GetValue(WEIGHTED_GAP);
                        const double AugmentedNormalPressure = k * CondGeometry[itNode].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) + epsilon * gn;     
                        
                        if (AugmentedNormalPressure < 0.0) // NOTE: This could be conflictive (< or <=)
                        {
                            CondGeometry[itNode].Set(ACTIVE, true);
                            DeactivateCondition = false; // The condition is active
                        }
                        else
                        {
                            CondGeometry[itNode].Set(ACTIVE, false);
//                             CondGeometry[itNode].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS, 0) = 0.0;
                        }
                        
//                         // Debug 
//                         std::cout << CondGeometry[itNode].Id() << " Gap: " << gn  << " Pressure: " << AugmentedNormalPressure << " Active: " << CondGeometry[itNode].Is(ACTIVE) << std::endl;
                        
                        CondGeometry[itNode].Set(VISITED, true);
                    }
                }
                
                // We deactivate the condition if necessary
                if (DeactivateCondition == true)
                {
                    itCond->Set(ACTIVE, false);
                }
            }
        }
    }
    
    /**
     * It resets the visited status in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the nodes changed
     */
    
    static inline void ResetVisited(ModelPart & rModelPart)
    {
        NodesArrayType& pNode = rModelPart.GetSubModelPart("Contact").Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            itNode->Set(VISITED, false);
        }
    }
    
    /**
     * It resets the value of the weighted gap and slip 
     * @param ModelPart: The model part to compute
     * @return The modelparts with the nodes changed
     */
    
    static inline void ResetWeightedValues(ModelPart & rModelPart)
    {
        NodesArrayType& pNode = rModelPart.GetSubModelPart("Contact").Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            if (itNode->Is(ACTIVE) == false)
            {
                itNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0) = ZeroVector(3);
            }
            itNode->GetValue(WEIGHTED_GAP)      = 0.0;
        }
    }
    
    /**
     * It resets the value of the weighted gap 
     * @param ModelPart: The model part to compute
     * @return The modelparts with the nodes changed
     */
    
    static inline void ResetWeightedALMFrictionlessValues(ModelPart & rModelPart)
    {
        NodesArrayType& pNode = rModelPart.GetSubModelPart("Contact").Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            if (itNode->Is(ACTIVE) == false)
            {
                itNode->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS, 0) = 0.0;
            }
            
            itNode->GetValue(WEIGHTED_GAP) = 0.0;
        }
    }
    
    /**
     * It calculates the matrix of coordinates of a geometry
     * @param nodes: The geometry to calculate
     * @param current: If we calculate the current coordinates or the initial ones
     * @return Coordinates: The matrix containing the coordinates of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> GetCoordinates(
        const GeometryType& nodes,
        const bool current = true
        )
    {
        /* DEFINITIONS */            
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> Coordinates;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            array_1d<double, 3> coord = nodes[iNode].Coordinates();
            
            if (current == false)
            {
                coord -= nodes[iNode].FastGetSolutionStepValue(DISPLACEMENT, 0);
            }

            for (unsigned int iDof = 0; iDof < TDim; iDof++)
            {
                Coordinates(iNode, iDof) = coord[iDof];
            }
        }
        
        return Coordinates;
    }

    /**
     * It calculates the vector of a variable of a geometry
     * @param nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return VarVector: The vector containing the variables of the geometry
     */
    
    template< unsigned int TNumNodes >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& nodes,
        const Variable<double>& rVarName,
        unsigned int step
        )
    {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> VarVector;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            VarVector[iNode] = nodes[iNode].FastGetSolutionStepValue(rVarName, step);
        }
        
        return VarVector;
    }
    
    template< unsigned int TNumNodes >
    static inline boost::numeric::ublas::bounded_matrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& nodes,
        const Variable<double>& rVarName,
        unsigned int step
        )
    {
        /* DEFINITIONS */        
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, 1> VarVector;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            VarVector(iNode, 0) = nodes[iNode].FastGetSolutionStepValue(rVarName, step);
        }
        
        return VarVector;
    }

    template< unsigned int TNumNodes >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& nodes,
        const Variable<double>& rVarName
        )
    {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> VarVector;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            VarVector[iNode] = nodes[iNode].GetValue(rVarName);
        }
        
        return VarVector;
    }
    
    template< unsigned int TNumNodes >
    static inline boost::numeric::ublas::bounded_matrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& nodes,
        const Variable<double>& rVarName
        )
    {
        /* DEFINITIONS */        
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, 1> VarVector;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            VarVector(iNode, 0) = nodes[iNode].GetValue(rVarName);
        }
        
        return VarVector;
    }
    
    /**
     * It calculates the matrix of a variable of a geometry
     * @param nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return VarMatrix: The matrix containing the variables of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& nodes,
        const Variable<array_1d<double,3> >& rVarName,
        unsigned int step
        )
    {
        /* DEFINITIONS */        
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> VarMatrix;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            const array_1d<double, 3> Value = nodes[iNode].FastGetSolutionStepValue(rVarName, step);
            for (unsigned int iDof = 0; iDof < TDim; iDof++)
            {
                VarMatrix(iNode, iDof) = Value[iDof];
            }
        }
        
        return VarMatrix;
    }

    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& nodes,
        const Variable<array_1d<double,3> >& rVarName
        )
    {
        /* DEFINITIONS */        
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> VarMatrix;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            const array_1d<double, 3> Value = nodes[iNode].GetValue(rVarName);
            for (unsigned int iDof = 0; iDof < TDim; iDof++)
            {
                VarMatrix(iNode, iDof) = Value[iDof];
            }
        }
        
        return VarMatrix;
    }
    
private:
};// class ContactUtilities
}
#endif /* KRATOS_CONTACT_UTILITIES defined */
