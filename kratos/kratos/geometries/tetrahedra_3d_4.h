//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//



#if !defined(KRATOS_TETRAHEDRA_3D_4_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_3D_4_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/triangle_3d_3.h"
#include "integration/tetrahedron_gauss_legendre_integration_points.h"

#include "geometries/plane.h"

namespace Kratos
{
/**
 * An eight node hexahedra geometry with linear shape functions
 */
template<class TPointType> class Tetrahedra3D4 : public Geometry<TPointType>
{
public:
    /**
     * Type Definitions
     */

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge and face geometries
     */
    typedef Line3D2<TPointType> EdgeType;
    typedef Triangle3D3<TPointType> FaceType;

    /**
     * Pointer definition of Tetrahedra3D4
     */
    KRATOS_CLASS_POINTER_DEFINITION(Tetrahedra3D4);

    /**
     * Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /**
     * A Vector of counted pointers to Geometries. Used for
     * returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /**
     * Type used for indexing in geometry class.std::size_t used for indexing
     * point or integration point access methods and also all other
     * methods which need point or integration point index.
     */
    typedef typename BaseType::IndexType IndexType;


    /**
     * This typed used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Array of counted pointers to point. This type used to hold
     * geometry's points.
     */
    typedef typename BaseType::PointsArrayType PointsArrayType;

    /**
     * This type used for representing an integration point in
     * geometry. This integration point is a point with an
     * additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A Vector of IntegrationPointType which used to hold
     * integration points related to an integration
     * method. IntegrationPoints functions used this type to return
     * their results.
     */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * A Vector of IntegrationPointsArrayType which used to hold
     * integration points related to different integration method
     * implemented in geometry.
     */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /**
     * A third order tensor used as shape functions' values
     * container.
     */
    typedef typename BaseType::ShapeFunctionsValuesContainerType
    ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType
    ShapeFunctionsLocalGradientsContainerType;

    /**
     * A third order tensor to hold jacobian matrices evaluated at
     * integration points. Jacobian and InverseOfJacobian functions
     * return this type as their result.
     */
    typedef typename BaseType::JacobiansType JacobiansType;

    /**
     * A third order tensor to hold shape functions' local
     * gradients. ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * Type of Matrix
     */
    typedef Matrix MatrixType;


    /**
     * Life Cycle
     */

    Tetrahedra3D4( const PointType& Point1, const PointType& Point2,
                   const PointType& Point3, const PointType& Point4 )
        : BaseType(PointsArrayType(), &msGeometryData)
    {
        this->Points().reserve(4);
        this->Points().push_back(typename PointType::Pointer(new PointType(Point1)));
        this->Points().push_back(typename PointType::Pointer(new PointType(Point2)));
        this->Points().push_back(typename PointType::Pointer(new PointType(Point3)));
        this->Points().push_back(typename PointType::Pointer(new PointType(Point4)));
    }

    Tetrahedra3D4( typename PointType::Pointer pPoint1,
                   typename PointType::Pointer pPoint2,
                   typename PointType::Pointer pPoint3,
                   typename PointType::Pointer pPoint4 )
        : BaseType(PointsArrayType(), &msGeometryData)
    {
        this->Points().reserve(4);
        this->Points().push_back(pPoint1);
        this->Points().push_back(pPoint2);
        this->Points().push_back(pPoint3);
        this->Points().push_back(pPoint4);
    }

    Tetrahedra3D4( const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
        if( this->PointsNumber() != 4)
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "Invalid points number. Expected 4, given " ,
                               this->PointsNumber());
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Tetrahedra3D4(Tetrahedra3D4 const& rOther)
        : BaseType(rOther)
    {
    }

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Tetrahedra3D4(Tetrahedra3D4<TOtherPointType> const& rOther)
        : BaseType(rOther)
    {
    }

    /// Destructor. Does nothing!!!
    virtual ~Tetrahedra3D4() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_Tetrahedra;
    }
    GeometryData::KratosGeometryType GetGeometryType()
    {
        return GeometryData::Kratos_Tetrahedra3D4;
    }

    /**
     * Operators
     */

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    Tetrahedra3D4& operator=(const Tetrahedra3D4& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    Tetrahedra3D4& operator=(Tetrahedra3D4<TOtherPointType> const & rOther)
    {
        BaseType::operator=(rOther);

        return *this;
    }


    /**
     * Operations
     */

    typename BaseType::Pointer Create(PointsArrayType const& ThisPoints) const
    {
        return typename BaseType::Pointer(new Tetrahedra3D4(ThisPoints));
    }

    virtual Geometry< Point<3> >::Pointer Clone() const
    {
        Geometry< Point<3> >::PointsArrayType NewPoints;
        //making a copy of the nodes TO POINTS (not Nodes!!!)
        for(IndexType i = 0 ; i < this->Points().size() ; i++)
            NewPoints.push_back(this->Points()[i]);
        //creating a geometry with the new points
        Geometry< Point<3> >::Pointer
        p_clone(new Tetrahedra3D4< Point<3> >(NewPoints));
        p_clone->ClonePoints();

        return p_clone;
    }

    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors(Vector& rResult) const
    {
        if(rResult.size() != 4)
            rResult.resize(4, false);
        std::fill(rResult.begin(), rResult.end(), 1.00 / 4.00);
        return rResult;
    }


    /**
     * Informations
     */

    /**
     * This method calculates and returns Length or charactereistic
     * length of this geometry depending on it's dimension. For one
     * dimensional geometry for example Line it returns length of it
     * and for the other geometries it gives Characteristic length
     * otherwise.
     *
     * @return double value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     *
     * :TODO: might be necessary to reimplement
     */
    virtual double Length() const
    {
        const double param = 2.0396489026555;  //12/raiz(2);
        return  param * pow(Volume(), 0.33333333333333); //sqrt(fabs( DeterminantOfJacobian(PointType())));
    }

    /**
     * This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     *
     * @return double value contains area or surface area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     *
     * :TODO: might be necessary to reimplement
     */
    virtual double Area() const
    {
        return this->Volume();
    }


    /**
     * This method calculate and return length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     *
     * :TODO: might be necessary to reimplement
     */

    virtual double Volume() const //closed formula for the linear triangle
    {
        const double onesixth = 1.0/6.0;

        const CoordinatesArrayType& rP0 = this->Points()[0].Coordinates();
        const CoordinatesArrayType& rP1 = this->Points()[1].Coordinates();
        const CoordinatesArrayType& rP2 = this->Points()[2].Coordinates();
        const CoordinatesArrayType& rP3 = this->Points()[3].Coordinates();

        double x10 = rP1[0] - rP0[0];
        double y10 = rP1[1] - rP0[1];
        double z10 = rP1[2] - rP0[2];

        double x20 = rP2[0] - rP0[0];
        double y20 = rP2[1] - rP0[1];
        double z20 = rP2[2] - rP0[2];

        double x30 = rP3[0] - rP0[0];
        double y30 = rP3[1] - rP0[1];
        double z30 = rP3[2] - rP0[2];

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return  detJ*onesixth;
    }


    virtual double DomainSize() const
    {

        return Volume();
    }


    /**
    * Returns a matrix of the local coordinates of all points
    * @param rResult a Matrix that will be overwritten by the results
    * @return the coordinates of all points of the current geometry
    */
    virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
    {
        if(rResult.size1()!= 4 || rResult.size2()!= 3)
            rResult.resize(4,3);

        rResult(0,0)=0.0;
        rResult(0,1)=0.0;
        rResult(0,2)=0.0;
        rResult(1,0)=1.0;
        rResult(1,1)=0.0;
        rResult(1,2)=0.0;
        rResult(2,0)=0.0;
        rResult(2,1)=1.0;
        rResult(2,2)=0.0;
        rResult(3,0)=0.0;
        rResult(3,1)=0.0;
        rResult(3,2)=1.0;

        return rResult;
    }
    /**
     * Returns whether given arbitrary point is inside the Geometry
     */
    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        this->PointLocalCoordinates( rResult, rPoint );
        if( rResult[0] >= 0.0-1.0e-8 && rResult[0] <= 1.0+1.0e-8 )
            if( rResult[1] >= 0.0-1.0e-8 && rResult[1] <= 1.0 +1.0e-8)
                if( rResult[2] >= 0.0-1.0e-8 && rResult[2] <= 1.0+1.0e-8 )
                    if( ((1.0-(rResult[0] + rResult[1] + rResult[2])) >= 0.0-1.0e-8)&&((1.0-(rResult[0] + rResult[1] + rResult[2])) <= 1.0+1.0e-8) )
                        return true;
        return false;
    }



    /** This method gives you number of all edges of this
    geometry.
    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
     */
    // will be used by refinement algorithm, thus uncommented. janosch.
    virtual SizeType EdgesNumber() const
    {
        return 6;
    }

    virtual SizeType FacesNumber() const
    {
        return 4;
    }

    /** This method gives you all edges of this geometry.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
     */
    virtual GeometriesArrayType Edges(void)
    {
        GeometriesArrayType edges = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer EdgePointerType;
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(1))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(1),
                                             this->pGetPoint(2))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(2),
                                             this->pGetPoint(0))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(3))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(1),
                                             this->pGetPoint(3))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(2),
                                             this->pGetPoint(3))) );

        return edges;
    }

    virtual GeometriesArrayType Faces(void)
    {
        GeometriesArrayType faces = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer FacePointerType;
        faces.push_back( FacePointerType(new FaceType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(2),
                                             this->pGetPoint(1))) );
        faces.push_back( FacePointerType(new FaceType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(3),
                                             this->pGetPoint(2))) );
        faces.push_back( FacePointerType(new FaceType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(1),
                                             this->pGetPoint(3))) );
        faces.push_back( FacePointerType(new FaceType(
                                             this->pGetPoint(2),
                                             this->pGetPoint(3),
                                             this->pGetPoint(1))) );
        return faces;
    }


    //Connectivities of faces required
    virtual void NumberNodesInFaces (boost::numeric::ublas::vector<unsigned int>& NumberNodesInFaces) const
    {
        NumberNodesInFaces.resize(4);
        // Linear Tetrahedra have elements of 3 nodes as faces
        NumberNodesInFaces[0]=3;
        NumberNodesInFaces[1]=3;
        NumberNodesInFaces[2]=3;
        NumberNodesInFaces[4]=3;

    }

    virtual void NodesInFaces (boost::numeric::ublas::matrix<unsigned int>& NodesInFaces) const
    {
        NodesInFaces.resize(4,4);
        NodesInFaces(0,0)=1;
        NodesInFaces(1,0)=2;
        NodesInFaces(2,0)=3;
        NodesInFaces(3,0)=0;//other node

        NodesInFaces(0,1)=2;
        NodesInFaces(1,1)=0;
        NodesInFaces(2,1)=3;
        NodesInFaces(3,1)=1;//other node

        NodesInFaces(0,2)=0;
        NodesInFaces(1,2)=1;
        NodesInFaces(2,2)=3;
        NodesInFaces(3,2)=2;//other node

        NodesInFaces(0,3)=0;
        NodesInFaces(1,3)=2;
        NodesInFaces(2,3)=1;
        NodesInFaces(3,3)=3;//other node
    }



    /**
     * Shape Function
     */

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     *
     * @return the value of the shape function at the given point
     * TODO: TO BE VERIFIED
     */
    virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint) const
    {
        switch( ShapeFunctionIndex )
        {
        case 0:
            return( 1.0-(rPoint[0]+rPoint[1]+rPoint[2]) );
        case 1:
            return( rPoint[0] );
        case 2:
            return( rPoint[1] );
        case 3:
            return( rPoint[2] );
        default:
            KRATOS_THROW_ERROR(std::logic_error,
                               "Wrong index of shape function!" , *this);
        }
        return 0;
    }
    
    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
    {
        if(rResult.size1() != this->PointsNumber() || rResult.size2() != this->LocalSpaceDimension())
            rResult.resize(this->PointsNumber(),this->LocalSpaceDimension(),false);
        
        CalculateShapeFunctionsLocalGradients(rResult, rPoint);
        
        return rResult;
    }

    /**
     * Calculates the Gradients of the shape functions.
     * Calculates the gradients of the shape functions with regard to the global
     * coordinates in all
     * integration points (\f$ \frac{\partial N^i}{\partial X_j} \f$)
     *
     * @param rResult a container which takes the calculated gradients
     * @param ThisMethod the given IntegrationMethod
     * @return the gradients of all shape functions with regard to the global coordinates
     *
     * KLUDGE: method call only works with explicit JacobiansType rather than creating
     * JacobiansType within argument list
     *
     * :TODO: TESTING!!!
     */
    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod) const
    {
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber(ThisMethod);
        if(integration_points_number == 0)
            KRATOS_THROW_ERROR(std::logic_error,
                               "This integration method is not supported" , *this);

        boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
        double x10 = this->Points()[1].X() - this->Points()[0].X();
        double y10 = this->Points()[1].Y() - this->Points()[0].Y();
        double z10 = this->Points()[1].Z() - this->Points()[0].Z();

        double x20 = this->Points()[2].X() - this->Points()[0].X();
        double y20 = this->Points()[2].Y() - this->Points()[0].Y();
        double z20 = this->Points()[2].Z() - this->Points()[0].Z();

        double x30 = this->Points()[3].X() - this->Points()[0].X();
        double y30 = this->Points()[3].Y() - this->Points()[0].Y();
        double z30 = this->Points()[3].Z() - this->Points()[0].Z();

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

        DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
        DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
        DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
        DN_DX(1,0) = y20 * z30 - y30 * z20;
        DN_DX(1,1) = z20 * x30 - x20 * z30;
        DN_DX(1,2) = x20 * y30 - y20 * x30;
        DN_DX(2,0) = -y10 * z30 + z10 * y30;
        DN_DX(2,1) = x10 * z30 - z10 * x30;
        DN_DX(2,2) = -x10 * y30 + y10 * x30;
        DN_DX(3,0) = y10 * z20 - z10 * y20;
        DN_DX(3,1) = -x10 * z20 + z10 * x20;
        DN_DX(3,2) = x10 * y20 - y10 * x20;

        DN_DX /= detJ;
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }

        for(unsigned int i=0; i<integration_points_number; i++)
                rResult[i] = DN_DX;
        

        return rResult;
    }

    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult
        , Vector& determinants_of_jacobian
        , IntegrationMethod ThisMethod) const
    {
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber(ThisMethod);
        if(integration_points_number == 0)
            KRATOS_THROW_ERROR(std::logic_error,
                               "This integration method is not supported" , *this);

        boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
        double x10 = this->Points()[1].X() - this->Points()[0].X();
        double y10 = this->Points()[1].Y() - this->Points()[0].Y();
        double z10 = this->Points()[1].Z() - this->Points()[0].Z();

        double x20 = this->Points()[2].X() - this->Points()[0].X();
        double y20 = this->Points()[2].Y() - this->Points()[0].Y();
        double z20 = this->Points()[2].Z() - this->Points()[0].Z();

        double x30 = this->Points()[3].X() - this->Points()[0].X();
        double y30 = this->Points()[3].Y() - this->Points()[0].Y();
        double z30 = this->Points()[3].Z() - this->Points()[0].Z();

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

        DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
        DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
        DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
        DN_DX(1,0) = y20 * z30 - y30 * z20;
        DN_DX(1,1) = z20 * x30 - x20 * z30;
        DN_DX(1,2) = x20 * y30 - y20 * x30;
        DN_DX(2,0) = -y10 * z30 + z10 * y30;
        DN_DX(2,1) = x10 * z30 - z10 * x30;
        DN_DX(2,2) = -x10 * y30 + y10 * x30;
        DN_DX(3,0) = y10 * z20 - z10 * y20;
        DN_DX(3,1) = -x10 * z20 + z10 * x20;
        DN_DX(3,2) = x10 * y20 - y10 * x20;

        DN_DX /= detJ;

        if(determinants_of_jacobian.size() != integration_points_number )
            determinants_of_jacobian.resize(integration_points_number,false);

        for(unsigned int i=0; i<integration_points_number; i++)
            determinants_of_jacobian[i] = detJ;

//        Volume = detJ*0.1666666666666666666667;

        //workaround by riccardo
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }
        for(unsigned int i=0; i<integration_points_number; i++)
                rResult[i] = DN_DX;
        

        return rResult;
    }




    /// detect if two tetrahedra are intersected
    virtual bool HasIntersection( const BaseType& rThisGeometry)
    {

        array_1d<Plane, 4>  plane;
        std::vector<BaseType> Intersection;

        //const BaseType& geom_1 = *this;
        const BaseType& geom_2 = rThisGeometry;

        GetPlanes(plane);
        Intersection.push_back(geom_2);
        for (unsigned int i = 0; i < 4; ++i)
        {
            std::vector<BaseType> inside;
            for (unsigned int j = 0; j < Intersection.size(); ++j)
            {
                SplitAndDecompose(Intersection[j], plane[i], inside);
            }
            Intersection = inside;
        }

        return bool (Intersection.size() > 0);
    }


    virtual bool HasIntersection(const Point<3, double>& rLowPoint, const Point<3, double>& rHighPoint)
    {
        return true;
    }


    void SplitAndDecompose(
        const BaseType& tetra, Plane& plane,
        std::vector<BaseType>& inside)
    {

        // Determine on which side of the plane the points of the tetrahedron lie.
        int i = 0;
        int positive = 0, negative = 0, zero = 0;
        array_1d<double,4> C;
        array_1d<int,4>  pos, neg, zer;

        noalias(C)   = ZeroVector(4);

        pos[0] = 0;
        neg[0] = 0;
        zer[0] = 0;
        pos[1] = 0;
        neg[1] = 0;
        zer[1] = 0;
        pos[2] = 0;
        neg[2] = 0;
        zer[2] = 0;
        pos[3] = 0;
        neg[3] = 0;
        zer[3] = 0;

        for (i = 0; i < 4; ++i)
        {
            const array_1d<double,3>& p = tetra[i].Coordinates();
            C[i] = plane.DistanceTo(p);
            if (C[i] > 0.00)
                pos[positive++] = i;
            else if (C[i] < 0.00)
                neg[negative++] = i;
            else
                zer[zero++] = i;
        }


        // For a split to occur, one of the c_i must be positive and one must
        // be negative.
        if (negative == 0)
        {
            // Tetrahedron is completely on the positive side of plane, full clip.
            return;
        }

        if (positive == 0)
        {
            // Tetrahedron is completely on the negative side of plane.
            inside.push_back(tetra);
            return;
        }

        double w0, w1, invCDiff;
        array_1d<array_1d<double, 3 >, 4> intp;
        array_1d<array_1d<double, 3 >, 4> V;

        if (positive == 3)
        {
            // +++-
            for (i = 0; i < positive; ++i)
            {
                invCDiff = 1.00/(C[pos[i]] - C[neg[0]]);
                w0        = -C[neg[0]]*invCDiff;
                w1        = +C[pos[i]]*invCDiff;
                V[pos[i]] = w0*tetra[pos[i]].Coordinates() +  w1*tetra[neg[0]].Coordinates();
            }
            inside.push_back(tetra);
        }

        else if (positive == 2)
        {
            if (negative == 2)
            {
                // ++--
                for (i = 0; i < positive; ++i)
                {
                    invCDiff = (1.00)/(C[pos[i]] - C[neg[0]]);
                    w0 = -C[neg[0]]*invCDiff;
                    w1 = +C[pos[i]]*invCDiff;
                    intp[i] = w0*tetra[pos[i]].Coordinates() + w1*tetra[neg[0]].Coordinates();
                }
                for (i = 0; i < negative; ++i)
                {
                    invCDiff = (1.00)/(C[pos[i]] - C[neg[1]]);
                    w0 = -C[neg[1]]*invCDiff;
                    w1 = +C[pos[i]]*invCDiff;
                    intp[i+2] = w0*tetra[pos[i]].Coordinates() + w1*tetra[neg[1]].Coordinates();
                }

                V[pos[0]] = intp[2];
                V[pos[1]] = intp[1];
                inside.push_back(tetra);
            }
            else
            {
                // ++-0
                for (i = 0; i < positive; ++i)
                {
                    invCDiff = (1.00)/(C[pos[i]] - C[neg[0]]);
                    w0 = -C[neg[0]]*invCDiff;
                    w1 = +C[pos[i]]*invCDiff;
                    V[pos[i]] = w0*tetra[pos[i]].Coordinates() + w1*tetra[neg[0]].Coordinates();
                }
                inside.push_back(tetra);
            }
        }


        else if (positive == 1)
        {
            if (negative == 3)
            {
                // +---
                for (i = 0; i < negative; ++i)
                {
                    invCDiff = (1.00)/(C[pos[0]] - C[neg[i]]);
                    w0       = -C[neg[i]]*invCDiff;
                    w1       = +C[pos[0]]*invCDiff;
                    intp[i]  = w0*tetra[pos[0]] + w1*tetra[neg[i]];
                }

                V[pos[0]] = intp[0];
                inside.push_back(tetra);
            }
            else if (negative == 2)
            {
                // +--0
                for (i = 0; i < negative; ++i)
                {
                    invCDiff = (1.00)/(C[pos[0]] - C[neg[i]]);
                    w0 = -C[neg[i]]*invCDiff;
                    w1 = +C[pos[0]]*invCDiff;
                    intp[i] = w0*tetra[pos[0]] + w1*tetra[neg[i]];
                }

                V[pos[0]] = intp[0];
                inside.push_back(tetra);
            }
            else
            {
                // +-00
                invCDiff        = (1.00)/(C[pos[0]] - C[neg[0]]);
                w0              = -C[neg[0]]*invCDiff;
                w1              = +C[pos[0]]*invCDiff;
                V[pos[0]]       = w0*tetra[pos[0]] + w1*tetra[neg[0]];
                inside.push_back(tetra);
            }
        }
    }


    void GetPlanes(array_1d<Plane, 4>& plane)
    {
        const BaseType& geom_1 = *this;
        array_1d<double, 3> edge10 = geom_1[1].Coordinates() - geom_1[0].Coordinates();
        array_1d<double, 3> edge20 = geom_1[2].Coordinates() - geom_1[0].Coordinates();
        array_1d<double, 3> edge30 = geom_1[3].Coordinates() - geom_1[0].Coordinates();
        array_1d<double, 3> edge21 = geom_1[2].Coordinates() - geom_1[1].Coordinates();
        array_1d<double, 3> edge31 = geom_1[3].Coordinates() - geom_1[1].Coordinates();


        plane[0].mNormal = MathUtils<double>::UnitCrossProduct(edge20,edge10);  // <v0,v2,v1>
        plane[1].mNormal = MathUtils<double>::UnitCrossProduct(edge10,edge30);  // <v0,v1,v3>
        plane[2].mNormal = MathUtils<double>::UnitCrossProduct(edge30,edge20);  // <v0,v3,v2>
        plane[3].mNormal = MathUtils<double>::UnitCrossProduct(edge21,edge31);  // <v1,v2,v3>


        double det = inner_prod(edge10, plane[3].mNormal);
        if (det < 0.00)
        {
            // The normals are inner pointing, reverse their directions.
            for (int i = 0; i < 4; ++i)
            {
                plane[i].mNormal = -plane[i].mNormal;
            }
        }

        for (int i = 0; i < 4; ++i)
        {
            plane[i].mConstant = inner_prod(geom_1[i].Coordinates(), plane[i].mNormal);
        }

    }



    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    virtual std::string Info() const
    {
        return "3 dimensional tetrahedra with four nodes in 3D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "3 dimensional tetrahedra with four nodes in 3D space";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points by the order they stored in the geometry
     * and then center point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
        std::cout << std::endl;
        rOStream << "    in Tetrahedra3D4 PrintData\t : " << std::endl;
        Matrix jacobian(3,3);
        this->Jacobian(jacobian, PointType());
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

protected:

    /**
     * there are no protected class members
     */

private:

    /**
     * Static Member Variables
     */
    static const GeometryData msGeometryData;


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
    }

    // serialization needs the default constructor
    Tetrahedra3D4(): BaseType(PointsArrayType(), &msGeometryData) {}

    /**
     * Private Operations
     */

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */

    static Matrix& CalculateShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint )
    {
        rResult(0,0) = -1.0;
        rResult(0,1) = -1.0;
        rResult(0,2) = -1.0;
        rResult(1,0) =  1.0;
        rResult(1,1) =  0.0;
        rResult(1,2) =  0.0;
        rResult(2,0) =  0.0;
        rResult(2,1) =  1.0;
        rResult(2,2) =  0.0;
        rResult(3,0) =  0.0;
        rResult(3,1) =  0.0;
        rResult(3,2) =  1.0;
        return rResult;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(
        typename BaseType::IntegrationMethod ThisMethod)
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 4;
        //setting up return matrix
        Matrix shape_function_values(integration_points_number, points_number );
        //loop over all integration points
        for(int pnt = 0; pnt < integration_points_number; pnt++)
        {
            shape_function_values(pnt,0) = ( 1.0
                                             -integration_points[pnt].X()
                                             -integration_points[pnt].Y()
                                             -integration_points[pnt].Z() );
            shape_function_values(pnt,1) = integration_points[pnt].X();
            shape_function_values(pnt,2) = integration_points[pnt].Y();
            shape_function_values(pnt,3) = integration_points[pnt].Z();
        }
        return shape_function_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the local gradients of all shape functions in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     *
     */
    static ShapeFunctionsGradientsType
    CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod)
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values(integration_points_number);
        //initialising container
        //loop over all integration points
        for( int pnt=0; pnt<integration_points_number; pnt++ )
        {
            Matrix result = ZeroMatrix(4,3);
            result(0,0) = -1.0;
            result(0,1) = -1.0;
            result(0,2) = -1.0;
            result(1,0) =  1.0;
            result(1,1) =  0.0;
            result(1,2) =  0.0;
            result(2,0) =  0.0;
            result(2,1) =  1.0;
            result(2,2) =  0.0;
            result(3,0) =  0.0;
            result(3,1) =  0.0;
            result(3,2) =  1.0;
            d_shape_f_values[pnt] = result;
        }
        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature<TetrahedronGaussLegendreIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TetrahedronGaussLegendreIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TetrahedronGaussLegendreIntegrationPoints3,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TetrahedronGaussLegendreIntegrationPoints4,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TetrahedronGaussLegendreIntegrationPoints5,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_3),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_4),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_5)
            }
        };
        return shape_functions_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    static const ShapeFunctionsLocalGradientsContainerType
    AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_1),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_2),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_3),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_4),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_5)
            }
        };
        return shape_functions_local_gradients;
    }


    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Tetrahedra3D4;


    /**
     * Un accessible methods
     */


};// Class Tetrahedra3D4


/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream, Tetrahedra3D4<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream, const Tetrahedra3D4<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


template<class TPointType> const
GeometryData Tetrahedra3D4<TPointType>::msGeometryData(
    3, 3, 3, GeometryData::GI_GAUSS_1,
    Tetrahedra3D4<TPointType>::AllIntegrationPoints(),
    Tetrahedra3D4<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

}// namespace Kratos.

#endif // KRATOS_TETRAHEDRA_3D_4_H_INCLUDED  defined 
