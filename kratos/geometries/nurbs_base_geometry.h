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

#if !defined(KRATOS_NURBS_BASE_GEOMETRY_H_INCLUDED) //modified by Matthias
#define  KRATOS_NURBS_BASE_GEOMETRY_H_INCLUDED  //modified by Matthias

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"

//Constants needed by OpenNURBS-Library
#define ON_SQRT_EPSILON 1.490116119385000000e-8

/////////////////////////////////////////////////////////////////
//INCLUDES IN THIS SECTION ARE TO AVOID INCLUSION OF OPENNURBS.h"
#include <memory.h>
#if defined(ON_COMPILER_XCODE)
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#if defined(ON_COMPILER_IRIX)
#include <alloca.h>
#endif
/////////////////////////////////////////////////////////////////

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
 * A four node quadrilateral geometry. While the shape functions are only defined in
 * 2D it is possible to define an arbitrary orientation in space. Thus it can be used for
 * defining surfaces on 3D elements.
 */

template<class TPointType> class NurbsPatchGeometry  //modified by Matthias
    : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;


    /**
     * Pointer definition of NurbsSurfaceGeometry3D
     */
    KRATOS_CLASS_POINTER_DEFINITION( NurbsPatchGeometry ); //modified by Matthias

    /**
     * Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /**
     * A Vector of counted pointers to Geometries.
     * Used for returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /**
     * Type used for indexing in geometry class.
     * std::size_t used for indexing
     * point or integration point access methods and also all other
     * methods which need point or integration point index.
     */
    typedef typename BaseType::IndexType IndexType;

    /**
     * This type is used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Array of counted pointers to point.
     * This type used to hold geometry's points.
     */
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    /**
     * Array of coordinates. Can be Nodes, Points or IntegrationPoints
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * This type used for representing an integration point in geometry.
     * This integration point is a point with an additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A Vector of IntegrationPointType which used to hold
     * integration points related to an integration
     * method.
     * IntegrationPoints functions used this type to return
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
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

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
     * A third order tensor to hold shape functions' local second derivatives.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType
    ShapeFunctionsSecondDerivativesType;

    /**
    * A third order tensor to hold shape functions' local third derivatives.
    * ShapefunctionsLocalGradients function return this
    * type as its result.
    */
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType
    ShapeFunctionsThirdDerivativesType;

    /**
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * A standard constructor which will create an empty NURBS-surface
    */

    NurbsPatchGeometry(GeometryData const* msGeometryData):
        BaseType( PointsArrayType(), msGeometryData)
    {}




    NurbsPatchGeometry( const PointerVector< Node<3> > ControlPoints,GeometryData const* msGeometryData)
           : BaseType( PointsArrayType(ControlPoints), msGeometryData )
    {
    }




    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     *
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> NurbsPatchGeometry( NurbsPatchGeometry<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~NurbsPatchGeometry() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const
    {
        return GeometryData::Kratos_Quadrilateral;
    }

    GeometryData::KratosGeometryType GetGeometryType() const
    {
        return GeometryData::Kratos_Triangle3D3;
    }

    ///@}
    ///@name Operators
    ///@{




        /**
     * DefineGeometries: In this function the both knot vectors in Xi and Eta direction
     * will be analysed. For each knot span combination (Xi x Eta) there will be generated
     * one new geometry of the type NurbsPatchGeometry. Each geometry will contain
     * exactly the information which is needed to compute the shape functions and its
     * derivatives in a knot span. Furthermore the upper and lower limits for the Xi- and
     * Eta-values are stored.
     *
     * @param GeometryContainer: a container which holds all created geometries.
     */
    virtual void DefineGeometries(std::vector< Geometry< Node<3> >::Pointer>& GeometryContainer)
    {
        KRATOS_ERROR << "Calling base class DefineGeometries method instead of drived class one. Please check the definition of derived class." << *this << std::endl;
    }


    virtual int GeometryNumber()
    {
        KRATOS_ERROR << "Calling base class GeometryNumber method instead of drived class one. Please check the definition of derived class." << *this << std::endl;
    return 0;
    }

    virtual int FindGeometryId(double Xi, double Eta)
    {

        return 0;

    }

protected:



private:
    ///@name Static Member Variables
    ///@{
//    static const GeometryData msGeometryData;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }



    ///@}
    ///@name Member Variables
    ///@{

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
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class NurbsPatchGeometry;

    ///@}
    ///@name Un accessible methods
    ///@{



    ///@}
}; // Class Geometry


//template<class TPointType> const
//GeometryData NurbsPatchGeometry<TPointType>::msGeometryData(
//    2, 3, 2,
//    GeometryData::GI_GAUSS_1,
//    NurbsPatchGeometry<TPointType>::AllIntegrationPoints(),
//    NurbsPatchGeometry<TPointType>::AllShapeFunctionsValues(),
//    AllShapeFunctionsLocalGradients()
//);

}// namespace Kratos.

#endif // KRATOS_QUADRILATERAL_3D_4_H_INCLUDED  defined 

