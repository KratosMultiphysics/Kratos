//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#if !defined (KRATOS_PYRAMID_3D_4_H_INCLUDED)
#define KRATOS_PYRAMID_3D_4_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"


namespace Kratos {

///@name Kratos Classes
///@{

/**
 */
template<class TPointType>
class Pyramid3D4 : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;

    /// Pointer definition of Pyramid3D4
    KRATOS_CLASS_POINTER_DEFINITION(Pyramid3D4);

    /** Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /** A Vector of counted pointers to Geometries. Used for
     returning edges of the geometry.
    */
    typedef BaseType::GeometryArrayType GeometryArrayType;

    /** Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /** Type used for indexing in geometry class.std::size_t used for indexing
     point or integration point access methods and also all other
    methods which need point or integration point index.
    */
    typedef BaseType::IndexType IndexType;


    /** This typed used to return size or dimension in
     geometry. Dimension, WorkingDimension, PointsNumber and
    ... return this type as their results.
    */
    typedef BaseType::SizeType SizeType;

    /** Array of counted pointers to point. This type used to hold
     geometry's points.
    */
    typedef  BaseType::PointsArrayType PointsArrayType;

    /** This type used for representing an integration point in
     geometry. This integration point is a point with an
    additional weight component.
    */
    typedef BaseType::IntegrationPointType IntegrationPointType;

    /** A Vector of IntegrationPointType which used to hold
     integration points related to an integration
    method. IntegrationPoints functions used this type to return
    their results.
    */
    typedef BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /** A Vector of IntegrationPointsArrayType which used to hold
     integration points related to different integration method
    implemented in geometry.
    */
    typedef BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
     continer.
    */
    typedef BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
     gradients container in geometry.
    */
    typedef BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold jacobian matrices evaluated at
     integration points. Jacobian and InverseOfJacobian functions
    return this type as their result.
    */
    typedef BaseType::JacobiansType JacobiansType;

    /** A third order tensor to hold shape functions' local
     gradients. ShapefunctionsLocalGradients function return this
    type as its result.
    */
    typedef BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /** Type of the normal vector used for normal to edges in geomety.
     */
    typedef BaseType::NormalType NormalType;

    ///@}
    ///@name Life Cycle
    ///@{

    Pyramid3D4(const PointsArrayType& ThisPoints)
    : BaseType(ThisPoints, msGeometryData)
    {
    }

    /** Copy constructor.
     Construct this geometry as a copy of given geometry.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    Pyramid3D4(Pyramid3D4 const& rOther)
    : BaseType(rOther)
    {
    }

    /** Copy constructor from a geometry with other point type.
     Construct this geometry as a copy of given geometry which
    has different type of points. The given goemetry's
    TOtherPointType* must be implicity convertible to this
    geometry PointType.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    template<class TOtherPointType> Pyramid3D4(Pyramid3D4<TOtherPointType> const& rOther)
    : BaseType(rOther)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /** Assignment operator.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    Pyramid3D4& operator=(const Pyramid3D4& rOther)
    {
        BaseType::operator=(rOther);

        return *this;
    }

    /** Assignment operator for geometries with different point type.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    template<class TOtherPointType>
    Pyramid3D4& operator=(Pyramid3D4<TOtherPointType> const & rOther)
    {
        BaseType::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual Kratos::shared_ptr< Geometry< Point > > Clone() const
    {
        Geometry< Point >::PointsArrayType NewPoints;

        //making a copy of the nodes TO POINTS (not Nodes!!!)
        for(IndexType i = 0 ; i < mPoints.size() ; i++)
        NewPoints.push_back(mPoints[i]);

        //creating a geometry with the new points
        Kratos::shared_ptr< Geometry< Point > > p_clone(new Pyramid3D4< Point >(NewPoints));
        p_clone->ClonePoints();

        return p_clone;
    }

    ///@}
    ///@name Informations
    ///@{

    /** This method calculate and return volume of this
     geometry. For one and two dimensional geometry it returns
    zero and for three dimensional it gives volume of geometry.

    @return double value contains volume.
    @see Length()
    @see Area()
    @see DomainSize()
    */
    virtual double Volume()
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }


    /** Calculates center of this geometry by a simple averaging algorithm.
     Each center point component calculated using:
    \f[
    c_i = \sum_j^n(x_i^j) / n
    \f]

    where \f$ c_i \f$ is component i of center point and \f$
    X_i^j \f$ is component i of j'th point of geometry and n is
    number of the points in this geometry.

    @return PointType which is the calculated center of this geometry.
    */
    virtual PointType Center() const
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    ///@}
    ///@name Shape Function
    ///@{

    /** This method gives value of given shape function evaluated in given
     point.

    @param rPoint Point of evaluation of the shape
    function. This point must be in local coordinate.

    @param ShapeFunctionIndex index of node which correspounding
    shape function evaluated in given integration point.

    @return Value of given shape function in given point.

    @see ShapeFunctionsValues
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    virtual double ShapeFunctionValue(IndexType ShapeFunctionIndex, const PointType& rPoint)
    {
        KRATOS_THROW_ERROR(std::logic_error,
                "Calling base class DeterminantOfJacobian method instead of drived class one. Please check the definition of derived class." , *this);

        return 0;
    }



    virtual void CalculateShapeFunctionsIntegrationPointsGradients(ShapeFunctionsGradientsType& rResult)
    {
        ShapeFunctionsIntegrationPointsGradients(mGeometryData.DefaultIntegrationMethod(), rResult);
    }

    virtual void ShapeFunctionsIntegrationPointsGradients(IntegrationMethod ThisMethod, ShapeFunctionsGradientsType& rResult)
    {
        KRATOS_THROW_ERROR(std::logic_error,
                "Calling base class ShapeFunctionsGaussPointsGradients method instead of derived class one. Please check the definition of derived class." , *this);
    }

    ///@}
    ///@name Input and output
    ///@{

    /** Turn back information as a string.

    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    virtual std::string Info() const
    {
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    /** Print geometry's data into given stream. Prints it's points
     by the order they stored in the geometry and then center
    point of geometry.

    @param rOStream Stream to print into it.
    @see PrintInfo()
    @see Info()
    */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static ShapeFunctionsValuesContainerType CalculateShapeFunctionsIntegrationPointsValues()
    {
    }

    static ShapeFunctionsLocalGradientsContainerType CalculateShapeFunctionsIntegrationPointsLocalGradients()
    {
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class Pyramid3D4;

    ///@}
    ///@name Un accessible methods
    ///@{

    Pyramid3D4();



    ///@}

}; // Class Geometry

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<class TPointType>
inline std::istream& operator >> (std::istream& rIStream,
                    Pyramid3D4<TPointType>& rThis);

/// output stream function
template<class TPointType>
inline std::ostream& operator << (std::ostream& rOStream,
                    const Pyramid3D4<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

template<class TPointType>
typename Triangle2D<TPointType>::IntegrationPointsContainerType Triangle2D<TPointType>::msIntegrationPoints = {
    Quadrature<TriangleGaussLegendreIntegrationPoints<1>, 2, IntegrationPoint<3> >::IntegrationPoints(),
    Quadrature<TriangleGaussLegendreIntegrationPoints<2>, 2, IntegrationPoint<3> >::IntegrationPoints()
};


template<class TPointType>
typename Triangle2D<TPointType>::ShapeFunctionsValuesContainerType
Triangle2D<TPointType>::msShapeFunctionsValues = {
    Triangle2D<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryType::GI_GAUSS_1),
    Triangle2D<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryType::GI_GAUSS_2)
};


template<class TPointType>
typename Triangle2D<TPointType>::ShapeFunctionsLocalGradientsContainerType
Triangle2D<TPointType>::msShapeFunctionsLocalGradients = {
    Triangle2D<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryType::GI_GAUSS_1),
    Triangle2D<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryType::GI_GAUSS_2)
};
template<class TPointType>
typename Pyramid3D4<TPointType>::IntegrationPointsContainerType Pyramid3D4<TPointType>::msIntegrationPoints = {}


template<class TPointType>
typename Pyramid3D4<TPointType>::ShapeFunctionsValuesContainerType
Pyramid3D4<TPointType>::msShapeFunctionsValues = {
//	  Pyramid3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryType::GI_GAUSS_1),
//	  Pyramid3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryType::GI_GAUSS_2)
};

template<class TPointType>
typename Pyramid3D4<TPointType>::ShapeFunctionsLocalGradientsContainerType
Pyramid3D4<TPointType>::msShapeFunctionsLocalGradients( = {
//	  Pyramid3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryType::GI_GAUSS_1),
//	  Pyramid3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryType::GI_GAUSS_2)
};

template<class TPointType>
const GeometryData Pyramid3D4<TPointType>::msGeometryData(/*Dimension*/,
                            /*WorkingSpaceDimension*/,
                            /*LocalSpaceDimension*/,
                            msIntegrationPoints,
                            msShapeFunctionsValues,
                            msShapeFunctionsLocalGradients);

}  // namespace Kratos.

#endif // KRATOS_PYRAMID_3D_4_H_INCLUDED defined
