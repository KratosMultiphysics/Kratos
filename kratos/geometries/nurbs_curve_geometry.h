//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//                   Tobias Teschemacher
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED

// System includes
#include <stdexcept>
#include <vector>

#include "includes/ublas_interface.h"

// External includes

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

#include "containers/array_1d.h"



namespace Kratos {

template <int TWorkingSpaceDimension, class TPointType>
class NurbsCurveGeometry : public Geometry<TPointType>
{
public:
    ///@name Type Definitions
    ///@{
    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;
    typedef NurbsCurveGeometry<TWorkingSpaceDimension, TPointType> GeometryType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    /** Array of counted pointers to point. This type used to hold
        geometry's points.*/
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    /// Counted pointer of NurbsSurfaceShapeFunction
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NurbsCurveGeometry);
    ///@}
    ///@name Life Cycle
    ///@{

    /* Conctructor for B-Spline curves. */
    NurbsCurveGeometry(
        const PointsArrayType& rThisPoints,
        const int PolynomialDegree,
        const Vector& rKnots)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mPolynomialDegree(PolynomialDegree)
        , mKnots(rKnots)
        , mIsRational(false)
    {
        KRATOS_DEBUG_ERROR_IF(rKnots.size() != NurbsUtilities::GetNbKnots(PolynomialDegree, rThisPoints.size()))
            << "Number of knots and control points do not match!" << std::endl;
    }

    /* Conctructor for NURBS curves. */
    NurbsCurveGeometry(
        const PointsArrayType& rThisPoints,
        const int PolynomialDegree,
        const Vector& rKnots,
        const Vector& rWeights)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mPolynomialDegree(PolynomialDegree)
        , mKnots(rKnots)
        , mIsRational(true)
        , mWeights(rWeights)
    {
        KRATOS_DEBUG_ERROR_IF(rKnots.size() != NurbsUtilities::GetNumberOfKnots(PolynomialDegree, rThisPoints.size()))
            << "Number of knots and control points do not match!" << std::endl;

        KRATOS_DEBUG_ERROR_IF(rWeights.size() != rThisPoints.size())
            << "Number of control points and weights do not match!" << std::endl;
    }

    explicit NurbsCurveGeometry(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /* Copy constructor.*/
    NurbsCurveGeometry(NurbsCurveGeometry const& rOther)
        : BaseType(rOther)
    {
    }

    /* Copy constructor from a geometry with different point type.*/
    template<class TOtherPointType> NurbsCurveGeometry(
        NurbsCurveGeometry<TWorkingSpaceDimension, TOtherPointType> const& rOther)
        : BaseType(rOther)
    {
    }

    /* Destructor.*/
    ~NurbsCurveGeometry() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_generic_family;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_generic_type;
    }

    ///@}
    ///@name Get and Set functions
    ///@{

    /* Checks if shape functions are rational or not.
    @return true if NURBS, false if B-Splines only (all weights are considered as 1)
    */
    bool IsRational() const
    {
        return mWeights.size() != 0;
    }

    /* Get Knot vector. This vector is defined to have a multiplicity of p
    at the beginning and end (NOT: p + 1).
    @return knot vector.
    */
    const Vector& Knots() const
    {
        return mKnots;
    }

    /* @return Gives the size of the knot vector.
    */
    SizeType NumberOfKnots() const
    {
        return mKnots.size();
    }

    /* Get Weights vector. All values are 1.0 for B-Splines, for NURBS those can be unequal 1.0.
    @return weights vector.
    */
    const Vector& Weights() const
    {
        return mWeights;
    }

    /* Provides the natural boundaries of the NURBS/B-Spline curve.
    @return domain interval.
    */
    Interval DomainInterval() const
    {
        return Interval(mKnots[mPolynomialDegree - 1], mKnots[NumberOfKnots() - mPolynomialDegree]);
    }

    /* Provides all knot span intervals of the curve.
    @return vector of knot span intervals.
    */
    std::vector<Interval> KnotSpanIntervals() const
    {
        const SizeType first_span = mPolynomialDegree - 1;
        const SizeType last_span = NumberOfKnots() - mPolynomialDegree - 1;

        const SizeType number_of_spans = last_span - first_span + 1;

        std::vector<Interval> result(number_of_spans);

        for (int i = 0; i < number_of_spans; i++) {
            const double t0 = mKnots[first_span + i];
            const double t1 = mKnots[first_span + i + 1];

            result[i] = Interval(t0, t1);
        }

        return result;
    }

    ///@}
    ///@name Operations
    ///@{

    /** This method maps from local space to working space and computes the
    * number of derivatives at the local space parameter in the dimension of the object.
    * @param LocalCoordinates The local coordinates in dimension space
    * @param Derivative Number of computed derivatives
    * @return std::vector<array_1d<double, 3>> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    std::vector<CoordinatesArrayType> GlobalDerivatives(
        const CoordinatesArrayType& rCoordinates,
        const int DerivativeOrder) const
    {
        NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, DerivativeOrder);

        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rCoordinates[0]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rCoordinates[0]);
        }

        std::vector<array_1d<double, 3>> derivatives(shape_function_container.NumberOfShapeFunctionRows());

        for (int order = 0; order < shape_function_container.NumberOfShapeFunctionRows(); order++) {
            int index_0 = shape_function_container.GetFirstNonzeroControlPoint();
            derivatives[order] = (*this)[index_0] * shape_function_container(order, 0);

            for (int i = 1; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
                int index = shape_function_container.GetFirstNonzeroControlPoint() + i;

                derivatives[order] += (*this)[index] * shape_function_container(order, i);
            }
        }

        return derivatives;
    }

    /** This method maps from dimension space to working space.
    * @param rResult array_1d<double, 3> with the coordinates in working space
    * @param LocalCoordinates The local coordinates in dimension space
    * @return array_1d<double, 3> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates
    ) const override
    {
        NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, 0);

        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rLocalCoordinates[0]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rLocalCoordinates[0]);
        }

        noalias(rResult) = ZeroVector(3);
        for (int i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            const int index = shape_function_container.GetFirstNonzeroControlPoint() + i;

            rResult += (*this)[index] * shape_function_container(0, i);
        }
        return rResult;
    }

    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, 0);

        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rCoordinates[0]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rCoordinates[0]);
        }

        if (rResult.size() != shape_function_container.NumberOfNonzeroControlPoints())
            rResult.resize(shape_function_container.NumberOfNonzeroControlPoints());

        for (int i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            rResult[i] = shape_function_container(0, i);
        }

        return rResult;
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, 1);

        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rCoordinates[0]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rCoordinates[0]);
        }

        if (rResult.size1() != 1
            && rResult.size2() != shape_function_container.NumberOfNonzeroControlPoints())
            rResult.resize(1, shape_function_container.NumberOfNonzeroControlPoints());

        for (int i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            rResult(0, i) = shape_function_container(1, i);
        }

        return rResult;
    }

    ///@}
    ///@name Information
    ///@{
    std::string Info() const override
    {
        return TWorkingSpaceDimension + " dimensional nurbs curve.";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << TWorkingSpaceDimension << " dimensional nurbs curve.";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }
    ///@}

protected:

private:
    ///@name Private Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    ///@}
    ///@name Private Member Variables
    ///@{

    const int mPolynomialDegree;
    Vector mKnots;
    bool mIsRational;
    Vector mWeights;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("PolynomialDegree", mPolynomialDegree);
        rSerializer.save("Knots", mKnots);
        rSerializer.save("IsRational", mIsRational);
        rSerializer.save("Weights", mWeights);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("PolynomialDegree", mPolynomialDegree);
        rSerializer.load("Knots", mKnots);
        rSerializer.load("IsRational", mIsRational);
        rSerializer.load("Weights", mWeights);
    }

    NurbsCurveGeometry() : BaseType(PointsArrayType(), &msGeometryData) {};

    ///@}
    ///@name Private Friends
    ///@{

    template<int TWorkingSpaceDimension, class TOtherPointType> friend class NurbsCurveGeometry;

    ///@}
}; // class NurbsCurveGeometry

template<int TWorkingSpaceDimension, class TPointType>
const GeometryData NurbsCurveGeometry<TWorkingSpaceDimension, TPointType>::msGeometryData(
    2,
    TWorkingSpaceDimension,
    2,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED defined