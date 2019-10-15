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
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"


namespace Kratos {

template <int TWorkingSpaceDimension, class TContainerPointType>
class NurbsCurveGeometry : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<typename TContainerPointType::value_type> BaseType;
    typedef NurbsCurveGeometry<TWorkingSpaceDimension, TContainerPointType> GeometryType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    /** Array of counted pointers to point. This type used to hold
        geometry's points.*/
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    typedef  typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /// Counted pointer of NurbsCurveGeometry
    KRATOS_CLASS_POINTER_DEFINITION(NurbsCurveGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Conctructor for B-Spline curves
    NurbsCurveGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegree,
        const Vector& rKnots)
        : BaseType(rThisPoints, &msGeometryData)
        , mPolynomialDegree(PolynomialDegree)
        , mKnots(rKnots)
    {
        KRATOS_ERROR_IF(rKnots.size() != NurbsUtilities::GetNumberOfKnots(PolynomialDegree, rThisPoints.size()))
            << "Number of knots and control points do not match!" << std::endl;
    }

    /// Conctructor for NURBS curves
    NurbsCurveGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegree,
        const Vector& rKnots,
        const Vector& rWeights)
        : BaseType(rThisPoints, &msGeometryData)
        , mPolynomialDegree(PolynomialDegree)
        , mKnots(rKnots)
        , mWeights(rWeights)
    {
        KRATOS_ERROR_IF(rKnots.size() != NurbsUtilities::GetNumberOfKnots(PolynomialDegree, rThisPoints.size()))
            << "Number of knots and control points do not match!" << std::endl;

        KRATOS_ERROR_IF(rWeights.size() != rThisPoints.size())
            << "Number of control points and weights do not match!" << std::endl;
    }

    explicit NurbsCurveGeometry(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    NurbsCurveGeometry(NurbsCurveGeometry<TWorkingSpaceDimension, TContainerPointType>  const& rOther)
        : BaseType(rOther)
        , mPolynomialDegree(rOther.mPolynomialDegree)
        , mKnots(rOther.mKnots)
        , mWeights(rOther.mWeights)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherContainerPointType> NurbsCurveGeometry(
        NurbsCurveGeometry<TWorkingSpaceDimension, TOtherContainerPointType> const& rOther)
        : BaseType(rOther, &msGeometryData)
        , mPolynomialDegree(rOther.mPolynomialDegree)
        , mKnots(rOther.mKnots)
        , mWeights(rOther.mWeights)
    {
    }

    /// Destructor
    ~NurbsCurveGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

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
    NurbsCurveGeometry& operator=(const NurbsCurveGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mPolynomialDegree = rOther.mPolynomialDegree;
        mKnots = rOther.mKnots;
        mWeights = rOther.mWeights;
        return *this;
    }

    /**
     * @brief Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherContainerPointType>
    NurbsCurveGeometry& operator=(
        NurbsCurveGeometry<TWorkingSpaceDimension, TOtherContainerPointType> const & rOther)
    {
        BaseType::operator=(rOther);
        mPolynomialDegree = rOther.mPolynomialDegree;
        mKnots = rOther.mKnots;
        mWeights = rOther.mWeights;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        PointsArrayType const& ThisPoints) const override
    {
        return Kratos::make_shared<NurbsCurveGeometry>(ThisPoints);
    }

    ///@}
    ///@name Get and Set functions
    ///@{

    /*
    * @brief Polynomial degree of this curve.
    * @return the polynomial degree.
    */
    SizeType PolynomialDegree() const
    {
        return mPolynomialDegree;
    }

    /* 
    * @brief Knot vector is defined to have a multiplicity of p
    *        at the beginning and end (NOT: p + 1).
    * @return knot vector.
    */
    const Vector& Knots() const
    {
        return mKnots;
    }

    /*
    * @brief The number of knots within the knot vector.
    * @return the size of the knot vector.
    */
    SizeType NumberOfKnots() const
    {
        return mKnots.size();
    }

    /*
    * @brief The number of nonzero control points for one point in the curve
    *        is given by p+1.
    * @return the number of nonzero control points.
    */
    SizeType NumberOfNonzeroControlPoints() const
    {
        return PolynomialDegree() + 1;
    }

    /*
    * @brief Checks if shape functions are rational or not.
    * @return true if NURBS,
    *         false if B-Splines only (all weights are considered as 1.0)
    */
    bool IsRational() const
    {
        return mWeights.size() != 0;
    }

    /*
    * @brief Provides weights vector. All values are 1.0 for B-Splines,
    *        for NURBS those can be unequal 1.0.
    *        Has size 0 if B-Spline.
    * @return weights vector.
    */
    const Vector& Weights() const
    {
        return mWeights;
    }

    /* 
    * @brief Provides the natural boundaries of the NURBS/B-Spline curve.
    * @return domain interval.
    */
    NurbsInterval DomainInterval() const
    {
        return NurbsInterval(mKnots[mPolynomialDegree - 1], mKnots[NumberOfKnots() - mPolynomialDegree]);
    }

    /*
    * @brief Provides all knot intervals within one curve.
    * @return vector of domain intervals.
    */
    std::vector<NurbsInterval> KnotSpanIntervals() const
    {
        const IndexType first_span = mPolynomialDegree - 1;
        const IndexType last_span = NumberOfKnots() - mPolynomialDegree - 1;

        const IndexType number_of_spans = last_span - first_span + 1;

        std::vector<NurbsInterval> result(number_of_spans);

        for (IndexType i = 0; i < number_of_spans; i++) {
            const double t0 = mKnots[first_span + i];
            const double t1 = mKnots[first_span + i + 1];

            result[i] = NurbsInterval(t0, t1);
        }

        return result;
    }

    ///@}
    ///@name Operation within Global Space
    ///@{

    /*
    * @brief This method maps from dimension space to working space.
    * From Piegl and Tiller, The NURBS Book, Algorithm A3.1/ A4.1
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
        for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            const IndexType index = shape_function_container.GetFirstNonzeroControlPoint() + i;

            rResult += (*this)[index] * shape_function_container(i, 0);
        }
        return rResult;
    }

    /** 
    * @brief This method maps from dimension space to working space and computes the
    *        number of derivatives at the dimension parameter.
    * From Piegl and Tiller, The NURBS Book, Algorithm A3.2/ A4.2
    * @param LocalCoordinates The local coordinates in dimension space
    * @param Derivative Number of computed derivatives
    * @return std::vector<array_1d<double, 3>> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const override
    {
        NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, DerivativeOrder);

        if (this->IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rLocalCoordinates[0]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rLocalCoordinates[0]);
        }

        if (rGlobalSpaceDerivatives.size() != DerivativeOrder + 1) {
            rGlobalSpaceDerivatives.resize(DerivativeOrder + 1);
        }


        for (IndexType order = 0; order < shape_function_container.NumberOfShapeFunctionRows(); order++) {
            IndexType index_0 = shape_function_container.GetFirstNonzeroControlPoint();
            rGlobalSpaceDerivatives[order] = (*this)[index_0] * shape_function_container(0, order);
            for (IndexType u = 1; u < shape_function_container.NumberOfNonzeroControlPoints(); u++) {
                IndexType index = shape_function_container.GetFirstNonzeroControlPoint() + u;

                rGlobalSpaceDerivatives[order] += (*this)[index] * shape_function_container(u, order);
            }
        }
    }

    ///@}
    ///@name Shape Function
    ///@{

    /*
    * @brief This function computes the shape functions at a certain point.
    * From Piegl and Tiller, The NURBS Book, Algorithm A3.1/ A4.1
    * @param rResult the given Vector, which will be overwritten by the solution
    * @param rCoordinates the given local coordinates, with the coordinates u and v.
    * @return vector of the shape functions at rCoordinates
    */
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

        for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            rResult[i] = shape_function_container(i, 0);
        }

        return rResult;
    }

    /*
    * @brief This function computes the first derivatives at a certain point.
    * From Piegl and Tiller, The NURBS Book, Algorithm A3.2/ A4.2
    * @param rResult the given Matrix which will be overwritten by the solution
    * @param rCoordinates the given local coordinates, with the coordinates u and v.
    * @return matrix of derivatives at rCoordinates.
    *         (0,i): dN/du, (1,i): dN/dv
    */
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

        for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            rResult(0, i) = shape_function_container(i, 1);
        }

        return rResult;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return TWorkingSpaceDimension + " dimensional nurbs curve.";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << TWorkingSpaceDimension << " dimensional nurbs curve.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Private Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Private Member Variables
    ///@{

    SizeType mPolynomialDegree;
    Vector mKnots;
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
        rSerializer.save("Weights", mWeights);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("PolynomialDegree", mPolynomialDegree);
        rSerializer.load("Knots", mKnots);
        rSerializer.load("Weights", mWeights);
    }

    NurbsCurveGeometry() : BaseType(PointsArrayType(), &msGeometryData) {};

    ///@}

}; // class NurbsCurveGeometry

template<int TWorkingSpaceDimension, class TContainerPointType>
const GeometryData NurbsCurveGeometry<TWorkingSpaceDimension, TContainerPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TContainerPointType>
const GeometryDimension NurbsCurveGeometry<TWorkingSpaceDimension, TContainerPointType>::msGeometryDimension(
    1, TWorkingSpaceDimension, 1);

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED defined
