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
//                   Andreas Apostolatos
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"


namespace Kratos {

template <int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
class NurbsCurveOnSurfaceGeometry : public NurbsCurveGeometry<2, TCurveContainerPointType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    // typedef Geometry<typename TCurveContainerPointType::value_type> BaseType;
    typedef NurbsCurveGeometry<2, TCurveContainerPointType> BaseType;
    typedef NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TSurfaceContainerPointType> GeometryType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef NurbsSurfaceGeometry<3, TSurfaceContainerPointType> NurbsSurfaceType;
    typedef NurbsCurveGeometry<2, TCurveContainerPointType> NurbsCurveType;

    /** Array of counted pointers to point. This type used to hold
        geometry's points.*/
    // typedef  typename TCurveContainerPointType PointsArrayType;
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    typedef  typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /// Counted pointer of NurbsCurveOnSurfaceGeometry
    KRATOS_CLASS_POINTER_DEFINITION(NurbsCurveOnSurfaceGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Conctructor for curve on surface
    /*NurbsCurveOnSurfaceGeometry(
        typename NurbsSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve)
        : BaseType(rThisPoints, &msGeometryData)
        , mpNurbsSurface(pSurface)
        , mpNurbsCurve(pCurve)
    {
    }*/

    NurbsCurveOnSurfaceGeometry(
        typename NurbsSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve)
        : BaseType(
            pCurve->Points(), 
            pCurve->PolynomialDegree(), 
            pCurve->Knots(), 
            pCurve->Weights()
        )
        , mpNurbsSurface(pSurface)
        , mpNurbsCurve(pCurve)
    {
    }

    explicit NurbsCurveOnSurfaceGeometry(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    NurbsCurveOnSurfaceGeometry(NurbsCurveOnSurfaceGeometry const& rOther)
        : BaseType(rOther)
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mpNurbsCurve(rOther.mpNurbsCurve)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherCurveContainerPointType, class TOtherSurfaceContainerPointType> NurbsCurveOnSurfaceGeometry(
        NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TOtherSurfaceContainerPointType> const& rOther)
        : BaseType(rOther)
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mpNurbsCurve(rOther.mpNurbsCurve)
    {
    }

    /// Destructor
    ~NurbsCurveOnSurfaceGeometry() override = default;

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
    NurbsCurveOnSurfaceGeometry& operator=(const NurbsCurveOnSurfaceGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsSurface = rOther.mpNurbsSurface;
        mpNurbsCurve = rOther.mpNurbsCurve;
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
    /*template<class TOtherContainerPointType>
    NurbsCurveOnSurfaceGeometry& operator=(
        NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TOtherContainerPointType> const & rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsSurface = rOther.mpNurbsSurface;
        mpNurbsCurve = rOther.mpNurbsCurve;
        return *this;
    }*/

    ///@}
    ///@name Operations
    ///@{

    /*typename BaseType::Pointer Create(
        PointsArrayType const& ThisPoints) const override
    {
        return Kratos::make_shared<NurbsCurveOnSurfaceGeometry>();
    }*/

    ///@}
    ///@name Get and Set functions
    ///@{

    /*
    * @brief Polynomial degree of this curve.
    * @return the polynomial degree.
    */
    SizeType PolynomialDegree() const
    {
        return mpNurbsCurve->PolynomialDegree();
    }

    /* 
    * @brief Knot vector is defined to have a multiplicity of p
    *        at the beginning and end (NOT: p + 1).
    * @return knot vector.
    */
    /*const Vector& Knots() const
    {
        return mpNurbsCurve->Knots();
    }*/

    ///*
    //* @brief The number of knots within the knot vector.
    //* @return the size of the knot vector.
    //*/
    //SizeType NumberOfKnots() const
    //{
    //    return mKnots.size();
    //}

    ///*
    //* @brief The number of nonzero control points for one point in the curve
    //*        is given by p+1.
    //* @return the number of nonzero control points.
    //*/
    //SizeType NumberOfNonzeroControlPoints() const
    //{
    //    return PolynomialDegree() + 1;
    //}

    ///*
    //* @brief Checks if shape functions are rational or not.
    //* @return true if NURBS,
    //*         false if B-Splines only (all weights are considered as 1.0)
    //*/
    //bool IsRational() const
    //{
    //    return mWeights.size() != 0;
    //}

    ///*
    //* @brief Provides weights vector. All values are 1.0 for B-Splines,
    //*        for NURBS those can be unequal 1.0.
    //*        Has size 0 if B-Spline.
    //* @return weights vector.
    //*/
    //const Vector& Weights() const
    //{
    //    //return mWeights;
    //}

    ///* 
    //* @brief Provides the natural boundaries of the NURBS/B-Spline curve.
    //* @return domain interval.
    //*/
    //Interval DomainInterval() const
    //{
    //    //return Interval(mKnots[mPolynomialDegree - 1], mKnots[NumberOfKnots() - mPolynomialDegree]);
    //}

    ///*
    //* @brief Provides all knot intervals within one curve.
    //* @return vector of domain intervals.
    //*/
    //std::vector<Interval> KnotSpanIntervals() const
    //{
    //    //const IndexType first_span = mPolynomialDegree - 1;
    //    //const IndexType last_span = NumberOfKnots() - mPolynomialDegree - 1;

    //    //const IndexType number_of_spans = last_span - first_span + 1;

    //    //std::vector<Interval> result(number_of_spans);

    //    //for (IndexType i = 0; i < number_of_spans; i++) {
    //    //    const double t0 = mKnots[first_span + i];
    //    //    const double t1 = mKnots[first_span + i + 1];

    //    //    result[i] = Interval(t0, t1);
    //    //}

    //    //return result;
    //}

    ///@}
    ///@name Operations
    ///@{

    /** 
    * @brief This method maps from dimension space to working space and computes the
    *        number of derivatives at the dimension parameter.
    * From ANurbs library (https://github.com/oberbichler/ANurbs)
    * @param LocalCoordinates The local coordinates in dimension space
    * @param Derivative Number of computed derivatives
    * @return std::vector<array_1d<double, 3>> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    std::vector<CoordinatesArrayType> GlobalDerivatives(
        const CoordinatesArrayType& rCoordinates,
        const SizeType DerivativeOrder) const
    {
        // derivatives of base geometries
        auto curve_derivatives = mpNurbsCurve->GlobalDerivatives(rCoordinates, DerivativeOrder);

        const CoordinatesArrayType surface_coordinates = {curve_derivatives[0][0], curve_derivatives[1][0]};

        auto surface_derivatives = mpNurbsSurface->GlobalDerivatives(surface_coordinates, DerivativeOrder);

        // compute combined derivatives

        std::vector<array_1d<double, 3>> derivatives(DerivativeOrder + 1);

        std::function<Vector(int, int, int)> c;

        c = [&](int DerivativeOrder, int i, int j) -> Vector {
            if (DerivativeOrder > 0) {
                Vector result = ZeroVector();

                for (int a = 1; a <= DerivativeOrder; a++) {
                    result += (
                        c(DerivativeOrder - a, i + 1, j) * curve_derivatives[0][a] +
                        c(DerivativeOrder - a, i, j + 1) * curve_derivatives[1][a]
                        ) * NurbsUtilities::GetBinomCoefficient(DerivativeOrder - 1, a - 1);
                }

                return result;
            }
            else {
                const int index = surface_derivatives->IndexOfShapeFunctionRow(i, j);
                return surface_derivatives[index];
            }
        };

        for (int i = 0; i <= DerivativeOrder; i++) {
            derivatives[i] = c(i, 0, 0);
        }

        return derivatives;
    }
        //NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, std::min(DerivativeOrder, PolynomialDegree()));

        //if (IsRational()) {
        //    shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rCoordinates[0]);
        //}
        //else {
        //    shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rCoordinates[0]);
        //}

        //std::vector<array_1d<double, 3>> derivatives(DerivativeOrder + 1);

        //for (IndexType order = 0; order < shape_function_container.NumberOfShapeFunctionRows(); order++) {
        //    IndexType index_0 = shape_function_container.GetFirstNonzeroControlPoint();
        //    derivatives[order] = (*this)[index_0] * shape_function_container(0, order);
        //    for (IndexType u = 1; u < shape_function_container.NumberOfNonzeroControlPoints(); u++) {
        //        IndexType index = shape_function_container.GetFirstNonzeroControlPoint() + u;

        //        derivatives[order] += (*this)[index] * shape_function_container(u, order);
        //    }
        //}

        ////fill up the vector with zeros
        //for (IndexType order = shape_function_container.NumberOfShapeFunctionRows(); order <= DerivativeOrder; order++) {
        //    IndexType index_0 = shape_function_container.GetFirstNonzeroControlPoint();
        //    derivatives[order] = (*this)[index_0] * 0.0;
        //}
        //return derivatives;
    // }

    /*
    * @brief This method maps from dimension space to working space.
    * From Piegl and Tiller, The NURBS Book, Algorithm A3.1/ A4.1
    * @param rResult array_1d<double, 3> with the coordinates in working space
    * @param LocalCoordinates The local coordinates in dimension space
    * @return array_1d<double, 3> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    /*CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates
    ) const override
    {
        const auto uv = m_curve_geometry->point_at(t);

        const auto point = m_surface_geometry->point_at(uv[0], uv[1]);

        return point;
        //NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, 0);

        //if (IsRational()) {
        //    shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rLocalCoordinates[0]);
        //}
        //else {
        //    shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rLocalCoordinates[0]);
        //}

        //noalias(rResult) = ZeroVector(3);
        //for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
        //    const IndexType index = shape_function_container.GetFirstNonzeroControlPoint() + i;

        //    rResult += (*this)[index] * shape_function_container(i, 0);
        //}
        //return rResult;
    }*/

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
        //NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, 0);

        //if (IsRational()) {
        //    shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rCoordinates[0]);
        //}
        //else {
        //    shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rCoordinates[0]);
        //}

        //if (rResult.size() != shape_function_container.NumberOfNonzeroControlPoints())
        //    rResult.resize(shape_function_container.NumberOfNonzeroControlPoints());

        //for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
        //    rResult[i] = shape_function_container(i, 0);
        //}

        //return rResult;
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
        //NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, 1);

        //if (IsRational()) {
        //    shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rCoordinates[0]);
        //}
        //else {
        //    shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rCoordinates[0]);
        //}

        //if (rResult.size1() != 1
        //    && rResult.size2() != shape_function_container.NumberOfNonzeroControlPoints())
        //    rResult.resize(1, shape_function_container.NumberOfNonzeroControlPoints());

        //for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
        //    rResult(0, i) = shape_function_container(i, 1);
        //}

        //return rResult;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "2 dimensional nurbs curve on 3D surface.";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "2 dimensional nurbs curve on 3D surface.";
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

    ///@}
    ///@name Private Member Variables
    ///@{

    typename NurbsSurfaceType::Pointer mpNurbsSurface;
    typename NurbsCurveType::Pointer mpNurbsCurve;

    ///@}
    ///@name Private Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("pNurbsSurface", mpNurbsSurface);
        rSerializer.save("pNurbsCurve", mpNurbsCurve);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("pNurbsSurface", mpNurbsSurface);
        rSerializer.load("pNurbsCurve", mpNurbsCurve);
    }

    ///@}

}; // class NurbsCurveOnSurfaceGeometry

/*template<int TWorkingSpaceDimension, class TContainerPointType>TPointType
const GeometryData NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimTPointTypeension, TContainerPointType>::msGeometryData(
    1,
    TWorkingSpaceDimension,
    2,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});*/

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED defined