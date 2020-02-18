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

#if !defined(KRATOS_NURBS_SURFACE_GEOMETRY_H_INCLUDED )
#define  KRATOS_NURBS_SURFACE_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"



namespace Kratos {

template <int TWorkingSpaceDimension, class TContainerPointType>
class NurbsSurfaceGeometry : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<typename TContainerPointType::value_type> BaseType;
    typedef NurbsSurfaceGeometry<TWorkingSpaceDimension, TContainerPointType> GeometryType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    /** Array of counted pointers to point. This type used to hold
        geometry's points.*/
    typedef  typename BaseType::PointsArrayType PointsArrayType;
    typedef  typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /// Counted pointer of NurbsSurfaceGeometry
    KRATOS_CLASS_POINTER_DEFINITION(NurbsSurfaceGeometry);
    ///@}
    ///@name Life Cycle
    ///@{

    /// Conctructor for B-Spline surfaces
    NurbsSurfaceGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const Vector& rKnotsU,
        const Vector& rKnotsV)
        : BaseType(rThisPoints, &msGeometryData)
        , mPolynomialDegreeU(PolynomialDegreeU)
        , mPolynomialDegreeV(PolynomialDegreeV)
        , mKnotsU(rKnotsU)
        , mKnotsV(rKnotsV)
    {
        CheckAndFitKnotVectors();
    }

    /// Conctructor for NURBS surfaces
    NurbsSurfaceGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rWeights)
        : BaseType(rThisPoints, &msGeometryData)
        , mPolynomialDegreeU(PolynomialDegreeU)
        , mPolynomialDegreeV(PolynomialDegreeV)
        , mKnotsU(rKnotsU)
        , mKnotsV(rKnotsV)
        , mWeights(rWeights)
    {
        CheckAndFitKnotVectors();

        KRATOS_ERROR_IF(rWeights.size() != rThisPoints.size())
            << "Number of control points and weights do not match!" << std::endl;
    }

    explicit NurbsSurfaceGeometry(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    NurbsSurfaceGeometry(NurbsSurfaceGeometry<TWorkingSpaceDimension, TContainerPointType> const& rOther)
        : BaseType(rOther)
        , mPolynomialDegreeU(rOther.mPolynomialDegreeU)
        , mPolynomialDegreeV(rOther.mPolynomialDegreeV)
        , mKnotsU(rOther.mKnotsU)
        , mKnotsV(rOther.mKnotsV)
        , mWeights(rOther.mWeights)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherContainerPointType> NurbsSurfaceGeometry(
        NurbsSurfaceGeometry<TWorkingSpaceDimension, TOtherContainerPointType> const& rOther)
        : BaseType(rOther, &msGeometryData)
        , mPolynomialDegreeU(rOther.mPolynomialDegreeU)
        , mPolynomialDegreeV(rOther.mPolynomialDegreeV)
        , mKnotsU(rOther.mKnotsU)
        , mKnotsV(rOther.mKnotsV)
        , mWeights(rOther.mWeights)
    {
    }

    /// Destructor.
    ~NurbsSurfaceGeometry() override = default;

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
    NurbsSurfaceGeometry& operator=(const NurbsSurfaceGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mPolynomialDegreeU = rOther.mPolynomialDegreeU;
        mPolynomialDegreeV = rOther.mPolynomialDegreeV;
        mKnotsU = rOther.mKnotsU;
        mKnotsV = rOther.mKnotsV;
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
    NurbsSurfaceGeometry& operator=(
        NurbsSurfaceGeometry<TWorkingSpaceDimension, TOtherContainerPointType> const & rOther)
    {
        BaseType::operator=(rOther);
        mPolynomialDegreeU = rOther.mPolynomialDegreeU;
        mPolynomialDegreeV = rOther.mPolynomialDegreeV;
        mKnotsU = rOther.mKnotsU;
        mKnotsV = rOther.mKnotsV;
        mWeights = rOther.mWeights;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        PointsArrayType const& ThisPoints) const override
    {
        return Kratos::make_shared<NurbsSurfaceGeometry>(ThisPoints);
    }

    ///@}
    ///@name Get and Set functions
    ///@{

    /* @return returns the polynomial degree 'p' in u direction. */
    SizeType PolynomialDegreeU() const
    {
        return mPolynomialDegreeU;
    }

    /* @return returns the polynomial degree 'p' in u direction. */
    SizeType PolynomialDegreeV() const
    {
        return mPolynomialDegreeV;
    }

    /* Get Knot vector in u-direction. This vector is defined to have
    a multiplicity of p at the beginning and end (NOT: p + 1).
    @return knot vector. */
    const Vector& KnotsU() const
    {
        return mKnotsU;
    }

    /* Get Knot vector in v-direction. This vector is defined to have
    a multiplicity of p at the beginning and end (NOT: p + 1).
    @return knot vector. */
    const Vector& KnotsV() const
    {
        return mKnotsV;
    }

    /* @return Gives the size of the knot vector in u-direction.
    */
    SizeType NumberOfKnotsU() const
    {
        return mKnotsU.size();
    }

    /* @return Gives the size of the knot vector in v-direction.
    */
    SizeType NumberOfKnotsV() const
    {
        return mKnotsV.size();
    }

    /* Checks if shape functions are rational or not.
    @return true if NURBS, false if B-Splines only (all weights are considered as 1) */
    bool IsRational() const
    {
        return mWeights.size() != 0;
    }

    /* Get Weights vector. All values are 1.0 for B-Splines, for NURBS those can be unequal 1.0.
    @return weights vector.
    */
    const Vector& Weights() const
    {
        return mWeights;
    }

    SizeType NumberOfControlPointsU() const
    {
        return NumberOfKnotsU() - PolynomialDegreeU() + 1;
    }

    SizeType NumberOfControlPointsV() const
    {
        return NumberOfKnotsV() - PolynomialDegreeV() + 1;
    }

    /* Provides the natural boundaries of the NURBS/B-Spline surface.
    @return domain interval.
    */
    NurbsInterval DomainIntervalU() const
    {
        return NurbsInterval(
            mKnotsU[mPolynomialDegreeU - 1],
            mKnotsU[NumberOfKnotsU() - mPolynomialDegreeU]);
    }

    /* Provides the natural boundaries of the NURBS/B-Spline surface.
    @return domain interval.
    */
    NurbsInterval DomainIntervalV() const
    {
        return NurbsInterval(
            mKnotsV[mPolynomialDegreeV - 1],
            mKnotsV[NumberOfKnotsV() - mPolynomialDegreeV]);
    }

    /* Provides all knot span intervals of the surface in u-direction.
    @return vector of knot span intervals.
    */
    std::vector<NurbsInterval> KnotSpanIntervalsU() const
    {
        const SizeType first_span = mPolynomialDegreeU - 1;
        const SizeType last_span = NumberOfKnotsU() - mPolynomialDegreeU - 1;

        const SizeType number_of_spans = last_span - first_span + 1;

        std::vector<NurbsInterval> result(number_of_spans);

        for (int i = 0; i < number_of_spans; i++) {
            const double t0 = mKnotsU[first_span + i];
            const double t1 = mKnotsU[first_span + i + 1];

            result[i] = NurbsInterval(t0, t1);
        }

        return result;
    }

    /* Provides all knot span intervals of the surface in u-direction.
    @return vector of knot span intervals.
    */
    std::vector<NurbsInterval> KnotSpanIntervalsV() const
    {
        const SizeType first_span = mPolynomialDegreeV - 1;
        const SizeType last_span = NumberOfKnotsV() - mPolynomialDegreeV - 1;

        const SizeType number_of_spans = last_span - first_span + 1;

        std::vector<NurbsInterval> result(number_of_spans);

        for (int i = 0; i < number_of_spans; i++) {
            const double t0 = mKnotsV[first_span + i];
            const double t1 = mKnotsV[first_span + i + 1];

            result[i] = NurbsInterval(t0, t1);
        }

        return result;
    }

    ///@}
    ///@name Operations
    ///@{

    /** This method maps from dimension space to working space.
    * @param rResult array_1d<double, 3> with the coordinates in working space
    * @param LocalCoordinates The local coordinates in dimension space,
                              here only the first 2 entries are taken into account.
    * @return array_1d<double, 3> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates
    ) const override
    {
        NurbsSurfaceShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, 0);

        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(
                mKnotsU, mKnotsV, mWeights, rLocalCoordinates[0], rLocalCoordinates[1]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(
                mKnotsU, mKnotsV, rLocalCoordinates[0], rLocalCoordinates[1]);
        }

        noalias(rResult) = ZeroVector(3);
        for (IndexType u = 0; u <= PolynomialDegreeU(); u++) {
            for (IndexType v = 0; v <= PolynomialDegreeV(); v++) {
                IndexType cp_index_u = shape_function_container.GetFirstNonzeroControlPointU() + u;
                IndexType cp_index_v = shape_function_container.GetFirstNonzeroControlPointV() + v;

                const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                    NumberOfControlPointsU(), NumberOfControlPointsV(), cp_index_u, cp_index_v);

                rResult += (*this)[index] * shape_function_container(u, v, 0);
            }
        }

        return rResult;
    }

    /** This method maps from local space to working space and computes the
    *   number of derivatives at the local space parameter in the dimension of the object.
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
        NurbsSurfaceShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, DerivativeOrder);

        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(
                mKnotsU, mKnotsV, mWeights, rLocalCoordinates[0], rLocalCoordinates[1]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(
                mKnotsU, mKnotsV, rLocalCoordinates[0], rLocalCoordinates[1]);
        }

        if (rGlobalSpaceDerivatives.size() != shape_function_container.NumberOfShapeFunctionRows()) {
            rGlobalSpaceDerivatives.resize(shape_function_container.NumberOfShapeFunctionRows());
        }

        for (IndexType shape_function_row_i = 0;
            shape_function_row_i < shape_function_container.NumberOfShapeFunctionRows();
            shape_function_row_i++) {
            for (IndexType u = 0; u <= PolynomialDegreeU(); u++) {
                for (IndexType v = 0; v <= PolynomialDegreeV(); v++) {
                    IndexType cp_index_u = shape_function_container.GetFirstNonzeroControlPointU() + u;
                    IndexType cp_index_v = shape_function_container.GetFirstNonzeroControlPointV() + v;

                    const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        NumberOfControlPointsU(), NumberOfControlPointsV(), cp_index_u, cp_index_v);

                    if (u == 0 && v==0)
                        rGlobalSpaceDerivatives[shape_function_row_i] =
                        (*this)[index] * shape_function_container(u, v, shape_function_row_i);
                    else
                        rGlobalSpaceDerivatives[shape_function_row_i] +=
                        (*this)[index] * shape_function_container(u, v, shape_function_row_i);
                }
            }
        }
    }

    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        NurbsSurfaceShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, 0);

        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(mKnotsU, mKnotsV, mWeights, rCoordinates[0], rCoordinates[1]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(mKnotsU, mKnotsV, rCoordinates[0], rCoordinates[1]);
        }

        if (rResult.size() != shape_function_container.NumberOfNonzeroControlPoints())
            rResult.resize(shape_function_container.NumberOfNonzeroControlPoints());

        for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            rResult[i] = shape_function_container(i, 0);
        }

        return rResult;
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        NurbsSurfaceShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, 0);

        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(mKnotsU, mKnotsV, mWeights, rCoordinates[0], rCoordinates[1]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(mKnotsU, mKnotsV, rCoordinates[0], rCoordinates[1]);
        }

        if (rResult.size1() != 2
            && rResult.size2() != shape_function_container.NumberOfNonzeroControlPoints())
            rResult.resize(2, shape_function_container.NumberOfNonzeroControlPoints());

        for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            rResult(0, i) = shape_function_container(i, 1);
            rResult(1, i) = shape_function_container(i, 2);
        }

        return rResult;
    }

    ///@}
    ///@name Information
    ///@{
    std::string Info() const override
    {
        return std::to_string(TWorkingSpaceDimension) + " dimensional nurbs surface.";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << TWorkingSpaceDimension << " dimensional nurbs surface.";
    }

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

    SizeType mPolynomialDegreeU;
    SizeType mPolynomialDegreeV;
    Vector mKnotsU;
    Vector mKnotsV;
    Vector mWeights;

    ///@}
    ///@name Private Operations
    ///@{

    /*
    * @brief Checks if the knot vector is coinciding with the number of
    *        control points and the polynomial degree. If the knot vectors
    *        have a multiplicity of p+1 in the beginning, it is reduced to p.
    */
    void CheckAndFitKnotVectors()
    {
        SizeType num_control_points = this->size();

        if (num_control_points !=
            (NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeU, mKnotsU.size())
                * NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeV, mKnotsV.size()))) {
            if (num_control_points ==
                (NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeU, mKnotsU.size() - 2)
                    * NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeV, mKnotsV.size() - 2))) {
                Vector KnotsU = ZeroVector(mKnotsU.size() - 2);
                for (SizeType i = 0; i < mKnotsU.size() - 2; ++i) {
                    KnotsU[i] = mKnotsU[i + 1];
                }
                mKnotsU = KnotsU;

                Vector KnotsV = ZeroVector(mKnotsV.size() - 2);
                for (SizeType i = 0; i < mKnotsV.size() - 2; ++i) {
                    KnotsV[i] = mKnotsV[i + 1];
                }
                mKnotsV = KnotsV;
            } else {
                KRATOS_ERROR
                    << "Number of controls points and polynomial degrees and number of knots do not match! " << std::endl
                    << " P: " << mPolynomialDegreeU << ", Q: " << mPolynomialDegreeV
                    << ", number of knots u: " << mKnotsU.size() << ", number of knots v: " << mKnotsV.size()
                    << ", number of control points: " << num_control_points << std::endl
                    << "Following condition must be achieved: ControlPoints.size() = (KnotsU.size() - P + 1) * (KnotsV.size() - Q + 1)"
                    << std::endl;
            }
        }
    }

    ///@}
    ///@name Private Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("PolynomialDegreeU", mPolynomialDegreeU);
        rSerializer.save("PolynomialDegreeV", mPolynomialDegreeV);
        rSerializer.save("KnotsU", mKnotsU);
        rSerializer.save("KnotsV", mKnotsV);
        rSerializer.save("Weights", mWeights);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("PolynomialDegreeU", mPolynomialDegreeU);
        rSerializer.load("PolynomialDegreeV", mPolynomialDegreeV);
        rSerializer.load("KnotsU", mKnotsU);
        rSerializer.load("KnotsV", mKnotsV);
        rSerializer.load("Weights", mWeights);
    }

    NurbsSurfaceGeometry() : BaseType(PointsArrayType(), &msGeometryData) {};

    ///@}

}; // class NurbsSurfaceGeometry

template<int TWorkingSpaceDimension, class TPointType>
const GeometryData NurbsSurfaceGeometry<TWorkingSpaceDimension, TPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TPointType>
const GeometryDimension NurbsSurfaceGeometry<TWorkingSpaceDimension, TPointType>::msGeometryDimension(
    2, TWorkingSpaceDimension, 2);

} // namespace Kratos

#endif // KRATOS_NURBS_SURFACE_GEOMETRY_H_INCLUDED defined
