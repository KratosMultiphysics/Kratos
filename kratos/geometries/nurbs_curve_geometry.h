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

#include "utilities/quadrature_points_utility.h"

#include "integration/integration_point_utilities.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TContainerPointType>
class NurbsCurveGeometry : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename TContainerPointType::value_type NodeType;

    typedef Geometry<NodeType> BaseType;
    typedef NurbsCurveGeometry<TWorkingSpaceDimension, TContainerPointType> GeometryType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    // using base class functionalities.
    using BaseType::CreateQuadraturePointGeometries;
    using BaseType::pGetPoint;
    using BaseType::GetPoint;

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
        CheckAndFitKnotVectors();
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
        CheckAndFitKnotVectors();

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
    ///@name Geometrical Information
    ///@{

    /// Returns number of points per direction.
    SizeType PointsNumberInDirection(IndexType LocalDirectionIndex) const override
    {
        if (LocalDirectionIndex == 0) {
            return this->PointsNumber();
        }
        KRATOS_ERROR << "Possible direction index in NurbsCurveGeometry is 0. Given direction index: "
            << LocalDirectionIndex << std::endl;
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the curve
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        KRATOS_DEBUG_ERROR_IF(LocalDirectionIndex != 0)
            << "Trying to access polynomial degree in direction " << LocalDirectionIndex
            << " from NurbsCurveGeometry #" << this->Id() << ". However, nurbs curves have only one direction."
            << std::endl;

        return mPolynomialDegree;
    }

    ///@}
    ///@name Get and Set functions
    ///@{

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
        return mPolynomialDegree + 1;
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

    /// Returns the number of spans.
    SizeType NumberOfKnotSpans(IndexType DirectionIndex = 0) const
    {
        SizeType knot_span_counter = 0;
        for (IndexType i = 0; i < mKnots.size() - 1; i++) {
            if (std::abs(mKnots[i] - mKnots[i + 1]) > 1e-6) {
                knot_span_counter++;
            }
        }
        return knot_span_counter;
    }

    /* @brief Provides knot spans of this nurbs curve.
     * @param resulting vector of span intervals.
     * @param index of chosen direction, for curves always 0.
     */
    void Spans(std::vector<double>& rSpans, IndexType DirectionIndex = 0) const
    {
        rSpans.resize(this->NumberOfKnotSpans(DirectionIndex) + 1);

        rSpans[0] = mKnots[0];

        IndexType counter = 1;
        for (IndexType i = 0; i < mKnots.size() - 1; i++) {
            if (std::abs(mKnots[i] - mKnots[i + 1]) > 1e-6) {
                rSpans[counter] = mKnots[i + 1];
                counter++;
            }
        }
    }

    ///@}
    ///@name Geometrical Informations
    ///@{

    /// Computes the length of a nurbs curve
    double Length() const override
    {
        IntegrationPointsArrayType integration_points;
        CreateIntegrationPoints(integration_points);

        double length = 0.0;
        for (IndexType i = 0; i < integration_points.size(); ++i) {
            double determinant_jacobian = this->DeterminantOfJacobian(integration_points[i]);
            length += integration_points[i].Weight() * determinant_jacobian;
        }
        return length;
    }

    ///@}
    ///@name Jacobian
    ///@{

    /// Computes the Jacobian with coordinates
    Matrix& Jacobian(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        const SizeType working_space_dimension = this->WorkingSpaceDimension();
        if (rResult.size1() != working_space_dimension || rResult.size2() != 1) {
            rResult.resize(working_space_dimension, 1, false);
        }
        rResult.clear();

        /// Compute shape functions
        NurbsCurveShapeFunction shape_function_container(mPolynomialDegree, 1);
        if (IsRational()) {
            shape_function_container.ComputeNurbsShapeFunctionValues(mKnots, mWeights, rCoordinates[0]);
        }
        else {
            shape_function_container.ComputeBSplineShapeFunctionValues(mKnots, rCoordinates[0]);
        }

        /// Compute Jacobian
        for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); i++) {
            const IndexType index = shape_function_container.GetFirstNonzeroControlPoint() + i;

            const array_1d<double, 3>& r_coordinates = (*this)[index].Coordinates();
            for (IndexType k = 0; k < working_space_dimension; ++k) {
                rResult(k, 0) += r_coordinates[k] * shape_function_container(i, 1);
            }
        }
        return rResult;
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points according to its the polynomial degrees.
     * @param result integration points.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints) const override
    {
        const SizeType points_per_span = PolynomialDegree(0) + 1;

        std::vector<double> spans;
        Spans(spans);

        this->CreateIntegrationPoints(
            rIntegrationPoints, spans, points_per_span);
    }

    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpanIntervals,
        SizeType IntegrationPointsPerSpan) const
    {
        const SizeType num_spans = rSpanIntervals.size() - 1;
        const SizeType number_of_integration_points =
            num_spans * IntegrationPointsPerSpan;

        if (rIntegrationPoints.size() != number_of_integration_points)
            rIntegrationPoints.resize(number_of_integration_points);

        typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();

        for (IndexType i = 0; i < num_spans; ++i)
        {
            IntegrationPointUtilities::IntegrationPoints1D(
                integration_point_iterator,
                IntegrationPointsPerSpan,
                rSpanIntervals[i], rSpanIntervals[i + 1]);
        }
    }

    ///@}
    ///@name Quadrature Points
    ///@{

    /* @brief creates a list of quadrature point geometries
     *        from a list of integration points.
     *
     * @param rResultGeometries list of quadrature point geometries.
     * @param rIntegrationPoints list of provided integration points.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints) override
    {
        // shape function container.
        NurbsCurveShapeFunction shape_function_container(
            mPolynomialDegree, NumberOfShapeFunctionDerivatives);

        // Resize containers.
        if (rResultGeometries.size() != rIntegrationPoints.size())
            rResultGeometries.resize(rIntegrationPoints.size());

        auto default_method = this->GetDefaultIntegrationMethod();
        SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();

        Matrix N(1, num_nonzero_cps);
        DenseVector<Matrix> shape_function_derivatives(NumberOfShapeFunctionDerivatives - 1);
        for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; i++) {
            shape_function_derivatives[i].resize(num_nonzero_cps, 1);
        }

        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i)
        {
            if (IsRational()) {
                shape_function_container.ComputeNurbsShapeFunctionValues(
                    mKnots, mWeights, rIntegrationPoints[i][0]);
            }
            else {
                shape_function_container.ComputeBSplineShapeFunctionValues(
                    mKnots, rIntegrationPoints[i][0]);
            }

            /// Get List of Control Points
            PointsArrayType nonzero_control_points(num_nonzero_cps);
            auto first_cp_index = shape_function_container.GetFirstNonzeroControlPoint();
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                nonzero_control_points(j) = pGetPoint(first_cp_index + j);
            }
            /// Get Shape Functions N
            if (NumberOfShapeFunctionDerivatives >= 0) {
                for (IndexType j = 0; j < num_nonzero_cps; j++) {
                    N(0, j) = shape_function_container(j, 0);
                }
            }

            /// Get Shape Function Derivatives DN_De, ...
            if (NumberOfShapeFunctionDerivatives > 0) {
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; n++) {
                    for (IndexType j = 0; j < num_nonzero_cps; j++) {
                        shape_function_derivatives[n](j, 0) = shape_function_container(j, n + 1);
                    }
                }
            }

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, rIntegrationPoints[i],
                N, shape_function_derivatives);

            rResultGeometries(i) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePoint(
                this->WorkingSpaceDimension(), 2, data_container, nonzero_control_points);
        }
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
        return std::to_string(TWorkingSpaceDimension) + " dimensional nurbs curve.";
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


    /*
    * @brief Checks if the knot vector is coinciding with the number of
    *        control points and the polynomial degree. If the knot vector
    *        has a multiplicity of p+1 in the beginning, it is reduced to p.
    */
    void CheckAndFitKnotVectors()
    {
        SizeType num_control_points = this->size();

        if (mKnots.size() != NurbsUtilities::GetNumberOfKnots(mPolynomialDegree, num_control_points)) {
            if ((mKnots.size() - 2) == NurbsUtilities::GetNumberOfKnots(mPolynomialDegree, num_control_points)) {
                Vector Knots = ZeroVector(mKnots.size() - 2);
                for (SizeType i = 0; i < mKnots.size() - 2; ++i) {
                    Knots[i] = mKnots[i + 1];
                }
                mKnots = Knots;
            } else {
                KRATOS_ERROR
                    << "Number of controls points, polynomial degree and number of knots do not match! " << std::endl
                    << " P: " << mPolynomialDegree << ", size of knot vector: " << mKnots.size()
                    << ", number of control points: " << num_control_points << "." << std::endl
                    << "Following condition must be achieved: Knots.size() = (ControlPoints.size() + PolynomialDegree - 1)."
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
