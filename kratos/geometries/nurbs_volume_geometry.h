//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_NURBS_VOLUME_GEOMETRY_H_INCLUDED )
#define  KRATOS_NURBS_VOLUME_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"

#include "utilities/quadrature_points_utility.h"

#include "integration/integration_point_utilities.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class NurbsVolumeGeometry
 * @ingroup KratosCore
 * @brief A volume geometry based on a full 3-dimensional BSpline tensor product.
 * @details Weights are not yet implemented. Therefore no NURBS, but only BSplines are constructed.
 *
 * @author Manuel Messmer
 **/
template <class TContainerPointType>
class NurbsVolumeGeometry : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename TContainerPointType::value_type NodeType;

    typedef Geometry<NodeType> BaseType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    // Using base class functionalities.
    using BaseType::pGetPoint;
    using BaseType::GetPoint;

    /// Pointer of NurbsVolumeGeometry
    KRATOS_CLASS_POINTER_DEFINITION(NurbsVolumeGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Conctructor for B-Spline volumes
    NurbsVolumeGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const SizeType PolynomialDegreeW,
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rKnotsW)
        : BaseType(rThisPoints, &msGeometryData)
        , mPolynomialDegreeU(PolynomialDegreeU)
        , mPolynomialDegreeV(PolynomialDegreeV)
        , mPolynomialDegreeW(PolynomialDegreeW)
        , mKnotsU(rKnotsU)
        , mKnotsV(rKnotsV)
        , mKnotsW(rKnotsW)
    {
        CheckAndFitKnotVectors();
    }

    /// Attention: Weigths are not yet implemented!
    /// Conctructor for NURBS volumes
    // NurbsVolumeGeometry(
    //     const PointsArrayType& rThisPoints,
    //     const SizeType PolynomialDegreeU,
    //     const SizeType PolynomialDegreeV,
    //     const SizeType PolynomialDegreeW,
    //     const Vector& rKnotsU,
    //     const Vector& rKnotsV,
    //     const Vector& rKnotsW,
    //     const Vector& rWeights)
    //     : BaseType(rThisPoints, &msGeometryData)
    //     , mPolynomialDegreeU(PolynomialDegreeU)
    //     , mPolynomialDegreeV(PolynomialDegreeV)
    //     , mPolynomialDegreeW(PolynomialDegreeW)
    //     , mKnotsU(rKnotsU)
    //     , mKnotsV(rKnotsV)
    //     , mKnotsW(rKnotsW)
    //     , mWeights(rWeights)
    // {
    //     CheckAndFitKnotVectors();

    //     KRATOS_ERROR_IF(rWeights.size() != rThisPoints.size())
    //         << "Number of control points and weights do not match!" << std::endl;
    // }

    explicit NurbsVolumeGeometry(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    NurbsVolumeGeometry(NurbsVolumeGeometry<TContainerPointType> const& rOther)
        : BaseType(rOther)
        , mPolynomialDegreeU(rOther.mPolynomialDegreeU)
        , mPolynomialDegreeV(rOther.mPolynomialDegreeV)
        , mPolynomialDegreeW(rOther.mPolynomialDegreeW)
        , mKnotsU(rOther.mKnotsU)
        , mKnotsV(rOther.mKnotsV)
        , mKnotsW(rOther.mKnotsW)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherContainerPointType> NurbsVolumeGeometry(
        NurbsVolumeGeometry<TOtherContainerPointType> const& rOther)
        : BaseType(rOther, &msGeometryData)
        , mPolynomialDegreeU(rOther.mPolynomialDegreeU)
        , mPolynomialDegreeV(rOther.mPolynomialDegreeV)
        , mPolynomialDegreeW(rOther.mPolynomialDegreeW)
        , mKnotsU(rOther.mKnotsU)
        , mKnotsV(rOther.mKnotsV)
        , mKnotsW(rOther.mKnotsW)
    {
    }

    /// Destructor.
    ~NurbsVolumeGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    NurbsVolumeGeometry& operator=(const NurbsVolumeGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mPolynomialDegreeU = rOther.mPolynomialDegreeU;
        mPolynomialDegreeV = rOther.mPolynomialDegreeV;
        mPolynomialDegreeW = rOther.mPolynomialDegreeW;
        mKnotsU = rOther.mKnotsU;
        mKnotsV = rOther.mKnotsV;
        mKnotsW = rOther.mKnotsW;
        return *this;
    }

    /// Assignment operator for geometries with different point type
    template<class TOtherContainerPointType>
    NurbsVolumeGeometry& operator=(
        NurbsVolumeGeometry<TOtherContainerPointType> const & rOther)
    {
        BaseType::operator=(rOther);
        mPolynomialDegreeU = rOther.mPolynomialDegreeU;
        mPolynomialDegreeV = rOther.mPolynomialDegreeV;
        mPolynomialDegreeW = rOther.mPolynomialDegreeW;
        mKnotsU = rOther.mKnotsU;
        mKnotsV = rOther.mKnotsV;
        mKnotsW = rOther.mKnotsW;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        PointsArrayType const& ThisPoints) const override
    {
        return Kratos::make_shared<NurbsVolumeGeometry>(ThisPoints);
    }

    ///@}
    ///@name Geometrical Information
    ///@{

    /// Returns number of points per direction.
    SizeType PointsNumberInDirection(IndexType LocalDirectionIndex) const override
    {
        if (LocalDirectionIndex == 0) {
            return this->NumberOfControlPointsU();
        }
        else if (LocalDirectionIndex == 1) {
            return this->NumberOfControlPointsV();
        }
        else if (LocalDirectionIndex == 2) {
            return this->NumberOfControlPointsW();
        }
        KRATOS_ERROR << "Possible direction index in NurbsVolumeGeometry reaches from 0-2. Given direction index: "
            << LocalDirectionIndex << std::endl;
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the volume in direction 0, 1, 2
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        KRATOS_DEBUG_ERROR_IF(LocalDirectionIndex > 2)
            << "Trying to access polynomial degree in direction " << LocalDirectionIndex
            << " from NurbsVolumeGeometry #" << this->Id() << ". Nurbs volume have only three directions."
            << std::endl;

        if (LocalDirectionIndex == 0) {
            return mPolynomialDegreeU;
        }
        else if(LocalDirectionIndex == 1){
            return mPolynomialDegreeV;
        }
        else {
            return mPolynomialDegreeW;
        }
    }

    ///@}
    ///@name Get and Set functions
    ///@{

    void SetInternals(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const SizeType PolynomialDegreeW,
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rKnotsW)
    {
        this->Points() = rThisPoints;
        mPolynomialDegreeU = PolynomialDegreeU;
        mPolynomialDegreeV = PolynomialDegreeV;
        mPolynomialDegreeW = PolynomialDegreeW;
        mKnotsU = rKnotsU;
        mKnotsV = rKnotsV;
        mKnotsW = rKnotsW;

        CheckAndFitKnotVectors();
    }

    /**
     * @return returns the polynomial degree 'p' in u direction.
     **/
    SizeType PolynomialDegreeU() const
    {
        return mPolynomialDegreeU;
    }

    /**
     *  @return returns the polynomial degree 'p' in v direction.
     **/
    SizeType PolynomialDegreeV() const
    {
        return mPolynomialDegreeV;
    }

    /**
     *  @return returns the polynomial degree 'p' in w direction.
     **/
    SizeType PolynomialDegreeW() const
    {
        return mPolynomialDegreeW;
    }

    /**
     * @brief Get Knot vector in u-direction.
     * @details This vector is defined to have a multiplicity of p at the beginning and end (NOT: p + 1).
     * @return knot vector.
     **/
    const Vector& KnotsU() const
    {
        return mKnotsU;
    }

    /**
     * @brief Get Knot vector in v-direction.
     * @details This vector is defined to have a multiplicity of p at the beginning and end (NOT: p + 1).
     * @return knot vector.
     **/
    const Vector& KnotsV() const
    {
        return mKnotsV;
    }

    /**
     * @brief Get Knot vector in w-direction.
     * @details This vector is defined to have a multiplicity of p at the beginning and end (NOT: p + 1).
     * @return knot vector.
     **/
    const Vector& KnotsW() const
    {
        return mKnotsW;
    }

    /**
     * @return Gives the size of the knot vector in u-direction.
     **/
    SizeType NumberOfKnotsU() const
    {
        return mKnotsU.size();
    }

    /**
     * @return Gives the size of the knot vector in v-direction.
     **/
    SizeType NumberOfKnotsV() const
    {
        return mKnotsV.size();
    }

    /**
     * @return Gives the size of the knot vector in w-direction.
     **/
    SizeType NumberOfKnotsW() const
    {
        return mKnotsW.size();
    }

    /**
     * @brief Checks if shape functions are rational or not.
     * @return False. Weights are not yet considered.
     **/
    bool IsRational() const
    {
        return false;
    }

    /// Attention weights are not yet implemented.
    /* Get Weights vector. All values are 1.0 for B-Splines, for NURBS those can be unequal 1.0.
    @return weights vector.
    */
    // const Vector& Weights() const
    // {
    //     return mWeights;
    // }

    /**
     * @return Gives the number of control points in u-direction.
     **/
    SizeType NumberOfControlPointsU() const
    {
        return NumberOfKnotsU() - PolynomialDegreeU() + 1;
    }

    /**
     * @return Gives the number of control points in v-direction.
     **/
    SizeType NumberOfControlPointsV() const
    {
        return NumberOfKnotsV() - PolynomialDegreeV() + 1;
    }

    /**
     * @return Gives the number of control points in w-direction.
     **/
    SizeType NumberOfControlPointsW() const
    {
        return NumberOfKnotsW() - PolynomialDegreeW() + 1;
    }

    /**
     * @brief Returns the number of spans in DirectionIndex=0:U, DirectionIndex=1:V and DirectionIndex=2
     * (which are larger than 0).
     **/
    SizeType NumberOfKnotSpans(IndexType DirectionIndex) const
    {
        SizeType knot_span_counter = 0;
        if (DirectionIndex == 0) {
            for (IndexType i = 0; i < mKnotsU.size() - 1; ++i) {
                if (std::abs(mKnotsU[i] - mKnotsU[i + 1]) > 1e-6) {
                    knot_span_counter++;
                }
            }
        }
        else if (DirectionIndex == 1) {
            for (IndexType i = 0; i < mKnotsV.size() - 1; ++i) {
                if (std::abs(mKnotsV[i] - mKnotsV[i + 1]) > 1e-6) {
                    knot_span_counter++;
                }
            }
        }
        else if (DirectionIndex == 2){
            for( IndexType i = 0; i < mKnotsW.size() - 1; ++i){
                if( std::abs(mKnotsW[i] - mKnotsW[i+1]) > 1e-6) {
                    knot_span_counter++;
                }
            }

        } else {
            KRATOS_ERROR << "NurbsVolumeGeometry::NumberOfKnotSpans: Direction index: "
                << DirectionIndex << " not available. Options are: 0, 1 and 2." << std::endl;
        }
        return knot_span_counter;
    }

    /**
     * @brief Provides all knot spans along the given direction.
     * @param return vector of span intervals.
     * @param index of direction: 0: U; 1: V; 2: W;
     **/
    void Spans(std::vector<double>& rSpans, IndexType DirectionIndex) const
    {
        rSpans.resize(this->NumberOfKnotSpans(DirectionIndex) + 1);

        if (DirectionIndex == 0) {
            rSpans[0] = mKnotsU[0];

            IndexType counter = 1;
            for (IndexType i = 0; i < mKnotsU.size() - 1; ++i) {
                if (std::abs(mKnotsU[i] - mKnotsU[i + 1]) > 1e-6) {
                    rSpans[counter] = mKnotsU[i + 1];
                    counter++;
                }
            }
        }
        else if (DirectionIndex == 1) {
            rSpans[0] = mKnotsV[0];

            IndexType counter = 1;
            for (IndexType i = 0; i < mKnotsV.size() - 1; ++i) {
                if (std::abs(mKnotsV[i] - mKnotsV[i + 1]) > 1e-6) {
                    rSpans[counter] = mKnotsV[i + 1];
                    counter++;
                }
            }
        }
        else if (DirectionIndex == 2) {
            rSpans[0] = mKnotsW[0];

            IndexType counter = 1;
            for (IndexType i = 0; i < mKnotsW.size() - 1; ++i) {
                if (std::abs(mKnotsW[i] - mKnotsW[i + 1]) > 1e-6) {
                    rSpans[counter] = mKnotsW[i + 1];
                    counter++;
                }
            }
        } else {
            KRATOS_ERROR << "NurbsVolumeGeometry::Spans: Direction index: "
                << DirectionIndex << " not available. Options are: 0, 1 and 2." << std::endl;
        }
    }

    /**
     * @brief Provides the natural boundaries of the NURBS/B-Spline volume.
     * @return Domain interval.
     **/
    NurbsInterval DomainIntervalU() const
    {
        return NurbsInterval(
            mKnotsU[mPolynomialDegreeU - 1],
            mKnotsU[NumberOfKnotsU() - mPolynomialDegreeU]);
    }

    /**
     * @brief Provides the natural boundaries of the NURBS/B-Spline volume.
     * @return Domain interval.
     **/
    NurbsInterval DomainIntervalV() const
    {
        return NurbsInterval(
            mKnotsV[mPolynomialDegreeV - 1],
            mKnotsV[NumberOfKnotsV() - mPolynomialDegreeV]);
    }

    /**
     * @brief Provides the natural boundaries of the NURBS/B-Spline volume.
     * @return Domain interval.
     **/
    NurbsInterval DomainIntervalW() const
    {
        return NurbsInterval(
            mKnotsW[mPolynomialDegreeW - 1],
            mKnotsW[NumberOfKnotsW() - mPolynomialDegreeW]);
    }

    /**
     * @brief Provides all knot span intervals of the volume in u-direction.
     * @return Vector of knot span intervals.
     **/
    std::vector<NurbsInterval> KnotSpanIntervalsU() const
    {
        const SizeType first_span = mPolynomialDegreeU - 1;
        const SizeType last_span = NumberOfKnotsU() - mPolynomialDegreeU - 1;

        const SizeType number_of_spans = last_span - first_span + 1;

        std::vector<NurbsInterval> result(number_of_spans);

        for (IndexType i = 0; i < number_of_spans; i++) {
            const double t0 = mKnotsU[first_span + i];
            const double t1 = mKnotsU[first_span + i + 1];

            result[i] = NurbsInterval(t0, t1);
        }

        return result;
    }

    /**
     * @brief Provides all knot span intervals of the volume in v-direction.
     * @return Vector of knot span intervals.
     **/
    std::vector<NurbsInterval> KnotSpanIntervalsV() const
    {
        const SizeType first_span = mPolynomialDegreeV - 1;
        const SizeType last_span = NumberOfKnotsV() - mPolynomialDegreeV - 1;

        const SizeType number_of_spans = last_span - first_span + 1;

        std::vector<NurbsInterval> result(number_of_spans);

        for (IndexType i = 0; i < number_of_spans; i++) {
            const double t0 = mKnotsV[first_span + i];
            const double t1 = mKnotsV[first_span + i + 1];

            result[i] = NurbsInterval(t0, t1);
        }

        return result;
    }

    /**
     * @brief Provides all knot span intervals of the volume in w-direction.
     * @return Vector of knot span intervals.
     **/
    std::vector<NurbsInterval> KnotSpanIntervalsW() const
    {
        const SizeType first_span = mPolynomialDegreeW - 1;
        const SizeType last_span = NumberOfKnotsW() - mPolynomialDegreeW - 1;

        const SizeType number_of_spans = last_span - first_span + 1;

        std::vector<NurbsInterval> result(number_of_spans);

        for (IndexType i = 0; i < number_of_spans; i++) {
            const double t0 = mKnotsW[first_span + i];
            const double t1 = mKnotsW[first_span + i + 1];

            result[i] = NurbsInterval(t0, t1);
        }

        return result;
    }

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree of the underlying surface.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        return IntegrationInfo(
            { PolynomialDegreeU() + 1, PolynomialDegreeV() + 1, PolynomialDegreeW() + 1 },
            { IntegrationInfo::QuadratureMethod::GAUSS, IntegrationInfo::QuadratureMethod::GAUSS, IntegrationInfo::QuadratureMethod::GAUSS });
    }

    ///@}
    ///@name Integration Points
    ///@{

    /**
     * @brief Creates integration points according to the polynomial degrees.
     * @return integration points.
     **/
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        const SizeType points_in_u = PolynomialDegreeU() + 1;
        const SizeType points_in_v = PolynomialDegreeV() + 1;
        const SizeType points_in_w = PolynomialDegreeW() + 1;

        CreateIntegrationPoints(
            rIntegrationPoints, points_in_u, points_in_v, points_in_w);
    }

    /**
     * @brief Creates integration points according to the input parameter.
     * @param return integration points.
     * @param points Number of points in u direction per span.
     * @param points Number of points in v direction per span.
     * @param points Number of points in w direction per span.
     **/
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        SizeType NumPointsPerSpanU,
        SizeType NumPointsPerSpanV,
        SizeType NumPointsPerSpanW) const
    {
        auto knot_span_intervals_u = KnotSpanIntervalsU();
        auto knot_span_intervals_v = KnotSpanIntervalsV();
        auto knot_span_intervals_w = KnotSpanIntervalsW();

        const SizeType number_of_integration_points =
            knot_span_intervals_u.size() * knot_span_intervals_v.size() * knot_span_intervals_w.size()
            * NumPointsPerSpanU * NumPointsPerSpanV * NumPointsPerSpanW;

        if (rIntegrationPoints.size() != number_of_integration_points) {
            rIntegrationPoints.resize(number_of_integration_points);
        }

        typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();

        for (IndexType i = 0; i < knot_span_intervals_u.size(); ++i) {
            for (IndexType j = 0; j < knot_span_intervals_v.size(); ++j) {
                for (IndexType k = 0; k < knot_span_intervals_w.size(); ++k) {
                    IntegrationPointUtilities::IntegrationPoints3D(
                        integration_point_iterator,
                        NumPointsPerSpanU, NumPointsPerSpanV, NumPointsPerSpanW,
                        knot_span_intervals_u[i].GetT0(), knot_span_intervals_u[i].GetT1(),
                        knot_span_intervals_v[j].GetT0(), knot_span_intervals_v[j].GetT1(),
                        knot_span_intervals_w[k].GetT0(), knot_span_intervals_w[k].GetT1());
                }
            }
        }
    }

    /**
     * @brief Computes jacobian matrix at the given coordinates.
     * @param rCoordinates Coordinates to be evaluated.
     * @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.
     * @todo Refactor such that addional 'ComputeBSplineShapeFunctionValues'-call can be omitted. Here it is only called to
     *       find the correct knotspans and to set the shape_function_member variables 'mFirstNonzeroControlPointU,-V,-W'.
     * @note This function is only required to compute e.g. the volume of the geometry. During an IGA-Analysis the corresponding function
     *       of the base class is called.
     **/
    Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rCoordinates ) const override
    {
        const SizeType working_space_dimension = this->WorkingSpaceDimension();
        const SizeType local_space_dimension = this->LocalSpaceDimension();
        const SizeType points_number = this->PointsNumber();

        rResult = ZeroMatrix(working_space_dimension,local_space_dimension);

        Matrix shape_functions_gradients(points_number, local_space_dimension);
        ShapeFunctionsLocalGradients( shape_functions_gradients, rCoordinates );

        // Get control point indices
        NurbsVolumeShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, mPolynomialDegreeW, 0);
        shape_function_container.ComputeBSplineShapeFunctionValues(mKnotsU,mKnotsV,mKnotsW,
            rCoordinates[0], rCoordinates[1], rCoordinates[2]);

        SizeType number_cp_u = NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeU, mKnotsU.size());
        SizeType number_cp_v = NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeV, mKnotsV.size());
        SizeType number_cp_w = NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeW, mKnotsW.size());

        std::vector<int> cp_indices = shape_function_container.ControlPointIndices(number_cp_u, number_cp_v, number_cp_w);
        SizeType number_of_cp = shape_function_container.NumberOfNonzeroControlPoints();

        for (IndexType i = 0; i < number_of_cp; ++i ) {
            const array_1d<double, 3>& r_coordinates = (*this)[cp_indices[i]].Coordinates();
            for(IndexType k = 0; k< working_space_dimension; ++k) {
                const double value = r_coordinates[k];
                for(IndexType m = 0; m < local_space_dimension; ++m) {
                    rResult(k,m) += value * shape_functions_gradients(i,m);
                }
            }
        }

        return rResult;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a list of quadrature point geometries
     *        from a list of integration points.
     * @param rResultGeometries List of quadrature point geometries.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     * @param rIntegrationInfo.
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        IntegrationInfo& rIntegrationInfo) override
    {
        IntegrationPointsArrayType IntegrationPoints;
        CreateIntegrationPoints(IntegrationPoints, rIntegrationInfo);

        // Makes sure we use assembly Option 2 in CreateQuadraturePointGeometries().
        IntegrationInfo integration_info(
            { PolynomialDegreeU() + 1, PolynomialDegreeV() + 1, PolynomialDegreeW() + 1 },
            { IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS });

        this->CreateQuadraturePointGeometries(
            rResultGeometries,
            NumberOfShapeFunctionDerivatives,
            IntegrationPoints,
            integration_info);
    }

    /**
     * @brief Creates a list of quadrature point geometries.
     *        from a list of integration points.
     * @param rResultGeometries List of quadrature point geometries.
     * @param rIntegrationPoints List of provided integration points.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     * @param rIntegrationInfo.
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        // Option 1: One QuadraturePointGeometry is created containing all integration points. This should be used
        // when all rIntegrationPoints are located inside the same element.
        if(  IntegrationInfo::QuadratureMethod::GAUSS == rIntegrationInfo.GetQuadratureMethod(0) ){
            KRATOS_ERROR_IF(NumberOfShapeFunctionDerivatives != 2) << "NumberOfShapeFunctionDerivatives must be 2.\n";
            // Shape function container.
            NurbsVolumeShapeFunction shape_function_container(
                mPolynomialDegreeU, mPolynomialDegreeV, mPolynomialDegreeW, NumberOfShapeFunctionDerivatives);

            // Resize containers.
            if (rResultGeometries.size() != 1){
                rResultGeometries.resize(1);
            }

            const auto default_method = this->GetDefaultIntegrationMethod();

            const SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();
            const SizeType num_points = rIntegrationPoints.size();
            KRATOS_ERROR_IF(num_points < 1) << "List of integration points is empty.\n";

            // Initialize containers.
            IntegrationPointsContainerType integration_points;
            ShapeFunctionsValuesContainerType shape_function_values;
            ShapeFunctionsLocalGradientsContainerType shape_function_gradients;

            integration_points[0] = rIntegrationPoints;
            shape_function_gradients[0].resize(rIntegrationPoints.size());
            shape_function_values[0].resize(rIntegrationPoints.size(), num_nonzero_cps);

            for( IndexType i_point = 0; i_point < num_points; ++i_point){
                shape_function_gradients[0][i_point].resize(num_nonzero_cps, 3);
            }

            // Centroid of points. This will be used to identify knot span.
            // Single point might be located on the boundary of two knot spans.
            array_1d<double, 3> centroid(3, 0.0);

            // Fill containers
            for (IndexType i_point = 0; i_point < num_points; ++i_point)
            {
                // Compute centroid.
                centroid += rIntegrationPoints[i_point].Coordinates();

                shape_function_container.ComputeBSplineShapeFunctionValues(
                    mKnotsU, mKnotsV, mKnotsW,
                    rIntegrationPoints[i_point][0], rIntegrationPoints[i_point][1], rIntegrationPoints[i_point][2]);

                /// Get Shape Functions.
                for (IndexType j = 0; j < num_nonzero_cps; ++j) {
                    shape_function_values[0](i_point, j) = shape_function_container(j, 0);
                }

                // Get Shape Function Derivatives DN_De, ...
                for (IndexType k = 0; k < 3; ++k) {
                    for (IndexType j = 0; j < num_nonzero_cps; ++j) {
                        shape_function_gradients[0][i_point](j, k) = shape_function_container(j, 1 + k);
                    }
                }
            }

            /// Get List of Control Points.
            PointsArrayType nonzero_control_points(num_nonzero_cps);
            centroid /= num_points;
            shape_function_container.ComputeBSplineShapeFunctionValues(
                    mKnotsU, mKnotsV, mKnotsW,
                    centroid[0], centroid[1], centroid[2]);

            auto cp_indices = shape_function_container.ControlPointIndices(
                NumberOfControlPointsU(), NumberOfControlPointsV(), NumberOfControlPointsW());

            for (IndexType j = 0; j < num_nonzero_cps; ++j) {
                nonzero_control_points(j) = pGetPoint(cp_indices[j]);
            }

            // Instantiate shape function container.
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, integration_points,
                shape_function_values, shape_function_gradients);

            // Create quadrature point geometry.
            rResultGeometries(0) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePoint(
                this->WorkingSpaceDimension(), 3, data_container, nonzero_control_points, this);
        }
        // Option 2: A list of QuadraturePointGeometry is created, one for each integration points.
        else if ( IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS == rIntegrationInfo.GetQuadratureMethod(0) ) {
            // Shape function container.
            NurbsVolumeShapeFunction shape_function_container(
                mPolynomialDegreeU, mPolynomialDegreeV, mPolynomialDegreeW, NumberOfShapeFunctionDerivatives);

            // Resize containers.
            if (rResultGeometries.size() != rIntegrationPoints.size())
                rResultGeometries.resize(rIntegrationPoints.size());

            auto default_method = this->GetDefaultIntegrationMethod();
            SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();

            Matrix N(1, num_nonzero_cps);
            DenseVector<Matrix> shape_function_derivatives(NumberOfShapeFunctionDerivatives - 1);

            for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; ++i) {
                const IndexType num_derivatives = (2 + i) * (3 + i) / 2;
                shape_function_derivatives[i].resize(num_nonzero_cps, num_derivatives);
            }

            for (IndexType i = 0; i < rIntegrationPoints.size(); ++i)
            {
                shape_function_container.ComputeBSplineShapeFunctionValues(
                    mKnotsU, mKnotsV, mKnotsW, rIntegrationPoints[i][0], rIntegrationPoints[i][1], rIntegrationPoints[i][2]);


                /// Get List of Control Points
                PointsArrayType nonzero_control_points(num_nonzero_cps);
                auto cp_indices = shape_function_container.ControlPointIndices(
                    NumberOfControlPointsU(), NumberOfControlPointsV(), NumberOfControlPointsW());
                for (IndexType j = 0; j < num_nonzero_cps; ++j) {
                    nonzero_control_points(j) = pGetPoint(cp_indices[j]);
                }
                /// Get Shape Functions N
                for (IndexType j = 0; j < num_nonzero_cps; ++j) {
                    N(0, j) = shape_function_container(j, 0);
                }

                /// Get Shape Function Derivatives DN_De, ...
                if (NumberOfShapeFunctionDerivatives > 0) {
                    IndexType shape_derivative_index = 1;
                    for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n) {
                        const IndexType num_derivatives = (2 + n) * (3 + n) / 2;
                        for (IndexType k = 0; k < num_derivatives; ++k) {
                            for (IndexType j = 0; j < num_nonzero_cps; ++j) {
                                shape_function_derivatives[n](j, k) = shape_function_container(j, shape_derivative_index + k);
                            }
                        }
                        shape_derivative_index += num_derivatives;
                    }
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                    default_method, rIntegrationPoints[i],
                    N, shape_function_derivatives);

                rResultGeometries(i) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePoint(
                    this->WorkingSpaceDimension(), 3, data_container, nonzero_control_points, this);
            }
        } else {
            KRATOS_ERROR << "Integration method not available.\n";
        }
    }

    ///@}
    ///@name Operations
    ///@{

/**
    * @brief Returns local coordinates in return for a point in physical coordinates.
    *        This function assumes the initial geometry to be a regular grid. Only linear mapping is applied.
    *
    * @param rPointGlobalCoordinates the point to which the
    *        projection has to be found.
    * @param rProjectedPointLocalCoordinates the location of the
    *        projection in local coordinates.
    * @param Tolerance not used.
    * @return This functions always returns 1 (successful projection).
    *         Since no Netwon-Raphson is performed, the projection is always considered as successful.
    */
    int ProjectionPointGlobalToLocalSpace(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const override
    {
        const CoordinatesArrayType& lower_point = this->begin()->GetInitialPosition();
        const CoordinatesArrayType& upper_point = (this->end()-1)->GetInitialPosition();

        const Vector& knots_u = this->KnotsU();
        const SizeType nb_knots_u = this->NumberOfKnotsU();
        const Vector& knots_v = this->KnotsV();
        const SizeType nb_knots_v = this->NumberOfKnotsV();
        const Vector& knots_w = this->KnotsW();
        const SizeType nb_knots_w = this->NumberOfKnotsW();

        rProjectedPointLocalCoordinates[0] = ((rPointGlobalCoordinates[0] - lower_point[0]) / std::abs( lower_point[0] - upper_point[0])
                        * std::abs(knots_u[nb_knots_u-1] - knots_u[0])) + knots_u[0];
        rProjectedPointLocalCoordinates[1] = ((rPointGlobalCoordinates[1] - lower_point[1]) / std::abs( lower_point[1] - upper_point[1])
                        * std::abs(knots_v[nb_knots_v-1] - knots_v[0])) + knots_v[0];
        rProjectedPointLocalCoordinates[2] = ((rPointGlobalCoordinates[2] - lower_point[2]) / std::abs( lower_point[2] - upper_point[2])
                        * std::abs(knots_w[nb_knots_w-1] - knots_w[0])) + knots_w[0];
        return 1;
    }

    /**
     * @brief This method maps from local space to working space.
     * @param LocalCoordinates The local coordinates in dimension space.
     * @return array_1d<double, 3> with the coordinates in working space.
     **/
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates
    ) const override
    {
        NurbsVolumeShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, mPolynomialDegreeW, 0);

        // Attention: Weights are not yet implemented.
        // if (IsRational()) {
        //     shape_function_container.ComputeNurbsShapeFunctionValues(
        //         mKnotsU, mKnotsV, mKnotsW, mWeights, rLocalCoordinates[0], rLocalCoordinates[1], rLocalCoordinates[2]);
        // }

        shape_function_container.ComputeBSplineShapeFunctionValues(
            mKnotsU, mKnotsV, mKnotsW, rLocalCoordinates[0], rLocalCoordinates[1], rLocalCoordinates[2]);


        noalias(rResult) = ZeroVector(3);
        for (IndexType u = 0; u <= PolynomialDegreeU(); ++u) {
            for (IndexType v = 0; v <= PolynomialDegreeV(); ++v) {
                for(IndexType w = 0; w <= PolynomialDegreeW(); ++w) {
                    IndexType cp_index_u = shape_function_container.GetFirstNonzeroControlPointU() + u;
                    IndexType cp_index_v = shape_function_container.GetFirstNonzeroControlPointV() + v;
                    IndexType cp_index_w = shape_function_container.GetFirstNonzeroControlPointW() + w;
                    const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        NumberOfControlPointsU(), NumberOfControlPointsV(), NumberOfControlPointsW(),
                            cp_index_u, cp_index_v, cp_index_w);
                    rResult += (*this)[index] * shape_function_container(u, v, w, 0);
                }
            }
        }

        return rResult;
    }

    /**
     * @brief This method maps from local space to working space and computes the
     *        derivatives at the local space parameter in the dimension of the object.
     * @param LocalCoordinates The local coordinates to be evaluated.
     * @param Derivative Number of computed derivatives.
     * @return std::vector<array_1d<double, 3>> with the coordinates and derivatives in working space.
     **/
    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const override
    {
        NurbsVolumeShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, mPolynomialDegreeW, DerivativeOrder);

        // Attention: Weights are not yet implemented.
        // if (IsRational()) {
        //     shape_function_container.ComputeNurbsShapeFunctionValues(
        //         mKnotsU, mKnotsV, mKnotsW, mWeights, rLocalCoordinates[0], rLocalCoordinates[1], rLocalCoordinates[2]);
        // }

        shape_function_container.ComputeBSplineShapeFunctionValues(
            mKnotsU, mKnotsV, mKnotsW, rLocalCoordinates[0], rLocalCoordinates[1], rLocalCoordinates[2]);


        if (rGlobalSpaceDerivatives.size() != shape_function_container.NumberOfShapeFunctionRows()) {
            rGlobalSpaceDerivatives.resize(shape_function_container.NumberOfShapeFunctionRows());
        }

        for (IndexType shape_function_row_i = 0;
            shape_function_row_i < shape_function_container.NumberOfShapeFunctionRows();
            ++shape_function_row_i) {
            for (IndexType u = 0; u <= PolynomialDegreeU(); ++u) {
                for (IndexType v = 0; v <= PolynomialDegreeV(); ++v) {
                    for (IndexType w = 0; w <= PolynomialDegreeW(); ++w) {
                        IndexType cp_index_u = shape_function_container.GetFirstNonzeroControlPointU() + u;
                        IndexType cp_index_v = shape_function_container.GetFirstNonzeroControlPointV() + v;
                        IndexType cp_index_w = shape_function_container.GetFirstNonzeroControlPointW() + w;

                        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            NumberOfControlPointsU(), NumberOfControlPointsV(), NumberOfControlPointsW(),
                                cp_index_u, cp_index_v, cp_index_w);

                        if (u == 0 && v==0 && w==0)
                            rGlobalSpaceDerivatives[shape_function_row_i] =
                            (*this)[index] * shape_function_container(u, v, w, shape_function_row_i);
                        else
                            rGlobalSpaceDerivatives[shape_function_row_i] +=
                            (*this)[index] * shape_function_container(u, v, w, shape_function_row_i);
                    }
                }
            }
        }
    }

    ///@}
    ///@name Shape Function
    ///@{

    /**
     * @brief Computes the shape function values in the local parameter space.
     * @details The values are only computed for nonzero control points.
     * @param rCoordinates Coordinates to be evaluated.
     * @return Vector that holds the values. Vector(CP-Index).
     **/
    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        NurbsVolumeShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, mPolynomialDegreeW, 0);

        // Attention: Weights are not yet implemented.
        // if (IsRational()) {
        //     shape_function_container.ComputeNurbsShapeFunctionValues(mKnotsU, mKnotsV, mKnotsW, mWeights,
        //         rCoordinates[0], rCoordinates[1], rCoordinates[2]);
        // }

        shape_function_container.ComputeBSplineShapeFunctionValues(mKnotsU, mKnotsV, mKnotsW,
            rCoordinates[0], rCoordinates[1], rCoordinates[2]);


        if (rResult.size() != shape_function_container.NumberOfNonzeroControlPoints())
            rResult.resize(shape_function_container.NumberOfNonzeroControlPoints());

        for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); ++i) {
            rResult[i] = shape_function_container(i, 0);
        }

        return rResult;
    }

    /**
     * @brief Computes the first derivatives in the local parameter space
     * @details The derivatives are only computed for nonzero control points.
     * @param rCoordinates Coordinates to be evaluated.
     * @return Matrix that holds the derivatives. Matrix(CP-Index,SpaceDirection).
     **/
    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        NurbsVolumeShapeFunction shape_function_container(mPolynomialDegreeU, mPolynomialDegreeV, mPolynomialDegreeW, 1);

        // Attention: Weights are not yet implemented.
        // if (IsRational()) {
        //     shape_function_container.ComputeNurbsShapeFunctionValues(mKnotsU, mKnotsV, mKnotsW, mWeights,
        //         rCoordinates[0], rCoordinates[1], rCoordinates[2]);
        // }

        shape_function_container.ComputeBSplineShapeFunctionValues(mKnotsU, mKnotsV, mKnotsW,
            rCoordinates[0], rCoordinates[1], rCoordinates[2]);


        if (rResult.size1() != shape_function_container.NumberOfNonzeroControlPoints()
            && rResult.size2() != 3)
            rResult.resize(shape_function_container.NumberOfNonzeroControlPoints(), 3);

        for (IndexType i = 0; i < shape_function_container.NumberOfNonzeroControlPoints(); ++i) {
            rResult(i, 0) = shape_function_container(i, 1);
            rResult(i, 1) = shape_function_container(i, 2);
            rResult(i, 2) = shape_function_container(i, 3);
        }

        return rResult;
    }

    ///@}
    ///@name Geometry Family
    ///@{

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Nurbs;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Nurbs_Volume;
    }

    ///@}
    ///@name Information
    ///@{

    /**
     * @brief Turn back information as a string.
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     **/
    std::string Info() const override
    {
        return "3 dimensional nurbs geometry.";
    }

    /**
     * @brief Print information about this object.
     * @param rOStream Stream to print into.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "3 dimensional nurbs geometry.";
    }

    /**
     * @brief Print geometry's data into given stream.
     * @details Prints Polynomial degree and number of knots in each direction.
     * @param rOStream Stream to print into.
     * @see PrintInfo()
     * @see Info()
     **/
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "PolynomialDegreeU: " << mPolynomialDegreeU << "." << std::endl;
        rOStream << "PolynomialDegreeV: " << mPolynomialDegreeV << "." << std::endl;
        rOStream << "PolynomialDegreeW: " << mPolynomialDegreeW << "." << std::endl;
        rOStream << "Number of Knots in u-direction: " << mKnotsU.size() << "." << std::endl;
        rOStream << "Number of Knots in v-direction: " << mKnotsV.size() << "." << std::endl;
        rOStream << "Number of Knots in w-direction: " << mKnotsW.size() << "." << std::endl;
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
    SizeType mPolynomialDegreeW;
    Vector mKnotsU;
    Vector mKnotsV;
    Vector mKnotsW;
    // Vector mWeights;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief   Checks if the knot vector is coinciding with the number of
     *          control points and the polynomial degree.
     * @details If the knot vectors have a multiplicity of p+1 in the beginning
     *          and the end, it is reduced to p.
     **/
    void CheckAndFitKnotVectors()
    {
        SizeType num_control_points = this->size();

        if (num_control_points !=
            (NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeU, mKnotsU.size())
                * NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeV, mKnotsV.size())
                  * NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeW, mKnotsW.size()))) {
            if (num_control_points ==
                (NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeU, mKnotsU.size() - 2)
                    * NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeV, mKnotsV.size() - 2)
                      * NurbsUtilities::GetNumberOfControlPoints(mPolynomialDegreeW, mKnotsW.size() - 2))) {
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

                Vector KnotsW = ZeroVector(mKnotsW.size() - 2);
                for (SizeType i = 0; i < mKnotsW.size() - 2; ++i) {
                    KnotsW[i] = mKnotsW[i + 1];
                }
                mKnotsW = KnotsW;

            } else {
                KRATOS_ERROR
                    << "Number of controls points and polynomial degrees and number of knots do not match! " << std::endl
                    << " P: " << mPolynomialDegreeU << ", Q: " << mPolynomialDegreeV << ", D: " << mPolynomialDegreeW
                    << ", number of knots u: " << mKnotsU.size() << ", number of knots v: " << mKnotsV.size() << ", number of knots w: " << mKnotsW.size()
                    << ", number of control points: " << num_control_points << std::endl
                    << "Following condition must be achieved: ControlPoints.size() = (KnotsU.size() - P + 1) * (KnotsV.size() - Q + 1) * (KnotsW.size() - D + 1)"
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
        rSerializer.save("PolynomialDegreeW", mPolynomialDegreeW);
        rSerializer.save("KnotsU", mKnotsU);
        rSerializer.save("KnotsV", mKnotsV);
        rSerializer.save("KnotsW", mKnotsW);
        // rSerializer.save("Weights", mWeights);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("PolynomialDegreeU", mPolynomialDegreeU);
        rSerializer.load("PolynomialDegreeV", mPolynomialDegreeV);
        rSerializer.load("PolynomialDegreeW", mPolynomialDegreeW);
        rSerializer.load("KnotsU", mKnotsU);
        rSerializer.load("KnotsV", mKnotsV);
        rSerializer.load("KnotsW", mKnotsW);
        // rSerializer.load("Weights", mWeights);
    }

    NurbsVolumeGeometry() : BaseType(PointsArrayType(), &msGeometryData) {};

    ///@}

}; // class NurbsVolumeGeometry

///@}
///@name Input and output
///@{

/// input stream function
template<class TPointType>
inline std::istream& operator >> ( std::istream& rIStream,
                                   NurbsVolumeGeometry<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const NurbsVolumeGeometry<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}
///@}

template<class TPointType>
const GeometryData NurbsVolumeGeometry<TPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TPointType>
const GeometryDimension NurbsVolumeGeometry<TPointType>::msGeometryDimension(3, 3);

} // namespace Kratos

#endif // KRATOS_NURBS_VOLUME_GEOMETRY_H_INCLUDED defined
