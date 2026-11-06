//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Wataru Fukuda
//

#pragma once

// System includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_local_refined_surface.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"
#include "utilities/geometry_utilities/local_refined_brep_trimming_utilities.h"

namespace Kratos
{

/**
 * @class LocalRefinedBrepSurface
 * @brief Topology wrapper for a locally-refined surface (e.g. THBSurfaceGeometry).
 *
 * @tparam TContainerPointType         Container of control-point nodes.
 * @tparam TLocalRefinedSurfaceType    Concrete locally-refined surface
 *                                     (THBSurfaceGeometry, …).
 * @tparam TShiftedBoundary            Shifted-boundary flag (default false).
 * @tparam TContainerPointEmbeddedType Container for embedded (trimming) nodes.
 */
template<
    class TContainerPointType,
    class TLocalRefinedSurfaceType,
    bool  TShiftedBoundary = false,
    class TContainerPointEmbeddedType = TContainerPointType>
class LocalRefinedBrepSurface
    : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LocalRefinedBrepSurface);

    typedef typename TContainerPointType::value_type PointType;

    typedef Geometry<PointType> BaseType;
    typedef Geometry<PointType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointer;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef TLocalRefinedSurfaceType LocalRefinedSurfaceType;

    typedef LocalRefinedBrepTrimmingUtilities<TShiftedBoundary>
        BrepTrimmingUtilitiesType;

    typedef BrepCurveOnLocalRefinedSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>
        BrepCurveOnLocalRefinedSurfaceType;
    typedef DenseVector<typename BrepCurveOnLocalRefinedSurfaceType::Pointer>
        BrepCurveOnLocalRefinedSurfaceLoopType;
    typedef DenseVector<DenseVector<typename BrepCurveOnLocalRefinedSurfaceType::Pointer>>
        BrepCurveOnLocalRefinedSurfaceLoopArrayType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for untrimmed locally-refined patch.
    LocalRefinedBrepSurface(
        typename LocalRefinedSurfaceType::Pointer pLocalRefinedSurface)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpLocalRefinedSurface(pLocalRefinedSurface)
        , mIsTrimmed(false)
    {
    }

    /// Constructor for trimmed locally-refined patch.
    LocalRefinedBrepSurface(
        typename LocalRefinedSurfaceType::Pointer pLocalRefinedSurface,
        const BrepCurveOnLocalRefinedSurfaceLoopArrayType& rOuterLoopArray,
        const BrepCurveOnLocalRefinedSurfaceLoopArrayType& rInnerLoopArray,
        bool IsTrimmed = true)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpLocalRefinedSurface(pLocalRefinedSurface)
        , mOuterLoopArray(rOuterLoopArray)
        , mInnerLoopArray(rInnerLoopArray)
        , mIsTrimmed(IsTrimmed)
    {
    }

    /// Copy constructor.
    LocalRefinedBrepSurface(const LocalRefinedBrepSurface& rOther)
        : BaseType(rOther)
        , mpLocalRefinedSurface(rOther.mpLocalRefinedSurface)
        , mOuterLoopArray(rOther.mOuterLoopArray)
        , mInnerLoopArray(rOther.mInnerLoopArray)
        , mIsTrimmed(rOther.mIsTrimmed)
    {
    }

    explicit LocalRefinedBrepSurface(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Destructor.
    ~LocalRefinedBrepSurface() override = default;

    ///@}
    ///@name Operators
    ///@{

    LocalRefinedBrepSurface& operator=(const LocalRefinedBrepSurface& rOther)
    {
        BaseType::operator=(rOther);
        mpLocalRefinedSurface = rOther.mpLocalRefinedSurface;
        mOuterLoopArray = rOther.mOuterLoopArray;
        mInnerLoopArray = rOther.mInnerLoopArray;
        mIsTrimmed = rOther.mIsTrimmed;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(PointsArrayType const& ThisPoints) const override
    {
        return typename BaseType::Pointer(new LocalRefinedBrepSurface(ThisPoints));
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    GeometryPointer pGetGeometryPart(const IndexType Index) override
    {
        const auto& const_this = *this;
        return std::const_pointer_cast<GeometryType>(const_this.pGetGeometryPart(Index));
    }

    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpLocalRefinedSurface;

        for (IndexType i = 0; i < mOuterLoopArray.size(); ++i)
            for (IndexType j = 0; j < mOuterLoopArray[i].size(); ++j)
                if (mOuterLoopArray[i][j]->Id() == Index)
                    return mOuterLoopArray[i][j];

        for (IndexType i = 0; i < mInnerLoopArray.size(); ++i)
            for (IndexType j = 0; j < mInnerLoopArray[i].size(); ++j)
                if (mInnerLoopArray[i][j]->Id() == Index)
                    return mInnerLoopArray[i][j];

        KRATOS_ERROR << "Index " << Index
            << " not existing in LocalRefinedBrepSurface: " << this->Id() << std::endl;
    }

    bool HasGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return true;

        for (IndexType i = 0; i < mOuterLoopArray.size(); ++i)
            for (IndexType j = 0; j < mOuterLoopArray[i].size(); ++j)
                if (mOuterLoopArray[i][j]->Id() == Index)
                    return true;

        for (IndexType i = 0; i < mInnerLoopArray.size(); ++i)
            for (IndexType j = 0; j < mInnerLoopArray[i].size(); ++j)
                if (mInnerLoopArray[i][j]->Id() == Index)
                    return true;

        return false;
    }

    typename LocalRefinedSurfaceType::Pointer pGetLocalRefinedSurface() const
    {
        return mpLocalRefinedSurface;
    }

    const BrepCurveOnLocalRefinedSurfaceLoopArrayType& GetOuterLoops() const { return mOuterLoopArray; }
    const BrepCurveOnLocalRefinedSurfaceLoopArrayType& GetInnerLoops() const { return mInnerLoopArray; }

    ///@}
    ///@name Mathematical Informations
    ///@{

    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpLocalRefinedSurface->PolynomialDegree(LocalDirectionIndex);
    }

    SizeType PointsNumberInDirection(IndexType DirectionIndex) const override
    {
        return mpLocalRefinedSurface->PointsNumberInDirection(DirectionIndex);
    }

    ///@}
    ///@name Information
    ///@{

    bool IsTrimmed() const { return mIsTrimmed; }

    ///@}
    ///@name Geometrical Operations
    ///@{

    Point Center() const override
    {
        return mpLocalRefinedSurface->Center();
    }

    int ProjectionPointGlobalToLocalSpace(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const override
    {
        return mpLocalRefinedSurface->ProjectionPointGlobalToLocalSpace(
            rPointGlobalCoordinates, rProjectedPointLocalCoordinates, Tolerance);
    }

    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates) const override
    {
        mpLocalRefinedSurface->GlobalCoordinates(rResult, rLocalCoordinates);
        return rResult;
    }

    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput) const override
    {
        if (rVariable == CHARACTERISTIC_GEOMETRY_LENGTH)
            mpLocalRefinedSurface->Calculate(rVariable, rOutput);
    }

    ///@}
    ///@name Integration Info
    ///@{

    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        return mpLocalRefinedSurface->GetDefaultIntegrationInfo();
    }

    ///@}
    ///@name Integration Points
    ///@{

    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        if (!mIsTrimmed)
        {
            mpLocalRefinedSurface->CreateIntegrationPoints(
                rIntegrationPoints, rIntegrationInfo);
        }
        else
        {
            std::vector<double> spans_u, spans_v;
            mpLocalRefinedSurface->SpansLocalSpace(spans_u, 0);
            mpLocalRefinedSurface->SpansLocalSpace(spans_v, 1);

            BrepTrimmingUtilitiesType::CreateBrepSurfaceTrimmingIntegrationPoints(
                rIntegrationPoints,
                mOuterLoopArray, mInnerLoopArray,
                spans_u, spans_v,
                rIntegrationInfo);
        }
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        mpLocalRefinedSurface->CreateQuadraturePointGeometries(
            rResultGeometries,
            NumberOfShapeFunctionDerivatives,
            rIntegrationPoints,
            rIntegrationInfo);

        for (IndexType i = 0; i < rResultGeometries.size(); ++i)
            rResultGeometries(i)->SetGeometryParent(this);
    }

    ///@}
    ///@name Shape Functions
    ///@{

    Vector& ShapeFunctionsValues(
        Vector& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpLocalRefinedSurface->ShapeFunctionsValues(rResult, rCoordinates);
        return rResult;
    }

    void ShapeFunctionsValuesAndCPIndices(
        const CoordinatesArrayType& rCoordinates,
        std::vector<IndexType>& rControlPointIndices,
        Vector& rShapeFunctionsValues,
        const IndexType DerivativeOrder = 0,
        DenseVector<Matrix>* pShapeFunctionDerivatives = nullptr) const
    {
        mpLocalRefinedSurface->ShapeFunctionsValuesAndCPIndices(
            rCoordinates, rControlPointIndices, rShapeFunctionsValues,
            DerivativeOrder, pShapeFunctionDerivatives);
    }

    void ShapeFunctionsValuesAndCPNodes(
        const CoordinatesArrayType& rCoordinates,
        PointsArrayType& rNonzeroControlPoints,
        Vector& rShapeFunctionsValues,
        const IndexType DerivativeOrder = 0,
        DenseVector<Matrix>* pShapeFunctionDerivatives = nullptr) const override
    {
        mpLocalRefinedSurface->ShapeFunctionsValuesAndCPNodes(
            rCoordinates, rNonzeroControlPoints, rShapeFunctionsValues,
            DerivativeOrder, pShapeFunctionDerivatives);
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpLocalRefinedSurface->ShapeFunctionsLocalGradients(rResult, rCoordinates);
        return rResult;
    }

    ///@}
    ///@name Geometry Family
    ///@{

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Brep;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Local_Refined_Brep_Surface;
    }

    ///@}
    ///@name Knot Span Accessors
    ///@{

    void SpansLocalSpace(std::vector<double>& rSpans, IndexType LocalDirectionIndex = 0) const override
    {
        mpLocalRefinedSurface->SpansLocalSpace(rSpans, LocalDirectionIndex);
    }

    void SpansLocalSpace(
        std::vector<double>& rSpans,
        IndexType DirectionIndex,
        double PositionOtherDirection) const override
    {
        mpLocalRefinedSurface->SpansLocalSpace(rSpans, DirectionIndex, PositionOtherDirection);
    }

    std::vector<double> KnotsU() const
    {
        std::vector<double> spans;
        mpLocalRefinedSurface->SpansLocalSpace(spans, 0);
        return spans;
    }

    std::vector<double> KnotsV() const
    {
        std::vector<double> spans;
        mpLocalRefinedSurface->SpansLocalSpace(spans, 1);
        return spans;
    }

    std::vector<std::array<double,4>> GetActiveCells() const
    {
        return mpLocalRefinedSurface->GetActiveCells();
    }

    bool ComputeSpanTriangulationLocalSpace(
        const double u0,
        const double u1,
        const double v0,
        const double v1,
        std::vector<Matrix>& rTrianglesLocalSpace) const
    {
        bool is_trimmed = false;
        BrepTrimmingUtilitiesType::ComputeSpanTriangulation(
            mOuterLoopArray,
            mInnerLoopArray,
            u0, u1, v0, v1,
            is_trimmed,
            rTrianglesLocalSpace);
        return is_trimmed;
    }

    ///@}
    ///@name Information
    ///@{

    std::string Info() const override { return "LocalRefinedBrepSurface"; }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LocalRefinedBrepSurface";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        rOStream << std::endl << "    LocalRefinedBrepSurface" << std::endl;
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;
    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Member Variables
    ///@{

    typename LocalRefinedSurfaceType::Pointer mpLocalRefinedSurface;

    BrepCurveOnLocalRefinedSurfaceLoopArrayType mOuterLoopArray;
    BrepCurveOnLocalRefinedSurfaceLoopArrayType mInnerLoopArray;

    bool mIsTrimmed = false;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("LocalRefinedSurface", mpLocalRefinedSurface);
        rSerializer.save("OuterLoopArray", mOuterLoopArray);
        rSerializer.save("InnerLoopArray", mInnerLoopArray);
        rSerializer.save("IsTrimmed", mIsTrimmed);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("LocalRefinedSurface", mpLocalRefinedSurface);
        rSerializer.load("OuterLoopArray", mOuterLoopArray);
        rSerializer.load("InnerLoopArray", mInnerLoopArray);
        rSerializer.load("IsTrimmed", mIsTrimmed);
    }

    LocalRefinedBrepSurface()
        : BaseType(PointsArrayType(), &msGeometryData)
    {}

    ///@}
};

///@name Static Type Declarations
///@{

template<class TContainerPointType, class TLocalRefinedSurfaceType, bool TShiftedBoundary, class TContainerPointEmbeddedType>
const GeometryData LocalRefinedBrepSurface<TContainerPointType, TLocalRefinedSurfaceType, TShiftedBoundary, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, class TLocalRefinedSurfaceType, bool TShiftedBoundary, class TContainerPointEmbeddedType>
const GeometryDimension LocalRefinedBrepSurface<TContainerPointType, TLocalRefinedSurfaceType, TShiftedBoundary, TContainerPointEmbeddedType>::msGeometryDimension(3, 2);

///@}

} // namespace Kratos
