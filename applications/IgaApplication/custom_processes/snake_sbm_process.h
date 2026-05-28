//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//                   Andrea Gorgi
//

#pragma once

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "spatial_containers/bins_dynamic.h"
#include "processes/process.h"
#include "geometries/nurbs_curve_geometry.h"

#include <queue>


namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class SnakeSbmProcess
 * @brief Process class for implementing Snake-based Surrogate Boundary Method (SBM).
 * This class provides various functions to create and manipulate surrogate boundaries
 * for Iga models in Kratos.
 */
class KRATOS_API(IGA_APPLICATION) SnakeSbmProcess
    : public Process
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SnakeSbmProcess
    KRATOS_CLASS_POINTER_DEFINITION(SnakeSbmProcess);

    class SearchPointType : public Point
    {
    public:
        using BaseType = Point;
        using ObjectType = Condition;

        KRATOS_CLASS_POINTER_DEFINITION(SearchPointType);

        SearchPointType() = default;

        SearchPointType(
            const double X,
            const double Y,
            const double Z)
            : BaseType(X, Y, Z)
        {
        }

        SearchPointType(
            Condition::Pointer pCondition,
            const double X,
            const double Y,
            const double Z)
            : BaseType(X, Y, Z)
            , mpCondition(std::move(pCondition))
        {
        }

        Condition::Pointer pGetObject() const
        {
            return mpCondition;
        }

        std::size_t Id() const
        {
            KRATOS_DEBUG_ERROR_IF(mpCondition.get() == nullptr)
                << "SearchPointType does not hold a condition." << std::endl;
            return mpCondition->Id();
        }

    private:
        Condition::Pointer mpCondition = nullptr;
    };

    class KnotSpanGrid2D
    {
    public:
        using StateType = std::int8_t;

        class RowProxy
        {
        public:
            RowProxy(StateType* pData, const int NumberOfColumns)
                : mpData(pData)
                , mNumberOfColumns(NumberOfColumns)
            {
            }

            StateType& operator[](const int Column)
            {
                return mpData[Column];
            }

            StateType operator[](const int Column) const
            {
                return mpData[Column];
            }

            std::size_t size() const
            {
                return static_cast<std::size_t>(mNumberOfColumns);
            }

        private:
            StateType* mpData = nullptr;
            int mNumberOfColumns = 0;
        };

        class ConstRowProxy
        {
        public:
            ConstRowProxy(const StateType* pData, const int NumberOfColumns)
                : mpData(pData)
                , mNumberOfColumns(NumberOfColumns)
            {
            }

            StateType operator[](const int Column) const
            {
                return mpData[Column];
            }

            std::size_t size() const
            {
                return static_cast<std::size_t>(mNumberOfColumns);
            }

        private:
            const StateType* mpData = nullptr;
            int mNumberOfColumns = 0;
        };

        KnotSpanGrid2D() = default;

        KnotSpanGrid2D(
            const int NumberOfRows,
            const int NumberOfColumns,
            const StateType InitialValue = 0)
            : mNumberOfRows(NumberOfRows)
            , mNumberOfColumns(NumberOfColumns)
            , mData(static_cast<std::size_t>(NumberOfRows) * static_cast<std::size_t>(NumberOfColumns), InitialValue)
        {
        }

        StateType& operator()(const int Row, const int Column)
        {
            return mData[FlattenIndex(Row, Column)];
        }

        StateType operator()(const int Row, const int Column) const
        {
            return mData[FlattenIndex(Row, Column)];
        }

        RowProxy operator[](const int Row)
        {
            return RowProxy(mData.data() + FlattenIndex(Row, 0), mNumberOfColumns);
        }

        ConstRowProxy operator[](const int Row) const
        {
            return ConstRowProxy(mData.data() + FlattenIndex(Row, 0), mNumberOfColumns);
        }

        int Rows() const
        {
            return mNumberOfRows;
        }

        int Columns() const
        {
            return mNumberOfColumns;
        }

        std::size_t Size() const
        {
            return mData.size();
        }

        std::size_t size() const
        {
            return static_cast<std::size_t>(mNumberOfRows);
        }

        std::size_t MemoryUsageInBytes() const
        {
            return mData.size() * sizeof(StateType);
        }

    private:
        std::size_t FlattenIndex(const int Row, const int Column) const
        {
            return static_cast<std::size_t>(Row) * static_cast<std::size_t>(mNumberOfColumns)
                + static_cast<std::size_t>(Column);
        }

        int mNumberOfRows = 0;
        int mNumberOfColumns = 0;
        std::vector<StateType> mData;
    };

    class KnotSpanGrid3D
    {
    public:
        using StateType = std::int8_t;

        class RowProxy
        {
        public:
            RowProxy(StateType* pData, const int NumberOfColumns)
                : mpData(pData)
                , mNumberOfColumns(NumberOfColumns)
            {
            }

            StateType& operator[](const int Column)
            {
                return mpData[Column];
            }

            StateType operator[](const int Column) const
            {
                return mpData[Column];
            }

            std::size_t size() const
            {
                return static_cast<std::size_t>(mNumberOfColumns);
            }

        private:
            StateType* mpData = nullptr;
            int mNumberOfColumns = 0;
        };

        class ConstRowProxy
        {
        public:
            ConstRowProxy(const StateType* pData, const int NumberOfColumns)
                : mpData(pData)
                , mNumberOfColumns(NumberOfColumns)
            {
            }

            StateType operator[](const int Column) const
            {
                return mpData[Column];
            }

            std::size_t size() const
            {
                return static_cast<std::size_t>(mNumberOfColumns);
            }

        private:
            const StateType* mpData = nullptr;
            int mNumberOfColumns = 0;
        };

        class LayerProxy
        {
        public:
            LayerProxy(StateType* pData, const int NumberOfRows, const int NumberOfColumns)
                : mpData(pData)
                , mNumberOfRows(NumberOfRows)
                , mNumberOfColumns(NumberOfColumns)
            {
            }

            RowProxy operator[](const int Row)
            {
                return RowProxy(mpData + static_cast<std::size_t>(Row) * static_cast<std::size_t>(mNumberOfColumns), mNumberOfColumns);
            }

            ConstRowProxy operator[](const int Row) const
            {
                return ConstRowProxy(mpData + static_cast<std::size_t>(Row) * static_cast<std::size_t>(mNumberOfColumns), mNumberOfColumns);
            }

            std::size_t size() const
            {
                return static_cast<std::size_t>(mNumberOfRows);
            }

        private:
            StateType* mpData = nullptr;
            int mNumberOfRows = 0;
            int mNumberOfColumns = 0;
        };

        class ConstLayerProxy
        {
        public:
            ConstLayerProxy(const StateType* pData, const int NumberOfRows, const int NumberOfColumns)
                : mpData(pData)
                , mNumberOfRows(NumberOfRows)
                , mNumberOfColumns(NumberOfColumns)
            {
            }

            ConstRowProxy operator[](const int Row) const
            {
                return ConstRowProxy(mpData + static_cast<std::size_t>(Row) * static_cast<std::size_t>(mNumberOfColumns), mNumberOfColumns);
            }

            std::size_t size() const
            {
                return static_cast<std::size_t>(mNumberOfRows);
            }

        private:
            const StateType* mpData = nullptr;
            int mNumberOfRows = 0;
            int mNumberOfColumns = 0;
        };

        KnotSpanGrid3D() = default;

        KnotSpanGrid3D(
            const int NumberOfLayers,
            const int NumberOfRows,
            const int NumberOfColumns,
            const StateType InitialValue = 0)
            : mNumberOfLayers(NumberOfLayers)
            , mNumberOfRows(NumberOfRows)
            , mNumberOfColumns(NumberOfColumns)
            , mData(
                static_cast<std::size_t>(NumberOfLayers)
                    * static_cast<std::size_t>(NumberOfRows)
                    * static_cast<std::size_t>(NumberOfColumns),
                InitialValue)
        {
        }

        StateType& operator()(const int Layer, const int Row, const int Column)
        {
            return mData[FlattenIndex(Layer, Row, Column)];
        }

        StateType operator()(const int Layer, const int Row, const int Column) const
        {
            return mData[FlattenIndex(Layer, Row, Column)];
        }

        LayerProxy operator[](const int Layer)
        {
            return LayerProxy(
                mData.data() + FlattenIndex(Layer, 0, 0),
                mNumberOfRows,
                mNumberOfColumns);
        }

        ConstLayerProxy operator[](const int Layer) const
        {
            return ConstLayerProxy(
                mData.data() + FlattenIndex(Layer, 0, 0),
                mNumberOfRows,
                mNumberOfColumns);
        }

        int Layers() const
        {
            return mNumberOfLayers;
        }

        int Rows() const
        {
            return mNumberOfRows;
        }

        int Columns() const
        {
            return mNumberOfColumns;
        }

        std::size_t Size() const
        {
            return mData.size();
        }

        std::size_t size() const
        {
            return static_cast<std::size_t>(mNumberOfLayers);
        }

        std::size_t MemoryUsageInBytes() const
        {
            return mData.size() * sizeof(StateType);
        }

    private:
        std::size_t FlattenIndex(const int Layer, const int Row, const int Column) const
        {
            return (static_cast<std::size_t>(Layer) * static_cast<std::size_t>(mNumberOfRows)
                    + static_cast<std::size_t>(Row))
                    * static_cast<std::size_t>(mNumberOfColumns)
                + static_cast<std::size_t>(Column);
        }

        int mNumberOfLayers = 0;
        int mNumberOfRows = 0;
        int mNumberOfColumns = 0;
        std::vector<StateType> mData;
    };

    ///@name Life Cycle
    ///@{

    /// Constructor
    SnakeSbmProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~SnakeSbmProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        CreateTheSnakeCoordinates();
    };
    
    void ExecuteInitialize() override
    {};

    void ExecuteInitializeSolutionStep() override
    {};

    void ExecuteBeforeSolutionLoop() override
    {};

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "",
            "skin_model_part_inner_initial_name" : "SkinModelPartInnerInitial",
            "skin_model_part_outer_initial_name" : "SkinModelPartOuterInitial",
            "skin_model_part_name" : "SkinModelPart",
            "create_surr_outer_from_surr_inner": false,
            "echo_level" : 0,
            "lambda_inner" : 0.5,
            "lambda_outer" : 0.5,
            "number_of_inner_loops": 0,
            "number_initial_points_if_importing_nurbs": 5000
        })" );

        return default_parameters;
    }

protected:
    /**
    * @brief Creates the initial snake coordinates for 2D skin.
    */
    void CreateTheSnakeCoordinates(bool RemoveIslands = false);

    // Helper to generate the outer-initial skin from the surrogate inner
    void GenerateOuterInitialFromSurrogateInner();

    /**
     * @brief Create a The Snake Coordinates object
     * 
     * @tparam TIsInnerLoop 
     * @param rIgaModelPart 
     * @param rSkinModelPartInitial 
     * @param rSkinModelPart 
     * @param NumberOfInnerLoops 
     * @param Lambda 
     */
    template <bool TIsInnerLoop>
    static void CreateTheSnakeCoordinates(
        const ModelPart& rSkinModelPartInitial,
        const std::size_t NumberOfLoops,
        const double Lambda,
        IndexType EchoLevel,
        ModelPart& rIgaModelPart,
        ModelPart& rSkinModelPart,
        const int NumberInitialPointsIfImportingNurbs,
        bool RemoveIslands = false
        );


    Model* mpModel = nullptr;
    Parameters mThisParameters;
    IndexType mEchoLevel;
    double mLambdaInner;
    double mLambdaOuter;
    std::size_t mNumberOfInnerLoops;
    int mNumberInitialPointsIfImportingNurbs;
    bool mCreateSurrOuterFromSurrInner;
    ModelPart* mpIgaModelPart = nullptr; 
    ModelPart* mpSkinModelPartInnerInitial = nullptr; 
    ModelPart* mpSkinModelPartOuterInitial = nullptr; 
    ModelPart* mpSkinModelPart = nullptr; 

    using PointType = SearchPointType;
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointTypePointer>;
    using PointIterator = PointVector::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;
    using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    using DynamicBinsPointerType = DynamicBins::PointerType;
    using NurbsCurveGeometryPointerType = NurbsCurveGeometry<2, PointerVector<Node>>::Pointer;
    using CoordinatesArrayType = Geometry<Node>::CoordinatesArrayType;

    /**
     * @brief Create a The Snake Coordinates 2 D object
     * 
     */
   void CreateTheSnakeCoordinates2D(bool RemoveIslands);

    /**
     * @brief 
     * 
     * @param rPointP 
     * @param rPointQ 
     * @param rPointR 
     * @return double 
     */
    static double Orientation(
        const Node& rPointP,
        const Node& rPointQ,
        const Node& rPointR)
    {
        return (rPointQ.X() - rPointP.X()) * (rPointR.Y() - rPointP.Y()) -
               (rPointQ.Y() - rPointP.Y()) * (rPointR.X() - rPointP.X());
    }

    /**
     * @brief 
     * 
     * @param rPointP 
     * @param rPointQ 
     * @param rPointR 
     * @return true 
     * @return false 
     */
    static bool OnSegment(
        const Node& rPointP,
        const Node& rPointQ,
        const Node& rPointR)
    {
        return (std::min(rPointP.X(), rPointR.X()) <= rPointQ.X() && rPointQ.X() <= std::max(rPointP.X(), rPointR.X()) &&
                std::min(rPointP.Y(), rPointR.Y()) <= rPointQ.Y() && rPointQ.Y() <= std::max(rPointP.Y(), rPointR.Y()));
    }

    /**
     * @brief 
     * 
     * @param rPointA 
     * @param rPointB 
     * @param rPointC 
     * @param rPointD 
     * @return true 
     * @return false 
     */
    static bool SegmentsIntersect(
        const Node& rPointA,
        const Node& rPointB,
        const Node& rPointC,
        const Node& rPointD)
    {
        const double orientation_1 = Orientation(rPointA, rPointB, rPointC);
        const double orientation_2 = Orientation(rPointA, rPointB, rPointD);
        const double orientation_3 = Orientation(rPointC, rPointD, rPointA);
        const double orientation_4 = Orientation(rPointC, rPointD, rPointB);

        if (orientation_1 * orientation_2 < 0.0 && orientation_3 * orientation_4 < 0.0) {
            return true;
        }

        if (std::abs(orientation_1) < 1e-14 && OnSegment(rPointA, rPointC, rPointB)) return true;
        if (std::abs(orientation_2) < 1e-14 && OnSegment(rPointA, rPointD, rPointB)) return true;
        if (std::abs(orientation_3) < 1e-14 && OnSegment(rPointC, rPointA, rPointD)) return true;
        if (std::abs(orientation_4) < 1e-14 && OnSegment(rPointC, rPointB, rPointD)) return true;

        return false;
    }

    /**
     * @brief Retrieves or creates a node in the specified model part.
     * 
     * @param rModelPart The model part to search in.
     * @param NodeId The ID of the node to retrieve or create.
     * @param NodeI The I coordinate of the node.
     * @param NodeJ The J coordinate of the node.
     * @param rKnotVectorU The U knot vector.
     * @param rKnotVectorV The V knot vector.
     */
    static void RetrieveOrCreateNodeInModelPart(
        ModelPart& rModelPart,
        const IndexType NodeId,
        const int NodeI,
        const int NodeJ,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV);

    /**
     * @brief Retrieves or creates a 3D node in the specified model part.
     *
     * @param rModelPart The model part to search in.
     * @param NodeId The ID of the node to retrieve or create.
     * @param NodeI The I coordinate of the node.
     * @param NodeJ The J coordinate of the node.
     * @param NodeK The K coordinate of the node.
     * @param rKnotVectorU The U knot vector.
     * @param rKnotVectorV The V knot vector.
     * @param rKnotVectorW The W knot vector.
     */
    static void RetrieveOrCreateNodeInModelPart(
        ModelPart& rModelPart,
        const IndexType NodeId,
        const int NodeI,
        const int NodeJ,
        const int NodeK,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rKnotVectorW);

private:

    /**
     * @brief Performs a single step in the snake algorithm for 2D models.
     * 
     * @param IdMatrixKnotSpansAvailable ID of the matrix tracking available knot spans.
     * @param rKnotSpansUV Knot spans in UV coordinates.
     * @param rConditionCoords XY coordinates of control points.
     * @param KnotStepUV Step size in UV space.
     * @param StartingPositionUV Starting position in UV space.
     * @param rSkinModelPart The skin model part to be updated.
     * @param rKnotSpansAvailable Knot spans available for the snake.
     */
    static void SnakeStep(
        const int IdMatrixKnotSpansAvailable, 
        const std::vector<std::vector<int>>& rKnotSpansUV, 
        const std::vector<std::vector<double>>& rConditionCoordinates, 
        const array_1d<double, 2>& rKnotStepUV, 
        const array_1d<double, 2>& rStartingPositionUV,
        ModelPart& rSkinModelPart, 
        std::vector<KnotSpanGrid2D>& rKnotSpansAvailable
        );
    
    /**
     * @brief Performs a single step in the snake algorithm for 2D models with NURBS geometries.
     * 
     * @param IdMatrixKnotSpansAvailable ID of the matrix tracking available knot spans.
     * @param rKnotSpansUV  Knot spans in UV coordinates.
     * @param rConditionCoord  XY coordinates of control points.
     * @param rKnotStepUV  Step size in UV space.
     * @param rStartingPosition Starting position in UV space.
     * @param rLocalCoords Local Coords of the point in the NURBS curve geometry.
     * @param rpCurve Pointer of the NURBS curve geometry.
     * @param rSkinModelPart The skin model part to be updated.
     * @param rKnotSpansAvailable  Knot spans available for the snake.
     */
    static void SnakeStepNurbs(
            const int IdMatrixKnotSpansAvailable, 
            const std::vector<std::vector<int>>& rKnotSpansUV, 
            const std::vector<std::vector<double>>& rConditionCoord, 
            const array_1d<double, 2>& rKnotStepUV, 
            const array_1d<double, 2>& rStartingPosition,
            const std::vector<double>& rLocalCoords,
            const NurbsCurveGeometryPointerType& rpCurve,
            ModelPart& rSkinModelPart, 
            std::vector<KnotSpanGrid2D>& rKnotSpansAvailable
            );

    /**
     * @brief Checks if a point is inside the given skin boundary.
     * @param rPointToSearch The point to be checked.
     * @param rPointsBin Spatial bins to search for nearby points.
     * @param rSkinModelPartIn Input skin model part.
     * @return True if the point is inside the boundary, false otherwise.
     */
    static bool IsPointInsideSkinBoundary(
        const Point& rPointToSearch, 
        DynamicBins& rPointsBin,
        const ModelPart& rSkinModelPartIn
        );

    /**
     * @brief Marks the knot spans as available for the snake algorithm.
     * @param IdMatrix ID of the matrix being updated.
     * @param rPointsBin Spatial bins to search for nearby points.
     * @param rSkinModelPart Skin model part being analyzed.
     * @param Lambda Relaxation parameter.
     * @param rNumberKnotSpansUV Number of knot spans in UV directions.
     * @param rKnotStepUV Step size in UV space.
     * @param rStartingPositionUV Starting position in UV space.
     * @param rKnotSpansAvailable Reference to the knot spans availability matrix.
     */
    static void MarkKnotSpansAvailable(
        const int IdMatrix,
        DynamicBins& rPointsBin, 
        const ModelPart& rSkinModelPart, 
        const double Lambda, 
        const std::vector<int>& rNumberKnotSpansUV, 
        const array_1d<double, 2>& rKnotStepUV, 
        const Vector& rStartingPositionUV,
        std::vector<KnotSpanGrid2D>& rKnotSpansAvailable
        );

    /**
     * @brief Create a Surrogate Buondary From Snake Inner object
     * 
     * @param IdMatrix 
     * @param rSkinModelPartInner 
     * @param rPointsBinInner 
     * @param rNumberKnotSpans 
     * @param rKnotVectorU 
     * @param rKnotVectorV 
     * @param rStartingPositionUV
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartInner 
     */
    static void CreateSurrogateBuondaryFromSnakeInner(
        const int IdMatrix,
        const ModelPart& rSkinModelPartInner,
        DynamicBins& rPointsBinInner,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rStartingPositionUV,
        std::vector<KnotSpanGrid2D>& rKnotSpansAvailable,
        ModelPart& rSurrogateModelPartInner
        );

    /**
     * @brief Create a Surrogate Buondary From Snake Outer object
     * 
     * @param IdMatrix 
     * @param rSkinModelPartOuter 
     * @param rPointsBinOuter 
     * @param rNumberKnotSpans 
     * @param rKnotVectorU 
     * @param rKnotVectorV 
     * @param rStartingPositionUV 
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartOuter 
     */
    static void CreateSurrogateBuondaryFromSnakeOuter(
        const int IdMatrix,
        const ModelPart& rSkinModelPartOuter,
        DynamicBins& rPointsBinOuter,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rStartingPositionUV,
        std::vector<KnotSpanGrid2D>& rKnotSpansAvailable,
        ModelPart& rSurrogateModelPartOuter
        );
    
    /**
     * @brief Checks if the knot span is at the border of the parameter sapce
     * 
     * @param rKnotSpanUV 
     * @param rNumberKnotSpansUV 
     * @return true 
     * @return false 
     */
    static bool IsInside(
        const std::vector<std::vector<int>> & rKnotSpanUV,
        const std::vector<int> &rNumberKnotSpansUV
        ); 

    
    /**
     * @brief Remove the islands surrogate domain disconnected from the main one
     * 
     * @tparam TIsInnerLoop 
     * @param grid 
     */
    template <bool TIsInnerLoop>
    static void KeepLargestZeroIsland(KnotSpanGrid2D& rGrid);

    /**
     * @brief Create a The Snake Coordinates 3 D object
     * 
     */
   void CreateTheSnakeCoordinates3D();


   /**
     * @brief Create a The Snake Coordinates object
     * 
     * @tparam TIsInnerLoop 
     * @param rIgaModelPart 
     * @param rSkinModelPartInitial 
     * @param rSkinModelPart 
     * @param NumberOfInnerLoops 
     * @param Lambda 
     */
    template <bool TIsInnerLoop>
    static void CreateTheSnakeCoordinates3D(
        const ModelPart& rSkinModelPartInitial,
        const std::size_t NumberOfInnerLoops,
        const double Lambda,
        IndexType EchoLevel,
        ModelPart& rIgaModelPart,
        ModelPart& rSkinModelPart 
        );
    
    /**
     * @brief Checks if the knot span is at the border of the parameter sapce
     * 
     * @param rKnotSpanUV 
     * @param rNumberKnotSpansUV 
     * @return true 
     * @return false 
     */
    static bool IsInside3D(
        const std::vector<std::vector<int>> & rKnotSpanUV,
        const std::vector<int> &rNumberKnotSpansUV
        ); 

    /**
     * @brief Removes isolated active cells from the 3D knot-span grids.
     * @param rKnotSpansAvailable 3D knot-span availability matrices for all loops.
     */
    static void RemoveIslands3D(
        std::vector<KnotSpanGrid3D>& rKnotSpansAvailable
        );

    
    /**
     * @brief Performs a single step in the snake algorithm for 3D models.
     * 
     * @param IdMatrixKnotSpansAvailable ID of the matrix tracking available knot spans.
     * @param rKnotSpansUV Knot spans in UV coordinates.
     * @param rConditionCoords XY coordinates of control points.
     * @param KnotStepUV Step size in UV space.
     * @param StartingPositionUV Starting position in UV space.
     * @param rSkinModelPart The skin model part to be updated.
     * @param rKnotSpansAvailable Knot spans available for the snake.
     */
    static void SnakeStep3D(
        const int IdMatrixKnotSpansAvailable, 
        const std::vector<std::vector<int>>& rKnotSpansUV, 
        const std::vector<std::vector<double>>& rConditionCoords, 
        const Vector KnotStepUV, 
        const Vector StartingPositionUV,
        ModelPart& rSkinModelPart, 
        std::vector<KnotSpanGrid3D>& rKnotSpansAvailable,
        array_1d<IndexType, 3>& ordered_ids
        );

    /**
     * @brief Marks the knot spans as available for the snake algorithm.
     * @param IdMatrix ID of the matrix being updated.
     * @param rPointsBin Spatial bins to search for nearby points.
     * @param rSkinModelPart Skin model part being analyzed.
     * @param Lambda Relaxation parameter.
     * @param rNumberKnotSpansUV Number of knot spans in UV directions.
     * @param rKnotStepUV Step size in UV space.
     * @param rStartingPositionUV Starting position in UV space.
     * @param rKnotSpansAvailable Reference to the knot spans availability matrix.
     */
    static void MarkKnotSpansAvailable3D(
        const int IdMatrix,
        DynamicBins& rPointsBin, 
        const ModelPart& rSkinModelPart, 
        const double Lambda, 
        const std::vector<int>& rNumberKnotSpansUVW, 
        const array_1d<double, 3>& rKnotStepUVW, 
        const Vector& rStartingPositionUV,
        std::vector<KnotSpanGrid3D>& rKnotSpansAvailable
        );

    /**
     * @brief Checks if a point is inside the given skin boundary.
     * @param rPointToSearch The point to be checked.
     * @param rPointsBin Spatial bins to search for nearby points.
     * @param rSkinModelPartIn Input skin model part.
     * @return True if the point is inside the boundary, false otherwise.
     */
    static bool IsPointInsideSkinBoundary3D(
        const Point& rPointToSearch, 
        DynamicBins& rPointsBin,
        const ModelPart& rSkinModelPartIn
        );

    
    /**
     * @brief Create a Surrogate Buondary From Snake Inner 3 D object
     * 
     * @param IdMatrix 
     * @param rSkinModelPartInner 
     * @param rPointsBinInner 
     * @param rNumberKnotSpans 
     * @param rKnotVectorU 
     * @param rKnotVectorV 
     * @param rKnotVectorW 
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartInner 
     */
    static void CreateSurrogateBuondaryFromSnakeInner3D(
        const int IdMatrix,
        const ModelPart& rSkinModelPartInner,
        DynamicBins& rPointsBinInner,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rKnotVectorW,
        const Vector& rStartingPositionUVW,
        std::vector<KnotSpanGrid3D>& rKnotSpansAvailable,
        ModelPart& rSurrogateModelPartInner
        );

    /**
     * @brief Create a Surrogate Buondary From Snake Outer 3 D object
     * 
     * @param IdMatrix 
     * @param rSkinModelPartOuter 
     * @param rPointsBinOuter 
     * @param rNumberKnotSpans 
     * @param rKnotVectorU 
     * @param rKnotVectorV 
     * @param rKnotVectorW 
     * @param rStartingPositionUVW 
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartOuter 
     */
    static void CreateSurrogateBuondaryFromSnakeOuter3D(
        const int IdMatrix,
        const ModelPart& rSkinModelPartOuter,
        DynamicBins& rPointsBinOuter,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rKnotVectorW,
        const Vector& rStartingPositionUVW,
        std::vector<KnotSpanGrid3D>& rKnotSpansAvailable,
        ModelPart& rSurrogateModelPartOuter
        );

}; // Class SnakeSbmProcess

}  // namespace Kratos.
