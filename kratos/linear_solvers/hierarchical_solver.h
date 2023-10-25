//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#pragma once

#include "geometries/register_kratos_components_for_geometry.h"
#include "linear_solvers/iterative_solver.h"
#include "factories/linear_solver_factory.h"
#include "utilities/sparse_matrix_multiplication_utility.h"

namespace Kratos {

namespace Internals {

// Helper class to hold the coarse problem Gauss-Seidel iteration as a solver pointer
template<
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>
>
class GaussSeidelIteration: public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GaussSeidelIteration);

    using SparseMatrixType = typename TSparseSpaceType::MatrixType;
    using VectorType = typename TSparseSpaceType::VectorType;

    explicit GaussSeidelIteration(Parameters Settings):
        LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>()
    {
        Settings.ValidateAndAssignDefaults(GetDefaultParameters());
        int smoothing_iterations = Settings["smoothing_iterations"].GetInt();
        KRATOS_ERROR_IF(smoothing_iterations < 0)
            << "smoothing_iterations must be a positive integer, got " << smoothing_iterations << std::endl;
        mSmoothingIterations = smoothing_iterations;
    }

    GaussSeidelIteration(const GaussSeidelIteration& rOther) = delete;

    ~GaussSeidelIteration() override = default;

    GaussSeidelIteration& operator=(const GaussSeidelIteration& rOther) = delete;

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        for (IndexType i = 0; i < mSmoothingIterations; ++i) {
            const auto& r_values = rA.value_data();
            const auto& r_row_indices = rA.index1_data();
            const auto& r_col_indices = rA.index2_data();
            for (IndexType row = 0; row < rX.size(); ++row) {
                IndexType row_start = r_row_indices[row];
                IndexType row_end = r_row_indices[row+1];

                double value = rB[row];
                double diag = 1.0;
                for (IndexType i = row_start; i < row_end; ++i) {
                    IndexType col = r_col_indices[i];
                    if (col == row) {
                        diag = r_values[i];
                    }
                    else {
                        value -= r_values[i] * rX[col];
                    }
                }
                rX[row] = value / diag;
            }
        }

        return true;
    }

private:

    std::size_t mSmoothingIterations = 0;

    Parameters GetDefaultParameters() const {
        return Parameters(R"({
            "solver_type": "gauss_seidel_smoothing",
            "smoothing_iterations": 5
        })");
    }
};

}

template<
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>
>
class HierarchicalSolver: public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HierarchicalSolver);

    using SparseSpace = TSparseSpaceType;
    using DenseSpace = TDenseSpaceType;
    using PreconditionerType = TPreconditionerType;
    using ReordererType = TReordererType;

    using BaseType = IterativeSolver<SparseSpace, DenseSpace, PreconditionerType, ReordererType>;

    using SparseMatrixType = typename SparseSpace::MatrixType;
    using SparseMatrixPointerType = typename SparseSpace::MatrixPointerType;
    using SizeType = typename SparseSpace::SizeType;
    using IndexType = typename SparseSpace::IndexType;

    using VectorType = typename BaseType::VectorType;
    using VectorPointerType = typename SparseSpace::VectorPointerType;

    using DofsArrayType = typename ModelPart::DofsArrayType;

    using MatrixMapType = std::map<std::pair<IndexType, IndexType>, double>;


    explicit HierarchicalSolver(Parameters Settings):
        BaseType()
    {
        KRATOS_TRY;

        Settings.ValidateAndAssignDefaults(GetDefaultParameters());
        BaseType::SetTolerance(Settings["tolerance"].GetDouble());
        BaseType::SetMaxIterationsNumber(Settings["max_iteration"].GetInt());

        auto factory = LinearSolverFactory<SparseSpace, DenseSpace>();
        mpCoarseSolver = factory.Create(Settings["coarse_solver_settings"]);

        KRATOS_ERROR_IF(mpCoarseSolver->AdditionalPhysicalDataIsNeeded())
            << "Solvers that require physical data are not supported as coarse solver" << std::endl;

        if (Settings["fine_solver_settings"]["solver_type"].GetString() == "gauss_seidel_smoothing") {
            mpFineSolver = Kratos::make_shared<
                Internals::GaussSeidelIteration<SparseSpace, DenseSpace, ReordererType>
            >(Settings["fine_solver_settings"]);
        }
        else {
            mpFineSolver = factory.Create(Settings["fine_solver_settings"]);
        }

        KRATOS_CATCH("");
    }


    HierarchicalSolver(const HierarchicalSolver& rOther) = delete;


    ~HierarchicalSolver() override = default;


    HierarchicalSolver& operator=(const HierarchicalSolver& rOther) = delete;


    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }


    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        DofsArrayType& rDofSet,
        ModelPart& rModelPart) override
    {
        if (mLinearIds.size() == 0) {
            ClassifyDofs(rModelPart);
        }

        if (mpInterpolationMatrix == nullptr || mpRestrictionMatrix == nullptr) {
            // Set up scale transfer operators
            MatrixMapType matrix_map{};
            BuildInterpolationMatrix(rModelPart, matrix_map);
            mpInterpolationMatrix = AsSparseMatrix(matrix_map, rDofSet.size(), mLinearIds.size());
            mpRestrictionMatrix = SparseSpace::CreateEmptyMatrixPointer();
            SparseMatrixMultiplicationUtility::TransposeMatrix(*mpRestrictionMatrix, *mpInterpolationMatrix);
        }

        if (mpCoarseA == nullptr) {
            mpCoarseA = SparseSpace::CreateEmptyMatrixPointer();
        }

        SparseMatrixType restricted_a;
        SparseMatrixMultiplicationUtility::MatrixMultiplication(*mpRestrictionMatrix, rA, restricted_a);
        SparseMatrixMultiplicationUtility::MatrixMultiplication(restricted_a, *mpInterpolationMatrix, *mpCoarseA);

        if (mpCoarseB == nullptr) {
            mpCoarseB = Kratos::make_shared<VectorType>(ZeroVector(mLinearIds.size()));
        }

        if (mpCoarseX == nullptr) {
            mpCoarseX = Kratos::make_shared<VectorType>(ZeroVector(mLinearIds.size()));
        }

        if (mpCoarseSolver->AdditionalPhysicalDataIsNeeded()) {
            mpCoarseSolver->ProvideAdditionalData(*mpCoarseA, *mpCoarseX, *mpCoarseB, GetCoarseDofSet(rDofSet), rModelPart);
        }
    }


    void Clear() override
    {
        mLinearIds.clear();
        mLinearRowMap.clear();

        mpCoarseA.reset();
        mpCoarseX.reset();
        mpCoarseB.reset();

        mpInterpolationMatrix.reset();
        mpRestrictionMatrix.reset();

        mCoarseDofSet.clear();
    }


    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Hierarchical Linear Solver";
        return  buffer.str();
    }


    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        IndexType iter = 0;
        this->SetIterationsNumber(iter);

        VectorType residual = VectorType(rB);
        VectorType delta = VectorType(rX.size());
        VectorType fine_rhs = VectorType(rX.size());
        VectorType aux = VectorType(rX.size());

        SparseMatrixType& r_coarse_a = *mpCoarseA;
        VectorType& r_coarse_b = *mpCoarseB;
        VectorType& r_coarse_x = *mpCoarseX;

        KRATOS_WATCH(rA.size1());
        KRATOS_WATCH(r_coarse_a.size1());

        constexpr bool init = true;
        constexpr bool no_init = false;

        this->mBNorm = SparseSpace::TwoNorm(rB);
        this->SetResidualNorm(this->mBNorm);

        while (this->IterationNeeded()) {
            this->SetIterationsNumber(++iter);

            // New approximation using fine solver/smoothing
            SparseSpace::SetToZero(delta);
            noalias(fine_rhs) = residual;
            mpFineSolver->Solve(rA, delta, fine_rhs);

            // Calculate coarse residual dl = rl - All*deltal - Alq*deltaq
            noalias(aux) = residual;
            parallel_axpy(rA, -delta, aux, no_init);
            parallel_axpy(*mpRestrictionMatrix, aux, r_coarse_b, init);

            // solve coarse problem
            mpCoarseSolver->Solve(r_coarse_a, r_coarse_x, r_coarse_b);

            // update the solution
            parallel_axpy(*mpInterpolationMatrix, r_coarse_x, delta, no_init);
            noalias(rX) += delta;

            // update the residual
            parallel_axpy(rA, -delta, residual, no_init);
            this->SetResidualNorm(SparseSpace::TwoNorm(residual));

            KRATOS_WATCH(this->GetResidualNorm());

            // check convergence
            if (this->IsConverged()) {
                return true;
            }
        }

        //return is_solved;
        return false;
    }

private:

    typename LinearSolver<SparseSpace, DenseSpace, ReordererType>::Pointer mpCoarseSolver;
    typename LinearSolver<SparseSpace, DenseSpace, ReordererType>::Pointer mpFineSolver;

    std::vector<IndexType> mLinearIds{};

    std::unordered_map<IndexType, IndexType> mLinearRowMap{};

    SparseMatrixPointerType mpCoarseA;
    VectorPointerType mpCoarseX;
    VectorPointerType mpCoarseB;

    SparseMatrixPointerType mpInterpolationMatrix = nullptr;
    SparseMatrixPointerType mpRestrictionMatrix = nullptr;

    std::vector<const Geometry<Node>*> mReferenceGeometries{};

    DofsArrayType mCoarseDofSet{};


    Parameters GetDefaultParameters() const {
        return Parameters(R"({
            "solver_type": "hierarchical_solver",
            "tolerance" : 1.0e-6,
            "max_iteration" : 200,
            "coarse_solver_settings": {},
            "fine_solver_settings": {
                "solver_type": "gauss_seidel_smoothing",
                "smoothing_iterations": 5
            }
        })");
    }


    //TODO: Is unordered_set -> sorted vector faster than sorted set insert?
    void ClassifyDofs(const ModelPart& rModelPart)
    {
        std::unordered_set<IndexType> linear_ids{};

        Element::EquationIdVectorType equation_ids;
        const auto& r_process_info = rModelPart.GetProcessInfo();

        for (const auto& r_elem: rModelPart.Elements()) {
            const Geometry<Node>& r_geom = r_elem.GetGeometry();
            r_elem.EquationIdVector(equation_ids, r_process_info);

            const SizeType block_size = equation_ids.size() / r_geom.PointsNumber();
            const IndexType high_order_offset = NumberOfLinearNodes(r_geom) * block_size;

            linear_ids.insert(
                equation_ids.begin(),
                equation_ids.begin() + high_order_offset);
        }

        for (const auto& r_cond: rModelPart.Conditions()) {
            const Geometry<Node>& r_geom = r_cond.GetGeometry();
            r_cond.EquationIdVector(equation_ids, r_process_info);

            const SizeType block_size = equation_ids.size() / r_geom.PointsNumber();
            const IndexType high_order_offset = NumberOfLinearNodes(r_geom) * block_size;

            linear_ids.insert(
                equation_ids.begin(),
                equation_ids.begin() + high_order_offset);
        }

        mLinearIds = std::vector<IndexType>(linear_ids.begin(), linear_ids.end());
        std::sort(mLinearIds.begin(), mLinearIds.end());

        mLinearRowMap.clear();
        for (std::size_t i = 0; i < mLinearIds.size(); i++) {
            mLinearRowMap.insert({mLinearIds[i], i});
        }
    }


    DofsArrayType& GetCoarseDofSet(DofsArrayType& rFineDofSet)
    {
        if (!mCoarseDofSet.size()) {
            mCoarseDofSet.reserve(mLinearIds.size());
            for (auto id: mLinearIds) {
                auto i_fine_dof = rFineDofSet.ptr_begin() + id;
                KRATOS_DEBUG_ERROR_IF((*i_fine_dof)->EquationId() != id) << "Unexpected dofset ordering" << std::endl;
                mCoarseDofSet.push_back(*i_fine_dof);
            }
        }

        return mCoarseDofSet;
    }


    SizeType NumberOfLinearNodes(const Geometry<Node>& rGeom) const
    {
        KRATOS_TRY;
        using Family = GeometryData::KratosGeometryFamily;
        switch (rGeom.GetGeometryFamily())
        {
            case Family::Kratos_Linear: return 2;
            case Family::Kratos_Triangle: return 3;
            case Family::Kratos_Quadrilateral: return 4;
            case Family::Kratos_Tetrahedra: return 4;
            case Family::Kratos_Hexahedra: return 8;
            case Family::Kratos_Prism: return 6;
            default: KRATOS_ERROR << "Unsupported geometry family " << (int)rGeom.GetGeometryFamily() << std::endl;
        }
        KRATOS_CATCH("");
    }

    const Geometry<Node>& LinearEquivalent(const Geometry<Node>& rGeom)
    {
        KRATOS_TRY;
        constexpr int offset = (int)GeometryData::KratosGeometryType::NumberOfGeometryTypes;
        std::size_t index = (3 - rGeom.WorkingSpaceDimension()) * offset + (int)rGeom.GetGeometryFamily();
        //TODO: is double checked locking necessary here? Can I initialize this in a safer way (not in parallel)?
        if (index >= mReferenceGeometries.size() || mReferenceGeometries[index] == nullptr) {
            KRATOS_CRITICAL_SECTION {
                if (index >= mReferenceGeometries.size()) mReferenceGeometries.resize(index+1);
                if (mReferenceGeometries[index] == nullptr) {
                    mReferenceGeometries[index] = &(ReferenceLinearEquivalent(rGeom));
                }
            }

        }
        return *(mReferenceGeometries[index]);
        KRATOS_CATCH("");
    }

    const Geometry<Node>& ReferenceLinearEquivalent(const Geometry<Node>& rGeometry) const
    {
        using GeoType = GeometryData::KratosGeometryType;
        switch (rGeometry.GetGeometryType())
        {
            case GeoType::Kratos_Hexahedra3D20:
            case GeoType::Kratos_Hexahedra3D27:
            case GeoType::Kratos_Hexahedra3D8:
                return KratosComponents<Geometry<Node>>::Get("Hexahedra3D8");
            case GeoType::Kratos_Prism3D15:
            case GeoType::Kratos_Prism3D6:
                return KratosComponents<Geometry<Node>>::Get("Prism3D6");
            case GeoType::Kratos_Quadrilateral2D4:
            case GeoType::Kratos_Quadrilateral2D8:
            case GeoType::Kratos_Quadrilateral2D9:
                return KratosComponents<Geometry<Node>>::Get("Quadrilateral2D4");
            case GeoType::Kratos_Quadrilateral3D4:
            case GeoType::Kratos_Quadrilateral3D8:
            case GeoType::Kratos_Quadrilateral3D9:
                return KratosComponents<Geometry<Node>>::Get("Quadrilateral3D4");
            case GeoType::Kratos_Tetrahedra3D10:
            case GeoType::Kratos_Tetrahedra3D4:
                return KratosComponents<Geometry<Node>>::Get("Tetrahedra3D4");
            case GeoType::Kratos_Triangle2D3:
            case GeoType::Kratos_Triangle2D6:
            case GeoType::Kratos_Triangle2D10:
            case GeoType::Kratos_Triangle2D15:
                return KratosComponents<Geometry<Node>>::Get("Triangle2D3");
            case GeoType::Kratos_Triangle3D3:
            case GeoType::Kratos_Triangle3D6:
                return KratosComponents<Geometry<Node>>::Get("Triangle3D3");
            case GeoType::Kratos_Line2D2:
            case GeoType::Kratos_Line2D3:
            case GeoType::Kratos_Line2D4:
            case GeoType::Kratos_Line2D5:
                return KratosComponents<Geometry<Node>>::Get("Line2D2");
            case GeoType::Kratos_Line3D2:
            case GeoType::Kratos_Line3D3:
                return KratosComponents<Geometry<Node>>::Get("Line3D2");

            default: KRATOS_ERROR << "Unsupported geometry type " << (int)rGeometry.GetGeometryType() << std::endl;
        }
        return KratosComponents<Geometry<Node>>::Get("DefaultOption");
    }


    void BuildInterpolationMatrix(const ModelPart& rModelPart, MatrixMapType& rMatrixMap)
    {
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        Element::EquationIdVectorType equation_ids;
        Matrix local_coordinates;
        Vector shape_functions;
        array_1d<double,3> coordinates(3, 0.0);

        for (const Element& r_element: rModelPart.Elements()) {
            r_element.EquationIdVector(equation_ids, r_process_info);
            const Geometry<Node>& r_geom = r_element.GetGeometry();
            r_geom.PointsLocalCoordinates(local_coordinates);

            const SizeType nodes_number = local_coordinates.size1();
            const SizeType local_dimension = local_coordinates.size2();
            const SizeType block_size = equation_ids.size() / nodes_number;

            const Geometry<Node>& r_linear_geometry = LinearEquivalent(r_geom);

            for (IndexType n = 0; n < r_linear_geometry.PointsNumber(); ++n) {
                for (IndexType b = 0; b < block_size; ++b) {
                    IndexType fine_eq_id = equation_ids[n*block_size+b];
                    rMatrixMap.insert({{fine_eq_id, mLinearRowMap[fine_eq_id]}, 1.0});
                }
            }

            for (IndexType n = r_linear_geometry.PointsNumber(); n < nodes_number; ++n) {
                for (IndexType d = 0; d < local_dimension; ++d) {
                    coordinates[d] = local_coordinates(n, d);
                }
                r_linear_geometry.ShapeFunctionsValues(shape_functions, coordinates);

                for (IndexType b = 0; b < block_size; ++b) {
                    IndexType fine_eq_id = equation_ids[n*block_size+b];
                    for (IndexType i = 0; i < shape_functions.size(); ++i) {
                        IndexType coarse_eq_id = equation_ids[i*block_size+b];
                        rMatrixMap.insert({{fine_eq_id, mLinearRowMap[coarse_eq_id]}, shape_functions[i]});
                    }
                }
            }
        }
    }


    SparseMatrixPointerType AsSparseMatrix(
        const MatrixMapType& rMatrixMap,
        SizeType FineSize,
        SizeType CoarseSize)
    {
        SparseMatrixPointerType p_matrix = Kratos::make_shared<SparseMatrixType>(
            FineSize, CoarseSize, rMatrixMap.size()); // rows, cols, non-zeros

        auto& r_values = p_matrix->value_data();
        auto& r_rows = p_matrix->index1_data();
        auto& r_cols = p_matrix->index2_data();
        IndexType current_row = 0;
        IndexType position = 0;
        r_rows[0] = 0;

        for (const auto& r_item: rMatrixMap) {
            const auto row = r_item.first.first;
            const auto col = r_item.first.second;

            while (row > current_row) {
                // pre-increment because we actually write the row END position
                r_rows[++current_row] = position;
            }

            r_cols[position] = col;
            r_values[position] = r_item.second;
            ++position;
        }

        r_rows[FineSize] = position;

        p_matrix->set_filled(r_rows.size(), r_values.size());

        return p_matrix;
    }


    void parallel_axpy(
        const SparseMatrixType& rA,
        const VectorType& rX,
        VectorType& rY,
        bool ResetY)
    {
        const auto& r_values = rA.value_data();
        const auto& r_rows = rA.index1_data();
        const auto& r_cols = rA.index2_data();

        IndexPartition<IndexType>(rY.size()).for_each([&](IndexType row) {
            double total = ResetY ? 0.0 : rY[row];
            for (IndexType i = r_rows[row]; i < r_rows[row+1]; i++) {
                IndexType col = r_cols[i];
                total += r_values[i] * rX[col];
            }
            rY[row] = total;
        });
    }

};

}