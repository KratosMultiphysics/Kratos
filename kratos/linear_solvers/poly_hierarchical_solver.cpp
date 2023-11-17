//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "poly_hierarchical_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/profiler.h"
#include "factories/linear_solver_factory.h"
#include "utilities/proxies.h"

// STL includes
#include <type_traits> // remove_reference_t, remove_cv_t
#include <sstream> // stringstream
#include <optional> // optional


namespace Kratos {


namespace detail {


class ParametersWithDefaults
{
public:
    ParametersWithDefaults(Parameters Settings,
                           Parameters Defaults)
        : mSettings(Settings),
          mDefaults(Defaults)
    {}

    ParametersWithDefaults operator[](const std::string& rEntryName)
    {
        if (mSettings.Has(rEntryName)) {
            if (mDefaults.Has(rEntryName)) {
                return ParametersWithDefaults(mSettings[rEntryName], mDefaults[rEntryName]);
            } else {
                KRATOS_ERROR << "Input parameters has entry '" << rEntryName
                             << "' but the defaults do not. Input parameters:\n"
                             << mSettings << "\nDefault parameters:\n"
                             << mDefaults;
            }
        } else {
            if (mDefaults.Has(rEntryName)) {
                Parameters defaults = mDefaults[rEntryName];
                mSettings.AddValue(rEntryName, defaults);
                return ParametersWithDefaults(mSettings[rEntryName], defaults);
            } else {
                KRATOS_ERROR << "Neither input parameters, nor defaults have an entry for "
                             << "'" << rEntryName << "'. Input parameters:\n"
                             << mSettings << "\nDefault parameters:\n"
                             << mDefaults;
            } // else (mDefaults.Has(rEntryName))
        } // else (mSettings.Has(rEntryName))
    }

    template <class T>
    auto Get()
    {
        using TValue = std::remove_cv_t<std::remove_reference_t<T>>;

        if (!mSettings.IsNull()) {
            return mSettings.Get<TValue>();
        }

        if (!mDefaults.IsNull()) {
            auto output = mDefaults.Get<TValue>();
            mSettings.Set<TValue>(output);
            return output;
        }

        KRATOS_ERROR << "Both the input parameters and the defaults are empty.\n";
        return TValue();
    }

    Parameters GetInput()
    {
        return mSettings;
    }

    Parameters GetDefaults()
    {
        return mDefaults;
    }


private:
    Parameters mSettings;

    Parameters mDefaults;
};


/** @brief Generate four submatrices from a square matrix.
 *
 *  @details Let @f(A@f) be the square input matrix, and
 *           @f(A_{ll}@f), @f(A_{lq}@f), @f(A_{ql}@f), @f(A_{qq}@f) the
 *           output matrices in this exact order. Then this function resizes
 *           and populates the output matrices such that
 *           @f[
 *              A = \begin{bmatrix}
 *                  A_{ll} & A_{lq} \\
 *                  A_{ql} & A_{qq}
 *              \end{bmatrix}
 *           @f]
 *
 *          @f(A_{ll}@f) is defined by a boolean mask, represented by a @c char
 *          array consisting of 1s at row indices where @f(A_{ll}@f) has
 *          components, and 0s everywhere else.
 *
 *  @throws if the mask @a rMask contains items other than 0 or 1.
 */
void MakeSubblocks(const TUblasSparseSpace<double>::MatrixType& rRootMatrix,
                   const std::vector<bool>& rMask,
                   const std::array<TUblasSparseSpace<double>::MatrixType*,4>& rOutput)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF_NOT(rMask.size() == rRootMatrix.size1())
        << "Mask size mismatch: mask of size " << rMask.size()
        << " provided for matrix of size " << rRootMatrix.size1() << "x" << rRootMatrix.size2();

    KRATOS_ERROR_IF_NOT(rRootMatrix.size1() == rRootMatrix.size2())
        << "Expecting a square matrix, but got " << rRootMatrix.size1() << "x" << rRootMatrix.size2();

    KRATOS_ERROR_IF_NOT(std::count(rOutput.begin(), rOutput.end(), nullptr) == 0)
        << "Output pointers must point to existing matrices";

    const std::size_t system_size = rMask.size();

    // Compute the index of each component local to the domain it belongs to.
    std::vector<std::size_t> block_local_indices(system_size);
    std::array<std::size_t,2> dof_count {0ul, 0ul}; // {masked_count, unmasked_count}
    for (std::size_t i_mask=0; i_mask<system_size; ++i_mask) {
        block_local_indices[i_mask] = dof_count[rMask[i_mask]]++;
    }

    // Resize submatrices.
    KRATOS_TRY
    rOutput[0]->resize(dof_count[0], dof_count[0], false);
    rOutput[1]->resize(dof_count[0], dof_count[1], false);
    rOutput[2]->resize(dof_count[1], dof_count[0], false);
    rOutput[3]->resize(dof_count[1], dof_count[1], false);
    KRATOS_CATCH("")

    // Fill submatrices
    // Note: the matrices are filled in proper row-major order,
    //       so each insertion should be O(1) (unless a reallocation
    //       is triggered).
    KRATOS_TRY
    for (auto it_row=rRootMatrix.begin1(); it_row!=rRootMatrix.end1(); ++it_row) {
        const std::size_t i_row = it_row.index1();
        const char row_mask = rMask[i_row];
        for (auto it_entry=it_row.begin(); it_entry!=it_row.end(); ++it_entry) {
            const std::size_t i_column = it_entry.index2();
            const char column_mask = rMask[i_column];
            const char i_block = (~(column_mask | (row_mask << 1))) & 0b11; // <== flattened index of the block
            rOutput[i_block]->insert_element(block_local_indices[i_row],
                                             block_local_indices[i_column],
                                             *it_entry);
        }
    }
    KRATOS_CATCH("")
}


} // namespace detail



template <class TSparseSpace,
          class TDenseSpace,
          class TReorderer>
struct PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    double mTolerance;

    int mVerbosity;

    std::size_t mMaxIterations;

    typename LinearSolver<TSparseSpace,TDenseSpace/*reorderer is omitted on purpose*/>::Pointer mpFineSolver;

    typename TSparseSpace::MatrixType mContractionOperator;

    struct CoarseSystem {
        typename TSparseSpace::MatrixType mA;

        Vector mX;

        Vector mB;

        typename LinearSolver<TSparseSpace,TDenseSpace/*reorderer is omitted on purpose*/>::Pointer mpSolver;
    } mCoarseSystem;
}; // struct PolyHierarchicalSolver::Impl



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PolyHierarchicalSolver(Parameters parameters)
    : mpImpl(new Impl)
{
    KRATOS_TRY
    Parameters default_parameters = this->GetDefaultParameters();
    parameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(parameters["solver_type"].GetString() == "poly_hierarchical")
        << "Requested a(n) '" << parameters["solver_type"].GetString() << "' solver,"
        << " but constructing an PolyHierarchicalSolver";

    detail::ParametersWithDefaults safe_parameters(parameters, default_parameters);
    mpImpl->mTolerance = safe_parameters["tolerance"].Get<double>();
    mpImpl->mVerbosity = safe_parameters["verbosity"].Get<int>();
    mpImpl->mMaxIterations = safe_parameters["max_iterations"].Get<int>();

    // Construct subsolvers
    KRATOS_TRY
    LinearSolverFactory<TSparseSpace,TDenseSpace> solver_factory;
    mpImpl->mCoarseSystem.mpSolver = solver_factory.Create(safe_parameters["coarse_settings"].GetInput());
    mpImpl->mpFineSolver = solver_factory.Create(safe_parameters["fine_settings"].GetInput());
    KRATOS_CATCH("")

    KRATOS_CATCH("")
}



// Necessary for PIMPL
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::~PolyHierarchicalSolver()
{
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                        Vector& rX,
                                                                        Vector& rB)
{
    //KRATOS_ERROR_IF_NOT(mpImpl->mCoarseMask.has_value())
    //    << "PolyHierarchicalSolver::Solve called before PolyHierarchicalSolver::ProvideAdditionalData";
    //auto& r_coarse_mask = mpImpl->mCoarseMask.value();
    //KRATOS_ERROR_IF_NOT(r_coarse_mask.size() == TSparseSpace::Size1(rA))
    //    << "DoF mask size " << r_coarse_mask.size()
    //    << " does not match system size " << TSparseSpace::Size1(rA);

    KRATOS_TRY

    // Dump the generated matrices to disk f requested.
    if (4 <= mpImpl->mVerbosity) {
        KRATOS_ERROR << "verbosity >=4 prints system matrices and terminates\n";
    }

    // Declarations
    const std::size_t fine_size = TSparseSpace::Size(rX);
    const std::size_t coarse_size = TSparseSpace::Size(mpImpl->mCoarseSystem.mX);
    Vector fine_delta(fine_size);
    Vector fine_residual(fine_size);
    Vector fine_tmp(fine_size);
    Vector coarse_delta(coarse_size);
    Vector coarse_residual(coarse_size);
    Vector coarse_tmp(coarse_size);

    const double rhs_norm = TSparseSpace::TwoNorm(rB);

    // Initializations
    // r = b - A * x
    TSparseSpace::Mult(rA, rX, fine_residual);
    fine_residual = rB - fine_residual;

    double residual_norm = 1.0;
    const std::size_t max_iterations = mpImpl->mMaxIterations;
    for (std::size_t i_iteration=0ul; i_iteration<max_iterations; ++i_iteration) {
        // Perform a few smoothing passes on the full system.
        KRATOS_TRY
        mpImpl->mpFineSolver->Solve(rA, fine_delta, fine_residual);
        KRATOS_CATCH("")

        // Contract the smoothed delta and residual to the coarse level.
        TSparseSpace::Mult(mpImpl->mContractionOperator, fine_delta, coarse_delta);
        TSparseSpace::Mult(mpImpl->mContractionOperator, fine_residual, coarse_residual);

        // Compute the restricted residual:
        // r_l -= A_ll * x_l + A_lh * x_h
        TSparseSpace::Mult(mpImpl->mCoarseSystem.mA, coarse_delta, coarse_tmp);
        coarse_residual = coarse_residual - coarse_tmp;

        // Solve coarse system.
        KRATOS_TRY
        mpImpl->mCoarseSystem.mpSolver->Solve(mpImpl->mCoarseSystem.mA, coarse_delta, coarse_residual);
        KRATOS_CATCH("")

        // Update solution
        noalias(fine_delta) += prod(trans(mpImpl->mContractionOperator), coarse_delta);
        rX += fine_delta;

        // Compute the fine_residual
        TSparseSpace::Mult(rA, rX, fine_tmp);
        noalias(fine_residual) = rB - fine_tmp;

        // Update status
        residual_norm = TSparseSpace::TwoNorm(fine_residual) / rhs_norm;
        KRATOS_INFO_IF("PolyHierarchicalSolver", 1 <= mpImpl->mVerbosity)
            << "iteration " << i_iteration
            << " residual " << residual_norm << "\n";

        // Early exit if converged
        if (residual_norm < mpImpl->mTolerance) {
            break;
        }
    }

    return residual_norm < mpImpl->mTolerance;
    KRATOS_CATCH("")
}


template <class TOutputIterator>
void GetLowerOrderNodeIndices(GeometryData::KratosGeometryFamily GeometryFamily,
                              TOutputIterator itOutput)
{
    using Family = GeometryData::KratosGeometryFamily;
    switch(GeometryFamily) {
        case Family::Kratos_Point: {
            *itOutput++ = 0ul;
            break;
        }
        case Family::Kratos_Linear: {
            *itOutput++ = 0ul;
            *itOutput   = 1ul;
            break;
        }
        case Family::Kratos_Triangle: {
            *itOutput++ = 0ul;
            *itOutput++ = 1ul;
            *itOutput   = 2ul;
            break;
        }
        case Family::Kratos_Quadrilateral: {
            *itOutput++ = 0ul;
            *itOutput++ = 1ul;
            *itOutput++ = 2ul;
            *itOutput   = 3ul;
            break;
        }
        case Family::Kratos_Tetrahedra: {
            *itOutput++ = 0ul;
            *itOutput++ = 1ul;
            *itOutput++ = 2ul;
            *itOutput   = 3ul;
            break;
        }
        case Family::Kratos_Hexahedra: {
            *itOutput++ = 0ul;
            *itOutput++ = 1ul;
            *itOutput++ = 2ul;
            *itOutput++ = 3ul;
            *itOutput++ = 4ul;
            *itOutput   = 5ul;
            break;
        }
        case Family::Kratos_Pyramid: {
            *itOutput++ = 0ul;
            *itOutput++ = 1ul;
            *itOutput++ = 2ul;
            *itOutput++ = 3ul;
            *itOutput   = 4ul;
            break;
        }
        default: {
            KRATOS_ERROR << "unsupported geometry family";
        }
    } // switch (GeometryFamily)
}


const Geometry<Node> GetLowerOrderGeometry(const Geometry<Node>& rGeometry)
{
    using GeometryEnum = GeometryData::KratosGeometryType;
    using GeometryRegistry = KratosComponents<Geometry<Node>>;
    const GeometryEnum geometry_type = rGeometry.GetGeometryType();
    switch (geometry_type) {
        case GeometryEnum::Kratos_Hexahedra3D8:     return GeometryRegistry::Get("Hexahedra3D8");
        case GeometryEnum::Kratos_Hexahedra3D20:    return GeometryRegistry::Get("Hexahedra3D8");
        case GeometryEnum::Kratos_Hexahedra3D27:    return GeometryRegistry::Get("Hexahedra3D8");
        case GeometryEnum::Kratos_Line2D2:          return GeometryRegistry::Get("Line2D2");
        case GeometryEnum::Kratos_Line2D3:          return GeometryRegistry::Get("Line2D2");
        case GeometryEnum::Kratos_Line2D4:          return GeometryRegistry::Get("Line2D2");
        case GeometryEnum::Kratos_Line2D5:          return GeometryRegistry::Get("Line2D2");
        case GeometryEnum::Kratos_Line3D2:          return GeometryRegistry::Get("Line3D2");
        case GeometryEnum::Kratos_Line3D3:          return GeometryRegistry::Get("Line3D2");
        case GeometryEnum::Kratos_Point2D:          return GeometryRegistry::Get("Point2D");
        case GeometryEnum::Kratos_Point3D:          return GeometryRegistry::Get("Point3D");
        case GeometryEnum::Kratos_Prism3D6:         return GeometryRegistry::Get("Prism3D6");
        case GeometryEnum::Kratos_Prism3D15:        return GeometryRegistry::Get("Prism3D6");
        case GeometryEnum::Kratos_Pyramid3D5:       return GeometryRegistry::Get("Pyramid3D5");
        case GeometryEnum::Kratos_Pyramid3D13:      return GeometryRegistry::Get("Pyramid3D5");
        case GeometryEnum::Kratos_Quadrilateral2D4: return GeometryRegistry::Get("Quadrilateral2D4");
        case GeometryEnum::Kratos_Quadrilateral2D8: return GeometryRegistry::Get("Quadrilateral2D4");
        case GeometryEnum::Kratos_Quadrilateral2D9: return GeometryRegistry::Get("Quadrilateral2D4");
        case GeometryEnum::Kratos_Quadrilateral3D4: return GeometryRegistry::Get("Quadrilateral3D4");
        case GeometryEnum::Kratos_Quadrilateral3D8: return GeometryRegistry::Get("Quadrilateral3D4");
        case GeometryEnum::Kratos_Quadrilateral3D9: return GeometryRegistry::Get("Quadrilateral3D4");
        case GeometryEnum::Kratos_Tetrahedra3D4:    return GeometryRegistry::Get("Tetrahedra3D4");
        case GeometryEnum::Kratos_Tetrahedra3D10:   return GeometryRegistry::Get("Tetrahedra3D4");
        case GeometryEnum::Kratos_Triangle2D3:      return GeometryRegistry::Get("Triangle2D3");
        case GeometryEnum::Kratos_Triangle2D6:      return GeometryRegistry::Get("Triangle2D3");
        case GeometryEnum::Kratos_Triangle2D10:     return GeometryRegistry::Get("Triangle2D3");
        case GeometryEnum::Kratos_Triangle2D15:     return GeometryRegistry::Get("Triangle2D3");
        case GeometryEnum::Kratos_Triangle3D3:      return GeometryRegistry::Get("Triangle3D3");
        case GeometryEnum::Kratos_Triangle3D6:      return GeometryRegistry::Get("Triangle3D3");
        default: {
            KRATOS_ERROR << "unsupported geometry type (" << (int)geometry_type << ")\n";
        }
    } // switch (rGeometry.GetGeometryType())
}


template <Globals::DataLocation TEntityLocation, class TSparseSpace>
void GetLowerOrderDofs(typename TSparseSpace::MatrixType& rA,
                       const ModelPart& rModelPart,
                       std::vector<bool>& rCoarseMask)
{
    KRATOS_TRY

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    using ThreadLocalStorage = std::pair<
        Element::DofsVectorType,    // <== all DoFs of the element/condition
        std::vector<std::size_t>    // <== indices of linear DoFs
    >;
    ThreadLocalStorage rTls; // <== defined only for the serial loop

    // Loop through elements and collect the DoFs' IDs
    // that are related to corner vertices.
    // @todo the parallel version of this loop relies on undefined behaviour,
    //       because all threads write to the same coarse mask at once => implement
    //       an ANY reduction. @matekelemen
    for (auto proxy : MakeProxy<TEntityLocation>(rModelPart)) {
    //block_for_each(MakeProxy<TEntityLocation>(rModelPart),
    //               ThreadLocalStorage(),
    //              [&r_process_info,&rA,&rCoarseMask](auto proxy, ThreadLocalStorage& rTls) {
        const auto& r_entity = proxy.GetEntity();
        if (r_entity.IsActive()) {
            rTls.first.clear();
            rTls.second.clear();
            const auto& r_geometry = r_entity.GetGeometry();
            const GeometryData::KratosGeometryFamily geometry_family = r_geometry.GetGeometryFamily();

            r_entity.GetDofList(rTls.first, r_process_info);
            GetLowerOrderNodeIndices(geometry_family,
                                     std::back_inserter(rTls.second));

            // The number of DoFs must be a multiple of the number of
            // (lower order) nodes, because each node should have the same DoFs.
            KRATOS_DEBUG_ERROR_IF(rTls.first.size() % rTls.second.size())
                << "DoF count mismatch in "
                << (std::is_same_v<std::remove_const_t<std::remove_reference_t<decltype(r_entity)>>,Element> ? "element" : "condition")
                << " " << r_entity.Id() << ": "
                << "number of DoFs (" << rTls.first.size() << ") "
                << "% number of nodes (" << rTls.second.size() << ")"
                << " != 0\n";

            const std::size_t dofs_per_node = rTls.first.size() / r_geometry.size();
            const std::size_t lower_order_dof_count = rTls.second.size() * dofs_per_node;

            for (std::size_t i_dof=0ul; i_dof<lower_order_dof_count; ++i_dof) {
                if (rTls.first[i_dof]->IsFree()) {
                    const std::size_t equation_id = rTls.first[i_dof]->EquationId();
                    KRATOS_DEBUG_ERROR_IF_NOT(equation_id < TSparseSpace::Size1(rA));
                    rCoarseMask[equation_id] = true;
                }
            }
        } // if (r_entity.IsActive())
    //});
    }

    KRATOS_CATCH("")
}


namespace Detail {


struct IndexPairTraits
{
    struct IsLess {
        bool operator()(std::pair<std::size_t,std::size_t> Left,
                        std::pair<std::size_t,std::size_t> Right) const noexcept
        {
            if (Left.first < Right.first) {
                return true;
            } else if (Left.first == Right.first) {
                return Left.second < Right.second;
            } else {
                return false;
            }
        }
    }; // struct IsLess
}; // struct IndexPairTraits


} // namespace Detail



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::ProvideAdditionalData(SparseMatrix& rA,
                                                                                        Vector& rX,
                                                                                        Vector& rB,
                                                                                        ModelPart::DofsArrayType& rDofs,
                                                                                        ModelPart& rModelPart)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

    const std::size_t system_size = TSparseSpace::Size1(rA);

    // Construct contraction operator
    KRATOS_INFO_IF("PolyhierarchicalSolver", 3 <= mpImpl->mVerbosity)
        << ((!TSparseSpace::Size1(mpImpl->mContractionOperator) || !TSparseSpace::Size2(mpImpl->mContractionOperator)) ? "construct" : "reconstruct")
        << " contraction operator\n";

    std::vector<bool> coarse_mask(TSparseSpace::Size1(rA), false); // <== yes, an intentional vector<bool>
    GetLowerOrderDofs<Globals::DataLocation::Element,TSparseSpace>(rA,
                                                                   rModelPart,
                                                                   coarse_mask);
    GetLowerOrderDofs<Globals::DataLocation::Condition,TSparseSpace>(rA,
                                                                     rModelPart,
                                                                     coarse_mask);

    // Compute an index map associating fine DoF indices with coarse ones.
    std::vector<std::size_t> coarse_dof_indices(system_size);
    std::size_t masked_count = 0;
    for (std::size_t i_mask=0; i_mask<system_size; ++i_mask) {
        const auto mask_value = coarse_mask[i_mask];
        coarse_dof_indices[i_mask] = (masked_count += mask_value);
    }

    typename TSparseSpace::MatrixType& r_contraction_operator = mpImpl->mContractionOperator;
    TSparseSpace::Resize(r_contraction_operator,
                         masked_count,
                         system_size);

    {
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        std::map<
            std::pair<std::size_t,std::size_t>,
            double,
            Detail::IndexPairTraits::IsLess
        > sparse_matrix_map;

        for (const auto& r_entity : rModelPart.Elements()) {
            if (r_entity.IsActive()) {
                const Geometry<Node>& r_fine_geometry = r_entity.GetGeometry();
                const Geometry<Node>& r_coarse_geometry = GetLowerOrderGeometry(r_fine_geometry);

                Element::EquationIdVectorType equation_ids;
                r_entity.EquationIdVector(equation_ids, r_process_info);

                const std::size_t fine_vertex_count = r_fine_geometry.PointsNumber();
                const std::size_t coarse_vertex_count = r_coarse_geometry.PointsNumber();
                const std::size_t dofs_per_node = equation_ids.size() / fine_vertex_count;

                for (std::size_t i_local_node_dof=0ul; i_local_node_dof<dofs_per_node; ++i_local_node_dof) {
                    for (std::size_t i_coarse_vertex=0ul; i_coarse_vertex<coarse_vertex_count; ++i_coarse_vertex) {
                        const std::size_t i_local_coarse_dof = i_coarse_vertex * dofs_per_node + i_local_node_dof;
                        const std::size_t i_coarse_dof = coarse_dof_indices[equation_ids[i_local_coarse_dof]];

                        for (std::size_t i_fine_vertex=0ul; i_fine_vertex<fine_vertex_count; ++i_fine_vertex) {
                            const std::size_t i_local_fine_dof = i_fine_vertex * dofs_per_node + i_local_node_dof;
                            const std::size_t i_fine_dof = equation_ids[i_local_fine_dof];
                            sparse_matrix_map.emplace(std::make_pair(i_coarse_dof, i_fine_dof), 0.5);
                        }

                        sparse_matrix_map[std::make_pair(i_coarse_dof, i_coarse_dof)] = 1.5;
                    }
                }
            } // if (r_entity.IsActive())
        }

        for (const auto& pair : sparse_matrix_map) {
            const std::size_t i_row = pair.first.first;
            const std::size_t i_column = pair.first.second;
            const double value = pair.second;
            r_contraction_operator.insert_element(i_row, i_column, value);
        }
    }

    KRATOS_ERROR_IF_NOT(TSparseSpace::Size1(mpImpl->mContractionOperator) == masked_count);
    KRATOS_ERROR_IF_NOT(TSparseSpace::Size2(mpImpl->mContractionOperator) == system_size);

    // Compute coarse system components
    // A_ll = C * A * C^T
    {
        typename TSparseSpace::MatrixType tmp = prod(mpImpl->mContractionOperator, rA);
        noalias(mpImpl->mCoarseSystem.mA) = prod(tmp, trans(mpImpl->mContractionOperator));
    }

    // x_l = C * x
    noalias(mpImpl->mCoarseSystem.mX) = prod(mpImpl->mContractionOperator, rX);

    // b_l = C * b
    noalias(mpImpl->mCoarseSystem.mB) = prod(mpImpl->mContractionOperator, rB);

    if (mpImpl->mCoarseSystem.mpSolver->AdditionalPhysicalDataIsNeeded()) {
        ModelPart::DofsArrayType coarse_dofs;
        for (auto it_dof=rDofs.ptr_begin(); it_dof!=rDofs.ptr_end(); ++it_dof) {
            if (coarse_mask[(*it_dof)->EquationId()]) {
                coarse_dofs.push_back(*it_dof);
            }
        }

        mpImpl->mCoarseSystem.mpSolver->ProvideAdditionalData(mpImpl->mCoarseSystem.mA,
                                                              mpImpl->mCoarseSystem.mX,
                                                              mpImpl->mCoarseSystem.mB,
                                                              coarse_dofs,
                                                              rModelPart);
    }

    if (mpImpl->mpFineSolver->AdditionalPhysicalDataIsNeeded()) {
        mpImpl->mpFineSolver->ProvideAdditionalData(rA,
                                                    rX,
                                                    rB,
                                                    rDofs,
                                                    rModelPart);
    }

    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
Parameters
PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::GetDefaultParameters()
{
    return Parameters(R"(
{
    "solver_type" : "poly_hierarchical",
    "verbosity" : 0,
    "max_iterations" : 50,
    "tolerance" : 1e-6,
    "coarse_settings" : {},
    "fine_settings" : {}
}
    )");
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                        DenseMatrix& rX,
                                                                        DenseMatrix& rB)
{
    return false;
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintInfo(std::ostream& rStream) const
{
    rStream << "PolyHierarchicalSolver";
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintData(std::ostream& rStream) const
{
    rStream
        << "tolerance     : " << mpImpl->mTolerance << "\n"
        << "max iterations: " << mpImpl->mMaxIterations << "\n"
        << "verbosity     : " << mpImpl->mVerbosity << "\n"
        ;
}



template
class PolyHierarchicalSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>,
    Reorderer<
        TUblasSparseSpace<double>,
        TUblasDenseSpace<double>
    >
>;


} // namespace Kratos
