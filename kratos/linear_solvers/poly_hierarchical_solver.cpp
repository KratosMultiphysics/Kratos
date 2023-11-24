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

// External includes
#include "eigen3/Eigen/Sparse"

// Project includes
#include "poly_hierarchical_solver.h"
#include "includes/code_location.h"
#include "includes/global_variables.h"
#include "spaces/ublas_space.h"
#include "utilities/profiler.h"
#include "factories/linear_solver_factory.h"
#include "utilities/proxies.h"
#include "utilities/sparse_matrix_multiplication_utility.h"

// STL includes
#include <boost/numeric/ublas/matrix_sparse.hpp>
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

    typename TSparseSpace::MatrixType mRestrictionOperator;

    typename TSparseSpace::MatrixType mInterpolationOperator;

    struct CoarseSystem {
        typename TSparseSpace::MatrixType mA;

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
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

    if (4 <= mpImpl->mVerbosity) {
        KRATOS_ERROR << "verbosity >=4 prints system matrices and terminates\n";
    }

    // Declarations
    const std::size_t coarse_size = TSparseSpace::Size1(mpImpl->mRestrictionOperator);
    const std::size_t fine_size = TSparseSpace::Size2(mpImpl->mRestrictionOperator);
    Vector fine_delta(fine_size);
    Vector fine_residual(fine_size);
    Vector fine_tmp(fine_size);
    Vector coarse_delta(coarse_size);
    Vector coarse_residual(coarse_size);

    const double rhs_norm = TSparseSpace::TwoNorm(rB);

    // Initializations
    // r = b - A * x
    TSparseSpace::Mult(rA, rX, fine_tmp);
    noalias(fine_residual) = rB - fine_tmp;

    double residual_norm = 1.0;
    const std::size_t max_iterations = mpImpl->mMaxIterations;
    for (std::size_t i_iteration=0ul; i_iteration<max_iterations; ++i_iteration) {
        // Perform a few smoothing passes on the full system.
        KRATOS_TRY
        TSparseSpace::SetToZero(fine_delta);
        mpImpl->mpFineSolver->Solve(rA, fine_delta, fine_residual);
        KRATOS_CATCH("")

        // Update residual
        TSparseSpace::Mult(rA, fine_delta, fine_tmp);
        fine_residual -= fine_tmp;

        // Restrict the updated residual
        TSparseSpace::Mult(mpImpl->mRestrictionOperator,
                           fine_residual,
                           coarse_residual);

        // Solve coarse system
        KRATOS_TRY
        TSparseSpace::SetToZero(coarse_delta);
        mpImpl->mCoarseSystem.mpSolver->Solve(mpImpl->mCoarseSystem.mA,
                                              coarse_delta,
                                              coarse_residual);
        KRATOS_CATCH("")

        // Update solution
        TSparseSpace::Mult(mpImpl->mInterpolationOperator, coarse_delta, fine_tmp);
        noalias(rX) += fine_delta + fine_tmp;

        // Update the fine residual
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


std::size_t GetLinearGeometryVertexCount(GeometryData::KratosGeometryFamily GeometryFamily)
{
    using Family = GeometryData::KratosGeometryFamily;
    switch(GeometryFamily) {
        case Family::Kratos_Point:          return 1;
        case Family::Kratos_Linear:         return 2;
        case Family::Kratos_Triangle:       return 3;
        case Family::Kratos_Quadrilateral:  return 4;
        case Family::Kratos_Tetrahedra:     return 4;
        case Family::Kratos_Pyramid:        return 5;
        case Family::Kratos_Hexahedra:      return 6;
        default: {
            KRATOS_ERROR << "unsupported geometry family";
        }
    } // switch (GeometryFamily)
}


template <Globals::DataLocation TEntityLocation, class TSparseSpace>
void GetLowerOrderDofs(typename TSparseSpace::MatrixType& rA,
                       const ModelPart& rModelPart,
                       std::vector<bool>& rCoarseMask)
{
    KRATOS_TRY

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    using ThreadLocalStorage = Element::DofsVectorType;
    ThreadLocalStorage rTls; // <== defined only for the serial loop

    // Loop through elements and collect the DoFs' IDs
    // that are related to corner vertices.
    // @todo the parallel version of this loop relies on undefined behaviour,
    //       because all threads write to the same coarse mask at once => implement
    //       an ANY reduction. @matekelemen
    for (auto proxy : MakeProxy<TEntityLocation>(rModelPart)) {
        const auto& r_entity = proxy.GetEntity();
        const Geometry<Node>& r_geometry = r_entity.GetGeometry();

        if (r_entity.IsActive() && r_geometry.PointsNumber()) {
            rTls.clear();
            const GeometryData::KratosGeometryFamily geometry_family = r_geometry.GetGeometryFamily();

            r_entity.GetDofList(rTls, r_process_info);

            // The number of DoFs must be a multiple of the number of
            // fine nodes, because each node should have the same DoFs.
            KRATOS_DEBUG_ERROR_IF(rTls.size() % r_geometry.PointsNumber())
                << "DoF count mismatch in "
                << (std::is_same_v<std::remove_const_t<std::remove_reference_t<decltype(r_entity)>>,Element> ? "element" : "condition")
                << " " << r_entity.Id() << ": "
                << "number of DoFs (" << rTls.size() << ") "
                << "% number of nodes (" << r_geometry.PointsNumber() << ")"
                << " != 0\n";

            const std::size_t dofs_per_node = rTls.size() / r_geometry.PointsNumber();
            const std::size_t coarse_node_count = GetLinearGeometryVertexCount(geometry_family);

            KRATOS_ERROR_IF_NOT(coarse_node_count <= r_geometry.PointsNumber())
                << "coarse node count " << coarse_node_count
                << " is higher than the number of fine vertices " << r_geometry.PointsNumber();

            for (std::size_t i_node=0ul; i_node<coarse_node_count; ++i_node) {
                if (r_geometry[i_node].IsNot(SLAVE)) {
                    const std::size_t i_dof_begin = i_node * dofs_per_node;
                    for (std::size_t i_local_dof=0ul; i_local_dof<dofs_per_node; ++i_local_dof) {
                        const std::size_t i_equation = rTls[i_dof_begin + i_local_dof]->EquationId();
                        KRATOS_DEBUG_ERROR_IF_NOT(i_equation < TSparseSpace::Size1(rA));
                        rCoarseMask[i_equation] = true;
                    }
                } // if node is not SLAVE
            }
        } // if (r_entity.IsActive())
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



template <class TOutputIterator>
void GetNonzeroQuadraticCoefficients(GeometryData::KratosGeometryType Geometry,
                                     TOutputIterator itOutput)
{
    using G = GeometryData::KratosGeometryType;
    using Pair = std::pair<
        std::size_t,    // index of the quadratic shape function
        double          // coefficient of the quadratic shape function
    >;
    using PairVector = std::vector<Pair>;
    switch (Geometry) {
        case G::Kratos_Point2D: {
            *itOutput = std::make_pair(0ul, PairVector {{0ul, 1.0}});
            break;
        }
        case G::Kratos_Point3D: {
            *itOutput = std::make_pair(0ul, PairVector {{0ul, 1.0}});
            break;
        }
        case G::Kratos_Line2D2: {
            *itOutput = std::make_pair(0ul, PairVector {{0ul, 1.0}});
            *itOutput = std::make_pair(1ul, PairVector {{1ul, 1.0}});
            break;
        }
        case G::Kratos_Line2D3: {
            *itOutput = std::make_pair(0ul, PairVector {{0ul, 1.0},
                                                        {2ul, 0.5}});
            *itOutput = std::make_pair(1ul, PairVector {{1ul, 1.0},
                                                        {2ul, 0.5}});
            break;
        }
        case G::Kratos_Line3D2: {
            *itOutput = std::make_pair(0ul, PairVector {{0ul, 1.0}});
            *itOutput = std::make_pair(1ul, PairVector {{1ul, 1.0}});
            break;
        }
        case G::Kratos_Line3D3: {
            *itOutput = std::make_pair(0ul, PairVector {{0ul, 1.0},
                                                        {2ul, 0.5}});
            *itOutput = std::make_pair(1ul, PairVector {{1ul, 1.0},
                                                        {2ul, 0.5}});
            break;
        }
        case G::Kratos_Triangle2D6: {
            *itOutput++ = std::make_pair(0ul,
                                         PairVector {
                                            {0ul, 1.0},
                                            {3ul, 0.5},
                                            {5ul, 0.5}
                                         });
            *itOutput++ = std::make_pair(1ul,
                                         PairVector {
                                            {1ul, 1.0},
                                            {3ul, 0.5},
                                            {4ul, 0.5}
                                         });
            *itOutput++ = std::make_pair(2ul,
                                         PairVector {
                                            {2ul, 1.0},
                                            {4ul, 0.5},
                                            {5ul, 0.5}
                                         });
            break;
        }
        case G::Kratos_Triangle3D6: {
            *itOutput++ = std::make_pair(0ul,
                                         PairVector {
                                            {0ul, 1.0},
                                            {3ul, 0.5},
                                            {5ul, 0.5}
                                         });
            *itOutput++ = std::make_pair(1ul,
                                         PairVector {
                                            {1ul, 1.0},
                                            {3ul, 0.5},
                                            {4ul, 0.5}
                                         });
            *itOutput++ = std::make_pair(2ul,
                                         PairVector {
                                            {2ul, 1.0},
                                            {4ul, 0.5},
                                            {5ul, 0.5}
                                         });
            break;
        }
        case G::Kratos_Tetrahedra3D10: {
            *itOutput++ = std::make_pair(0ul,
                                         PairVector {
                                            {0ul, 1.0},
                                            {4ul, 0.5},
                                            {6ul, 0.5},
                                            {7ul, 0.5}
                                         });
            *itOutput++ = std::make_pair(1ul,
                                         PairVector {
                                            {1ul, 1.0},
                                            {4ul, 0.5},
                                            {5ul, 0.5},
                                            {8ul, 0.5}
                                         });
            *itOutput++ = std::make_pair(2ul,
                                         PairVector {
                                            {2ul, 1.0},
                                            {5ul, 0.5},
                                            {6ul, 0.5},
                                            {9ul, 0.5}
                                         });
            *itOutput++ = std::make_pair(3ul,
                                         PairVector {
                                            {3ul, 1.0},
                                            {7ul, 0.5},
                                            {8ul, 0.5},
                                            {9ul, 0.5}
                                         });
            break;
        }
        default: KRATOS_ERROR << "unsupported geometry type " << (int) Geometry << "\n";
    }
}



template <Globals::DataLocation TEntityLocation>
void MapHigherToLowerOrder(const ModelPart& rModelPart,
                           const std::vector<bool>& rCoarseDofMask,
                           const std::vector<std::size_t>& rCoarseDofMap,
                           std::map<
                                std::pair<std::size_t,std::size_t>,
                                double,
                                Detail::IndexPairTraits::IsLess
                           >& rRestrictionMap,
                           std::map<
                                std::pair<std::size_t,std::size_t>,
                                double,
                                Detail::IndexPairTraits::IsLess
                           >& rInterpolationMap)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    KRATOS_TRY

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    for (auto proxy : MakeProxy<TEntityLocation>(rModelPart)) {
        const auto& r_entity = proxy.GetEntity();
        if (r_entity.IsActive()) {
            // Get the geometry of the current element, and its
            // associated lower order geometry.
            const Geometry<Node>& r_fine_geometry = r_entity.GetGeometry();
            const std::size_t fine_vertex_count = r_fine_geometry.PointsNumber();

            // Collect the DoF indices defined on the fine system.
            Element::EquationIdVectorType local_to_fine_dof_map;
            r_entity.EquationIdVector(local_to_fine_dof_map, r_process_info);

            KRATOS_DEBUG_ERROR_IF_NOT(fine_vertex_count <= local_to_fine_dof_map.size());
            KRATOS_DEBUG_ERROR_IF(local_to_fine_dof_map.size() % fine_vertex_count);
            const std::size_t dofs_per_node = local_to_fine_dof_map.size() / fine_vertex_count;

            std::vector<std::pair<
                std::size_t,        // <== coarse DoF index
                std::vector<std::pair<
                    std::size_t,    // <== fine DoF index
                    double          // <== fine DoF coefficient
                >>
            >> restriction_coefficients;
            restriction_coefficients.reserve(4ul); // <== reserve assuming a tetrahedron (linear vertex count)
            GetNonzeroQuadraticCoefficients(r_fine_geometry.GetGeometryType(),
                                            std::back_inserter(restriction_coefficients));

            for (const auto& r_restriction_pair : restriction_coefficients) {
                const std::size_t i_local_coarse_vertex = r_restriction_pair.first;
                const std::size_t local_coarse_dof_begin = i_local_coarse_vertex * dofs_per_node;

                for (const auto& r_fine_coefficient_pair : r_restriction_pair.second) {
                    const std::size_t i_local_fine_vertex = r_fine_coefficient_pair.first;
                    const std::size_t local_fine_dof_begin = i_local_fine_vertex * dofs_per_node;

                    for (std::size_t i_node_dof=0ul; i_node_dof<dofs_per_node; ++i_node_dof) {
                        const std::size_t i_coarse_dof_in_fine_system = local_to_fine_dof_map[local_coarse_dof_begin + i_node_dof];

                        if (rCoarseDofMask[i_coarse_dof_in_fine_system]) {
                            const std::size_t i_coarse_dof = rCoarseDofMap[i_coarse_dof_in_fine_system];
                            const std::size_t i_fine_dof = local_to_fine_dof_map[local_fine_dof_begin + i_node_dof];
                            rRestrictionMap.emplace(std::make_pair(i_coarse_dof, i_fine_dof), r_fine_coefficient_pair.second);
                            rInterpolationMap.emplace(std::make_pair(i_fine_dof, i_coarse_dof), r_fine_coefficient_pair.second);
                        }
                    }
                }
            }
        } // if (r_entity.IsActive())
    }
    KRATOS_CATCH("")
}



template <class TValue>
Eigen::SparseMatrix<TValue> MakeEigenMatrix(const typename TUblasSparseSpace<TValue>::MatrixType& rInput)
{
    Eigen::SparseMatrix<TValue,Eigen::RowMajor> output(rInput.size1(), rInput.size2());
    output.reserve(rInput.value_data().size());

    const auto it_row_end = rInput.end1();
    for (auto it_row=rInput.begin1(); it_row!=it_row_end; ++it_row) {
        const auto it_entry_end = it_row.end();
        for (auto it_entry=it_row.begin(); it_entry!=it_entry_end; ++it_entry) {
            const std::size_t i_row = it_entry.index1();
            const std::size_t i_column = it_entry.index2();
            const TValue value = *it_entry;
            output.insert(i_row, i_column) = value;
        }
    }

    return output;
}



template <class TValue>
void Assign(const Eigen::SparseMatrix<TValue>& rSource,
            typename TUblasSparseSpace<TValue>::MatrixType& rTarget)
{
    KRATOS_ERROR_IF_NOT(rSource.isCompressed());
    using EigenMatrix = Eigen::SparseMatrix<TValue>;
    rTarget.resize(rSource.rows(), rSource.cols());
    for (long i_outer=0l; i_outer<rSource.outerSize(); ++i_outer) {
        for (auto it_entry=typename EigenMatrix::InnerIterator(rSource,i_outer); it_entry; ++it_entry) {
            rTarget.insert_element(it_entry.row(), it_entry.col(), it_entry.value());
        }
    }
}



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

    if (4 <= mpImpl->mVerbosity) {
        TSparseSpace::WriteMatrixMarketMatrix("system_matrix.mm", rA, false);
    }

    const std::size_t system_size = TSparseSpace::Size1(rA);

    // Construct restriction operator
    KRATOS_INFO_IF("PolyhierarchicalSolver", 3 <= mpImpl->mVerbosity)
        << ((!TSparseSpace::Size1(mpImpl->mRestrictionOperator) || !TSparseSpace::Size2(mpImpl->mRestrictionOperator)) ? "construct" : "reconstruct")
        << " restriction operator\n";

    std::vector<bool> coarse_mask(system_size, false); // <== yes, an intentional vector<bool>
    GetLowerOrderDofs<Globals::DataLocation::Element,TSparseSpace>(rA,
                                                                   rModelPart,
                                                                   coarse_mask);

    // Conditions shouldn't introduce coarse DoFs, but I'll leave this
    // here for the time being in case I'm missing some exocotic boundary
    // condition implementations @matekelemen.
    //GetLowerOrderDofs<Globals::DataLocation::Condition,TSparseSpace>(rA,
    //                                                                 rModelPart,
    //                                                                 coarse_mask);

    // Compute an index map associating fine DoF indices with coarse ones
    // and collect coarse DoFs.
    // Indices at fine DoFs should not be accessed
    // (doing so will result in a segfault instead of silent miscalculations)
    std::vector<std::size_t> coarse_dof_indices(system_size,
                                                std::numeric_limits<std::size_t>::max());
    ModelPart::DofsArrayType coarse_dofs;
    coarse_dofs.reserve(std::size_t(0.4 * system_size)); // <== reserve assuming quadratic tets
    std::size_t masked_count = 0ul;

    for (std::size_t i_mask=0; i_mask<system_size; ++i_mask) {
        if (coarse_mask[i_mask]) {
            coarse_dof_indices[i_mask] = masked_count++;
            coarse_dofs.push_back(*(rDofs.ptr_begin() + i_mask));
        }
    }

    typename TSparseSpace::MatrixType& r_restriction_operator = mpImpl->mRestrictionOperator;
    typename TSparseSpace::MatrixType& r_interpolation_operator = mpImpl->mInterpolationOperator;

    {
        std::map<
            std::pair<std::size_t,std::size_t>,
            double,
            Detail::IndexPairTraits::IsLess
        > restriction_map, interpolation_map;
        MapHigherToLowerOrder<Globals::DataLocation::Element>(rModelPart,
                                                              coarse_mask,
                                                              coarse_dof_indices,
                                                              restriction_map,
                                                              interpolation_map);

        // Conditions shouldn't introduce coarse DoFs, but I'll leave this
        // here for the time being in case I'm missing some exotic boundary
        // condition implementations @matekelemen.
        //MapHigherToLowerOrder<Globals::DataLocation::Condition>(rModelPart,
        //                                                        coarse_dof_indices,
        //                                                        restriction_map,
        //                                                        interpolation_map);

        {
            KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

            TSparseSpace::Resize(r_restriction_operator,
                                 masked_count,
                                 system_size);
            for (auto pair : restriction_map) {
                const std::size_t i_row = pair.first.first;
                const std::size_t i_column = pair.first.second;
                KRATOS_ERROR_IF_NOT(i_row < masked_count);
                KRATOS_ERROR_IF_NOT(i_column < system_size);
                const double value = pair.second;
                r_restriction_operator.insert_element(i_row, i_column, value);
            }

            TSparseSpace::Resize(r_interpolation_operator,
                                 system_size,
                                 masked_count);
            for (auto pair : interpolation_map) {
                const std::size_t i_row = pair.first.first;
                const std::size_t i_column = pair.first.second;
                const double value = pair.second;
                r_interpolation_operator.insert_element(i_row, i_column, value);
            }
        } // end profiling
    } // destruct restriction_map and interpolation_map

    if (4 <= mpImpl->mVerbosity) {
        TSparseSpace::WriteMatrixMarketMatrix("restriction_operator.mm", mpImpl->mRestrictionOperator, false);
        TSparseSpace::WriteMatrixMarketMatrix("interpolation_operator.mm", mpImpl->mInterpolationOperator, false);
    }

    KRATOS_ERROR_IF_NOT(TSparseSpace::Size1(mpImpl->mRestrictionOperator) == masked_count);
    KRATOS_ERROR_IF_NOT(TSparseSpace::Size2(mpImpl->mRestrictionOperator) == system_size);
    KRATOS_ERROR_IF_NOT(TSparseSpace::Size1(mpImpl->mInterpolationOperator) == system_size);
    KRATOS_ERROR_IF_NOT(TSparseSpace::Size2(mpImpl->mInterpolationOperator) == masked_count);

    // Compute coarse system components
    // coarse_system = restriction_operator * fine_system * restriction_operator^T
    {
        KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
        typename TSparseSpace::MatrixType tmp;
        TSparseSpace::Resize(tmp, masked_count, system_size);
        TSparseSpace::Resize(mpImpl->mCoarseSystem.mA, masked_count, masked_count);

        KRATOS_INFO("PolyHierarchicalSolver")
                << "first multiplication begin "
                << "(" << mpImpl->mRestrictionOperator.size1() << ", " << mpImpl->mRestrictionOperator.size2() << ")" << " [" << mpImpl->mRestrictionOperator.nnz() << "]"
                << " x "
                << "(" << rA.size1() << ", " << rA.size2() << ")" << " [" << rA.nnz() << "]"
                << " = "
                << "(" << tmp.size1() << ", " << tmp.size2() << ")" << " [" << tmp.nnz() << "]"
                << "\n";
        auto eigen_restriction_operator = MakeEigenMatrix<typename TSparseSpace::DataType>(mpImpl->mRestrictionOperator);
        auto eigen_system_matrix = MakeEigenMatrix<typename TSparseSpace::DataType>(rA);
        //SparseMatrixMultiplicationUtility::MatrixMultiplication(mpImpl->mRestrictionOperator,
        //                                                        rA,
        //                                                        tmp);
        KRATOS_INFO("PolyHierarchicalSolver") << "first multiplication end\n";

        {
            // The interpolation operator can easily be converted into column-major format
            // since its transpose is already available (the restriction operator)
            // => create a temporary column-major interpolation operator for the duration
            //    of the matrix multiplication to perform it in a more cache-friendly manner.
            //using UblasCompressedMatrix = boost::numeric::ublas::compressed_matrix<typename TSparseSpace::DataType>;
            //if constexpr (std::is_same_v<typename TSparseSpace::MatrixType,UblasCompressedMatrix>)

            KRATOS_INFO("PolyHierarchicalSolver")
                << "first multiplication begin "
                << "(" << tmp.size1() << ", " << tmp.size2() << ")" << " [" << tmp.nnz() << "]"
                << " x "
                << "(" << mpImpl->mInterpolationOperator.size1() << ", " << mpImpl->mInterpolationOperator.size2() << ")" << " [" << mpImpl->mRestrictionOperator.nnz() << "]"
                << " = "
                << "(" << mpImpl->mCoarseSystem.mA.size1() << ", " << mpImpl->mCoarseSystem.mA.size2() << ")" << " [" << mpImpl->mCoarseSystem.mA.nnz() << "]"
                << "\n";
            auto eigen_interpolation_operator = MakeEigenMatrix<typename TSparseSpace::DataType>(mpImpl->mInterpolationOperator);
            decltype(eigen_system_matrix) eigen_coarse_system = eigen_restriction_operator * eigen_system_matrix * eigen_interpolation_operator;

            eigen_system_matrix.makeCompressed();
            Assign<typename TSparseSpace::DataType>(eigen_coarse_system, mpImpl->mCoarseSystem.mA);

            KRATOS_INFO("PolyHierarchicalSolver")
                << eigen_coarse_system.rows() << " "
                << eigen_coarse_system.cols() << " "
                << eigen_coarse_system.nonZeros() << " "
                << mpImpl->mCoarseSystem.mA.index1_data().size() << " "
                << mpImpl->mCoarseSystem.mA.index2_data().size() << " "
                << mpImpl->mCoarseSystem.mA.value_data().size() << " "
                << "\n";
            //SparseMatrixMultiplicationUtility::MatrixMultiplication(tmp,
            //                                                        mpImpl->mInterpolationOperator,
            //                                                        mpImpl->mCoarseSystem.mA);
            KRATOS_INFO("PolyHierarchicalSolver") << "second multiplication end\n";
        }

        if (4 <= mpImpl->mVerbosity) {
            TSparseSpace::WriteMatrixMarketMatrix("coarse_system_matrix.mm", mpImpl->mCoarseSystem.mA, false);
        }
    } // end profiling


    if (mpImpl->mCoarseSystem.mpSolver->AdditionalPhysicalDataIsNeeded()) {
        // Restrict solution and RHS vectors
        Vector coarse_x(masked_count), coarse_b(masked_count);
        TSparseSpace::Mult(mpImpl->mRestrictionOperator,
                           rX,
                           coarse_x);
        TSparseSpace::Mult(mpImpl->mRestrictionOperator,
                           rB,
                           coarse_b);

        // Provide additional physical data for the coarse system
        mpImpl->mCoarseSystem.mpSolver->ProvideAdditionalData(mpImpl->mCoarseSystem.mA,
                                                              coarse_x,
                                                              coarse_b,
                                                              coarse_dofs,
                                                              rModelPart);
    }

    if (mpImpl->mpFineSolver->AdditionalPhysicalDataIsNeeded()) {
        // Provide additional physical data for the fine system
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
