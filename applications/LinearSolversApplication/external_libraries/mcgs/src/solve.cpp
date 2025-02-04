#define MCGS_INTERNAL

// --- External Includes ---
#ifdef MCGS_OPENMP
#include <omp.h> // omp_get_num_threads
#endif

// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::solve, mcgs::CSRAdaptor
#include "partition.hpp" // mcgs::Partition
#include "defineMacros.hpp" // MCGS_EXPORT_SYMBOL

// --- STL Includes ---
#include <cstddef> // std::size_t
#include <vector> // std::vector
#include <algorithm> // std::copy, std::clamp
#include <cmath> // std::sqrt
#include <iostream> // std::cout, std::cerr


namespace mcgs {


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
TValue residual(const CSRAdaptor<TIndex,TValue>& rMatrix,
                const TValue* pSolution,
                const TValue* pRHS) noexcept
{
    TValue residual = 0;

    #ifdef MCGS_OPENMP
    #pragma omp parallel for reduction(+: residual)
    #endif
    for (int iRow=0; iRow<static_cast<int>(rMatrix.rowCount); ++iRow) {
        TValue residualComponent = pRHS[iRow];
        const TIndex iRowBegin = rMatrix.pRowExtents[iRow];
        const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];

        for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
            const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
            residualComponent -= rMatrix.pEntries[iEntry] * pSolution[iColumn];
        } // for iEntry in range(iRowBegin, iRowEnd)

        residual += residualComponent * residualComponent;
    } // for iRow in range(0, rowCount)

    return std::sqrt(residual);
}


/// @brief Do one Gauss-Seidel iteration in serial.
template <class TIndex, class TValue>
int sweep(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const unsigned long iRowBegin,
          const unsigned long iRowEnd,
          const SolveSettings<TIndex,TValue> settings)
{
    for (TIndex iRow=iRowBegin; iRow<static_cast<TIndex>(iRowEnd); ++iRow) {
        TValue value = pRHS[iRow];
        TValue diagonal = 1;

        const TIndex iEntryBegin = rMatrix.pRowExtents[iRow];
        const TIndex iEntryEnd = rMatrix.pRowExtents[iRow + 1];

        for (TIndex iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
            const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
            const TValue nonzero = rMatrix.pEntries[iEntry];

            if (iRow == iColumn) diagonal = nonzero;
            else value -= nonzero * pSolution[iColumn];
        } /*for iEntry in range(iEntryBegin, iEntryEnd)*/

        pSolution[iRow] += settings.relaxation * (value / diagonal - pSolution[iRow]);
    }

    return MCGS_SUCCESS;
}


/// @brief Do one Gauss-Seidel iteration, distributing jobs row-wise.
template <class TIndex, class TValue>
int rowWiseSweep(TValue* pSolution,
                 const TValue* pSolutionBuffer,
                 const CSRAdaptor<TIndex,TValue>& rMatrix,
                 const TValue* pRHS,
                 const SolveSettings<TIndex,TValue> settings,
                 const TIndex iRowBegin,
                 const TIndex iRowEnd,
                 [[maybe_unused]] const int threadCount)
{
    #ifdef MCGS_OPENMP
    #pragma omp parallel for num_threads(threadCount)
    #endif
    for (int iRow=iRowBegin; iRow<static_cast<int>(iRowEnd); ++iRow) {
        TValue value = pRHS[iRow];
        TValue diagonal = 1;

        const TIndex iEntryBegin = rMatrix.pRowExtents[iRow];
        const TIndex iEntryEnd = rMatrix.pRowExtents[iRow + 1];

        for (TIndex iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
            const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
            const TValue nonzero = rMatrix.pEntries[iEntry];

            if (iColumn < static_cast<TIndex>(iRow)) value -= nonzero * pSolution[iColumn];
            else if (static_cast<TIndex>(iRow) < iColumn) value -= nonzero * pSolutionBuffer[iColumn];
            else diagonal = nonzero;
        } /*for iEntry in range(iEntryBegin, iEntryEnd)*/

        pSolution[iRow] += settings.relaxation * (value / diagonal - pSolution[iRow]);
    } // omp parallel for

    return MCGS_SUCCESS;
}


/// @brief Do one Gauss-Seidel iteration, distributing jobs by chunks of nonzeros in the matrix.
template <class TIndex, class TValue>
int entrywiseSweep(TValue* pSolution,
                   const TValue* pSolutionBuffer,
                   const CSRAdaptor<TIndex,TValue>& rMatrix,
                   const TValue* pRHS,
                   const SolveSettings<TIndex,TValue> settings,
                   const TIndex iRowBegin,
                   const TIndex iRowEnd,
                   const int threadCount)
{
    const TIndex partitionRowCount = iRowEnd - iRowBegin;
    const auto itEntryBegin = rMatrix.pRowExtents + iRowBegin;
    const auto itEntryEnd = rMatrix.pRowExtents + iRowEnd;
    const TIndex iEntryBegin = *itEntryBegin;
    const TIndex iEntryEnd = *itEntryEnd;
    const TIndex entryCount = iEntryEnd - iEntryBegin;

    std::vector<TValue> diagonals(partitionRowCount);
    std::vector<TValue> updates(partitionRowCount);
    std::copy(pRHS + iRowBegin, pRHS + iRowEnd, updates.data());

    std::vector<TIndex> threadEntryExtents(threadCount + 1);
    threadEntryExtents.front() = iEntryBegin;
    {
        const TIndex chunkSize = entryCount / threadCount + (entryCount % threadCount ? 1 : 0);
        for (TIndex iEnd=1; iEnd<static_cast<TIndex>(threadCount) + 1; ++iEnd) {
            threadEntryExtents[iEnd] = std::min(
                iEntryEnd,
                threadEntryExtents[iEnd - 1] + chunkSize
            );
        }
    }

    #ifdef MCGS_OPENMP
    #pragma omp parallel num_threads(threadCount)
    #endif
    {
        std::vector<TValue> localUpdates(partitionRowCount, static_cast<TValue>(0));

        #ifdef MCGS_OPENMP
        const TIndex iThread = omp_get_thread_num();
        #else
        const TIndex iThread = 0;
        #endif

        const TIndex iLocalEntryBegin = threadEntryExtents[iThread];
        const TIndex iLocalEntryEnd = threadEntryExtents[iThread + 1];

        const auto itLocalEntryBegin = std::max(std::upper_bound(itEntryBegin, itEntryEnd, iLocalEntryBegin) - 1,
                                                itEntryBegin);
        TIndex iRow = std::distance(rMatrix.pRowExtents, itLocalEntryBegin);
        TIndex iLocalRow = iRow - iRowBegin;
        const TIndex* itRowEnd = rMatrix.pRowExtents + iRow + 1;

        for (TIndex iEntry=iLocalEntryBegin; iEntry<iLocalEntryEnd; ++iEntry) {
            while (*itRowEnd <= iEntry) {
                ++iRow;
                ++iLocalRow;
                ++itRowEnd;
            }

            const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
            const TValue nonzero = rMatrix.pEntries[iEntry];

            if (iColumn < iRow) {
                localUpdates[iLocalRow] -= nonzero * pSolution[iColumn];
            } else if (iRow < iColumn) {
                localUpdates[iLocalRow] -= nonzero * pSolutionBuffer[iColumn];
            } else {
                diagonals[iLocalRow] = nonzero;
            }
        } // for iEntry in range(iLocalEntryBegin, iLocalEntryEnd)

        #ifdef MCGS_OPENMP
        #pragma omp critical
        #endif
        {
            for (TIndex iLocal=0; iLocal<static_cast<TIndex>(updates.size()); ++iLocal) {
                updates[iLocal] += localUpdates[iLocal];
            }
        } // omp critical
    } // omp parallel

    #ifdef MCGS_OPENMP
    #pragma omp parallel for num_threads(threadCount)
    #endif
    for (int iRow=iRowBegin; iRow<static_cast<int>(iRowEnd); ++iRow) {
        const TIndex iLocalRow = iRow - iRowBegin;
        pSolution[iRow] += settings.relaxation * (updates[iLocalRow] / diagonals[iLocalRow] - pSolution[iRow]);
    }

    return MCGS_SUCCESS;
}


/// @brief Decide which parallelization strategy to use and perform a single Gauss-Seidel iteration.
template <class TIndex, class TValue>
int dispatchSweep(TValue* pSolution,
                  const TValue* pSolutionBuffer,
                  const CSRAdaptor<TIndex,TValue>& rMatrix,
                  const TValue* pRHS,
                  const SolveSettings<TIndex,TValue> settings,
                  const TIndex iRowBegin,
                  const TIndex iRowEnd,
                  const int threadCount)
{
    if (iRowEnd < iRowBegin) {
        if (1 <= settings.verbosity) {
            std::cerr << "mcgs: error: invalid range [" << iRowBegin << ", " << iRowEnd << "[\n";
        }
        return MCGS_FAILURE;
    }

    if (settings.parallelization == Parallelization::RowWise) {
        return rowWiseSweep(pSolution,
                            pSolutionBuffer,
                            rMatrix,
                            pRHS,
                            settings,
                            iRowBegin,
                            iRowEnd,
                            threadCount);
    } /*if settings.parallelization == RowWise*/ else if (settings.parallelization == Parallelization::EntryWise) {
        return entrywiseSweep(pSolution,
                                pSolutionBuffer,
                                rMatrix,
                                pRHS,
                                settings,
                                iRowBegin,
                                iRowEnd,
                                threadCount);
    } /*if settings.parallelization == Parallelization::EntryWise*/ else if (settings.parallelization == Parallelization::None) {
        return sweep(pSolution,
                     rMatrix,
                     pRHS,
                     iRowBegin,
                     iRowEnd,
                     settings);
    } // /*if settings.parallelization == Parallelization::None*/

    return MCGS_SUCCESS;
}


/// @brief Perform Gauss-Seidel iterations in serial.
template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const SolveSettings<TIndex,TValue> settings)
{
    std::vector<TValue> buffer(rMatrix.columnCount);
    const TValue initialResidual = 3 <= settings.verbosity ?
                                   residual(rMatrix, pSolution, pRHS) :
                                   static_cast<TValue>(1);

    for (TIndex iIteration=0; iIteration<settings.maxIterations; ++iIteration) {
        if (sweep(pSolution,
                  rMatrix,
                  pRHS,
                  TIndex(0),
                  rMatrix.rowCount,
                  settings) != MCGS_SUCCESS) {
            if (1 <= settings.verbosity) {
                std::cerr << "mcgs: error: serial Gauss-Seidel failed at iteration "
                          << iIteration
                          << '\n';
            }
            return MCGS_FAILURE;
        }

        if (3 <= settings.verbosity) {
            std::cout << "iteration " << iIteration
                      << " residual: "
                      << residual(rMatrix, pSolution, pRHS) / initialResidual
                      << "\n";
        }
    }

    return MCGS_SUCCESS;
}


/// @brief Perform Gauss-Seidel iterations in parallel.
template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const Partition<TIndex>* pPartition,
          const SolveSettings<TIndex,TValue> settings)
{
    // Query the maximum number of threads MCGS is allowed to use.
    #ifdef MCGS_OPENMP
    const int maxThreadCount = omp_get_max_threads();
    #else
    const int maxThreadCount = 1;
    #endif

    if (maxThreadCount < 1) {
        if (1 <= settings.verbosity) {
            std::cerr << "mcgs: error: maximum number of allowed threads must be at least 1\n";
        }
        return MCGS_FAILURE;
    }

    // Run serial smoothing if no parallelization is allowed
    // or only a single thread is available.
    if (settings.parallelization == Parallelization::None || maxThreadCount == 1) {
        return solve(pSolution,
                     rMatrix,
                     pRHS,
                     settings);
    }

    // Run parallel smoothing otherwise.
    else {
        // Collect how many threads should execute each partition.
        std::vector<int> threadCounts(pPartition->size());

        #ifdef MCGS_OPENMP
        #pragma omp parallel for
        #endif
        for (int iPartition=0; iPartition<static_cast<int>(pPartition->size()); ++iPartition) {
            // @todo Find a dynamic way of approximating the optimal load of a single thread.
            //std::size_t nonzeroCount = 0ul;
            //for (auto itPartition=pPartition->begin(iPartition); itPartition!=pPartition->end(iPartition); ++itPartition) {
            //    const TIndex iRow = *itPartition;
            //    nonzeroCount += rMatrix.pRowExtents[iRow + 1] - rMatrix.pRowExtents[iRow];
            //} // for itPartition in pPartition[iPartition]
            //
            //threadCounts[iPartition] = std::clamp(nonzeroCount / 1024,
            //                                      static_cast<std::size_t>(1),
            //                                      std::min(pPartition->size(iPartition),
            //                                               static_cast<std::size_t>(maxThreadCount)));
            threadCounts[iPartition] = maxThreadCount;
        }

        std::vector<TValue> buffer(rMatrix.columnCount);
        const TValue initialResidual = 3 <= settings.verbosity ?
                                    residual(rMatrix, pSolution, pRHS) :
                                    static_cast<TValue>(1);

        for (TIndex iIteration=0; iIteration<settings.maxIterations; ++iIteration) {
            std::copy(pSolution, pSolution + rMatrix.rowCount, buffer.data());

            for (typename Partition<TIndex>::size_type iPartition=0; iPartition<pPartition->size(); ++iPartition) {
                const auto threadCount = threadCounts[iPartition];
                if (pPartition->isContiguous()) {
                    if (dispatchSweep(pSolution,
                                      buffer.data(),
                                      rMatrix,
                                      pRHS,
                                      settings,
                                      *pPartition->begin(iPartition),
                                      *pPartition->end(iPartition),
                                      threadCount) != MCGS_SUCCESS) {
                        if (1 <= settings.verbosity) {
                            std::cerr << "mcgs: error: parallel Gauss-Seidel failed at iteration "
                                      << iIteration
                                      << " on partition "
                                      << iPartition
                                      << '\n';
                        }
                        return MCGS_FAILURE;
                    }
                } /*if pPartition->isContiguous()*/ else {
                    if (1 <= settings.verbosity) {
                        std::cerr << "mcgs: error: parallel Gauss-Seidel requires a reordered system\n";
                    }
                    return MCGS_FAILURE;
                }
            } // for iPartition in range(partitionCount)

            if (3 <= settings.verbosity) {
                std::cout << "mcgs: iteration " << iIteration
                        << ", residual "
                        << residual(rMatrix, pSolution, pRHS) / initialResidual
                        << '\n';
            } // if 3 <= settings.verbosity
        } // for iIteration in range(settings.maxIterations)
    }

    return MCGS_SUCCESS;
}



#define MCGS_INSTANTIATE_SOLVE(TIndex, TValue)                                                 \
    template MCGS_EXPORT_SYMBOL TValue residual(const CSRAdaptor<TIndex,TValue>&,              \
                                                const TValue*,                                 \
                                                const TValue*) noexcept;                       \
    template MCGS_EXPORT_SYMBOL int solve<TIndex,TValue>(TValue*,                              \
                                                         const CSRAdaptor<TIndex,TValue>&,     \
                                                         const TValue*,                        \
                                                         const SolveSettings<TIndex,TValue>);  \
    template MCGS_EXPORT_SYMBOL int solve<TIndex,TValue>(TValue*,                              \
                                                         const CSRAdaptor<TIndex,TValue>&,     \
                                                         const TValue*,                        \
                                                         const Partition<TIndex>*,             \
                                                         const SolveSettings<TIndex,TValue>)

MCGS_INSTANTIATE_SOLVE(int, double);
MCGS_INSTANTIATE_SOLVE(long, double);
MCGS_INSTANTIATE_SOLVE(unsigned, double);
MCGS_INSTANTIATE_SOLVE(unsigned long, double);
MCGS_INSTANTIATE_SOLVE(int, float);
MCGS_INSTANTIATE_SOLVE(long, float);
MCGS_INSTANTIATE_SOLVE(unsigned, float);
MCGS_INSTANTIATE_SOLVE(unsigned long, float);

#undef MCGS_INSTANTIATE_SOLVE
#undef MCGS_INTERNAL
#include "undefineMacros.hpp"


} // namespace mcgs
