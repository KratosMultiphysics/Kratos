#define MCGS_INTERNAL

// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::color
#include "multithreading.hpp" // MCGS_MUTEX_ARRAY
#include "defineMacros.hpp" // MCGS_EXPORT_SYMBOL

// --- STL Includes ---
#include <vector> // std::vector
#include <cstddef> // std::size_t
#include <algorithm> // std::min, std::max, std::equal_range
#include <numeric> // std::iota
#include <random> // std::mt19937, std::uniform_int_distribution
#include <iostream> // std::cout, std::cerr


namespace mcgs {


namespace detail {


template <class TIndex>
struct IndexPairTraits
{
    using value_type = std::pair<TIndex,TIndex>;

    struct Less
    {
        bool operator()(value_type left, value_type right) const noexcept
        {
            if (left.first < right.first) {
                return true;
            } else if (left.first == right.first) {
                return left.second < right.second;
            } else {
                return false;
            }
        }
    }; // struct Less
}; // struct IndexPairTraits


} // namespace detail


template <class TIndex>
using NeighborSet = std::vector<TIndex>;


/// @brief Collect all edges of an undirected graph.
template <class TIndex, class TValue>
std::vector<NeighborSet<TIndex>> collectNeighbors(const CSRAdaptor<TIndex,TValue>& rMatrix,
                                                  const ColorSettings<TValue> settings,
                                                  [[maybe_unused]] MCGS_MUTEX_ARRAY& rMutexes)
{
    std::vector<NeighborSet<TIndex>> neighbors(rMatrix.columnCount);

    #ifdef MCGS_OPENMP
    #pragma omp parallel for
    #endif
    for (int iRow=0; iRow<static_cast<int>(rMatrix.rowCount); ++iRow) {
        const TIndex iRowBegin = rMatrix.pRowExtents[iRow];
        const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];

        for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
            const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
            const TValue value = rMatrix.pEntries[iEntry];
            if (settings.tolerance <= std::abs(value) && static_cast<TIndex>(iRow) != iColumn) {
                {
                    MCGS_ACQUIRE_MUTEX(rMutexes[iRow]);
                    const auto [itBegin, itEnd] = std::equal_range(neighbors[iRow].begin(),
                                                                   neighbors[iRow].end(),
                                                                   iColumn);
                    if (itBegin == itEnd) neighbors[iRow].insert(itBegin, iColumn);
                    MCGS_RELEASE_MUTEX(rMutexes[iRow]);
                }

                {
                    MCGS_ACQUIRE_MUTEX(rMutexes[iColumn]);
                    const auto [itBegin, itEnd] = std::equal_range(neighbors[iColumn].begin(),
                                                                   neighbors[iColumn].end(),
                                                                   static_cast<TIndex>(iRow));
                    if (itBegin == itEnd) neighbors[iColumn].insert(itBegin, iRow);
                    MCGS_RELEASE_MUTEX(rMutexes[iColumn]);
                }
            }
        } // for iEntry in range(iRowBegin, iRowEnd)
    } // for iRow in range(rowCount)

    return neighbors;
}


using Mask = std::vector<char>;


template <class TIndex, class TColor>
bool isColored(const TIndex iVertex,
               const NeighborSet<TIndex>* pNeighborMap,
               const TColor* pColors,
               const Mask& rColoredMask)
{
    const TColor currentColor = pColors[iVertex];

    // If there's a conflict, keep the coloring of the vertex with the higher index.
    bool colored = true;

    for (const TIndex iNeighbor : pNeighborMap[iVertex]) {
        const TColor neighborColor = pColors[iNeighbor];
        if (neighborColor == currentColor) {
            if (iVertex < iNeighbor) {
                // The current vertex has a lower index than the neighbor
                // it's in conflict with => give up on this vertex.
                colored = false;
                break;
            } else if (rColoredMask[iNeighbor]) {
                // Although the current vertex would win a tiebreaker against
                // its neighbor it's in conflict with, the neighbor's color is
                // already set and cannot be changed => give up on this vertex.
                colored = false;
                break;
            }

            // Otherwise, the current vertex has a higher index
            // than its conflicting neighbor, who is still waiting
            // to be colored, winning the tiebreaker
            // => hang on to this vertex.
        }
    } // for iNeighbor in neighbors[iVertex]

    return colored;
}


template <class TColor>
struct Palette
{
    TColor maxColor;

    std::vector<TColor> palette;
}; // struct Palette


template <class TColor>
void extendPalette(Palette<TColor>& rPalette, const TColor color)
{
    rPalette.palette.push_back(color);
    rPalette.maxColor = std::max(color, rPalette.maxColor);
}


template <class TColor>
void extendPalette(Palette<TColor>& rPalette)
{
    extendPalette(rPalette, rPalette.maxColor + 1);
}


template <class TColor>
void removeFromPalette(const TColor color,
                       Palette<TColor>& rPalette)
{
    if (rPalette.maxColor < color) return;

    // The palette's colors are assumed to be sorted
    const auto [itBegin, itEnd] = std::equal_range(rPalette.palette.begin(),
                                                   rPalette.palette.end(),
                                                   color);
    rPalette.palette.erase(itBegin, itEnd);
}


template <class TIndex, class TValue, class TColor>
MCGS_EXPORT_SYMBOL
int color(TColor* pColors,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const ColorSettings<TValue> settings)
{
    // Cheap sanity checks
    if (rMatrix.rowCount < 0) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: invalid number of rows " << rMatrix.rowCount << "\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.columnCount < 0) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: invalid number of columns " << rMatrix.columnCount << "\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.entryCount < 0) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: invalid number of nonzeros " << rMatrix.entryCount << "\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pRowExtents) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: missing row data\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pColumnIndices) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: missing column data\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pEntries) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: missing nonzeros\n";
        return MCGS_FAILURE;
    }

    if (!pColors) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: missing output array\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.rowCount != rMatrix.columnCount) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: expecting a square matrix, but got "
                                              << rMatrix.rowCount << "x" << rMatrix.columnCount << "\n";
        return MCGS_FAILURE;
    }

    MCGS_MUTEX_ARRAY mutexes(rMatrix.rowCount);
    for ([[maybe_unused]] MCGS_MUTEX& rMutex : mutexes) MCGS_INITIALIZE_MUTEX(rMutex);

    // Collect all edges of the graph
    // (symmetric version of the input matrix)
    const auto neighbors = collectNeighbors(rMatrix, settings, mutexes);

    // Find the minimum and maximum vertex degrees.
    TIndex minDegree = std::numeric_limits<TIndex>::max();
    TIndex maxDegree = 0;

    #if defined(MCGS_OPENMP) && 201107 <= _OPENMP
    #pragma omp parallel for reduction(min: minDegree) reduction(max: maxDegree)
    #endif
    for (int iRow=0; iRow<static_cast<int>(rMatrix.rowCount); ++iRow) {
        const TIndex degree = static_cast<TIndex>(neighbors[iRow].size());
        minDegree = std::min(minDegree, degree);
        maxDegree = std::max(maxDegree, degree);
    } // for iRow in range(rowCount)

    if (2 <= settings.verbosity) {
        std::cout << "mcgs: max vertex degree is " << maxDegree << '\n'
                  << "mcgs: min vertex degree is " << minDegree << '\n'
                  ;
    }

    // Allocate the palette of every vertex to the max possible (maximum vertex degree).
    // An extra entry at the end of each palette indicates the palette's actual size.
    const TIndex shrinkingFactor = 0 < settings.shrinkingFactor ?
                                   static_cast<TIndex>(settings.shrinkingFactor) :
                                   std::max(TIndex(1), minDegree);

    const TIndex initialPaletteSize = std::max(
        TIndex(1),
        TIndex(double(maxDegree) / double(shrinkingFactor)));

    if (3 <= settings.verbosity) {
        std::cout << "mcgs: initial palette size is " << initialPaletteSize << std::endl;
    }

    // Initialize the palette of all vertices to a shrunk set
    std::vector<Palette<TColor>> palettes(rMatrix.columnCount);

    #ifdef MCGS_OPENMP
    #pragma omp parallel for
    #endif
    for (int iVertex=0; iVertex<static_cast<int>(rMatrix.columnCount); ++iVertex) {
        palettes[iVertex].palette.resize(initialPaletteSize);
        std::iota(palettes[iVertex].palette.begin(), palettes[iVertex].palette.end(), TColor(0));
        palettes[iVertex].maxColor = palettes[iVertex].palette.back();
    } // for iVertex in range(columnCount)

    // Track vertices that need to be colored.
    Mask coloredMask(rMatrix.columnCount, false);
    std::vector<TIndex> uncolored(rMatrix.rowCount);
    std::iota(uncolored.begin(), uncolored.end(), TIndex(0));

    #ifdef MCGS_OPENMP
    [[maybe_unused]] const int threadCount = omp_get_max_threads();
    #else
    [[maybe_unused]] const int threadCount = 1;
    #endif

    // Keep coloring until all vertices are colored.
    std::size_t iterationCount = 0ul;
    int stallCounter = 0;

    while (!uncolored.empty()) {
        const TIndex uncoloredCount = uncolored.size();
        if (3 <= settings.verbosity) {
            std::cout << "mcgs: coloring iteration " << iterationCount++
                      << " (" << uncoloredCount << "/" << rMatrix.columnCount
                      << " left to color)\n";
        }

        // Assign random colors to each remaining vertex from their palette.
        #ifdef MCGS_OPENMP
        #pragma omp parallel
        #endif
        {

            {
                #ifdef MCGS_OPENMP
                    const int iThread = omp_get_thread_num();
                #else
                    const int iThread = 0;
                #endif
                constexpr int seed = 0;
                std::mt19937 randomGenerator(seed);

                // To produce reproducible results, the random number
                // generator must be consistent for all possible number
                // of threads. To emulate what we'd get from a single-threaded
                // run, we need to partition the work manually and discard
                // the sequence of random numbers generated by other threads
                // before the current one.
                const TIndex uiThread = static_cast<TIndex>(iThread);
                const TIndex minChunkSize = uncoloredCount / static_cast<TIndex>(threadCount);
                const TIndex chunkLeftovers = uncoloredCount % static_cast<TIndex>(threadCount);
                const TIndex iChunkBegin = uiThread * minChunkSize + std::min(chunkLeftovers, uiThread);
                const TIndex iChunkEnd = iChunkBegin + minChunkSize + (chunkLeftovers <= uiThread ? 0 : 1);

                if (iChunkBegin != iChunkEnd) {
                    randomGenerator.discard(iChunkBegin);
                }

                for (TIndex iVisit=iChunkBegin; iVisit<iChunkEnd; ++iVisit) {
                    const TIndex iVertex = uncolored[iVisit];
                    const TColor paletteSize = palettes[iVertex].palette.size();
                    //const TColor iColorMax = paletteSize ? paletteSize - 1 : static_cast<TColor>(0);
                    //const TColor iColor = std::uniform_int_distribution<TColor>(TColor(0), iColorMax)(randomGenerator);
                    const TColor iColor = randomGenerator() % paletteSize;
                    pColors[iVertex] = palettes[iVertex].palette[iColor];
                } // for iVisit in range(iChunkBegin, iChunkEnd)

                #ifdef MCGS_OPENMP
                #pragma omp barrier
                #endif
            }

            #ifdef MCGS_OPENMP
            #pragma omp for
            #endif
            for (int iVisit=0; iVisit<static_cast<int>(uncoloredCount); ++iVisit) {
                const TIndex iVertex = uncolored[iVisit];
                const bool colored = isColored(iVertex, neighbors.data(), pColors, coloredMask);

                // If the current vertex has a valid color, remove
                // its color from the palettes of its neighbors.
                if (colored) {
                    MCGS_ACQUIRE_MUTEX(mutexes[iVertex]);
                    coloredMask[iVertex] = true;
                    palettes[iVertex].palette = std::vector<TColor> {};
                    MCGS_RELEASE_MUTEX(mutexes[iVertex]);


                    for (TIndex iNeighbor : neighbors[iVertex]) {
                        if (!coloredMask[iNeighbor]) {
                            MCGS_ACQUIRE_MUTEX(mutexes[iNeighbor]);
                            removeFromPalette(pColors[iVertex], palettes[iNeighbor]);
                            MCGS_RELEASE_MUTEX(mutexes[iNeighbor]);
                        } // if !coloredMask[iNeighbor]
                    } // for iNeighbor in neighbors[iVertex]
                } // if colored
            } // for iVertex in uncolored

            #ifdef MCGS_OPENMP
            #pragma omp for
            #endif
            for (int iVisit=0; iVisit<static_cast<int>(uncoloredCount); ++iVisit) {
                const TIndex iVertex = uncolored[iVisit];
                const bool needsExtension = !coloredMask[iVertex] && palettes[iVertex].palette.empty();
                if (needsExtension) {
                    TColor extension = palettes[iVertex].maxColor + 1;
                    extendPalette(palettes[iVertex], extension);
                } // if needsExtension
            } // for iVisit in range(uncoloredCount)
        } // omp parallel

        // Update remaining vertices
        {
            std::vector<TIndex> swap;
            swap.reserve(uncolored.size());

            for (TIndex iVisit=0; iVisit<static_cast<TIndex>(uncolored.size()); ++iVisit) {
                const TIndex iVertex = uncolored[iVisit];
                if (!coloredMask[iVertex]) swap.push_back(iVertex);
            }

            uncolored.swap(swap);
        }

        if (static_cast<TIndex>(uncolored.size()) == uncoloredCount) {
            // Failed to color any vertices => extend the palette of some random vertices
            const TIndex maxExtensions = std::max(TIndex(1), TIndex(25 * uncolored.size() / 100));
            TIndex extensionCounter = 0;

            // @todo parallelize
            for (TIndex iUncolored=0; iUncolored<static_cast<TIndex>(uncolored.size()) && extensionCounter < maxExtensions; ++iUncolored) {
                extendPalette(palettes[uncolored[iUncolored]]);
                ++extensionCounter;
            }

            if (!extensionCounter) {
                ++stallCounter;
                if (0 <= settings.maxStallCount && settings.maxStallCount <= stallCounter) {
                    if (1 <= settings.verbosity) std::cerr << "mcgs: error: reached stall limit (" << settings.maxStallCount << ")\n";
                    return MCGS_FAILURE;
                }
            } else {
                stallCounter = 0;
            }
        } else {
            stallCounter = 0;
        }
    };

    for ([[maybe_unused]] MCGS_MUTEX& rMutex : mutexes) MCGS_DEINITIALIZE_MUTEX(rMutex);

    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_COLOR(TIndex, TValue, TColor)                                 \
    template MCGS_EXPORT_SYMBOL int color(TColor* pColors,                             \
                                          const CSRAdaptor<TIndex,TValue>& rMatrix,    \
                                          const ColorSettings<TValue> settings)

MCGS_INSTANTIATE_COLOR(int, double, unsigned);
MCGS_INSTANTIATE_COLOR(long, double, unsigned);
MCGS_INSTANTIATE_COLOR(unsigned, double, unsigned);
MCGS_INSTANTIATE_COLOR(unsigned long, double, unsigned);

#undef MCGS_INSTANTIATE_COLOR
#undef MCGS_INTERNAL
#include "undefineMacros.hpp"


} // namespace mcgs
