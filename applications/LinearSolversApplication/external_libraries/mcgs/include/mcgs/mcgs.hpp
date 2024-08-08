#pragma once

/// @def MCGS_EXPORT_SYMBOL
/// @brief Exposes the symbol in the shared object.
#if _WIN32
    #ifdef MCGS_INTERNAL
        #define MCGS_EXPORT_SYMBOL __declspec(dllexport)
    #else
        #define MCGS_EXPORT_SYMBOL __declspec(dllimport)
    #endif
#else
    #define MCGS_EXPORT_SYMBOL __attribute__((visibility ("default")))
#endif


namespace mcgs {


/// @def MCGS_SUCCESS
/// @brief Return value indicating success.
#define MCGS_SUCCESS 0


/// @def MCGS_FAILURE
/// @brief Return value indicating failure.
#define MCGS_FAILURE 1


/// @brief Adaptor for matrices stored in the compressed sparse row format.
///
/// @tparam TIndex Integer type of stored indices (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries (for now, only @p double is supported).
///
/// @details This class is a hollow adaptor that does not own any of the arrays associated
///          with the matrix, and does not have mutable access to them. Its purpose is restricted
///          to providing read access to the matrix.
///
/// @see <a href="https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)">CSR format</a>.
template <class TIndex, class TValue>
struct CSRAdaptor
{
    unsigned long rowCount;         ///< @brief Number of rows in the matrix.
    unsigned long columnCount;      ///< @brief Number of columns in the matrix.
    unsigned long entryCount;       ///< @brief Number of stored entries in the matrix.
    const TIndex* pRowExtents;      ///< @brief Index array defining entries related to each row within the @ref CSRAdaptor::pEntries "array of entries". Size at least @ref rowCount + 1.
    const TIndex* pColumnIndices;   ///< @brief Index array assigning column indices to each entry in the @ref CSRAdaptor::pEntries "array of entries". Size at least @ref entryCount.
    const TValue* pEntries;         ///< @brief Array of stored entries. Size at least @ref entryCount.

    /// @brief Default constructor creating an invalid object that must be initialized by the user.
    CSRAdaptor() noexcept
        : rowCount(0ul),
          columnCount(0ul),
          entryCount(0ul),
          pRowExtents(nullptr),
          pColumnIndices(nullptr),
          pEntries(nullptr)
    {}
};


/// @brief Settings for the coloring algorithm.
/// @tparam TValue Number type of stored entries in the related matrix (for now, only @p double is supported).
template <class TValue>
struct ColorSettings
{
    /// @brief Parameter controlling the initial palette size.
    /// @details This is the only algorithmic parameter, but it has a significant
    ///          influence on the quality of the resulting coloring as well as
    ///          the runtime of the algorithm. In general, a shrinking factor that
    ///          leads to a good coloring also has shorter runtimes.
    ///
    ///          The initial palette size @f$p@f$ depends on the maximum vertex degree
    ///          @f$v_{max}@f$ and shrinking factor @f$s@f$:
    ///          @f[ p = \frac{v_{max}}{s} @f]
    ///
    /// @note
    /// - If a nonpositive value is passed, the minimum vertex degree is used instead.
    /// - Defaults to @p -1.
    int shrinkingFactor;

    /// @brief Controls how many consequent failed coloring iterations to allow before terminating.
    /// @details Since the algorithm randomly picks colors for each vertex from their remaining
    ///          palettes, it is possible to end up with iterations that fail to assign valid colors
    ///          to any of the remaining vertices. This is characteristic of the final stages of the
    ///          coloring process for runs with too tight initial palettes. Allowing a longer chain
    ///          of failed iterations makes the algorithm more robust in these cases.
    ///
    /// @note
    /// - Negative values for @p maxStallCount allow infinite stalling.
    /// - Defaults to @p 1e3.
    int maxStallCount;

    /// @brief Minimum absolute value for matrix entries to treat as nonzero.
    /// @details Some software store actual zeros in their sparse matrices. For example, this can
    ///          happen during the naive enforcement of Dirichlet conditions or linear constraints
    ///          in FEM software. Topology complexity can be reduced in these cases by defining
    ///          a small tolerance, under which entries in the sparse matrix are ignored.
    /// @note
    /// - Nonpositive values for @p tolerance will respect the matrix topology as is.
    /// - Defaults to @p 0.
    TValue tolerance;

    /// @brief Verbosity level controlling the volume of output to @p std::cerr and @p std::cout.
    /// @details Each level is a strict superset of lower levels.
    ///          - @p 0: no output is written to @p std::cerr or @p std::cout.
    ///          - @p 1: error messages are written to @p std::cerr.
    ///          - @p 2: information outside iterations is written to @p std::cout.
    ///          - @p 3: information from iterations is written to @p std::cout.
    ///
    /// @note Defaults to @p 1.
    unsigned short verbosity;

    /// @brief Constructor initializing settings to their default values.
    ColorSettings() noexcept
        : shrinkingFactor(-1),
          maxStallCount(1e3),
          tolerance(0),
          verbosity(1)
    {}
};


/// @brief Compute an approximate coloring of a graph.
///
/// @tparam TIndex Integer type of stored indices in the matrix (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries in the matrix (for now, only @p double is supported).
/// @tparam TColor Integer type of vertex colors (@p unsigned or <c>unsigned long</c>).
///
/// @details @par Output
///          The output is written to @p pColors, which stores the color of each row at the matching index.
///
/// @details @par Algorithm
///          The algorithm randomly assigns colors to each uncolored vertex from their palette of remaining
///          colors. After each round of assignment, vertices with valid colors are frozen and a basic
///          conflict resolution step is carried out. The color of newly frozen vertices is removed from their
///          neighbors' palettes.
///
///          The quality of the resulting coloring, as well as the number of iterations to find it strongly
///          depend on the size of the initial palettes. This can be adjusted via the
///          @ref ColorSettings::shrinkingFactor "shrinking factor", which is the only algorithmic parameter.
///
/// @details @par Parallelism
///          Coloring is performed in parallel on the CPU in a shared memory model using @a OpenMP, if
///          @a MCGS was compiled with @a OpenMP support. The maximum number of used threads is controlled
///          by @a OpenMP.
///
/// @details @par Stalling
///          Due to the random nature of assignments, iterations that fail to freeze the color of any remaining
///          vertex can occur. Such an iteration is labelled a stalled iteration, and the function is set to
///          fail after a given number of consequent stalls. The length of tolerated consequent stalls is
///          controlled by @ref ColorSettings::maxStallCount.
///
/// @details @par Errors
///          The algorithm can fail if
///          - the number of consequent stalls reaches @ref ColorSettings::maxStallCount
///          - the number of rows, columns, or nonzeros in the input matrix is negative
///          - the input matrix is not square
///          - @ref CSRAdaptor::pRowExtents, @ref CSRAdaptor::pColumnIndices, or @ref CSRAdaptor::pEntries is @p nullptr.
///          - @p pColors is @p nullptr
///          If any of the above is detected, the function returns @ref MCGS_FAILURE.
///
///          No exceptions are thrown by the function directly, but exceptions from the C++ standard library
///          can propagate to the user (provided that exceptions are enabled by the compiler).
///
///          Segmentation faults can occur if the user manages the memory of the input matrix
///          or output colors incorrectly. Requirements on matrix storage are specified by @ref CSRAdaptor.
///          @p pColors must be at least of length @ref CSRAdaptor::rowCount.
///
/// @param pColors output array of vertex colors. The memory must be managed by the user, and the array
///                must have a length at least as the number of columns in the input matrix @p rMatrix.
/// @param rMatrix @ref CSRAdaptor "CSR matrix adaptor" representing the graph to be colored.
/// @param settings @ref ColorSettings "algorithmic parameters and other settings" to apply during coloring.
///
/// @return @ref MCGS_SUCCESS if successful, otherwise @ref MCGS_FAILURE.
///
/// @note This function is a loose implementation of the algorithm described in
///       <i>Fast Distributed Algorithms for Brooks-Vizing Colorings</i>, doi:<c>10.1006/jagm.2000.1097</c>.
template <class TIndex, class TValue, class TColor>
MCGS_EXPORT_SYMBOL
int color(TColor* pColors,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const ColorSettings<TValue> settings);


/// @brief Partition of a graph with respect to a coloring.
/// @details This class is meant for internal use. It stores
///          an ordering of the representing matrix' rows, as
///          well as the extents of each group within the partition.
///
///          Users can create partitions with @ref makePartition and
///          @b must destroy them with @ref destroyPartition.
///
/// @see makePartition
/// @see destroyPartition
template <class TIndex>
class Partition;


/// @brief Construct the @ref Partition of a graph with respect to a coloring.
///
/// @tparam TIndex Integer type of stored indices in the matrix (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TColor Integer type of vertex colors (@p unsigned or <c>unsigned long</c>).
///
/// @param pColors Array of vertex colors.
/// @param rowCount Number of vertices (must match the number of rows of the matrix representation of the graph).
///
/// @return a valid pointer to a partition if successful, otherwise @p nullptr.
///
/// @note
/// - successful partitionings @b must be destroyed explicitly by the user with @ref destroyPartition.
/// - if partitioning fails, a @p nullptr is returned. In such cases, the user must not invoke
///   @ref destroyPartition with the returned pointer.
///
/// @see destroyPartition
template <class TIndex, class TColor>
MCGS_EXPORT_SYMBOL
[[nodiscard]] Partition<TIndex>* makePartition(const TColor* pColors,
                                               const TIndex rowCount);


/// @brief Destroy a partition that was successfully constructed by @ref makePartition.
///
/// @tparam TIndex Integer type of stored indices in the matrix (@p int, @p long, @p unsigned or <c>unsigned long</c>).
///
/// @param pPartition Pointer to the @ref Partition to destroy.
///
/// @note successful calls to @ref makePartition must be followed by a matching call to @p destroyPartition,
///       once the constructed @ref Partition is no longer used. Memory allocated by the partition will otherwise
///       not be released.
///
/// @see makePartition
template <class TIndex>
MCGS_EXPORT_SYMBOL
void destroyPartition(Partition<TIndex>* pPartition);


/// @brief Reorder rows and columns of a CSR matrix as well as a matching dense vector, with respect to a coloring.
///
/// @tparam TIndex Integer type of stored indices in the matrix (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries in the matrix (for now, only @p double is supported).
///
/// @details This function reorders the rows of a matrix stored in the compressed sparse row format such that
///          rows with identical colors become contiguous. Columns of the matrix, as well as the associated
///          right hand side vector @p pRHS are reordered accordingly.
///
///          A new partition is returned in which rows of the same color are contiguous. Such a partition is
///          a prerequisite for running @ref solve(TValue*,const CSRAdaptor&,const TValue*,const Partition*,const SolveSettings)
///          "parallel versions of Gauss-Seidel smoothing".
///
/// @param rowCount Number of rows in the matrix.
/// @param columnCount Number of columns in the matrix.
/// @param entryCount Number of entries stored in the compressed matrix.
/// @param pRowExtents Index array defining entries in the matrix related to each row in @p pEntries.
///                    Size must be at least <c>rowCount + 1</c>.
/// @param pColumnIndices Index array assigning column indices to each entry in @p pEntries.
///                       Size must be at least @p entryCount.
/// @param pEntries Array of stored entries in the matrix. Size must be at least @p entryCount.
/// @param pSolution Initial solution vector stored in a dense contiguous array. Size must be at least @p columnCount.
///                  Pass @p nullptr if you don't want to reorder the initial solution vector.
/// @param pRHS Right hand side vector stored in a dense contiguous array. Size must be at least @p rowCount.
///             Pass @p nullptr if you don't want to reorder the right hand side vector.
/// @param pPartition Pointer to the partition that encodes the coloring of the input matrix.
///
/// @return a pointer to a contiguous partition of the reordered system if the reordering is successful,
///         otherwise @p nullptr.
///
/// @note
/// - The function allocates a duplicate of the input matrix, vector and partition during execution.
/// - The returned partition must be destroyed explicitly by the user with @ref destroyPartition
///   (provided the reordering was successful).
///
/// @see revertReorder(TValue*,const TIndex,const Partition*)
/// @see revertReorder(const TIndex,const TIndex,const TIndex,TIndex*,TIndex*,TValue*,TValue*,TValue*,const Partition*)
template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
[[nodiscard]] Partition<TIndex>* reorder(const unsigned long rowCount, const unsigned long columnCount, const unsigned long entryCount,
                                         TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pEntries,
                                         TValue* pSolution, TValue* pRHS,
                                         const Partition<TIndex>* pPartition);


/// @brief Restore the original order of a @ref reorder "reordered" vector.
///
/// @tparam TIndex Integer type of the partition (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries in the vector (for now, only @p double is supported).
///
/// @param pRHS Dense vector to be reordered.
/// @param rowCount Size of the dense vector (should be equal to the number of columns in the associated matrix).
/// @param pPartition Pointer to the original @ref Partition the vector was reordered by.
///
/// @return @ref MCGS_SUCCESS if successful, otherwise @ref MCGS_FAILURE.
///
/// @fn revertReorder(TValue*,const unsigned long,const Partition*)
template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int revertReorder(TValue* pRHS,
                  const unsigned long rowCount,
                  const Partition<TIndex>* pPartition);


/// @brief Restore the original order of a @ref reorder "reordered" CSR matrix and associated right hand side vector.
///
/// @tparam TIndex Integer type of stored indices in the matrix (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries in the matrix (for now, only @p double is supported).
///
/// @param rowCount Number of rows in the matrix.
/// @param columnCount Number of columns in the matrix.
/// @param entryCount Number of entries stored in the compressed matrix.
/// @param pRowExtents Index array defining entries in the matrix related to each row in @p pEntries.
///                    Size must be at least <c>rowCount + 1</c>.
/// @param pColumnIndices Index array assigning column indices to each entry in @p pEntries.
///                       Size must be at least @p entryCount.
/// @param pEntries Array of stored entries in the matrix. Size must be at least @p entryCount.
/// @param pSolution Initial solution vector stored in a dense contiguous array. Size must be at least @p columnCount.
///                  Pass @p nullptr if you don't want to reorder the initial solution vector.
/// @param pRHS Right hand side vector stored in a dense contiguous array. Size must be at least @p rowCount.
///             Pass @p nullptr if you don't want to reorder the right hand side vector.
/// @param pPartition Pointer to the original @ref Partition the matrix and vector were reordered by.
///
/// @return @ref MCGS_SUCCESS if successful, otherwise @ref MCGS_FAILURE.
///
/// @fn revertReorder(const unsigned long,const unsigned long,const unsigned long,TIndex*,TIndex*,TValue*,TValue*,TValue*,const Partition*)
template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int revertReorder(const unsigned long rowCount, const unsigned long columnCount, const unsigned long entryCount,
                  TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pEntries,
                  TValue* pSolution, TValue* pRHS,
                  const Partition<TIndex>* pPartition);


/// @brief Compute the 2-norm of the residual of a linear system's approximate solution.
///
/// @tparam TIndex Integer type of stored indices in the matrix (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries in the matrix (for now, only @p double is supported).
///
/// @details @f[ r = \sqrt{\sum_i{(b_i - a_{ij}x_j)^2}} @f]
///
/// @param rMatrix Left hand side square matrix in compressed sparse row format (@f$a_{ij}@f$).
/// @param pSolution Approximate solution vector (@f$x_j@f$).
///                  Size must be at least equal to the number of @ref CSRAdaptor::columnCount "columns" in the matrix.
/// @param pRHS Right hand side vector (@f$b_i@f$).
///             Size must be at least equal to the number of @ref CSRAdaptor::rowCount "rows" in the matrix.
///
/// @return The residual's 2-norm.
template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
TValue residual(const CSRAdaptor<TIndex,TValue>& rMatrix,
                const TValue* pSolution,
                const TValue* pRHS) noexcept;


/// @brief Enum for all supported Gauss-Seidel parallelization strategies.
enum struct Parallelization
{
    None        = 1,    ///< @brief Perform Gauss-Seidel iterations in serial.
    RowWise     = 2,    ///< @brief Distribute work assuming each row has the same number of entries.
    EntryWise   = 4     ///< @brief Distribute work along equal chunks of entries.
};


/// @brief Settings for Gauss-Seidel relaxation.
///
/// @tparam TIndex Integer type of stored indices in the matrix (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries in the matrix (for now, only @p double is supported).
template <class TIndex, class TValue>
struct SolveSettings
{
    /// @brief Absolute tolerance used for convergence checks.
    /// @note
    /// - If @ref residualRelativeTolerance "both" tolerances are non-positive,
    ///   residuals are not computed and relaxation ends when the number of iterations
    ///   reaches @ref maxIterations.
    /// - Defaults to @p -1.
    TValue residualAbsoluteTolerance;

    /// @brief Relative tolerance used for convergence checks.
    /// @note
    /// - If @ref residualAbsoluteTolerance "both" tolerances are non-positive,
    ///   residuals are not computed and relaxation ends when the number of iterations
    ///   reaches @ref maxIterations.
    /// - Defaults to @p -1.
    TValue residualRelativeTolerance;

    /// @brief Maximum number of Gauss-Seidel iterations to allow.
    /// @note
    /// - Defaults to @p 1.
    TIndex maxIterations;

    /// @brief Overrelaxation parameter.
    /// @note Defaults to @p 1.0.
    /// @see <a href="https://en.wikipedia.org/wiki/Successive_over-relaxation">Successive Over-relaxation</a>
    TValue relaxation;

    /// @brief Parallelization strategy.
    /// @note Defaults to @ref Parallelization::RowWise "row-wise" parallelization.
    /// @see Parallelization
    Parallelization parallelization;

    /// @brief Verbosity level controlling the volume of output to @p std::cerr and @p std::cout.
    /// @details Each level is a strict superset of lower levels.
    ///          - @p 0: no output is written to @p std::cerr or @p std::cout.
    ///          - @p 1: error messages are written to @p std::cerr.
    ///          - @p 2: information outside iterations is written to @p std::cout.
    ///          - @p 3: information from iterations is written to @p std::cout.
    ///
    /// @note Defaults to @p 1.
    unsigned short verbosity;

    /// @brief Constructor initializing settings to their default values.
    SolveSettings() noexcept
        : residualAbsoluteTolerance(-1),
          residualRelativeTolerance(-1),
          maxIterations(1),
          relaxation(1),
          parallelization(Parallelization::RowWise),
          verbosity(1)
    {}
};


/// @brief Perform Gauss-Seidel relaxation in serial.
///
/// @tparam TIndex Integer type of stored indices (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries (for now, only @p double is supported).
///
/// @details @par Algorithm
///          Given the linear system @f$ a_{ij} x_j^0 = b_i @f$ and the decomposition of
///          @f$a_{ij}@f$ into the sum of a strictly lower triangular matrix @f$l_{ij}@f$,
///          a diagonal matrix @f$d_{ij}@f$ and strictly upper triangular matrix @f$u_{ij}@f$
///          (@f$ a_{ij} = l_{ij} + d_{ij} + u_{ij} @f$), this function computes
///          @f[
///             x_j^{k+1} = (1 - \omega)x_j^k + \frac{\omega}{a_{jj}} (b_j - \sum_{i<j}{a_{ji} x_i^{k+1}} - \sum_{i>j}{a_{ji} x_i^k})
///          @f]
///          where @f$\omega@f$ is a user-selected relaxation parameter.
///
/// @param pSolution Solution array @f$x_j@f$. Size must be at least equal to the number of @ref CSRAdaptor::columnCount "columns" in @p rMatrix.
/// @param rMatrix Left hand side square matrix in compressed sparse row format @f$a_{ij}@f$.
/// @param pRHS Right hand side vector @f$b_i@f$. Size must be at least equal to the number of @ref CSRAdaptor::rowCount "rows" in @p rMatrix.
/// @param settings @ref SolveSettings "algorithmic parameters and other settings" to apply during relaxation.
///
/// @return @ref MCGS_SUCCESS if successful, otherwise @ref MCGS_FAILURE.
///
/// @see solve(TValue*,const CSRAdaptor&,const TValue*,const Partition*,const SolveSettings)
/// @fn solve(TValue*,const CSRAdaptor&,const TValue*,const SolveSettings)
template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const SolveSettings<TIndex,TValue> settings);


/// @brief Perform Gauss-Seidel relaxation on a reordered system in parallel.
///
/// @tparam TIndex Integer type of stored indices (@p int, @p long, @p unsigned or <c>unsigned long</c>).
/// @tparam TValue Number type of stored entries (for now, only @p double is supported).
///
/// @details @par Algorithm
///          Given the linear system @f$ a_{ij} x_j^0 = b_i @f$ and the decomposition of
///          @f$a_{ij}@f$ into the sum of a strictly lower triangular matrix @f$l_{ij}@f$,
///          a diagonal matrix @f$d_{ij}@f$ and strictly upper triangular matrix @f$u_{ij}@f$
///          (@f$ a_{ij} = l_{ij} + d_{ij} + u_{ij} @f$), this function computes
///          @f[
///             x_j^{k+1} = (1 - \omega)x_j^k + \frac{\omega}{a_{jj}} (b_j - \sum_{i<j}{a_{ji} x_i^{k+1}} - \sum_{i<j}{a_{ji} x_i^k})
///          @f]
///          where @f$\omega@f$ is a user-selected relaxation parameter.
///
/// @details @par Prerequisites
///          The input system must be reordered with respect to a coloring of the system matrix' graph.
///          This ensures that the algorithm doesn't run into race conditions and performs Gauss-Seidel
///          iterations instead of some aliased product.
///
/// @param pSolution Solution array @f$x_j@f$. Size must be at least equal to the number of @ref CSRAdaptor::columnCount "columns" in @p rMatrix.
/// @param rMatrix Left hand side square matrix in compressed sparse row format @f$a_{ij}@f$.
/// @param pRHS Right hand side vector @f$b_i@f$. Size must be at least equal to the number of @ref CSRAdaptor::rowCount "rows" in @p rMatrix.
/// @param pPartition Contiguous partition of the input matrix.
/// @param settings @ref SolveSettings "algorithmic parameters and other settings" to apply during relaxation.
///
/// @return @ref MCGS_SUCCESS if successful, otherwise @ref MCGS_FAILURE.
///
/// @see solve(TValue*,const CSRAdaptor&,const TValue*,const SolveSettings)
/// @see color
/// @see reorder
/// @fn solve(TValue*,const CSRAdaptor&,const TValue*,const Partition*,const SolveSettings)
template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const Partition<TIndex>* pPartition,
          const SolveSettings<TIndex,TValue> settings);


} // namespace mcgs


#undef MCGS_EXPORT_SYMBOL
