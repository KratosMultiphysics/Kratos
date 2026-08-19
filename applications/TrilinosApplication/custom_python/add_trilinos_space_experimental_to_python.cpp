//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "trilinos_space_experimental.h"
#include "custom_python/add_trilinos_space_experimental_to_python.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos::Python
{
namespace py = pybind11;

#ifdef HAVE_TPETRA

using ExperimentalTrilinosSparseSpaceType = TrilinosSpaceExperimental<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;

/**
 * @brief Thin RAII wrapper around a Tpetra FECrsMatrix shared pointer.
 * @details Pybind11 cannot bind Teuchos::RCP directly, so this wrapper
 *          exposes GetPointer() and GetReference() to Python.
 */
class ExperimentalAuxiliaryMatrixWrapper
{
public:
    using TrilinosMatrixType = typename ExperimentalTrilinosSparseSpaceType::MatrixType;
    using TrilinosMatrixPointerType = typename ExperimentalTrilinosSparseSpaceType::MatrixPointerType;

    ExperimentalAuxiliaryMatrixWrapper(TrilinosMatrixPointerType pMatrix) : mpMatrix(pMatrix) {}
    virtual ~ExperimentalAuxiliaryMatrixWrapper() = default;
    TrilinosMatrixPointerType& GetPointer() { return mpMatrix; }
    TrilinosMatrixType& GetReference() { return *mpMatrix; }

private:
    TrilinosMatrixPointerType mpMatrix;
};

/**
 * @brief Thin RAII wrapper around a Tpetra FEMultiVector shared pointer.
 * @details Pybind11 cannot bind Teuchos::RCP directly, so this wrapper
 *          exposes GetPointer() and GetReference() to Python.
 */
class ExperimentalAuxiliaryVectorWrapper
{
public:
    using TrilinosVectorType = typename ExperimentalTrilinosSparseSpaceType::VectorType;
    using TrilinosVectorPointerType = typename ExperimentalTrilinosSparseSpaceType::VectorPointerType;

    ExperimentalAuxiliaryVectorWrapper(TrilinosVectorPointerType pVector) : mpVector(pVector) {}
    virtual ~ExperimentalAuxiliaryVectorWrapper() = default;
    TrilinosVectorPointerType& GetPointer() { return mpVector; }
    TrilinosVectorType& GetReference() { return *mpVector; }

private:
    TrilinosVectorPointerType mpVector;
};

namespace // anonymous: internal linkage for pybind11 helper trampolines
{

// --- Wrapper accessor helpers ---

ExperimentalTrilinosSparseSpaceType::MatrixType& GetMatRef(ExperimentalAuxiliaryMatrixWrapper& rWrapper)
{
    return rWrapper.GetReference();
}

ExperimentalTrilinosSparseSpaceType::VectorType& GetVecRef(ExperimentalAuxiliaryVectorWrapper& rWrapper)
{
    return rWrapper.GetReference();
}

// --- Pointer factories ---

ExperimentalAuxiliaryMatrixWrapper CreateEmptyMatrixPointer(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    Teuchos::MpiComm<int>& rComm
    )
{
    return ExperimentalAuxiliaryMatrixWrapper(rDummy.CreateEmptyMatrixPointer(rComm));
}

ExperimentalAuxiliaryVectorWrapper CreateEmptyVectorPointer(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    Teuchos::MpiComm<int>& rComm
    )
{
    return ExperimentalAuxiliaryVectorWrapper(rDummy.CreateEmptyVectorPointer(rComm));
}

ExperimentalAuxiliaryVectorWrapper CreateVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MapPointerType& pMap
    )
{
    return ExperimentalAuxiliaryVectorWrapper(rDummy.CreateVector(pMap));
}

ExperimentalAuxiliaryMatrixWrapper CreateMatrix(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::GraphPointerType& pGraph
    )
{
    return ExperimentalAuxiliaryMatrixWrapper(rDummy.CreateMatrix(pGraph));
}

ExperimentalAuxiliaryVectorWrapper CreateVectorCopy(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::VectorType& rV
    )
{
    return ExperimentalAuxiliaryVectorWrapper(rDummy.CreateVectorCopy(rV));
}

ExperimentalAuxiliaryMatrixWrapper CreateMatrixCopy(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rA
    )
{
    return ExperimentalAuxiliaryMatrixWrapper(rDummy.CreateMatrixCopy(rA));
}

// --- Clear / resize ---

void ClearMatrix(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixPointerType& pA
    )
{
    rDummy.Clear(pA);
}

void ClearVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorPointerType& pX
    )
{
    rDummy.Clear(pX);
}

void ResizeVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorPointerType& pX,
    ExperimentalTrilinosSparseSpaceType::SizeType N
    )
{
    rDummy.Resize(pX, N);
}

// --- Zero-fill ---

void SetToZeroMatrix(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rA
    )
{
    rDummy.SetToZero(rA);
}

void SetToZeroVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX
    )
{
    rDummy.SetToZero(rX);
}

// --- Norms ---

double TwoNorm(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX
    )
{
    return rDummy.TwoNorm(rX);
}

double GetDiagonalNorm(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rA
    )
{
    return rDummy.GetDiagonalNorm(rA);
}

double GetAveragevalueDiagonal(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rA
    )
{
    return rDummy.GetAveragevalueDiagonal(rA);
}

double GetMaxDiagonal(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rA
    )
{
    return rDummy.GetMaxDiagonal(rA);
}

double GetMinDiagonal(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rA
    )
{
    return rDummy.GetMinDiagonal(rA);
}

// --- Dot product ---

double Dot(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    ExperimentalTrilinosSparseSpaceType::VectorType& rY
    )
{
    return rDummy.Dot(rX, rY);
}

// --- Vector / matrix algebra ---

void ScaleAndAdd(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const double A,
    const ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    const double B,
    ExperimentalTrilinosSparseSpaceType::VectorType& rY
    )
{
    rDummy.ScaleAndAdd(A, rX, B, rY);
}

void ScaleAndAddMatrix(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const double A,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rX,
    const double B,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rY
    )
{
    rDummy.ScaleAndAdd(A, rX, B, rY);
}

void UnaliasedAdd(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    const double A,
    const ExperimentalTrilinosSparseSpaceType::VectorType& rY
    )
{
    rDummy.UnaliasedAdd(rX, A, rY);
}

void Assign(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    const double A,
    const ExperimentalTrilinosSparseSpaceType::VectorType& rY
    )
{
    rDummy.Assign(rX, A, rY);
}

void InplaceMult(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    const double A
    )
{
    rDummy.InplaceMult(rX, A);
}

void SetValue(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    const double Value
    )
{
    rDummy.Set(rX, Value);
}

// --- Matrix-vector / matrix-matrix products ---

void Mult(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rA,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    ExperimentalTrilinosSparseSpaceType::VectorType& rY
    )
{
    rDummy.Mult(rA, rX, rY);
}

void TransposeMult(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rA,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    ExperimentalTrilinosSparseSpaceType::VectorType& rY
    )
{
    rDummy.TransposeMult(rA, rX, rY);
}

// --- Size queries ---

ExperimentalTrilinosSparseSpaceType::IndexType Size(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType const& rV
    )
{
    return rDummy.Size(rV);
}

ExperimentalTrilinosSparseSpaceType::IndexType Size1(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType const& rM
    )
{
    return rDummy.Size1(rM);
}

ExperimentalTrilinosSparseSpaceType::IndexType Size2(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType const& rM
    )
{
    return rDummy.Size2(rM);
}

// --- Copy ---

void CopyVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    ExperimentalTrilinosSparseSpaceType::VectorType& rY
    )
{
    rDummy.Copy(rX, rY);
}

void CopyMatrix(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rX,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rY
    )
{
    rDummy.Copy(rX, rY);
}

void CopyMatrixValues(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rA,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rB
    )
{
    rDummy.CopyMatrixValues(rA, rB);
}

// --- Element access ---

double GetValue(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    std::size_t I
    )
{
    return rDummy.GetValue(rX, I);
}

void SetValueVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    ExperimentalTrilinosSparseSpaceType::IndexType I,
    const double Value
    )
{
    rDummy.SetValue(rX, I, Value);
}

void SetValueMatrix(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rA,
    ExperimentalTrilinosSparseSpaceType::IndexType I,
    ExperimentalTrilinosSparseSpaceType::IndexType J,
    const double Value
    )
{
    rDummy.SetValue(rA, I, J, Value);
}

py::object GatherValues(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rX,
    const std::vector<int>& rIndexArray
    )
{
    std::vector<double> values(rIndexArray.size());
    rDummy.GatherValues(rX, rIndexArray, values.data());
    return py::cast(values);
}

// --- Assembly ---

void GlobalAssembleMatrix(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rA
    )
{
    rDummy.GlobalAssemble(rA);
}

void GlobalAssembleVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::VectorType& rV
    )
{
    rDummy.GlobalAssemble(rV);
}

void ManualFinalize(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType& rA
    )
{
    rDummy.ManualFinalize(rA);
}

// --- Matrix Market I/O ---

ExperimentalAuxiliaryMatrixWrapper ReadMatrixMarket(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const std::string& rFileName,
    Teuchos::MpiComm<int>& rComm
    )
{
    return ExperimentalAuxiliaryMatrixWrapper(rDummy.ReadMatrixMarket(rFileName, rComm));
}

ExperimentalAuxiliaryVectorWrapper ReadMatrixMarketVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const std::string& rFileName,
    ExperimentalTrilinosSparseSpaceType::CommunicatorPointerType pComm,
    const int N
    )
{
    return ExperimentalAuxiliaryVectorWrapper(rDummy.ReadMatrixMarketVector(rFileName, pComm, N));
}

void WriteMatrixMarketMatrix(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const char* pFileName,
    const ExperimentalTrilinosSparseSpaceType::MatrixType& rA,
    const bool Symmetric
    )
{
    rDummy.WriteMatrixMarketMatrix(pFileName, rA, Symmetric);
}

void WriteMatrixMarketVector(
    ExperimentalTrilinosSparseSpaceType& rDummy,
    const char* pFileName,
    const ExperimentalTrilinosSparseSpaceType::VectorType& rV
    )
{
    rDummy.WriteMatrixMarketVector(pFileName, rV);
}

} // anonymous namespace

#endif

void AddBasicOperationsExperimental(py::module& m)
{
#ifdef HAVE_TPETRA

    // Expose underlying Tpetra/Teuchos types so Python can receive or pass them
    py::class_<Teuchos::MpiComm<int>>(m, "Experimental_TeuchosMpiComm");
    py::class_<Tpetra::FECrsMatrix<>>(m, "Experimental_FECrsMatrix");
    py::class_<Tpetra::FEMultiVector<>>(m, "Experimental_FEMultiVector");

    // Smart-pointer wrappers (Teuchos::RCP cannot be bound directly)
    py::class_<ExperimentalAuxiliaryMatrixWrapper>(m, "ExperimentalTrilinosMatrixPointer")
        .def("GetReference", GetMatRef, py::return_value_policy::reference_internal)
        ;

    py::class_<ExperimentalAuxiliaryVectorWrapper>(m, "ExperimentalTrilinosVectorPointer")
        .def("GetReference", GetVecRef, py::return_value_policy::reference_internal)
        ;

    py::class_<ExperimentalTrilinosSparseSpaceType>(m, "ExperimentalTrilinosSparseSpace")
        .def(py::init<>())
        // --- Pointer factories ---
        .def("CreateEmptyMatrixPointer", CreateEmptyMatrixPointer)
        .def("CreateEmptyVectorPointer", CreateEmptyVectorPointer)
        .def("CreateVector",             CreateVector)
        .def("CreateMatrix",             CreateMatrix)
        .def("CreateVectorCopy",         CreateVectorCopy)
        .def("CreateMatrixCopy",         CreateMatrixCopy)
        // --- Clear / resize ---
        .def("ClearMatrix",  ClearMatrix)
        .def("ClearVector",  ClearVector)
        .def("ResizeVector", ResizeVector)
        // --- Zero-fill ---
        .def("SetToZeroMatrix", SetToZeroMatrix)
        .def("SetToZeroVector", SetToZeroVector)
        // --- Norms ---
        .def("TwoNorm",                 TwoNorm)
        .def("GetDiagonalNorm",         GetDiagonalNorm)
        .def("GetAveragevalueDiagonal", GetAveragevalueDiagonal)
        .def("GetMaxDiagonal",          GetMaxDiagonal)
        .def("GetMinDiagonal",          GetMinDiagonal)
        // --- Dot product ---
        .def("Dot", Dot)
        // --- Vector / matrix algebra ---
        .def("UnaliasedAdd",      UnaliasedAdd)
        .def("ScaleAndAdd",       ScaleAndAdd)
        .def("ScaleAndAddMatrix", ScaleAndAddMatrix)
        .def("Assign",            Assign)
        .def("InplaceMult",       InplaceMult)
        .def("Set",               SetValue)
        // --- Matrix-vector / matrix-matrix products ---
        .def("Mult",          Mult)
        .def("TransposeMult", TransposeMult)
        // --- Size queries ---
        .def("Size",  Size)
        .def("Size1", Size1)
        .def("Size2", Size2)
        // --- Copy ---
        .def("CopyVector",       CopyVector)
        .def("CopyMatrix",       CopyMatrix)
        .def("CopyMatrixValues", CopyMatrixValues)
        // --- Element access ---
        .def("GetValue",       GetValue)
        .def("SetValueVector", SetValueVector)
        .def("SetValueMatrix", SetValueMatrix)
        .def("GatherValues",   GatherValues)
        // --- Assembly ---
        .def("GlobalAssembleMatrix", GlobalAssembleMatrix)
        .def("GlobalAssembleVector", GlobalAssembleVector)
        .def("ManualFinalize",       ManualFinalize)
        // --- Matrix Market I/O ---
        .def("ReadMatrixMarket",        ReadMatrixMarket)
        .def("ReadMatrixMarketVector",  ReadMatrixMarketVector)
        .def("WriteMatrixMarketMatrix", WriteMatrixMarketMatrix)
        .def("WriteMatrixMarketVector", WriteMatrixMarketVector)
        // --- Static queries ---
        .def_static("IsDistributed",           &ExperimentalTrilinosSparseSpaceType::IsDistributed)
        .def_static("IsDistributedSpace",      &ExperimentalTrilinosSparseSpaceType::IsDistributedSpace)
        .def_static("FastestDirectSolverList", &ExperimentalTrilinosSparseSpaceType::FastestDirectSolverList)
        ;

#endif
}

} // namespace Kratos::Python