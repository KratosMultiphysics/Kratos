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

namespace Kratos::Python {
namespace py = pybind11;

#ifdef HAVE_TPETRA

using ExperimentalTrilinosSparseSpaceType = TrilinosSpaceExperimental<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;

class ExperimentalAuxiliaryMatrixWrapper {
public:
  typedef typename ExperimentalTrilinosSparseSpaceType::MatrixType
      TrilinosMatrixType;
  typedef typename ExperimentalTrilinosSparseSpaceType::MatrixPointerType
      TrilinosMatrixPointerType;

  ExperimentalAuxiliaryMatrixWrapper(TrilinosMatrixPointerType p) : mp(p) {};
  virtual ~ExperimentalAuxiliaryMatrixWrapper() {}
  TrilinosMatrixPointerType &GetPointer() { return mp; }
  TrilinosMatrixType &GetReference() { return *mp; }

private:
  TrilinosMatrixPointerType mp;
};

class ExperimentalAuxiliaryVectorWrapper {
public:
  typedef typename ExperimentalTrilinosSparseSpaceType::VectorType
      TrilinosVectorType;
  typedef typename ExperimentalTrilinosSparseSpaceType::VectorPointerType
      TrilinosVectorPointerType;

  ExperimentalAuxiliaryVectorWrapper(TrilinosVectorPointerType p) : mp(p) {};
  virtual ~ExperimentalAuxiliaryVectorWrapper() {}
  TrilinosVectorPointerType &GetPointer() { return mp; }
  TrilinosVectorType &GetReference() { return *mp; }

private:
  TrilinosVectorPointerType mp;
};

namespace {

double ExperimentalDot(ExperimentalTrilinosSparseSpaceType &dummy,
                       ExperimentalTrilinosSparseSpaceType::VectorType &rX,
                       ExperimentalTrilinosSparseSpaceType::VectorType &rY) {
  return dummy.Dot(rX, rY);
}
void ExperimentalScaleAndAdd(
    ExperimentalTrilinosSparseSpaceType &dummy, const double A,
    const ExperimentalTrilinosSparseSpaceType::VectorType &rX, const double B,
    ExperimentalTrilinosSparseSpaceType::VectorType &rY) {
  dummy.ScaleAndAdd(A, rX, B, rY);
}
void ExperimentalScaleAndAddMatrix(
    ExperimentalTrilinosSparseSpaceType &dummy, const double A,
    const ExperimentalTrilinosSparseSpaceType::MatrixType &rX, const double B,
    ExperimentalTrilinosSparseSpaceType::MatrixType &rY) {
  dummy.ScaleAndAdd(A, rX, B, rY);
}
void ExperimentalMult(ExperimentalTrilinosSparseSpaceType &dummy,
                      ExperimentalTrilinosSparseSpaceType::MatrixType &rA,
                      ExperimentalTrilinosSparseSpaceType::VectorType &rX,
                      ExperimentalTrilinosSparseSpaceType::VectorType &rY) {
  dummy.Mult(rA, rX, rY);
}
void ExperimentalTransposeMult(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType &rA,
    ExperimentalTrilinosSparseSpaceType::VectorType &rX,
    ExperimentalTrilinosSparseSpaceType::VectorType &rY) {
  dummy.TransposeMult(rA, rX, rY);
}
ExperimentalTrilinosSparseSpaceType::IndexType
ExperimentalSize(ExperimentalTrilinosSparseSpaceType &dummy,
                 ExperimentalTrilinosSparseSpaceType::VectorType const &rV) {
  return dummy.Size(rV);
}
ExperimentalTrilinosSparseSpaceType::IndexType
ExperimentalSize1(ExperimentalTrilinosSparseSpaceType &dummy,
                  ExperimentalTrilinosSparseSpaceType::MatrixType const &rM) {
  return dummy.Size1(rM);
}
ExperimentalTrilinosSparseSpaceType::IndexType
ExperimentalSize2(ExperimentalTrilinosSparseSpaceType &dummy,
                  ExperimentalTrilinosSparseSpaceType::MatrixType const &rM) {
  return dummy.Size2(rM);
}
void ExperimentalClearMatrix(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::MatrixPointerType &pA) {
  dummy.Clear(pA);
}
void ExperimentalClearVector(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorPointerType &pX) {
  dummy.Clear(pX);
}
void ExperimentalResizeVector(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorPointerType &pX,
    ExperimentalTrilinosSparseSpaceType::SizeType n) {
  dummy.Resize(pX, n);
}
void ExperimentalSetToZeroMatrix(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType &A) {
  dummy.SetToZero(A);
}
void ExperimentalSetToZeroVector(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorType &x) {
  dummy.SetToZero(x);
}
double ExperimentalTwoNorm(ExperimentalTrilinosSparseSpaceType &dummy,
                           ExperimentalTrilinosSparseSpaceType::VectorType &x) {
  return dummy.TwoNorm(x);
}
void ExperimentalUnaliasedAdd(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorType &x, const double A,
    const ExperimentalTrilinosSparseSpaceType::VectorType &y) {
  dummy.UnaliasedAdd(x, A, y);
}
ExperimentalAuxiliaryMatrixWrapper
ExperimentalCreateEmptyMatrixPointer(ExperimentalTrilinosSparseSpaceType &dummy,
                                     Teuchos::MpiComm<int> &rComm) {
  return ExperimentalAuxiliaryMatrixWrapper(
      dummy.CreateEmptyMatrixPointer(rComm));
}
ExperimentalAuxiliaryVectorWrapper
ExperimentalCreateEmptyVectorPointer(ExperimentalTrilinosSparseSpaceType &dummy,
                                     Teuchos::MpiComm<int> &rComm) {
  return ExperimentalAuxiliaryVectorWrapper(
      dummy.CreateEmptyVectorPointer(rComm));
}
ExperimentalTrilinosSparseSpaceType::MatrixType &
ExperimentalGetMatRef(ExperimentalAuxiliaryMatrixWrapper &dummy) {
  return dummy.GetReference();
}
ExperimentalTrilinosSparseSpaceType::VectorType &
ExperimentalGetVecRef(ExperimentalAuxiliaryVectorWrapper &dummy) {
  return dummy.GetReference();
}

// Additional bindings

void ExperimentalCopyVector(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::VectorType &rX,
    ExperimentalTrilinosSparseSpaceType::VectorType &rY) {
  dummy.Copy(rX, rY);
}
void ExperimentalCopyMatrix(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType &rX,
    ExperimentalTrilinosSparseSpaceType::MatrixType &rY) {
  dummy.Copy(rX, rY);
}
void ExperimentalCopyMatrixValues(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType &rA,
    const ExperimentalTrilinosSparseSpaceType::MatrixType &rB) {
  dummy.CopyMatrixValues(rA, rB);
}
void ExperimentalAssign(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorType &rX, const double A,
    const ExperimentalTrilinosSparseSpaceType::VectorType &rY) {
  dummy.Assign(rX, A, rY);
}
void ExperimentalInplaceMult(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorType &rX, const double A) {
  dummy.InplaceMult(rX, A);
}
void ExperimentalSetV(ExperimentalTrilinosSparseSpaceType &dummy,
                      ExperimentalTrilinosSparseSpaceType::VectorType &rX,
                      double value) {
  dummy.Set(rX, value);
}
double ExperimentalGetValue(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::VectorType &rX, std::size_t I) {
  return dummy.GetValue(rX, I);
}
void ExperimentalSetValueVector(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorType &rX,
    ExperimentalTrilinosSparseSpaceType::IndexType i, double value) {
  dummy.SetValue(rX, i, value);
}
void ExperimentalSetValueMatrix(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType &rA,
    ExperimentalTrilinosSparseSpaceType::IndexType i,
    ExperimentalTrilinosSparseSpaceType::IndexType j, double value) {
  dummy.SetValue(rA, i, j, value);
}
pybind11::object ExperimentalGatherValues(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorType &rX,
    const std::vector<int> &IndexArray) {
  std::vector<double> values(IndexArray.size());
  dummy.GatherValues(rX, IndexArray, values.data());
  return pybind11::cast(values);
}
void ExperimentalGlobalAssembleMatrix(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType &rA) {
  dummy.GlobalAssemble(rA);
}
void ExperimentalGlobalAssembleVector(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::VectorType &rV) {
  dummy.GlobalAssemble(rV);
}
void ExperimentalManualFinalize(
    ExperimentalTrilinosSparseSpaceType &dummy,
    ExperimentalTrilinosSparseSpaceType::MatrixType &rA) {
  dummy.ManualFinalize(rA);
}
double ExperimentalGetDiagonalNorm(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType &rA) {
  return dummy.GetDiagonalNorm(rA);
}
double ExperimentalGetAveragevalueDiagonal(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType &rA) {
  return dummy.GetAveragevalueDiagonal(rA);
}
double ExperimentalGetMaxDiagonal(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType &rA) {
  return dummy.GetMaxDiagonal(rA);
}
double ExperimentalGetMinDiagonal(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType &rA) {
  return dummy.GetMinDiagonal(rA);
}
ExperimentalAuxiliaryVectorWrapper ExperimentalCreateVectorCopy(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::VectorType &rV) {
  return ExperimentalAuxiliaryVectorWrapper(dummy.CreateVectorCopy(rV));
}
ExperimentalAuxiliaryMatrixWrapper ExperimentalCreateMatrixCopy(
    ExperimentalTrilinosSparseSpaceType &dummy,
    const ExperimentalTrilinosSparseSpaceType::MatrixType &rA) {
  return ExperimentalAuxiliaryMatrixWrapper(dummy.CreateMatrixCopy(rA));
}

} // namespace

#endif

void AddBasicOperationsExperimental(pybind11::module &m) {
#ifdef HAVE_TPETRA

  py::class_<Teuchos::MpiComm<int>>(m, "Experimental_TeuchosMpiComm")
      //.def(py::init< Teuchos::MpiComm<int>& >())
      ;

  py::class_<Tpetra::FECrsMatrix<>>(m, "Experimental_FECrsMatrix")
      //.def(py::init< Tpetra::FECrsMatrix<>& >())
      ;

  py::class_<Tpetra::FEMultiVector<>>(m, "Experimental_FEMultiVector")
      //.def(py::init< Tpetra::FEMultiVector<>& >())
      ;

  py::class_<ExperimentalAuxiliaryMatrixWrapper>(
      m, "ExperimentalTrilinosMatrixPointer")
      .def("GetReference", ExperimentalGetMatRef,
           py::return_value_policy::reference_internal);

  py::class_<ExperimentalAuxiliaryVectorWrapper>(
      m, "ExperimentalTrilinosVectorPointer")
      .def("GetReference", ExperimentalGetVecRef,
           py::return_value_policy::reference_internal);

  py::class_<ExperimentalTrilinosSparseSpaceType>(
      m, "ExperimentalTrilinosSparseSpace")
      .def(py::init<>())
      // --- Clear / Resize ---
      .def("ClearMatrix", ExperimentalClearMatrix)
      .def("ClearVector", ExperimentalClearVector)
      .def("ResizeVector", ExperimentalResizeVector)
      // --- Zero-fill ---
      .def("SetToZeroMatrix", ExperimentalSetToZeroMatrix)
      .def("SetToZeroVector", ExperimentalSetToZeroVector)
      // --- Norms ---
      .def("TwoNorm", ExperimentalTwoNorm)
      .def("GetDiagonalNorm", ExperimentalGetDiagonalNorm)
      .def("GetAveragevalueDiagonal", ExperimentalGetAveragevalueDiagonal)
      .def("GetMaxDiagonal", ExperimentalGetMaxDiagonal)
      .def("GetMinDiagonal", ExperimentalGetMinDiagonal)
      // --- Dot / vector algebra ---
      .def("Dot", ExperimentalDot)
      .def("UnaliasedAdd", ExperimentalUnaliasedAdd)
      .def("ScaleAndAdd", ExperimentalScaleAndAdd)
      .def("ScaleAndAddMatrix", ExperimentalScaleAndAddMatrix)
      .def("Assign", ExperimentalAssign)
      .def("InplaceMult", ExperimentalInplaceMult)
      // --- Matrix-vector / matrix-matrix products ---
      .def("Mult", ExperimentalMult)
      .def("TransposeMult", ExperimentalTransposeMult)
      // --- Sizes ---
      .def("Size", ExperimentalSize)
      .def("Size1", ExperimentalSize1)
      .def("Size2", ExperimentalSize2)
      // --- Copy ---
      .def("CopyVector", ExperimentalCopyVector)
      .def("CopyMatrix", ExperimentalCopyMatrix)
      .def("CopyMatrixValues", ExperimentalCopyMatrixValues)
      // --- Pointer factories ---
      .def("CreateEmptyMatrixPointer", ExperimentalCreateEmptyMatrixPointer)
      .def("CreateEmptyVectorPointer", ExperimentalCreateEmptyVectorPointer)
      .def("CreateVectorCopy", ExperimentalCreateVectorCopy)
      .def("CreateMatrixCopy", ExperimentalCreateMatrixCopy)
      // --- Element access ---
      .def("GetValue", ExperimentalGetValue)
      .def("SetValueVector", ExperimentalSetValueVector)
      .def("SetValueMatrix", ExperimentalSetValueMatrix)
      .def("GatherValues", ExperimentalGatherValues)
      // --- Set scalar ---
      .def("Set", ExperimentalSetV)
      // --- Assembly ---
      .def("GlobalAssembleMatrix", ExperimentalGlobalAssembleMatrix)
      .def("GlobalAssembleVector", ExperimentalGlobalAssembleVector)
      .def("ManualFinalize", ExperimentalManualFinalize)
      // --- Static queries ---
      .def_static("IsDistributed",
                  &ExperimentalTrilinosSparseSpaceType::IsDistributed)
      .def_static("IsDistributedSpace",
                  &ExperimentalTrilinosSparseSpaceType::IsDistributedSpace)
      .def_static(
          "FastestDirectSolverList",
          &ExperimentalTrilinosSparseSpaceType::FastestDirectSolverList);

#endif
}

} // namespace Kratos::Python.
