//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

//Trilinos includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

// Project includes
#include "includes/define.h"
#include "trilinos_application.h"
#include "trilinos_space.h"
#include "custom_python/trilinos_pointer_wrapper.h"
#include "custom_python/add_trilinos_space_to_python.h"
// #include "spaces/ublas_space.h"
// #include "add_trilinos_linear_solvers_to_python.h"
#include "includes/model_part.h"

// Teuchos parameter list
#include "Teuchos_ParameterList.hpp"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;
typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
//typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

typedef Epetra_FECrsMatrix FECrsMatrix;

void EraseAll(std::string& ThisString, std::string ToBeRemoved)
{
    int position;
    while ((position = ThisString.find_first_of(ToBeRemoved)) >= 0)
    {
        ThisString.erase(position, ToBeRemoved.size());
    }

}

std::string ErrorCleaner(std::string const& Input)
{
    std::string output(Input);

    EraseAll(output, "boost::numeric::");

    return output;
}

double Dot(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& rX, TrilinosSparseSpaceType::VectorType& rY)
{
    return dummy.Dot(rX, rY);
}

void ScaleAndAdd(TrilinosSparseSpaceType& dummy, const double A, const TrilinosSparseSpaceType::VectorType& rX, const double B, TrilinosSparseSpaceType::VectorType& rY)
// rY = (A * rX) + (B * rY)
{
    dummy.ScaleAndAdd(A, rX, B, rY);
}

void Mult(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType& rA, TrilinosSparseSpaceType::VectorType& rX, TrilinosSparseSpaceType::VectorType& rY)
//rY=A*rX (the product is stored inside the rY)
{
    dummy.Mult(rA, rX, rY);
}

void TransposeMult(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType& rA, TrilinosSparseSpaceType::VectorType& rX, TrilinosSparseSpaceType::VectorType& rY)
{
    dummy.TransposeMult(rA, rX, rY);
}

TrilinosSparseSpaceType::IndexType Size(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType const& rV)
{
    return dummy.Size(rV);
}

TrilinosSparseSpaceType::IndexType Size1(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType const& rM)
{
    return dummy.Size1(rM);
}

TrilinosSparseSpaceType::IndexType Size2(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType const& rM)
{
    return dummy.Size2(rM);
}

// void ResizeMatrix(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType& A, unsigned int i1, unsigned int i2)
// {
//     dummy.Resize(A, i1, i2);
// }

void ResizeVector(TrilinosSparseSpaceType& dummy, AuxiliaryVectorWrapper& px, unsigned int i1)
{
    dummy.Resize(px.GetPointer(), i1);
}

void SetToZeroMatrix(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType& A)
{
    dummy.SetToZero(A);
}

void SetToZeroVector(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x)
{
    dummy.SetToZero(x);
}

void ClearMatrix(TrilinosSparseSpaceType& dummy, AuxiliaryMatrixWrapper& pA)
{
    dummy.Clear(pA.GetPointer());
}

void ClearVector(TrilinosSparseSpaceType& dummy, AuxiliaryVectorWrapper& px)
{
    dummy.Clear(px.GetPointer());
}

double TwoNorm(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x)
{
    return dummy.TwoNorm(x);
}

void UnaliasedAdd(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x, const double A, const TrilinosSparseSpaceType::VectorType& rY) // x+= a*Y
{
    dummy.UnaliasedAdd(x, A, rY);
}

// bool (TrilinosLinearSolverType::*pointer_to_solve)(TrilinosLinearSolverType::SparseMatrixType& rA, TrilinosLinearSolverType::VectorType& rX, TrilinosLinearSolverType::VectorType& rB) = &TrilinosLinearSolverType::Solve;

//************************************************************************************************

Epetra_MpiComm CreateCommunicator()
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    return comm;
}

//************************************************************************************************


//************************************************************************************************

AuxiliaryMatrixWrapper CreateEmptyMatrixPointer(TrilinosSparseSpaceType& dummy, Epetra_MpiComm& Comm)
{
    return AuxiliaryMatrixWrapper(dummy.CreateEmptyMatrixPointer(Comm));
}

AuxiliaryVectorWrapper CreateEmptyVectorPointer(TrilinosSparseSpaceType& dummy, Epetra_MpiComm& Comm)
{
    return AuxiliaryVectorWrapper(dummy.CreateEmptyVectorPointer(Comm));
}

AuxiliaryMatrixWrapper ReadMatrixMarketMatrix(TrilinosSparseSpaceType& dummy, const std::string FileName,Epetra_MpiComm& Comm)
{
    return AuxiliaryMatrixWrapper(dummy.ReadMatrixMarket(FileName, Comm));
}

Epetra_FECrsMatrix& GetMatRef(AuxiliaryMatrixWrapper& dummy)
{
    return dummy.GetReference();
}

Epetra_FEVector& GetVecRef(AuxiliaryVectorWrapper& dummy)
{
    return dummy.GetReference();
}

void SetValue(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x, std::size_t i, double value)
{
    dummy.SetValue(x,i,value);
}



//************************************************************************************************
//************************************************************************************************
//teuchos paramter list

  void SetDoubleValue(Teuchos::ParameterList& dummy, const std::string &name, double value)
{
    dummy.set(name, value);
}

  void SetIntValue(Teuchos::ParameterList& dummy, const std::string &name, int value)
{
    dummy.set(name, value);
}

  void SetCharValue(Teuchos::ParameterList& dummy, const std::string &name, const char value[])
{
    dummy.set(name, value);
}

  void SetBoolValue(Teuchos::ParameterList& dummy, const std::string &name, int value)
{
    if(value == 0)
        dummy.set(name, false);
    else
        dummy.set(name, true);
}

  void SetSublistIntValue(Teuchos::ParameterList& dummy, const std::string &sublist_name, const std::string &name, int value)
{
    dummy.sublist(sublist_name).set(name, value);
}

  void SetSublistDoubleValue(Teuchos::ParameterList& dummy, const std::string &sublist_name, const std::string &name, double value)
{
    dummy.sublist(sublist_name).set(name, value);
}

  void SetSublistCharValue(Teuchos::ParameterList& dummy, const std::string &sublist_name, const std::string &name, const char value[])
{
    dummy.sublist(sublist_name).set(name, value);
}

  void SetSublistBoolValue(Teuchos::ParameterList& dummy, const std::string &sublist_name, const std::string &name, double value)
{
    dummy.sublist(sublist_name).set(name, value);
}

void  AddBasicOperations(pybind11::module& m)
{
    py::class_< Epetra_MpiComm > (m,"Epetra_MpiComm")
    .def(py::init< Epetra_MpiComm& >())
    .def("MyPID",&Epetra_MpiComm::MyPID)
    .def("NumProc",&Epetra_MpiComm::NumProc)
    ;

    //NOTE: deliberatly avoiding defining a Pointer handler, to make it incompatible with the Kratos. all uses should pass through the AuxiliaryMatrixWrapper
    py::class_< Epetra_FECrsMatrix  > (m,"Epetra_FECrsMatrix")
    .def(py::init< Epetra_FECrsMatrix& >())
    .def("__str__", PrintObject<Epetra_FECrsMatrix>)
    ;

    //NOTE: deliberatly avoiding defining a Pointer handler, to make it incompatible with the Kratos. all uses should pass through the AuxiliaryVectorWrapper
    py::class_< Epetra_FEVector > (m,"Epetra_FEVector")
    .def(py::init< Epetra_FEVector& >())
    .def("SetValue", SetValue)
    .def("__str__", PrintObject<Epetra_FEVector>)
    ;

    py::class_< AuxiliaryMatrixWrapper > (m,"TrilinosMatrixPointer")//.def(py::init< TrilinosSparseSpaceType::MatrixPointerType > ())
    .def("GetReference", GetMatRef, py::return_value_policy::reference_internal)
    ;

    py::class_< AuxiliaryVectorWrapper > (m,"TrilinosVectorPointer")//.def(py::init< TrilinosSparseSpaceType::VectorPointerType > ())
    .def("GetReference", GetVecRef, py::return_value_policy::reference_internal)
    ;

    //typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;
    //typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
    //typedef TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;

    //********************************************************************
    //********************************************************************


    py::class_< TrilinosSparseSpaceType> (m,"TrilinosSparseSpace")
    .def(py::init<>())
    .def("ClearMatrix", ClearMatrix)
    .def("ClearVector", ClearVector)
//     .def("ResizeMatrix", ResizeMatrix)
    .def("ResizeVector", ResizeVector)
    .def("SetToZeroMatrix", SetToZeroMatrix)
    .def("SetToZeroVector", SetToZeroVector)
    .def("TwoNorm", TwoNorm)
    //the dot product of two vectors
    .def("Dot", Dot)
    //the matrix-vector multiplication
    .def("Mult", Mult)
    //          .def("TransposeMult", TransposeMult)
    .def("Size", Size)
    .def("Size1", Size1)
    .def("Size2", Size2)
    .def("UnaliasedAdd", UnaliasedAdd)
    .def("ScaleAndAdd", ScaleAndAdd)
    .def("CreateEmptyMatrixPointer", CreateEmptyMatrixPointer)
    .def("CreateEmptyVectorPointer", CreateEmptyVectorPointer)
    .def("ReadMatrixMarketMatrix", ReadMatrixMarketMatrix)
    .def("SetValue", SetValue)
    ;

    m.def("CreateCommunicator", CreateCommunicator);
    m.def("ErrorCleaner", ErrorCleaner);

    //********************************************************************
    //********************************************************************
    py::class_< Teuchos::ParameterList > (m,"ParameterList").def(py::init<>())
    .def("set", SetDoubleValue)
    .def("set", SetIntValue)
    .def("set", SetCharValue)
    .def("setboolvalue", SetBoolValue)
    .def("SetSublistIntValue", SetSublistIntValue)
    .def("SetSublistDoubleValue", SetSublistDoubleValue)
    .def("SetSublistCharValue", SetSublistCharValue)
    .def("SetSublistBoolValue", SetSublistBoolValue)
//     .def(self_ns::str(self))
    ;
}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
