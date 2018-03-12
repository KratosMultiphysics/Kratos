//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
// #include "spaces/ublas_space.h"
// #include "add_trilinos_linear_solvers_to_python.h"
#include "includes/model_part.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

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



typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
//typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

typedef Epetra_FECrsMatrix FECrsMatrix;

void prova(TrilinosSparseSpaceType& dummy, FECrsMatrix& rX)
{
    rX.PutScalar(0.0);
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

void ResizeVector(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorPointerType& px, unsigned int i1)
{
    KRATOS_WATCH("Within ResizeVector")
    dummy.Resize(px, i1);
}

void SetToZeroMatrix(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType& A)
{
    dummy.SetToZero(A);
}

void SetToZeroVector(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x)
{
    dummy.SetToZero(x);
}

void ClearMatrix(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixPointerType& pA)
{
    dummy.Clear(pA);
}

void ClearVector(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorPointerType& px)
{
    dummy.Clear(px);
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

TrilinosSparseSpaceType::MatrixPointerType CreateEmptyMatrixPointer(TrilinosSparseSpaceType& dummy, Epetra_MpiComm& Comm)
{
    return dummy.CreateEmptyMatrixPointer(Comm);
}

TrilinosSparseSpaceType::VectorPointerType CreateEmptyVectorPointer(TrilinosSparseSpaceType& dummy, Epetra_MpiComm& Comm)
{
    return dummy.CreateEmptyVectorPointer(Comm);
}

Epetra_FECrsMatrix& GetMatRef(TrilinosSparseSpaceType::MatrixPointerType& dummy)
{
    KRATOS_WATCH("insdie GetMatRef")
    KRATOS_WATCH(dummy)
    KRATOS_WATCH(&dummy)
    return *(dummy.get());
}

Epetra_FEVector& GetVecRef(TrilinosSparseSpaceType::VectorPointerType& dummy)
{
    return *(dummy.get());
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

    class_< Epetra_MpiComm > (m,"Epetra_MpiComm")
    .def(init< Epetra_MpiComm& >())
    .def("MyPID",&Epetra_MpiComm::MyPID)
    .def("NumProc",&Epetra_MpiComm::NumProc)
    ;

    class_< Epetra_FECrsMatrix > (m,"Epetra_FECrsMatrix")
    .def(init< Epetra_FECrsMatrix& >())
    .def("__repr__",[](const Epetra_FECrsMatrix& self){
            std::stringstream ss;
            ss << self;
            return ss.str();
        })
    ;

    class_< Epetra_FEVector > (m,"Epetra_FEVector").def(init< Epetra_FEVector& >())
    .def("SetValue", SetValue)
    .def("__repr__",[](const Epetra_FEVector& self){
            std::stringstream ss;
            ss << self;
            return ss.str();
        })
    ;

    class_< TrilinosSparseSpaceType::MatrixPointerType > (m,"TrilinosMatrixPointer")//.def(init< TrilinosSparseSpaceType::MatrixPointerType > ())
    .def("GetReference", GetMatRef, return_value_policy::reference)
    ;

    class_< TrilinosSparseSpaceType::VectorPointerType > (m,"TrilinosVectorPointer")//.def(init< TrilinosSparseSpaceType::VectorPointerType > ())
    .def("GetReference", GetVecRef, return_value_policy::reference)
    ;

    //typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;
    //typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
    //typedef TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;

    //********************************************************************
    //********************************************************************


    class_< TrilinosSparseSpaceType> (m,"TrilinosSparseSpace")
    .def(init<>())
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
    // 		 .def("TransposeMult", TransposeMult)
    .def("Size", Size)
    .def("Size1", Size1)
    .def("Size2", Size2)
    .def("UnaliasedAdd", UnaliasedAdd)
    .def("ScaleAndAdd", ScaleAndAdd)
    .def("CreateEmptyMatrixPointer", CreateEmptyMatrixPointer)
    .def("CreateEmptyVectorPointer", CreateEmptyVectorPointer)
    .def("ReadMatrixMarketMatrix", &TrilinosSparseSpaceType::ReadMatrixMarket)
    .def("SetValue", SetValue)
    .def("__repr__",[](const TrilinosSparseSpaceType& self){
            std::stringstream ss;
            self.PrintInfo(ss);
            return ss.str();
        })
    ;


    m.def("CreateCommunicator", CreateCommunicator);
    m.def("ErrorCleaner", ErrorCleaner);

    //********************************************************************
    //********************************************************************
    class_< Teuchos::ParameterList > (m,"ParameterList").def(init<>())
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
