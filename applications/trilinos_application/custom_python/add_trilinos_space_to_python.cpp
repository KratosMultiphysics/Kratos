//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-12-09 20:20:55 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>

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

//strategies
// #include "solving_strategies/strategies/solving_strategy.h"
// #include "solving_strategies/strategies/residualbased_linear_strategy.h"
// #include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

//schemes
// #include "solving_strategies/schemes/scheme.h"
// #include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
// #include "custom_strategies/schemes/trilinos_residualbased_lagrangian_monolithic_scheme.h"
// #include "../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
// #include "custom_strategies/schemes/trilinos_predictorcorrector_velocity_bossak_scheme.h"

//convergence criterias
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"
// #include "solving_strategies/convergencecriterias/displacement_criteria.h"
// 
// //Builder And Solver
// // #include "solving_strategies/builder_and_solvers/builder_and_solver.h"
// #include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver.h"
// #include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
// #include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_vec.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_mixed.h"

//linear solvers
// #include "linear_solvers/linear_solver.h"

//utilities
#include "python/pointer_vector_set_python_interface.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

// #include "external_includes/aztec_solver.h"
// #include "external_includes/amesos_solver.h"
// #include "external_includes/ml_solver.h"

//configuration files
// #include "../../incompressible_fluid_application/custom_strategies/strategies/solver_configuration.h"
// #include "custom_strategies/strategies/trilinos_fractionalstep_configuration.h"
// #include "../../incompressible_fluid_application/custom_strategies/strategies/fractional_step_strategy.h"
// #include "../../incompressible_fluid_application/incompressible_fluid_application.h"



namespace Kratos {

    namespace Python {

        using namespace boost::python;

        void EraseAll(std::string& ThisString, std::string ToBeRemoved) {
            int position;
            while ((position = ThisString.find_first_of(ToBeRemoved)) >= 0) {
                ThisString.erase(position, ToBeRemoved.size());
            }

        }

        std::string ErrorCleaner(std::string const& Input) {
            std::string output(Input);

            EraseAll(output, "boost::numeric::");

            return output;
        }



        typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
      //typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

        typedef Epetra_FECrsMatrix FECrsMatrix;

        void prova(TrilinosSparseSpaceType& dummy, FECrsMatrix& rX) {
            rX.PutScalar(0.0);
        }

        double Dot(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& rX, TrilinosSparseSpaceType::VectorType& rY) {
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

        void TransposeMult(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType& rA, TrilinosSparseSpaceType::VectorType& rX, TrilinosSparseSpaceType::VectorType& rY) {
            dummy.TransposeMult(rA, rX, rY);
        }

        TrilinosSparseSpaceType::IndexType Size(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType const& rV) {
            return dummy.Size(rV);
        }

        void ResizeMatrix(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType& A, unsigned int i1, unsigned int i2) {
            dummy.Resize(A, i1, i2);
        }

        void ResizeVector(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x, unsigned int i1) {
            dummy.Resize(x, i1);
        }

        void SetToZeroMatrix(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixType& A) {
            dummy.SetToZero(A);
        }

        void SetToZeroVector(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x) {
            dummy.SetToZero(x);
        }

        void ClearMatrix(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::MatrixPointerType& pA) {
            dummy.Clear(pA);
        }

        void ClearVector(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorPointerType& px) {
            dummy.Clear(px);
        }

        double TwoNorm(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x) {
            return dummy.TwoNorm(x);
        }

        void UnaliasedAdd(TrilinosSparseSpaceType& dummy, TrilinosSparseSpaceType::VectorType& x, const double A, const TrilinosSparseSpaceType::VectorType& rY) // x+= a*Y
        {
            dummy.UnaliasedAdd(x, A, rY);
        }

        // bool (TrilinosLinearSolverType::*pointer_to_solve)(TrilinosLinearSolverType::SparseMatrixType& rA, TrilinosLinearSolverType::VectorType& rX, TrilinosLinearSolverType::VectorType& rB) = &TrilinosLinearSolverType::Solve;

        //************************************************************************************************

        Epetra_MpiComm CreateCommunicator() {
            Epetra_MpiComm comm(MPI_COMM_WORLD);
            return comm;
        }

        //************************************************************************************************


        //************************************************************************************************

        TrilinosSparseSpaceType::MatrixPointerType CreateEmptyMatrixPointer(TrilinosSparseSpaceType& dummy, Epetra_MpiComm& Comm) {
            return dummy.CreateEmptyMatrixPointer(Comm);
        }

        TrilinosSparseSpaceType::VectorPointerType CreateEmptyVectorPointer(TrilinosSparseSpaceType& dummy, Epetra_MpiComm& Comm) {
            return dummy.CreateEmptyVectorPointer(Comm);
        }

        Epetra_FECrsMatrix& GetMatRef(TrilinosSparseSpaceType::MatrixPointerType& dummy) {
            return *dummy;
        }

        Epetra_FEVector& GetVecRef(TrilinosSparseSpaceType::VectorPointerType& dummy) {
            return *dummy;
        }



        //************************************************************************************************
        //************************************************************************************************
        //teuchos paramter list

        void SetDoubleValue(Teuchos::ParameterList& dummy, const string &name, double value) {
            dummy.set(name, value);
        }

        void SetIntValue(Teuchos::ParameterList& dummy, const string &name, int value) {
            dummy.set(name, value);
        }

        void SetCharValue(Teuchos::ParameterList& dummy, const string &name, const char value[]) {
            dummy.set(name, value);
        }

        void SetBoolValue(Teuchos::ParameterList& dummy, const string &name, int value) {
            if(value == 0)
                dummy.set(name, false);
            else
                dummy.set(name, true);
        }

	void  AddBasicOperations()
        {

            class_< Epetra_MpiComm > ("Epetra_MpiComm", init< Epetra_MpiComm& >())
                    ;

            class_< Epetra_FECrsMatrix > ("Epetra_FECrsMatrix", init< Epetra_FECrsMatrix& >())
                    ;

            class_< Epetra_FEVector > ("Epetra_FEVector", init< Epetra_FEVector& >())
                    ;

            class_< TrilinosSparseSpaceType::MatrixPointerType > ("TrilinosMatrixPointer", init< TrilinosSparseSpaceType::MatrixPointerType > ())
                    .def("GetReference", GetMatRef, return_value_policy<reference_existing_object > ())
                    ;

            class_< TrilinosSparseSpaceType::VectorPointerType > ("TrilinosVectorPointer", init< TrilinosSparseSpaceType::VectorPointerType > ())
                    .def("GetReference", GetVecRef, return_value_policy<reference_existing_object > ())
                    ;

            //typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;
            //typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
            //typedef TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;

            //********************************************************************
            //********************************************************************


            class_< TrilinosSparseSpaceType, boost::noncopyable > ("TrilinosSparseSpace", init<>())
                    .def("ClearMatrix", ClearMatrix)
                    .def("ClearVector", ClearVector)
                    .def("ResizeMatrix", ResizeMatrix)
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
                    .def("UnaliasedAdd", UnaliasedAdd)
                    .def("ScaleAndAdd", ScaleAndAdd)
                    .def("prova", prova)
                    .def("CreateEmptyMatrixPointer", CreateEmptyMatrixPointer)
                    .def("CreateEmptyVectorPointer", CreateEmptyVectorPointer)
                    ;


            def("CreateCommunicator", CreateCommunicator);
            def("ErrorCleaner", ErrorCleaner);

            //********************************************************************
            //********************************************************************
            class_< Teuchos::ParameterList, boost::noncopyable > ("ParameterList", init<>())
                    .def("set", SetDoubleValue)
                    .def("set", SetIntValue)
                    .def("set", SetCharValue)
                    .def("setboolvalue", SetBoolValue)
                    .def(self_ns::str(self))
                    ;

 

        }


    } // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
