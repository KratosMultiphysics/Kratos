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

#include "custom_python/add_trilinos_convergence_criterias_to_python.h" 

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
#include "spaces/ublas_space.h"
#include "add_trilinos_linear_solvers_to_python.h"
#include "includes/model_part.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
//#include "solving_strategies/convergencecriterias/displacement_criteria.h"
// 
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"



namespace Kratos {

    namespace Python {

        using namespace boost::python;


	void  AddConvergenceCriterias()
        {
        typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
      //typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
	typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

        typedef Epetra_FECrsMatrix FECrsMatrix;


            //********************************************************************
            //********************************************************************
            //convergence criteria base class
            typedef ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosConvergenceCriteria;
            class_< TrilinosConvergenceCriteria, boost::noncopyable > ("TrilinosConvergenceCriteria", init<>())
                    .def("SetActualizeRHSFlag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::SetActualizeRHSFlag)
                    .def("GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag)
                    .def("PreCriteria", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::PreCriteria)
                    .def("PostCriteria", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::PostCriteria)
                    .def("Initialize", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Initialize)
                    .def("InitializeSolutionStep", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::FinalizeSolutionStep)
                    ;

            class_< TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
                    bases< TrilinosConvergenceCriteria >,
                    boost::noncopyable >
                    ("TrilinosDisplacementCriteria", init< double, double, Epetra_MpiComm& >());

            class_< TrilinosUPCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
                    bases< TrilinosConvergenceCriteria >,
                    boost::noncopyable >
                    ("TrilinosUPCriteria", init< double, double, double, double, Epetra_MpiComm& >());

        }


    } // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
