// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"
#include "utilities/process_factory_utility.h"
// Strategies

// Schemes

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#ifdef INCLUDE_MMG
    #include "custom_strategies/custom_convergencecriterias/error_mesh_criteria.h"
#endif

// Builders and solvers

// Linear solvers

namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef boost::shared_ptr<ProcessFactoryUtility> ProcessesListType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    
    // Base types
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
        
    // Custom strategy types
    
    // Custom scheme types

    // Custom convergence criterion types
#ifdef INCLUDE_MMG
    typedef ErrorMeshCriteria< SparseSpaceType,  LocalSpaceType > ErrorMeshCriteriaType;
#endif

    // Custom builder and solvers types
    
    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************
             
    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************

    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

#ifdef INCLUDE_MMG
    // Displacement Convergence Criterion
    class_< ErrorMeshCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "ErrorMeshCriteria", 
            init<ModelPart&, Parameters>())
            .def(init<ModelPart&, Parameters, ProcessesListType>())
            .def("SetEchoLevel", &ErrorMeshCriteriaType::SetEchoLevel)
            ;
#endif         
            
    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

}

}  // namespace Python.

} // Namespace Kratos

