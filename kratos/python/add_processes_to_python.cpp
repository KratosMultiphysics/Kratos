// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-22 17:13:57 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

#include "processes/process.h"
#include "python/add_processes_to_python.h"
#include "processes/find_nodal_h_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_conditions_neighbours_process.h"
#include "processes/find_elements_neighbours_process.h"
#include "processes/calculate_nodal_area_process.h"
#include "processes/node_erase_process.h"
#include "processes/element_erase_process.h"
#include "processes/condition_erase_process.h"
#include "processes/eliminate_isolated_nodes_process.h"
#include "processes/calculate_signed_distance_to_3d_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h" 
#include "processes/translation_operation.h"
#include "processes/rotation_operation.h"
#include "processes/tetrahedral_mesh_orientation_check.h"
#include "processes/compute_bdfcoefficients_process.h"
#include "processes/variational_distance_calculation_process.h"
#include "processes/levelset_convection_process.h"
#include "processes/apply_constant_scalarvalue_process.h"
#include "processes/apply_constant_vectorvalue_process.h"
#include "processes/check_skin_process.h"
#include "processes/replace_elements_and_condition_process.h"
#include "includes/node.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{
    



void  AddProcessesToPython()
{
    using namespace boost::python;
    

    
    class_<Process>("Process")
            .def("Execute",&Process::Execute)
            .def("ExecuteInitialize",&Process::ExecuteInitialize)
            .def("ExecuteBeforeSolutionLoop",&Process::ExecuteBeforeSolutionLoop)
            .def("ExecuteInitializeSolutionStep",&Process::ExecuteInitializeSolutionStep)
            .def("ExecuteFinalizeSolutionStep",&Process::ExecuteFinalizeSolutionStep)
            .def("ExecuteBeforeOutputStep",&Process::ExecuteBeforeOutputStep)
            .def("ExecuteAfterOutputStep",&Process::ExecuteAfterOutputStep)
            .def("ExecuteFinalize",&Process::ExecuteFinalize)
            .def(self_ns::str(self))
    ;

    class_<FindNodalHProcess, bases<Process> >("FindNodalHProcess",init<ModelPart&>())
    .def("Execute",&FindNodalHProcess::Execute)
    ;
    class_<FindNodalNeighboursProcess, bases<Process> >("FindNodalNeighboursProcess",
            init<ModelPart&, int, int>())
    .def("ClearNeighbours",&FindNodalNeighboursProcess::ClearNeighbours)
    ;

    class_<FindConditionsNeighboursProcess, bases<Process> >("FindConditionsNeighboursProcess",
            init<ModelPart&, int, int>())
    .def("ClearNeighbours",&FindConditionsNeighboursProcess::ClearNeighbours)
    ;

    class_<FindElementalNeighboursProcess, bases<Process> >("FindElementalNeighboursProcess",
            init<ModelPart&, int, int>())
    .def("ClearNeighbours",&FindElementalNeighboursProcess::ClearNeighbours)
    ;

    class_<CalculateNodalAreaProcess, bases<Process> >("CalculateNodalAreaProcess",
            init<ModelPart&, unsigned int>())
    ;

    class_<NodeEraseProcess, bases<Process> >("NodeEraseProcess",
            init<ModelPart&>())
    ;
    
    class_<ElementEraseProcess, bases<Process> >("ElementEraseProcess",
            init<ModelPart&>())
    ;
    
    class_<ConditionEraseProcess, bases<Process> >("ConditionEraseProcess",
            init<ModelPart&>())
    ;

    class_<EliminateIsolatedNodesProcess, bases<Process> >("EliminateIsolatedNodesProcess",
            init<ModelPart&>())
    ;

    class_<CalculateSignedDistanceTo3DSkinProcess, bases<Process>, boost::noncopyable >("CalculateSignedDistanceTo3DSkinProcess",
            init<ModelPart&, ModelPart&>())
    .def("GenerateSkinModelPart",&CalculateSignedDistanceTo3DSkinProcess::GenerateSkinModelPart)
    .def("MappingPressureToStructure",&CalculateSignedDistanceTo3DSkinProcess::MappingPressureToStructure)
    ;

   class_<CalculateSignedDistanceTo3DConditionSkinProcess, bases<Process> >("CalculateSignedDistanceTo3DConditionSkinProcess",
            init<ModelPart&, ModelPart&>())
    ;
    
    class_<TranslationOperation, bases<Process> >("TranslationOperation",
            init<ModelPart&, boost::numeric::ublas::vector<int> ,boost::numeric::ublas::vector<int> ,unsigned int>())
    ;

    class_<RotationOperation, bases<Process> >("RotationOperation",
            init<ModelPart&, boost::numeric::ublas::vector<int> ,boost::numeric::ublas::vector<int> ,unsigned int>())
    ;

    class_<TetrahedralMeshOrientationCheck, bases<Process>, boost::noncopyable >("TetrahedralMeshOrientationCheck",
            init<ModelPart&, bool>())
            .def("SwapAll",&TetrahedralMeshOrientationCheck::SwapAll)
            .def("SwapNegativeElements",&TetrahedralMeshOrientationCheck::SwapNegativeElements)
    ;
    
    class_<ComputeBDFCoefficientsProcess, bases<Process>, boost::noncopyable >("ComputeBDFCoefficientsProcess",
            init<ModelPart&, const unsigned int>())
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    class_<VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType > , bases<Process>, boost::noncopyable >("VariationalDistanceCalculationProcess2D",
            init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;
    class_<VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType > , bases<Process>, boost::noncopyable >("VariationalDistanceCalculationProcess3D",
            init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;    
    
    class_<LevelSetConvectionProcess<2> , bases<Process>, boost::noncopyable >("LevelSetConvectionProcess2D",
            init<Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double >())
    ;
    class_<LevelSetConvectionProcess<3> , bases<Process>, boost::noncopyable >("LevelSetConvectionProcess3D",
            init<Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double>())
    ;   


    class_<ApplyConstantScalarValueProcess , bases<Process>, boost::noncopyable >("ApplyConstantScalarValueProcess",
            init<ModelPart&, const Variable<double>&, double, std::size_t, Flags>())
            .def(init< ModelPart&, Parameters& >())
            .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, double, std::size_t, Flags>())
            .def(init<ModelPart&, const Variable<int>&, int, std::size_t, Flags>())
            .def(init<ModelPart&, const Variable<bool>&, bool, std::size_t, Flags>())
            .def_readonly("VARIABLE_IS_FIXED", &ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)
    ; 

    class_<ApplyConstantVectorValueProcess , bases<Process>, boost::noncopyable >("ApplyConstantVectorValueProcess",
            init<ModelPart&, const Variable<array_1d<double, 3 > >& , const double, const Vector , std::size_t, Flags>())
            .def(init< ModelPart&, Parameters& >())
            .def_readonly("X_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::X_COMPONENT_FIXED)
            .def_readonly("Y_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Y_COMPONENT_FIXED)
            .def_readonly("Z_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Z_COMPONENT_FIXED)
    ; 

    class_<CheckSkinProcess , bases<Process>, boost::noncopyable >("CheckSkinProcess",
            init<ModelPart&, Flags>())
    ; 
        
    class_<ReplaceElementsAndConditionsProcess , bases<Process>, boost::noncopyable >("ReplaceElementsAndConditionsProcess",
            init<ModelPart&, Parameters>())
    ;
    
    //typedef PointerVectorSet<Node<3>, IndexedObject> NodesContainerType;
    //typedef PointerVectorSet<Dof<double>, IndexedObject> DofsContainerType;

    //class_<AddDofsNodalProcess<Variable<double> >, bases<Process> >("AddDoubleDofsNodalProcess")
    // .def(init<Variable<double>, NodesContainerType&, DofsContainerType&>())
    // ;
    //class_<AddDofsNodalProcess<VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >, bases<Process> >("AddArrayComponentDofsNodalProcess")
    // ;
}

}  // namespace Python.

} // Namespace Kratos

