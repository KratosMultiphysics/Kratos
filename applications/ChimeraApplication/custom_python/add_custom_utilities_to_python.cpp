//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================
// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Processes

#include "custom_utilities/vtk_output.hpp"
#include "custom_utilities/quadtree_binary_cell.h"
#include "custom_utilities/quadtree_binary.h"
//#include "custom_utilities/multipoint_constraint_data.hpp"

#include "custom_utilities/interpolation_utility.h"

namespace Kratos
{

template< unsigned int TDim >
void AuxCreateFluidBoundaryFaces(InterpolationUtility<TDim>& ThisUtility, ModelPart& rBackground)
{
    std::string ConditionName("ChimeraFluidCouplingCondition2D");
    const Condition& rReferenceCondition = KratosComponents<Condition>::Get(ConditionName);
    ThisUtility.CreateBoundaryFaces(rBackground,rReferenceCondition);
}


template< unsigned int TDim >
void AuxCreateThermalBoundaryFaces(InterpolationUtility<TDim>& ThisUtility, ModelPart& rBackground)
{
    std::string ConditionName("ChimeraThermalCouplingCondition2D");
    const Condition& rReferenceCondition = KratosComponents<Condition>::Get(ConditionName);
    ThisUtility.CreateBoundaryFaces(rBackground,rReferenceCondition);
}

namespace Python
{


  void  AddCustomUtilitiesToPython()
  {
    using namespace boost::python;


      //typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
      //typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
      //typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


          /// Processes
      /*class_<ApplyMultipointConstraintsProcessChimera, boost::noncopyable, bases<Process>>("ApplyMultipointConstraintsProcessChimera", init<ModelPart&>())
      .def(init< ModelPart&, Parameters& >())
      .def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcessChimera::AddMasterSlaveRelation)
      .def("PrintData", &ApplyMultipointConstraintsProcessChimera::PrintData);*/

      class_<VtkOutput, boost::noncopyable>("VtkOutput", init< ModelPart&, std::string, Parameters >())
      .def("PrintOutput", &VtkOutput::PrintOutput)
      .def("PrintOutput", &VtkOutput::PrintOutputSubModelPart);      

      class_ < InterpolationUtility<2>, boost::noncopyable >( "InterpolationUtility2D", init<>() )
              .def("Test", &InterpolationUtility<2>::Test )
              //.def("DirichletBoundaryCalculation", &InterpolationUtility < 2 > ::DirichletBoundaryCalculation<double>)
              .def("Projection_Flux", &InterpolationUtility < 2 > ::Projection_Flux<double>)
              //.def("DirichletBoundaryCalculationKreis", &InterpolationUtility < 2 > ::DirichletBoundaryCalculationKreis<double>)
              //.def("NeumannBoundaryCalculationRechteck", &InterpolationUtility < 2 > ::NeumannBoundaryCalculationRechteck<double>)
              //.def("DirichletBoundaryCalculationRechteck", &InterpolationUtility < 2 > ::DirichletBoundaryCalculationRechteck<double>)
              .def("CreateFluidBoundaryFaces", AuxCreateFluidBoundaryFaces<2>)
              .def("CreateThermalBoundaryFaces", AuxCreateThermalBoundaryFaces<2>)
              //.def("DirichletBoundaryCalculationFlag", &InterpolationUtility < 2 > ::DirichletBoundaryCalculationFlag<double>)
              .def("Projection_Vector", &InterpolationUtility < 2 > ::Projection_Vector< array_1d<double,3> >)
              .def("Projection_Scalar", &InterpolationUtility < 2 > ::Projection_Vector< double >)
              .def("Projection_Traction", &InterpolationUtility < 2 > ::Projection_Traction)
              .def("Projection_Traction_Cavity", &InterpolationUtility < 2 > ::Projection_Traction_Cavity)
              .def("Projection_Traction_Nodal", &InterpolationUtility < 2 > ::Projection_Traction_Nodal)
              .def("Projection_Flux_Gauss", &InterpolationUtility < 2 > :: Projection_Flux_Gauss)
              ;

      class_ < InterpolationUtility<3>, boost::noncopyable >( "InterpolationUtility3D", init<>() )
              .def("Test", &InterpolationUtility<3>::Test )
              //.def("DirichletBoundaryCalculation", &InterpolationUtility < 2 > ::DirichletBoundaryCalculation<double>)

              //.def("DirichletBoundaryCalculationKreis", &InterpolationUtility < 2 > ::DirichletBoundaryCalculationKreis<double>)
              //.def("NeumannBoundaryCalculationRechteck", &InterpolationUtility < 2 > ::NeumannBoundaryCalculationRechteck<double>)
              //.def("DirichletBoundaryCalculationRechteck", &InterpolationUtility < 2 > ::DirichletBoundaryCalculationRechteck<double>)
              //.def("CreateBoundaryFaces", &InterpolationUtility < 3 > ::CreateBoundaryFaces)
              //.def("DirichletBoundaryCalculationFlag", &InterpolationUtility < 2 > ::DirichletBoundaryCalculationFlag<double>)
              .def("Projection_Vector", &InterpolationUtility < 3 > ::Projection_Vector< array_1d<double,3> >)
              .def("Projection_Scalar", &InterpolationUtility < 3 > ::Projection_Vector< double >)
              .def("Projection_Traction", &InterpolationUtility < 3 > ::Projection_Traction)
              .def("ExtractBoundaryFaces",&InterpolationUtility < 3 >::ExtractBoundaryFaces)
              ;


  }





}  // namespace Python.

} // Namespace Kratos
