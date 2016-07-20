/*
==============================================================================
KratosTestApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//   
//   Project Name:        Kratos       
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "includes/model_part.h"

//#include "custom_utilities/CreateSPHParticle.h"
//#include "custom_utilities/CreateMLSParticle.h"
#include "custom_utilities/CreateMLSParticleGauss.h"
#include "custom_utilities/CreateLMEParticleGauss.h"

#include "custom_utilities/neighbours_calculator_SPH.h"
#include "custom_utilities/neighbours_calculator_MLS.h"
#include "custom_utilities/neighbours_calculator_LME.h"
//#include "custom_processes/node_and_element_erase_process.h"
#include "custom_utilities/nodal_values_utility.h"
#include "custom_utilities/gauss_coordiantes_update.h"
#include "custom_utilities/recompute_neighbours.h"

namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;


//     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//     typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

//     typedef Geometry<Node<3> >::GeometryType GeometryType;



	    
    class_<CreateMLSParticleGauss > ("CreateMLSParticleGauss", init<ModelPart& ,unsigned int >())  //the input parameters is a model part
            .def("Execute", &CreateMLSParticleGauss::Execute);
    class_<CreateLMEParticleGauss > ("CreateLMEParticleGauss", init<ModelPart& ,unsigned int >())  //the input parameters is a model part
            .def("Execute", &CreateLMEParticleGauss::Execute);


    //class_<CreateMLSParticle > ("CreateMLSParticle", init<ModelPart& ,unsigned int, float  >())  //the input parameters is a model part
            //.def("Execute", &CreateMLSParticle::Execute);
	    
    //class_<Neighbours_Calculator_MLS<LinearMLSKernel> > ("Neighbours_Calculator_MLS", init<>())  //the input parameters is a model part
            //.def("Search_Neighbours", &Neighbours_Calculator_MLS<LinearMLSKernel>::Search_Neighbours);
	    


//    class_<WCSPHUtils>("WCSPHUtils", init<>())
//            .def("UpdateDensities",&WCSPHUtils::UpdateDensities)
//            .def("UpdatePositions",&WCSPHUtils::UpdatePositions)
//            .def("UpdateVelocities",&WCSPHUtils::UpdateVelocities)
//            .def("ApplyFixedBoundary",&WCSPHUtils::ApplyFixedBoundary)
//            .def("ApplyFreeSurfaceBoundary",&WCSPHUtils::ApplyFreeSurfaceBoundary)
//            .def("UpdatePressures",&WCSPHUtils::UpdatePressures)
//            .def("CheckTimeStep",&WCSPHUtils::CheckTimeStep)
//            .def("ViscousTerm",&WCSPHUtils::ViscousTerm)
//            .def("GradientOfPressure",&WCSPHUtils::GradientOfPressure)
//            .def("RemoveBoundaryForce",&WCSPHUtils::RemoveBoundaryForce)
//            .def("DivergenceOfVelocity",&WCSPHUtils::DivergenceOfVelocity)
//            .def("GetPressureAcceleration",&WCSPHUtils::GetPressureAcceleration)
//            .def("GetViscousAcceleration",&WCSPHUtils::GetViscousAcceleration)
//            .def("GetNeighbours",&WCSPHUtils::GetNeighbours)
//            .def("Printer",&WCSPHUtils::Printer)
//            .def("ApplyInitialConditions",&WCSPHUtils::ApplyInitialConditions)
//            .def("MarkOuterNodes",&WCSPHUtils::MarkOuterNodes)
//            .def("GetBoundaryAcceleration",&WCSPHUtils::GetBoundaryAcceleration);


//    class_<PCISPHUtils>("PCISPHUtils", init<>())
//            .def("UpdatePositions",&PCISPHUtils::UpdatePositions)
//            .def("UpdateVelocities",&PCISPHUtils::UpdateVelocities)
//            .def("ApplyFixedBoundary",&PCISPHUtils::ApplyFixedBoundary)
//            .def("ApplyFreeSurfaceBoundary",&PCISPHUtils::ApplyFreeSurfaceBoundary)
//            .def("CheckTimeStep",&PCISPHUtils::CheckTimeStep)
//            .def("ViscousTerm",&PCISPHUtils::ViscousTerm)
//            .def("GradientOfPressure",&PCISPHUtils::GradientOfPressure)
//            .def("RemoveBoundaryForce",&PCISPHUtils::RemoveBoundaryForce)
//            .def("GetViscousAcceleration",&PCISPHUtils::GetViscousAcceleration)
//            .def("DivergenceOfVelocity",&PCISPHUtils::DivergenceOfVelocity)
//            .def("GetNeighbours",&PCISPHUtils::GetNeighbours)
//            .def("Printer",&PCISPHUtils::Printer)
//            .def("Iterate",&PCISPHUtils::Iterate)
//            .def("MarkOuterNodes",&PCISPHUtils::MarkOuterNodes)
//            .def("SetPressureToZero",&PCISPHUtils::SetPressureToZero)
//            .def("GetBoundaryAcceleration",&PCISPHUtils::GetBoundaryAcceleration);

//    class_<ISPHUtils>("ISPHUtils", init<>())
//            .def("UpdatePositions",&ISPHUtils::UpdatePositions)
//            .def("UpdateToFinalVelocities",&ISPHUtils::UpdateToFinalVelocities)
//            .def("UpdateToTemporaryVelocities",&ISPHUtils::UpdateToTemporaryVelocities)
//            .def("ApplyFixedBoundary",&ISPHUtils::ApplyFixedBoundary)
//            .def("ApplyFreeSurfaceBoundary",&ISPHUtils::ApplyFreeSurfaceBoundary)
//            .def("CheckTimeStep",&ISPHUtils::CheckTimeStep)
//            .def("ViscousTerm",&ISPHUtils::ViscousTerm)
//            .def("GradientOfPressure",&ISPHUtils::GradientOfPressure)
//            .def("RemoveBoundaryForce",&ISPHUtils::RemoveBoundaryForce)
//            .def("GetViscousAcceleration",&ISPHUtils::GetViscousAcceleration)
//            .def("GetNeighbours",&ISPHUtils::GetNeighbours)
//            .def("Printer",&ISPHUtils::Printer)
//            .def("MarkOuterNodes",&ISPHUtils::MarkOuterNodes)
//            .def("SetPressureToZero",&ISPHUtils::SetPressureToZero)
//            .def("SolvePressurePoissonEquation",&ISPHUtils::SolvePressurePoissonEquation)
//            .def("GetBoundaryAcceleration",&ISPHUtils::GetBoundaryAcceleration);

    //class_<NodeAndElementEraseProcess>("NodeAndElementEraseProcess", init < ModelPart& >())
            //.def("Execute", &NodeAndElementEraseProcess::Execute);


    class_<Nodal_Values_Utility>
    ( "Nodal_Values_Utility", init<ModelPart&, int >() )
    .def( "CalculateNodalStress", &Nodal_Values_Utility::CalculateNodalStress )
    .def( "CalculateNodalArea", &Nodal_Values_Utility::CalculateNodalArea )
    ;

    class_<Gauss_Coordinates_Update_Utility>
    ( "Gauss_Coordinates_Update_Utility", init<ModelPart& >() )
    .def( "UpdateGaussCoordinates", &Gauss_Coordinates_Update_Utility::UpdateGaussCoordinates )
    ;

    class_<Recompute_Neighbours_Utility>
    ( "RecomputeNeighbours", init<ModelPart& >() )
    .def( "Recompute_Neighbours", &Recompute_Neighbours_Utility::Recompute_Neighbours )
    ;



}





}  // namespace Python.

} // Namespace Kratos

