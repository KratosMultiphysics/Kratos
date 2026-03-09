/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-23 12:50:01 $
//   Revision:            $Revision: 1.10 $
//
//


// System includes

// External includes



// Project includes
// #include "includes/define.h
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_meshers_to_python.h"
#include "includes/model_part.h"

#include "external_includes/tetgen_pfem_refine.h"
#include "external_includes/tetgen_mesh_suite_optimized.h"
#include "external_includes/trigen_mesh_suite.h"
#include "external_includes/trigen_refine.h"


namespace Kratos
{

namespace Python
{

/*	void TetRegenerate(TetGenModeler& Mesher,ModelPart& model_part, double alpha_shape)
	{
		Mesher.ReGenerateMesh(model_part,
			KratosComponents<Element>::Get("Fluid3D"),
			KratosComponents<Condition>::Get("SurfaceCondition3D3N"),alpha_shape	);
	}

	void TetRegenerateLagrangian(TetGenModeler& Mesher,ModelPart& model_part, double alpha_shape)
	{
		Mesher.ReGenerateMesh(model_part,
			KratosComponents<Element>::Get("TotalLagrangianFLuid"),
			KratosComponents<Condition>::Get("SurfaceCondition3D3N"),alpha_shape	);
	}
*/
void TetRegeneratePfemUlf3D(TetGenPfemModeler& Mesher,ModelPart& model_part, double alpha_shape)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("UpdatedLagrangianFluid3D"),
                          KratosComponents<Condition>::Get("SurfaceCondition3D3N"),alpha_shape	);
}
void TetRegeneratePfemUlf3DInc(TetGenPfemModeler& Mesher,ModelPart& model_part, double alpha_shape)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("UpdatedLagrangianFluid3Dinc"),
                          KratosComponents<Condition>::Get("SurfaceCondition3D3N"),alpha_shape	);
}
void TetRegeneratePfem3DInc(TetGenPfemModeler& Mesher,ModelPart& model_part, double alpha_shape)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("Fluid3D"),
                          KratosComponents<Condition>::Get("SurfaceCondition3D3N"),alpha_shape	);
}

void TriRefinePFEM(TriGenCDTrefine & Mesher,ModelPart& model_part,bool refine)
{
    Mesher.RefineCDT(model_part,
                     refine,
                     KratosComponents<Element>::Get("Fluid2D"));
}

void TriRegenerate(TriGenModeler& Mesher,ModelPart& model_part,double alpha_shape)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("Fluid2D"),
                          KratosComponents<Condition>::Get("LineCondition2D2N"),alpha_shape	);
}

void TriRegenerateCoupled(TriGenModeler& Mesher,ModelPart& model_part,double alpha_shape)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("Fluid2DCoupled"),
                          KratosComponents<Condition>::Get("LineCondition2D2N"),alpha_shape	);
}

void TriRegenerateUpdatedLagrangian(TriGenModeler& Mesher,ModelPart& model_part,double alpha_shape)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("UpdatedLagrangianFluid2D"),
                          KratosComponents<Condition>::Get("LineCondition2D2N"),alpha_shape	);
}
void TriRegenerateUpdatedLagrangianTest(TriGenModeler& Mesher,ModelPart& model_part,double alpha_shape)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("UpdatedLagrangianFluid2Dinc"),
                          KratosComponents<Condition>::Get("LineCondition2D2N"),alpha_shape	);
}
void TetRegenerateUpdatedLagrangian(TetGenModeler& Mesher,ModelPart& model_part,double alpha_shape)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("UpdatedLagrangianFluid3D"),
                          KratosComponents<Condition>::Get("SurfaceCondition3D3N"),alpha_shape	);
}

void TetRegenerateUpdatedLagrangianInc(TetGenModeler& Mesher,ModelPart& model_part,double alpha_shape)
{
    //KRATOS_WATCH("AAAAAAAAAKKKKKKKKKKKKKKKKKKKK")
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get("UpdatedLagrangianFluid3Dinc"),
                          KratosComponents<Condition>::Get("SurfaceCondition3D3N"),alpha_shape	);
}
/*
void TriRegenerateulf_pressure(TriGenModeler& Mesher,ModelPart& model_part,double alpha_shape)
{
	Mesher.ReGenerateMesh(model_part,
		KratosComponents<Element>::Get("ulf_pressure2D"),
		KratosComponents<Condition>::Get("LineCondition2D2N"),alpha_shape	);
}
*/
void  AddMeshersToPython(pybind11::module& m)
{

    namespace py = pybind11;

    py::class_<TetGenModeler >(m,"TetGenModeler")
    .def(py::init< >())
    // .def("ReGenerateMesh",TetRegenerate)
    // .def("ReGenerateMesh_Lagrangian",TetRegenerateLagrangian)
    .def("ReGenerateUpdatedLagrangian3D",TetRegenerateUpdatedLagrangian)
    .def("ReGenerateUpdatedLagrangian3Dinc",TetRegenerateUpdatedLagrangianInc)
    ;

    py::class_<TetGenPfemModeler >(m,"TetGenPfemModeler")
    .def(py::init< >())
    .def("ReGenerateMeshPfemUlf3D",TetRegeneratePfemUlf3D)
    .def("ReGenerateMeshPfemUlf3Dinc",TetRegeneratePfemUlf3DInc)
    .def("ReGenerateMeshPfem3Dinc",TetRegeneratePfem3DInc)
    ;

    py::class_<TriGenModeler >(m,"TriGenModeler")
    .def(py::init< >())
    .def("ReGenerateMesh",TriRegenerate)
    .def("ReGenerateMeshCoupled",TriRegenerateCoupled)
    .def("ReGenerateUpdatedLagrangian",TriRegenerateUpdatedLagrangian)
    .def("RegenerateUpdatedLagrangian2Dinc",TriRegenerateUpdatedLagrangianTest)
    // .def("ReGenerateulf_pressure",TriRegenerateulf_pressure)
    ;
    py::class_<TriGenCDTrefine >(m,"TriRefine")
    .def(py::init< >())
    .def("RefineMesh",RefineCDT)
    ;

}

}  // namespace Python.

} // Namespace Kratos

