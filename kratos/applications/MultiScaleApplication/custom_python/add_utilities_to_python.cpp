/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

// System includes
#include <iostream>

// External includes
#include <boost/python.hpp>
#include "includes/constitutive_law.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

// Project includes
#include "add_utilities_to_python.h"
#include "custom_utilities/load_function.h"

#include "custom_utilities/rve_macroscale_status.h"
#include "custom_utilities/rve_boundary.h"
#include "custom_utilities/rve_boundary_2D.h"
#include "custom_utilities/rve_boundary_3D.h"
#include "custom_utilities/rve_utilities_model_part.h"
#include "custom_utilities/rve_utilities_element_info.h"
#include "custom_utilities/rve_adapter.h"
#include "constitutive_laws/rve_constitutive_law.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{
namespace Python
{

using namespace boost::python;

void TestOpenMP(int omp_max_num_threads, int nested)
{
#ifdef _OPENMP
	std::stringstream ss;
	ss << "SETTING OPENMP PARAMS\n";
	ss << "before:\n";
	ss << " - num: " << omp_get_num_threads() << std::endl;
	ss << " - nested: " << omp_get_nested() << std::endl;
	ss << "SETTING: " << omp_max_num_threads << ", " << nested << std::endl;
	omp_set_num_threads(omp_max_num_threads);
	omp_set_nested(nested);
	ss << " - num: " << omp_get_num_threads() << std::endl;
	ss << " - nested: " << omp_get_nested() << std::endl;
	std::cout << ss.str();
#endif
}

template<class Tobj, class T>
T ConstitutiveLawGetValue(Tobj& obj, const Variable<T>& a, T& b)
{
	return obj.GetValue(a, b);
}

void AddUtilitiesToPython()
{

	// ===============================================================================
	//
	// load functions
	//
	// ===============================================================================

	typedef LoadFunction<double> LoadFunctionbaseType;
	enum_<LoadFunctionbaseType::ProlongationType>("ProlongationType")
		.value("Zero", LoadFunctionbaseType::Zero)
		.value("Constant", LoadFunctionbaseType::Constant)
		.value("Linear", LoadFunctionbaseType::Linear)
		;

	class_<LoadFunctionbaseType, LoadFunctionbaseType::Pointer, boost::noncopyable >(
		"LoadFunction", 
		init<>())
		;
		
	typedef PieceWiseLoadFunction<double> PieceWiseLoadFunctionType;
	class_<PieceWiseLoadFunctionType, PieceWiseLoadFunctionType::Pointer, bases< LoadFunctionbaseType >, boost::noncopyable >(
		"PieceWiseLoadFunction",
		init<const boost::python::list &>())
		.def(init<const boost::python::list &, LoadFunctionbaseType::ProlongationType>())
		.def(init<const boost::python::list &, const boost::python::list &, LoadFunctionbaseType::ProlongationType>())
		;

	// ===============================================================================
	//
	// RVE -  These should be moved to .. AddRveToPython
	//
	// ===============================================================================


	class_<RveBoundary, RveBoundary::Pointer, boost::noncopyable>(
		"RveBoundary",
		init<>())
		.def(self_ns::str(self))
		;

	class_<RveBoundary2D, RveBoundary2D::Pointer, bases< RveBoundary >, boost::noncopyable>(
		"RveBoundary2D",
		init<ModelPart&>())
		.def(init<ModelPart&, double>())
		.def("AddConditions", &RveBoundary2D::AddConditions)
		;

	class_<RveBoundary3D, RveBoundary3D::Pointer, bases< RveBoundary >, boost::noncopyable>(
		"RveBoundary3D",
		init<ModelPart&>())
		.def(init<ModelPart&, double>())
		.def("AddConditions", &RveBoundary3D::AddConditions)
		;

	def("RveCloneModelPart", &RveUtilities::CloneModelPart);

	def("TestOpenMP", &TestOpenMP);

	class_<RveMacroscaleStatus, RveMacroscaleStatus::Pointer, boost::noncopyable>(
		"RveMacroscaleStatus",
		init<>())
		.def(self_ns::str(self))
		.add_property("StrainVector", 
			make_function(&RveMacroscaleStatus::GetStrainVector, return_value_policy<copy_const_reference>()), 
			&RveMacroscaleStatus::SetStrainVector)
		;

	typedef RveUtilities::SolidElementInfo    SolidElementInfoType;
	typedef SolidElementInfoType::IndexType   ElementInfoIndexType;
	typedef RveUtilities::ShellElementInfo    ShellElementInfoType;

	class_<SolidElementInfoType, boost::noncopyable>(
		"SolidElementInfo",
		init<>())
		.def(init<ElementInfoIndexType, ElementInfoIndexType>())
		.def("GetStringExtension", &SolidElementInfoType::GetStringExtension)
		.def(self_ns::str(self))
		.add_property("ElementID", &SolidElementInfoType::GetElementID, &SolidElementInfoType::SetElementID)
		.add_property("GaussPointID", &SolidElementInfoType::GetGaussPointID, &SolidElementInfoType::SetGaussPointID)
		;

	class_<ShellElementInfoType, bases<SolidElementInfoType>, boost::noncopyable>(
		"ShellElementInfo",
		init<>())
		.def(init<ElementInfoIndexType, ElementInfoIndexType>())
		.def(init<ElementInfoIndexType, ElementInfoIndexType, ElementInfoIndexType, ElementInfoIndexType>())
		.add_property("PlyID", &ShellElementInfoType::GetPlyID, &ShellElementInfoType::SetPlyID)
		.add_property("PlyIntegrationPointID", &ShellElementInfoType::GetPlyIntegrationPointID, &ShellElementInfoType::SetPlyIntegrationPointID)
		;



	typedef UblasSpace<double, CompressedMatrix, Vector>   SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector>			   LocalSpaceType;
	typedef LinearSolver<SparseSpaceType, LocalSpaceType>  LinearSolverBaseType;

	/*typedef RveAdapter<SparseSpaceType, LocalSpaceType, LinearSolverBaseType> RveAdapterType;
	class_<RveAdapterType, RveAdapterType::Pointer, boost::noncopyable>(
		"RveAdapter",
		init<>())
		.def("SetRveData", &RveAdapterType::SetRveData)
		.def("RveGenerated", &RveAdapterType::RveGenerated)
		.def("RveGenerationRequested", &RveAdapterType::RveGenerationRequested)
		.def("WorkingSpaceDimension", &RveAdapterType::WorkingSpaceDimension)
		.def("GetStrainSize", &RveAdapterType::GetStrainSize)
		.def(self_ns::str(self))
		.def("TestMaterialResponse", &RveAdapterType::TestMaterialResponse)
		;

	typedef ConstitutiveLawAdapter<RveAdapterType> RveConstitutiveLawBaseType;
	class_<RveConstitutiveLawBaseType, RveConstitutiveLawBaseType::Pointer, bases<ConstitutiveLaw>, boost::noncopyable>(
		"RveConstitutiveLawBaseType", no_init)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawBaseType, double>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawBaseType, Vector>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawBaseType, Matrix>)
        .def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawBaseType, array_1d<double,3> >)
        .def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawBaseType, array_1d<double,6> >)
		;

	typedef RveConstitutiveLaw<RveAdapterType> RveConstitutiveLawType;
	class_<RveConstitutiveLawType, RveConstitutiveLawType::Pointer, 
		   bases<RveConstitutiveLawBaseType>, 
		   boost::noncopyable>(
		"RveConstitutiveLaw",
		init<const RveAdapterType::Pointer&>())
		.def("GetModelPart", &RveConstitutiveLawType::GetModelPart)
		.def("TestMaterialResponse", &RveConstitutiveLawType::TestMaterialResponse)
		;*/

	// TYPES

	typedef RveAdapterSettings<3> RveAdapterSettings_PlaneStress;
	typedef RveAdapterSettings<6> RveAdapterSettings_3D;

	// PLANE STRESS

	typedef RveAdapter<SparseSpaceType, LocalSpaceType, LinearSolverBaseType, RveAdapterSettings_PlaneStress> RvePlaneStressAdapterType;
	class_<RvePlaneStressAdapterType, RvePlaneStressAdapterType::Pointer, boost::noncopyable>(
		"RvePlaneStressAdapter",
		init<>())
		.def("SetRveData", &RvePlaneStressAdapterType::SetRveData)
		.def("RveGenerated", &RvePlaneStressAdapterType::RveGenerated)
		.def("RveGenerationRequested", &RvePlaneStressAdapterType::RveGenerationRequested)
		.def("WorkingSpaceDimension", &RvePlaneStressAdapterType::WorkingSpaceDimension)
		.def("GetStrainSize", &RvePlaneStressAdapterType::GetStrainSize)
		.def(self_ns::str(self))
		.def("TestMaterialResponse", &RvePlaneStressAdapterType::TestMaterialResponse)
		;

	typedef ConstitutiveLawAdapter<RvePlaneStressAdapterType> RveConstitutiveLawPlaneStressBaseType;
	class_<RveConstitutiveLawPlaneStressBaseType, RveConstitutiveLawPlaneStressBaseType::Pointer, bases<ConstitutiveLaw>, boost::noncopyable>(
		"RveConstitutiveLawPlaneStressBaseType", no_init)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawPlaneStressBaseType, double>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawPlaneStressBaseType, Vector>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawPlaneStressBaseType, Matrix>)
        .def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawPlaneStressBaseType, array_1d<double,3> >)
        .def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawPlaneStressBaseType, array_1d<double,6> >)
		;

	typedef RveConstitutiveLaw<RvePlaneStressAdapterType> RveConstitutiveLawPlaneStressType;
	class_<RveConstitutiveLawPlaneStressType, RveConstitutiveLawPlaneStressType::Pointer, 
		   bases<RveConstitutiveLawPlaneStressBaseType>, 
		   boost::noncopyable>(
		"RveConstitutiveLawPlaneStress",
		init<const RvePlaneStressAdapterType::Pointer&>())
		.def("GetModelPart", &RveConstitutiveLawPlaneStressType::GetModelPart)
		.def("TestMaterialResponse", &RveConstitutiveLawPlaneStressType::TestMaterialResponse)
		;

	// 3D

	typedef RveAdapter<SparseSpaceType, LocalSpaceType, LinearSolverBaseType, RveAdapterSettings_3D> Rve3DAdapterType;
	class_<Rve3DAdapterType, Rve3DAdapterType::Pointer, boost::noncopyable>(
		"Rve3DAdapter",
		init<>())
		.def("SetRveData", &Rve3DAdapterType::SetRveData)
		.def("RveGenerated", &Rve3DAdapterType::RveGenerated)
		.def("RveGenerationRequested", &Rve3DAdapterType::RveGenerationRequested)
		.def("WorkingSpaceDimension", &Rve3DAdapterType::WorkingSpaceDimension)
		.def("GetStrainSize", &Rve3DAdapterType::GetStrainSize)
		.def(self_ns::str(self))
		.def("TestMaterialResponse", &Rve3DAdapterType::TestMaterialResponse)
		;

	typedef ConstitutiveLawAdapter<Rve3DAdapterType> RveConstitutiveLaw3DBaseType;
	class_<RveConstitutiveLaw3DBaseType, RveConstitutiveLaw3DBaseType::Pointer, bases<ConstitutiveLaw>, boost::noncopyable>(
		"RveConstitutiveLaw3DBaseType", no_init)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLaw3DBaseType, double>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLaw3DBaseType, Vector>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLaw3DBaseType, Matrix>)
        .def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLaw3DBaseType, array_1d<double,3> >)
        .def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLaw3DBaseType, array_1d<double,6> >)
		;

	typedef RveConstitutiveLaw<Rve3DAdapterType> RveConstitutiveLaw3DType;
	class_<RveConstitutiveLaw3DType, RveConstitutiveLaw3DType::Pointer, 
		   bases<RveConstitutiveLaw3DBaseType>, 
		   boost::noncopyable>(
		"RveConstitutiveLaw3D",
		init<const Rve3DAdapterType::Pointer&>())
		.def("GetModelPart", &RveConstitutiveLaw3DType::GetModelPart)
		.def("TestMaterialResponse", &RveConstitutiveLaw3DType::TestMaterialResponse)
		;

}


}

}
