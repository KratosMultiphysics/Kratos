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

#include "custom_utilities/rve_utilities_model_part.h"
#include "custom_utilities/rve_utilities_element_info.h"

#include "custom_utilities/rve_adapter_v2.h"
#include "constitutive_laws/rve_constitutive_law.h"

#include "custom_utilities/rve_macroscale_data.h"
#include "custom_utilities/rve_geometry_descriptor.h"
#include "custom_utilities/rve_linear_system_of_equations.h"
#include "custom_utilities/rve_constraint_handler.h"
#include "custom_utilities/rve_constraint_handler_zbf_sd.h"
#include "custom_utilities/rve_constraint_handler_pbf_sd.h"
#include "custom_utilities/rve_constraint_handler_wpbf_sd.h"
#include "custom_utilities/rve_constraint_handler_zbf_sd_thick_shell.h"
#include "custom_utilities/rve_constraint_handler_pbf_sd_thick_shell.h"
#include "custom_utilities/rve_constraint_handler_pbfzu_sd_thick_shell.h"
#include "custom_utilities/rve_constraint_handler_pbfzr_sd_thick_shell.h"
#include "custom_utilities/rve_constraint_handler_pbfwts_sd_thick_shell.h"
#include "custom_utilities/rve_constraint_handler_zbf_sd_thermal.h"
#include "custom_utilities/rve_constraint_handler_pbf_sd_thermal.h"
#include "custom_utilities/rve_constraint_handler_wpbf_sd_thermal.h"
#include "custom_utilities/rve_homogenizer.h"
#include "custom_utilities/rve_homogenizer_thermal.h"
#include "custom_utilities/rve_homogenizer_thick_shell.h"

#include "custom_utilities/rve_predictor_calculator.h"
#include "custom_conditions/periodic_condition_lm_2D2N.h"
#include "geometries/line_2d_2.h"

#include "custom_utilities/rve_material_database.h"
#include "custom_utilities/rve_material_database_3D.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "includes/table.h"

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

void RveGeometryDescriptor_SetUserCornerNodes_Helper(RveGeometryDescriptor& self, const boost::python::list & values)
{
	size_t n = len(values);
	RveGeometryDescriptor::IndexContainerType ids(n);
	for(size_t i = 0; i < n; i++) 
		ids[i] = boost::python::extract<double>(values[i]);
	self.SetUserCornerNodes(ids);
}




void Table_example_cpp_01()
{
	// define a table of double(key) and double(values)
	// here we use DISTANCE and PRESSURE as an example
	typedef Table<double,double> TableType;
	TableType aTable;

	// let's fill the table.
	// we will create a gaussian function, centered at DISTANCE=0.0,
	// with a half-width=20.0, and a peak-value=50.0
	// since the table is not a real function, we need to discretize the function.
	
	double a = 50.0; // height of the bell (peak-value)
	double b = 0.0; // abscissa of the peak-value
	double c = 20.0; // half width of the bell
	
	unsigned int num_samples = 10;
	double distance_increment = c / double(num_samples);
	double current_distance = 0.0;

	// add the first point
	aTable.PushBack(current_distance, a);

	// add the other points until we reach the half-distance of 20.0 units
	for(unsigned int i = 0; i < num_samples; i++)
	{
		current_distance += distance_increment;
		double val = a*std::exp(-std::pow(current_distance-b, 2)/(2.0*c));
		aTable.PushBack(current_distance, val);
	}

	// this function doesn't acutally reach the Zero-value. And since the table
	// will interpolate the results if we are at an Abscissa greater than the last one,
	// we will obtain (very small) negative values.
	// to avoid this let's add two extra points with a Zero-value.
	current_distance += distance_increment;
	aTable.PushBack(current_distance, 0.0);
	current_distance += distance_increment;
	aTable.PushBack(current_distance, 0.0);

	// one the table has been created, we can access it as follows:
	// note: we want to use a linear interpolation in this example, so we use the function GetValue(X).
	// If one doesn't want to interpolate the results, the function GetNearestValue(X) can be used.
	unsigned int num_samples_for_output = 20;
	distance_increment = c / double(num_samples_for_output);
	current_distance = 0.0;

	std::stringstream outbuffer;
	// get the first value
	outbuffer << "DISTANCE\tPRESSURE\n";
	outbuffer << current_distance << "\t" << aTable.GetValue(current_distance) << std::endl;
	// get the other points until we reach the half-distance of 20.0 units
	for(unsigned int i = 0; i < num_samples_for_output; i++)
	{
		current_distance += distance_increment;
		outbuffer << current_distance << "\t" << aTable.GetValue(current_distance) << std::endl;
	}
	// get some more points just to make sure the interpolation after we reach the bell with gives always Zero
	current_distance += distance_increment;
	outbuffer << current_distance << "\t" << aTable.GetValue(current_distance) << std::endl;
	current_distance += distance_increment;
	outbuffer << current_distance << "\t" << aTable.GetValue(current_distance) << std::endl;
	// print
	std::cout << outbuffer.str();
}



void AddPeriodicCondition(ModelPart& mp, int slave_node, int master_node)
{
	unsigned int max_cond_id = 0;
	for(ModelPart::ConditionIterator it = mp.ConditionsBegin(); it != mp.ConditionsEnd(); ++it)
		if(max_cond_id < it->GetId()) max_cond_id = it->GetId();
	max_cond_id++;
	ModelPart::NodeType::Pointer n1 = mp.pGetNode(slave_node);
	ModelPart::NodeType::Pointer n2 = mp.pGetNode(master_node);
	Condition::GeometryType::PointsArrayType points;
	points.push_back(n1);
	points.push_back(n2);
	Condition::GeometryType::Pointer new_geom( new Line2D2<Condition::GeometryType::PointType>( points ) );
	PeriodicConditionLM2D2N::Pointer new_cond( new PeriodicConditionLM2D2N(max_cond_id, new_geom) );
	mp.AddCondition( new_cond );
}



void AddUtilitiesToPython()
{

	def("Table_example_cpp_01", &Table_example_cpp_01);
	def("AddPeriodicCondition", &AddPeriodicCondition);

	def("RveCloneModelPart", &RveUtilities::CloneModelPart);
	//def("RveCloneModelPart2Physics", &RveUtilities::CloneModelPart2Physics);
	def("ReorientNarrowQuads", &RveUtilities::ReorientNarrowQuads);
	def("ReorientNarrowQuadsReplaceWithInterface", &RveUtilities::ReorientNarrowQuadsReplaceWithInterface);
	def("ReorientQuadsX", &RveUtilities::ReorientQuadsX);
	def("ReorientQuadsY", &RveUtilities::ReorientQuadsY);
	def("ReorientQuadsCCWBL", &RveUtilities::ReorientQuadsCCWBL);
	def("TestOpenMP", &TestOpenMP);

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


	typedef UblasSpace<double, CompressedMatrix, Vector>   SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector>			   LocalSpaceType;
	typedef LinearSolver<SparseSpaceType, LocalSpaceType>  LinearSolverBaseType;
	
	class_<RveMaterialDatabase, RveMaterialDatabase::Pointer,
		boost::noncopyable>(
		"RveMaterialDatabase",
		init<std::string, std::string, std::string>())
		.def(self_ns::str(self))
		;

	class_<RveMacroscaleData, RveMacroscaleData::Pointer, 
		boost::noncopyable>(
		"RveMacroscaleData",
		init<>())
		.def(self_ns::str(self))
		;

	class_<RveGeometryDescriptor, RveGeometryDescriptor::Pointer, 
		boost::noncopyable>(
		"RveGeometryDescriptor",
		init<>())
		.def("Build", &RveGeometryDescriptor::Build)
		.def("SetUserCornerNodes", &RveGeometryDescriptor_SetUserCornerNodes_Helper)
		.def(self_ns::str(self))
		;

	class_<RvePredictorCalculator, RvePredictorCalculator::Pointer,
		boost::noncopyable>(
			"RvePredictorCalculator",
			init<std::string, std::string, std::string>())
		;

	typedef RveConstraintHandler<SparseSpaceType, LocalSpaceType> RveConstraintHandlerBaseType;
	class_<RveConstraintHandlerBaseType, RveConstraintHandlerBaseType::Pointer, 
		boost::noncopyable>(
		"RveConstraintHandler",
		init<>())
		;
	typedef RveConstraintHandler_ZBF_SD<SparseSpaceType, LocalSpaceType> RveConstraintHandler_ZBF_SD_Type;
	class_<RveConstraintHandler_ZBF_SD_Type, RveConstraintHandler_ZBF_SD_Type::Pointer, 
		bases<RveConstraintHandlerBaseType>, 
		boost::noncopyable>(
		"RveConstraintHandler_ZBF_SD",
		init<>())
		;
	typedef RveConstraintHandler_PBF_SD<SparseSpaceType, LocalSpaceType> RveConstraintHandler_PBF_SD_Type;
	class_<RveConstraintHandler_PBF_SD_Type, RveConstraintHandler_PBF_SD_Type::Pointer, 
		bases<RveConstraintHandlerBaseType>, 
		boost::noncopyable>(
		"RveConstraintHandler_PBF_SD",
		init<>())
		;
	typedef RveConstraintHandler_WPBF_SD<SparseSpaceType, LocalSpaceType> RveConstraintHandler_WPBF_SD_Type;
	class_<RveConstraintHandler_WPBF_SD_Type, RveConstraintHandler_WPBF_SD_Type::Pointer, 
		bases<RveConstraintHandlerBaseType>, 
		boost::noncopyable>(
		"RveConstraintHandler_WPBF_SD",
		init<>())
		;
	typedef RveConstraintHandler_ZBF_SD_ThickShell<SparseSpaceType, LocalSpaceType> RveConstraintHandler_ZBF_SD_ThicShell_Type;
	class_<RveConstraintHandler_ZBF_SD_ThicShell_Type, RveConstraintHandler_ZBF_SD_ThicShell_Type::Pointer, 
		bases<RveConstraintHandlerBaseType>, 
		boost::noncopyable>(
		"RveConstraintHandler_ZBF_SD_ThickShell",
		init<>())
		;
	typedef RveConstraintHandler_PBF_SD_ThickShell<SparseSpaceType, LocalSpaceType> RveConstraintHandler_PBF_SD_ThicShell_Type;
	class_<RveConstraintHandler_PBF_SD_ThicShell_Type, RveConstraintHandler_PBF_SD_ThicShell_Type::Pointer, 
		bases<RveConstraintHandlerBaseType>, 
		boost::noncopyable>(
		"RveConstraintHandler_PBF_SD_ThickShell",
		init<>())
		;
	typedef RveConstraintHandler_PBFZU_SD_ThickShell<SparseSpaceType, LocalSpaceType> RveConstraintHandler_PBFZU_SD_ThicShell_Type;
	class_<RveConstraintHandler_PBFZU_SD_ThicShell_Type, RveConstraintHandler_PBFZU_SD_ThicShell_Type::Pointer, 
		bases<RveConstraintHandlerBaseType>, 
		boost::noncopyable>(
		"RveConstraintHandler_PBFZU_SD_ThickShell",
		init<>())
		;
	typedef RveConstraintHandler_PBFZR_SD_ThickShell<SparseSpaceType, LocalSpaceType> RveConstraintHandler_PBFZR_SD_ThicShell_Type;
	class_<RveConstraintHandler_PBFZR_SD_ThicShell_Type, RveConstraintHandler_PBFZR_SD_ThicShell_Type::Pointer, 
		bases<RveConstraintHandlerBaseType>, 
		boost::noncopyable>(
		"RveConstraintHandler_PBFZR_SD_ThickShell",
		init<>())
		;
	typedef RveConstraintHandler_PBFWTS_SD_ThickShell<SparseSpaceType, LocalSpaceType> RveConstraintHandler_PBFWTS_SD_ThicShell_Type;
	class_<RveConstraintHandler_PBFWTS_SD_ThicShell_Type, RveConstraintHandler_PBFWTS_SD_ThicShell_Type::Pointer, 
		bases<RveConstraintHandlerBaseType>, 
		boost::noncopyable>(
		"RveConstraintHandler_PBFWTS_SD_ThickShell",
		init<>())
		;

	typedef RveLinearSystemOfEquations<SparseSpaceType, LocalSpaceType> RveLinearSystemOfEquationsType;
	class_<RveLinearSystemOfEquationsType, RveLinearSystemOfEquationsType::Pointer, 
		boost::noncopyable>(
		"RveLinearSystemOfEquations",
		init<LinearSolverBaseType::Pointer>())
		;

	typedef RveHomogenizer<SparseSpaceType, LocalSpaceType> RveHomogenizerType;
	class_<RveHomogenizerType, RveHomogenizerType::Pointer, 
		boost::noncopyable>(
		"RveHomogenizer",
		init<>())
		;

	typedef RveHomogenizerThickShell<SparseSpaceType, LocalSpaceType> RveHomogenizerThickShellType;
	class_<RveHomogenizerThickShellType, RveHomogenizerThickShellType::Pointer, bases<RveHomogenizerType>,
		boost::noncopyable>(
		"RveHomogenizerThickShell",
		init<>())
		;

	// PLANE STRESS

	typedef RveAdapterV2<SparseSpaceType, LocalSpaceType, RveAdapterSettings_PlaneStress> RvePlaneStressAdapterV2Type;
	class_<RvePlaneStressAdapterV2Type, RvePlaneStressAdapterV2Type::Pointer, boost::noncopyable>(
		"RvePlaneStressAdapterV2",
		init<>())
		.def("SetPredictorData", &RvePlaneStressAdapterV2Type::SetPredictorData)
		.def("SetRveDataAfterPredictor", &RvePlaneStressAdapterV2Type::SetRveDataAfterPredictor)
		.def("SetRveData", &RvePlaneStressAdapterV2Type::SetRveData)
		.def("RveGenerated", &RvePlaneStressAdapterV2Type::RveGenerated)
		.def("RveGenerationRequested", &RvePlaneStressAdapterV2Type::RveGenerationRequested)
		.def("WorkingSpaceDimension", &RvePlaneStressAdapterV2Type::WorkingSpaceDimension)
		.def("GetStrainSize", &RvePlaneStressAdapterV2Type::GetStrainSize)
		.def(self_ns::str(self))
		;

	typedef ConstitutiveLawAdapter<RvePlaneStressAdapterV2Type> RveConstitutiveLawV2PlaneStressBaseType;
	class_<RveConstitutiveLawV2PlaneStressBaseType, RveConstitutiveLawV2PlaneStressBaseType::Pointer, bases<ConstitutiveLaw>, boost::noncopyable>(
		"RveConstitutiveLawV2PlaneStressBase", no_init)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2PlaneStressBaseType, double>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2PlaneStressBaseType, Vector>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2PlaneStressBaseType, Matrix>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2PlaneStressBaseType, array_1d<double,3> >)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2PlaneStressBaseType, array_1d<double,6> >)
		;

	typedef RveConstitutiveLaw<RvePlaneStressAdapterV2Type> RveConstitutiveLawV2PlaneStressType;
	class_<RveConstitutiveLawV2PlaneStressType, RveConstitutiveLawV2PlaneStressType::Pointer, 
		   bases<RveConstitutiveLawV2PlaneStressBaseType>, 
		   boost::noncopyable>(
		"RveConstitutiveLawV2PlaneStress",
		init<const RvePlaneStressAdapterV2Type::Pointer&>())
		.def("GetModelPart", &RveConstitutiveLawV2PlaneStressType::GetModelPart)
		;

	// 3D

	typedef RveAdapterV2<SparseSpaceType, LocalSpaceType, RveAdapterSettings_3D> Rve3DAdapterV2Type;
	class_<Rve3DAdapterV2Type, Rve3DAdapterV2Type::Pointer, boost::noncopyable>(
		"Rve3DAdapterV2",
		init<>())
		.def("SetPredictorData", &Rve3DAdapterV2Type::SetPredictorData)
		.def("SetRveDataAfterPredictor", &Rve3DAdapterV2Type::SetRveDataAfterPredictor)
		.def("SetRveData", &Rve3DAdapterV2Type::SetRveData)
		.def("RveGenerated", &Rve3DAdapterV2Type::RveGenerated)
		.def("RveGenerationRequested", &Rve3DAdapterV2Type::RveGenerationRequested)
		.def("WorkingSpaceDimension", &Rve3DAdapterV2Type::WorkingSpaceDimension)
		.def("GetStrainSize", &Rve3DAdapterV2Type::GetStrainSize)
		.def(self_ns::str(self))
		;

	typedef ConstitutiveLawAdapter<Rve3DAdapterV2Type> RveConstitutiveLawV23DBaseType;
	class_<RveConstitutiveLawV23DBaseType, RveConstitutiveLawV23DBaseType::Pointer, bases<ConstitutiveLaw>, boost::noncopyable>(
		"RveConstitutiveLawV23DBase", no_init)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV23DBaseType, double>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV23DBaseType, Vector>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV23DBaseType, Matrix>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV23DBaseType, array_1d<double,3> >)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV23DBaseType, array_1d<double,6> >)
		;

	typedef RveConstitutiveLaw<Rve3DAdapterV2Type> RveConstitutiveLawV23DType;
	class_<RveConstitutiveLawV23DType, RveConstitutiveLawV23DType::Pointer, 
		   bases<RveConstitutiveLawV23DBaseType>, 
		   boost::noncopyable>(
		"RveConstitutiveLawV23D",
		init<const Rve3DAdapterV2Type::Pointer&>())
		.def("GetModelPart", &RveConstitutiveLawV23DType::GetModelPart)
		;

	// ThickShell

	typedef RveAdapterV2<SparseSpaceType, LocalSpaceType, RveAdapterSettings_ThickShell> RveThickShellAdapterV2Type;
	class_<RveThickShellAdapterV2Type, RveThickShellAdapterV2Type::Pointer, boost::noncopyable>(
		"RveThickShellAdapterV2",
		init<>())
		.def("SetRveData", &RveThickShellAdapterV2Type::SetRveData)
		.def("RveGenerated", &RveThickShellAdapterV2Type::RveGenerated)
		.def("RveGenerationRequested", &RveThickShellAdapterV2Type::RveGenerationRequested)
		.def("WorkingSpaceDimension", &RveThickShellAdapterV2Type::WorkingSpaceDimension)
		.def("GetStrainSize", &RveThickShellAdapterV2Type::GetStrainSize)
		.def(self_ns::str(self))
		;

	typedef ConstitutiveLawAdapter<RveThickShellAdapterV2Type> RveConstitutiveLawV2ThickShellBaseType;
	class_<RveConstitutiveLawV2ThickShellBaseType, RveConstitutiveLawV2ThickShellBaseType::Pointer, bases<ConstitutiveLaw>, boost::noncopyable>(
		"RveConstitutiveLawV2ThickShellBase", no_init)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThickShellBaseType, double>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThickShellBaseType, Vector>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThickShellBaseType, Matrix>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThickShellBaseType, array_1d<double,3> >)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThickShellBaseType, array_1d<double,6> >)
		;

	typedef RveConstitutiveLaw<RveThickShellAdapterV2Type> RveConstitutiveLawV2ThickShellType;
	class_<RveConstitutiveLawV2ThickShellType, RveConstitutiveLawV2ThickShellType::Pointer, 
		   bases<RveConstitutiveLawV2ThickShellBaseType>, 
		   boost::noncopyable>(
		"RveConstitutiveLawV2ThickShell",
		init<const RveThickShellAdapterV2Type::Pointer&>())
		.def("GetModelPart", &RveConstitutiveLawV2ThickShellType::GetModelPart)
		;

	// ThinShell

	typedef RveAdapterV2<SparseSpaceType, LocalSpaceType, RveAdapterSettings_ThinShell> RveThinShellAdapterV2Type;
	class_<RveThinShellAdapterV2Type, RveThinShellAdapterV2Type::Pointer, boost::noncopyable>(
		"RveThinShellAdapterV2",
		init<>())
		.def("SetRveData", &RveThinShellAdapterV2Type::SetRveData)
		.def("RveGenerated", &RveThinShellAdapterV2Type::RveGenerated)
		.def("RveGenerationRequested", &RveThinShellAdapterV2Type::RveGenerationRequested)
		.def("WorkingSpaceDimension", &RveThinShellAdapterV2Type::WorkingSpaceDimension)
		.def("GetStrainSize", &RveThinShellAdapterV2Type::GetStrainSize)
		.def(self_ns::str(self))
		;

	typedef ConstitutiveLawAdapter<RveThinShellAdapterV2Type> RveConstitutiveLawV2ThinShellBaseType;
	class_<RveConstitutiveLawV2ThinShellBaseType, RveConstitutiveLawV2ThinShellBaseType::Pointer, bases<ConstitutiveLaw>, boost::noncopyable>(
		"RveConstitutiveLawV2ThinShellBase", no_init)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThinShellBaseType, double>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThinShellBaseType, Vector>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThinShellBaseType, Matrix>)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThinShellBaseType, array_1d<double,3> >)
		.def("GetValue", ConstitutiveLawGetValue<RveConstitutiveLawV2ThinShellBaseType, array_1d<double,6> >)
		;

	typedef RveConstitutiveLaw<RveThinShellAdapterV2Type> RveConstitutiveLawV2ThinShellType;
	class_<RveConstitutiveLawV2ThinShellType, RveConstitutiveLawV2ThinShellType::Pointer, 
		   bases<RveConstitutiveLawV2ThinShellBaseType>, 
		   boost::noncopyable>(
		"RveConstitutiveLawV2ThinShell",
		init<const RveThinShellAdapterV2Type::Pointer&>())
		.def("GetModelPart", &RveConstitutiveLawV2ThinShellType::GetModelPart)
		;
}


}

}
