//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "add_constitutive_law_to_python.h"
#include "includes/properties.h"
#include "includes/define_python.h"
#include "includes/constitutive_law.h"
#include "containers/variable.h"

namespace Kratos::Python
{
namespace py = pybind11;

// Helper to cast for overloaded functions
using namespace pybind11::literals;

void AddConstitutiveLawToPython(pybind11::module& m)
{
    py::enum_<ConstitutiveLaw::StrainMeasure>(m,"StrainMeasure")
        .value("StrainMeasure_Infinitesimal", ConstitutiveLaw::StrainMeasure_Infinitesimal)
        .value("StrainMeasure_GreenLagrange", ConstitutiveLaw::StrainMeasure_GreenLagrange)
        .value("StrainMeasure_Hencky_Material",ConstitutiveLaw::StrainMeasure_Hencky_Material)
        .value("StrainMeasure_Hencky_Spatial",ConstitutiveLaw::StrainMeasure_Hencky_Spatial)
        .value("StrainMeasure_Deformation_Gradient",ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
        .value("StrainMeasure_Right_CauchyGreen",ConstitutiveLaw::StrainMeasure_Right_CauchyGreen)
        .value("StrainMeasure_Left_CauchyGreen",ConstitutiveLaw::StrainMeasure_Left_CauchyGreen)
        .export_values();

    py::enum_<ConstitutiveLaw::StressMeasure>(m,"StressMeasure")
        .value("StressMeasure_PK1",ConstitutiveLaw::StressMeasure_PK1)
        .value("StressMeasure_PK2", ConstitutiveLaw::StressMeasure_PK2)
        .value("StressMeasure_Kirchhoff",ConstitutiveLaw::StressMeasure_Kirchhoff)
        .value("StressMeasure_Cauchy",ConstitutiveLaw::StressMeasure_Cauchy)
        .export_values();

    py::class_< ConstitutiveLaw::Features, ConstitutiveLaw::Features::Pointer>(m,"ConstitutiveLawFeatures")
        .def(py::init<>() )
        .def("SetOptions",&ConstitutiveLaw::Features::SetOptions)
        .def("SetStrainSize",&ConstitutiveLaw::Features::SetStrainSize)
        .def("SetSpaceDimension",&ConstitutiveLaw::Features::SetSpaceDimension)
        .def("SetStrainMeasures",&ConstitutiveLaw::Features::SetStrainMeasures)
        .def("SetStrainMeasure",&ConstitutiveLaw::Features::SetStrainMeasure)
        .def("GetOptions", &ConstitutiveLaw::Features::GetOptions)
        .def("GetStrainSize", &ConstitutiveLaw::Features::GetStrainSize)
        .def("GetSpaceDimension", &ConstitutiveLaw::Features::GetSpaceDimension)
        .def("GetStrainMeasures", &ConstitutiveLaw::Features::GetStrainMeasures, py::return_value_policy::reference_internal)
        ;

py::class_< ConstitutiveLaw::Parameters, ConstitutiveLaw::Parameters::Pointer>(m,"ConstitutiveLawParameters")
        .def(py::init<>())
        .def(py::init<const ConstitutiveLaw::GeometryType&, const Properties&, const ConstitutiveLaw::ProcessInfoType&>())
        .def("CheckAllParameters", &ConstitutiveLaw::Parameters::CheckAllParameters)
        .def("CheckMechanicalVariables", &ConstitutiveLaw::Parameters::CheckMechanicalVariables)
        .def("CheckShapeFunctions", &ConstitutiveLaw::Parameters::CheckShapeFunctions)
        .def("CheckInfoMaterialGeometry", &ConstitutiveLaw::Parameters::CheckInfoMaterialGeometry)
        .def("Set", &ConstitutiveLaw::Parameters::Set)
        .def("Reset", &ConstitutiveLaw::Parameters::Reset)
        .def("SetOptions", &ConstitutiveLaw::Parameters::SetOptions)
        .def("SetDeterminantF", &ConstitutiveLaw::Parameters::SetDeterminantF)
        .def("SetShapeFunctionsValues", &ConstitutiveLaw::Parameters::SetShapeFunctionsValues)
        .def("SetShapeFunctionsDerivatives", &ConstitutiveLaw::Parameters::SetShapeFunctionsDerivatives)
        .def("SetStrainVector", &ConstitutiveLaw::Parameters::SetStrainVector)
        .def("SetStressVector", &ConstitutiveLaw::Parameters::SetStressVector)
        .def("SetConstitutiveMatrix", &ConstitutiveLaw::Parameters::SetConstitutiveMatrix)
        .def("SetProcessInfo", &ConstitutiveLaw::Parameters::SetProcessInfo)
        .def("SetMaterialProperties", &ConstitutiveLaw::Parameters::SetMaterialProperties)
        .def("SetElementGeometry", &ConstitutiveLaw::Parameters::SetElementGeometry)
        .def("SetDeformationGradientF", &ConstitutiveLaw::Parameters::SetDeformationGradientF)
        .def("GetOptions", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetOptions), py::return_value_policy::reference_internal)
        .def("GetDeterminantF", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetDeterminantF))
        .def("GetDeformationGradientF", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetDeformationGradientF), py::return_value_policy::reference_internal)
        .def("GetStrainVector", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetStrainVector), py::return_value_policy::reference_internal)
        .def("GetStressVector", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetStressVector), py::return_value_policy::reference_internal)
        .def("GetConstitutiveMatrix", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetConstitutiveMatrix), py::return_value_policy::reference_internal)
        .def("GetShapeFunctionsValues", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetShapeFunctionsValues), py::return_value_policy::reference_internal)
        .def("GetProcessInfo", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetProcessInfo), py::return_value_policy::reference_internal)
        .def("GetMaterialProperties", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetMaterialProperties), py::return_value_policy::reference_internal)
        .def("GetElementGeometry", py::overload_cast<>(&ConstitutiveLaw::Parameters::GetElementGeometry), py::return_value_policy::reference_internal)
        ;


    py::class_< ConstitutiveLaw, ConstitutiveLaw::Pointer , Flags >(m,"ConstitutiveLaw")
        .def(py::init<>())
        .def("Clone", &ConstitutiveLaw::Clone)
        .def("Create", static_cast<ConstitutiveLaw::Pointer(ConstitutiveLaw::*)(Kratos::Parameters) const>(&ConstitutiveLaw::Create))
        .def("Create", static_cast<ConstitutiveLaw::Pointer(ConstitutiveLaw::*)(Kratos::Parameters, const Properties&) const>(&ConstitutiveLaw::Create))
        .def("WorkingSpaceDimension", &ConstitutiveLaw::WorkingSpaceDimension)
        .def("GetStrainSize", &ConstitutiveLaw::GetStrainSize)
        .def("GetStressMeasure", &ConstitutiveLaw::GetStressMeasure)
        .def("IsIncremental", &ConstitutiveLaw::IsIncremental)
        .def("Has", static_cast<bool (ConstitutiveLaw::*)(const Variable<bool>&)>(&ConstitutiveLaw::Has))
        .def("Has", static_cast<bool (ConstitutiveLaw::*)(const Variable<int>&)>(&ConstitutiveLaw::Has))
        .def("Has", static_cast<bool (ConstitutiveLaw::*)(const Variable<double>&)>(&ConstitutiveLaw::Has))
        .def("Has", static_cast<bool (ConstitutiveLaw::*)(const Variable<array_1d<double, 3>>&)>(&ConstitutiveLaw::Has))
        .def("Has", static_cast<bool (ConstitutiveLaw::*)(const Variable<Vector>&)>(&ConstitutiveLaw::Has))
        .def("Has", static_cast<bool (ConstitutiveLaw::*)(const Variable<Matrix>&)>(&ConstitutiveLaw::Has))
        .def("GetValue", [](ConstitutiveLaw& rCL, const Variable<bool>& rVar) { bool v; rCL.GetValue(rVar, v); return v; })
        .def("GetValue", [](ConstitutiveLaw& rCL, const Variable<int>& rVar) { int v; rCL.GetValue(rVar, v); return v; })
        .def("GetValue", [](ConstitutiveLaw& rCL, const Variable<double>& rVar) { double v; rCL.GetValue(rVar, v); return v; })
        .def("GetValue", [](ConstitutiveLaw& rCL, const Variable<array_1d<double, 3>>& rVar) { array_1d<double, 3> v; rCL.GetValue(rVar, v); return v; })
        .def("GetValue", [](ConstitutiveLaw& rCL, const Variable<Vector>& rVar) { Vector v; rCL.GetValue(rVar, v); return v; })
        .def("GetValue", [](ConstitutiveLaw& rCL, const Variable<Matrix>& rVar) { Matrix v; rCL.GetValue(rVar, v); return v; })
        .def("SetValue", py::overload_cast<const Variable<int>&, const int&, const ProcessInfo&>(&ConstitutiveLaw::SetValue))
        .def("SetValue", py::overload_cast<const Variable<double>&, const double&, const ProcessInfo&>(&ConstitutiveLaw::SetValue))
        .def("SetValue", py::overload_cast<const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&, const ProcessInfo&>(&ConstitutiveLaw::SetValue))
        .def("SetValue", py::overload_cast<const Variable<Vector>&, const Vector&, const ProcessInfo&>(&ConstitutiveLaw::SetValue))
        .def("SetValue", py::overload_cast<const Variable<Matrix>&, const Matrix&, const ProcessInfo&>(&ConstitutiveLaw::SetValue))
        .def("CalculateValue", [](ConstitutiveLaw& rCL, ConstitutiveLaw::Parameters& rP, const Variable<bool>& rVar) { bool v; rCL.CalculateValue(rP, rVar, v); return v; })
        .def("CalculateValue", [](ConstitutiveLaw& rCL, ConstitutiveLaw::Parameters& rP, const Variable<int>& rVar) { int v; rCL.CalculateValue(rP, rVar, v); return v; })
        .def("CalculateValue", [](ConstitutiveLaw& rCL, ConstitutiveLaw::Parameters& rP, const Variable<double>& rVar) { double v; rCL.CalculateValue(rP, rVar, v); return v; })
        .def("CalculateValue", [](ConstitutiveLaw& rCL, ConstitutiveLaw::Parameters& rP, const Variable<array_1d<double, 3>>& rVar) { array_1d<double, 3> v; rCL.CalculateValue(rP, rVar, v); return v; })
        .def("CalculateValue", [](ConstitutiveLaw& rCL, ConstitutiveLaw::Parameters& rP, const Variable<Vector>& rVar) { Vector v; rCL.CalculateValue(rP, rVar, v); return v; })
        .def("CalculateValue", [](ConstitutiveLaw& rCL, ConstitutiveLaw::Parameters& rP, const Variable<Matrix>& rVar) { Matrix v; rCL.CalculateValue(rP, rVar, v); return v; })
        .def("CalculateMaterialResponse", &ConstitutiveLaw::CalculateMaterialResponse)
        .def("CalculateMaterialResponsePK1", &ConstitutiveLaw::CalculateMaterialResponsePK1)
        .def("CalculateMaterialResponsePK2", &ConstitutiveLaw::CalculateMaterialResponsePK2)
        .def("CalculateMaterialResponseKirchhoff", &ConstitutiveLaw::CalculateMaterialResponseKirchhoff)
        .def("CalculateMaterialResponseCauchy", &ConstitutiveLaw::CalculateMaterialResponseCauchy)
        .def("InitializeMaterialResponse", &ConstitutiveLaw::InitializeMaterialResponse)
        .def("InitializeMaterialResponsePK1", &ConstitutiveLaw::InitializeMaterialResponsePK1)
        .def("InitializeMaterialResponsePK2", &ConstitutiveLaw::InitializeMaterialResponsePK2)
        .def("InitializeMaterialResponseKirchhoff", &ConstitutiveLaw::InitializeMaterialResponseKirchhoff)
        .def("InitializeMaterialResponseCauchy", &ConstitutiveLaw::InitializeMaterialResponseCauchy)
        .def("FinalizeMaterialResponse", &ConstitutiveLaw::FinalizeMaterialResponse)
        .def("FinalizeMaterialResponsePK1", &ConstitutiveLaw::FinalizeMaterialResponsePK1)
        .def("FinalizeMaterialResponsePK2", &ConstitutiveLaw::FinalizeMaterialResponsePK2)
        .def("FinalizeMaterialResponseKirchhoff", &ConstitutiveLaw::FinalizeMaterialResponseKirchhoff)
        .def("FinalizeMaterialResponseCauchy", &ConstitutiveLaw::FinalizeMaterialResponseCauchy)
        .def("InitializeMaterial", &ConstitutiveLaw::InitializeMaterial)
        .def("ResetMaterial", &ConstitutiveLaw::ResetMaterial)
        .def("TransformStrains", &ConstitutiveLaw::TransformStrains, py::return_value_policy::reference_internal)
        .def("TransformPK1Stresses", &ConstitutiveLaw::TransformPK1Stresses, py::return_value_policy::reference_internal)
        .def("TransformPK2Stresses", &ConstitutiveLaw::TransformPK2Stresses, py::return_value_policy::reference_internal)
        .def("TransformKirchhoffStresses", &ConstitutiveLaw::TransformKirchhoffStresses, py::return_value_policy::reference_internal)
        .def("TransformCauchyStresses", &ConstitutiveLaw::TransformCauchyStresses, py::return_value_policy::reference_internal)
        .def("PullBackConstitutiveMatrix", &ConstitutiveLaw::PullBackConstitutiveMatrix)
        .def("PushForwardConstitutiveMatrix", &ConstitutiveLaw::PushForwardConstitutiveMatrix)
        .def("Check", &ConstitutiveLaw::Check)
        .def("GetLawFeatures", &ConstitutiveLaw::GetLawFeatures)
        .def("Info", [](const ConstitutiveLaw& rThisConstitutiveLaw){return rThisConstitutiveLaw.Info();})
        .def_readonly_static("USE_ELEMENT_PROVIDED_STRAIN", &ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)
        .def_readonly_static("COMPUTE_STRESS", &ConstitutiveLaw::COMPUTE_STRESS)
        .def_readonly_static("COMPUTE_CONSTITUTIVE_TENSOR", &ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)
        .def_readonly_static("COMPUTE_STRAIN_ENERGY", &ConstitutiveLaw::COMPUTE_STRAIN_ENERGY)
        .def_readonly_static("ISOCHORIC_TENSOR_ONLY", &ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)
        .def_readonly_static("VOLUMETRIC_TENSOR_ONLY", &ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)
        .def_readonly_static("MECHANICAL_RESPONSE_ONLY", &ConstitutiveLaw::MECHANICAL_RESPONSE_ONLY)
        .def_readonly_static("THERMAL_RESPONSE_ONLY", &ConstitutiveLaw::THERMAL_RESPONSE_ONLY)
        .def_readonly_static("INCREMENTAL_STRAIN_MEASURE", &ConstitutiveLaw::INCREMENTAL_STRAIN_MEASURE)
        .def_readonly_static("INITIALIZE_MATERIAL_RESPONSE", &ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE)
        .def_readonly_static("FINALIZE_MATERIAL_RESPONSE", &ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE)
        .def_readonly_static("FINITE_STRAINS", &ConstitutiveLaw::FINITE_STRAINS)
        .def_readonly_static("INFINITESIMAL_STRAINS", &ConstitutiveLaw::INFINITESIMAL_STRAINS)
        .def_readonly_static("THREE_DIMENSIONAL_LAW", &ConstitutiveLaw::THREE_DIMENSIONAL_LAW)
        .def_readonly_static("PLANE_STRAIN_LAW", &ConstitutiveLaw::PLANE_STRAIN_LAW)
        .def_readonly_static("PLANE_STRESS_LAW", &ConstitutiveLaw::PLANE_STRESS_LAW)
        .def_readonly_static("AXISYMMETRIC_LAW", &ConstitutiveLaw::AXISYMMETRIC_LAW)
        .def_readonly_static("U_P_LAW", &ConstitutiveLaw::U_P_LAW)
        .def_readonly_static("ISOTROPIC", &ConstitutiveLaw::ISOTROPIC)
        .def_readonly_static("ANISOTROPIC", &ConstitutiveLaw::ANISOTROPIC)
        ;

}
}  // namespace Kratos::Python.