//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
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

void  AddConstitutiveLawToPython(pybind11::module& m)
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
    ;

    py::enum_<ConstitutiveLaw::StressMeasure>(m,"StressMeasure")
        .value("StressMeasure_PK1",ConstitutiveLaw::StressMeasure_PK1)
        .value("StressMeasure_PK2", ConstitutiveLaw::StressMeasure_PK2)
        .value("StressMeasure_Kirchhoff",ConstitutiveLaw::StressMeasure_Kirchhoff)
        .value("StressMeasure_Cauchy",ConstitutiveLaw::StressMeasure_Cauchy)
        .export_values();
    ;

    py::class_< ConstitutiveLaw::Features, ConstitutiveLaw::Features::Pointer>(m,"ConstitutiveLawFeatures")
        .def(py::init<>() )
      .def("SetOptions",&ConstitutiveLaw::Features::SetOptions)
      .def("SetStrainSize",&ConstitutiveLaw::Features::SetStrainSize)
      .def("SetSpaceDimension",&ConstitutiveLaw::Features::SetSpaceDimension)
      .def("SetStrainMeasures",&ConstitutiveLaw::Features::SetStrainMeasures)
      .def("SetStrainMeasure",&ConstitutiveLaw::Features::SetStrainMeasure)
      .def("GetOptions",[](ConstitutiveLaw::Features& rThisFeatures){ return rThisFeatures.GetOptions();})
      .def("GetStrainSize",[](ConstitutiveLaw::Features& rThisFeatures){ return rThisFeatures.GetStrainSize();})
      .def("GetSpaceDimension",[](ConstitutiveLaw::Features& rThisFeatures){ return rThisFeatures.GetSpaceDimension();})
      .def("GetStrainMeasures",[](ConstitutiveLaw::Features& rThisFeatures){ return rThisFeatures.GetStrainMeasures(); }, py::return_value_policy::reference_internal)
      ;

    py::class_< ConstitutiveLaw::Parameters, ConstitutiveLaw::Parameters::Pointer>(m,"ConstitutiveLawParameters")
        .def(py::init<>() )
        .def(py::init< const ConstitutiveLaw::GeometryType& ,const Properties&, const ConstitutiveLaw::ProcessInfoType& >() )
        .def("CheckAllParameters",&ConstitutiveLaw::Parameters::CheckAllParameters)
        .def("CheckMechanicalVariables",&ConstitutiveLaw::Parameters::CheckMechanicalVariables)
        .def("CheckShapeFunctions",&ConstitutiveLaw::Parameters::CheckShapeFunctions)
        .def("CheckInfoMaterialGeometry",&ConstitutiveLaw::Parameters::CheckInfoMaterialGeometry)
        .def("Set",&ConstitutiveLaw::Parameters::Set)
        .def("Reset",&ConstitutiveLaw::Parameters::Reset)
        .def("SetOptions",&ConstitutiveLaw::Parameters::SetOptions)
        .def("SetDeterminantF",&ConstitutiveLaw::Parameters::SetDeterminantF)
        .def("SetShapeFunctionsValues",&ConstitutiveLaw::Parameters::SetShapeFunctionsValues)
        .def("SetShapeFunctionsDerivatives",&ConstitutiveLaw::Parameters::SetShapeFunctionsDerivatives)
        .def("SetStrainVector",&ConstitutiveLaw::Parameters::SetStrainVector)
        .def("SetStressVector",&ConstitutiveLaw::Parameters::SetStressVector)
        .def("SetConstitutiveMatrix",&ConstitutiveLaw::Parameters::SetConstitutiveMatrix)
        .def("SetProcessInfo",&ConstitutiveLaw::Parameters::SetProcessInfo)
        .def("SetMaterialProperties",&ConstitutiveLaw::Parameters::SetMaterialProperties)
        .def("SetElementGeometry",&ConstitutiveLaw::Parameters::SetElementGeometry)
        .def("SetDeformationGradientF",&ConstitutiveLaw::Parameters::SetDeformationGradientF)
        .def("GetOptions",[](ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetOptions();})
        .def("GetDeterminantF",[](ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetDeterminantF();})
//         .def("GetShapeFunctionsValues",&ConstitutiveLaw::Parameters::GetShapeFunctionsValues)
//         .def("GetShapeFunctionsDerivatives",&ConstitutiveLaw::Parameters::GetShapeFunctionsDerivatives)
        .def("GetDeformationGradientF",[](ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetDeformationGradientF();}, py::return_value_policy::reference_internal)
        .def("GetDeformationGradientF",[](ConstitutiveLaw::Parameters& rThisParameters, ConstitutiveLaw::DeformationGradientMatrixType& F){ return rThisParameters.GetDeformationGradientF(F);}, py::return_value_policy::reference_internal)
        .def("GetStrainVector",[](ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetStrainVector();}, py::return_value_policy::reference_internal)
        .def("GetStrainVector",[](ConstitutiveLaw::Parameters& rThisParameters, ConstitutiveLaw::StrainVectorType& strain){ return rThisParameters.GetStrainVector(strain);}, py::return_value_policy::reference_internal)
        .def("GetStressVector",[](ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetStressVector();}, py::return_value_policy::reference_internal)
        .def("GetStressVector",[](ConstitutiveLaw::Parameters& rThisParameters, ConstitutiveLaw::StressVectorType& stress){ return rThisParameters.GetStressVector(stress);}, py::return_value_policy::reference_internal)
        .def("GetConstitutiveMatrix",[](ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetConstitutiveMatrix();}, py::return_value_policy::reference_internal)
        .def("GetConstitutiveMatrix",[](ConstitutiveLaw::Parameters& rThisParameters, ConstitutiveLaw::VoigtSizeMatrixType& C){ return rThisParameters.GetConstitutiveMatrix(C);}, py::return_value_policy::reference_internal)
        .def("GetShapeFunctionsValues",&ConstitutiveLaw::Parameters::GetShapeFunctionsValues, py::return_value_policy::reference_internal)
        .def("GetProcessInfo",&ConstitutiveLaw::Parameters::GetProcessInfo, py::return_value_policy::reference_internal)
        .def("GetMaterialProperties",&ConstitutiveLaw::Parameters::GetMaterialProperties, py::return_value_policy::reference_internal)
        .def("GetElementGeometry",&ConstitutiveLaw::Parameters::GetElementGeometry, py::return_value_policy::reference_internal)
    ;


    py::class_< ConstitutiveLaw, ConstitutiveLaw::Pointer , Flags >(m,"ConstitutiveLaw")
    .def(py::init<>() )
    .def("Create",[](ConstitutiveLaw& rThisConstitutiveLaw, Kratos::Parameters NewParameters){ return rThisConstitutiveLaw.Create(NewParameters);})
    .def("Create",[](ConstitutiveLaw& rThisConstitutiveLaw, Kratos::Parameters NewParameters, const Properties& rProperties){ return rThisConstitutiveLaw.Create(NewParameters, rProperties);})
    .def("Clone",&ConstitutiveLaw::Clone)
    .def("WorkingSpaceDimension",&ConstitutiveLaw::WorkingSpaceDimension)
    .def("GetStrainSize",&ConstitutiveLaw::GetStrainSize)
    .def("GetStressMeasure",&ConstitutiveLaw::GetStressMeasure)
    .def("IsIncremental",&ConstitutiveLaw::IsIncremental)
    .def("WorkingSpaceDimension",&ConstitutiveLaw::WorkingSpaceDimension)
    .def("Has", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<bool>& rThisVariable){ return rThisConstitutiveLaw.Has(rThisVariable); })
    .def("Has", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<int>& rThisVariable){ return rThisConstitutiveLaw.Has(rThisVariable); })
    .def("Has", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<double>& rThisVariable){ return rThisConstitutiveLaw.Has(rThisVariable); })
    .def("Has", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<array_1d<double,3>>& rThisVariable){ return rThisConstitutiveLaw.Has(rThisVariable); })
    .def("Has", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<Vector>& rThisVariable){ return rThisConstitutiveLaw.Has(rThisVariable); })
    .def("Has", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<Matrix>& rThisVariable){ return rThisConstitutiveLaw.Has(rThisVariable); })
    .def("GetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<bool>& rThisVariable, bool& value ){ bool tmp = rThisConstitutiveLaw.GetValue(rThisVariable, value); return tmp;} )
    .def("GetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<int>& rThisVariable, int& value ){ int tmp = rThisConstitutiveLaw.GetValue(rThisVariable, value); return tmp;} )
    .def("GetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<double>& rThisVariable, double& value ){ double tmp = rThisConstitutiveLaw.GetValue(rThisVariable, value); return tmp;} )
    .def("GetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<array_1d<double,3>>& rThisVariable, array_1d<double,3>& value ){ array_1d<double,3> tmp = rThisConstitutiveLaw.GetValue(rThisVariable, value); return tmp;}  )
    .def("GetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<Vector>& rThisVariable, Vector& value ){ Vector tmp = rThisConstitutiveLaw.GetValue(rThisVariable, value); return tmp;} )
    .def("GetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<Matrix>& rThisVariable, Matrix& value ){ Matrix tmp = rThisConstitutiveLaw.GetValue(rThisVariable, value); return tmp;} )
    .def("SetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<int>& rThisVariable, const int& value, const ProcessInfo& rCurrentProcessInfo){ rThisConstitutiveLaw.SetValue(rThisVariable, value, rCurrentProcessInfo); })
    .def("SetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<double>& rThisVariable, const double& value, const ProcessInfo& rCurrentProcessInfo){ rThisConstitutiveLaw.SetValue(rThisVariable, value, rCurrentProcessInfo); })
    .def("SetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<array_1d<double,3>>& rThisVariable, const array_1d<double,3>& value, const ProcessInfo& rCurrentProcessInfo){ rThisConstitutiveLaw.SetValue(rThisVariable, value, rCurrentProcessInfo); }  )
    .def("SetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<Vector>& rThisVariable, const Vector& value, const ProcessInfo& rCurrentProcessInfo){ rThisConstitutiveLaw.SetValue(rThisVariable, value, rCurrentProcessInfo); })
    .def("SetValue", [](ConstitutiveLaw& rThisConstitutiveLaw, const Variable<Matrix>& rThisVariable, const Matrix& value, const ProcessInfo& rCurrentProcessInfo){ rThisConstitutiveLaw.SetValue(rThisVariable, value, rCurrentProcessInfo); })
    .def("CalculateValue", [](ConstitutiveLaw& rThisConstitutiveLaw, ConstitutiveLaw::Parameters& rValues, const Variable<bool>& rThisVariable, bool& value){ bool tmp = rThisConstitutiveLaw.CalculateValue(rValues, rThisVariable, value); return tmp;} )
    .def("CalculateValue", [](ConstitutiveLaw& rThisConstitutiveLaw, ConstitutiveLaw::Parameters& rValues, const Variable<int>& rThisVariable, int& value){ int tmp = rThisConstitutiveLaw.CalculateValue(rValues, rThisVariable, value); return tmp;} )
    .def("CalculateValue", [](ConstitutiveLaw& rThisConstitutiveLaw, ConstitutiveLaw::Parameters& rValues, const Variable<double>& rThisVariable, double& value){ double tmp = rThisConstitutiveLaw.CalculateValue(rValues, rThisVariable, value); return tmp;} )
    .def("CalculateValue", [](ConstitutiveLaw& rThisConstitutiveLaw, ConstitutiveLaw::Parameters& rValues, const Variable<array_1d<double,3>>& rThisVariable, array_1d<double,3>& value){ array_1d<double,3> tmp = rThisConstitutiveLaw.CalculateValue(rValues, rThisVariable, value); return tmp;}  )
    .def("CalculateValue", [](ConstitutiveLaw& rThisConstitutiveLaw, ConstitutiveLaw::Parameters& rValues, const Variable<Vector>& rThisVariable, Vector& value){ Vector tmp = rThisConstitutiveLaw.CalculateValue(rValues, rThisVariable, value); return tmp;} )
    .def("CalculateValue", [](ConstitutiveLaw& rThisConstitutiveLaw, ConstitutiveLaw::Parameters& rValues, const Variable<Matrix>& rThisVariable, Matrix& value){ Matrix tmp = rThisConstitutiveLaw.CalculateValue(rValues, rThisVariable, value); return tmp;} )
    .def("CalculateMaterialResponse",[](ConstitutiveLaw& rThisConstitutiveLaw, ConstitutiveLaw::Parameters& rValues,const ConstitutiveLaw::StressMeasure& rStressMeasure){rThisConstitutiveLaw.CalculateMaterialResponse (rValues,rStressMeasure);})
    .def("CalculateMaterialResponsePK1",&ConstitutiveLaw::CalculateMaterialResponsePK1)
    .def("CalculateMaterialResponsePK2",&ConstitutiveLaw::CalculateMaterialResponsePK2)
    .def("CalculateMaterialResponseKirchhoff",&ConstitutiveLaw::CalculateMaterialResponseKirchhoff)
    .def("CalculateMaterialResponseCauchy",&ConstitutiveLaw::CalculateMaterialResponseCauchy)
    .def("InitializeMaterialResponse",&ConstitutiveLaw::InitializeMaterialResponse)
    .def("InitializeMaterialResponsePK1",&ConstitutiveLaw::InitializeMaterialResponsePK1)
    .def("InitializeMaterialResponsePK2",&ConstitutiveLaw::InitializeMaterialResponsePK2)
    .def("InitializeMaterialResponseKirchhoff",&ConstitutiveLaw::InitializeMaterialResponseKirchhoff)
    .def("InitializeMaterialResponseCauchy",&ConstitutiveLaw::InitializeMaterialResponseCauchy)
    .def("FinalizeMaterialResponse",&ConstitutiveLaw::FinalizeMaterialResponse)
    .def("FinalizeMaterialResponsePK1",&ConstitutiveLaw::FinalizeMaterialResponsePK1)
    .def("FinalizeMaterialResponsePK2",&ConstitutiveLaw::FinalizeMaterialResponsePK2)
    .def("FinalizeMaterialResponseKirchhoff",&ConstitutiveLaw::FinalizeMaterialResponseKirchhoff)
    .def("FinalizeMaterialResponseCauchy",&ConstitutiveLaw::FinalizeMaterialResponseCauchy)
    .def("InitializeMaterial",&ConstitutiveLaw::InitializeMaterial)
    .def("ResetMaterial",&ConstitutiveLaw::ResetMaterial)
    .def("TransformStrains",&ConstitutiveLaw::TransformStrains, py::return_value_policy::reference_internal)
//     .def("TransformStresses",&ConstitutiveLaw::TransformStresses)
//     .def("TransformStresses",&ConstitutiveLaw::TransformStresses)
    .def("TransformPK1Stresses",&ConstitutiveLaw::TransformPK1Stresses, py::return_value_policy::reference_internal)
    .def("TransformPK2Stresses",&ConstitutiveLaw::TransformPK2Stresses, py::return_value_policy::reference_internal)
    .def("TransformKirchhoffStresses",&ConstitutiveLaw::TransformKirchhoffStresses, py::return_value_policy::reference_internal)
    .def("TransformCauchyStresses",&ConstitutiveLaw::TransformCauchyStresses, py::return_value_policy::reference_internal)
    .def("PullBackConstitutiveMatrix",&ConstitutiveLaw::PullBackConstitutiveMatrix)
    .def("PushForwardConstitutiveMatrix",&ConstitutiveLaw::PushForwardConstitutiveMatrix)
    .def("Check",&ConstitutiveLaw::Check)
    .def("GetLawFeatures",&ConstitutiveLaw::GetLawFeatures)
    .def("Info", [](ConstitutiveLaw& rThisConstitutiveLaw){return rThisConstitutiveLaw.Info();})
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