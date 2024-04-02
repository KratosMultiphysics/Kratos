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

using SizeType = std::size_t;
typedef ConstitutiveLaw ConstitutiveLawBaseType;
template<class TVariableType> bool ConstitutiveLawHas(ConstitutiveLaw& rThisConstitutiveLaw, TVariableType const& rThisVariable) { return rThisConstitutiveLaw.Has(rThisVariable); }

//dirty trick. give back a copy instead of a reference
template<class TDataType> const TDataType ConstitutiveLawGetValue(ConstitutiveLaw& rThisConstitutiveLaw, const Variable<TDataType >& rThisVariable, TDataType& value )
{
    TDataType tmp = rThisConstitutiveLaw.GetValue(rThisVariable, value);
    return tmp;
}

// Function to export CalculateValue(...).
// Returns a copy instead of a reference, as GetValue does.
template<class TDataType> const TDataType ConstitutiveLawCalculateValue(
        ConstitutiveLaw& rThisConstitutiveLaw,
        ConstitutiveLaw::Parameters& rValues,
        const Variable<TDataType >& rThisVariable,
        TDataType& value)
{
    TDataType tmp = rThisConstitutiveLaw.CalculateValue(rValues, rThisVariable, value);
    return tmp;
}

template<class TDataType> void ConstitutiveLawSetValue(ConstitutiveLaw& rThisConstitutiveLaw, const Variable<TDataType>& rThisVariable, TDataType& value, const ProcessInfo& rCurrentProcessInfo)
{ rThisConstitutiveLaw.SetValue(rThisVariable, value, rCurrentProcessInfo); }

//only exporting functions with the new interface - considering deprecated the old interface
void NewInterfaceCalculateMaterialResponse(ConstitutiveLaw& rThisConstitutiveLaw, ConstitutiveLaw::Parameters& rValues,const ConstitutiveLaw::StressMeasure& rStressMeasure)
{rThisConstitutiveLaw.CalculateMaterialResponse (rValues,rStressMeasure);}

Flags GetFeaturesOptions(ConstitutiveLaw::Features& rThisFeatures){ return rThisFeatures.GetOptions();}
SizeType GetStrainSizeFeatures(ConstitutiveLaw::Features& rThisFeatures){ return rThisFeatures.GetStrainSize();}
SizeType GetSpaceDimensionFeatures(ConstitutiveLaw::Features& rThisFeatures){ return rThisFeatures.GetSpaceDimension();}
std::vector<ConstitutiveLaw::StrainMeasure>& GetStrainMeasuresFeatures(ConstitutiveLaw::Features& rThisFeatures){ return rThisFeatures.GetStrainMeasures(); }

Flags GetLawOptions(ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetOptions();}

double GetDeterminantF1(ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetDeterminantF();}

ConstitutiveLaw::StrainVectorType& GetStrainVector1(ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetStrainVector();}
ConstitutiveLaw::StrainVectorType& GetStrainVector2(ConstitutiveLaw::Parameters& rThisParameters, ConstitutiveLaw::StrainVectorType& strain){ return rThisParameters.GetStrainVector(strain);}
ConstitutiveLaw::StressVectorType& GetStressVector1(ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetStressVector();}
ConstitutiveLaw::StressVectorType& GetStressVector2(ConstitutiveLaw::Parameters& rThisParameters, ConstitutiveLaw::StressVectorType& stress){ return rThisParameters.GetStressVector(stress);}
ConstitutiveLaw::VoigtSizeMatrixType& GetConstitutiveMatrix1(ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetConstitutiveMatrix();}
ConstitutiveLaw::VoigtSizeMatrixType& GetConstitutiveMatrix2(ConstitutiveLaw::Parameters& rThisParameters, ConstitutiveLaw::VoigtSizeMatrixType& C){ return rThisParameters.GetConstitutiveMatrix(C);}
const ConstitutiveLaw::DeformationGradientMatrixType& GetDeformationGradientF1(ConstitutiveLaw::Parameters& rThisParameters){ return rThisParameters.GetDeformationGradientF();}
ConstitutiveLaw::DeformationGradientMatrixType& GetDeformationGradientF2(ConstitutiveLaw::Parameters& rThisParameters, ConstitutiveLaw::DeformationGradientMatrixType& F){ return rThisParameters.GetDeformationGradientF(F);}

ConstitutiveLaw::Pointer CreateWithoutProperties(ConstitutiveLaw& rThisConstitutiveLaw, Kratos::Parameters NewParameters){ return rThisConstitutiveLaw.Create(NewParameters);}
ConstitutiveLaw::Pointer CreateWithProperties(ConstitutiveLaw& rThisConstitutiveLaw, Kratos::Parameters NewParameters, const Properties& rProperties){ return rThisConstitutiveLaw.Create(NewParameters, rProperties);}

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
      .def("GetOptions",GetFeaturesOptions)
      .def("GetStrainSize",GetStrainSizeFeatures)
      .def("GetSpaceDimension",GetSpaceDimensionFeatures)
      .def("GetStrainMeasures",&GetStrainMeasuresFeatures, py::return_value_policy::reference_internal)
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
        .def("GetOptions",GetLawOptions)
        .def("GetDeterminantF",GetDeterminantF1)
//         .def("GetShapeFunctionsValues",&ConstitutiveLaw::Parameters::GetShapeFunctionsValues)
//         .def("GetShapeFunctionsDerivatives",&ConstitutiveLaw::Parameters::GetShapeFunctionsDerivatives)
        .def("GetDeformationGradientF",&GetDeformationGradientF1, py::return_value_policy::reference_internal)
        .def("GetDeformationGradientF",&GetDeformationGradientF2, py::return_value_policy::reference_internal)
        .def("GetStrainVector",&GetStrainVector1, py::return_value_policy::reference_internal)
        .def("GetStrainVector",&GetStrainVector2, py::return_value_policy::reference_internal)
        .def("GetStressVector",&GetStressVector1, py::return_value_policy::reference_internal)
        .def("GetStressVector",&GetStressVector2, py::return_value_policy::reference_internal)
        .def("GetConstitutiveMatrix",&GetConstitutiveMatrix1, py::return_value_policy::reference_internal)
        .def("GetConstitutiveMatrix",&GetConstitutiveMatrix2, py::return_value_policy::reference_internal)
        .def("GetShapeFunctionsValues",&ConstitutiveLaw::Parameters::GetShapeFunctionsValues, py::return_value_policy::reference_internal)
        .def("GetProcessInfo",&ConstitutiveLaw::Parameters::GetProcessInfo, py::return_value_policy::reference_internal)
        .def("GetMaterialProperties",&ConstitutiveLaw::Parameters::GetMaterialProperties, py::return_value_policy::reference_internal)
        .def("GetElementGeometry",&ConstitutiveLaw::Parameters::GetElementGeometry, py::return_value_policy::reference_internal)
    ;


    py::class_< ConstitutiveLaw, ConstitutiveLaw::Pointer , Flags >(m,"ConstitutiveLaw")
    .def(py::init<>() )
    .def("Create",CreateWithoutProperties)
    .def("Create",CreateWithProperties)
    .def("Clone",&ConstitutiveLaw::Clone)
    .def("WorkingSpaceDimension",&ConstitutiveLaw::WorkingSpaceDimension)
    .def("GetStrainSize",&ConstitutiveLaw::GetStrainSize)
    .def("GetStressMeasure",&ConstitutiveLaw::GetStressMeasure)
    .def("IsIncremental",&ConstitutiveLaw::IsIncremental)
    .def("WorkingSpaceDimension",&ConstitutiveLaw::WorkingSpaceDimension)
    .def("Has", &ConstitutiveLawHas< Variable<bool> >)
    .def("Has", &ConstitutiveLawHas< Variable<int> >)
    .def("Has", &ConstitutiveLawHas< Variable<double> >)
    .def("Has", &ConstitutiveLawHas< Variable<array_1d<double,3> > >)
    .def("Has", &ConstitutiveLawHas< Variable<Vector> >)
    .def("Has", &ConstitutiveLawHas< Variable<Matrix> >)
    .def("GetValue", &ConstitutiveLawGetValue<bool> )
    .def("GetValue", &ConstitutiveLawGetValue<int> )
    .def("GetValue", &ConstitutiveLawGetValue<double> )
    .def("GetValue", &ConstitutiveLawGetValue<array_1d<double,3>  >)
    .def("GetValue", &ConstitutiveLawGetValue<Vector >)
    .def("GetValue", &ConstitutiveLawGetValue<Matrix >)
    .def("SetValue", &ConstitutiveLawSetValue<int> )
    .def("SetValue", &ConstitutiveLawSetValue<double> )
    .def("SetValue", &ConstitutiveLawSetValue<array_1d<double,3>  >)
    .def("SetValue", &ConstitutiveLawSetValue<Vector >)
    .def("SetValue", &ConstitutiveLawSetValue<Matrix >)
    .def("CalculateValue", &ConstitutiveLawCalculateValue<bool> )
    .def("CalculateValue", &ConstitutiveLawCalculateValue<int> )
    .def("CalculateValue", &ConstitutiveLawCalculateValue<double> )
    .def("CalculateValue", &ConstitutiveLawCalculateValue<array_1d<double,3>  >)
    .def("CalculateValue", &ConstitutiveLawCalculateValue<Vector >)
    .def("CalculateValue", &ConstitutiveLawCalculateValue<Matrix >)
    .def("CalculateMaterialResponse",&NewInterfaceCalculateMaterialResponse)
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
