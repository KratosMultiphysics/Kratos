/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-04-24 10:30:22 $
//   Revision:            $Revision: 1.3 $
//
//


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "python/add_constitutive_law_to_python.h"

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "add_constitutive_law_to_python.h"
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
//#include "python/add_mesh_to_python.h"
//#include "python/pointer_vector_set_python_interface.h"
//#include "python/variable_indexing_python.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

typedef ConstitutiveLaw ConstitutiveLawBaseType;
typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;

template<class TVariableType> bool ConstitutiveLawHas(ConstitutiveLaw& this_constitutive_law, TVariableType const& rThisVariable) { return this_constitutive_law.Has(rThisVariable); }

//dirty trick. give back a copy instead of a reference
template<class TDataType> const TDataType ConstitutiveLawGetValue(ConstitutiveLaw& this_constitutive_law, const Variable<TDataType >& rThisVariable, TDataType& value ) 
{ 
    TDataType tmp = this_constitutive_law.GetValue(rThisVariable, value);
    return tmp;
}
template<class TDataType> void ConstitutiveLawSetValue(ConstitutiveLaw& this_constitutive_law, const Variable<TDataType>& rThisVariable, TDataType& value, const ProcessInfo& rCurrentProcessInfo)
{ this_constitutive_law.SetValue(rThisVariable, value, rCurrentProcessInfo); }

//only exporting functions with the new interface - considering deprecated the old interface
void NewInterfaceCalculateMaterialResponse(ConstitutiveLaw& this_constitutive_law, ConstitutiveLaw::Parameters& rValues,const ConstitutiveLaw::StressMeasure& rStressMeasure)
{this_constitutive_law.CalculateMaterialResponse (rValues,rStressMeasure);}

Vector& GetStrainVector1(ConstitutiveLaw::Parameters& this_params){ return this_params.GetStrainVector();}
Vector& GetStrainVector2(ConstitutiveLaw::Parameters& this_params, Vector& strain){ return this_params.GetStrainVector(strain);}
Vector& GetStressVector1(ConstitutiveLaw::Parameters& this_params){ return this_params.GetStrainVector();}
Vector& GetStressVector2(ConstitutiveLaw::Parameters& this_params, Vector& stress){ return this_params.GetStressVector(stress);}
Matrix& GetConstitutiveMatrix1(ConstitutiveLaw::Parameters& this_params){ return this_params.GetConstitutiveMatrix();}
Matrix& GetConstitutiveMatrix2(ConstitutiveLaw::Parameters& this_params, Matrix& C){ return this_params.GetConstitutiveMatrix(C);}
const Matrix& GetDeformationGradientF1(ConstitutiveLaw::Parameters& this_params){ return this_params.GetDeformationGradientF();}
Matrix& GetDeformationGradientF2(ConstitutiveLaw::Parameters& this_params, Matrix& F){ return this_params.GetDeformationGradientF(F);}



void  AddConstitutiveLawToPython()
{
    
    enum_<ConstitutiveLaw::StrainMeasure>("StrainMeasure")
        .value("StrainMeasure_Infinitesimal", ConstitutiveLaw::StrainMeasure_Infinitesimal)
        .value("StrainMeasure_GreenLagrange", ConstitutiveLaw::StrainMeasure_GreenLagrange)
        .value("StrainMeasure_Hencky_Material",ConstitutiveLaw::StrainMeasure_Hencky_Material)
        .value("StrainMeasure_Hencky_Spatial",ConstitutiveLaw::StrainMeasure_Hencky_Spatial)
        .value("StrainMeasure_Deformation_Gradient",ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
        .value("StrainMeasure_Right_CauchyGreen",ConstitutiveLaw::StrainMeasure_Right_CauchyGreen)
        .value("StrainMeasure_Left_CauchyGreen",ConstitutiveLaw::StrainMeasure_Left_CauchyGreen)
    ;

    enum_<ConstitutiveLaw::StressMeasure>("StressMeasure")
        .value("StressMeasure_PK1",ConstitutiveLaw::StressMeasure_PK1)
        .value("StressMeasure_PK2", ConstitutiveLaw::StressMeasure_PK2)
        .value("StressMeasure_Kirchhoff",ConstitutiveLaw::StressMeasure_Kirchhoff)
        .value("StressMeasure_Cauchy",ConstitutiveLaw::StressMeasure_Cauchy)
    ;


    class_< ConstitutiveLaw::Parameters, ConstitutiveLaw::Parameters::Pointer, boost::noncopyable>
    ("ConstitutiveLawParameters", init<>() )
        .def(init< const ConstitutiveLaw::GeometryType& ,const Properties&, const ConstitutiveLaw::ProcessInfoType& >() )
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
        .def("SetShapeFunctionsValues",&ConstitutiveLaw::Parameters::SetShapeFunctionsValues)
        .def("SetProcessInfo",&ConstitutiveLaw::Parameters::SetProcessInfo)
        .def("SetMaterialProperties",&ConstitutiveLaw::Parameters::SetMaterialProperties)
        .def("SetElementGeometry",&ConstitutiveLaw::Parameters::SetElementGeometry)
//         .def("GetDeterminantF" , &ConstitutiveLaw::Parameters::GetDeterminantF)
//         .def("GetShapeFunctionsValues",&ConstitutiveLaw::Parameters::GetShapeFunctionsValues)
//         .def("GetShapeFunctionsDerivatives",&ConstitutiveLaw::Parameters::GetShapeFunctionsDerivatives)
        .def("GetDeformationGradientF",&GetDeformationGradientF1, return_internal_reference<>())
        .def("GetDeformationGradientF",&GetDeformationGradientF2, return_internal_reference<>())
        .def("GetStrainVector",&GetStrainVector1, return_internal_reference<>())
        .def("GetStrainVector",&GetStrainVector2, return_internal_reference<>())
        .def("GetStressVector",&GetStressVector1, return_internal_reference<>())
        .def("GetStressVector",&GetStressVector2, return_internal_reference<>())
        .def("GetConstitutiveMatrix",&GetConstitutiveMatrix1, return_internal_reference<>())
        .def("GetConstitutiveMatrix",&GetConstitutiveMatrix2, return_internal_reference<>())
        .def("GetShapeFunctionsValues",&ConstitutiveLaw::Parameters::GetShapeFunctionsValues, return_internal_reference<>())
        .def("GetProcessInfo",&ConstitutiveLaw::Parameters::GetProcessInfo, return_internal_reference<>())
        .def("GetMaterialProperties",&ConstitutiveLaw::Parameters::GetMaterialProperties, return_internal_reference<>())
        .def("GetElementGeometry",&ConstitutiveLaw::Parameters::GetElementGeometry, return_internal_reference<>())    
    ;

    
    class_< ConstitutiveLaw, ConstitutiveLaw::Pointer , bases<Flags>, boost::noncopyable >
    ("ConstitutiveLaw", init<>() )
    .def("Clone",&ConstitutiveLaw::Clone)
    .def("WorkingSpaceDimension",&ConstitutiveLaw::WorkingSpaceDimension)
    .def("GetStrainSize",&ConstitutiveLaw::GetStrainSize)
    .def("GetStressMeasure",&ConstitutiveLaw::GetStressMeasure)
    .def("IsIncremental",&ConstitutiveLaw::IsIncremental)
    .def("Has", &ConstitutiveLawHas< Variable<int> >)
    .def("Has", &ConstitutiveLawHas< Variable<double> >)
    .def("Has", &ConstitutiveLawHas< Variable<array_1d<double,3> > >)
    .def("Has", &ConstitutiveLawHas< Variable<Vector> >)
    .def("Has", &ConstitutiveLawHas< Variable<Matrix> >)
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
    .def("CalculateMaterialResponse",&NewInterfaceCalculateMaterialResponse)
    .def("CalculateMaterialResponsePK1",&ConstitutiveLaw::CalculateMaterialResponsePK1)
    .def("CalculateMaterialResponsePK2",&ConstitutiveLaw::CalculateMaterialResponsePK2)
    .def("CalculateMaterialResponseKirchhoff",&ConstitutiveLaw::CalculateMaterialResponseKirchhoff)
    .def("CalculateMaterialResponseCauchy",&ConstitutiveLaw::CalculateMaterialResponseCauchy)
    .def("FinalizeMaterialResponse",&ConstitutiveLaw::FinalizeMaterialResponse)
    .def("FinalizeMaterialResponsePK1",&ConstitutiveLaw::FinalizeMaterialResponsePK1)
    .def("FinalizeMaterialResponsePK2",&ConstitutiveLaw::FinalizeMaterialResponsePK2)
    .def("FinalizeMaterialResponseKirchhoff",&ConstitutiveLaw::FinalizeMaterialResponseKirchhoff)
    .def("FinalizeMaterialResponseCauchy",&ConstitutiveLaw::FinalizeMaterialResponseCauchy)
    .def("ResetMaterial",&ConstitutiveLaw::ResetMaterial)
    .def("TransformStrains",&ConstitutiveLaw::TransformStrains, return_internal_reference<>())
//     .def("TransformStresses",&ConstitutiveLaw::TransformStresses)
//     .def("TransformStresses",&ConstitutiveLaw::TransformStresses)
    .def("TransformPK1Stresses",&ConstitutiveLaw::TransformPK1Stresses, return_internal_reference<>())
    .def("TransformPK2Stresses",&ConstitutiveLaw::TransformPK2Stresses, return_internal_reference<>())
    .def("TransformKirchhoffStresses",&ConstitutiveLaw::TransformKirchhoffStresses, return_internal_reference<>())
    .def("TransformCauchyStresses",&ConstitutiveLaw::TransformCauchyStresses, return_internal_reference<>())
    .def("PullBackConstitutiveMatrix",&ConstitutiveLaw::PullBackConstitutiveMatrix)
    .def("PushForwardConstitutiveMatrix",&ConstitutiveLaw::PushForwardConstitutiveMatrix)
//     .def("GetLawFeatures",&ConstitutiveLaw::GetLawFeatures)
    .def_readonly("COMPUTE_STRAIN", &ConstitutiveLaw::COMPUTE_STRAIN)
    .def_readonly("COMPUTE_STRESS", &ConstitutiveLaw::COMPUTE_STRESS)
    .def_readonly("COMPUTE_CONSTITUTIVE_TENSOR", &ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)
    .def_readonly("COMPUTE_STRAIN_ENERGY", &ConstitutiveLaw::COMPUTE_STRAIN_ENERGY)
    .def_readonly("ISOCHORIC_TENSOR_ONLY", &ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)
    .def_readonly("VOLUMETRIC_TENSOR_ONLY", &ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)
    .def_readonly("TOTAL_TENSOR", &ConstitutiveLaw::TOTAL_TENSOR)
    .def_readonly("FINALIZE_MATERIAL_RESPONSE", &ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE)
    .def_readonly("FINITE_STRAINS", &ConstitutiveLaw::FINITE_STRAINS)
    .def_readonly("INFINITESIMAL_STRAINS", &ConstitutiveLaw::INFINITESIMAL_STRAINS)
    .def_readonly("PLANE_STRAIN_LAW", &ConstitutiveLaw::PLANE_STRAIN_LAW)
    .def_readonly("PLANE_STRESS_LAW", &ConstitutiveLaw::PLANE_STRESS_LAW)
    .def_readonly("AXISYMMETRIC_LAW", &ConstitutiveLaw::AXISYMMETRIC_LAW)
    .def_readonly("U_P_LAW", &ConstitutiveLaw::U_P_LAW)
    .def_readonly("ISOTROPIC", &ConstitutiveLaw::ISOTROPIC)
    .def_readonly("ANISOTROPIC", &ConstitutiveLaw::ANISOTROPIC)
    ;

}
}  // namespace Python.
}  // namespace Kratos.
