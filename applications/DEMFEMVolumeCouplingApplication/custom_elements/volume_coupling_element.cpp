// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    
//

// System includes

// External includes

// Project includes
// #include "utilities/math_utils.h"
// #include "utilities/geometry_utilities.h"
// #include "utilities/atomic_utilities.h"

// Application includes
// #include "custom_elements/base_solid_element.h"
// #include "structural_mechanics_application_variables.h"
// #include "custom_utilities/structural_mechanics_element_utilities.h"
// #include "custom_utilities/constitutive_law_utilities.h"
#include "../DEMFEM_volume_coupling_application.h"
#include "volume_coupling_element.h"


// 
namespace Kratos
{

VolumeCouplingElement::VolumeCouplingElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : SmallDisplacement(NewId, pGeometry)
{
    // Custom constructor code here
}

VolumeCouplingElement::VolumeCouplingElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SmallDisplacement(NewId, pGeometry, pProperties)
{
    // Custom constructor code here
}
VolumeCouplingElement::VolumeCouplingElement()
{
    // Custom implementation code here
}

Element::Pointer VolumeCouplingElement::Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<VolumeCouplingElement>(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

Element::Pointer VolumeCouplingElement::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<VolumeCouplingElement>(NewId, pGeom, pProperties);
}



Element::Pointer VolumeCouplingElement::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
{
    VolumeCouplingElement::Pointer p_new_elem = Kratos::make_intrusive<VolumeCouplingElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Dynamic casting to the correct type
    BaseSolidElement* p_base_solid_elem = dynamic_cast<BaseSolidElement*>(p_new_elem.get());
    if (p_base_solid_elem != nullptr) {
        p_base_solid_elem->SetIntegrationMethod(this->GetIntegrationMethod());
        p_base_solid_elem->SetConstitutiveLawVector(this->GetConstitutiveLawVector());
    }

    return p_new_elem;
}


VolumeCouplingElement::~VolumeCouplingElement()
{
    // Custom destructor code here
}

double VolumeCouplingElement::GetIntegrationWeight(const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
                                                   const IndexType PointNumber,
                                                   const double detJ) const 
{
    //auto P = this->GetIntegrationMethod();
    // auto N = row(this->GetGeometry().ShapeFunctionsValues(this->GetIntegrationMethod()), PointNumber);
    // double interpolated_coupling_weight_at_int_point=0; 
    // //std::cout<<"shape function at int point: "<<N;
    // //std::cout<<"geometry size: "<<this->GetGeometry().size();
    // for (int i=0; i < this->GetGeometry().size(); i++)
    // {
    //    interpolated_coupling_weight_at_int_point += N[i]*this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT);
    //    //KRATOS_WATCH(this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT));
    //    //KRATOS_WATCH( N[i]);
       
    // }
    //KRATOS_WATCH(interpolated_coupling_weight_at_int_point);
    //function to check if interpolated_coupling_weight_at_int_point is not zero
    // if (interpolated_coupling_weight_at_int_point!=0)
    // {
    //     return (0.5) * SmallDisplacement::GetIntegrationWeight(rThisIntegrationPoints,PointNumber,detJ);
    // }
    // else
    // {
    //     return SmallDisplacement::GetIntegrationWeight(rThisIntegrationPoints,PointNumber,detJ);
    // }
    return SmallDisplacement::GetIntegrationWeight(rThisIntegrationPoints,PointNumber,detJ);
    // interpolated_coupling_weight_at_int_point = 0; // incase weidght dem + weight fem !=1, this must be used for hybrid region to have 0 dem-weight => 1 fem weight.
    // return (1-interpolated_coupling_weight_at_int_point) * SmallDisplacement::GetIntegrationWeight(rThisIntegrationPoints,PointNumber,detJ);

}

} // Namespace Kratos


// namespace Kratos
// {

// VolumeCouplingElement::VolumeCouplingElement( IndexType NewId, GeometryType::Pointer pGeometry )
//     : SmallDisplacement( NewId, pGeometry )
// {
    
// }

// VolumeCouplingElement::VolumeCouplingElement( ) // Default constructor needed for serialization
//     : SmallDisplacement( )
// {
    
// }

// double VolumeCouplingElement::GetIntegrationWeight(
//     const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
//     const IndexType PointNumber,
//     const double detJ
//     ) const 
// {

//     auto N = row(this->GetGeometry().ShapeFunctionsValues(this->GetIntegrationMethod()), PointNumber);
//     double interpolated_coupling_weight_at_int_point=0; 
//     for (int i=0; i < this->GetGeometry().size(); i++)
//     {
//        interpolated_coupling_weight_at_int_point += N[i]*this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT);
//     }
//     return (1-interpolated_coupling_weight_at_int_point) * SmallDisplacement::GetIntegrationWeight(rThisIntegrationPoints,PointNumber,detJ);
//     KRATOS_WATCH(interpolated_coupling_weight_at_int_point);
// }


// } // Namespace Kratos


// Element::Pointer VolumeCouplingElement::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
// {
//     VolumeCouplingElement::Pointer p_new_elem = Kratos::make_intrusive<VolumeCouplingElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
//     p_new_elem->SetData(this->GetData());
//     p_new_elem->Set(Flags(*this));

//     // Setting integration method and constitutive law vector using the accessor methods
//     p_new_elem->mThisIntegrationMethod=mThisIntegrationMethod;
//     p_new_elem->mConstitutiveLawVector=mConstitutiveLawVector;

//     return p_new_elem;
// }
// Element::Pointer VolumeCouplingElement::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
// {
//     VolumeCouplingElement::Pointer p_new_elem = Kratos::make_intrusive<VolumeCouplingElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
//     p_new_elem->SetData(this->GetData());
//     p_new_elem->Set(Flags(*this));

//     // Dynamic casting to the correct type
//     BaseSolidElement* p_base_solid_elem = dynamic_cast<BaseSolidElement*>(p_new_elem.get());
//     if (p_base_solid_elem != nullptr) {
//         p_base_solid_elem->SetIntegrationMethod(mThisIntegrationMethod);
//         p_base_solid_elem->SetConstitutiveLawVector(mConstitutiveLawVector);
//     }

//     return p_new_elem;
// }