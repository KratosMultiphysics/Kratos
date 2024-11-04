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
#include <iostream>

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
    auto P = this->GetIntegrationMethod();
    auto N = row(this->GetGeometry().ShapeFunctionsValues(this->GetIntegrationMethod()), PointNumber);

    double interpolated_coupling_weight_at_int_point=0; 
    //std::cout<<"shape function at int point: "<<N;
    //std::cout<<"geometry size: "<<this->GetGeometry().size();
    // if condition to check if the integration pointy cordinate is less than 0.16
    // if (this->GetGeometry().Center()[1] < 0.16 ) // this line is important to get expected max stress 

    // {

   
    for (int i=0; i < this->GetGeometry().size(); i++)
    {
         
        //   if (this->GetGeometry().Center()[1] < 0.16) // this line is important to get expected max stress  // check y cordinate of integration point to be less than 0.16
                //  {
                   interpolated_coupling_weight_at_int_point += N[i]*this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT);
                //  }  // interpolated_coupling_weight_at_int_point += this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT)/8; // same weight on all integration points in the element  
       //std::cout<<" nodal coupling weight ="<<this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT)<<std::endl;
       //KRATOS_WATCH( N[i]);
       
    }

    // }

    // interpolated_coupling_weight_at_int_point = 0.5;
    return (1-interpolated_coupling_weight_at_int_point) * SmallDisplacement::GetIntegrationWeight(rThisIntegrationPoints,PointNumber,detJ);
   

}



void VolumeCouplingElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
) {
    // Call base class function to fill rOutput with stress values
    SmallDisplacement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    // print element number
    // std::cout<<"element number: "<<this->Id()<<std::endl;

    // Check if we are dealing with stress vectors and if the custom condition is met
    if ((rVariable == CAUCHY_STRESS_VECTOR) && (rCurrentProcessInfo[ACTIVATION_LEVEL] > 0.0)) {
        // Access integration points
        //print activation level with std::cout
        // std::cout<<"activation level: "<<rCurrentProcessInfo[ACTIVATION_LEVEL]<<std::endl;
        // std::cout<<"For stress:"<<std::endl;
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints(this->GetIntegrationMethod());
        const SizeType number_of_integration_points = integration_points.size();
        array_1d<double, 3> global_coordinates;
        

        // if (this->GetGeometry().Center()[1] < 0.16) // this line is important to get expected max stress 
        //     {
            //std::cout<<"##############################################Element centre location y location : "<<this->GetGeometry().Center()[1]<<std::endl;
            // Loop over each integration point
            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {

            
            // Calculate the coupling weight for this integration point
            auto N = row(this->GetGeometry().ShapeFunctionsValues(this->GetIntegrationMethod()), point_number);
            double coupling_weight_on_this_integration_point = 0.0;
            // check y cordinate of integration point to be less than 0.16

            for (SizeType i = 0; i < this->GetGeometry().size(); ++i) 
            {
                //  if (this->GetGeometry().Center()[1] < 0.16) // this line is important to get expected max stress  // check y cordinate of integration point to be less than 0.16
                //  {
                    coupling_weight_on_this_integration_point += N[i] * this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT); 
                //  }  //    coupling_weight_on_this_integration_point += this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT)/8;  // same weight on all integration points in the element  
            }
            

            // for experiements:
            // coupling_weight_on_this_integration_point =0.5;
            // Scale the stress vector at this integration point by the calculated coupling weight
            
            rOutput[point_number] *= (1-coupling_weight_on_this_integration_point);

            std::cout<<"integration point weight "<<(1-coupling_weight_on_this_integration_point)<<std::endl;
            // print the location of the integration point
            std::cout<<"integration point location: "<<integration_points[point_number].Coordinates()<<std::endl;
            // print the location of the element
            std::cout<<"element location: "<<this->GetGeometry().Center()<<std::endl;
            // print the 'CAUCHY_STRESS_VECTOR' for this element
            std::cout<<"CAUCHY_STRESS_VECTOR: "<<rOutput[point_number]<<std::endl;
            // get cordinates of integration point using globalCoordinates function
            //auto global_coordinates = ZeroVector(3);
            
            global_coordinates = this->GetGeometry().GlobalCoordinates(global_coordinates,integration_points[point_number].Coordinates());
            // print the global coordinates of the integration point
            std::cout<<"global coordinates: "<<global_coordinates<<std::endl;

  
            }
            // }
            // else
            // {       
                
            //     for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number)
            //      {      
            //             double coupling_weight_on_this_integration_point = 0.0;
            //             std::cout<<"integration point weight "<<(1-coupling_weight_on_this_integration_point)<<std::endl;
            //             // print the location of the integration point
            //             std::cout<<"integration point location: "<<integration_points[point_number].Coordinates()<<std::endl;
            //             // print the location of the element
            //             std::cout<<"element location: "<<this->GetGeometry().Center()<<std::endl;
            //             // print the 'CAUCHY_STRESS_VECTOR' for this element
            //             std::cout<<"CAUCHY_STRESS_VECTOR: "<<rOutput[point_number]<<std::endl;
            //             // get cordinates of integration point using globalCoordinates function
            //             //auto global_coordinates = ZeroVector(3);
                        
            //             global_coordinates = this->GetGeometry().GlobalCoordinates(global_coordinates,integration_points[point_number].Coordinates());
            //             // print the global coordinates of the integration point
            //             std::cout<<"global coordinates: "<<global_coordinates<<std::endl;
            //     }
            // }
        
    }
    if ((rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) && (rCurrentProcessInfo[ACTIVATION_LEVEL] > 0.0)) {
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints(this->GetIntegrationMethod());
        const SizeType number_of_integration_points = integration_points.size();
        array_1d<double, 3> global_coordinates;
        
        // std::cout<<"For strain:"<<std::endl;
        // if (this->GetGeometry().Center()[1] < 0.16) // this line is important to get expected max stress 
        //     {
        //     //std::cout<<"##############################################Element centre location y location : "<<this->GetGeometry().Center()[1]<<std::endl;
        //     // Loop over each integration point
            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {

            
            // Calculate the coupling weight for this integration point
            auto N = row(this->GetGeometry().ShapeFunctionsValues(this->GetIntegrationMethod()), point_number);
            double coupling_weight_on_this_integration_point = 0.0;
            

            for (SizeType i = 0; i < this->GetGeometry().size(); ++i) 
            {       
                //  if (this->GetGeometry().Center()[1] < 0.16) // this line is important to get expected max stress  // check y cordinate of integration point to be less than 0.16
                //  {
                    coupling_weight_on_this_integration_point += N[i] * this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT); 
                //  } 
                //    coupling_weight_on_this_integration_point += this->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT)/8;  // same weight on all integration points in the element  
            }
            


            
            std::cout<<"integration point weight "<<(1-coupling_weight_on_this_integration_point)<<std::endl;
            // print the location of the integration point
            std::cout<<"integration point location: "<<integration_points[point_number].Coordinates()<<std::endl;
            // print the location of the element
            std::cout<<"element location: "<<this->GetGeometry().Center()<<std::endl;
            // print the 'CAUCHY_STRAIN_VECTOR' for this element
            std::cout<<"GREEN_LAGRANGE_STRAIN_VECTOR: "<<rOutput[point_number]<<std::endl;
            // get cordinates of integration point using globalCoordinates function
            //auto global_coordinates = ZeroVector(3);
            
            global_coordinates = this->GetGeometry().GlobalCoordinates(global_coordinates,integration_points[point_number].Coordinates());
            // print the global coordinates of the integration point
            std::cout<<"global coordinates: "<<global_coordinates<<std::endl;

  
            }
         }
    //         else
    //         {       
                
    //             for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number)
    //              {      
    //                     double coupling_weight_on_this_integration_point = 0.0;
    //                     std::cout<<"integration point weight "<<(1-coupling_weight_on_this_integration_point)<<std::endl;
    //                     // print the location of the integration point
    //                     std::cout<<"integration point location: "<<integration_points[point_number].Coordinates()<<std::endl;
    //                     // print the location of the element
    //                     std::cout<<"element location: "<<this->GetGeometry().Center()<<std::endl;
    //                     // print the 'CAUCHY_STRESS_VECTOR' for this element
    //                     std::cout<<"GREEN_LAGRANGE_STRAIN_VECTOR: "<<rOutput[point_number]<<std::endl;
    //                     // get cordinates of integration point using globalCoordinates function
    //                     //auto global_coordinates = ZeroVector(3);
                        
    //                     global_coordinates = this->GetGeometry().GlobalCoordinates(global_coordinates,integration_points[point_number].Coordinates());
    //                     // print the global coordinates of the integration point
    //                     std::cout<<"global coordinates: "<<global_coordinates<<std::endl;
    //             }
    //         }
        
    // }
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