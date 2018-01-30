//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2017-01-09 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

// External includes

// Project includes
#include "rigid_body_element.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/variables.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    RigidBodyElement::RigidBodyElement() : Element() {}
            
    RigidBodyElement::RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {
        mpIntegrationScheme = NULL;
    }
      
    RigidBodyElement::RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {
        mpIntegrationScheme = NULL;
    }
      
    RigidBodyElement::RigidBodyElement(IndexType NewId, NodesArrayType const& ThisNodes)
    : Element(NewId, ThisNodes) {
        mpIntegrationScheme = NULL;
    }
    
    Element::Pointer RigidBodyElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
        return Element::Pointer(new RigidBodyElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }      
    
    // Destructor
    RigidBodyElement::~RigidBodyElement() {

    }
      
    void RigidBodyElement::Initialize(ProcessInfo& r_process_info) {

    }   
    
    void RigidBodyElement::SetIntegrationScheme(DEMIntegrationScheme::Pointer& integration_scheme){

    }
    
    void RigidBodyElement::CustomInitialize(ProcessInfo& r_process_info) {

    }    
    
    void RigidBodyElement::SetOrientation(const Quaternion<double> Orientation) {

    }

    void RigidBodyElement::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info)
    {

    } //Calculate

    void RigidBodyElement::UpdateLinearDisplacementAndVelocityOfNodes() {
               
    }    
    
    void RigidBodyElement::UpdatePositionOfNodes() {
                               
    }   
    
    void RigidBodyElement::CollectForcesAndTorquesFromNodes() {
           
    }   
    
    void RigidBodyElement::GetRigidBodyElementForce(const array_1d<double,3>& gravity) {
        
    }
    
    void RigidBodyElement::ComputeAdditionalForces(const array_1d<double,3>& gravity) {

    }   

    double RigidBodyElement::SlowGetDensity() { return GetProperties()[PARTICLE_DENSITY];}
    
    int RigidBodyElement::SlowGetRigidBodyElementMaterial() { return GetProperties()[PARTICLE_MATERIAL];}

    void RigidBodyElement::Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag ) {

    }   
}  // namespace Kratos.
