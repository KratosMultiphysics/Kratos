/*
==============================================================================
KratosStructuralApplication
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
//   Last modified by:    $Author: virginia $
//   Date:                $Date: 2009-01-23 14:39:59 $
//   Revision:            $Revision: 1.27 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"
#include "DEM_application.h"

//#include <omp.h>

namespace Kratos
{
//************************************************************************************
//************************************************************************************
Particle_Contact_Element::Particle_Contact_Element( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); //fa falta???
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
Particle_Contact_Element::Particle_Contact_Element( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();//fa falta???
}


//create contact elements instances.

Element::Pointer Particle_Contact_Element::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new Particle_Contact_Element( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Particle_Contact_Element::~Particle_Contact_Element()
{
}

std::string Particle_Contact_Element::Info() const
{
    std::stringstream buffer;
    buffer << "Particle Contact Element" << std::endl;
    return buffer.str();
}



void Particle_Contact_Element::Initialize() {
    KRATOS_TRY
    
    mMeanContactArea = 0.0;    
    mFailureCriterionState = 0.0;    
    mLocalContactAreaLow  = 0.0;
    mLocalContactAreaHigh = 0.0;
    mLocalContactForce[0] = 0.0;
    mLocalContactForce[1] = 0.0;
    mLocalContactForce[2] = 0.0;
    mUnidimendionalDamage = 0.0;
    mContactFailure = 0.0;
    mContactSigma = 0.0;
    mContactTau = 0.0;
    
    
    this->GetValue(LOW_POISSON_FORCE) = 0.0;  
    this->GetValue(HIGH_POISSON_FORCE) = 0.0; 
    this->GetValue(LOCAL_CONTACT_FORCE)[0] = 0.0;
    this->GetValue(LOCAL_CONTACT_FORCE)[1] = 0.0;
    this->GetValue(LOCAL_CONTACT_FORCE)[2] = 0.0;
    this->GetValue(CONTACT_SIGMA)          = 0.0;
    this->GetValue(CONTACT_TAU)            = 0.0;
    this->GetValue(CONTACT_FAILURE)        = 0.0;
    this->GetValue(FAILURE_CRITERION_STATE) = 0.0;
    this->GetValue(UNIDIMENSIONAL_DAMAGE)  = 0.0; 
   
    KRATOS_CATCH( "" )
}

void Particle_Contact_Element::PrepareForPrinting() {
    KRATOS_TRY
            
    this->GetValue(LOCAL_CONTACT_FORCE)[0] = mLocalContactForce[0];
    this->GetValue(LOCAL_CONTACT_FORCE)[1] = mLocalContactForce[1];
    this->GetValue(LOCAL_CONTACT_FORCE)[2] = mLocalContactForce[2];
    this->GetValue(CONTACT_SIGMA)          = mContactSigma;    
    this->GetValue(CONTACT_TAU)            = mContactTau;
    this->GetValue(CONTACT_FAILURE)        = mContactFailure;
    this->GetValue(FAILURE_CRITERION_STATE)= mFailureCriterionState;
    this->GetValue(UNIDIMENSIONAL_DAMAGE)  = mUnidimendionalDamage; 
    
    KRATOS_CATCH( "" )    
}

void Particle_Contact_Element::CalculateMeanContactArea(const bool has_mpi)
{
    KRATOS_TRY
            
    if(!has_mpi) mMeanContactArea = 0.5 * ( mLocalContactAreaLow + mLocalContactAreaHigh );
    else mMeanContactArea = mLocalContactAreaLow;
    
    KRATOS_CATCH( "" )
    
}


void Particle_Contact_Element::GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
        
    //if(rVariable == LOCAL_CONTACT_FORCE) {  //3D VARIABLE WITH COMPONENTS                   
        rOutput.resize(1);
        const Particle_Contact_Element* const_this = dynamic_cast< const Particle_Contact_Element* >(this); //To ensure we don't set the value here
        rOutput[0][0] = const_this->GetValue(rVariable)[0];
        rOutput[0][1] = const_this->GetValue(rVariable)[1];
        rOutput[0][2] = const_this->GetValue(rVariable)[2];        
    //}    
}
 
  void Particle_Contact_Element::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo) {        
        Output.resize(1);
        const Particle_Contact_Element* const_this = dynamic_cast< const Particle_Contact_Element* >(this); //To ensure we don't set the value here                
        Output[0] = double(const_this->GetValue(rVariable)); 
  } 
    


//************************************************************************************
//************************************************************************************
void Particle_Contact_Element::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
  
}

////************************************************************************************
////************************************************************************************

void Particle_Contact_Element::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    mContactTau           = 0.0;
    mContactSigma         = 0.0;    
    mLocalContactForce[0] = 0.0;
    mLocalContactForce[1] = 0.0;
    mLocalContactForce[2] = 0.0;
    if (mFailureCriterionState<1.0) {
        mFailureCriterionState = 0.0;
    } // else we keep it at 1.0.

}

////************************************************************************************
////************************************************************************************
void Particle_Contact_Element::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{

 
}

void Particle_Contact_Element::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
{ 
 
  if (rVariable == DUMMY_DEBUG_DOUBLE)
  {
    
    /*
    if(this->GetValue(CONTACT_SIGMA) != 0.0)
    {
     //KRATOS_WATCH( ( this->GetValue(CONTACT_SIGMA) - this->GetValue(CONTACT_TAU) ) /(this->GetValue(CONTACT_SIGMA) ) )
    }
    */
   // KRATOS_WATCH( (this->GetValue(LOCAL_CONTACT_FORCE_HIGH)[2]-this->GetValue(LOCAL_CONTACT_FORCE_LOW)[2])/this->GetValue(LOCAL_CONTACT_FORCE_HIGH)[2]  )
    
  }
  
}










} // Namespace Kratos


