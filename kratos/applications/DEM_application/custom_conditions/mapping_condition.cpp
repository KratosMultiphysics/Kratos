//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "includes/define.h"
#include "custom_conditions/mapping_condition.h"
#include "DEM_application.h"

#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{
	using namespace GeometryFunctions;

//***********************************************************************************
//***********************************************************************************


// Constructor

MAPcond::MAPcond()
{
}

// Constructor

MAPcond::MAPcond(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

// Constructor

MAPcond::MAPcond(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
    //setting up the nodal degrees of freedom
}

//***********************************************************************************
//***********************************************************************************

Condition::Pointer MAPcond::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new MAPcond(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//***********************************************************************************
//***********************************************************************************
// Destructor

MAPcond::~MAPcond()
{
}

//***********************************************************************************
//***********************************************************************************

void MAPcond::Initialize()
{
  
}

//***********************************************************************************
//***********************************************************************************

void MAPcond::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    
 }

  
 void MAPcond::AddExplicitContribution(const VectorType& rRHS,
                         const Variable<VectorType>& rRHSVariable,
                         Variable<array_1d<double,3> >& rDestinationVariable,
                         const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    KRATOS_CATCH( "" )
}

 
 

void MAPcond::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
{
    
}

void MAPcond::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)   
{  
  
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
