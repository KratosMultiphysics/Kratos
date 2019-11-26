//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "custom_conditions/mapping_condition.h"
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

void MAPcond::Initialize(const ProcessInfo& rCurrentProcessInfo)
{

}

//***********************************************************************************
//***********************************************************************************

void MAPcond::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& r_process_info)
{

}


 void MAPcond::AddExplicitContribution(const VectorType& rRHS,
                         const Variable<VectorType>& rRHSVariable,
                         Variable<array_1d<double,3> >& rDestinationVariable,
                         const ProcessInfo& r_process_info)
{
    KRATOS_TRY


    KRATOS_CATCH( "" )
}

} // Namespace Kratos.
