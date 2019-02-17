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

/* **************************************************************************************
 *
 *   Last Modified by:    $Author: rrossi $
 *   Date:                $Date: 2008-10-13 07:00:53 $
 *   Revision:            $Revision: 1.12 $
 *
 * ***************************************************************************************/


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/total_lagrangian_vel.h"
#include "includes/constitutive_law.h"
#include "structural_application.h"

namespace Kratos
{



//***********************************************************************************
//***********************************************************************************
// -------- //
//  PUBLIC  //
// -------- //

// Constructor

TotalLagrangianVelocityBased::TotalLagrangianVelocityBased(IndexType NewId, GeometryType::Pointer pGeometry)
    : TotalLagrangian(NewId, pGeometry)
{
}

// Constructor

TotalLagrangianVelocityBased::TotalLagrangianVelocityBased(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : TotalLagrangian(NewId, pGeometry, pProperties)
{
}

//***********************************************************************************
//***********************************************************************************

Element::Pointer TotalLagrangianVelocityBased::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{

    return Element::Pointer(new TotalLagrangianVelocityBased(NewId, GetGeometry().Create(ThisNodes), pProperties));

}

//***********************************************************************************
//***********************************************************************************
// Destructor

TotalLagrangianVelocityBased::~TotalLagrangianVelocityBased()
{
}

//***********************************************************************************
//***********************************************************************************

void TotalLagrangianVelocityBased::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        int number_of_nodes = GetGeometry().size();
        int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int dim2 = number_of_nodes * dim;

        if ( rResult.size() != dim2 )
            rResult.resize( dim2, false );

        for ( int i = 0; i < number_of_nodes; i++ )
        {
            int index = i * dim;
            rResult[index] = GetGeometry()[i].GetDof( VELOCITY_X ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof( VELOCITY_Y ).EquationId();

            if ( dim == 3 )
                rResult[index + 2] = GetGeometry()[i].GetDof( VELOCITY_Z ).EquationId();
        }
		
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void TotalLagrangianVelocityBased::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Y ) );

            if ( GetGeometry().WorkingSpaceDimension() == 3 )
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Z ) );
            }
        }
		
    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************



} // Namespace Kratos.
