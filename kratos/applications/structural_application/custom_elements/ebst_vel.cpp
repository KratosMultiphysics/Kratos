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
#include "custom_elements/ebst_vel.h"
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

    EbstVel::EbstVel()
    {
    }

    // Constructor

    EbstVel::EbstVel(IndexType NewId, GeometryType::Pointer pGeometry)
    : Ebst(NewId, pGeometry)
    {
    }

    // Constructor

    EbstVel::EbstVel(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Ebst(NewId, pGeometry, pProperties)
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    Element::Pointer EbstVel::Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties) const
    {

        return Element::Pointer(new EbstVel(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    //***********************************************************************************
    //***********************************************************************************
    // Destructor

    EbstVel::~EbstVel()
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    void EbstVel::EquationIdVector(
            EquationIdVectorType& rResult,
            ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);

        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);
        unsigned int dim = number_of_nodes * 3;

        if (rResult.size() != dim)
            rResult.resize(dim, false);

        //nodes of the central element
        for (int i = 0; i < 3; i++)
        {
            int index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
        }

        //adding the ids ofthe neighbouring nodes
        int index = 9;
        for (unsigned int i = 0; i < 3; i++)
        {
            if (HasNeighbour(i, neigb[i]))
            {
                rResult[index] = neigb[i].GetDof(VELOCITY_X).EquationId();
                rResult[index + 1] = neigb[i].GetDof(VELOCITY_Y).EquationId();
                rResult[index + 2] = neigb[i].GetDof(VELOCITY_Z).EquationId();
                index += 3;
            }
        }

        
        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void EbstVel::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);
        ElementalDofList.resize(0);

        //nodes of the central element
        for (unsigned int i = 0; i < GetGeometry().size(); i++)
        {
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
        }

        //adding the dofs ofthe neighbouring nodes
        for (unsigned int i = 0; i < 3; i++)
        {
            if (HasNeighbour(i, neigb[i]))
            {
                ElementalDofList.push_back(neigb[i].pGetDof(VELOCITY_X));
                ElementalDofList.push_back(neigb[i].pGetDof(VELOCITY_Y));
                ElementalDofList.push_back(neigb[i].pGetDof(VELOCITY_Z));
            }
        }
        KRATOS_CATCH("")
    }



} // Namespace Kratos.
