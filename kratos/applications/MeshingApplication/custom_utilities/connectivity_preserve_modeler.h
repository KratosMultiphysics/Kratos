// Kratos Multi-Physics
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 	
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-11-19 15:38:01 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_GENERATE_MODEL_PART_MODELER_INCLUDED )
#define  KRATOS_GENERATE_MODEL_PART_MODELER_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "modeler/modeler.h"

#include "meshing_application.h"
namespace Kratos
{

class ConnectivityPreserveModeler : public Modeler
{
public:

    ConnectivityPreserveModeler()
    {
    }

    /// Destructor.

    virtual ~ConnectivityPreserveModeler()
    {
    }

    //**********************************************************************************************
    //**********************************************************************************************
    ///This function fills the @param DestinationModelPart using the data obtained from @param  OriginModelPart
    ///the elements and conditions of the DestinationModelPart part use the same connectivity (and id) as the
    ///OriginModelPart but their type is determined by @param rReferenceElement and @param rReferenceBoundaryCondition

    virtual void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition
    )
    {
        KRATOS_TRY;

        DestinationModelPart.Nodes().clear();


        DestinationModelPart.Conditions().clear();
        DestinationModelPart.Elements().clear();




        //assigning ProcessInfo
        DestinationModelPart.SetProcessInfo(  OriginModelPart.pGetProcessInfo() );
        DestinationModelPart.SetBufferSize(OriginModelPart.GetBufferSize());
//             DestinationModelPart.SetProcessInfo( OriginModelPart.GetProcessInfo() );

        //assigning Properties
        DestinationModelPart.SetProperties(OriginModelPart.pProperties());

        //assigning the nodes to the new model part

        DestinationModelPart.Nodes() = OriginModelPart.Nodes();

        //generating the elements
        for (ModelPart::ElementsContainerType::iterator iii = OriginModelPart.ElementsBegin(); iii != OriginModelPart.ElementsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();
            Element::Pointer p_element = rReferenceElement.Create(iii->Id(), iii->GetGeometry(), properties);
            
            //actually use the geometry of the old element, so that memory is saved!!
            p_element->pGetGeometry() = iii->pGetGeometry();
            
            DestinationModelPart.Elements().push_back(p_element);
        }

        //generating the conditions
        for (ModelPart::ConditionsContainerType::iterator iii = OriginModelPart.ConditionsBegin(); iii != OriginModelPart.ConditionsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();

            Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(iii->Id(), iii->GetGeometry(), properties);
            
            //assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_condition->pGetGeometry() = iii->pGetGeometry();
            
            DestinationModelPart.Conditions().push_back(p_condition);
        }
        
        //generating tables
	DestinationModelPart.Tables() = OriginModelPart.Tables();


        Communicator::Pointer pComm = OriginModelPart.GetCommunicator().Create();
        DestinationModelPart.SetCommunicator(pComm);

        KRATOS_CATCH("");
    }


    //**********************************************************************************************
    //**********************************************************************************************
private:

    ConnectivityPreserveModeler & operator=(ConnectivityPreserveModeler const& rOther);

    /// Copy constructor.
    ConnectivityPreserveModeler(ConnectivityPreserveModeler const& rOther);


};

} // namespace Kratos.

#endif //KRATOS_GENERATE_MODEL_PART_MODELER_INCLUDED  defined


