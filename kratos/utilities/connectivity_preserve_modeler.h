//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
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
    ) override
    {
        KRATOS_TRY;

        //TODO: do this better, with the remove function
//         DestinationModelPart.Nodes().clear();
//         DestinationModelPart.Conditions().clear();
//         DestinationModelPart.Elements().clear();
        
        for(auto it = DestinationModelPart.NodesBegin(); it != DestinationModelPart.NodesEnd(); it++)
            it->Set(TO_ERASE);
        DestinationModelPart.RemoveNodes(TO_ERASE);

        for(auto it = DestinationModelPart.ElementsBegin(); it != DestinationModelPart.ElementsEnd(); it++)
            it->Set(TO_ERASE);
        DestinationModelPart.RemoveElements(TO_ERASE);
        
        for(auto it = DestinationModelPart.ElementsBegin(); it != DestinationModelPart.ElementsEnd(); it++)
            it->Set(TO_ERASE);
        DestinationModelPart.RemoveConditions(TO_ERASE);


        //assigning ProcessInfo
        DestinationModelPart.SetProcessInfo(  OriginModelPart.pGetProcessInfo() );
        
        if(DestinationModelPart.IsSubModelPart())
        {
            if( DestinationModelPart.GetBufferSize()!=OriginModelPart.GetBufferSize())
            {
                KRATOS_THROW_ERROR(std::runtime_error,"DestinationModelPart  is a SubModelPart and its buffer size does not coincide with the one of the original model. Setting of the buffer size is not possible","")
            }
        }
        else
            DestinationModelPart.SetBufferSize(OriginModelPart.GetBufferSize());

        //assigning Properties
        DestinationModelPart.SetProperties(OriginModelPart.pProperties());

        //assigning the nodes to the new model part
        DestinationModelPart.AddNodes(OriginModelPart.NodesBegin(), OriginModelPart.NodesEnd()); // = OriginModelPart.Nodes();

        //generating the elements
        ModelPart::ElementsContainerType temp_elements;
        for (ModelPart::ElementsContainerType::iterator iii = OriginModelPart.ElementsBegin(); iii != OriginModelPart.ElementsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();
            Element::Pointer p_element = rReferenceElement.Create(iii->Id(), iii->GetGeometry(), properties);
            
            //actually use the geometry of the old element, so that memory is saved!!
            p_element->pGetGeometry() = iii->pGetGeometry();
            
            temp_elements.push_back(p_element);
        }
        DestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());

        //generating the conditions
        ModelPart::ConditionsContainerType temp_conditions;
        for (ModelPart::ConditionsContainerType::iterator iii = OriginModelPart.ConditionsBegin(); iii != OriginModelPart.ConditionsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();

            Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(iii->Id(), iii->GetGeometry(), properties);
            
            //assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_condition->pGetGeometry() = iii->pGetGeometry();
            
            temp_conditions.push_back(p_condition);
        }
        DestinationModelPart.AddConditions(temp_conditions.begin(), temp_conditions.end());
        
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


