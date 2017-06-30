//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala 
//
//


#ifndef DUPLICATE_INTERFACE_PROCESS_H
#define DUPLICATE_INTERFACE_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "processes/find_nodal_neighbours_process.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{

class DuplicateInterfaceProcess : public Process
{
public:

    /// Pointer definition of DuplicateInterfaceProcess
    KRATOS_CLASS_POINTER_DEFINITION(DuplicateInterfaceProcess);

    /// Constructor.
    DuplicateInterfaceProcess(  ModelPart& main_model_part ) : mr_main_model_part(main_model_part)
    {   
        IsNeighbourFormulated = false;
        CalculateNeighbourInformation();
    }

    /// Destructor.
    virtual ~DuplicateInterfaceProcess(){

    }


    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_TRY;



        KRATOS_CATCH("");
    }


    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        
        KRATOS_CATCH("");
    }

    void DuplicateInterface(ModelPart& model_part_1, ModelPart& model_part_2, ModelPart& interface_model_part){
        long int currentNumNodes = mr_main_model_part.NumberOfNodes();
        std::cout<<"interface_model_part name :: "<<interface_model_part.Name()<<std::endl;
        if(!IsNeighbourFormulated){
            CalculateNeighbourInformation();
        } else {

            for(ModelPart::NodesContainerType::iterator interface_node = interface_model_part.NodesBegin();
                                                        interface_node != interface_model_part.NodesEnd(); interface_node++){
                // Create a duplicat of the current interface node
                unsigned int oldNodeId = interface_node->Id();                
                mr_main_model_part.CreateNewNode(currentNumNodes+1, interface_node->X(), interface_node->Y(), interface_node->Z());
                currentNumNodes++;
                unsigned int newNodeId = currentNumNodes;
                
                // Access the element neighbours of the current interface node
                WeakPointerVector< Element >& rneigh_el = interface_node->GetValue(NEIGHBOUR_ELEMENTS);
                for(auto nei_elem : rneigh_el){
                    const unsigned int number_of_nodes = nei_elem.GetGeometry().PointsNumber();
                    if(model_part_1.GetMesh().HasElement(nei_elem.Id())){ // If the model_part_1 has the element replace this element's node
                        for (unsigned int j = 0; j < number_of_nodes; j++)
                        {
                            if (nei_elem.GetGeometry()[j].Id() == oldNodeId)
                            { 
                                nei_elem.GetGeometry()[j] = mr_main_model_part.Nodes()[newNodeId];
                            }
                        }
                    }
                }
            }
        }
    }

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DuplicateInterfaceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "DuplicateInterfaceProcess";}

    /// Print object's data.
    void PrintData() {
        
    }



protected:
    ///@name Protected static Member Variables
    ///@{
        ModelPart&                                 mr_main_model_part;
    ///@}
    ///@name Protected member Variables
    ///@{
        bool IsNeighbourFormulated;
private:

    /// Assignment operator.
    DuplicateInterfaceProcess& operator=(DuplicateInterfaceProcess const& rOther){return *this;}

    void CalculateNeighbourInformation(){
       FindNodalNeighboursProcess neighbour_finder = FindNodalNeighboursProcess(mr_main_model_part, 10, 25); // 10 and 25 just approximate values
       neighbour_finder.Execute();
       IsNeighbourFormulated = true;
    }

}; // Class MoveRotorProcess

};  // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
