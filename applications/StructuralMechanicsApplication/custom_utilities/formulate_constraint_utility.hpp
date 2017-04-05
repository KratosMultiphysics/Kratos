// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Aditya Ghantasala
// 


#if !defined(FORMULATE_CONSTRAINT_UTILITY_H_INCLUDED )
#define  FORMULATE_CONSTRAINT_UTILITY_H_INCLUDED

// System includes
#include <iostream>
#include <assert.h>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "structural_mechanics_application_variables.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/dof.h"

#include "constraint_slave.hpp"

namespace Kratos
{
///@name Kratos Globals
///@{
    
///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{
  
/** \brief InterfacePreprocessCondition 
 * Creates Model Parts containing the interface
 */

class AddMultiPointConstraint
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef ModelPart::NodesContainerType                   NodesArrayType;
    typedef ModelPart::ElementsContainerType             ElementsArrayType;
    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;
    typedef typename Element::DofsVectorType                DofsVectorType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > VariableComponentType;
    typedef Dof<double>::Pointer DofPointerType;
    typedef Dof<double> DofType;


    ///@}
    ///@name Life Cycle
    ///@{
    
    ///@}
    ///@name Operators
    ///@{
    AddMultiPointConstraint(){} //Empty Constructor
    ~AddMultiPointConstraint(){} // Empty Destructor
    void ApplyConstraint(Node<3> &MasterNode, VariableComponentType& MasterVariable, Node<3> &SlaveNode, VariableComponentType& SlaveVariable, float weight){

    	// Master dofs and weights vectors for the slave nodes.
    	MpcData &mpcData       = SlaveNode.GetValue(SLAVES);
    	bool isSlave = SlaveNode.GetValue(IS_SLAVE);

    	if(!isSlave){
    		std::cout<<"The node passed as slave is not flagged slave ... "<< isSlave<< std::endl;;
    		return;
    	} else {
    		//std::cout<<"Formulating constraint with :: "<<MasterVariable<<" as Master and :: "<<SlaveVariable<<" as Slave"<<std::endl;
    	}

    	// Adding information about the slave DOFS to the Slave node
    	DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
    	DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);

    	mpcData.AddSlave(pointerSlaveDOF);
    	mpcData.AddMaster(pointerSlaveDOF, pointerMasterDOF, weight);

    	//mpcData.GetInfo();
    }
    ///@}
    ///@name Operations
    ///@{


    
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}
private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class InterfacePreprocessCondition
}

#endif
