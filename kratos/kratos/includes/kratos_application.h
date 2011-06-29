/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_KRATOS_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"
#include "includes/element.h"
#include "includes/condition.h"

namespace Kratos
{
  ///@name Kratos Classes
  ///@{

  /// This class defines the interface with kernel for all applications in Kratos.
  /** The application class defines the interface necessary for providing the information
      needed by Kernel in order to configure the whole sistem correctly.

  */

  class KratosApplication
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of KratosApplication
      KRATOS_CLASS_POINTER_DEFINITION(KratosApplication);

  
      ///@}
      ///@name Life Cycle 
      ///@{ 
     
      /// Default constructor.
	KratosApplication();


      /// Copy constructor.
	KratosApplication(KratosApplication const& rOther) : 
			mpVariableData(rOther.mpVariableData),
			mpDoubleVariables(rOther.mpDoubleVariables),
			mpArray1DVariables(rOther.mpArray1DVariables),
			mpVectorVariables(rOther.mpVectorVariables),
			mpMatrixVariables(rOther.mpMatrixVariables),
			mpArray1DVariableComponents(rOther.mpArray1DVariableComponents),
			mpElements(rOther.mpElements), 
			mpConditions(rOther.mpConditions){}



      /// Destructor.
      virtual ~KratosApplication(){}

      

      ///@}
      ///@name Operations
      ///@{

      virtual void Register()

	  {

		  RegisterVariables();

	  }

      
     void RegisterVariables();

      ///@}
      ///@name Access
      ///@{



//	template<class TComponentType>
//		typename KratosComponents<TComponentType>::ComponentsContainerType& GetComponents(TComponentType const& rComponentType)
//	{
//		return KratosComponents<TComponentType>::GetComponents();
//	}



	

	

	// I have to see why the above version is not working for multi thread ... 
	// Anyway its working with these functions.Pooyan.
	KratosComponents<Variable<double> >::ComponentsContainerType& GetComponents(Variable<double> const& rComponentType)
	{
		return *mpDoubleVariables;
	}

	KratosComponents<Variable<array_1d<double, 3> > >::ComponentsContainerType& GetComponents(Variable<array_1d<double, 3> >  const& rComponentType)
	{
		return *mpArray1DVariables;
	}

	KratosComponents<Variable<Vector> >::ComponentsContainerType& GetComponents(Variable<Vector> const& rComponentType)
	{
		return *mpVectorVariables;
	}

	KratosComponents<Variable<Matrix> >::ComponentsContainerType& GetComponents(Variable<Matrix>  const& rComponentType)
	{
		return *mpMatrixVariables;
	}

	KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::ComponentsContainerType& GetComponents(VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > const& rComponentType)
	{
		return *mpArray1DVariableComponents;
	}
	
	KratosComponents<VariableData>::ComponentsContainerType& GetVariables()
	{
	  return *mpVariableData;
	  
	}
	
	KratosComponents<Element>::ComponentsContainerType& GetElements()
	{
	  return *mpElements;
	}
	
	KratosComponents<Condition>::ComponentsContainerType& GetConditions()
	{
	  return *mpConditions;
	}
	
	void SetComponents(KratosComponents<VariableData>::ComponentsContainerType const& VariableDataComponents)
	
	{
	  for(KratosComponents<VariableData>::ComponentsContainerType::iterator i = mpVariableData->begin() ;
	  
	  i != mpVariableData->end() ; i++)
	  
	  {
	    std::string const& variable_name = i->second->Name();
	    KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = VariableDataComponents.find(variable_name);
	    
	    if(i_variable == VariableDataComponents.end()) 
	      
	      KRATOS_ERROR(std::logic_error, "This variable is not registered in Kernel : ",   *(i_variable->second));
	    
	    unsigned int variable_key = i_variable->second->Key();
	    
	    if(variable_key == 0)
	      
	      KRATOS_ERROR(std::logic_error, "This variable is not initialized in Kernel : ",   *(i_variable->second));
	    
	    
	    
	    //			KRATOS_WATCH(i_variable->second.get());
	    
	    //			KRATOS_WATCH(i->second.get().Key());
	    
	    //			KRATOS_WATCH(variable_key);
	    
	    i->second->SetKey(variable_key);
	    
	  }
	  
	  //			KRATOS_WATCH("!!!!!!!!!!!!!!!!!!!!! END SETTING COMPONENETS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	  
	}
	
	
	
	
	void SetComponents(KratosComponents<Element>::ComponentsContainerType const& ElementComponents)
	
	{
	  // It's better to make a loop over new components and add them if they are NOT already exist in application. Or make an ERROR for incompatibility between applications.
	  
	  mpElements->insert(ElementComponents.begin(), ElementComponents.end());
	  
	}
	
	
	
	void SetComponents(KratosComponents<Condition>::ComponentsContainerType const& ConditionComponents)
	
	{
	  
	  mpConditions->insert(ConditionComponents.begin(), ConditionComponents.end());
	  
	}
	

	Serializer::RegisteredObjectsContainerType& GetRegisteredObjects()
        {
            return *mpRegisteredObjects;
        }

	Serializer::RegisteredObjectsNameContainerType& GetRegisteredObjectsName()
        {
            return *mpRegisteredObjectsName;
        }
      

      ///@}

      ///@name Inquiry

      ///@{

      

      

      ///@}      

      ///@name Input and output

      ///@{



      /// Turn back information as a string.

      virtual std::string Info() const

      {

	return "KratosApplication";

      }

      

      /// Print information about this object.

      virtual void PrintInfo(std::ostream& rOStream) const

      {

	rOStream << Info();

      }



      /// Print object's data.

      virtual void PrintData(std::ostream& rOStream) const

      {

		rOStream << "Variables:" << std::endl;

		KratosComponents<VariableData>().PrintData(rOStream);

		rOStream << std::endl;

		rOStream << "Elements:" << std::endl;

		KratosComponents<Element>().PrintData(rOStream);

		rOStream << std::endl;

		rOStream << "Conditions:" << std::endl;

		KratosComponents<Condition>().PrintData(rOStream);

      }

      

            

      ///@}      

      ///@name Friends

      ///@{

      

            

      ///@}

      

    protected:

      ///@name Protected static Member Variables 

      ///@{ 

       

      ///@} 

      ///@name Protected member Variables 

      ///@{ 
		const Condition  mCondition3D; 

		const Condition  mCondition2D;        

        

		KratosComponents<VariableData>::ComponentsContainerType* mpVariableData;

		KratosComponents<Variable<double> >::ComponentsContainerType* mpDoubleVariables;

		KratosComponents<Variable<array_1d<double, 3> > >::ComponentsContainerType* mpArray1DVariables;

		KratosComponents<Variable<Vector> >::ComponentsContainerType* mpVectorVariables;

		KratosComponents<Variable<Matrix> >::ComponentsContainerType* mpMatrixVariables;
		
		KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::ComponentsContainerType* mpArray1DVariableComponents; 

		KratosComponents<Element>::ComponentsContainerType* mpElements;

		KratosComponents<Condition>::ComponentsContainerType* mpConditions;

                Serializer::RegisteredObjectsContainerType* mpRegisteredObjects;

                Serializer::RegisteredObjectsNameContainerType* mpRegisteredObjectsName;

        

        

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

      ///@name Private Inquiry 

      ///@{ 

        

        

      ///@}    

      ///@name Un accessible methods 

      ///@{ 

      

      /// Assignment operator.

      KratosApplication& operator=(KratosApplication const& rOther);



        

      ///@}    

        

    }; // Class KratosApplication 



  ///@} 

  

  ///@name Type Definitions       

  ///@{ 

  

  

  ///@} 

  ///@name Input and output 

  ///@{ 

        

 

  /// input stream function

  inline std::istream& operator >> (std::istream& rIStream, 

				    KratosApplication& rThis);



  /// output stream function

  inline std::ostream& operator << (std::ostream& rOStream, 

				    const KratosApplication& rThis)

    {

      rThis.PrintInfo(rOStream);

      rOStream << std::endl;

      rThis.PrintData(rOStream);



      return rOStream;

    }

  ///@} 

  

  

}  // namespace Kratos.



#endif // KRATOS_KRATOS_APPLICATION_H_INCLUDED  defined 





