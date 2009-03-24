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
//   Date:                $Date: 2008-01-30 10:10:35 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_KERNEL_H_INCLUDED )
#define  KRATOS_KERNEL_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"


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
  
  /// Short class definition.
  /** Detail class definition.
  */
  class Kernel
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Kernel
      KRATOS_CLASS_POINTER_DEFINITION(Kernel);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Kernel()
	  {
		  mKratosApplication.RegisterVariables();
//		  		KratosApplication::RegisterVariables();
	  }

      /// Copy constructor.
	  Kernel(Kernel const& rOther){}


      /// Destructor.
      virtual ~Kernel(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      void Initialize()
	{
	  unsigned int j = 0;
	  for(KratosComponents<VariableData>::ComponentsContainerType::iterator i = KratosComponents<VariableData>::GetComponents().begin() ;
		  i != KratosComponents<VariableData>::GetComponents().end() ; i++)
		  //const_cast<VariableData&>(i->second.get()).SetKey(++j);
		  i->second.get().SetKey(++j);
	}
      
      void AddApplication(KratosApplication& NewApplication)
	{
	typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

		//NewApplication.SetComponents(KratosComponents<VariableData>::pGetComponents(), 
		//	KratosComponents<Element>::pGetComponents(), 
		//	KratosComponents<Condition>::pGetComponents());
			KRATOS_WATCH("kerneal entered in AddApplication");
		NewApplication.Register();
			KRATOS_WATCH("Application Registered");
		KratosComponents<VariableData>::GetComponents().insert(NewApplication.GetVariables().begin(),
															   NewApplication.GetVariables().end());
			KRATOS_WATCH("Variables Registered");

		KratosComponents<Variable<double> >::GetComponents().insert(NewApplication.GetComponents(Variable<double>("NONE")).begin(),
															   NewApplication.GetComponents(Variable<double>("NONE")).end());
		KratosComponents<Variable<array_1d<double, 3> > >::GetComponents().insert(NewApplication.GetComponents(Variable<array_1d<double, 3> >("NONE")).begin(),
															   NewApplication.GetComponents(Variable<array_1d<double, 3> >("NONE")).end());
		KratosComponents<Variable<Vector> >::GetComponents().insert(NewApplication.GetComponents(Variable<Vector>("NONE")).begin(),
															   NewApplication.GetComponents(Variable<Vector>("NONE")).end());
		KratosComponents<Variable<Matrix> >::GetComponents().insert(NewApplication.GetComponents(Variable<Matrix>("NONE")).begin(),
															   NewApplication.GetComponents(Variable<Matrix>("NONE")).end());
		Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > temp_adaptor(DISPLACEMENT, 0); // the displacement is not important, only an array_1d variable is needed!
		KratosComponents<array_1d_component_type>::GetComponents().insert(NewApplication.GetComponents(array_1d_component_type("NONE", temp_adaptor)).begin(),				   NewApplication.GetComponents(array_1d_component_type("NONE", temp_adaptor)).end());

	
		KratosComponents<Element>::GetComponents().insert(NewApplication.GetElements().begin(),
														  NewApplication.GetElements().end());
			KRATOS_WATCH("Elements Registered");
		KratosComponents<Condition>::GetComponents().insert(NewApplication.GetConditions().begin(),
															NewApplication.GetConditions().end());
			KRATOS_WATCH("Conditions Registered");

		KratosComponents<Variable<double> >::GetComponents().insert(NewApplication.GetComponents(Variable<double>("NONE")).begin(),
															   NewApplication.GetComponents(Variable<double>("NONE")).end());

	}
      
      void InitializeApplication(KratosApplication& NewApplication)
	{
		NewApplication.SetComponents(KratosComponents<VariableData>::GetComponents());
		NewApplication.SetComponents(KratosComponents<Element>::GetComponents());
		NewApplication.SetComponents(KratosComponents<Condition>::GetComponents());
	}
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  return "kernel";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "kernel";
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
		rOStream << "Echo Finished" << std::endl;
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
	  KratosApplication mKratosApplication;
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
      void RegisterVariables();
        
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
      Kernel& operator=(Kernel const& rOther);

        
      ///@}    
        
    }; // Class Kernel 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Kernel& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Kernel& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_KERNEL_H_INCLUDED  defined 


