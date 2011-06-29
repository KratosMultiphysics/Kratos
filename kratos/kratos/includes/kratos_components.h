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

//   Last Modified by:    $Author: pooyan $

//   Date:                $Date: 2008-10-29 14:31:24 $

//   Revision:            $Revision: 1.4 $

//

//





#if !defined(KRATOS_KRATOS_COMPONENTS_H_INCLUDED )

#define  KRATOS_KRATOS_COMPONENTS_H_INCLUDED







// System includes

#include <string>

#include <iostream>

#include <map> 





// External includes 

#include <boost/ref.hpp>



// Project includes

#include "includes/define.h"

#include "containers/variable_data.h"





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

  template<class TComponentType>

  class KratosComponents

    {

    public:

      ///@name Type Definitions

      ///@{

      

      /// Pointer definition of KratosComponents

      KRATOS_CLASS_POINTER_DEFINITION(KratosComponents);



      typedef std::map<std::string, const TComponentType* > ComponentsContainerType;

  

      typedef typename ComponentsContainerType::value_type ValueType;



      ///@}

      ///@name Life Cycle 

      ///@{ 

      

      /// Default constructor.

      KratosComponents(){}



//       KratosComponents(std::string const& Name, TComponentType const& ThisComponent)

//       {

//  	msComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));

//       }



      /// Destructor.

      virtual ~KratosComponents(){}

      



      ///@}

      ///@name Operators 

      ///@{

      

      

      ///@}

      ///@name Operations

      ///@{



      static void Add(std::string const& Name, TComponentType const& ThisComponent)

      {

	msComponents.insert(typename ComponentsContainerType::value_type(Name , &ThisComponent));

      }

      

 //     static void Add(std::string const& Name, TComponentType const& ThisComponent, ComponentsContainerType& ThisComponents)

 //     {

	////ThisComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));

	//msComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));

 //     }

      

      static TComponentType const& Get(std::string const& Name)

      {
	typename ComponentsContainerType::iterator i =  msComponents.find(Name);
/* 	if(i == msComponents.end()) */
/* 	   KRATOS_ERROR(std::invalid_argument, "The component is not registered!", Name); */
	return *(i->second);

      }

      

      static ComponentsContainerType & GetComponents()

      {

	return msComponents;

      }

      

	  static ComponentsContainerType * pGetComponents()

      {

	return &msComponents;

      }



      static void Register()

	{

	}

      

      

      ///@}

      ///@name Access

      ///@{ 

      

      

      ///@}

      ///@name Inquiry

      ///@{

      
      static bool Has(std::string const& Name)

      {

	return (msComponents.find(Name) != msComponents.end());

      }


      

      ///@}      

      ///@name Input and output

      ///@{



      /// Turn back information as a string.

      virtual std::string Info() const

      {

	return "Kratos components";

      }

      

      /// Print information about this object.

      virtual void PrintInfo(std::ostream& rOStream) const

      {

	rOStream << "Kratos components";

      }



      /// Print object's data.

      virtual void PrintData(std::ostream& rOStream) const

      {

	for(typename ComponentsContainerType::const_iterator i = msComponents.begin() ; i != msComponents.end() ; ++i)

	  rOStream << "    " << i->first << std::endl;

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



      static ComponentsContainerType msComponents;

        

        

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

      KratosComponents& operator=(KratosComponents const& rOther);



      /// Copy constructor.

      KratosComponents(KratosComponents const& rOther);



        

      ///@}    

        

    }; // Class KratosComponents 



  ///@} 



  template<>

  class KratosComponents<VariableData>

    {

    public:

      ///@name Type Definitions

      ///@{

      

      /// Pointer definition of KratosComponents

      KRATOS_CLASS_POINTER_DEFINITION(KratosComponents);



      typedef std::map<std::string, VariableData* > ComponentsContainerType;

  

      typedef ComponentsContainerType::value_type ValueType;



      ///@}

      ///@name Life Cycle 

      ///@{ 

      

      /// Default constructor.

      KratosComponents(){}



//       KratosComponents(std::string const& Name, TComponentType const& ThisComponent)

//       {

//  	msComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));

//       }



      /// Destructor.

      virtual ~KratosComponents(){}

      



      ///@}

      ///@name Operators 

      ///@{

      

      

      ///@}

      ///@name Operations

      ///@{



      static void Add(std::string const& Name, VariableData& ThisComponent)

      {

	msComponents.insert(ComponentsContainerType::value_type(Name ,&ThisComponent));

      }

      

	  static std::size_t Size() 

      {

		return msComponents.size();

      }

      

 //     static void Add(std::string const& Name, VariableData& ThisComponent, ComponentsContainerType& ThisComponents)

 //     {

	////ThisComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));

	//msComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::ref(ThisComponent)));

 //     }

      

      static VariableData & Get(std::string const& Name)

      {

	return *(msComponents.find(Name)->second);

      }

      static VariableData* pGet(std::string const& Name)

      {

	return (msComponents.find(Name)->second);

      }

      

      static ComponentsContainerType & GetComponents()

      {

	return msComponents;

      }

      

	  static ComponentsContainerType * pGetComponents()

      {

	return &msComponents;

      }



      static void Register()

	{

	}

      

      

      ///@}

      ///@name Access

      ///@{ 

      

      

      ///@}

      ///@name Inquiry

      ///@{

      
      static bool Has(std::string const& Name)
      {
	return (msComponents.find(Name) != msComponents.end());
      }


      

      ///@}      

      ///@name Input and output

      ///@{



      /// Turn back information as a string.

      virtual std::string Info() const

      {

	return "Kratos components";

      }

      

      /// Print information about this object.

      virtual void PrintInfo(std::ostream& rOStream) const

      {

	rOStream << "Kratos components";

      }



      /// Print object's data.

      virtual void PrintData(std::ostream& rOStream) const

      {

	for(ComponentsContainerType::const_iterator i = msComponents.begin() ; i != msComponents.end() ; ++i)

	  rOStream << "    " << *(i->second) << std::endl;

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



      static ComponentsContainerType msComponents;

        

        

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

      KratosComponents& operator=(KratosComponents const& rOther);



      /// Copy constructor.

      KratosComponents(KratosComponents const& rOther);



        

      ///@}    

        

    }; // Class KratosComponents 



  template<class TComponentType>

  typename KratosComponents<TComponentType>::ComponentsContainerType KratosComponents<TComponentType>::msComponents;

  

  

  ///@name Type Definitions       

  ///@{ 

  

  

  ///@} 

  ///@name Input and output 

  ///@{ 

        

 

  /// input stream function

//   template<class TComponentType>

//   inline std::istream& operator >> (std::istream& rIStream, 

// 				    KratosComponents<TComponentType>& rThis);



  /// output stream function

  template<class TComponentType>

  inline std::ostream& operator << (std::ostream& rOStream, 

				    const KratosComponents<TComponentType>& rThis)

    {

      rThis.PrintInfo(rOStream);

      rOStream << std::endl;

      rThis.PrintData(rOStream);



      return rOStream;

    }

  ///@} 



  template<class TComponentType> inline void AddComponent(std::string const& Name, TComponentType const& ThisComponent)

  {

    KratosComponents<TComponentType>::Add(Name, ThisComponent);

  }

  

}  // namespace Kratos.



#endif // KRATOS_KRATOS_COMPONENTS_H_INCLUDED  defined 





