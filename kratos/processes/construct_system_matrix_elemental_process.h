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
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_CONSTRUCT_SYSTEM_MATRIX_ELEMENTAL_PROCESS_H_INCLUDED )
#define  KRATOS_CONSTRUCT_SYSTEM_MATRIX_ELEMENTAL_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "containers/pointer_vector_set.h"


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
  template<class TSystemSpaceType, class TDirichletSpaceType = TSystemSpaceType>
    class ConstructSystemMatrixElementalProcess : public Process
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of ConstructSystemMatrixElementalProcess
      KRATOS_CLASS_POINTER_DEFINITION(ConstructSystemMatrixElementalProcess);

      typedef typename TSystemSpaceType::MatrixType SystemMatrixType;
      
      typedef typename TDirichletSpaceType::MatrixType DirichletMatrixType;
    
      typedef Element ElementType;

      typedef PointerVectorSet<ElementType, IndexedObject> ElementsContainerType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Constructor.
      ConstructSystemMatrixElementalProcess(ElementsContainerType & rElements, SystemMatrixType& rSystemMatrix, DirichletMatrixType& rDirichletMatrix)
      : mElements(rElements), mSystemMatrix(rSystemMatrix), mDirichletMatrix(rDirichletMatrix){}

      /// Destructor.
      virtual ~ConstructSystemMatrixElementalProcess(){}
      

      ///@}
      ///@name Operators 
      ///@{

      
      ///@}
      ///@name Operations
      ///@{

      virtual void Execute()
	{
	  std::size_t equation_size = TSystemSpaceType::Size1(mSystemMatrix);
	  std::vector<std::vector<std::size_t> > indices(equation_size);
	  std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));

	  ProcessInfo process_info;
	  for(ElementsContainerType::iterator i_element = mElements.begin() ; i_element != mElements.end() ; i_element++)
	  {
	    Element::EquationIdVectorType ids(3,0);
	    
	    i_element->EquationIdVector(ids, process_info);
	    for(std::size_t i = 0 ; i < 3 ; i++)
	      if(ids[i] < equation_size)
		for(std::size_t j = 0 ; j < 3 ; j++)
		  if(ids[j] < equation_size)
		    indices[ids[i]].push_back(ids[j]);
		  else
		    dirichlet_indices[ids[i]].push_back(ids[j]- equation_size);
	    
	  }
	  for(std::size_t i = 0 ; i < indices.size() ; i++)
	  {
	    std::vector<std::size_t>& row_indices = indices[i];
	    std::vector<std::size_t>& dirichlet_row_indices = dirichlet_indices[i];
	    std::sort(row_indices.begin(), row_indices.end());
	    std::unique(row_indices.begin(), row_indices.end());
	    std::sort(dirichlet_row_indices.begin(), dirichlet_row_indices.end());
	    std::unique(dirichlet_row_indices.begin(), dirichlet_row_indices.end());
	    for(std::size_t j = 0 ; j < row_indices.size() ; j++)
	    mSystemMatrix()(i,row_indices[j]) = 0.00;
	    for(std::size_t j = 0 ; j < dirichlet_row_indices.size() ; j++)
	    mDirichletMatrix()(i,dirichlet_row_indices[j]) = 0.00;
	  }
	  for(std::size_t i = 0 ; i < indices.size() ; i++) // To make put at least one element per each row. Pooyan.
	    mDirichletMatrix()(i,0) = 0.00;
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
	  return "ConstructSystemMatrixElementalProcess";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "ConstructSystemMatrixElementalProcess";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  rOStream << "    System matrix non zeros : " << mSystemMatrix().non_zeros() << std::endl;
	  rOStream << "    Dirichlet matrix non zeros : " << mDirichletMatrix().non_zeros() << std::endl;
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
        
      ElementsContainerType& mElements;
      SystemMatrixType& mSystemMatrix;
      DirichletMatrixType& mDirichletMatrix;
        
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
      ConstructSystemMatrixElementalProcess& operator=(ConstructSystemMatrixElementalProcess const& rOther);

      /// Copy constructor.
      ConstructSystemMatrixElementalProcess(ConstructSystemMatrixElementalProcess const& rOther);

        
      ///@}    
        
    }; // Class ConstructSystemMatrixElementalProcess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function

  template<class TSystemSpaceType, class TDirichletSpaceType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    ConstructSystemMatrixElementalProcess<TSystemSpaceType, TDirichletSpaceType>& rThis);

  /// output stream function
  template<class TSystemSpaceType, class TDirichletSpaceType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ConstructSystemMatrixElementalProcess<TSystemSpaceType, TDirichletSpaceType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_CONSTRUCT_SYSTEM_MATRIX_ELEMENTAL_PROCESS_H_INCLUDED  defined 


