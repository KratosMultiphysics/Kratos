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


#if !defined( KRATOS_PROCESS_INFO_H_INCLUDED )
#define  KRATOS_PROCESS_INFO_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/data_value_container.h"
#include "includes/variables.h"


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
  class ProcessInfo : public DataValueContainer
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of ProcessInfo
      KRATOS_CLASS_POINTER_DEFINITION(ProcessInfo);

      typedef DataValueContainer BaseType;

      typedef std::size_t SizeType;
  
      typedef std::size_t IndexType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ProcessInfo() :
	DataValueContainer(),
	mIsTimeStep(true),
	mSolutionStepIndex(),
	mpPreviousSolutionStepInfo(),
	mpPreviousTimeStepInfo() 
	{
	}

      /// Copy constructor.
      ProcessInfo(const ProcessInfo& Other) : 
	DataValueContainer(Other),
	mIsTimeStep(Other.mIsTimeStep),
	mSolutionStepIndex(Other.mSolutionStepIndex),
	mpPreviousSolutionStepInfo(Other.mpPreviousSolutionStepInfo),
	mpPreviousTimeStepInfo(Other.mpPreviousTimeStepInfo) 
	{
	}

      /// Destructor.
      virtual ~ProcessInfo(){}
      

      ///@}
      ///@name Operators 
      ///@{

      /// Assignment operator.
      ProcessInfo& operator=(const ProcessInfo& rOther)
	{
	  BaseType::operator=(rOther);
	  
	  mIsTimeStep = rOther.mIsTimeStep;
	  mSolutionStepIndex = rOther.mSolutionStepIndex;
	  mpPreviousSolutionStepInfo = rOther.mpPreviousSolutionStepInfo;
	  mpPreviousTimeStepInfo = rOther.mpPreviousTimeStepInfo;

	  return *this;
	}
     
      
      ///@}
      ///@name Time steps
      ///@{

      void CreateTimeStepInfo(IndexType SolutionStepIndex = 0)
	{
	  CreateSolutionStepInfo(SolutionStepIndex);
	  mIsTimeStep = true;
	}

      void CloneTimeStepInfo(IndexType SolutionStepIndex = 0)
	{
	  CreateSolutionStepInfo(SolutionStepIndex);
	  mIsTimeStep = true;
	}

/*       void CloneTimeStepInfo(IndexType SolutionStepIndex, IndexType SourceSolutionStepIndex) */
/* 	{ */
/* 	  CloneSolutionStepInfo(SolutionStepIndex, SourceSolutionStepIndex); */
/* 	  mIsTimeStep = true; */
/* 	} */

      void CloneTimeStepInfo(IndexType SolutionStepIndex, ProcessInfo const &  SourceSolutionStepInfo)
	{
	  CloneSolutionStepInfo(SolutionStepIndex, SourceSolutionStepInfo);
	  mIsTimeStep = true;
	}

      void CreateTimeStepInfo(double NewTime, IndexType SolutionStepIndex = 0)
	{
	  CreateTimeStepInfo(SolutionStepIndex);
	  SetCurrentTime(NewTime);
	}

      void CloneTimeStepInfo(double NewTime, IndexType SourceSolutionStepIndex = 0)
	{
	  CloneTimeStepInfo(SourceSolutionStepIndex);
	  SetCurrentTime(NewTime);
	}

/*       void CloneTimeStepInfo(double NewTime, IndexType SolutionStepIndex = 0) */
/* 	{ */
/* 	  CloneTimeStepInfo(SolutionStepIndex); */
/* 	  SetCurrentTime(NewTime); */
/* 	} */

/*       void CloneTimeStepInfo(double NewTime, IndexType SolutionStepIndex = 0, IndexType SourceSolutionStepIndex = 0) */
/* 	{ */
/* 	  CloneTimeStepInfo(SolutionStepIndex, SourceSolutionStepIndex); */
/* 	  SetCurrentTime(NewTime); */
/* 	} */

      void CloneTimeStepInfo(double NewTime, IndexType SolutionStepIndex, ProcessInfo const &  SourceSolutionStepInfo)
	{
	  CloneTimeStepInfo(SolutionStepIndex, SourceSolutionStepInfo);
	  SetCurrentTime(NewTime);
	}

      void SetAsTimeStepInfo()
	{
	  mIsTimeStep = true;
          SetCurrentTime((*this)(TIME));
	}

      void SetAsTimeStepInfo(double NewTime)
	{
	  mIsTimeStep = true;
	  SetCurrentTime(NewTime);
	}

      ProcessInfo::Pointer pGetPreviousTimeStepInfo(IndexType StepsBefore = 1)
	{
	  if(StepsBefore > 1)
	    return mpPreviousTimeStepInfo->pGetPreviousTimeStepInfo(--StepsBefore);
	  
	  if(StepsBefore == 0)
	    KRATOS_ERROR(std::invalid_argument, "Steps before = 0", "");
	    
	  if(!mpPreviousTimeStepInfo)
	    KRATOS_ERROR(std::invalid_argument, "No previous time step exist.", "");

	  return mpPreviousTimeStepInfo;
	    
	}
      
      const ProcessInfo::Pointer pGetPreviousTimeStepInfo(IndexType StepsBefore = 1) const
	{
	  if(StepsBefore > 1)
	    return mpPreviousTimeStepInfo->pGetPreviousTimeStepInfo(--StepsBefore);
	  
	  if(StepsBefore == 0)
	    KRATOS_ERROR(std::invalid_argument, "Steps before = 0", "");
	    
	  if(!mpPreviousTimeStepInfo)
	    KRATOS_ERROR(std::invalid_argument, "No previous time step exist.", "");

	  return mpPreviousTimeStepInfo;
	    
	}
      
      ProcessInfo& GetPreviousTimeStepInfo(IndexType StepsBefore = 1)
	{
	  return *pGetPreviousTimeStepInfo(StepsBefore);
	    
	}
      
      ProcessInfo const& GetPreviousTimeStepInfo(IndexType StepsBefore = 1) const
	{
	  return *pGetPreviousTimeStepInfo(StepsBefore);
	    
	}
      
      IndexType GetPreviousTimeStepIndex(IndexType StepsBefore = 1) const
	{
	  return GetPreviousTimeStepInfo(StepsBefore).GetSolutionStepIndex();
	}

      ///@}
      ///@name Operations
      ///@{

      void SetCurrentTime(double NewTime)
      {
	  (*this)(TIME) = NewTime;
	  if(!mpPreviousTimeStepInfo)
	    (*this)(DELTA_TIME) = NewTime;
	  else
	    (*this)(DELTA_TIME) = NewTime -  mpPreviousTimeStepInfo->GetValue(TIME);
      }

      void ReIndexBuffer(SizeType BufferSize, IndexType BaseIndex = 0)
	{
	  mSolutionStepIndex = BaseIndex;
	  if(BufferSize > 1)
	    if(mpPreviousSolutionStepInfo)
	      mpPreviousSolutionStepInfo->ReIndexBuffer(BufferSize - 1, BaseIndex + 1);
	} 
      
      
      ///@}
      ///@name Solution Step Data
      ///@{
      
      void CreateSolutionStepInfo(IndexType SolutionStepIndex = 0)
	{
	  mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
	  mSolutionStepIndex = SolutionStepIndex;
	  if(mIsTimeStep)
	    mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
	  mIsTimeStep = false;
	  Clear();
	}

//       void CloneSolutionStepInfo(IndexType SolutionStepIndex = 0)
// 	{
// 	  mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
// 	  mSolutionStepIndex = SolutionStepIndex;
// 	  if(mIsTimeStep)
// 	    mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
// 	  mIsTimeStep = false;
// 	}

//       void CloneSolutionStepInfo(IndexType SolutionStepIndex = 0, IndexType SourceSolutionStepIndex = 0)
// 	{
// 	  ProcessInfo& source_info = FindSolutionStepInfo(SourceSolutionStepIndex);
// 	  if(source_info.GetSolutionStepIndex() == SourceSolutionStepIndex)
// 	    {
// 	      mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
// 	      mSolutionStepIndex = SolutionStepIndex;
// 	      BaseType::operator=(source_info);
// 	      if(mIsTimeStep)
// 		mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
// 	      mIsTimeStep = false;
// 	    }
// 	  else
// 	    CreateSolutionStepInfo(SolutionStepIndex);
// 	}

      void CloneSolutionStepInfo()
	{
	  mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
	  mSolutionStepIndex = 0;
	  if(mIsTimeStep)
	    mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
	  mIsTimeStep = false;
	}

      void CloneSolutionStepInfo(IndexType SourceSolutionStepIndex)
	{
	  ProcessInfo& source_info = FindSolutionStepInfo(SourceSolutionStepIndex);
	  if(source_info.GetSolutionStepIndex() == SourceSolutionStepIndex)
	    {
	      mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
	      mSolutionStepIndex = 0;
	      BaseType::operator=(source_info);
	      if(mIsTimeStep)
		mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
	      mIsTimeStep = false;
	    }
	  else
	    CreateSolutionStepInfo(0);
	}

      void CloneSolutionStepInfo(IndexType SolutionStepIndex, ProcessInfo const &  SourceSolutionStepInfo)
	{
	      mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
	      mSolutionStepIndex = SolutionStepIndex;
	      BaseType::operator=(SourceSolutionStepInfo);
	      if(mIsTimeStep)
		mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
	      mIsTimeStep = false;
	}
      
      ProcessInfo& FindSolutionStepInfo(IndexType ThisIndex)
	{
	  if(mSolutionStepIndex == ThisIndex)
	    return *this;

	  if(!mpPreviousSolutionStepInfo)
	    return *this;
	  
	  return mpPreviousTimeStepInfo->FindSolutionStepInfo(ThisIndex);
	    
	}

      void RemoveSolutionStepInfo(IndexType SolutionStepIndex)
	{
	  if(!mpPreviousSolutionStepInfo)
	    return;
	  
	  if(mpPreviousSolutionStepInfo->GetSolutionStepIndex() == SolutionStepIndex)
	    mpPreviousSolutionStepInfo = mpPreviousSolutionStepInfo->pGetPreviousSolutionStepInfo();
	  else
	    mpPreviousSolutionStepInfo->RemoveSolutionStepInfo(SolutionStepIndex);
	}

      void ClearHistory(IndexType StepsBefore = 0)
	{
	  if(StepsBefore == 0)
	    {
	      mpPreviousTimeStepInfo = Pointer();
	      mpPreviousSolutionStepInfo = Pointer();
	    }
	  else 
	    {
	      if(mpPreviousTimeStepInfo)
		mpPreviousTimeStepInfo->ClearHistory(--StepsBefore);
	      if(mpPreviousSolutionStepInfo)
		mpPreviousSolutionStepInfo->ClearHistory(--StepsBefore);
	    }
	}

      ProcessInfo::Pointer pGetPreviousSolutionStepInfo(IndexType StepsBefore = 1)
	{
	  if(StepsBefore > 1)
	    return mpPreviousSolutionStepInfo->pGetPreviousSolutionStepInfo(--StepsBefore);
	  
	  if(StepsBefore == 0)
	    KRATOS_ERROR(std::invalid_argument, "Steps before = 0", "");
	    
	  if(!mpPreviousSolutionStepInfo)
	    KRATOS_ERROR(std::invalid_argument, "No previous time step exist.", "");

	  return mpPreviousSolutionStepInfo;
	    
	}
      
      const ProcessInfo::Pointer pGetPreviousSolutionStepInfo(IndexType StepsBefore = 1) const
	{
	  if(StepsBefore > 1)
	    return mpPreviousSolutionStepInfo->pGetPreviousSolutionStepInfo(--StepsBefore);
	  
	  if(StepsBefore == 0)
	    KRATOS_ERROR(std::invalid_argument, "Steps before = 0", "");
	    
	  if(!mpPreviousSolutionStepInfo)
	    KRATOS_ERROR(std::invalid_argument, "No previous time step exist.", "");

	  return mpPreviousSolutionStepInfo;
	    
	}
      
      ProcessInfo& GetPreviousSolutionStepInfo(IndexType StepsBefore = 1)
	{
	  return *pGetPreviousSolutionStepInfo(StepsBefore);
	    
	}
      
      ProcessInfo const& GetPreviousSolutionStepInfo(IndexType StepsBefore = 1) const
	{
	  return *pGetPreviousSolutionStepInfo(StepsBefore);
	    
	}
      
      IndexType GetPreviousSolutionStepIndex(IndexType StepsBefore = 1) const
	{
	  return GetPreviousSolutionStepInfo(StepsBefore).GetSolutionStepIndex();
	}

      ///@}
      ///@name Access
      ///@{ 
      
      IndexType GetSolutionStepIndex() const
	{
	  return mSolutionStepIndex;
	}
      
      void SetSolutionStepIndex(IndexType NewIndex)
	{
	  mSolutionStepIndex = NewIndex;
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
	  return "Process Info";
	}

      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info();
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  rOStream << "    Current solution step index : " << mSolutionStepIndex << std::endl;
	  BaseType::PrintData(rOStream);
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
        
      bool mIsTimeStep;
        
      IndexType mSolutionStepIndex;

      ProcessInfo::Pointer mpPreviousSolutionStepInfo;

      ProcessInfo::Pointer mpPreviousTimeStepInfo;
        
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
      
       
      ///@}    
        
    }; // Class ProcessInfo 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    ProcessInfo& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ProcessInfo& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_PROCESS_INFO_H_INCLUDED  defined 


