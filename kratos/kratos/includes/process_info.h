// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
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
#include "containers/flags.h"


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
class ProcessInfo : public DataValueContainer, public Flags
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
		BaseType(),
        Flags(),
        mIsTimeStep(true),
        mSolutionStepIndex(),
        mpPreviousSolutionStepInfo(),
        mpPreviousTimeStepInfo()
    {
    }

    /// Copy constructor.
    ProcessInfo(const ProcessInfo& Other) :
		BaseType(Other),
        Flags(Other),
        mIsTimeStep(Other.mIsTimeStep),
        mSolutionStepIndex(Other.mSolutionStepIndex),
        mpPreviousSolutionStepInfo(Other.mpPreviousSolutionStepInfo),
        mpPreviousTimeStepInfo(Other.mpPreviousTimeStepInfo)
    {
    }

    /// Destructor.
    virtual ~ProcessInfo() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ProcessInfo& operator=(const ProcessInfo& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);

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
            KRATOS_THROW_ERROR(std::invalid_argument, "Steps before = 0", "");

        if(!mpPreviousTimeStepInfo)
            KRATOS_THROW_ERROR(std::invalid_argument, "No previous time step exist.", "");

        return mpPreviousTimeStepInfo;

    }

    const ProcessInfo::Pointer pGetPreviousTimeStepInfo(IndexType StepsBefore = 1) const
    {
        if(StepsBefore > 1)
            return mpPreviousTimeStepInfo->pGetPreviousTimeStepInfo(--StepsBefore);

        if(StepsBefore == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "Steps before = 0", "");

        if(!mpPreviousTimeStepInfo)
            KRATOS_THROW_ERROR(std::invalid_argument, "No previous time step exist.", "");

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
		BaseType::Clear();
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
            KRATOS_THROW_ERROR(std::invalid_argument, "Steps before = 0", "");

        if(!mpPreviousSolutionStepInfo)
            KRATOS_THROW_ERROR(std::invalid_argument, "No previous time step exist.", "");

        return mpPreviousSolutionStepInfo;

    }

    const ProcessInfo::Pointer pGetPreviousSolutionStepInfo(IndexType StepsBefore = 1) const
    {
        if(StepsBefore > 1)
            return mpPreviousSolutionStepInfo->pGetPreviousSolutionStepInfo(--StepsBefore);

        if(StepsBefore == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "Steps before = 0", "");

        if(!mpPreviousSolutionStepInfo)
            KRATOS_THROW_ERROR(std::invalid_argument, "No previous time step exist.", "");

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
    ///@name Serialization
    ///@{

    friend class Serializer;


    virtual void save(Serializer& rSerializer) const
    {
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("Is Time Step",mIsTimeStep);
        rSerializer.save("Solution Step Index",mSolutionStepIndex);
        rSerializer.save("Previous Solution Step Info",mpPreviousSolutionStepInfo);
        rSerializer.save("Previous Time Step Info", mpPreviousTimeStepInfo);
    }

    virtual void load(Serializer& rSerializer)
    {
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
        rSerializer.load("Is Time Step",mIsTimeStep);
        rSerializer.load("Solution Step Index",mSolutionStepIndex);
        rSerializer.load("Previous Solution Step Info",mpPreviousSolutionStepInfo);
        rSerializer.load("Previous Time Step Info", mpPreviousTimeStepInfo);
    }



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


