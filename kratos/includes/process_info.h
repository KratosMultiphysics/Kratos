//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
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

/// ProcessInfo holds the current value of different solution parameters.
/**
 * ProcessInfo holds the current value of different solution parameters.
 * It can be used to keep variables like time, solution step, non linear step, or any other variable defined in Kratos.
 * Its variable base interface provides a clear and flexible access to these data.
*/
class KRATOS_API(KRATOS_CORE) ProcessInfo : public DataValueContainer, public Flags
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
    ~ProcessInfo() override {}


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

    void SetAsTimeStepInfo();

    void SetAsTimeStepInfo(double NewTime);

    ProcessInfo::Pointer pGetPreviousTimeStepInfo(IndexType StepsBefore = 1);

    const ProcessInfo::Pointer pGetPreviousTimeStepInfo(IndexType StepsBefore = 1) const;

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

    void SetCurrentTime(double NewTime);

    void ReIndexBuffer(SizeType BufferSize, IndexType BaseIndex = 0);


    ///@}
    ///@name Solution Step Data
    ///@{

    void CreateSolutionStepInfo(IndexType SolutionStepIndex = 0);

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

    void CloneSolutionStepInfo();

    void CloneSolutionStepInfo(IndexType SourceSolutionStepIndex);

    void CloneSolutionStepInfo(IndexType SolutionStepIndex, ProcessInfo const &  SourceSolutionStepInfo);

    ProcessInfo& FindSolutionStepInfo(IndexType ThisIndex);

    void RemoveSolutionStepInfo(IndexType SolutionStepIndex);

    void ClearHistory(IndexType StepsBefore = 0);

    ProcessInfo::Pointer pGetPreviousSolutionStepInfo(IndexType StepsBefore = 1);

    const ProcessInfo::Pointer pGetPreviousSolutionStepInfo(IndexType StepsBefore = 1) const;

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
    std::string Info() const override
    {
        return "Process Info";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

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


