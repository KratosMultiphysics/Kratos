//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//


// System includes


// External includes


// Project includes
#include "includes/process_info.h"
#include "includes/variables.h"

namespace Kratos
{

using IndexType = ProcessInfo::IndexType;

/// Default constructor.
ProcessInfo::ProcessInfo() :
    BaseType(),
    Flags(),
    mIsTimeStep(true),
    mSolutionStepIndex(),
    mpPreviousSolutionStepInfo(),
    mpPreviousTimeStepInfo()
{
}

/// Copy constructor.
ProcessInfo::ProcessInfo(const ProcessInfo& Other) :
    BaseType(Other),
    Flags(Other),
    mIsTimeStep(Other.mIsTimeStep),
    mSolutionStepIndex(Other.mSolutionStepIndex),
    mpPreviousSolutionStepInfo(Other.mpPreviousSolutionStepInfo),
    mpPreviousTimeStepInfo(Other.mpPreviousTimeStepInfo)
{
}

/// Destructor.
ProcessInfo::~ProcessInfo() {}

/// Assignment operator.
ProcessInfo& ProcessInfo::operator=(const ProcessInfo& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator=(rOther);

    mIsTimeStep = rOther.mIsTimeStep;
    mSolutionStepIndex = rOther.mSolutionStepIndex;
    mpPreviousSolutionStepInfo = rOther.mpPreviousSolutionStepInfo;
    mpPreviousTimeStepInfo = rOther.mpPreviousTimeStepInfo;

    return *this;
}

void ProcessInfo::CreateTimeStepInfo(IndexType SolutionStepIndex)
{
    CreateSolutionStepInfo(SolutionStepIndex);
    mIsTimeStep = true;
}

void ProcessInfo::CloneTimeStepInfo(IndexType SolutionStepIndex)
{
    CreateSolutionStepInfo(SolutionStepIndex);
    mIsTimeStep = true;
}

void ProcessInfo::CloneTimeStepInfo(IndexType SolutionStepIndex, ProcessInfo const &  SourceSolutionStepInfo)
{
    CloneSolutionStepInfo(SolutionStepIndex, SourceSolutionStepInfo);
    mIsTimeStep = true;
}

void ProcessInfo::CreateTimeStepInfo(double NewTime, IndexType SolutionStepIndex)
{
    CreateTimeStepInfo(SolutionStepIndex);
    SetCurrentTime(NewTime);
}

void ProcessInfo::CloneTimeStepInfo(double NewTime, IndexType SourceSolutionStepIndex)
{
    CloneTimeStepInfo(SourceSolutionStepIndex);
    SetCurrentTime(NewTime);
}

void ProcessInfo::CloneTimeStepInfo(double NewTime, IndexType SolutionStepIndex, ProcessInfo const &  SourceSolutionStepInfo)
{
    CloneTimeStepInfo(SolutionStepIndex, SourceSolutionStepInfo);
    SetCurrentTime(NewTime);
}

void ProcessInfo::SetAsTimeStepInfo()
{
    mIsTimeStep = true;
    SetCurrentTime((*this)(TIME));
}

void ProcessInfo::SetAsTimeStepInfo(double NewTime)
{
    mIsTimeStep = true;
    SetCurrentTime(NewTime);
}

ProcessInfo::Pointer ProcessInfo::pGetPreviousTimeStepInfo(IndexType StepsBefore)
{
    if(StepsBefore > 1)
        return mpPreviousTimeStepInfo->pGetPreviousTimeStepInfo(--StepsBefore);

    if(StepsBefore == 0)
        KRATOS_THROW_ERROR(std::invalid_argument, "Steps before = 0", "");

    if(!mpPreviousTimeStepInfo)
        KRATOS_THROW_ERROR(std::invalid_argument, "No previous time step exist.", "");

    return mpPreviousTimeStepInfo;
}

const ProcessInfo::Pointer ProcessInfo::pGetPreviousTimeStepInfo(IndexType StepsBefore) const
{
    if(StepsBefore > 1)
        return mpPreviousTimeStepInfo->pGetPreviousTimeStepInfo(--StepsBefore);

    if(StepsBefore == 0)
        KRATOS_THROW_ERROR(std::invalid_argument, "Steps before = 0", "");

    if(!mpPreviousTimeStepInfo)
        KRATOS_THROW_ERROR(std::invalid_argument, "No previous time step exist.", "");

    return mpPreviousTimeStepInfo;
}

void ProcessInfo::SetCurrentTime(double NewTime)
{
    (*this)(TIME) = NewTime;
    if(!mpPreviousTimeStepInfo)
        (*this)(DELTA_TIME) = NewTime;
    else
        (*this)(DELTA_TIME) = NewTime -  mpPreviousTimeStepInfo->GetValue(TIME);
}

void ProcessInfo::ReIndexBuffer(SizeType BufferSize, IndexType BaseIndex)
{
    mSolutionStepIndex = BaseIndex;
    if(BufferSize > 1)
        if(mpPreviousSolutionStepInfo)
            mpPreviousSolutionStepInfo->ReIndexBuffer(BufferSize - 1, BaseIndex + 1);
}

void ProcessInfo::CreateSolutionStepInfo(IndexType SolutionStepIndex)
{
    mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
    mSolutionStepIndex = SolutionStepIndex;
    if(mIsTimeStep)
        mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
    mIsTimeStep = false;
    BaseType::Clear();
}

void ProcessInfo::CloneSolutionStepInfo()
{
    mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
    mSolutionStepIndex = 0;
    if(mIsTimeStep)
        mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
    mIsTimeStep = false;
}

void ProcessInfo::CloneSolutionStepInfo(IndexType SourceSolutionStepIndex)
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

void ProcessInfo::CloneSolutionStepInfo(IndexType SolutionStepIndex, ProcessInfo const &  SourceSolutionStepInfo)
{
    mpPreviousSolutionStepInfo = Pointer(new ProcessInfo(*this));
    mSolutionStepIndex = SolutionStepIndex;
    BaseType::operator=(SourceSolutionStepInfo);
    if(mIsTimeStep)
        mpPreviousTimeStepInfo = mpPreviousSolutionStepInfo;
    mIsTimeStep = false;
}

ProcessInfo& ProcessInfo::FindSolutionStepInfo(IndexType ThisIndex)
{
    if(mSolutionStepIndex == ThisIndex)
        return *this;

    if(!mpPreviousSolutionStepInfo)
        return *this;

    return mpPreviousTimeStepInfo->FindSolutionStepInfo(ThisIndex);
}

void ProcessInfo::RemoveSolutionStepInfo(IndexType SolutionStepIndex)
{
    if(!mpPreviousSolutionStepInfo)
        return;

    if(mpPreviousSolutionStepInfo->GetSolutionStepIndex() == SolutionStepIndex)
        mpPreviousSolutionStepInfo = mpPreviousSolutionStepInfo->pGetPreviousSolutionStepInfo();
    else
        mpPreviousSolutionStepInfo->RemoveSolutionStepInfo(SolutionStepIndex);
}

void ProcessInfo::ClearHistory(IndexType StepsBefore)
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

ProcessInfo::Pointer ProcessInfo::pGetPreviousSolutionStepInfo(IndexType StepsBefore)
{
    if(StepsBefore > 1)
        return mpPreviousSolutionStepInfo->pGetPreviousSolutionStepInfo(--StepsBefore);

    if(StepsBefore == 0)
        KRATOS_THROW_ERROR(std::invalid_argument, "Steps before = 0", "");

    if(!mpPreviousSolutionStepInfo)
        KRATOS_THROW_ERROR(std::invalid_argument, "No previous time step exist.", "");

    return mpPreviousSolutionStepInfo;
}

const ProcessInfo::Pointer ProcessInfo::pGetPreviousSolutionStepInfo(IndexType StepsBefore) const
{
    if(StepsBefore > 1)
        return mpPreviousSolutionStepInfo->pGetPreviousSolutionStepInfo(--StepsBefore);

    if(StepsBefore == 0)
        KRATOS_THROW_ERROR(std::invalid_argument, "Steps before = 0", "");

    if(!mpPreviousSolutionStepInfo)
        KRATOS_THROW_ERROR(std::invalid_argument, "No previous time step exist.", "");

    return mpPreviousSolutionStepInfo;
}

void ProcessInfo::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    rSerializer.save("Is Time Step",mIsTimeStep);
    rSerializer.save("Solution Step Index",mSolutionStepIndex);
    rSerializer.save("Previous Solution Step Info",mpPreviousSolutionStepInfo);
    rSerializer.save("Previous Time Step Info", mpPreviousTimeStepInfo);
}

void ProcessInfo::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    rSerializer.load("Is Time Step",mIsTimeStep);
    rSerializer.load("Solution Step Index",mSolutionStepIndex);
    rSerializer.load("Previous Solution Step Info",mpPreviousSolutionStepInfo);
    rSerializer.load("Previous Time Step Info", mpPreviousTimeStepInfo);
}

}  // namespace Kratos.


