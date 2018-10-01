//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Author2 Fullname
//

#ifndef KRATOS_STATISTICAL_UTILITIES_H_INCLUDED
#define KRATOS_STATISTICAL_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <functional>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

class StatisticsSampler
{
public:

StatisticsSampler(unsigned int NumValues, unsigned int Offset):
    mNumValues(NumValues),
    mOffset(Offset)
{}

virtual ~StatisticsSampler() {}

virtual void SampleDataPoint(std::vector<double>::iterator& BufferIterator)
{}

unsigned int GetSize() const {
    return mNumValues;
}

unsigned int GetValueOffset() const
{
    KRATOS_ERROR << "Method not implemented" << std::endl;
}

unsigned int GetComponentOffset(unsigned int i) const
{
    KRATOS_ERROR << "Method not implemented" << std::endl;
}

unsigned int GetComponentOffset(unsigned int i, unsigned int j) const
{
    KRATOS_ERROR << "Method not implemented" << std::endl;
}

protected:

unsigned int GetOffset() const
{
    return mOffset;
}

private:

unsigned int mNumValues;

unsigned int mOffset;

};

class ScalarAverageSampler: public StatisticsSampler
{
public:

ScalarSampler(unsigned int Offset, std::function<double> Getter):
  StatisticsSampler(1, Offset),
  mGetter(Getter)
{}

~ScalarSampler() override {}

void SampleDataPoint(std::vector<double>::iterator& BufferIterator, const std::vector<double>& rData)
{
    *BufferIterator = mGetter();
    ++BufferIterator;
}

unsigned int GetValueOffset() const override {
    return this->GetOffset();
}

private:

std::function<double> mGetter;

};

template< class VectorType >
class VectorAverageSampler: public StatisticsSampler
{
public:

VectorSampler(std::function<VectorType> Getter, unsigned int VectorSize):
    StatisticsSampler(VectorSize),
    mGetter(Getter)
{}

~VectorSampler() override {}

void SampleDataPoint(std::vector<double>::iterator& BufferIterator) override {
    TVectorType result = mGetter();
    for (unsigned int i = 0; i < this->GetSize(); i++) {
        *BufferIterator = result[i];
        ++BufferIterator;
    }
}

private:

std::function<VectorType> mGetter;
};

/// Short class definition.
/** Detail class definition.
  */
template< class ValueContainterType >
class StatisticsUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StatisticsUtilities
    KRATOS_CLASS_POINTER_DEFINITION(StatisticsUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StatisticsUtilities() {}

    /// Destructor.
    virtual ~StatisticsUtilities() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void AddQuantity(StatisticsSampler::Pointer pAverage) {
        KRATOS_ERROR_IF(mInitialized) << "Trying to add more statistics quantities after initialization" << std::endl;
        mAverages.push_back(pAverage);
    }

    void AddQuantity(HigherOrderStatistic::Pointer pHigherOrderStatistic)  {
        KRATOS_ERROR_IF(mInitialized) << "Trying to add more statistics quantities after initialization" << std::endl;
        mHigherOrderStatistics.push_back(pHigherOrderStatistic);
    }

    void Initialize() {}

    void Update(const ValueContainerType& rMeasurement) {}

    void Combine(const StatisticsUtilities& rOther) {}

    ValueContainerType Output() {}

    void LoadFromFile() {}

    void SaveToFile() {}

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
        std::stringstream buffer;
        buffer << "StatisticsUtilities";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "StatisticsUtilities"; }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const {}

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

    std::vector<StatisticsSampler::Pointer> mAverages;

    std::vector<HigherOrderStatistic::Pointer> mHigherOrderStatistics;

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
    StatisticsUtilities &operator=(StatisticsUtilities const &rOther) {}

    /// Copy constructor.
    StatisticsUtilities(StatisticsUtilities const &rOther) {}

    ///@}

}; // Class StatisticsUtilities

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                StatisticsUtilities &rThis) {}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const StatisticsUtilities &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_STATISTICAL_UTILITIES_H_INCLUDED  defined
