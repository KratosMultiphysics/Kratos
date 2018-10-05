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
//

#ifndef KRATOS_STATISTICS_RECORD_H_INCLUDED
#define KRATOS_STATISTICS_RECORD_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "statistics_utilities.h"

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

/// Short class definition.
/** Detail class definition.
  */
class StatisticsRecord
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StatisticsRecord
    KRATOS_CLASS_POINTER_DEFINITION(StatisticsRecord);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StatisticsRecord()
    : mInitialized(false)
    , mDataBufferSize(0)
    , mRecordedSteps(0)
    {}

    /// Destructor.
    virtual ~StatisticsRecord()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void AddResult(StatisticsSampler::Pointer pResult);

    void InitializeStorage();

    void SampleIntegrationPointResults(ModelPart& rModelPart);

    void UpdateStatistics(Element* pElement);

    void FinalizeStatistics(ModelPart::ElementsContainerType& rElements);

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
        buffer << "StatisticsRecord";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "StatisticsRecord"; }

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

    static std::vector<double> mUpdateBuffer;
    #pragma omp threadprivate(mUpdateBuffer)

    static std::vector<double> mMeasurementBuffer;
    #pragma omp threadprivate(mMeasurementBuffer)

    ///@}
    ///@name Member Variables
    ///@{

    bool mInitialized;

    std::size_t mDataBufferSize;

    std::size_t mRecordedSteps;

    std::vector<StatisticsSampler::Pointer> mAverageData;

    std::vector<StatisticsSampler::Pointer> mHigherOrderData;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const {}

    void load(Serializer& rSerializer) {}

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
    StatisticsRecord &operator=(StatisticsRecord const &rOther) = delete;

    /// Copy constructor.
    StatisticsRecord(StatisticsRecord const &rOther) {}

    ///@}

}; // Class StatisticsRecord

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                StatisticsRecord &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const StatisticsRecord &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_STATISTICS_RECORD_H_INCLUDED  defined
