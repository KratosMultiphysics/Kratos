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

#ifndef KRATOS_STATISTICS_DATA_H_INCLUDED
#define KRATOS_STATISTICS_DATA_H_INCLUDED

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

/// Short class definition.
/** Detail class definition.
  */
template< class ValueContainterType >
class StatisticsData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StatisticsData
    KRATOS_CLASS_POINTER_DEFINITION(StatisticsData);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StatisticsData() {}

    /// Destructor.
    virtual ~StatisticsData() {}

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

    void Combine(const StatisticsData& rOther) {}

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
        buffer << "StatisticsData";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "StatisticsData"; }

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
    StatisticsData &operator=(StatisticsData const &rOther) {}

    /// Copy constructor.
    StatisticsData(StatisticsData const &rOther) {}

    ///@}

}; // Class StatisticsData

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                StatisticsData< std::vector<double> > &rThis) {}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const StatisticsData< std::vector<double> > &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_STATISTICS_DATA_H_INCLUDED  defined
