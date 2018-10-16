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
#include "containers/pointer_vector.h"
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

/// Main class for online statistics calculation.
/** This class manages the definition, update and output of statistics calculated during a simulation.
 *  It is intended to be stored as a ProcessInfo variable and called at specific points during the
 *  simulation procedure.
 *  @see IntegrationPointStatisticsProcess for an implementation using this to compute statistics at the mesh integration points.
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
    ///@name Operations
    ///@{

    /// Add one first order statistic to be tracked during the simulation.
    /** @param[in] pResult pointer to the statistic quantity.
     */
    void AddResult(StatisticsSampler::Pointer pResult);

    /// Add one second or higher order statistic to be tracked during the simulation.
    /** @param[in] pResult pointer to the statistic quantity.
     */
    void AddHigherOrderStatistic(StatisticsSampler::Pointer pResult);

    /// Initialize elemental storage for the tracked statistics.
    /** Note that all results should be added before calling this.
     *  @param rElements List of elements whose integration points will be used to record statistics.
     */
    void InitializeStorage(ModelPart::ElementsContainerType& rElements);

    /// Record a new sample point.
    /** This function expects InitializeStorage to have been called in advance.
     *  Normal usage will be to call this once for time step (or every time a new measurement is available).
     *  It is also assumed that the elements will call the UpdateStatistics method of this class internally.
     *  @see FluidElement for an example of use.
     *  @param rModelPart ModelPart containing the elements where statistics are recorded.
     */
    void SampleIntegrationPointResults(ModelPart& rModelPart);

    /// Update statistics for a single element.
    /** This function should be called by elements supporting statistics calculation on integration points.
     *  @see FluidElement for an example of use.
     *  @param rElement The element.
     */
    void UpdateStatistics(Element* pElement);

    /// Output recorded statistics to a comma-separated value file.
    /** @param[in] rModelPart The model part instance where statistics are recorded.
     *  @param[in] rOutputFileName Base name for the output file.
     *  The .csv extension (and in MPI the rank) will be automatically appended to the provided file name.
     */
    void PrintToFile(const ModelPart &rModelPart, const std::string& rOutputFileName) const;

    /// Output current values of all recorded quantities as a vector of doubles.
    /** Provided for debug and testing purposes ONLY. */
    std::vector<double> OutputForTest(ModelPart::ElementsContainerType& rElements) const;

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

private:
    ///@name Static Member Variables
    ///@{

    static std::vector<double> mUpdateBuffer;
    #pragma omp threadprivate(mUpdateBuffer)

    ///@}
    ///@name Member Variables
    ///@{

    bool mInitialized;

    std::size_t mDataBufferSize;

    std::size_t mRecordedSteps;

    PointerVector<StatisticsSampler> mAverageData;

    PointerVector<StatisticsSampler> mHigherOrderData;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const {}

    void load(Serializer& rSerializer) {}

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
