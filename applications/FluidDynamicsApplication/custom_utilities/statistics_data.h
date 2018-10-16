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
#include "includes/element.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Internal container for integration point statistitcs on a given element.
/** @see StatisticsRecord @see IntegrationPointStatisticsProcess.
  */
class StatisticsData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StatisticsData
    KRATOS_CLASS_POINTER_DEFINITION(StatisticsData);

    typedef std::vector<double> ValueContainerType;

    typedef Matrix::iterator1 IntegrationPointDataView;
    typedef Matrix::const_iterator1 IntegrationPointDataConstView;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StatisticsData() {}

    /// Copy constructor.
    StatisticsData(StatisticsData const &rOther) {
        mData = rOther.mData;
    }

    /// Destructor.
    virtual ~StatisticsData() {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    StatisticsData &operator=(StatisticsData const &rOther) {
        mData = rOther.mData;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Initialize internal storage for this container.
    /** @param rElement The element instance this container stores data for.
     *  @param MeasurementSize Total number of quantities that are recorded.
     *  For vector and matrix statistics, each component should be counted as a different quantity.
     */
    void InitializeStorage(Element& rElement, std::size_t MeasurementSize)
    {
        const GeometryData::IntegrationMethod integration_method = rElement.GetIntegrationMethod();
        std::size_t number_of_integration_points = rElement.GetGeometry().IntegrationPointsNumber(integration_method);

        mData.resize(number_of_integration_points, MeasurementSize, false);
        mData = ZeroMatrix(number_of_integration_points, MeasurementSize);
    }

    /// Record a new realization for the measured statistics.
    /** This function assumes that InitializeStorage has been called in advance.
     *  @param rElement The element instance this container stores data for.
     *  @param rRecordedStatistics List of first order statistics stored in this container.
     *  @param rHigherOrderStatistics List of higher order statistics stored in this container.
     *  @param NumberOfMeasurements How many times this function has been called so far (including this one).
     */
    void UpdateMeasurement(
        const Element& rElement,
        const PointerVector<StatisticsSampler>& rStatisticsSamplers,
        const PointerVector<StatisticsSampler>& rHigherOrderStatistics,
        ValueContainerType& rUpdate,
        std::size_t NumMeasurements)
    {
        KRATOS_DEBUG_ERROR_IF(NumMeasurements == 0)
        << "Trying to update statistics, but providied number of recorded steps is zero" << std::endl;

        const Geometry<Node<3>> &r_geometry = rElement.GetGeometry();
        const GeometryData::IntegrationMethod integration_method = rElement.GetIntegrationMethod();
        Matrix shape_functions;
        typename Geometry<Node<3>>::ShapeFunctionsGradientsType shape_gradients;
        this->CalculateGeometryData(r_geometry, integration_method, shape_functions, shape_gradients);

        for (unsigned int g = 0; g < shape_functions.size1(); g++)
        {
            auto N = row(shape_functions, g);
            auto &rDN_DN = shape_gradients[g];

            auto it_update_buffer = rUpdate.begin();
            for (auto it_sampler = rStatisticsSamplers.begin(); it_sampler != rStatisticsSamplers.end(); ++it_sampler)
            {
                it_sampler->SampleDataPoint(r_geometry, N, rDN_DN, it_update_buffer);
            }

            if (NumMeasurements > 1) { // Second order and higher statistics start from the second iteration
                // loop on higer order statistics. They require the already written rMeasurements, (const) mData and the number of steps as input
                for (auto it_statistic = rHigherOrderStatistics.begin(); it_statistic != rHigherOrderStatistics.end(); ++it_statistic)
                {
                    it_statistic->SampleDataPoint(it_update_buffer, IntegrationPointData(g), rUpdate, NumMeasurements);
                }
            }

            for (unsigned int i = 0; i < rUpdate.size(); i++)
            {
                mData(g,i) += rUpdate[i];
            }
        }
    }

    /// Write computed statistics to output file.
    /** @param rOutputStream output file stream.
     *  @param rElement Element instance this container stores data for.
     *  @param rRecordedStatistics List of first order statistics stored in this container.
     *  @param rHigherOrderStatistics List of higher order statistics stored in this container.
     *  @param NumberOfMeasurements How many realizations have been recorded up to this point.
     *  @param rSeparator Separator to use in the output file.
     */
    void WriteToCSVOutput(
        std::ofstream& rOutputStream,
        const Element& rElement,
        const PointerVector<StatisticsSampler>& rRecordedStatistics,
        const PointerVector<StatisticsSampler>& rHigherOrderStatistics,
        std::size_t NumberOfMeasurements,
        const std::string& rSeparator) const
    {
        const Geometry<Node<3>> &r_geometry = rElement.GetGeometry();
        const GeometryData::IntegrationMethod integration_method = rElement.GetIntegrationMethod();
        Matrix shape_functions;
        typename Geometry<Node<3>>::ShapeFunctionsGradientsType shape_gradients;
        this->CalculateGeometryData(r_geometry, integration_method, shape_functions, shape_gradients);

        for (unsigned int g = 0; g < shape_functions.size1(); g++)
        {
            // Print element ID and integration point index
            rOutputStream << rElement.Id() << rSeparator << g << rSeparator;

            // Print integration point coordinates
            array_1d<double,3> coordinates(3,0.0);
            for (unsigned int n = 0; n < shape_functions.size2(); n++)
                coordinates += shape_functions(g,n) * r_geometry[n].Coordinates();
            rOutputStream << coordinates[0] << rSeparator << coordinates[1] << rSeparator << coordinates[2];

            auto data_iterator = IntegrationPointData(g).begin();
            for (auto it_sampler = rRecordedStatistics.begin(); it_sampler != rRecordedStatistics.end(); ++it_sampler)
            {
                it_sampler->OutputResult(rOutputStream,data_iterator,NumberOfMeasurements,rSeparator);
            }

            for (auto it_sampler = rHigherOrderStatistics.begin(); it_sampler != rHigherOrderStatistics.end(); ++it_sampler)
            {
                it_sampler->OutputResult(rOutputStream,data_iterator,NumberOfMeasurements,rSeparator);
            }

            rOutputStream << "\n";
        }
    }

    ///@}
    ///@name Inquiry
    ///@{

    /// How many integration points this container holds data for.
    std::size_t NumberOfIntegrationPoints() const
    {
        return mData.size1();
    }

    /// How many quantities this container holds data for.
    /** When working with vector or matrix quantities, each individual component
     *  is counted as a different quantity.
     */
    std::size_t NumberOfStatisticalQuantities() const
    {
        return mData.size2();
    }

    ///@}
    ///@name Access
    ///@{

    /// Access internal data for a given integration point.
    /** @param IntegrationPointIndex index of the requested integration point.
     */
    IntegrationPointDataView IntegrationPointData(std::size_t IntegrationPointIndex)
    {
        KRATOS_DEBUG_ERROR_IF(IntegrationPointIndex >= mData.size1())
            << "Asking for integration point number " << IntegrationPointIndex
            << " but only " << mData.size1() << " points are recorded." << std::endl;
        return (mData.begin1() + IntegrationPointIndex);
    }

    /// Access internal data for a given integration point.
    /** @param IntegrationPointIndex index of the requested integration point.
     */
    IntegrationPointDataConstView IntegrationPointData(std::size_t IntegrationPointIndex) const
    {
        KRATOS_DEBUG_ERROR_IF(IntegrationPointIndex >= mData.size1())
            << "Asking for integration point number " << IntegrationPointIndex
            << " but only " << mData.size1() << " points are recorded." << std::endl;
        return (mData.begin1() + IntegrationPointIndex);
    }

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

protected:

    ///@name Protected Operations
    ///@{

    void CalculateGeometryData(
        const Geometry<Node<3>>& rGeometry,
        const GeometryData::IntegrationMethod IntegrationMethod,
        Matrix &rN,
        typename Geometry<Node<3>>::ShapeFunctionsGradientsType &rDN_DX) const
    {
        const unsigned int number_of_nodes = rGeometry.PointsNumber();
        const unsigned int number_of_gauss_points = rGeometry.IntegrationPointsNumber(IntegrationMethod);

        Vector det_J;
        rGeometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX,det_J,IntegrationMethod);

        rN.resize(number_of_gauss_points,number_of_nodes,false);
        rN = rGeometry.ShapeFunctionsValues(IntegrationMethod);
    }

    ///@}

private:

    ///@name Member Variables
    ///@{

    Matrix mData;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const {}

    void load(Serializer& rSerializer) {}

    ///@}

}; // Class StatisticsData

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                StatisticsData/*< std::vector<double> >*/ &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const StatisticsData/*< std::vector<double> >*/ &rThis)
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
