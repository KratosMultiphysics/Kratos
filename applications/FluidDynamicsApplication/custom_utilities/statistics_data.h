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
class StatisticsData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StatisticsData
    KRATOS_CLASS_POINTER_DEFINITION(StatisticsData);

    typedef std::vector<double> ValueContainerType;

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

    void UpdateMeasurement(
        /*const*/ Element* pElement,
        const std::vector<StatisticsSampler::Pointer>& rStatisticsSamplers,
        ValueContainerType& rMeasurement,
        ValueContainerType& rUpdate,
        std::size_t NumMeasurements)
    {
        if (NumMeasurements <= 1)
        {
            InitializeStorage(pElement, rMeasurement.size());
        }
        else
        {
            const Geometry<Node<3>> &r_geometry = pElement->GetGeometry();
            const GeometryData::IntegrationMethod integration_method = pElement->GetIntegrationMethod();
            Matrix shape_functions;
            typename Geometry<Node<3>>::ShapeFunctionsGradientsType shape_gradients;
            this->CalculateGeometryData(r_geometry, integration_method, shape_functions, shape_gradients);

            for (unsigned int g = 0; g < shape_functions.size1(); g++)
            {
                auto N = row(shape_functions, g);
                auto &rDN_DN = shape_gradients[g];

                auto it_measurement_buffer = rMeasurement.begin();
                for (auto it_sampler = rStatisticsSamplers.begin(); it_sampler != rStatisticsSamplers.end(); ++it_sampler)
                {
                    (**it_sampler).SampleDataPoint(r_geometry, N, rDN_DN, it_measurement_buffer);
                }

                // loop on higer order statistics. They require the already written rMeasurements, (const) mData and the number of steps as input

                for (unsigned int i = 0; i < rUpdate.size(); i++)
                {
                    mData(g,i) += rUpdate[i];
                }
            }
        }
    }

/*    void CalculateUpdateDelta(
        const Element* pElement,
        const std::vector<StatisticsSampler::Pointer>& rStatisticsSamplers,
        ValueContainerType& rMeasurements,
        std::size_t NumMeasurements) const
    {
        const Geometry< Node<3> >& r_geometry = pElement->GetGeometry();
        const GeometryData::IntegrationMethod integration_method = pElement->GetIntegrationMethod();
        Matrix shape_functions;
        typename Geometry<Node<3>>::ShapeFunctionsGradientsType shape_gradients;
        this->CalculateGeometryData(r_geometry,integration_method,shape_functions,shape_gradients);

        auto it_measurement_buffer = rMeasurements.begin();
        for (auto it_sampler = rStatisticsSamplers.begin(); it_sampler != rStatisticsSamplers.end(); ++it_sampler) {
            //(**it_sampler).SampleDataPoint(r_geometry,shape_functions,shape_gradients,it_measurement_buffer);
        }

        // loop on higer order statistics. They require the already written rMeasurements, (const) mData and the number of steps as input
    }

    void UpdateMeasurement(const ValueContainerType& rUpdate)
    {
        for (unsigned int i = 0; i < mData.size(); i++)
        {
            mData[i] += rUpdate[i];
        }
    }*/

    void Combine(const StatisticsData& rOther) {}

    /*ValueContainerType Output() {
        return mData;
    }*/

    void LoadFromFile() {}

    void SaveToFile() {}

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    std::size_t NumberOfIntegrationPoints()
    {
        return mData.size1();
    }

    typename Matrix::iterator2 DataIterator(std::size_t IntegrationPointIndex)
    {
        KRATOS_DEBUG_ERROR_IF(IntegrationPointIndex >= mData.size1())
            << "Asking for integration point number " << IntegrationPointIndex
            << " but only " << mData.size1() << " points are recorded." << std::endl;
        return (mData.begin1() + IntegrationPointIndex).begin();
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

    Matrix mData;

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

    void InitializeStorage(Element* pElement, std::size_t MeasurementSize)
    {
        const GeometryData::IntegrationMethod integration_method = pElement->GetIntegrationMethod();
        std::size_t number_of_integration_points = pElement->GetGeometry().IntegrationPointsNumber(integration_method);

        mData.resize(number_of_integration_points, MeasurementSize, false);
        mData = ZeroMatrix(number_of_integration_points, MeasurementSize);
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

}; // Class StatisticsData

///@}

///@name Type Definitions
///@{

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
