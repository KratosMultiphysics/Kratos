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
template< class ValueContainerType >
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

    void Initialize(std::size_t BufferSize)
    {
        mData.clear();
        mData.resize(BufferSize);
        for (auto iter = mData.begin(); iter != mData.end(); ++iter)
        {
            *iter = 0.0;
        }
    }

    void CalculateUpdateDelta(
        const Element* pElement,
        const std::vector<StatisticsSampler::Pointer>& rStatisticsSamplers,
        ValueContainerType& rMeasurements,
        std::size_t NumMeasurements)
    {
        const Geometry< Node<3> >& r_geometry = pElement->GetGeometry();
        const GeometryData::IntegrationMethod integration_method = pElement->GetIntegrationMethod();
        Matrix shape_functions;
        typename Geometry<Node<3>>::ShapeFunctionsGradientsType shape_gradients;
        this->CalculateGeometryData(r_geometry,integration_method,shape_functions,shape_gradients);

        auto it_measurement_buffer = rMeasurements.begin();
        for (auto it_sampler = rStatisticsSamplers.begin(); it_sampler != rStatisticsSamplers.end(); ++it_sampler) {
            (**it_sampler).SampleDataPoint(r_geometry,shape_functions,shape_gradients,it_measurement_buffer);
        }
    }

    void AddMeasurement(const ValueContainerType& rMeasurement, std::size_t NumSteps) {}

    void Combine(const StatisticsData& rOther) {}

    ValueContainerType Output() {
        return mData;
    }

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

    ValueContainerType mData;

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
                                StatisticsData< std::vector<double> > &rThis)
{
    return rIStream;
}

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
