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

#ifndef KRATOS_STATISTICS_UTILITIES_H_INCLUDED
#define KRATOS_STATISTICS_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <functional>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry.h"

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

KRATOS_CLASS_POINTER_DEFINITION(StatisticsSampler);

typedef Matrix::iterator1 IntegrationPointDataView;
typedef Matrix::const_iterator2 IntegrationPointDataViewIterator;

StatisticsSampler(unsigned int NumValues, unsigned int Offset = 0):
    mNumValues(NumValues),
    mOffset(Offset)
{}

virtual ~StatisticsSampler() {}

// For first-order statistics: read data directly
virtual void SampleDataPoint(
    const Geometry< Node<3> >& rGeometry,
    const Vector& rShapeFunctions,
    const Matrix& rShapeDerivatives,
    std::vector<double>::iterator& BufferIterator)
{}

// For higher-order statistics: operate on first order data
virtual void SampleDataPoint(
    std::vector<double>::iterator& BufferIterator,
    const StatisticsSampler::IntegrationPointDataView& rCurrentStatistics,
    const std::vector<double>& rNewMeasurement,
    const std::size_t NumberOfMeasurements)
{}

unsigned int GetSize() const {
    return mNumValues;
}

virtual unsigned int GetValueOffset() const
{
    KRATOS_ERROR << "Method not implemented" << std::endl;
}

virtual unsigned int GetComponentOffset(std::size_t i) const
{
    KRATOS_DEBUG_ERROR_IF(i >= mNumValues)
    << "Trying to access component index " << i << ", but only "
    << mNumValues << " components are stored." << std::endl;

    return mOffset + i;
}

virtual std::size_t GetComponentOffset(std::size_t i, std::size_t j) const
{
    KRATOS_ERROR << "Method not implemented" << std::endl;
}

virtual std::size_t ComponentIndex(std::size_t i, std::size_t j) const
{
    return i*mNumValues + j;
}

std::size_t GetOffset() const
{
    return mOffset;
}

void SetOffset(std::size_t Offset) {
    mOffset = Offset;
}

virtual void OutputResult(
    std::ofstream& rOutStream,
    IntegrationPointDataViewIterator& rDataBuffer,
    std::size_t SampleSize,
    std::string& rSeparator) const
{
    for (std::size_t i = 0; i < mNumValues; i++) {
        rOutStream << rSeparator << Finalize(*rDataBuffer, SampleSize);
        ++rDataBuffer;
    }
}

virtual double Finalize(double Value, std::size_t SampleSize) const
{
    return Value / SampleSize;
}

private:

const std::size_t mNumValues;

std::size_t mOffset;

friend class Serializer;

void save(Serializer& rSerializer) const {}

void load(Serializer& rSerializer) {}

StatisticsSampler():
    mNumValues(0),
    mOffset(0)
{}

};

class ScalarAverageSampler: public StatisticsSampler
{
public:

ScalarAverageSampler(std::function<double(const Geometry< Node<3> >&, const Vector&, const Matrix&)> Getter):
  StatisticsSampler(1),
  mGetter(Getter)
{}

~ScalarAverageSampler() override {}

void SampleDataPoint(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives, std::vector<double>::iterator& BufferIterator) override
{
    *BufferIterator = mGetter(rGeometry,rShapeFunctions,rShapeDerivatives);
    ++BufferIterator;
}

unsigned int GetValueOffset() const override {
    return this->GetOffset();
}

private:

std::function<double(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives)> mGetter;

};

template< class VectorType >
class VectorAverageSampler: public StatisticsSampler
{
public:

VectorAverageSampler(std::function<VectorType(const Geometry< Node<3> >&, const Vector&, const Matrix&)> Getter, unsigned int VectorSize):
    StatisticsSampler(VectorSize),
    mGetter(Getter)
{}

~VectorAverageSampler() override {}

void SampleDataPoint(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives, std::vector<double>::iterator& BufferIterator) override {
    VectorType result = mGetter(rGeometry,rShapeFunctions,rShapeDerivatives);
    for (unsigned int i = 0; i < this->GetSize(); i++) {
        *BufferIterator = result[i];
        ++BufferIterator;
    }
}

private:

std::function<VectorType(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives)> mGetter;
};

class VarianceSampler : public StatisticsSampler
{
public:

VarianceSampler(const StatisticsSampler::Pointer pQuantity1, const StatisticsSampler::Pointer pQuantity2):
    StatisticsSampler(pQuantity1->GetSize() * pQuantity2->GetSize()),
    mpQuantity1(pQuantity1),
    mpQuantity2(pQuantity2)
{}

void SampleDataPoint(
    std::vector<double>::iterator& BufferIterator,
    const StatisticsSampler::IntegrationPointDataView& rCurrentStatistics,
    const std::vector<double>& rNewMeasurement,
    const std::size_t NumberOfMeasurements) override
{
    const double update_factor = 1.0 / ((NumberOfMeasurements-1)*NumberOfMeasurements);
    for (std::size_t i = 0; i < mpQuantity1->GetSize(); i++)
    {
        double current_total_i = *(rCurrentStatistics.begin() + mpQuantity1->GetComponentOffset(i));
        double new_measurement_i = rNewMeasurement[mpQuantity1->GetComponentOffset(i)];
        double delta_i = (NumberOfMeasurements-1)*current_total_i - new_measurement_i;
        for (std::size_t j = 0; j < mpQuantity2->GetSize(); j++)
        {
            double current_total_j = *(rCurrentStatistics.begin() + mpQuantity2->GetComponentOffset(j));
            double new_measurement_j = rNewMeasurement[mpQuantity2->GetComponentOffset(j)];
            double delta_j = (NumberOfMeasurements-1)*current_total_j - new_measurement_j;
            (*BufferIterator) = update_factor * delta_i * delta_j;
            ++BufferIterator;
        }
    }

}

// Override the Finalize method to implement the unbiased variance (divide by n-1)
double Finalize(double Value, std::size_t SampleSize) const override
{
    return Value / (SampleSize - 1);
}

protected:

VarianceSampler(const StatisticsSampler::Pointer pQuantity1, const StatisticsSampler::Pointer pQuantity2, std::size_t DataSize):
    StatisticsSampler(DataSize),
    mpQuantity1(pQuantity1),
    mpQuantity2(pQuantity2)
{}

StatisticsSampler::Pointer GetQuantity1()
{
    return mpQuantity1;
}

private:

const StatisticsSampler::Pointer mpQuantity1;

const StatisticsSampler::Pointer mpQuantity2;

};

class SymmetricVarianceSampler: public VarianceSampler
{
public:

SymmetricVarianceSampler(const StatisticsSampler::Pointer pQuantity1):
    VarianceSampler(pQuantity1, pQuantity1, ((pQuantity1->GetSize()+1) * pQuantity1->GetSize()) / 2)
{}

std::size_t ComponentIndex(std::size_t i, std::size_t j) const override
{
    if (i <= j)
    {
        return i*this->GetSize() - (i*(i+1))/2 + j;
    }
    else
    {
        return j*this->GetSize() - (j*(j+1))/2 + i;
    }
}

void SampleDataPoint(
    std::vector<double>::iterator& BufferIterator,
    const StatisticsSampler::IntegrationPointDataView& rCurrentStatistics,
    const std::vector<double>& rNewMeasurement,
    const std::size_t NumberOfMeasurements) override
{
    const double update_factor = 1.0 / ((NumberOfMeasurements-1)*NumberOfMeasurements);
    for (std::size_t i = 0; i < GetQuantity1()->GetSize(); i++)
    {
        double current_total_i = *(rCurrentStatistics.begin() + GetQuantity1()->GetComponentOffset(i));
        double new_measurement_i = rNewMeasurement[GetQuantity1()->GetComponentOffset(i)];
        double delta_i = (NumberOfMeasurements-1)*current_total_i - new_measurement_i;
        for (std::size_t j = i; j < GetQuantity1()->GetSize(); j++)
        {
            double current_total_j = *(rCurrentStatistics.begin() + GetQuantity1()->GetComponentOffset(j));
            double new_measurement_j = rNewMeasurement[GetQuantity1()->GetComponentOffset(j)];
            double delta_j = (NumberOfMeasurements-1)*current_total_j - new_measurement_j;
            (*BufferIterator) = update_factor * delta_i * delta_j;
            ++BufferIterator;
        }
    }
}

};

class ComponentwiseVarianceSampler: public VarianceSampler
{
public:

ComponentwiseVarianceSampler(
    const StatisticsSampler::Pointer pQuantity1,
    std::size_t ComponentIndex1,
    const StatisticsSampler::Pointer pQuantity2,
    std::size_t ComponentIndex2)
    : VarianceSampler(pQuantity1,pQuantity2,1)
    , mComponent1(ComponentIndex1)
    , mComponent2(ComponentIndex2)
{}

void SampleDataPoint(
    std::vector<double>::iterator& BufferIterator,
    const StatisticsSampler::IntegrationPointDataView& rCurrentStatistics,
    const std::vector<double>& rNewMeasurement,
    const std::size_t NumberOfMeasurements) override
{
    const double update_factor = 1.0 / ((NumberOfMeasurements-1)*NumberOfMeasurements);

    double current_total_i = *(rCurrentStatistics.begin() + GetQuantity1()->GetComponentOffset(mComponent1));
    double new_measurement_i = rNewMeasurement[GetQuantity1()->GetComponentOffset(mComponent1)];
    double delta_i = (NumberOfMeasurements-1)*current_total_i - new_measurement_i;

    double current_total_j = *(rCurrentStatistics.begin() + GetQuantity1()->GetComponentOffset(mComponent2));
    double new_measurement_j = rNewMeasurement[GetQuantity1()->GetComponentOffset(mComponent2)];
    double delta_j = (NumberOfMeasurements-1)*current_total_j - new_measurement_j;
    (*BufferIterator) = update_factor * delta_i * delta_j;
    ++BufferIterator;
}

private:

std::size_t mComponent1;

std::size_t mComponent2;

};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                StatisticsSampler &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const StatisticsSampler &rThis)
{
    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_STATISTICS_UTILITIES_H_INCLUDED  defined
