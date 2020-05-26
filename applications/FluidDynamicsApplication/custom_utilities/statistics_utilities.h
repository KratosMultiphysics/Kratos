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

///@name Kratos Classes
///@{

/// Base class for statistical measurements.
class StatisticsSampler
{
public:

KRATOS_CLASS_POINTER_DEFINITION(StatisticsSampler);

typedef Matrix::iterator1 IntegrationPointDataView;
typedef Matrix::const_iterator2 IntegrationPointDataViewIterator;

/// Define a new StatisticsSampler instance.
/** Note that this is the base class, unless you are implementing a new statistic,
 *  you should work with derived classes directly.
 *  @param NumValues Number of different quantities (components) to be managed.
 */
StatisticsSampler(std::size_t NumValues):
    mNumValues(NumValues),
    mOffset(0)
{}

/// Destructor
virtual ~StatisticsSampler() {}

/// For first-order statistics: read data directly.
/** To be implemented in derived classes.
 *  Note that the BufferIterator argument will be advanced during the call.
 *  @param[in] rGeometry The geometry for the element that contains the point where statistics are sampled.
 *  @param[in] rShapeFunctions Local values of elemental shape functions at the sampling point.
 *  @param[in] rShapeDerivatives Local values of elemental shape function gradients at the sampling point.
 *  @param BufferIterator Iterator pointing to the start of the storage of the statistics managed by this class.
 */
virtual void SampleDataPoint(
    const Geometry< Node<3> >& rGeometry,
    const Vector& rShapeFunctions,
    const Matrix& rShapeDerivatives,
    std::vector<double>::iterator& BufferIterator)
{}

/// For higher-order statistics: operate on lower order data.
/** To be implemented in derived classes.
 *  Note that the BufferIterator argument will be advanced during the call.
 *  @param BufferIterator Iterator pointing to the start of the storage of the statistics managed by this class.
 *  @param[in] rCurrentStatistics Container of last known values for recorded quantities at this sampling point.
 *  @param[in] rNewMeasurement Container of new measurements for recorded quantities at this sampling point.
 *  @param[in] NumberOfMeasurements Total number of realizations measured at the sampling point (including the current one).
 */
virtual void SampleDataPoint(
    std::vector<double>::iterator& BufferIterator,
    const StatisticsSampler::IntegrationPointDataView& rCurrentStatistics,
    const std::vector<double>& rNewMeasurement,
    const std::size_t NumberOfMeasurements)
{}

/// Number of quantities managed by this statistic.
std::size_t GetSize() const {
    return mNumValues;
}

/// Offset (from the start of the space allocated to this statistic) for the storage of component i.
virtual std::size_t GetComponentOffset(std::size_t i) const
{
    KRATOS_DEBUG_ERROR_IF(i >= mNumValues)
    << "Trying to access component index " << i << ", but only "
    << mNumValues << " components are stored." << std::endl;

    return mOffset + i;
}

/// Helper returning the correct argument in calls to GetComponentOffset for matrix quantities.
virtual std::size_t ComponentIndex(std::size_t i, std::size_t j) const
{
    return i*mNumValues + j;
}

/// Offset (from the start of the statistics container) to the first component stored by this statistic.
std::size_t GetOffset() const
{
    return mOffset;
}

/// Assign a new offset (from the start of the statistics container) to the first component stored by this statistic.
void SetOffset(std::size_t Offset) {
    mOffset = Offset;
}

/// Write results managed by this class to output buffer.
/** This function is called when printing to a file.
 *  Note that the rDataBuffer argument will be advanced during the call.
 *  @param rOutStream stream for the output file.
 *  @param rDataBuffer Iterator pointing to the data managed by this statistic.
 *  @param SampleSize Number of measured realizations for the recorded quantity.
 *  @param rSeparator Separator for the output csv file (for example a comma or semicolon).
 */
virtual void OutputResult(
    std::ofstream& rOutStream,
    IntegrationPointDataViewIterator& rDataBuffer,
    std::size_t SampleSize,
    const std::string& rSeparator) const
{
    for (std::size_t i = 0; i < mNumValues; i++) {
        rOutStream << rSeparator << Finalize(*rDataBuffer, SampleSize);
        ++rDataBuffer;
    }
}

/// Write header for the output file.
/** This function is called when printing to a file.
 *  @param rOutStream stream for the output file.
 *  @param rSeparator Separator for the output csv file (for example a comma or semicolon).
 */
virtual void OutputHeader(
    std::ofstream& rOutStream,
    const std::string& rSeparator) const
{
    rOutStream << rSeparator;
}

/// Post-process internal data to produce the final value of the statistical result.
/** In most cases, this corresponds to dividing by sample size.
 */
virtual double Finalize(double Value, std::size_t SampleSize) const
{
    return Value / SampleSize;
}

/// Get the string associated one of the components of this statistic.
/** This is used to obtain the "name" for the statistic to be printed in the output file header.
 *  @param[in] ComponentIndex index of the required component.
 */
std::string GetTag(std::size_t ComponentIndex) const
{
    KRATOS_DEBUG_ERROR_IF(ComponentIndex >= mNumValues)
    << "Requesting tag for component " << ComponentIndex
    << " but only " << mNumValues << " components are stored." << std::endl;
    return mTags[ComponentIndex];
}

protected:

void AddTag(std::string Tag)
{
    mTags.push_back(Tag);
}

private:

const std::size_t mNumValues;

std::size_t mOffset;

std::vector<std::string> mTags;

friend class Serializer;

void save(Serializer& rSerializer) const {}

void load(Serializer& rSerializer) {}

StatisticsSampler():
    mNumValues(0),
    mOffset(0)
{}

};

/// This class manages the computation of the average of a scalar quantity.
class ScalarAverageSampler: public StatisticsSampler
{
public:

/// Initialize a sampler for a scalar average quantity.
/** @param Getter Helper function to evaluate the quantity of interest.
 *  @param Tag Name to identify the quantity on result files.
 *  @see Internals::MakeSamplerAtLocalCoordinate for the definition of the Getter argument.
 */
ScalarAverageSampler(
    std::function<double(const Geometry< Node<3> >&, const Vector&, const Matrix&)> Getter,
    const std::string& Tag)
    :
    StatisticsSampler(1),
    mGetter(Getter)
{
    AddTag(Tag);
}

~ScalarAverageSampler() override {}

void SampleDataPoint(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives, std::vector<double>::iterator& BufferIterator) override
{
    *BufferIterator = mGetter(rGeometry,rShapeFunctions,rShapeDerivatives);
    ++BufferIterator;
}

void OutputHeader(
    std::ofstream& rOutStream,
    const std::string& rSeparator) const override
{
    rOutStream << "<" << GetTag(0) << ">" << rSeparator;
}

private:

std::function<double(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives)> mGetter;

};


/// This class manages the computation of the average of a vector quantity.
template< class VectorType >
class VectorAverageSampler: public StatisticsSampler
{
public:

/// Initialize a sampler for a vector average quantity.
/** @param Gettr Helper function to evaluate the quantity of interest.
 *  @param VectorSize Number of components of the vector.
 *  @param Tags Names to identify each of the components on result files.
 *  @see Internals::MakeSamplerAtLocalCoordinate for the definition of the Getter argument.
 */
VectorAverageSampler(
    std::function<VectorType(const Geometry< Node<3> >&, const Vector&, const Matrix&)> Getter,
    std::size_t VectorSize,
    std::vector<std::string>& Tags)
    :
    StatisticsSampler(VectorSize),
    mGetter(Getter)
{
    for (auto it_tag = Tags.begin(); it_tag != Tags.end(); ++it_tag)
    {
        AddTag(*it_tag);
    }
}

~VectorAverageSampler() override {}

void SampleDataPoint(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives, std::vector<double>::iterator& BufferIterator) override {
    VectorType result = mGetter(rGeometry,rShapeFunctions,rShapeDerivatives);
    for (unsigned int i = 0; i < this->GetSize(); i++) {
        *BufferIterator = result[i];
        ++BufferIterator;
    }
}

void OutputHeader(
    std::ofstream& rOutStream,
    const std::string& rSeparator) const override
{
    for (std::size_t i = 0; i < GetSize(); i++)
    {
        rOutStream << "<" << GetTag(i) << ">" << rSeparator;
    }
}

private:

std::function<VectorType(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives)> mGetter;
};

/// This class manages the computation of the (co)variance between two given quantities (scalar or vector).
class VarianceSampler : public StatisticsSampler
{
public:

/// Initialize a sampler for the (co)variance between two given quantities.
/** Note that online recording of a variance requires also recording the averages of each quantity.
 *  @param pQuantity1 Sampler for the average of the first quantity of the variance.
 *  @param pQuantity2 Sampler for the average of the second quantity of the variance.
 */
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
        double delta_i = (NumberOfMeasurements-1)*new_measurement_i - current_total_i;
        for (std::size_t j = 0; j < mpQuantity2->GetSize(); j++)
        {
            double current_total_j = *(rCurrentStatistics.begin() + mpQuantity2->GetComponentOffset(j));
            double new_measurement_j = rNewMeasurement[mpQuantity2->GetComponentOffset(j)];
            double delta_j = (NumberOfMeasurements-1)*new_measurement_j - current_total_j;
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

const StatisticsSampler::Pointer GetQuantity1() const
{
    return mpQuantity1;
}

const StatisticsSampler::Pointer GetQuantity2() const
{
    return mpQuantity2;
}

void OutputHeader(
    std::ofstream& rOutStream,
    const std::string& rSeparator) const override
{
    for (std::size_t i = 0; i < mpQuantity1->GetSize(); i++)
    {
        for (std::size_t j = 0; j < mpQuantity2->GetSize(); j++)
        {
            rOutStream << "<" << mpQuantity1->GetTag(i) << "'" << mpQuantity2->GetTag(j) << "'>" << rSeparator;
        }
    }
}

private:

const StatisticsSampler::Pointer mpQuantity1;

const StatisticsSampler::Pointer mpQuantity2;

};

/// This class manages the computation of the variance for a given quantity.
/** It is essentially intended as a substitution of VarianceSampler for vector quantities.
 *  It takes advantage of the fact that variance(ui,uj) = variance(uj,ui), reducing
 *  the number of values to be tracked.
 */
class SymmetricVarianceSampler: public VarianceSampler
{
public:

/// Initialize a sampler for the variance of a vecor quantity.
/** Note that online recording of a variance requires also recording the average of the corresponding quantity.
 *  @param pQuantity1 Sampler for the average of the quantity of interest.
 */
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
        double delta_i = (NumberOfMeasurements-1)*new_measurement_i - current_total_i;
        for (std::size_t j = i; j < GetQuantity1()->GetSize(); j++)
        {
            double current_total_j = *(rCurrentStatistics.begin() + GetQuantity1()->GetComponentOffset(j));
            double new_measurement_j = rNewMeasurement[GetQuantity1()->GetComponentOffset(j)];
            double delta_j = (NumberOfMeasurements-1)*new_measurement_j - current_total_j;
            (*BufferIterator) = update_factor * delta_i * delta_j;
            ++BufferIterator;
        }
    }
}

void OutputHeader(
    std::ofstream& rOutStream,
    const std::string& rSeparator) const override
{
    for (std::size_t i = 0; i < GetQuantity1()->GetSize(); i++)
    {
        for (std::size_t j = i; j < GetQuantity1()->GetSize(); j++)
        {
            rOutStream << "<" << GetQuantity1()->GetTag(i) << "'" << GetQuantity1()->GetTag(j) << "'>" << rSeparator;
        }
    }
}


};

/// This class manages the computation of the (co)variance when one or both variables are Vector components.
/** Note that it is only prepared to accept scalar or vector component arguments, for full vectors, use
 *  the general VarianceSampler.
 */
class ComponentwiseVarianceSampler: public VarianceSampler
{
public:

/// Initialize a sampler for the (co)variance between two given quantities.
/** Note that online recording of a variance requires also recording the averages of each quantity.
 *  When using this class with a scalar quantity, pass 0 as its ComponentIndex.
 *  @param pQuantity1 Sampler for the average of the first quantity of the variance.
 *  @param ComponentIndex1 Index of the component of pQuantity1 to be used in the (co)variance.
 *  @param pQuantity2 Sampler for the average of the second quantity of the variance.
 *  @param ComponentIndex2 Index of the component of pQuantity2 to be used in the (co)variance.
 */
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
    double delta_i = (NumberOfMeasurements-1)*new_measurement_i - current_total_i;

    double current_total_j = *(rCurrentStatistics.begin() + GetQuantity2()->GetComponentOffset(mComponent2));
    double new_measurement_j = rNewMeasurement[GetQuantity2()->GetComponentOffset(mComponent2)];
    double delta_j = (NumberOfMeasurements-1)*new_measurement_j - current_total_j;
    (*BufferIterator) = update_factor * delta_i * delta_j;
    ++BufferIterator;
}


void OutputHeader(
    std::ofstream& rOutStream,
    const std::string& rSeparator) const override
{
    rOutStream << "<" << GetQuantity1()->GetTag(mComponent1) << "'" << GetQuantity2()->GetTag(mComponent2) << "'>" << rSeparator;
}

private:

std::size_t mComponent1;

std::size_t mComponent2;

};

/// This class manages the computation of third order moments of scalars or vector components.
/** Note that, since vector correlations quickly involve a lot of cross-components (up to 27 in 3D),
 *  support is only provided for single-valued quantities (that is, combinations of scalars or vector components).
 */
class ThirdOrderCorrelationSampler : public StatisticsSampler
{
public:

/// Initialize a sampler for the third order moment involving three given quantities.
/** Note that online recording of a third order statistic requires
 *  recording the averages of each quantity and the covariances for each pair.
 *  When using this class with a scalar quantity, pass 0 as its ComponentIndex.
 *  For variances, use pVarianceij->GetComponentIndex() to obtain the corresponding VarianceComponent argument.
 *  @param pQuantity1 Sampler for the average of the first quantity of the third order moment.
 *  @param ComponentIndex1 Index of the component of pQuantity1 to be used in calculation.
 *  @param pQuantity2 Sampler for the average of the second quantity of the third order moment.
 *  @param ComponentIndex2 Index of the component of pQuantity2 to be used in calculation.
 *  @param pQuantity3 Sampler for the average of the third quantity of the third order moment.
 *  @param ComponentIndex3 Index of the component of pQuantity3 to be used in calculation.
 *  @param pVariance12 Sampler for the covariance of the first and second quantities.
 *  @param VarianceComponent12 Index of the component of pVariance12 to be used in calculation.
 *  @param pVariance13 Sampler for the covariance of the first and third quantities.
 *  @param VarianceComponent13 Index of the component of pVariance13 to be used in calculation.
 *  @param pVariance23 Sampler for the covariance of the second and third quantities.
 *  @param VarianceComponent23 Index of the component of pVariance23 to be used in calculation.
 */
ThirdOrderCorrelationSampler(
    const StatisticsSampler::Pointer pQuantity1,
    const std::size_t QuantityComponent1,
    const StatisticsSampler::Pointer pQuantity2,
    const std::size_t QuantityComponent2,
    const StatisticsSampler::Pointer pQuantity3,
    const std::size_t QuantityComponent3,
    const StatisticsSampler::Pointer pVariance12,
    const std::size_t VarianceComponent12,
    const StatisticsSampler::Pointer pVariance13,
    const std::size_t VarianceComponent13,
    const StatisticsSampler::Pointer pVariance23,
    const std::size_t VarianceComponent23):
    StatisticsSampler(1),
    mpQuantity1(pQuantity1),
    mpQuantity2(pQuantity2),
    mpQuantity3(pQuantity3),
    mComponent1(QuantityComponent1),
    mComponent2(QuantityComponent2),
    mComponent3(QuantityComponent3),
    mpVariance12(pVariance12),
    mpVariance13(pVariance13),
    mpVariance23(pVariance23),
    mVarianceComponent12(VarianceComponent12),
    mVarianceComponent13(VarianceComponent13),
    mVarianceComponent23(VarianceComponent23)
{}

void SampleDataPoint(
    std::vector<double>::iterator& BufferIterator,
    const StatisticsSampler::IntegrationPointDataView& rCurrentStatistics,
    const std::vector<double>& rNewMeasurement,
    const std::size_t NumberOfMeasurements) override
{
    const double update_factor_1 = 1.0 / ((NumberOfMeasurements-1)*NumberOfMeasurements);
    const double update_factor_2 = (NumberOfMeasurements-2)*update_factor_1*update_factor_1;

    const std::size_t value_1_offset = mpQuantity1->GetComponentOffset(mComponent1);
    const std::size_t value_2_offset = mpQuantity2->GetComponentOffset(mComponent2);
    const std::size_t value_3_offset = mpQuantity3->GetComponentOffset(mComponent3);

    const std::size_t variance_12_offset = mpVariance12->GetComponentOffset(mVarianceComponent12);
    const std::size_t variance_13_offset = mpVariance13->GetComponentOffset(mVarianceComponent13);
    const std::size_t variance_23_offset = mpVariance23->GetComponentOffset(mVarianceComponent23);

    double current_total_1 = *(rCurrentStatistics.begin() + value_1_offset);
    double new_measurement_1 = rNewMeasurement[value_1_offset];
    double delta_1 = (NumberOfMeasurements-1)*new_measurement_1 - current_total_1;

    double current_total_2 = *(rCurrentStatistics.begin() + value_2_offset);
    double new_measurement_2 = rNewMeasurement[value_2_offset];
    double delta_2 = (NumberOfMeasurements-1)*new_measurement_2 - current_total_2;

    double current_total_3 = *(rCurrentStatistics.begin() + value_3_offset);
    double new_measurement_3 = rNewMeasurement[value_3_offset];
    double delta_3 = (NumberOfMeasurements-1)*new_measurement_3 - current_total_3;

    double current_variance_12 = *(rCurrentStatistics.begin() + variance_12_offset);
    double current_variance_13 = *(rCurrentStatistics.begin() + variance_13_offset);
    double current_variance_23 = *(rCurrentStatistics.begin() + variance_23_offset);

    double update = update_factor_2 * delta_1 * delta_2 * delta_3;
    update -= update_factor_1 * current_variance_12 * delta_3;
    update -= update_factor_1 * current_variance_13 * delta_2;
    update -= update_factor_1 * current_variance_23 * delta_1;

    (*BufferIterator) = update;
    ++BufferIterator;
}

virtual void OutputHeader(
    std::ofstream& rOutStream,
    const std::string& rSeparator) const override
{
    rOutStream << "<"
    << mpQuantity1->GetTag(mComponent1) << "'"
    << mpQuantity2->GetTag(mComponent2) << "'"
    << mpQuantity3->GetTag(mComponent3) << "'>" << rSeparator;
}

private:

const StatisticsSampler::Pointer mpQuantity1;

const StatisticsSampler::Pointer mpQuantity2;

const StatisticsSampler::Pointer mpQuantity3;

const std::size_t mComponent1;

const std::size_t mComponent2;

const std::size_t mComponent3;


const StatisticsSampler::Pointer mpVariance12;

const StatisticsSampler::Pointer mpVariance13;

const StatisticsSampler::Pointer mpVariance23;

const std::size_t mVarianceComponent12;

const std::size_t mVarianceComponent13;

const std::size_t mVarianceComponent23;


};

namespace Internals {

class MakeSamplerAtLocalCoordinate {
public:
    static std::function<double(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives)> ValueGetter(const Variable<double>& rVariable) {
        return [rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> double {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeFunctions.size()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            double value = 0.0;
            for (unsigned int i =  0; i < rGeometry.size(); i++) {
                value += rGeometry[i].FastGetSolutionStepValue(rVariable) * rShapeFunctions[i];
            }
            return value;
        };
    }


    static std::function< array_1d<double,3>(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) > ValueGetter(const Variable<array_1d<double,3>>& rVariable) {
        return [rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> array_1d<double,3> {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeFunctions.size()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            array_1d<double,3> value = ZeroVector(3);
            for (unsigned int i =  0; i < rGeometry.size(); i++) {
                value += rGeometry[i].FastGetSolutionStepValue(rVariable) * rShapeFunctions[i];
            }
            return value;
        };
    }

    static std::function< Matrix(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) > ValueGetter(const Variable<Matrix>& rVariable) {
        return [rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> Matrix {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeFunctions.size()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            Matrix value = ZeroMatrix(3,3);
            for (unsigned int i =  0; i < rGeometry.size(); i++) {
                value += rGeometry[i].FastGetSolutionStepValue(rVariable) * rShapeFunctions[i];
            }
            return value;
        };
    }


    static std::function< array_1d<double,3>(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) > GradientGetter(const Variable<double>& rVariable) {
        return [rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> array_1d<double,3> {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeDerivatives.size1()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            array_1d<double,3> gradient = ZeroVector(3);
            for (unsigned int n =  0; n < rGeometry.size(); n++) {
                const auto& value = rGeometry[n].FastGetSolutionStepValue(rVariable);
                for (unsigned int i = 0; i < rShapeDerivatives.size2(); i++)
                    gradient[i] += value * rShapeDerivatives(n,i);
            }
            return gradient;
        };
    }


    static std::function< Matrix(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) > GradientGetter(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, Variable<array_1d<double,3>>& rVariable) {
        return [rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> Matrix {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeDerivatives.size1()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            Matrix gradient = ZeroMatrix(3,3);
            for (unsigned int n =  0; n < rGeometry.size(); n++) {
                const auto& value = rGeometry[n].FastGetSolutionStepValue(rVariable);
                for (unsigned int i = 0; i < rShapeDerivatives.size2(); i++)
                {
                    for (unsigned int j = 0; j < rShapeDerivatives.size2(); j++)
                        gradient(i,j) += value[i] * rShapeDerivatives(n,j); //dui/dxj
                }
            }
            return gradient;
        };
    }
};

}

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
