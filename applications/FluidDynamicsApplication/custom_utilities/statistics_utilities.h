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
#include "geometries/geometry.h"
#include "includes/ublas_interface.h"

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

StatisticsSampler(unsigned int NumValues, unsigned int Offset = 0):
    mNumValues(NumValues),
    mOffset(Offset)
{}

KRATOS_CLASS_POINTER_DEFINITION(StatisticsSampler);

virtual ~StatisticsSampler() {}

virtual void SampleDataPoint(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives, std::vector<double>::iterator& BufferIterator)
{}

unsigned int GetSize() const {
    return mNumValues;
}

virtual unsigned int GetValueOffset() const
{
    KRATOS_ERROR << "Method not implemented" << std::endl;
}

virtual unsigned int GetComponentOffset(unsigned int i) const
{
    KRATOS_ERROR << "Method not implemented" << std::endl;
}

virtual unsigned int GetComponentOffset(unsigned int i, unsigned int j) const
{
    KRATOS_ERROR << "Method not implemented" << std::endl;
}

unsigned int GetOffset() const
{
    return mOffset;
}

void SetOffset(std::size_t Offset) {
    mOffset = Offset;
}

template<class TIteratorType>
void Finalize(TIteratorType& rBufferIterator, std::size_t NumberOfMeasurements)
{
    for (std::size_t i = 0; i < mNumValues; i++) {
        *rBufferIterator /= NumberOfMeasurements;
        KRATOS_WATCH(*rBufferIterator);
        ++rBufferIterator;
    }
}


private:

unsigned int mNumValues;

unsigned int mOffset;

};

class ScalarAverageSampler: public StatisticsSampler
{
public:

ScalarAverageSampler(unsigned int Offset, std::function<double(const Geometry< Node<3> >&, const Vector&, const Matrix&)> Getter):
  StatisticsSampler(1, Offset),
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

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_STATISTICS_UTILITIES_H_INCLUDED  defined
