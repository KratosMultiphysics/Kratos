// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kevin Braun, https://github.com/MFusseder
//

#ifndef VECTOR_MATH_H
#define VECTOR_MATH_H

// System includes

// External includes

// Project includes
#include "includes/process_info.h"

namespace Kratos
{

/** \brief DerivativeBuilder
*
* The purpose of this class is to perform varies calculations with vectors that are filled with certain data types (e.g. array_1d<double, 3>, Matrix etc.)
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) VectorMath
{
public:

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;


// These functions add two vectors 

template <typename TDataType>
static void Addition( std::vector<std::vector<TDataType>>& rOutput, const std::vector<std::vector<TDataType>>& rInput )
{
    KRATOS_ERROR_IF_NOT( rOutput.size() == rInput.size() ) << "Not Possible to add 2 vectors of different sizes" <<std::endl;
        
    for(IndexType i = 0; i < rOutput.size(); ++i)
        Addition(rOutput[i], rInput[i]);
}

template <typename TDataType>
static void Addition( std::vector<TDataType>& rOutput, const std::vector<TDataType>& rInput )
{
    KRATOS_ERROR_IF_NOT( rOutput.size() == rInput.size() ) << "Not Possible to add 2 vectors of different sizes" <<std::endl;
        
    for(IndexType i = 0; i < rOutput.size(); ++i)
        Addition(rOutput[i], rInput[i]);
}

static void Addition( array_1d<double, 3>& rOutput, const array_1d<double, 3>& rInput )
{   
    for (IndexType dir_it = 0; dir_it < 3; ++dir_it)
        rOutput[dir_it] += rInput[dir_it];
}

static void Addition( Matrix& rOutput, const Matrix& rInput )
{    
   for( IndexType i = 0; i < rOutput.size1(); ++i )
        for( IndexType j = 0; j < rOutput.size2(); ++j )
            rOutput(i,j) += rInput(i,j); 
}


// These functions set all entries of a vector to zero.

template <typename TDataType>
static void SetToZero( std::vector<std::vector<TDataType>>& rOutput )
{
    for (IndexType i = 0; i < rOutput.size(); ++i)
        SetToZero(rOutput[i]);
} 

template <typename TDataType>
static void SetToZero( std::vector<TDataType>& rOutput )
{
    for (IndexType i = 0; i < rOutput.size(); ++i)
        SetToZero(rOutput[i]);
}

static void SetToZero( array_1d<double, 3>& rOutput )
{
    for (IndexType dir_it = 0; dir_it < 3; ++dir_it)
        rOutput[dir_it] = 0;
}

static void SetToZero(Matrix& rOutput)
{    
   for( IndexType i = 0; i < rOutput.size1(); ++i )
        for( IndexType j = 0; j < rOutput.size2(); ++j )
            rOutput(i,j) = 0; 
}

// These functions size the components of a vector

template <typename TDataType>
static void SizeVectorComponents(std::vector<std::vector<TDataType>>& rVector)
{   
    for (IndexType i = 0; i < rVector.size(); ++i)
        SizeVectorComponents(rVector[i]);   
}


template <typename TDataType>
static void SizeVectorComponents(std::vector<TDataType>& rVector)
{   
    for (IndexType i = 0; i < rVector.size(); ++i)
        SizeVectorComponents(rVector[i]);   
}

static void SizeVectorComponents(array_1d<double,3>& rComponent){}

static void SizeVectorComponents(Matrix& rComponent)
{   
    if (rComponent.size1() != 3 || rComponent.size1() != 3 )
    rComponent.resize(3,3);
}

// These functions multiply a vector by a scalar factor

template <typename TDataType>
static void MultiplyByFactor( std::vector<std::vector<TDataType>>& rOutput ,const double& rFactor )
{
    for (IndexType i = 0; i < rOutput.size(); ++i)
        MultiplyByFactor(rOutput[i], rFactor);
}

template <typename TDataType>
static void MultiplyByFactor( std::vector<TDataType>& rOutput ,const double& rFactor )
{
    for (IndexType i = 0; i < rOutput.size(); ++i)
        MultiplyByFactor(rOutput[i], rFactor);
}

static void MultiplyByFactor(array_1d<double, 3>& rOutput , const double& rFactor)
{    
   for( IndexType dir_it = 0; dir_it < 3; ++dir_it )
        rOutput[dir_it] *= rFactor; 
}

static void MultiplyByFactor(Matrix& rOutput , const double& rFactor)
{    
   for( IndexType i = 0; i < rOutput.size1(); ++i )
        for( IndexType j = 0; j < rOutput.size2(); ++j )
            rOutput(i,j) *= rFactor; 
}

// These functions compute a scalar product of two vectors. One of the vectors is filled with double components, the other with
// an arbritrary type (e.g. array_1d<double, 3>, array_1d<double, 6>, Matrix etc.). Needed in direct sensitivity Postprocess. 

template <typename TDataType>
static void ScalarProduct(const std::vector<std::vector<TDataType>>& rFactor1,
                            const Vector& rFactor2,
                            std::vector<TDataType>& rScalarProduct)
{    
    SetToZero(rScalarProduct);    
    for(IndexType i = 0; i < rFactor1.size(); ++i)
        for( IndexType j = 0; j < rFactor1[0].size(); ++j )
            MultiplyFactorsForScalarProduct(rScalarProduct[j], rFactor1[i][j], rFactor2[i]);
}

template <typename TDataType>
static void ScalarProduct(const std::vector<TDataType>& rFactor1 ,
                            const Vector& rFactor2,
                            TDataType& rScalarProduct)
{    
    SetToZero(rScalarProduct);
    for( IndexType i = 0; i < rFactor1.size(); ++i )
        MultiplyFactorsForScalarProduct(rScalarProduct, rFactor1[i], rFactor2[i]);
} 

// These functions compute the product of a vector filled with templated data types (e.g. array_1d<double, 3>, Matrix etc.)
// and a component of TDataType. Needed in direct sensitivity Postprocess.

template <typename TDataType>
static void Product(const std::vector<TDataType>& rFactor1,
                            const TDataType& rFactor2,
                            TDataType& rScalarProduct)
{
    SetToZero(rScalarProduct);
    for(IndexType i = 0; i < rFactor1.size(); ++i)
        MultiplyFactorsForProduct(rScalarProduct, rFactor1[i], rFactor2);
}

    
private:

static void MultiplyFactorsForScalarProduct( array_1d<double, 3>& rScalarProduct ,const array_1d<double, 3>& rFactor1 , const double& rFactor2)
{    
   for( IndexType dir_it = 0; dir_it < 3; ++dir_it )
        rScalarProduct[dir_it] += rFactor1[dir_it] * rFactor2; 
}

static void MultiplyFactorsForScalarProduct( Matrix& rScalarProduct ,const Matrix& rFactor1 , const double& rFactor2)
{     
   for( IndexType i = 0; i < rFactor1.size1(); ++i )
        for( IndexType j = 0; j < rFactor1.size2(); ++j )
            rScalarProduct(i,j) += rFactor1(i,j) * rFactor2; 
}

static void MultiplyFactorsForProduct( array_1d<double, 3>& rScalarProduct ,const array_1d<double, 3>& rFactor1 , const array_1d<double, 3>& rFactor2)
{    
   for( IndexType dir_it = 0; dir_it < 3; ++dir_it )
        rScalarProduct[dir_it] += rFactor1[dir_it] * rFactor2[dir_it]; 
}

static void MultiplyFactorsForProduct( Matrix& rScalarProduct , const Matrix& rFactor1 , const Matrix& rFactor2)
{    
   KRATOS_ERROR << "MultiplyFactorsForProduct not yet implemented for matrix type." << std::endl;
}


}; // class DerivativeBuilder


} /* namespace Kratos.*/

#endif /* VECTOR_MATH_H defined */