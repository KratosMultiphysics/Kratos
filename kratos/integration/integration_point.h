//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                    
//

#if !defined(KRATOS_INTEGRATION_POINT_H_INCLUDED )
#define  KRATOS_INTEGRATION_POINT_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

#include <numeric>
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "geometries/point.h"


namespace Kratos
{

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
template<std::size_t TDimension, class TDataType = double, class TWeightType = double>
class IntegrationPoint : public Point
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IntegrationPoint
    KRATOS_CLASS_POINTER_DEFINITION(IntegrationPoint);

    typedef Point BaseType;

    typedef Point PointType;

    typedef typename Point::CoordinatesArrayType CoordinatesArrayType;

    typedef typename Point::IndexType IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IntegrationPoint() : BaseType(), mWeight()
    {
    }

    /// 1d constructor.
    explicit IntegrationPoint(TDataType const& NewX) : BaseType(NewX), mWeight()
    {
    }

    /// 1d constructor.
    IntegrationPoint(TDataType const& NewX, TDataType const& NewW) : BaseType(NewX), mWeight(NewW)
    {
    }

    /// 2d constructor.
    IntegrationPoint(TDataType const& NewX, TDataType const& NewY, TDataType const& NewW) : BaseType(NewX, NewY), mWeight(NewW)
    {
    }

    /// 3d constructor
    IntegrationPoint(TDataType const& NewX, TDataType const& NewY, TDataType const& NewZ, TDataType const& NewW) : BaseType(NewX, NewY, NewZ), mWeight(NewW)
    {
    }

    /** Copy constructor. Initialize this integration point with the coordinates
    of given integration point.*/
    IntegrationPoint(IntegrationPoint const& rOtherIntegrationPoint)
        : BaseType(rOtherIntegrationPoint), mWeight(rOtherIntegrationPoint.mWeight) {}

    /** Copy constructor from other dimension integration
    point. Initialize this integration point with the
    coordinates of given integration point.*/
    template<std::size_t TOtherDimension>
    IntegrationPoint(IntegrationPoint<TOtherDimension, TDataType, TWeightType> const& rOtherIntegrationPoint)
        : BaseType(rOtherIntegrationPoint), mWeight(rOtherIntegrationPoint.mWeight) {}

    /** Point constructor. Initialize this integration point with the coordinates
    of given point.*/
    explicit IntegrationPoint(PointType const& rOtherPoint)
        : BaseType(rOtherPoint), mWeight() {}

    /** Point constructor with weight. Initialize this integration point with the coordinates
    of given point and given weight*/
    IntegrationPoint(PointType const& rOtherPoint, TWeightType NewWeight)
        : BaseType(rOtherPoint), mWeight(NewWeight) {}

    /** Constructor using coordinates stored in given array. Initialize
    this integration point with the coordinates in the array. Integration wieght initializes to zero. */
    explicit IntegrationPoint(CoordinatesArrayType const& rOtherCoordinates)
        : BaseType(rOtherCoordinates), mWeight() {}

    /** Constructor using coordinates stored in given array. Initialize
    this integration point with the coordinates in the array and given wieght. */
    IntegrationPoint(CoordinatesArrayType const& rOtherCoordinates, TWeightType NewWeight)
        : BaseType(rOtherCoordinates), mWeight(NewWeight) {}

    /** Constructor using coordinates stored in given array. Initialize
    this integration point with the coordinates in the array. Integration wieght initializes to zero. */
    template<class TVectorType>
    explicit IntegrationPoint(vector_expression<TVectorType> const&  rOtherCoordinates)
        : BaseType(rOtherCoordinates), mWeight() {}

    /** Constructor using coordinates stored in given array. Initialize
    this integration point with the coordinates in the array and given wieght. */
    template<class TVectorType>
    IntegrationPoint(vector_expression<TVectorType> const&  rOtherCoordinates, TWeightType NewWeight)
        : BaseType(rOtherCoordinates), mWeight(NewWeight) {}

    /** Constructor using coordinates stored in given std::vector. Initialize
    this integration point with the coordinates in the array. Integration wieght initializes to zero. */
    explicit IntegrationPoint(std::vector<TDataType> const&  rOtherCoordinates)
        : BaseType(rOtherCoordinates), mWeight() {}

    /** Constructor using coordinates stored in given std::vector. Initialize
    this integration point with the coordinates in the array and given wieght. */
    IntegrationPoint(std::vector<TDataType> const&  rOtherCoordinates, TWeightType NewWeight)
        : BaseType(rOtherCoordinates), mWeight(NewWeight) {}

    /// Destructor.
    ~IntegrationPoint() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IntegrationPoint& operator=(const IntegrationPoint& rOther)
    {
        BaseType::operator =(rOther);
        mWeight = rOther.mWeight;

        return *this;
    }

    /// Point assignment operator.
    IntegrationPoint& operator=(const PointType& OtherPoint)
    {
        BaseType::operator =(OtherPoint);

        return *this;
    }

    bool operator==(const IntegrationPoint& rOther)
    {
        return (mWeight == rOther.mWeight) && (BaseType::operator ==(rOther));
    }

    /// Assignment operator with different dimension.
    template<std::size_t TOtherDimension>
    IntegrationPoint& operator=(const IntegrationPoint<TOtherDimension>& rOther)
    {
        BaseType::operator =(rOther);
        mWeight = rOther.Weight();

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    TWeightType Weight() const
    {
        return mWeight;
    }

    TWeightType& Weight()
    {
        return mWeight;
    }

    void SetWeight(TWeightType NewWeight)
    {
        mWeight = NewWeight;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << TDimension << " dimensional integration point";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << TDimension << " dimensional integration point";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        if(!TDimension)
            return;

        rOStream << "("  << this->operator[](0);

        for(IndexType i = 1 ; i < TDimension  ; i++)
            rOStream << " , " << this->operator[](i);

        rOStream << "), weight = " << mWeight;
    }



    ///@}
    ///@name Friends
    ///@{

    template<std::size_t TOtherDimension, class TOtherDataType, class TOtherWeightType> friend class IntegrationPoint;

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

    TWeightType mWeight;

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


    ///@}

}; // Class IntegrationPoint

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<std::size_t TDimension, class TDataType, class TWeightType>
inline std::istream& operator >> (std::istream& rIStream,
                                  IntegrationPoint<TDimension, TDataType, TWeightType>& rThis);

/// output stream function
template<std::size_t TDimension, class TDataType, class TWeightType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const IntegrationPoint<TDimension, TDataType, TWeightType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_POINT_H_INCLUDED  defined 


