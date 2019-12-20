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

#if !defined(KRATOS_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_BOUNDING_BOX_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Representing a bounding box by storing the min and max points
/** It stores the min and max points and have constructor for it construction with any container of points.
 *  TPointType should provide access operator [] to its coordinate and deep copy operator=
*/
template <typename TPointType>
class BoundingBox
{
	static constexpr int Dimension = 3;
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BoundingBox
    KRATOS_CLASS_POINTER_DEFINITION(BoundingBox);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    BoundingBox(){
        for (int i = 0; i < Dimension; i++)
        {
            mMinPoint[i] = 0.00;
            mMaxPoint[i] = 0.00;
        }
    };

	BoundingBox(TPointType const& MinPoint, TPointType const& MaxPoint) :
		mMinPoint(MinPoint), mMaxPoint(MaxPoint) {}

    /// Copy constructor
	BoundingBox( const BoundingBox &Other) :
		mMinPoint( Other.mMinPoint), mMaxPoint( Other.mMaxPoint) {}


    /// Construction with container of points.
	template<typename TIteratorType>
	BoundingBox(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd) {
        Set(PointsBegin, PointsEnd);
	}

    /// Destructor.
    virtual ~BoundingBox(){}

    ///@}
    ///@name Operators
    ///@{


    /// Assignment operator.
    BoundingBox& operator=(BoundingBox const& rOther){
        mMinPoint = rOther.mMinPoint;
        mMaxPoint = rOther.mMaxPoint;

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    template<typename TIteratorType>
    void Set(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd){
        if (PointsBegin == PointsEnd) {
            for (int i = 0; i < Dimension; i++) {
                mMinPoint[i] = 0.00;
                mMaxPoint[i] = 0.00;
            }
            return;
        }

        for (int i = 0; i < Dimension; i++) {
            mMinPoint[i] = (*PointsBegin)[i];
            mMaxPoint[i] = (*PointsBegin)[i];
        }

        Extend(PointsBegin, PointsEnd);
    }

    template<typename TIteratorType>
    void Extend(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd){

        for (TIteratorType i_point = PointsBegin; i_point != PointsEnd; i_point++){
            for (int i = 0; i < Dimension; i++)
            {
                if ((*i_point)[i] < mMinPoint[i]) mMinPoint[i] = (*i_point)[i];
                if ((*i_point)[i] > mMaxPoint[i]) mMaxPoint[i] = (*i_point)[i];
            }
        }
    }

    void Extend(double Margin){
        for (int i = 0; i < Dimension; i++){
            mMinPoint[i] -= Margin;
            mMaxPoint[i] += Margin;
        }

    }


    ///@}
    ///@name Access
    ///@{

	TPointType& GetMinPoint() { return mMinPoint; }
	TPointType const& GetMinPoint() const { return mMinPoint; }

	TPointType& GetMaxPoint() { return mMaxPoint; }
	TPointType const& GetMaxPoint() const { return mMaxPoint; }


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
        buffer << "BoundingBox" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "BoundingBox";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    
    TPointType mMinPoint;
	TPointType mMaxPoint;


    ///@}

}; // Class BoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template <typename TPointType>
inline std::istream& operator >> (std::istream& rIStream,
                BoundingBox<TPointType>& rThis){
                    return rIStream;
                }

/// output stream function
template <typename TPointType>
inline std::ostream& operator << (std::ostream& rOStream,
                const BoundingBox<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOUNDING_BOX_H_INCLUDED  defined


