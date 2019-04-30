//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//


#if !defined(KRATOS_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_BOUNDING_BOX_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>



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
template<class TPointType,  class TPointerType>
class BoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BoundingBox
    KRATOS_CLASS_POINTER_DEFINITION(BoundingBox);

    typedef TPointerType  PointerType;
    typedef BoundingBox<TPointType, PointerType > BoundingBoxType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BoundingBox() : mHighPoint(), mLowPoint()
    {
        std::cout<< "Calling empty constructor" <<std::endl;
    }

    explicit BoundingBox(const TPointType& Point) :  mLowPoint(Point), mHighPoint(Point)
    {
    }

    BoundingBox(const TPointType& lowpoint, const TPointType& highpoint ) :  mLowPoint(lowpoint), mHighPoint(highpoint)
    {
    }

    BoundingBox(const PointerType Object, const TPointType& lowpoint, const TPointType& highpoint) :  mLowPoint(lowpoint), mHighPoint(highpoint), mObject(Object)
    {
    }


    /// Destructor.
    virtual ~BoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    void Set(const TPointerType Object, const TPointType& lowpoint, const TPointType& highpoint)
    {
        mLowPoint  = lowpoint;
        mHighPoint = highpoint;
        mObject    = Object;

    }

    TPointType const& HighPoint() const
    {
        return mHighPoint;
    }


    TPointType& HighPoint()
    {
        return mHighPoint;
    }

    TPointType const& LowPoint() const
    {
        return mLowPoint;
    }


    TPointType& LowPoint()
    {
        return mLowPoint;
    }


    /// Assignment operator.
    BoundingBox& operator=(BoundingBox const& rOther)
    {
        mHighPoint = rOther.mHighPoint;
        mLowPoint  = rOther.mLowPoint;
        mObject     = rOther.mObject;
        return *this;
    }


    /// Copy constructor.
    BoundingBox(BoundingBox const& rOther) :
        mHighPoint(rOther.mHighPoint),
        mLowPoint(rOther.mLowPoint),
        mObject(rOther.mObject)
    {
    }



    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BoundingBox";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << mHighPoint << " , " << mLowPoint;
    }

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
    TPointType mHighPoint;
    TPointType mLowPoint;
    TPointerType  mObject;


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


}; // Class BoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TPointType, class TPointerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  BoundingBox<TPointType, TPointerType>& rThis);

/// output stream function
template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BoundingBox<TPointType, TPointerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


