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


#if !defined(KRATOS_SPATIAL_SEARCH_RESULT_H_INCLUDED )
#define  KRATOS_SPATIAL_SEARCH_RESULT_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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
template <typename TObjectType>
class SpatialSearchResult
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialSearchResult
    KRATOS_CLASS_POINTER_DEFINITION(SpatialSearchResult);
	
    using TPointerType = GlobalPointer<TObjectType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
   	SpatialSearchResult() : mpObject(nullptr), mDistance(0.00), mIsObjectFound(false), mIsDistanceCalculated(false) {}
       
	SpatialSearchResult(TObjectType* pObject) : mpObject(pObject), mDistance(0.00), mIsObjectFound(false), mIsDistanceCalculated(false) {
		if (mpObject.get() != nullptr)
			mIsObjectFound = true;
	}

	SpatialSearchResult(SpatialSearchResult const& /* Other */) = default;

	SpatialSearchResult(SpatialSearchResult&& /* Other */) = default;

    /// Destructor.
    virtual ~SpatialSearchResult(){}

    ///@}
    ///@name Operators
    ///@{

        SpatialSearchResult& operator=(SpatialSearchResult const& /*Other*/) = default;


    ///@}
    ///@name Operations
    ///@{

	void Reset() {
		mpObject = mpObject(nullptr);
		mDistance = 0.00;
		mIsObjectFound = false;
		mIsDistanceCalculated = false;
	}

    ///@}
    ///@name Access
    ///@{

	TPointerType Get() { return mpObject; }
	TPointerType const Get() const { return mpObject; }
	void Set(TObjectType* pObject) {
		mpObject = pObject;
		mIsObjectFound = true;
	}
	double GetDistance() const { return mDistance; }
	void SetDistance(double TheDistance) {
		mDistance = TheDistance;
		mIsDistanceCalculated = true;
	}

    ///@}
    ///@name Inquiry
    ///@{

	bool IsObjectFound() const { return mIsObjectFound; }

	bool IsDistanceCalculated() const { return mIsDistanceCalculated; }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
            buffer << "SpatialSearchResult" ;
            return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SpatialSearchResult";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

private:
    ///@name Member Variables
    ///@{

	TPointerType mpObject;
	double mDistance;
	bool mIsObjectFound;
	bool mIsDistanceCalculated;

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

}; // Class SpatialSearchResult

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template <typename TObjectType>
inline std::istream& operator >> (std::istream& rIStream,
                SpatialSearchResult<TObjectType>& rThis){
                    return rIStream;
                }

/// output stream function
template <typename TObjectType>
inline std::ostream& operator << (std::ostream& rOStream,
                const SpatialSearchResult<TObjectType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPATIAL_SEARCH_RESULT_H_INCLUDED  defined


