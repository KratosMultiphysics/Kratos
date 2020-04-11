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

///@name Kratos Classes
///@{

/// This class is the result of the spatial searches.
/** It provides:
 *  - A global pointer to the object found.
 *  - Distance to the object if IsDistanceCalculated() is true
 *  - IsObjectFound if for example search nearest fails or not
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

    /// Constructor with the resulted object   
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

    /// Returns the global pointer to the object
	TPointerType Get() { 
        return mpObject; 
    }
    
    /// Returns a const global pointer to the object
	TPointerType const Get() const { 
        return mpObject; 
    }

    /// Set the object to be pointed
	void Set(TObjectType* pObject) {
		mpObject = pObject;
		mIsObjectFound = true;
	}

    /// Getting the result distance
	double GetDistance() const { 
        return mDistance; 
    }

    /// Setting the result distance
	void SetDistance(double TheDistance) {
		mDistance = TheDistance;
		mIsDistanceCalculated = true;
	}

    ///@}
    ///@name Inquiry
    ///@{

    /// Returns true if the object is set 
	bool IsObjectFound() const { 
        return mIsObjectFound; 
    }

    /// Returns true if the distance is set
	bool IsDistanceCalculated() const { 
        return mIsDistanceCalculated; 
    }

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
    ///@name Private Operations
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


