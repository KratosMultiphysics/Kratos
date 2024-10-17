//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

// External includes

// Project includes

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

template< class  TConfigure>
class Cell
{
public:
    ///@name Type Definitions
    ///@{

    /// configure types
    typedef std::size_t  SizeType;
    typedef typename TConfigure::PointType               PointType;
    typedef typename TConfigure::PointerType             PointerType;
    typedef typename TConfigure::ContainerType           ContainerType;
    typedef typename TConfigure::IteratorType            IteratorType;
    typedef typename TConfigure::ResultContainerType     ResultContainerType;
    typedef typename TConfigure::ResultIteratorType      ResultIteratorType;
    typedef typename TConfigure::DistanceIteratorType	 DistanceIteratorType;

    typedef std::vector<PointerType>     LocalContainerType;
    typedef typename LocalContainerType::iterator LocalIteratorType;

    /// Pointer definition of Cell
    KRATOS_CLASS_POINTER_DEFINITION(Cell);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Cell()
    {
    }

    /// Destructor.
    virtual ~Cell() {}


    void Add(const PointerType& ThisObject)
    {
        mObjects.push_back(ThisObject);
    }

    void Remove(const PointerType& ThisObject)
    {
        mObjects.erase(std::remove(mObjects.begin(),mObjects.end(),ThisObject),mObjects.end());
    }

    // This is a zero based index of object in local array
    void Remove(const std::size_t Index)
    {
        std::swap(mObjects[Index], mObjects.back());
        mObjects.pop_back();
    }

    void Clear()
    {
        mObjects.clear();
    }


    void AllocateCell(const std::size_t size)
    {
        mObjects.reserve(size);
    }


    /// Assignment operator.
    Cell& operator=(Cell const& rOther)
    {
        mObjects   = rOther.mObjects;
        return *this;
    }


    /// Copy constructor.
    Cell(Cell const& rOther) :
        mObjects(rOther.mObjects)
    {
    }

    void SearchObjects(PointerType& rThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
    {
        for(LocalIteratorType i_object = Begin() ; i_object != End()  && NumberOfResults < MaxNumberOfResults ; i_object++)
        {
            if(TConfigure::Intersection(rThisObject, *i_object))
            {
                ResultIteratorType repeated_object = std::find(Result-NumberOfResults, Result, *i_object);
                if(repeated_object==Result)
                {
                    *Result   = *i_object;
                    Result++;
                    NumberOfResults++;
                }
            }
        }
    }

    void SearchObjects(PointerType& rThisObject, ResultContainerType& Result)
    {
        for(LocalIteratorType i_object = Begin() ; i_object != End(); i_object++)
        {
            if(TConfigure::Intersection(rThisObject, *i_object))
            {
                ResultIteratorType repeated_object = std::find(Result.begin(), Result.end(), *i_object);
                if(repeated_object==Result.end())
                {
                    Result.push_back(*i_object);
                }
            }
        }
    }

    void SearchObjectsExclusive(PointerType& rThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
    {
        for(LocalIteratorType i_object = Begin() ; i_object != End()  && NumberOfResults < MaxNumberOfResults ; i_object++)
        {
            if( rThisObject != *i_object )
            {
                if(TConfigure::Intersection(rThisObject, *i_object))
                {
                    ResultIteratorType repeated_object = std::find(Result-NumberOfResults, Result, *i_object);
                    if(repeated_object==Result)
                    {
                        *Result   = *i_object;
                        Result++;
                        NumberOfResults++;
                    }
                }
            }
        }
    }

    void SearchObjectsExclusive(PointerType& rThisObject, ResultContainerType& Result)
    {
        for(LocalIteratorType i_object = Begin() ; i_object != End(); i_object++)
        {
            if( rThisObject != *i_object )
            {
                if(TConfigure::Intersection(rThisObject, *i_object))
                {
                    ResultIteratorType repeated_object = std::find(Result.begin(), Result.end(), *i_object);
                    if(repeated_object==Result.end())
                    {
                        Result.push_back(*i_object);
                    }
                }
            }
        }
    }

    void SearchObjectsInRadius(PointerType& rThisObject, double const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
    {
        for(LocalIteratorType i_object = Begin() ; i_object != End()  && NumberOfResults < MaxNumberOfResults ; i_object++)
        {
            if(TConfigure::Intersection(rThisObject, *i_object, Radius))
            {
                ResultIteratorType repeated_object = std::find(Result-NumberOfResults, Result, *i_object);
                if(repeated_object==Result)
                {
                    *Result   = *i_object;
                    Result++;
                    NumberOfResults++;
                }
            }
        }
    }

    void SearchObjectsInRadiusExclusive(PointerType& rThisObject, double const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
    {
        for(LocalIteratorType i_object = Begin() ; i_object != End()  && NumberOfResults < MaxNumberOfResults ; i_object++)
        {
            if( rThisObject != *i_object )
            {
                if(TConfigure::Intersection(rThisObject, *i_object, Radius))
                {
                    ResultIteratorType repeated_object = std::find(Result-NumberOfResults, Result, *i_object);
                    if(repeated_object==Result)
                    {
                        *Result   = *i_object;
                        Result++;
                        NumberOfResults++;
                    }
                }
            }
        }
    }

    void SearchObjectsInRadius(PointerType& rThisObject, double const& Radius, ResultIteratorType& Result, DistanceIteratorType& Distances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
    {
        for(LocalIteratorType i_object = Begin() ; i_object != End()  && NumberOfResults < MaxNumberOfResults ; i_object++)
        {
            if(TConfigure::Intersection(rThisObject, *i_object, Radius))
            {
                ResultIteratorType repeated_object = std::find(Result-NumberOfResults, Result, *i_object);
                if(repeated_object==Result)
                {
                    double distance = 0;
                    TConfigure::Distance(rThisObject,*i_object,distance); // squared distance function
                    *Result   = *i_object;
                    Result++;
                    *Distances = distance;
                    Distances++;
                    NumberOfResults++;
                }
            }
        }
    }

    void SearchObjectsInRadiusExclusive(PointerType& rThisObject, double const& Radius, ResultIteratorType& Result, DistanceIteratorType& Distances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
    {
        for(LocalIteratorType i_object = Begin() ; i_object != End()  && NumberOfResults < MaxNumberOfResults ; i_object++)
        {
            if( rThisObject != *i_object )
            {
                if(TConfigure::Intersection(rThisObject, *i_object, Radius))
                {
                    ResultIteratorType repeated_object = std::find(Result-NumberOfResults, Result, *i_object);
                    if(repeated_object==Result)
                    {
                        double distance = 0;
                        TConfigure::Distance(rThisObject,*i_object,distance); // squared distance function
                        *Result   = *i_object;
                        Result++;
                        *Distances = distance;
                        Distances++;
                        NumberOfResults++;
                    }
                }
            }
        }
    }

    LocalIteratorType Begin()
    {
        return mObjects.begin();
    }

    LocalIteratorType End()
    {
        return mObjects.end();
    }

    SizeType Size()
    {
        return mObjects.size();
    }

    LocalIteratorType Begin() const
    {
        return mObjects.begin();
    }

    LocalIteratorType End() const
    {
        return mObjects.end();
    }

    SizeType Size() const
    {
        return mObjects.size();
    }

    PointerType GetObject(std::size_t Index)
    {
        return mObjects[Index];
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
        return "Cell Class ";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        return;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        return;
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

    std::vector<PointerType> mObjects;

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

}; // Class Cell

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class  TConfigure>
inline std::istream& operator >> (std::istream& rIStream,
                                  Cell<TConfigure>& rThis)
{
    return rIStream;
}

/// output stream function
template< class  TConfigure>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Cell<TConfigure>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.


