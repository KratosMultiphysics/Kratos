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


#if !defined(KRATOS_DATA_VALUE_CONTAINER_H_INCLUDED )
#define KRATOS_DATA_VALUE_CONTAINER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <vector>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "includes/kratos_components.h"
#include "includes/exception.h"

#ifdef KRATOS_DEBUG
#include "utilities/openmp_utils.h"
#endif


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
class KRATOS_API(KRATOS_CORE) DataValueContainer
{

public:
    ///@name Type Definitions
    ///@{
    KRATOS_DEFINE_LOCAL_FLAG(OVERWRITE_OLD_VALUES);

    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(DataValueContainer);

    /// Type of the container used for variables
    typedef std::pair<const VariableData*, void*> ValueType;

    /// Type of the container used for variables
    typedef std::vector<ValueType> ContainerType;

    /// Type of the container used for variables
    typedef std::vector<ValueType>::iterator iterator;

    /// Type of the container used for variables
    typedef std::vector<ValueType>::const_iterator const_iterator;

    /// Type of the container used for variables
    typedef std::vector<ValueType>::size_type SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataValueContainer() {}

    /// Copy constructor.
    DataValueContainer(DataValueContainer const& rOther)
    {
        for(const_iterator i = rOther.mData.begin() ; i != rOther.mData.end() ; ++i)
            mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
    }

    /// Destructor.
    virtual ~DataValueContainer()
    {
        for(iterator i = mData.begin() ; i != mData.end() ; ++i)
            i->first->Delete(i->second);
    }


    ///@}
    ///@name Operators
    ///@{

    template<class TDataType> const TDataType& operator()(const VariableData& rThisVariable) const
    {
        return GetValue<TDataType>(rThisVariable);
    }

    template<class TDataType> TDataType& operator()(const Variable<TDataType>& rThisVariable)
    {
        return GetValue<TDataType>(rThisVariable);
    }

    template<class TDataType> const TDataType& operator()(const Variable<TDataType>& rThisVariable) const
    {
        return GetValue<TDataType>(rThisVariable);
    }

    template<class TDataType> TDataType& operator[](const VariableData& rThisVariable)
    {
        return GetValue<TDataType>(rThisVariable);
    }

    template<class TDataType> const TDataType& operator[](const VariableData& rThisVariable) const
    {
        return GetValue<TDataType>(rThisVariable);
    }

    template<class TDataType> TDataType& operator[](const Variable<TDataType>& rThisVariable)
    {
        return GetValue<TDataType>(rThisVariable);
    }

    template<class TDataType> const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
    {
        return GetValue<TDataType>(rThisVariable);
    }

    iterator begin()
    {
        return mData.begin();
    }

    const_iterator begin() const
    {
        return mData.begin();
    }

    iterator end()
    {
        return mData.end();
    }

    const_iterator end() const
    {
        return mData.end();
    }

    /// Assignment operator.
    DataValueContainer& operator=(const DataValueContainer& rOther)
    {
        Clear();

        for(const_iterator i = rOther.mData.begin() ; i != rOther.mData.end() ; ++i)
            mData.push_back(ValueType(i->first, i->first->Clone(i->second)));

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    template<class TDataType> TDataType& GetValue(const Variable<TDataType>& rThisVariable)
    {
        typename ContainerType::iterator i;

        if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey())))  != mData.end())
            return *(static_cast<TDataType*>(i->second) + rThisVariable.GetComponentIndex());

#ifdef KRATOS_DEBUG
        if(OpenMPUtils::IsInParallel() != 0)
            KRATOS_ERROR << "attempting to do a GetValue for: " << rThisVariable << " unfortunately the variable is not in the database and the operations is not threadsafe (this function is being called from within a parallel region)" << std::endl;
#endif

        auto p_source_variable = &rThisVariable.GetSourceVariable();
        mData.push_back(ValueType(p_source_variable,p_source_variable->Clone(p_source_variable->pZero())));

        return *(static_cast<TDataType*>(mData.back().second) + rThisVariable.GetComponentIndex());
    }

    //TODO: make the variable of the constant version consistent with the one of the "classical" one
    template<class TDataType> const TDataType& GetValue(const Variable<TDataType>& rThisVariable) const
    {
        typename ContainerType::const_iterator i;

        if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey())))  != mData.end())
            return *(static_cast<const TDataType*>(i->second) + rThisVariable.GetComponentIndex());

        return rThisVariable.Zero();
    }

    SizeType Size()
    {
        return mData.size();
    }

    template<class TDataType> void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue)
    {
        typename ContainerType::iterator i;

        if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey())))  != mData.end()) {
            *(static_cast<TDataType*>(i->second) + rThisVariable.GetComponentIndex()) = rValue;
        } else {
            auto p_source_variable = &rThisVariable.GetSourceVariable();
            mData.push_back(ValueType(p_source_variable,p_source_variable->Clone(p_source_variable->pZero())));
            *(static_cast<TDataType*>(mData.back().second) + rThisVariable.GetComponentIndex()) = rValue;
        }
    }

    template<class TDataType> void Erase(const Variable<TDataType>& rThisVariable)
    {
        typename ContainerType::iterator i;

        if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey())))  != mData.end())
        {
            i->first->Delete(i->second);
            mData.erase(i);
        }
    }

    void Clear()
    {
        for(ContainerType::iterator i = mData.begin() ; i != mData.end() ; i++)
            i->first->Delete(i->second);

        mData.clear();
    }

    void Merge(const DataValueContainer& rOther, Flags Options);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return (std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey())) != mData.end());
    }

    bool IsEmpty() const
    {
        return mData.empty();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return std::string("data value container");
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "data value container";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        for(const_iterator i = mData.begin() ; i != mData.end() ; ++i)
        {
            rOStream <<"    ";
            i->first->Print(i->second, rOStream);
            rOStream << std::endl;
        }
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

    class IndexCheck
    {
        std::size_t mI;
    public:
        explicit IndexCheck(std::size_t I) : mI(I) {}
        bool operator()(const ValueType& I)
        {
            return I.first->SourceKey() == mI;
        }
    };

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ContainerType mData;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


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

}; // Class DataValueContainer


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DataValueContainer& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DataValueContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_DATA_VALUE_CONTAINER_H_INCLUDED  defined


