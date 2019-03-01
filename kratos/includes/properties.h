//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//


#if !defined(KRATOS_PROPERTIES_H_INCLUDED )
#define  KRATOS_PROPERTIES_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <unordered_map>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/data_value_container.h"
#include "includes/process_info.h"
#include "includes/table.h"


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

/// Properties encapsulates data shared by different Elements or Conditions
/**
 * Properties encapsulates data shared by different Elements or Conditions. It can store any type of data and provides a variable base access to them.
 *  These are all parameters that can be shared between Element. Usually material parameters are common for a set of element, so this category of data is referred as properties.
 * But in general it can be any common parameter for a group of Elements. Sharing these data as properties reduces the memory used by the application and also helps updating them if necessary.
 * As mentioned before Properties is a shared data container between Elements or Conditions. In finite element problems there are several parameters which are the same for a set of elements and conditions.
 * Thermal conductivity, elasticity of the material and viscosity of the fluid are examples of these parameters. Properties holds these data and is shared by elements or Conditions. This eliminates memory overhead due to redundant copies of these data for each element and Condition. Properties also can be used to access nodal data if it is necessary.
 * It is important to mention that accessing the nodal data via Properties is not the same as accessing it via Node. When user asks Properties for a variable data in a Node, the process starts with finding the variable in the Properties data container and if it does not exist then get it from Node.
 * This means that the priority of data is with the one stored in Properties and then in Node.
 */
class Properties : public IndexedObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Properties
    KRATOS_CLASS_POINTER_DEFINITION(Properties);

#ifdef  _WIN32 // work around for windows int64_t error
    typedef __int64 int64_t;
#endif
    typedef IndexedObject BaseType;

    typedef DataValueContainer ContainerType;

    typedef Node<3> NodeType;

    typedef NodeType::IndexType IndexType;

    typedef Table<double> TableType;

    typedef std::unordered_map<std::size_t, TableType> TablesContainerType; // This is a provisional implmentation and should be changed to hash. Pooyan.


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit Properties(IndexType NewId = 0) : BaseType(NewId), mData(), mTables() {}

    /// Copy constructor.
    Properties(const Properties& rOther) : BaseType(rOther), mData(rOther.mData), mTables(rOther.mTables) {}

    /// Destructor.
    ~Properties() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Properties& operator=(const Properties& rOther)
    {
        BaseType::operator=(rOther);
        mData = rOther.mData;
        return *this;
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rV)
    {
        return GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rV) const
    {
        return GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type& operator[](const TVariableType& rV)
    {
        return GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator[](const TVariableType& rV) const
    {
        return GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rV, NodeType& rThisNode)
    {
        return GetValue(rV, rThisNode);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rV, NodeType const& rThisNode) const
    {
        return GetValue(rV, rThisNode);
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rV, NodeType& rThisNode, IndexType SolutionStepIndex)
    {
        return GetValue(rV, rThisNode, SolutionStepIndex);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rV, NodeType const& rThisNode, IndexType SolutionStepIndex) const
    {
        return GetValue(rV, rThisNode, SolutionStepIndex);
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rV, NodeType& rThisNode, ProcessInfo const& rCurrentProcessInfo)
    {
        return GetValue(rV, rThisNode, rCurrentProcessInfo.GetSolutionStepIndex());
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rV, NodeType const& rThisNode, ProcessInfo const& rCurrentProcessInfo) const
    {
        return GetValue(rV, rThisNode, rCurrentProcessInfo.GetSolutionStepIndex());
    }

    ///@}
    ///@name Operations
    ///@{

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rV)
    {
        return mData.GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rV) const
    {
        return mData.GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rV, NodeType& rThisNode)
    {
        if(mData.Has(rV))
            return mData.GetValue(rV);

        return rThisNode.GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rV, NodeType const& rThisNode) const
    {
        if(mData.Has(rV))
            return mData.GetValue(rV);

        return rThisNode.GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rV, NodeType& rThisNode, IndexType SolutionStepIndex)
    {
        if(mData.Has(rV))
            return mData.GetValue(rV);

        return rThisNode.GetValue(rV, SolutionStepIndex);
    }

    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rV, NodeType const& rThisNode, IndexType SolutionStepIndex) const
    {
        if(mData.Has(rV))
            return mData.GetValue(rV);

        return rThisNode.GetValue(rV, SolutionStepIndex);
    }

    template<class TVariableType>
    void SetValue(TVariableType const& rV, typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rV, rValue);
    }

    bool HasVariables() const
    {
        return !mData.IsEmpty();
    }

    template<class TXVariableType, class TYVariableType>
    TableType& GetTable(const TXVariableType& XVariable, const TYVariableType& YVariable)
    {
        return mTables[Key(XVariable.Key(), YVariable.Key())];
    }

    template<class TXVariableType, class TYVariableType>
    TableType const& GetTable(const TXVariableType& XVariable, const TYVariableType& YVariable) const
    {
        return mTables.at(Key(XVariable.Key(), YVariable.Key()));
    }

    template<class TXVariableType, class TYVariableType>
    void SetTable(const TXVariableType& XVariable, const TYVariableType& YVariable, TableType const& rThisTable)
    {
        mTables[Key(XVariable.Key(), YVariable.Key())] = rThisTable;
    }

    bool HasTables() const
    {
        return !mTables.empty();
    }

    bool IsEmpty() const
    {
        return !( HasVariables() || HasTables() );
    }

    int64_t Key(std::size_t XKey, std::size_t YKey) const
    {
        int64_t result_key = XKey;
        result_key = result_key << 32;
        result_key |= YKey; // I know that the key is less than 2^32 so I don't need zeroing the upper part
        return result_key;
    }

    ///@}
    ///@name Access
    ///@{

    ContainerType& Data()
    {
        return mData;
    }

    ContainerType const& Data() const
    {
        return mData;
    }

    TablesContainerType& Tables()
    {
        return mTables;
    }

    TablesContainerType const& Tables() const
    {
        return mTables;
    }


    ///@}
    ///@name Inquiry
    ///@{

    template<class TVariableType>
    bool Has(TVariableType const& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    template<class TXVariableType, class TYVariableType>
    bool HasTable(const TXVariableType& XVariable, const TYVariableType& YVariable) const
    {
        return (mTables.find(Key(XVariable.Key(), YVariable.Key())) != mTables.end());
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Properties";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream <<  "Properties";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        mData.PrintData(rOStream);
        rOStream << "This properties contains " << mTables.size() << " tables";
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

    ContainerType mData;
    TablesContainerType mTables;

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
        rSerializer.save("Data", mData);
        rSerializer.save("Tables", mTables);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
        rSerializer.load("Data", mData);
        rSerializer.load("Tables", mTables);
    }

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

}; // Class Properties

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Properties& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Properties& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PROPERTIES_H_INCLUDED  defined


