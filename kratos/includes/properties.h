// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#if !defined(KRATOS_PROPERTIES_H_INCLUDED )
#define  KRATOS_PROPERTIES_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <map> // This is a provisional implmentation and should be changed to hash. Pooyan.


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/data_value_container.h"
//#include "containers/all_variables_data_value_container.h"
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

/// Short class definition.
/** Detail class definition.
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

	typedef std::map<int64_t, TableType> TablesContainerType; // This is a provisional implmentation and should be changed to hash. Pooyan.


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
	Properties(IndexType NewId = 0) : BaseType(NewId), mData(), mTables() {}

    /// Copy constructor.
    Properties(const Properties& rOther) : BaseType(rOther), mData(rOther.mData), mTables(rOther.mTables) {}

    /// Destructor.
    virtual ~Properties() {}


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
        mData.GetValue(rV) = rValue;
    }

    template<class XVariableType, class YVariableType>
    TableType& GetTable(const XVariableType& XVariable, const YVariableType& YVariable)
    {
		return mTables[Key(XVariable, YVariable)];
    }

    template<class XVariableType, class YVariableType>
    TableType const& GetTable(const XVariableType& XVariable, const YVariableType& YVariable) const
    {
		return mTables[Key(XVariable.Key(), YVariable.Key())];
    }

    template<class XVariableType, class YVariableType>
    void SetTable(const XVariableType& XVariable, const YVariableType& YVariable, TableType const& rThisTable)
    {
		mTables[Key(XVariable.Key(), YVariable.Key())] = rThisTable;
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



    ///@}
    ///@name Inquiry
    ///@{

    template<class TVariableType>
    bool Has(TVariableType const& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Properties";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream <<  "Properties";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
        rSerializer.save("Data", mData);
    }

    void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
        rSerializer.load("Data", mData);
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


