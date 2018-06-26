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


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/variable_data.h"
#include "input_output/logger.h"
#include "includes/fnv_1a_hash.h"


namespace Kratos
{

    /// Constructor.
	VariableData::VariableData(const std::string& NewName, std::size_t NewSize, bool Iscomponent, char ComponentIndex) : mName(NewName), mKey(0), mSize(NewSize), mIsComponent(Iscomponent) {
        mKey = GenerateKey(mName, mSize, Iscomponent, ComponentIndex);
    }

	/// Copy constructor
    VariableData::VariableData(const VariableData& rOtherVariable)
        : mName(rOtherVariable.mName), mKey(rOtherVariable.mKey), mSize(rOtherVariable.mSize), mIsComponent(rOtherVariable.mIsComponent) {}

	VariableData::KeyType VariableData::GenerateKey(const std::string& Name, std::size_t Size, bool IsComponent, char ComponentIndex)
	{
        std::uint64_t key = Size;
        key <<= 32;
        key += FNV1a32Hash::CalculateHash(Name.c_str());

        key <<= 1;
        key += IsComponent;
        key <<= 7;
        key += ComponentIndex;

		return key;

	}

    void* VariableData::Clone(const void* pSource) const
    {
        return 0;
    }

    void* VariableData::Copy(const void* pSource, void* pDestination) const
    {
        return 0;
    }

    void VariableData::Assign(const void* pSource, void* pDestination) const {}

    void VariableData::AssignZero(void* pDestination) const {}

    void VariableData::Destruct(void* pSource) const {}

    void VariableData::Delete(void* pSource) const {}

    void VariableData::Print(const void* pSource, std::ostream& rOStream) const {}

    void VariableData::PrintData(const void* pSource, std::ostream& rOStream) const {}

    void VariableData::Allocate(void** pData) const
    {
    }

    void VariableData::Save(Serializer& rSerializer, void* pData) const
    {
    }

    void VariableData::Load(Serializer& rSerializer, void* pData) const
    {
    }

    /// NOTE: This function is for internal use and not
    /// to change arbitrary any variable's key
    void VariableData::SetKey(KeyType NewKey)
    {
        mKey = NewKey;
    }


    /// Turn back information as a string.
    std::string VariableData::Info() const
    {
        std::stringstream buffer;
        buffer << mName << " variable data" <<" #" << static_cast<unsigned int>(mKey);
        return buffer.str();
    }

    /// Print information about this object.
    void VariableData::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << mName << " variable data";
    }

    /// Print object's data.
    void VariableData::PrintData(std::ostream& rOStream) const
    {
        rOStream <<" #" << static_cast<unsigned int>(mKey);
    }

    void VariableData::save(Serializer& rSerializer) const
    {
        rSerializer.save("Name",mName);
        rSerializer.save("Key",mKey);
		rSerializer.save("IsComponent", mIsComponent);
    }

    void VariableData::load(Serializer& rSerializer)
    {
        rSerializer.load("Name",mName);
        rSerializer.load("Key",mKey);
		rSerializer.load("IsComponent", mIsComponent);
    }

  
}  // namespace Kratos.


