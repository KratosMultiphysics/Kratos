//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
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
#include "utilities/logger.h"


namespace Kratos
{

    /// Constructor.
	VariableData::VariableData(const std::string& NewName, std::size_t NewSize, bool Iscomponent) : mName(NewName), mKey(0), mSize(NewSize), mIsComponent(Iscomponent) {}

	/// Copy constructor
    VariableData::VariableData(const VariableData& rOtherVariable)
        : mName(rOtherVariable.mName), mKey(rOtherVariable.mKey), mSize(rOtherVariable.mSize), mIsComponent(rOtherVariable.mIsComponent) {}

	VariableData::KeyType VariableData::GenerateKey(const std::string& Name, std::size_t Size, std::size_t ComponentIndex)
	{
		// For generation of 32bit hash key I use the 
		// Jenkins's one-at-a-time hash taken from wikipedia 
		// I could have used better hash functions like murmur or lookup3
		// But finally choose this for simplicity. Pooyan.
		// https://en.wikipedia.org/wiki/Jenkins_hash_function

		unsigned int hash, i;
		for(hash = i = 0; i < Name.size(); ++i)
		{
			hash += Name[i];
			hash += (hash << 10);
			hash ^= (hash >> 6);
		}
		hash += (hash << 3);
		hash ^= (hash >> 11);
		hash += (hash << 15);

		//if(Size > 127) // to be store in a char
		//	KRATOS_ERROR << "A variable of Kratos cannot be larger than 127 bytes and variable " << Name << " has sizeof " << Size;
//		KeyType key = hash;
//		key << 2; // This is for adding the 
		return hash;

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
        buffer << mName << " variable data";
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


