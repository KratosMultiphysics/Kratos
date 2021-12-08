//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_SERIALIZER_INCLUDED
#define CO_SIM_IO_SERIALIZER_INCLUDED

// System includes
#include <string>
#include <cstring>
#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <sstream>
#include <fstream>
#include <memory>
#include <array>
#include <vector>
#include <utility>

// Project includes
#include "define.hpp"

#define CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(type)                    \
    void load(std::string const & rTag, type& rValue)                \
    {                                                                \
        load_trace_point(rTag);                                      \
        read(rValue);                                                \
    }                                                                \
    void load(std::string const & rTag, type const& rValue)          \
    {                                                                \
        load_trace_point(rTag);                                      \
        read(const_cast<type&>(rValue));                             \
    }                                                                \
    void load_base(std::string const & rTag, type& rValue)           \
    {                                                                \
        load_trace_point(rTag);                                      \
        read(rValue);                                                \
    }

#define CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(type)                   \
    void save(std::string const & rTag, type const & rValue)        \
    {                                                               \
      save_trace_point(rTag);                                       \
      write(rValue);                                                \
    }                                                               \
    void save_base(std::string const & rTag, type const & rValue)   \
    {                                                               \
      save_trace_point(rTag);                                       \
      write(rValue);                                                \
    }

#define CO_SIM_IO_SERIALIZATION_DIRECT_CREATE(type)                 \
    void* create(std::string const & rTag, type* prototype)         \
    {                                                               \
      type* p_new = new type;                                       \
      load(rTag, *p_new);                                           \
      return p_new;                                                 \
    }

#define CO_SIM_IO_SERIALIZER_MODE_BINARY \
    if(!mTrace) {
#define CO_SIM_IO_SERIALIZER_MODE_ASCII \
    } else {
#define CO_SIM_IO_SERIALIZER_MODE_END \
    }

namespace CoSimIO {
namespace Internals {

///@name CoSimIO Classes
///@{

/**
 * @class Serializer
 *
 * @brief The serialization consists in storing the state of an object into a storage format like data file or memory buffer and also retrieving the object from such a media.
 *
 * @details The serialization consists in storing the state of an object into a storage format like data file or memory buffer and also retrieving the object from such a media.
 * The idea of serialization is based on saving all object's data consecutively in the file or buffer and then load it in the same order.
 *
 * @author Pooyan Dadvand
 */
class CO_SIM_IO_API Serializer
{
public:
    ///@name  Enum's
    ///@{

    enum PointerType {SP_INVALID_POINTER, SP_BASE_CLASS_POINTER, SP_DERIVED_CLASS_POINTER};
    enum TraceType {SERIALIZER_NO_TRACE=0, SERIALIZER_TRACE_ERROR=1, SERIALIZER_TRACE_ALL=2};

    ///@}
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;

    typedef void* (*ObjectFactoryType)();

    typedef std::map<void*, void*> LoadedPointersContainerType;

    typedef std::map<std::string, ObjectFactoryType> RegisteredObjectsContainerType;

    typedef std::map<std::string, std::string> RegisteredObjectsNameContainerType;

    typedef std::set<const void*> SavedPointersContainerType;

    typedef std::iostream BufferType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit Serializer(BufferType* pBuffer, TraceType const& rTrace=SERIALIZER_NO_TRACE) :
        mpBuffer(pBuffer), mTrace(rTrace), mNumberOfLines(0)
    {
    }

    /// Destructor.
    virtual ~Serializer()
    {
        delete mpBuffer;
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Serializer& operator=(Serializer const& rOther) = delete;

    /// Copy constructor.
    Serializer(Serializer const& rOther) = delete;

    /// Sets the Serializer in a state ready to be loaded
    /// Note: If the same object is loaded twice before deleting it from memory all its pointers will be duplicated.
    void SetLoadState();

    ///@}
    ///@name Operations
    ///@{
    /// This function returns the "trace type" used in initializing the serializer.
    /// Trace type is one of SERIALIZER_NO_TRACE,SERIALIZER_TRACE_ERROR,SERIALIZER_TRACE_ALL
    TraceType GetTraceType() const {return mTrace;}

    void SetBuffer(BufferType* pBuffer)
    {
        mpBuffer = pBuffer;
    }

    template<class TDataType>
    static void* Create()
    {
        return new TDataType;
    }

    template<class TDataType>
    static void Register(std::string const & rName, TDataType const& pPrototype)
    {
        msRegisteredObjects.insert(RegisteredObjectsContainerType::value_type(rName,Create<TDataType>));
        msRegisteredObjectsName.insert(RegisteredObjectsNameContainerType::value_type(typeid(TDataType).name(), rName));
    }

    template<class TDataType>
    void load(std::string const & rTag, TDataType& rObject)
    {
        load_trace_point(rTag);
        rObject.load(*this);
    }

    template<class TDataType>
    void load(std::string const & rTag, std::shared_ptr<TDataType>& pValue)
    {
        PointerType pointer_type = SP_INVALID_POINTER;
        void* p_pointer;
        read(pointer_type);

        if(pointer_type != SP_INVALID_POINTER)
        {
            read(p_pointer);
            LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
            if(i_pointer == mLoadedPointers.end())
            {
                if(pointer_type == SP_BASE_CLASS_POINTER)
                {
                    if(!pValue) {
                        pValue = std::shared_ptr<TDataType>(new TDataType);
                    }
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    CO_SIM_IO_ERROR_IF(i_prototype == msRegisteredObjects.end())
                        << "There is no object registered in CoSimIO with name : "
                        << object_name << std::endl;

                    if(!pValue) {
                        pValue = std::shared_ptr<TDataType>(static_cast<TDataType*>((i_prototype->second)()));
                    }
                }

                // Load the pointer address before loading the content
                mLoadedPointers[p_pointer]=&pValue;
                load(rTag, *pValue);
            }
            else
            {
                pValue = *static_cast<std::shared_ptr<TDataType>*>((i_pointer->second));
            }
        }
    }

    template<class TDataType>
    void load(std::string const & rTag, CoSimIO::intrusive_ptr<TDataType>& pValue)
    {
        PointerType pointer_type = SP_INVALID_POINTER;
        void* p_pointer;
        read(pointer_type);

        if(pointer_type != SP_INVALID_POINTER)
        {
            read(p_pointer);
            LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
            if(i_pointer == mLoadedPointers.end())
            {
                if(pointer_type == SP_BASE_CLASS_POINTER)
                {
                    if(!pValue) {
                        pValue = CoSimIO::intrusive_ptr<TDataType>(new TDataType);
                    }
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    CO_SIM_IO_ERROR_IF(i_prototype == msRegisteredObjects.end())
                        << "There is no object registered in CoSimIO with name : "
                        << object_name << std::endl;

                    if(!pValue) {
                        pValue = CoSimIO::intrusive_ptr<TDataType>(static_cast<TDataType*>((i_prototype->second)()));
                    }
                }

                // Load the pointer address before loading the content
                mLoadedPointers[p_pointer]=&pValue;
                load(rTag, *pValue);
            }
            else
            {
                pValue = *static_cast<CoSimIO::intrusive_ptr<TDataType>*>((i_pointer->second));
            }
        }
    }

    template<class TDataType>
    void load(std::string const & rTag, std::unique_ptr<TDataType>& pValue)
    {
        PointerType pointer_type = SP_INVALID_POINTER;
        void* p_pointer;
        read(pointer_type);

        if(pointer_type != SP_INVALID_POINTER)
        {
            read(p_pointer);
            LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
            if(i_pointer == mLoadedPointers.end())
            {
                if(pointer_type == SP_BASE_CLASS_POINTER)
                {
                    if(!pValue) {
                        pValue = std::unique_ptr<TDataType>(new TDataType);
                    }
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    CO_SIM_IO_ERROR_IF(i_prototype == msRegisteredObjects.end())
                        << "There is no object registered in CoSimIO with name : "
                        << object_name << std::endl;

                    if(!pValue) {
                        pValue = std::move(std::unique_ptr<TDataType>(static_cast<TDataType*>((i_prototype->second)())));
                    }
                }

                // Load the pointer address before loading the content
                mLoadedPointers[p_pointer]=pValue.get();
                load(rTag, *pValue);
            }
            else
            {
                pValue = std::move(std::unique_ptr<TDataType>(static_cast<TDataType*>((i_pointer->second))));
            }
        }
    }

    template<class TDataType>
    void load(std::string const & rTag, TDataType*& pValue)
    {
        PointerType pointer_type = SP_INVALID_POINTER;
        void* p_pointer;
        read(pointer_type);

        if(pointer_type != SP_INVALID_POINTER)
        {
            read(p_pointer);
            LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
            if(i_pointer == mLoadedPointers.end())
            {
                if(pointer_type == SP_BASE_CLASS_POINTER)
                {
                    if(!pValue) {
                        pValue = new TDataType;
                    }
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    CO_SIM_IO_ERROR_IF(i_prototype == msRegisteredObjects.end())
                        << "There is no object registered in CoSimIO with name : "
                        << object_name << std::endl;

                    if(!pValue) {
                        pValue = static_cast<TDataType*>((i_prototype->second)());
                    }

                }

                // Load the pointer address before loading the content
                mLoadedPointers[p_pointer]=&pValue;
                load(rTag, *pValue);
            }
            else
            {
                pValue = *static_cast<TDataType**>((i_pointer->second));
            }
        }
    }

    template<class TDataType, std::size_t TDataSize>
    void load(std::string const & rTag, std::array<TDataType, TDataSize>& rObject)
    {
        load_trace_point(rTag);
        for (SizeType i = 0; i < TDataSize; i++)
            load("E", rObject[i]);
    }

    template<class TDataType>
    void load(std::string const & rTag, std::vector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        SizeType size;

        load("size", size);

        rObject.resize(size);

        for (SizeType i = 0 ; i < size ; i++)
            load("E", rObject[i]);
    }

    template<class TKeyType, class TDataType>
    void load(std::string const & rTag, std::map<TKeyType, TDataType>& rObject)
    {
        load_map(rTag, rObject);
    }

    template<class TKeyType, class TDataType>
    void load(std::string const & rTag, std::unordered_map<TKeyType, TDataType>& rObject)
    {
        load_map(rTag, rObject);
    }


    template<class TFirstType, class TSecondType>
    void load(std::string const & rTag, std::pair<TFirstType, TSecondType>& rObject)
    {
        load_trace_point(rTag);
        load("First", rObject.first);
        load("Second", rObject.second);
    }

    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(bool)
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(int)
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(long)
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(double)
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(float)
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(unsigned long)
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(unsigned int)
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(std::string)
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(long long)
//#ifdef  _WIN32 // work around for windows int64_t error
//    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(__int64)
//#endif
#ifdef  _WIN64 // work around for windows size_t error in win64
    CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(std::size_t)
#endif
    // CO_SIM_IO_SERIALIZATION_DIRECT_LOAD(std::complex<double>)

    template<class TDataType, std::size_t TDataSize>
    void save(std::string const & rTag, std::array<TDataType, TDataSize> const& rObject)
    {
        save_trace_point(rTag);
        for (SizeType i = 0; i < TDataSize; i++)
            save("E", rObject[i]);
    }

    template<class TDataType>
    void save(std::string const & rTag, std::vector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        SizeType size = rObject.size();

        save("size", size);

        for (SizeType i = 0 ; i < size ; i++)
            save("E", rObject[i]);
    }

    template<class TKeyType, class TDataType>
    void save(std::string const & rTag, std::map<TKeyType, TDataType> const& rObject)
    {
        save_map(rTag, rObject);
    }

    template<class TKeyType, class TDataType>
    void save(std::string const & rTag, std::unordered_map<TKeyType, TDataType> const& rObject)
    {
        save_map(rTag, rObject);
    }

    template<class TDataType>
    void save(std::string const & rTag, TDataType const& rObject)
    {
        save_trace_point(rTag);
        rObject.save(*this);
    }

    template<class TDataType>
    void save(std::string const & rTag, std::shared_ptr<TDataType> pValue)
    {
        save(rTag, pValue.get());
    }

    template<class TDataType>
    void save(std::string const & rTag, CoSimIO::intrusive_ptr<TDataType> pValue)
    {
        save(rTag, pValue.get());
    }


    template<class TDataType>
    void save(std::string const & rTag, std::unique_ptr<TDataType> const& pValue)
    {
        save(rTag, pValue.get());
    }

    template<class TDataType>
    void save(std::string const & rTag, const TDataType * pValue)
    {
        if (pValue) {
            if(IsDerived(pValue)) {
                write(SP_DERIVED_CLASS_POINTER);
            } else {
                write(SP_BASE_CLASS_POINTER);
            }

            SavePointer(rTag,pValue);
        } else {
            write(SP_INVALID_POINTER);
        }
    }

    template<class TDataType>
    void save(std::string const & rTag, TDataType * pValue)
    {
        if (pValue) {
            if (IsDerived(pValue)) {
                write(SP_DERIVED_CLASS_POINTER);
            } else {
                write(SP_BASE_CLASS_POINTER);
            }

            SavePointer(rTag,pValue);
        } else {
            write(SP_INVALID_POINTER);
        }
    }

    template<class TDataType>
    bool IsDerived(const TDataType * pValue)
    {
        return typeid(TDataType).name() != typeid(*pValue).name();
    }

    void save(std::string const & rTag, const char * pValue)
    {
        save_trace_point(rTag);
        write(std::string(pValue));
    }

    template<class TFirstType, class TSecondType>
    void save(std::string const & rTag, const std::pair<TFirstType, TSecondType>& rObject)
    {
        save_trace_point(rTag);
        save("First", rObject.first);
        save("Second", rObject.second);
    }

    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(bool)
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(int)
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(long)
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(double)
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(float)
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(unsigned long)
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(unsigned int)
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(std::string)
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(long long)
//#ifdef  _WIN32 // work around for windows int64_t error
//    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(__int64)
//#endif
#ifdef  _WIN64 // work around for windows size_t error in win64
    CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(std::size_t)
#endif
    // CO_SIM_IO_SERIALIZATION_DIRECT_SAVE(std::complex<double>)


    template<class TDataType>
    void load_base(std::string const & rTag, TDataType& rObject)
    {
        load_trace_point(rTag);
        rObject.TDataType::load(*this);
    }

    template<class TDataType>
    void load_base(std::string const & rTag, std::vector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        load(rTag, rObject);
    }

    template<class TDataType>
    void save_base(std::string const & rTag, std::vector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        save(rTag, rObject);
    }

    template<class TDataType>
    void save_base(std::string const & rTag, TDataType const& rObject)
    {
        save_trace_point(rTag);
        rObject.TDataType::save(*this);
    }

    void save_trace_point(std::string const & rTag)
    {
        if(mTrace) {
            write(rTag);
        }
    }

    bool load_trace_point(std::string const & rTag)
    {
        if(mTrace == SERIALIZER_TRACE_ERROR) {// only reporting the errors
            std::string read_tag;
            read(read_tag);
            if(read_tag == rTag) {
                return true;
            } else {
                std::stringstream buffer;
                buffer << "In line " << mNumberOfLines;
                buffer << " the trace tag is not the expected one:" << std::endl;
                buffer << "    Tag found : " << read_tag << std::endl;
                buffer << "    Tag given : " << rTag << std::endl;
                CO_SIM_IO_ERROR << buffer.str() << std::endl;
            }
        } else if (mTrace == SERIALIZER_TRACE_ALL) {// also reporting matched tags.
            std::string read_tag;
            read(read_tag);
            if(read_tag == rTag) {
                CO_SIM_IO_INFO("Serializer") << "In line " << mNumberOfLines << " loading " << rTag << " as expected" << std::endl;
                return true;
            } else {
                std::stringstream buffer;
                buffer << "In line " << mNumberOfLines;
                buffer << " the trace tag is not the expected one:" << std::endl;
                buffer << "    Tag found : " << read_tag << std::endl;
                buffer << "    Tag given : " << rTag << std::endl;
                CO_SIM_IO_ERROR << buffer.str() << std::endl;
            }
        }
        return false;

    }

    ///@}
    ///@name Access
    ///@{

    BufferType* pGetBuffer()
    {
        return mpBuffer;
    }

    static RegisteredObjectsContainerType& GetRegisteredObjects()
    {
        return msRegisteredObjects;
    }

    static RegisteredObjectsNameContainerType& GetRegisteredObjectsName()
    {
        return msRegisteredObjectsName;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Serializer";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static RegisteredObjectsContainerType msRegisteredObjects;
    static RegisteredObjectsNameContainerType msRegisteredObjectsName;

    ///@}
    ///@name Member Variables
    ///@{

    BufferType* mpBuffer;
    TraceType mTrace;
    SizeType mNumberOfLines;

    SavedPointersContainerType mSavedPointers;
    LoadedPointersContainerType mLoadedPointers;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TDataType>
    void SavePointer(std::string const & rTag, const TDataType * pValue)
    {
        write(pValue);
        if (mSavedPointers.find(pValue) == mSavedPointers.end()) {
            mSavedPointers.insert(pValue);
            if (IsDerived(pValue)) {
                typename RegisteredObjectsNameContainerType::iterator i_name = msRegisteredObjectsName.find(typeid (*pValue).name());

                if (i_name == msRegisteredObjectsName.end()) {
                    CO_SIM_IO_ERROR << "There is no object registered in CoSimIO with type id : "
                                 << typeid (*pValue).name() << std::endl;
                } else {
                    write(i_name->second);
                }
            }

            save(rTag, *pValue);
        }
    }

    template<class TMapType>
    void load_map(std::string const & rTag, TMapType& rObject)
    {
        load_trace_point(rTag);
        SizeType size = rObject.size();

        load("size", size);

        for (SizeType i = 0 ; i < size ; i++) {
            typename TMapType::value_type temp;
            load("E", temp);
            rObject.insert(std::move(temp));
        }
    }

    template<class TMapType>
    void save_map(std::string const & rTag, TMapType const& rObject)
    {
        save_trace_point(rTag);
        SizeType size = rObject.size();

        save("size", size);

        for (auto& i : rObject) {
            save("E", i);
        }
    }

    void read(PointerType& rValue)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        int temp;
        mpBuffer->read((char *)(&temp),sizeof(PointerType));
        rValue = PointerType(temp);

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        int temp;
        *mpBuffer >> temp;
        rValue = PointerType(temp);
        mNumberOfLines++;

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    void write(PointerType const& rValue)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        int ptr = (int)rValue;
        const char * data = reinterpret_cast<const char*>(&ptr);
        mpBuffer->write(data,sizeof(PointerType));

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        *mpBuffer << int(rValue) << std::endl;

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    void read(std::string& rValue)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        SizeType size;
        mpBuffer->read((char *)(&size),sizeof(SizeType));
        char* c_binStream = new char [size];
        mpBuffer->read(c_binStream,size);
        std::string s_binStream(c_binStream,size);
        rValue = s_binStream;
        delete [] c_binStream;

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        // going to the first '"'
        std::getline( *mpBuffer,rValue, '\"');
        // reading the string itself until second '"'
        std::getline( *mpBuffer,rValue, '\"');
        mNumberOfLines++;

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    void write(std::string const& rValue)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        const char * data = rValue.c_str();
        SizeType rData_size = rValue.length() * sizeof(char);

        const char * data1 = reinterpret_cast<const char *>(&rData_size);

        mpBuffer->write(data1,sizeof(SizeType));
        mpBuffer->write(data,rData_size);

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        *mpBuffer << "\"" << rValue << "\"" << std::endl;

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void read(TDataType& rData)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        mpBuffer->read((char *)(&rData),sizeof(TDataType));

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        *mpBuffer >> rData;
        mNumberOfLines++;

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void write(TDataType const& rData)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        const char * data = reinterpret_cast<const char*>(&rData);
        mpBuffer->write(data,sizeof(TDataType));

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        *mpBuffer << rData << std::endl;

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void read(std::vector<TDataType>& rData)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        SizeType size;
        mpBuffer->read((char *)(&size),sizeof(SizeType));

        rData.resize(size);

        read(rData.begin(), rData.end(), sizeof(TDataType));

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        std::size_t size;
        *mpBuffer >> size;
        rData.resize(size);
        mNumberOfLines++;

        read(rData.begin(), rData.end());

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void write(std::vector<TDataType> const& rData)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        SizeType rData_size = rData.size();

        const char * data = reinterpret_cast<const char *>(&rData_size);
        mpBuffer->write(data,sizeof(SizeType));

        write(rData.begin(), rData.end(), sizeof(TDataType));

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        *mpBuffer << rData.size() << std::endl;
        write(rData.begin(), rData.end());

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    template<class TIteratorType>
    void read(TIteratorType First, TIteratorType Last, SizeType size)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        for (; First != Last ; First++) {
            mpBuffer->read((char *)First,sizeof(size));
        }

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        for (; First != Last ; First++) {
            *mpBuffer >> *First;
            mNumberOfLines++;
        }

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    template<class TIteratorType>
    void write(TIteratorType First, TIteratorType Last, SizeType size)
    {
        CO_SIM_IO_SERIALIZER_MODE_BINARY

        for (; First != Last ; First++) {
            const char * data = reinterpret_cast<const char *>(First);
            mpBuffer->write(data,sizeof(size));
        }

        CO_SIM_IO_SERIALIZER_MODE_ASCII

        for (; First != Last ; First++) {
            *mpBuffer << *First << std::endl;
        }

        CO_SIM_IO_SERIALIZER_MODE_END
    }

    ///@}

}; // Class Serializer

///@}

#define CO_SIM_IO_SERIALIZE_SAVE_BASE_CLASS(Serializer, BaseType) \
    Serializer.save_base("BaseClass",*static_cast<const BaseType *>(this));

#define CO_SIM_IO_SERIALIZE_LOAD_BASE_CLASS(Serializer, BaseType) \
    Serializer.load_base("BaseClass",*static_cast<BaseType *>(this));

} // namespace Internals
} // namespace CoSimIO

#undef CO_SIM_IO_SERIALIZER_MODE_BINARY
#undef CO_SIM_IO_SERIALIZER_MODE_ASCII
#undef CO_SIM_IO_SERIALIZER_MODE_END

#undef CO_SIM_IO_SERIALIZATION_DIRECT_LOAD
#undef CO_SIM_IO_SERIALIZATION_DIRECT_SAVE

#endif // CO_SIM_IO_SERIALIZER_INCLUDED
