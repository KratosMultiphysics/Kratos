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


#if !defined(KRATOS_SERIALIZER_H_INCLUDED )
#define  KRATOS_SERIALIZER_H_INCLUDED



// System includes
#include <string>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <fstream>



// External includes


// Project includes
#include "includes/define.h"
#include "containers/buffer.h"
#include "containers/weak_pointer_vector.h"
// #include "containers/variable.h"


#define KRATOS_SERIALIZATION_DIRECT_LOAD(type)                           \
    void load(std::string const & rTag, type& rValue)                \
    {                                                                \
        load_trace_point(rTag);                                      \
            read(rValue);                                             \
    }                                \
                                     \
    void load_base(std::string const & rTag, type& rValue)           \
    {                                                                \
          load_trace_point(rTag);                                      \
      read(rValue);                                             \
    }

#define KRATOS_SERIALIZATION_DIRECT_SAVE(type)                           \
    void save(std::string const & rTag, type const & rValue)         \
    {                                                                \
      save_trace_point(rTag);                                      \
      write(rValue);                                             \
    }                                    \
                                     \
    void save_base(std::string const & rTag, type const & rValue)    \
    {                                                                \
      save_trace_point(rTag);                                      \
      write(rValue);                                             \
    }

#define KRATOS_SERIALIZATION_DIRECT_CREATE(type)                         \
    void* create(std::string const & rTag, type* prototype)          \
    {                                                                \
      type* p_new = new type;                                        \
      load(rTag, *p_new);                                            \
      return p_new;                                                  \
    }

#define KRATOS_SERIALIZER_MODE_BINARY \
    if(!mTrace) {
#define KRATOS_SERIALIZER_MODE_ASCII \
    } else {
#define KRATOS_SERIALIZER_MODE_END \
    }
namespace Kratos
{

class ModelPart;
class VariableData;
template <class TDataType> class Variable;
//  template <class TDataType> class KratosComponents;


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
class KRATOS_API(KRATOS_CORE) Serializer
{
public:
    ///@name  Enum's
    ///@{

    enum PointerType {SP_INVALID_POINTER, SP_BASE_CLASS_POINTER, SP_DERIVED_CLASS_POINTER};
    enum TraceType {SERIALIZER_NO_TRACE=0, SERIALIZER_TRACE_ERROR=1, SERIALIZER_TRACE_ALL=2};

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Serializer
    KRATOS_CLASS_POINTER_DEFINITION(Serializer);

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
    Serializer(TraceType const& rTrace=SERIALIZER_NO_TRACE) : mpBuffer(new std::stringstream(std::ios::binary|std::ios::in|std::ios::out)), mTrace(rTrace), mNumberOfLines(0)
    {
    }

    Serializer(std::string const& Filename, TraceType const& rTrace=SERIALIZER_NO_TRACE) : mTrace(rTrace), mNumberOfLines(0)
    {
        std::fstream* p_file = new std::fstream(std::string(Filename+".rest").c_str(), std::ios::binary|std::ios::in|std::ios::out);
        if(!(*p_file))
        {
            delete p_file;
            p_file = new std::fstream(std::string(Filename+".rest").c_str(), std::ios::binary|std::ios::out);
        }
        mpBuffer = p_file;
        if(!(*mpBuffer))
            KRATOS_THROW_ERROR(std::invalid_argument, "Error opening input file : ", std::string(Filename+".rest"));
    }

    /// Destructor.
    virtual ~Serializer()
    {
        delete mpBuffer;
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
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
	//std::cout<<" REGISTERED OBJECT "<<rName<<" TypeID "<<typeid (TDataType).name()<<std::endl;
	//msRegisteredObjects.insert(RegisteredObjectsContainerType::value_type(rName,&pPrototype));
    }

    template<class TDataType>
    void load(std::string const & rTag, TDataType& rObject)
    {
        load_trace_point(rTag);
        rObject.load(*this);
    }

    template<class TDataType>
    void load(std::string const & rTag, boost::shared_ptr<TDataType>& pValue)
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
                    if(!pValue)
                        pValue = boost::shared_ptr<TDataType>(new TDataType);

                    load(rTag, *pValue);
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    if(i_prototype == msRegisteredObjects.end())
                        KRATOS_THROW_ERROR(std::runtime_error, "There is no object registered in Kratos with name : ", object_name)

                        if(!pValue)
                            pValue = boost::shared_ptr<TDataType>(static_cast<TDataType*>((i_prototype->second)()));

                    load(rTag, *pValue);

                }
                mLoadedPointers[p_pointer]=&pValue;
            }
            else
                pValue = *static_cast<boost::shared_ptr<TDataType>*>((i_pointer->second));
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
                    if(!pValue)
                        pValue = new TDataType;

                    load(rTag, *pValue);
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    if(i_prototype == msRegisteredObjects.end())
                        KRATOS_THROW_ERROR(std::runtime_error, "There is no object registered in Kratos with name : ", object_name)

                        if(!pValue)
                            pValue = static_cast<TDataType*>((i_prototype->second)());

                    load(rTag, *pValue);

                }
                mLoadedPointers[p_pointer]=&pValue;
            }
            else
            {
                pValue = *static_cast<TDataType**>((i_pointer->second));
            }
        }
    }

    template<class TDataType>
    void load(std::string const & rTag, boost::weak_ptr<TDataType>& pValue)
    {
        // This is for testing. I have to change it. Pooyan.
        //KRATOS_THROW_ERROR(std::logic_error, "The serialization for weak_ptrs is not implemented yet", "")
//    read(*pValue);
    }

    template<class TDataType>
    void load(std::string const & rTag, WeakPointerVector<TDataType>& pValue)
    {
        // This is for testing. I have to change it. Pooyan.
        //KRATOS_THROW_ERROR(std::logic_error, "The serialization for weak_ptrs is not implemented yet", "")
//    read(*pValue);
    }

    template<class TDataType>
    void load(std::string const & rTag, const Variable<TDataType>* pVariable)
    {
        load_trace_point(rTag);
        std::string name;
        read(name);

        pVariable = static_cast<const Variable<TDataType>*>(GetVariableData(name));
    }

    template<class TDataType>
    void load(std::string const & rTag, std::vector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        SizeType size;

        load("size", size);

        rObject.resize(size);

        for(SizeType i = 0 ; i < size ; i++)
            load("E", rObject[i]);
//    read(rObject);
    }

    template<class TDataType>
    void load(std::string const & rTag, boost::numeric::ublas::vector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        SizeType size;

        load("size", size);

        rObject.resize(size,false);

        for(SizeType i = 0 ; i < size ; i++)
            load("E", rObject[i]);
//    read(rObject);
    }

    template<class TDataType, std::size_t TDimension>
    void load(std::string const & rTag, array_1d<TDataType, TDimension>& rObject)
    {
        load_trace_point(rTag);
	//rObject = array_1d<TDataType, TDimension>(); //it generates a warnning --> commented 23/09/2015 <--
        for(SizeType i = 0 ; i < TDimension ; i++)
            load("E", rObject[i]);
//    read(rObject);
    }

    KRATOS_SERIALIZATION_DIRECT_LOAD(bool)
    KRATOS_SERIALIZATION_DIRECT_LOAD(int)
    KRATOS_SERIALIZATION_DIRECT_LOAD(long)
    KRATOS_SERIALIZATION_DIRECT_LOAD(double)
    KRATOS_SERIALIZATION_DIRECT_LOAD(unsigned long)
    KRATOS_SERIALIZATION_DIRECT_LOAD(unsigned int)
    KRATOS_SERIALIZATION_DIRECT_LOAD(std::string)
    KRATOS_SERIALIZATION_DIRECT_LOAD(Matrix)
    KRATOS_SERIALIZATION_DIRECT_LOAD(long long)
//#ifdef  _WIN32 // work around for windows int64_t error
//    KRATOS_SERIALIZATION_DIRECT_LOAD(__int64)
//#endif
#ifdef  _WIN64 // work around for windows size_t error in win64
    KRATOS_SERIALIZATION_DIRECT_LOAD(std::size_t)
#endif



    template<class TDataType>
    void save(std::string const & rTag, std::vector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        SizeType size = rObject.size();

        save("size", size);

        for(SizeType i = 0 ; i < size ; i++)
            save("E", rObject[i]);
//    write(rObject);
    }

    template<class TDataType>
    void save(std::string const & rTag, boost::numeric::ublas::vector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        SizeType size = rObject.size();

        save("size", size);

        for(SizeType i = 0 ; i < size ; i++)
            save("E", rObject[i]);
//    write(rObject);
    }

    template<class TDataType, std::size_t TDimension>
    void save(std::string const & rTag, array_1d<TDataType, TDimension> const& rObject)
    {
        save_trace_point(rTag);
        for(SizeType i = 0 ; i < TDimension ; i++)
            save("E", rObject[i]);

//    write(rObject);
    }

    template<class TDataType>
    void save(std::string const & rTag, TDataType const& rObject)
    {
        save_trace_point(rTag);
        rObject.save(*this);
    }

    template<class TDataType>
    void save(std::string const & rTag, const Variable<TDataType>* pVariable)
    {
        save_trace_point(rTag);
        write(pVariable->Name());
    }


    template<class TDataType>
    void save(std::string const & rTag, boost::shared_ptr<TDataType> pValue)
    {
        save(rTag, pValue.get());
    }

    template<class TDataType>
    void save(std::string const & rTag, const TDataType * pValue)
    {
        if(pValue)
        {
            if(IsDerived(pValue))
                write(SP_DERIVED_CLASS_POINTER);
            else
                write(SP_BASE_CLASS_POINTER);

            SavePointer(rTag,pValue);
        }
        else
        {
            write(SP_INVALID_POINTER);
        }
    }

    template<class TDataType>
    bool IsDerived(TDataType * pValue)
    {
      if (strcmp(typeid(TDataType).name(), typeid(*pValue).name()) != 0) {
	return true;
      }
      else {
	return false;
      }
      // bool is_derived = (typeid(TDataType) != typeid(*pValue));
//    std::cout << "for TDataType : " << typeid(TDataType).name() << " and *pValue type : " << typeid(*pValue).name() << " is derived : " << is_derived << std::endl;
      //return is_derived;
    }


    template<class TDataType>
    void save(std::string const & rTag, TDataType * pValue)
    {
        if(pValue)
        {
            if(IsDerived(pValue))
            {
                write(SP_DERIVED_CLASS_POINTER);
            }
            else
            {
                write(SP_BASE_CLASS_POINTER);
            }

            SavePointer(rTag,pValue);
        }
        else
        {
            write(SP_INVALID_POINTER);
        }
    }

    template<class TDataType>
    void save(std::string const & rTag, boost::weak_ptr<TDataType> pValue)
    {
        // This is for testing. I have to implement it. Pooyan.
        //KRATOS_THROW_ERROR(std::logic_error, "The serialization for weak_ptrs is not implemented yet", "")
//    write(*pValue);
    }

    template<class TDataType>
    void save(std::string const & rTag, Kratos::WeakPointerVector<TDataType> pValue)
    {
        // This is for testing. I have to implement it. Pooyan.
        //KRATOS_THROW_ERROR(std::logic_error, "The serialization for weak_ptrs is not implemented yet", "")
//    write(*pValue);
    }

    template<class TDataType>
    void save(std::string const & rTag, boost::shared_ptr<const TDataType> pValue)
    {
        // This is for testing. I have to change it. Pooyan.
//          save_trace_point(rTag);
//    write(*pValue);
        save(rTag, pValue.get());
    }

    void save(std::string const & rTag, const char * pValue)
    {
        save_trace_point(rTag);
        write(std::string(pValue));
    }

    KRATOS_SERIALIZATION_DIRECT_SAVE(bool)
    KRATOS_SERIALIZATION_DIRECT_SAVE(int)
    KRATOS_SERIALIZATION_DIRECT_SAVE(long)
    KRATOS_SERIALIZATION_DIRECT_SAVE(double)
    KRATOS_SERIALIZATION_DIRECT_SAVE(unsigned long)
    KRATOS_SERIALIZATION_DIRECT_SAVE(unsigned int)
    KRATOS_SERIALIZATION_DIRECT_SAVE(std::string)
//        KRATOS_SERIALIZATION_DIRECT_SAVE(Vector)
    KRATOS_SERIALIZATION_DIRECT_SAVE(Matrix)
    KRATOS_SERIALIZATION_DIRECT_SAVE(long long)
//#ifdef  _WIN32 // work around for windows int64_t error
//    KRATOS_SERIALIZATION_DIRECT_SAVE(__int64)
//#endif
#ifdef  _WIN64 // work around for windows size_t error in win64
    KRATOS_SERIALIZATION_DIRECT_SAVE(std::size_t)
#endif


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
    void load_base(std::string const & rTag, boost::numeric::ublas::vector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        load(rTag, rObject);
    }

    template<class TDataType, std::size_t TDimension>
    void load_base(std::string const & rTag, array_1d<TDataType, TDimension>& rObject)
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
    void save_base(std::string const & rTag, boost::numeric::ublas::vector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        save(rTag, rObject);
    }

    template<class TDataType, std::size_t TDimension>
    void save_base(std::string const & rTag, array_1d<TDataType, TDimension> const& rObject)
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
        if(mTrace)
        {
            write(rTag);
        }
    }

    bool load_trace_point(std::string const & rTag)
    {
        if(mTrace == SERIALIZER_TRACE_ERROR) // only reporting the errors
        {
            std::string read_tag;
            read(read_tag);
            if(read_tag == rTag)
                return true;
            else
            {
                std::stringstream buffer;
                buffer << "In line " << mNumberOfLines;
                buffer << " the trace tag is not the expected one:" << std::endl;
                buffer << "    Tag found : " << read_tag << std::endl;
                buffer << "    Tag given : " << rTag << std::endl;
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            }
        }
        else if(mTrace == SERIALIZER_TRACE_ALL) // also reporting matched tags.
        {
            std::string read_tag;
            read(read_tag);
            if(read_tag == rTag)
            {
                std::cout << "In line " << mNumberOfLines;
                std::cout << " loading " << rTag << " as expected" << std::endl;
                return true;
            }
            else
            {
                std::stringstream buffer;
                buffer << "In line " << mNumberOfLines;
                buffer << " the trace tag is not the expected one:" << std::endl;
                buffer << "    Tag found : " << read_tag << std::endl;
                buffer << "    Tag given : " << rTag << std::endl;
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
    ///@name Inquiry
    ///@{


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
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    template<class TDataType>
    void SavePointer(std::string const & rTag, const TDataType * pValue)
    {
        write(pValue);
        if (mSavedPointers.find(pValue) == mSavedPointers.end())
        {
            if (IsDerived(pValue))
            {
                typename RegisteredObjectsNameContainerType::iterator i_name = msRegisteredObjectsName.find(typeid (*pValue).name());

                if (i_name == msRegisteredObjectsName.end())
                    KRATOS_THROW_ERROR(std::runtime_error, "There is no object registered in Kratos with type id : ", typeid (*pValue).name())
                    else
                        write(i_name->second);


            }

            save(rTag, *pValue);
            //        pValue->save(*this);
            mSavedPointers.insert(pValue);
        }
    }

    VariableData* GetVariableData(std::string const & VariableName);

//        void read(bool& rData)
//        {
//            int temp;
//            read(temp);
//            rData = temp << 1;
//
//         }
//
//        void write(bool& rData)
//        {
//            int temp(rData);
//            write(temp);
//        }
//
//        void read(std::string& rValue)
//        {
//            *mpBuffer >> rValue;
//        }
//
//        void write(std::string& rValue)
//        {
//            *mpBuffer << rValue;
//        }

    void read(PointerType& rValue)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        int temp;
        mpBuffer->read((char *)(&temp),sizeof(PointerType));
        rValue = PointerType(temp);

        KRATOS_SERIALIZER_MODE_ASCII

        int temp;
        *mpBuffer >> temp;
        rValue = PointerType(temp);
        mNumberOfLines++;

        KRATOS_SERIALIZER_MODE_END
    }

    void write(PointerType const& rValue)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        int ptr = (int)rValue;
        const char * data = reinterpret_cast<const char*>(&ptr);
        mpBuffer->write(data,sizeof(PointerType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << int(rValue) << std::endl;

        KRATOS_SERIALIZER_MODE_END
    }

    void read(std::string& rValue)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType size;
        mpBuffer->read((char *)(&size),sizeof(SizeType));
        char* c_binStream = new char [size];
        mpBuffer->read(c_binStream,size);
        std::string s_binStream(c_binStream,size);
        rValue = s_binStream;
        delete [] c_binStream;

        KRATOS_SERIALIZER_MODE_ASCII

        // going to the first '"'
        std::getline( *mpBuffer,rValue, '\"');
        // reading the string itself until second '"'
        std::getline( *mpBuffer,rValue, '\"');
        mNumberOfLines++;

        KRATOS_SERIALIZER_MODE_END
    }

    void write(std::string const& rValue)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        const char * data = rValue.c_str();
        SizeType rData_size = rValue.length() * sizeof(char);

        const char * data1 = reinterpret_cast<const char *>(&rData_size);

        mpBuffer->write(data1,sizeof(SizeType));
        mpBuffer->write(data,rData_size);

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << "\"" << rValue << "\"" << std::endl;

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void read(TDataType& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        mpBuffer->read((char *)(&rData),sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer >> rData;
        mNumberOfLines++;

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void write(TDataType const& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        const char * data = reinterpret_cast<const char*>(&rData);
        mpBuffer->write(data,sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << rData << std::endl;

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void read(std::vector<TDataType>& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType size;
        mpBuffer->read((char *)(&size),sizeof(SizeType));

        rData.resize(size);

        read(rData.begin(), rData.end(), sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        std::size_t size;
        *mpBuffer >> size;
        rData.resize(size);
        mNumberOfLines++;

        read(rData.begin(), rData.end());

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void write(std::vector<TDataType> const& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType rData_size = rData.size();

        const char * data = reinterpret_cast<const char *>(&rData_size);
        mpBuffer->write(data,sizeof(SizeType));

        write(rData.begin(), rData.end(), sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << rData.size() << std::endl;
        write(rData.begin(), rData.end());

        KRATOS_SERIALIZER_MODE_END
    }

//        template<class TDataType, std::size_t TDimenasion>
//        void read(array_1d<TDataType, TDimenasion>& rData)
//        {
//            read(rData.begin(), rData.end());
//        }
//
//        template<class TDataType, std::size_t TDimension>
//        void write(array_1d<TDataType, TDimension> const& rData)
//        {
//           write(rData.begin(), rData.end());
//        }
//
//        template<class TDataType>
//        void read(boost::numeric::ublas::vector<TDataType>& rData)
//        {
//            std::size_t size;
//            *mpBuffer >> size;
//            rData.resize(size,false);
//
//            read(rData.begin(), rData.end());
//        }
//
//        template<class TDataType>
//        void write(boost::numeric::ublas::vector<TDataType> const& rData)
//        {
//            *mpBuffer << rData.size() << std::endl;
//            write(rData.begin(), rData.end());
//        }

    template<class TDataType>
    void read(boost::numeric::ublas::matrix<TDataType>& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType size1;
        SizeType size2;

        mpBuffer->read((char *)(&size1),sizeof(SizeType));
        mpBuffer->read((char *)(&size2),sizeof(SizeType));

        rData.resize(size1,size2);

        read(rData.data().begin(), rData.data().end(), sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        SizeType size1;
        SizeType size2;

        *mpBuffer >> size1;
        mNumberOfLines++;
        *mpBuffer >> size2;
        mNumberOfLines++;

        rData.resize(size1,size2);

        read(rData.data().begin(), rData.data().end(),0);

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void write(boost::numeric::ublas::matrix<TDataType> const& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType rData_size1 = rData.size1();
        SizeType rData_size2 = rData.size2();

        const char * data1 = reinterpret_cast<const char *>(&rData_size1);
        const char * data2 = reinterpret_cast<const char *>(&rData_size2);

        mpBuffer->write(data1,sizeof(SizeType));
        mpBuffer->write(data2,sizeof(SizeType));

        write(rData.data().begin(), rData.data().end(), sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << rData.size1() << std::endl;
        *mpBuffer << rData.size2() << std::endl;

        write(rData.data().begin(), rData.data().end(),0);

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TIteratorType>
    void read(TIteratorType First, TIteratorType Last, SizeType size)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        for(; First != Last ; First++)
        {
            mpBuffer->read((char *)First,sizeof(size));
        }

        KRATOS_SERIALIZER_MODE_ASCII

        for(; First != Last ; First++)
        {
            *mpBuffer >> *First;
            mNumberOfLines++;

        }

        KRATOS_SERIALIZER_MODE_END
    }
    template<class TIteratorType>
    void write(TIteratorType First, TIteratorType Last, SizeType size)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        for(; First != Last ; First++)
        {
            const char * data = reinterpret_cast<const char *>(First);
            mpBuffer->write(data,sizeof(size));
        }

        KRATOS_SERIALIZER_MODE_ASCII

        for(; First != Last ; First++)
            *mpBuffer << *First << std::endl;

        KRATOS_SERIALIZER_MODE_END
    }

    inline SizeType BlockCompatibleSize(SizeType rSize)
    {
        typedef char BlockType;
        const SizeType block_size = sizeof(BlockType);
        return static_cast<SizeType>(((block_size - 1) + rSize) / block_size);
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

    /// Assignment operator.
    Serializer& operator=(Serializer const& rOther);

    /// Copy constructor.
    Serializer(Serializer const& rOther);


    ///@}

}; // Class Serializer

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


//   template<class TDataType>
//   inline Serializer& operator >> (Serializer& rThis, TDataType& rObject)
//   {
//     rThis.load(rObject);

//     return rThis;
//   }


//   template<class TDataType>
//   inline Serializer& operator << (Serializer& rThis, TDataType& rObject)
//   {
//     rThis.save(rObject, KRATOS_VERSION);

//     return rThis;
//   }
/// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
//                  Serializer& rThis);

/// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
//                  const Serializer& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
///@}


}  // namespace Kratos.

#undef KRATOS_SERIALIZER_MODE_BINARY
#undef KRATOS_SERIALIZER_MODE_ASCII
#undef KRATOS_SERIALIZER_MODE_END

#undef KRATOS_SERIALIZATION_DIRECT_LOAD
#undef KRATOS_SERIALIZATION_DIRECT_SAVE

#endif // KRATOS_SERIALIZER_H_INCLUDED  defined
