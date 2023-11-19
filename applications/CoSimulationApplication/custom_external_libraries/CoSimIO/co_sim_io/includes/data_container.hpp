//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_DATA_CONTAINER_INCLUDED
#define CO_SIM_IO_DATA_CONTAINER_INCLUDED

// System includes
#include <vector>
#include <algorithm> // std::max
#include <ostream>

// Project includes
#include "define.hpp"
#include "serializer.hpp"

namespace CoSimIO {
namespace Internals {

template<typename TDataType>
class DataContainer
{
public:
    DataContainer() = default;
    virtual ~DataContainer() = default;

    DataContainer(const DataContainer&) = delete;
    DataContainer& operator=(const DataContainer&) = delete;

    virtual std::size_t size() const = 0;
    virtual void resize(const std::size_t NewSize) = 0;

    virtual TDataType* data() = 0;
    virtual const TDataType* data() const = 0;
    const TDataType* data_const() const
    {
        const auto& r_const_this = *this;
        return this->data();
    }

    const TDataType& operator[](const std::size_t Index) const
    {
        return this->data()[Index];
    }

    TDataType& operator[](const std::size_t Index)
    {
        return this->data()[Index];
    }

private:
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        CO_SIM_IO_TRY

        rSerializer.save("size", size());
        for (std::size_t i=0; i<size(); ++i) {
            rSerializer.save("v", data()[i]);
        }

        CO_SIM_IO_CATCH
    }

    virtual void load(Serializer& rSerializer)
    {
        CO_SIM_IO_TRY

        std::size_t new_size;
        rSerializer.load("size", new_size);
        if (size() != new_size) {
            resize(new_size);
        }

        for (std::size_t i=0; i<size(); ++i) {
            rSerializer.load("v", data()[i]);
        }

        CO_SIM_IO_CATCH
    }
};

/// output stream function
template<typename TDataType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DataContainer<TDataType>& rThis)
{
    const std::size_t size = rThis.size();

    rOStream << "[";
    if(size>0) rOStream << rThis[0];
    if(size>1) {
        for(std::size_t i=1; i<size; ++i)
            rOStream<<", "<<rThis[i];
    }
    rOStream << "]";

    return rOStream;
}

template<typename TDataType>
class DataContainerStdVector : public DataContainer<TDataType>
{
public:
    explicit DataContainerStdVector(std::vector<TDataType>& rVector)
        : mrVector(rVector) {}

    std::size_t size() const override {return mrVector.size();}
    void resize(const std::size_t NewSize) override {mrVector.resize(NewSize);} // resize does not change the capacity if resized to a smaller size
    const TDataType* data() const override {return mrVector.data();}
    TDataType* data() override {return mrVector.data();}

private:
    std::vector<TDataType>& mrVector;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        CO_SIM_IO_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DataContainer<TDataType>)
    }

    void load(Serializer& rSerializer) override
    {
        CO_SIM_IO_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DataContainer<TDataType>)
    }
};

template<typename TDataType>
class DataContainerStdVectorReadOnly : public DataContainer<TDataType>
{
public:
    explicit DataContainerStdVectorReadOnly(const std::vector<TDataType>& rVector)
        : mrVector(rVector) {}

    std::size_t size() const override {return mrVector.size();}
    void resize(const std::size_t NewSize) override {CO_SIM_IO_ERROR << "Resizing of readonly object is not possible!" << std::endl;}
    const TDataType* data() const override {return mrVector.data();}
    TDataType* data() override {CO_SIM_IO_ERROR << "Using non-const member of readonly object!" << std::endl; return nullptr;}

private:
    const std::vector<TDataType>& mrVector;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        CO_SIM_IO_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DataContainer<TDataType>)
    }

    void load(Serializer& rSerializer) override {CO_SIM_IO_ERROR << "Loading a readonly object is not possible!" << std::endl;}
};

template<typename TDataType>
class DataContainerRawMemory : public DataContainer<TDataType>
{
public:
    explicit DataContainerRawMemory(TDataType** ppData, const std::size_t Size)
        : mppData(ppData), mSize(Size), mCapacity(Size) {}

    std::size_t size() const override {return mSize;};
    void resize(const std::size_t NewSize) override
    {
        mSize = NewSize;
        if (NewSize > mCapacity) { // only increase the capacity if too small => same behavior as std::vector
            if (mCapacity == 0) {
                *mppData = nullptr; // initial allocation if not done outside
            }
            mCapacity = std::max(mCapacity*2, NewSize); // increase size beyond what is necessary to avoid frequent reallocations (like std::vector does)

            *mppData = (TDataType *)realloc(*mppData, (mCapacity)*sizeof(TDataType));

            CO_SIM_IO_ERROR_IF_NOT(*mppData) << "Memory reallocation failed!";
        }
    };

    const TDataType* data() const override {return *mppData;}
    TDataType* data() override {return *mppData;}

private:
    TDataType** mppData;
    std::size_t mSize;
    std::size_t mCapacity;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        CO_SIM_IO_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DataContainer<TDataType>)
    }

    void load(Serializer& rSerializer) override
    {
        CO_SIM_IO_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DataContainer<TDataType>)
    }
};

template<typename TDataType>
class DataContainerRawMemoryReadOnly : public DataContainer<TDataType>
{
public:
    explicit DataContainerRawMemoryReadOnly(const TDataType* pData, const std::size_t Size)
        : mpData(pData), mSize(Size) {}

    std::size_t size() const override {return mSize;};
    void resize(const std::size_t NewSize) override {CO_SIM_IO_ERROR << "Resizing of readonly object is not possible!" << std::endl;};
    const TDataType* data() const override {return mpData;}
    TDataType* data() override {CO_SIM_IO_ERROR << "Using non-const member of readonly object!" << std::endl; return nullptr;}

private:
    const TDataType* mpData;
    const std::size_t mSize;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        CO_SIM_IO_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DataContainer<TDataType>)
    }

    void load(Serializer& rSerializer) override {CO_SIM_IO_ERROR << "Loading a readonly object is not possible!" << std::endl;}
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_DATA_CONTAINER_INCLUDED
