#include "custom_io/hdf5_data_value_container_io.h"

#include "custom_utilities/registered_variable_lookup.h"

namespace Kratos
{

template <typename TVariable>
class ReadVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    HDF5::File& rFile,
                    std::string const& rPrefix,
                    DataValueContainer& rData)
    {
        typename TVariable::Type value;
        rFile.ReadAttribute(rPrefix + "/DataValues", rVariable.Name(), value);
        rData[rVariable] = value;
    }
};

template <typename TVariable>
class WriteVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    HDF5::File& rFile,
                    std::string const& rPrefix,
                    DataValueContainer const& rData)
    {
        rFile.WriteAttribute(rPrefix + "/DataValues", rVariable.Name(), rData[rVariable]);
    }
};

namespace HDF5
{
namespace Internals
{

void ReadDataValueContainer(File& rFile, std::string const& rPrefix, DataValueContainer& rData)
{
    KRATOS_TRY;
    std::vector<std::string> attr_names = rFile.GetAttributeNames(rPrefix + "/DataValues");
    for (const auto& r_name : attr_names)
        RegisteredVariableLookup<Variable<int>, Variable<double>,
                                 Variable<Vector<double>>, Variable<Matrix<double>>>(r_name)
            .Execute<ReadVariableFunctor>(rFile, rPrefix, rData);
    KRATOS_CATCH("Path: \"" + rPrefix + "/DataValues\".");
}

void WriteDataValueContainer(File& rFile, std::string const& rPrefix, DataValueContainer const& rData)
{
    KRATOS_TRY;
    rFile.AddPath(rPrefix + "/DataValues");
    for (auto it = rData.begin(); it != rData.end(); ++it)
        try
        {
            RegisteredVariableLookup<Variable<int>, Variable<double>, Variable<Vector<double>>,
                                     Variable<Matrix<double>>>(it->first->Name())
                .Execute<WriteVariableFunctor>(rFile, rPrefix, rData);
        }
        catch (Exception& e)
        {
        }
    KRATOS_CATCH("Path: \"" + rPrefix + "/DataValues\".");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
