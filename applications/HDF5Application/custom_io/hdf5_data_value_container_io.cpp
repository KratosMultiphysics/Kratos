//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "includes/constitutive_law.h"
#include "utilities/data_type_traits.h"

// Application includes

// Include base h
#include "custom_io/hdf5_data_value_container_io.h"

namespace Kratos
{
namespace HDF5Utilities
{

template<class TDataType>
bool ReadComponent(
    HDF5::File& rFile,
    const std::string& rVariableName,
    const std::string& rPrefix,
    DataValueContainer& rData)
{
    KRATOS_TRY

    if (KratosComponents<Variable<TDataType>>::Has(rVariableName)) {
        const auto& r_variable = KratosComponents<Variable<TDataType>>::Get(rVariableName);
        if constexpr(std::is_same_v<TDataType, ConstitutiveLaw::Pointer>) {
            // reading the constitutive law name
            std::string constitutive_law_name;
            rFile.ReadAttribute(rPrefix + "/DataValues", rVariableName, constitutive_law_name);
            KRATOS_ERROR_IF_NOT(KratosComponents<ConstitutiveLaw>::Has(constitutive_law_name)) << "Kratos components missing \"" << constitutive_law_name << "\"" << std::endl;

            // TODO: The constitutive law cannot be constructed because, the constitutive law does not retain
            //       the construction mechanism it used in creating it (especially in composites, where the composite
            //       the structure of composite is given in json).
            // auto p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get(constitutive_law_name).Create(cl_parameters, rProperty);
            // rData[r_variable] = p_constitutive_law;
        } else {
            TDataType value{};
            rFile.ReadAttribute(rPrefix + "/DataValues", rVariableName, value);
            rData[r_variable] = value;
        }
        return true;
    } else {
        return false;
    }

    KRATOS_CATCH("");
}

template<class TDataType>
bool WriteComponent(
    HDF5::File& rFile,
    const std::string& rVariableName,
    const std::string& rPrefix,
    const DataValueContainer& rData)
{
    KRATOS_TRY

    if (KratosComponents<Variable<TDataType>>::Has(rVariableName)) {
        const auto& r_variable = KratosComponents<Variable<TDataType>>::Get(rVariableName);
        const auto& r_value = rData[r_variable];
        if constexpr(std::is_same_v<TDataType, ConstitutiveLaw::Pointer>) {
            auto components_cl = KratosComponents<ConstitutiveLaw>::GetComponents();
            std::string cl_name = "";
            for (const auto& comp_cl : components_cl) {
                if (r_value->HasSameType(r_value.get(), comp_cl.second)) {
                    cl_name = comp_cl.first;
                    break;
                }
            }
            KRATOS_ERROR_IF(cl_name == "") << "Kratos components missing \"" << r_value << "\"" << std::endl;
            rFile.WriteAttribute(rPrefix + "/DataValues", rVariableName, cl_name);
        } else {
            rFile.WriteAttribute(rPrefix + "/DataValues", rVariableName, r_value);
        }
        return true;
    } else {
        return false;
    }

    KRATOS_CATCH("");
}

template<class... TDataTypes>
void Read(
    HDF5::File& rFile,
    const std::string& rVariableName,
    const std::string& rPrefix,
    DataValueContainer& rData)
{
    KRATOS_TRY

    const bool is_read = (... || ReadComponent<TDataTypes>(rFile, rVariableName, rPrefix, rData));

    KRATOS_ERROR_IF_NOT(is_read) << "The variable \"" << rVariableName << "\" not found in registered variables list.";

    KRATOS_CATCH("");
}

template<class... TDataTypes>
void Write(
    HDF5::File& rFile,
    const std::string& rVariableName,
    const std::string& rPrefix,
    const DataValueContainer& rData)
{
    KRATOS_TRY

    const bool is_written = (... || WriteComponent<TDataTypes>(rFile, rVariableName, rPrefix, rData));

    KRATOS_ERROR_IF_NOT(is_written) << "The variable \"" << rVariableName << "\" not found in registered variables list.";

    KRATOS_CATCH("");
}

} // namespace HDF5Utilities

namespace HDF5
{
namespace Internals
{

void ReadDataValueContainer(
    File& rFile,
    const std::string& rPrefix,
    DataValueContainer& rData)
{
    KRATOS_TRY;

    const auto& attr_names = rFile.GetAttributeNames(rPrefix + "/DataValues");

    for (const auto& r_name : attr_names) {
        HDF5Utilities  ::Read<
            ConstitutiveLaw::Pointer,
            bool,
            int,
            double,
            std::string,
            array_1d<double, 3>,
            array_1d<double, 4>,
            array_1d<double, 6>,
            array_1d<double, 9>,
            Kratos::Vector,
            Kratos::Matrix>(rFile, r_name, rPrefix, rData);
    }

    KRATOS_CATCH("Path: \"" + rPrefix + "/DataValues\".");
}

void WriteDataValueContainer(
    File& rFile,
    const std::string& rPrefix,
    const DataValueContainer& rData)
{
    KRATOS_TRY;

    rFile.AddPath(rPrefix + "/DataValues");

    for (auto it = rData.begin(); it != rData.end(); ++it) {
        HDF5Utilities  ::Write<
            ConstitutiveLaw::Pointer,
            bool,
            int,
            double,
            std::string,
            array_1d<double, 3>,
            array_1d<double, 4>,
            array_1d<double, 6>,
            array_1d<double, 9>,
            Kratos::Vector,
            Kratos::Matrix>(rFile, it->first->Name(), rPrefix, rData);
    }

    KRATOS_CATCH("Path: \"" + rPrefix + "/DataValues\".");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.