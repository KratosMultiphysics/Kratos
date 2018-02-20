#include "custom_io/hdf5_data_value_container_io.h"

#include "includes/kratos_components.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

void ReadDataValueContainer(File& rFile, std::string const& rPrefix, DataValueContainer& rData)
{
    KRATOS_TRY;

    std::vector<std::string> attr_names;
    rFile.GetAttributeNames(rPrefix + "/DataValues", attr_names);

    for (const auto& r_name : attr_names)
    {
        if (KratosComponents<Variable<int>>::Has(r_name))
        {
            int value;
            rFile.ReadAttribute(rPrefix + "/DataValues", r_name, value);
            const Variable<int>& rVARIABLE = KratosComponents<Variable<int>>::Get(r_name);
            rData[rVARIABLE] = value;
        }
        else if (KratosComponents<Variable<double>>::Has(r_name))
        {
            double value;
            rFile.ReadAttribute(rPrefix + "/DataValues", r_name, value);
            const Variable<double>& rVARIABLE =
                KratosComponents<Variable<double>>::Get(r_name);
            rData[rVARIABLE] = value;
        }
        else if (KratosComponents<Variable<Vector<double>>>::Has(r_name))
        {
            Vector<double> value;
            rFile.ReadAttribute(rPrefix + "/DataValues", r_name, value);
            const Variable<Vector<double>>& rVARIABLE =
                KratosComponents<Variable<Vector<double>>>::Get(r_name);
            rData[rVARIABLE] = value;
        }
        else if (KratosComponents<Variable<Matrix<double>>>::Has(r_name))
        {
            Matrix<double> value;
            rFile.ReadAttribute(rPrefix + "/DataValues", r_name, value);
            const Variable<Matrix<double>>& rVARIABLE =
                KratosComponents<Variable<Matrix<double>>>::Get(r_name);
            rData[rVARIABLE] = value;
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable type: " << r_name << std::endl;
        }
    }

    KRATOS_CATCH("Path: \"" + rPrefix + "/DataValues\".");
}

void WriteDataValueContainer(File& rFile, std::string const& rPrefix, DataValueContainer const& rData)
{
    KRATOS_TRY;

    rFile.AddPath(rPrefix + "/DataValues");

    for (auto it = rData.begin(); it != rData.end(); ++it)
    {
        const std::string& r_name = it->first->Name();
        if (KratosComponents<Variable<int>>::Has(r_name))
        {
            const Variable<int>& rVARIABLE = KratosComponents<Variable<int>>::Get(r_name);
            rFile.WriteAttribute(rPrefix + "/DataValues", r_name, rData[rVARIABLE]);
        }
        else if (KratosComponents<Variable<double>>::Has(r_name))
        {
            const Variable<double>& rVARIABLE = KratosComponents<Variable<double>>::Get(r_name);
            rFile.WriteAttribute(rPrefix + "/DataValues", r_name, rData[rVARIABLE]);
        }
        else if (KratosComponents<Variable<Vector<double>>>::Has(r_name))
        {
            const Variable<Vector<double>>& rVARIABLE = KratosComponents<Variable<Vector<double>>>::Get(r_name);
            rFile.WriteAttribute(rPrefix + "/DataValues", r_name, rData[rVARIABLE]);
        }
        else if (KratosComponents<Variable<Matrix<double>>>::Has(r_name))
        {
            const Variable<Matrix<double>>& rVARIABLE = KratosComponents<Variable<Matrix<double>>>::Get(r_name);
            rFile.WriteAttribute(rPrefix + "/DataValues", r_name, rData[rVARIABLE]);
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable type: " << r_name << std::endl;
        }
    }

    KRATOS_CATCH("Path: \"" + rPrefix + "/DataValues\".");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
