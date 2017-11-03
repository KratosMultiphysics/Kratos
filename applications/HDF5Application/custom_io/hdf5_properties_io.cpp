#include "custom_io/hdf5_properties_io.h"

#include <sstream>
#include "includes/kratos_components.h"

namespace Kratos
{
namespace HDF5
{
PropertiesIO::PropertiesIO(std::string Prefix, File::Pointer pFile)
    : mPrefix(Prefix), mpFile(pFile)
{
}

void PropertiesIO::ReadProperties(PropertiesContainerType& rProperties)
{
    KRATOS_TRY;

    Vector<int> prop_ids;
    mpFile->ReadAttribute(mPrefix + "/Properties", "Ids", prop_ids);
    std::vector<std::string> attr_names;

    for (unsigned i = 0; i < prop_ids.size(); ++i)
    {
        std::stringstream pstream;
        pstream << mPrefix << "/Properties(" << prop_ids[i] << ")";
        std::string path = pstream.str();
        mpFile->GetAttributeNames(path, attr_names);
        PropertiesType r_properties = rProperties[i];

        for (const auto& r_name : attr_names)
        {
            if (KratosComponents<Variable<int>>::Has(r_name))
            {
                int value;
                mpFile->ReadAttribute(path, r_name, value);
                const Variable<int>& rVARIABLE =
                    KratosComponents<Variable<int>>::Get(r_name);
                r_properties[rVARIABLE] = value;
            }
            else if (KratosComponents<Variable<double>>::Has(r_name))
            {
                double value;
                mpFile->ReadAttribute(path, r_name, value);
                const Variable<double>& rVARIABLE =
                    KratosComponents<Variable<double>>::Get(r_name);
                r_properties[rVARIABLE] = value;
            }
            else if (KratosComponents<Variable<Vector<double>>>::Has(r_name))
            {
                Vector<double> value;
                mpFile->ReadAttribute(path, r_name, value);
                const Variable<Vector<double>>& rVARIABLE =
                    KratosComponents<Variable<Vector<double>>>::Get(r_name);
                r_properties[rVARIABLE] = value;
            }
            else if (KratosComponents<Variable<Matrix<double>>>::Has(r_name))
            {
                Matrix<double> value;
                mpFile->ReadAttribute(path, r_name, value);
                const Variable<Matrix<double>>& rVARIABLE =
                    KratosComponents<Variable<Matrix<double>>>::Get(r_name);
                r_properties[rVARIABLE] = value;
            }
            else
            {
                KRATOS_ERROR << "Unsupported variable type: " << r_name << std::endl;
            }
        }
    }

    KRATOS_CATCH("");
}

void PropertiesIO::WriteProperties(Properties const& rProperties)
{
    KRATOS_TRY;

    std::stringstream pstream;
    pstream << mPrefix << "/Properties(" << rProperties.Id() << ")";
    std::string path = pstream.str();
    mpFile->AddPath(path);

    const PropertiesType::ContainerType& r_data = rProperties.Data();
    for (auto it = r_data.begin(); it != r_data.end(); ++it)
    {
        const std::string& r_name = it->first->Name();
        if (KratosComponents<Variable<int>>::Has(r_name))
        {
            const Variable<int>& rVARIABLE = KratosComponents<Variable<int>>::Get(r_name);
            mpFile->WriteAttribute(path, r_name, r_data[rVARIABLE]);
        }
        else if (KratosComponents<Variable<double>>::Has(r_name))
        {
            const Variable<double>& rVARIABLE = KratosComponents<Variable<double>>::Get(r_name);
            mpFile->WriteAttribute(path, r_name, r_data[rVARIABLE]);
        }
        else if (KratosComponents<Variable<Vector<double>>>::Has(r_name))
        {
            const Variable<Vector<double>>& rVARIABLE = KratosComponents<Variable<Vector<double>>>::Get(r_name);
            mpFile->WriteAttribute(path, r_name, r_data[rVARIABLE]);
        }
        else if (KratosComponents<Variable<Matrix<double>>>::Has(r_name))
        {
            const Variable<Matrix<double>>& rVARIABLE = KratosComponents<Variable<Matrix<double>>>::Get(r_name);
            mpFile->WriteAttribute(path, r_name, r_data[rVARIABLE]);
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable type: " << r_name << std::endl;
        }
    }

    KRATOS_CATCH("");
}

void PropertiesIO::WriteProperties(PropertiesContainerType const& rProperties)
{
    KRATOS_TRY;

    Vector<int> prop_ids(rProperties.size());
    unsigned i = 0;
    for (const PropertiesType& r_properties : rProperties)
    {
        prop_ids[i++] = r_properties.Id();
        WriteProperties(r_properties);
    }
    mpFile->WriteAttribute(mPrefix + "/Properties", "Ids", prop_ids);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
