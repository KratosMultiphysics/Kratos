// schema.cpp

#include "includes/schema.h"

namespace Kratos
{

std::unordered_map<std::string, std::string_view> Schema::mByName;
std::unordered_map<std::string, std::vector<std::string_view>> Schema::mByType;


void Schema::Register(const std::string& name, const std::string& type, std::string_view value)
{
    auto& by_name = GetByName();
    auto& by_type = GetByType();

    const auto it = by_name.find(name);

    if (it != by_name.end()) {
        // Same registration coming from another DSO.
        if (it->second == value) {
            return;
        }

        // Same name but different schema: real conflict.
        KRATOS_ERROR
            << "Schema '" << name
            << "' is already registered with a different value."
            << std::endl;
    }

    by_name.emplace(name, value);
    by_type[type].push_back(value);
}


std::string_view Schema::Get(const std::string& name)
{
    auto& by_name = Schema::GetByName();
    const auto it = by_name.find(name);

    if (it == by_name.end()) {
        KRATOS_ERROR
            << "Registry item with name "
            << name << " does not exist."
            << std::endl;
    }

    return it->second;
}

const std::vector<std::string_view>& Schema::FindInstances(const std::string& type)
{
    auto& by_type = Schema::GetByType();
    const auto it = by_type.find(type);

    if (it == by_type.end()) {
        KRATOS_ERROR << "Registry item with type "
                     << type << " does not exist." << std::endl;
    }

    return it->second;
}

Schema::ByNameMap& Schema::GetByName()
{
    static ByNameMap mByName;
    return mByName;
}

Schema::ByTypeMap& Schema::GetByType()
{
    static ByTypeMap mByType;
    return mByType;
}

}