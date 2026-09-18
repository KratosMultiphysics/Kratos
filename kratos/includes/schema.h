#pragma once

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <memory>
#include <typeinfo>
#include <iostream>
#include <ranges>

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) Schema
{
public:
    Schema() = delete;

    using ByNameMap = std::unordered_map<std::string, std::string_view>;
    using ByTypeMap = std::unordered_map<std::string, std::vector<std::string_view>>;

    static void Register(const std::string& name, const std::string& type, std::string_view value);

    static std::string_view Get(const std::string& name);
    static ByNameMap& GetByName();
    static ByTypeMap& GetByType();

    static const std::vector<std::string_view>& FindInstances(const std::string& type);

private:
    static ByNameMap mByName;
    static ByTypeMap mByType;
};

template<class T>
class SchemaRegistrar
{
public:
    SchemaRegistrar(std::string_view name, std::string_view type, std::string_view schema)
    {
        Schema::Register(std::string(name), std::string(type), schema);
    }
};

}