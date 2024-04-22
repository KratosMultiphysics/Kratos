
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#pragma once

// System includes

// External includes

// Project includes
#include "factories/factory.h"
#include "includes/kratos_components.h"
#include "includes/kratos_parameters.h"

#include "modeler/voxel_mesh_generator_modeler.h"

namespace Kratos {

namespace Internals {

    template<class BaseType>
    class BaseRegisteredComponent {
    public:
        using BaseComponentType = BaseType;
        explicit BaseRegisteredComponent() {}
        virtual ~BaseRegisteredComponent() = default;
        virtual typename BaseType::Pointer Create(VoxelMeshGeneratorModeler& rModeler, Parameters Settings) const = 0;
    };

    template<class BaseType, class RegisteredType>
    class RegisteredComponent: public BaseRegisteredComponent<BaseType> {
    public:
        explicit RegisteredComponent() {}
        ~RegisteredComponent() override = default;
        typename BaseType::Pointer Create(VoxelMeshGeneratorModeler& rModeler, Parameters Settings) const override
        {
            return Kratos::make_shared<RegisteredType>(rModeler, Settings);
        }
    };

}

template<class ComponentType>
class VoxelMesherComponentFactory: public FactoryBase {
    std::string mName;
public:
    using RegisteredType = Internals::BaseRegisteredComponent<ComponentType>;

    explicit VoxelMesherComponentFactory(): mName("voxel mesher component") {};

    explicit VoxelMesherComponentFactory(const std::string& rName): mName(rName) {};

    ~VoxelMesherComponentFactory() override = default;

    bool Has(const std::string& rComponentName) const override
    {
        return KratosComponents<RegisteredType>::Has(rComponentName);
    }

    typename ComponentType::Pointer Create(VoxelMeshGeneratorModeler& rModeler, Parameters Settings) const
    {
        const auto type = Settings["type"].GetString();
        KRATOS_ERROR_IF_NOT(Has(type))
            << "Trying to construct " << mName << " with unregistered type \"" << type << "\"" << std::endl
            << "The list of available options are: " << std::endl
            << KratosComponents<RegisteredType>() << std::endl;
        return KratosComponents<RegisteredType>::Get(type).Create(rModeler, Settings);
    }

    template<class SpecializedType>
    static void Register(const std::string& rName) {
        static const Internals::RegisteredComponent<ComponentType, SpecializedType> component{};
        KratosComponents<RegisteredType>::Add(rName, component);
    }

};

}
