//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_FACTORY_H_INCLUDED )
#define  KRATOS_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"

namespace Kratos
{
///@name Type Definitions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class FactoryBase
 * @ingroup KratosCore
 * @brief Here we define some common methods
 * @details Defines the base class factory methods
 * @author Vicente Mataix Ferrandiz
 */
class FactoryBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FactoryBase
    KRATOS_CLASS_POINTER_DEFINITION(FactoryBase);

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    explicit FactoryBase(){}

    /** Destructor.
     */
    virtual ~FactoryBase(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks if the linear class is registered
     * @param rClassName The nanme of the class
     * @return True if registered, false otherwise
     */
    virtual bool Has(const std::string& rClassName) const
    {
        KRATOS_ERROR << "Methods must be implemented in the derived class" << std::endl;
        return false;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "FactoryBase";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << this->Info() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/**
 * @class Factory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of class
 * @details Defines the base class factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TClass The class to create the factory
 */
template<typename TClass>
class Factory
    : public FactoryBase
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the class
    typedef TClass ClassType;

    /// Pointer definition of Factory
    KRATOS_CLASS_POINTER_DEFINITION(Factory);

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    explicit Factory(){}

    /** Destructor.
     */
    virtual ~Factory() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks if the linear class is registered
     * @param rClassName The nanme of the class
     * @return True if registered, false otherwise
     */
    bool Has(const std::string& rClassName) const override
    {
        return KratosComponents<ClassType>::Has( rClassName );
    }

    /**
     * @brief This method creates a new class
     * @param Arguments The arguments of the method
     * @return The pointer to the class of interest
     * @tparam TArgumentsType Variadic template arguments
     */
    template<typename... TArgumentsType >
    typename ClassType::Pointer Create(TArgumentsType&&... Arguments) const
    {
        std::tuple<TArgumentsType...> args(Arguments...);
        constexpr std::size_t args_size = std::tuple_size<std::tuple<TArgumentsType...>>::value;
        Kratos::Parameters Settings = std::get<args_size - 1>(args);
        const std::string& r_name = Settings["name"].GetString();
        KRATOS_ERROR_IF_NOT(Has( r_name )) << "Trying to construct a class with type name= " << r_name << "\n" <<
                                            "Which does not exist. The list of available options (for currently loaded applications) are: \n" <<
                                            KratosComponents<ClassType>() << std::endl;
        const ClassType& aux = KratosComponents<ClassType>::Get( r_name );
        return aux.Create(std::forward<TArgumentsType>(Arguments)...);
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Factory";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info() << std::endl;
        const auto factory_components = KratosComponents<ClassType>::GetComponents();
        for (const auto& r_comp : factory_components) {
            rOStream << "\t" << r_comp.first << std::endl;
        }
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
};

namespace Internals
{
    template <typename TBaseCategoryType>
    class RegisteredPrototypeBase {
        public:
        RegisteredPrototypeBase() = default;
    };

    template <typename TClassType, typename TBaseCategoryType>
    class RegisteredPrototype : public  RegisteredPrototypeBase<TBaseCategoryType> {
        public:
        explicit RegisteredPrototype(const std::string& rName, const TClassType& rPrototype)
        {
            KratosComponents<TBaseCategoryType>::Add(rName, rPrototype);
        }
    };

}

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template<class TClass>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Factory<TClass>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_FACTORY_H_INCLUDED  defined
