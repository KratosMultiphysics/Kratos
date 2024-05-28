// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef VARIANT_DATA_CONTAINER_INCLUDE_HPP
#define VARIANT_DATA_CONTAINER_INCLUDE_HPP

//// STL includes
#include <variant>
#include <iostream>
#include <string>
#include <algorithm>
#include <exception>
#include <type_traits>
#include <vector>
////Project includes
#include "includes/define.hpp"

namespace queso {

/**
 * @class  Component
 * @author Manuel Messmer
 * @brief  Container for all available Types of variant data container. Available Types are:
 *         PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType.
 * @details Stores Name (Key) and Value of Parameter.
**/
class Component {
public:
    ///@name Type Definitions
    ///@{
    typedef std::variant<PointType, Vector3i, bool, double, unsigned long, std::string, IntegrationMethodType> ComponentType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    Component(std::string Name, ComponentType NewComponent) : mName(Name), mComponents(NewComponent)
    {}

    ///@}
    ///@name Operations
    ///@{

    /// Returns Value (const)
    const ComponentType& Get() const{
        return mComponents;
    }

    /// Returns Value (non-const)
    ComponentType& Get(){
        return mComponents;
    }

    /// Returns Name (Key)
    const std::string& Name() const {
        return mName;
    }

private:
    ///@}
    ///@name Private Member variables
    ///@{

    std::string mName;
    ComponentType mComponents;

    ///@}

}; // End class Component


/**
 * @class  VariantDataContainer
 * @author Manuel Messmer
 * @brief  Pure base class for containers to store parameters of different types. Derived class must override
 *         GetDefaultComponents() and GetAvailableComponents().
 *         EnsureValidValues() may be overriden to set e.g. limits on certain parameters.
 * @see Component
**/
class VariantDataContainer {
public:

    typedef std::vector<Component> ComponentVectorType;
    typedef std::vector<std::pair<std::string, const std::type_info*>> AvailableComponentVectorType;

    ///@name std::visit Structs
    ///@{
    struct TypeVisit {
        const std::type_info* operator()(const PointType& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const Vector3i& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const unsigned long& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const double& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const std::string& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const bool& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const IntegrationMethodType& rValue){return &typeid(rValue); };
    };

    struct PrintVisit {
        PrintVisit(std::ostream& rOStream) : mOStream(rOStream){}
        void operator()(const PointType& rValue){mOStream << rValue; };
        void operator()(const Vector3i& rValue){mOStream << rValue;};
        void operator()(const unsigned long& rValue){ mOStream << rValue;};
        void operator()(const double& rValue){mOStream << rValue;};
        void operator()(const std::string& rValue){mOStream << rValue;};
        void operator()(const bool& rValue){mOStream << rValue; };
        void operator()(const IntegrationMethodType& rValue){ mOStream << rValue; };

    private:
        std::ostream& mOStream;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    VariantDataContainer() {
    }

    /// Constructor with list of components
    VariantDataContainer(ComponentVectorType& Component) : mComponents(Component) {
    }

    /// Destructor
    virtual ~VariantDataContainer() = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Set values to given component. If component is not stored yet, component is added.
    ///        Values that are supposed to be 'Double' or 'Vector3d' but are falsely given as 'int' or 'Vector3i'
    ///        are casted to their correct types. This allows the values in the 'QuESoParameters.json' to be provided as '0', '[0, 0, 0]',
    ///        which are supposed to be '0.0', '[0.0, 0.0, 0.0]'.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @param rValue New value of parameter.
    template<typename type>
    void Set(const std::string& rName, const type& rValue){
        auto p_value = pFind<type>(rName);
        if( p_value ){ // Component already exists -> Override.
            *p_value = rValue;
        }
        else {
            // Try to cast 0 -> 0.0 and [0, 0, 0] -> [0.0, 0.0, 0.0] when applicable.
            if( !CastAmbiguousTypesAndSet(rName, rValue) ) {
               // Component does not exist -> Add new one.
               mComponents.push_back(Component(rName, rValue));
            }
        }
        EnsureValidValues();
        CheckComponents();
    }

    /// @brief Returns true if parameter component exists.
    /// @tparam type
    /// @param rName Name (Key) of component.
    /// @return bool
    template<typename type>
    bool Has( const std::string& rName ) const {
        const auto p_value = pFind<type>(rName);
        if( p_value ){
            return true;
        }
        return false;
    }

    /// @brief Get parameter value. Throws an error if parameter is not available.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @return const type&
    template<typename type>
    const type& Get( const std::string& rName ) const {
        const auto p_value = pFind<type>(rName);
        if( p_value ){
            return *p_value;
        }

        QuESo_ERROR << "Component: '" + rName + "' not found.\n";
    }

    /// @brief Print all components (Name, Value).
    /// @param rOStream
    void PrintInfo(std::ostream& rOStream) const {
        for( auto& value : mComponents ){
            rOStream << value.Name() << ": ";
            std::visit(PrintVisit(rOStream), value.Get());
            rOStream << "\n";
        }
    }

    /// @brief Adds all default values (GetDefaultComponents()) to current set of parameters.
    ///        GetDefaultComponents() must be overriden in derived class.
    /// @see GetDefaultComponents()
    void AddDefaults(){
        const auto& r_default_components = GetDefaultComponents();
        for( auto& value : r_default_components){
            AddValueIfTypesMatch<PointType>(value);
            AddValueIfTypesMatch<Vector3i>(value);
            AddValueIfTypesMatch<bool>(value);
            AddValueIfTypesMatch<double>(value);
            AddValueIfTypesMatch<unsigned long>(value);
            AddValueIfTypesMatch<std::string>(value);
            AddValueIfTypesMatch<IntegrationMethodType>(value);
        }
    }

    /// @brief Checks if all components are part of GetAvailableComponents(). Also checks if all
    ///        components have the correct Types.
    ///        GetAvailableComponents() must be overriden in derived class.
    /// @see GetAvailableComponents()
    void CheckComponents() const {
        const auto& r_available_components = GetAvailableComponents();
        for( auto& r_components : mComponents ){
            const auto& current_name = r_components.Name();

            const auto p_pair_found = std::find_if( r_available_components.begin(), r_available_components.end(),
                [&current_name](const auto& rPair) { return rPair.first == current_name; } );

            if( p_pair_found != r_available_components.end() ){
                const auto p_current_type_info = std::visit(TypeVisit{}, r_components.Get());
                const auto p_ref_type_info = p_pair_found->second;
                if( !(*p_current_type_info ==  *p_ref_type_info) ){ // If type is wrong.
                    QuESo_ERROR << "Name: '" + current_name +
                        "' is not provided with correct Type.\n";
                }
            }
            else { // If name is part of r_available_components.
                QuESo_ERROR << "Name: '" + current_name +
                    "' is not a valid Parameter.\n";
            }
        }
    }

    /// @brief Interface to set limits of certain parameters.
    ///        May be overriden in derived class.
    virtual void EnsureValidValues() {};

    /// @brief Interface to define default parameters in derived class.
    ///        Must be overriden in derived class.
    /// @return ComponentVectorType
    virtual const ComponentVectorType& GetDefaultComponents() const = 0;

    /// @brief Interface to define available parameters in derived class.
    ///        Must be overriden in derived class.
    /// @return AvailableComponentVectorType
    virtual const AvailableComponentVectorType& GetAvailableComponents() const = 0;

private:

    ///@}
    ///@name Private operations
    ///@{

    /// @brief Adds Value to mComponents if Types of 'type' and 'rValues' match
    ///        and 'rValues' does not yet exist.
    /// @tparam type
    /// @param rValues New Component.
    template<typename type>
    void AddValueIfTypesMatch(const Component& rValues){
        auto p_type_id = std::visit(TypeVisit{}, rValues.Get());
        if( *p_type_id == typeid(type) ){
            auto p_value = pFind<type>(rValues.Name());
            if( !p_value ){
                mComponents.push_back( rValues );
            }
        }
    }

    /// @brief Casts values that are supposed to be 'Double' or 'Vector3d' but are falsely given as 'int' or 'Vector3i' to their correct types
    ///        and calls Set() again. Correct type is taken from GetAvailableComponents().
    /// @param rName Name (Key) of parameter.
    /// @param rValue New value of parameter.
    /// @param return bool. True if types were succesfully cast.
    bool CastAmbiguousTypesAndSet(const std::string& rName, const Component::ComponentType& rValue){
        const auto& r_available_components = GetAvailableComponents();
        const auto p_pair_found = std::find_if( r_available_components.begin(), r_available_components.end(),
                [&rName](const auto& rPair) { return rPair.first == rName; } );

        if( p_pair_found != r_available_components.end() ){
            const auto p_current_type_info = std::visit(TypeVisit{}, rValue);
            const auto p_ref_type_info = p_pair_found->second;
            // 0 -> 0.0 when required.
            if( *p_current_type_info == typeid(unsigned long) && *p_ref_type_info == typeid(double)) {
                auto p_value = const_cast<unsigned long*>(std::get_if<unsigned long>(&rValue));
                Set<double>(rName, static_cast<double>(*p_value));
                return true;
            }
            // [0, 0, 0] -> [0.0, 0.0, 0.0] when required.
            else if( *p_current_type_info == typeid(Vector3i) && *p_ref_type_info == typeid(Vector3d)) {
                auto p_value = const_cast<Vector3i*>(std::get_if<Vector3i>(&rValue));
                const auto& r_value = *p_value;
                Set<Vector3d>(rName, Vector3d{static_cast<double>(r_value[0]),
                                              static_cast<double>(r_value[1]),
                                              static_cast<double>(r_value[2])} );
                return true;
            }
        }

        return false;
    }

    /// @brief Returns Value of component if it exists. Returns nullptr otherwise. Const version.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @return const type*
    template<class type>
    const type* pFind( const std::string& rName ) const {
        for( auto& r_component : mComponents){
            const type* p_value = std::get_if<type>(&r_component.Get());
            if( p_value ){
                if( r_component.Name() == rName ){
                    return p_value;
                }
            }
        }
        return nullptr;
    }

    /// @brief Returns Value of component if it exists. Returns nullptr otherwise. Non-Const version.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @return const type*
    template<class type>
    type* pFind( const std::string& rName ) {
        for( auto& r_component : mComponents){
            type* p_value = std::get_if<type>(&r_component.Get());
            if( p_value ){
                if( r_component.Name() == rName ){
                    return p_value;
                }
            }
        }
        return nullptr;
    }

    ///@}
    ///@name Private member variables
    ///@{

    // List of components.
    ComponentVectorType mComponents;

}; // End of class VariantDataContainer
///@}
} // End Namespace queso
#endif // VARIANT_DATA_CONTAINER_INCLUDE_HPP