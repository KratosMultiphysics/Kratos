//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

# pragma once

// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/accessor.h"
#include "includes/node.h"
#include "includes/indexed_object.h"
#include "containers/data_value_container.h"
#include "includes/process_info.h"
#include "includes/table.h"
#include "utilities/string_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class Properties
 * @ingroup KratosCore
 * @brief Properties encapsulates data shared by different Elements or Conditions. It can store any type of data and provides a variable base access to them.
 * @details These are all parameters that can be shared between Element. Usually material parameters are common for a set of element, so this category of data is referred as properties.
 * But in general it can be any common parameter for a group of Elements. Sharing these data as properties reduces the memory used by the application and also helps updating them if necessary.
 * As mentioned before Properties is a shared data container between Elements or Conditions. In finite element problems there are several parameters which are the same for a set of elements and conditions.
 * Thermal conductivity, elasticity of the material and viscosity of the fluid are examples of these parameters. Properties holds these data and is shared by elements or Conditions. This eliminates memory overhead due to redundant copies of these data for each element and Condition. Properties also can be used to access nodal data if it is necessary.
 * It is important to mention that accessing the nodal data via Properties is not the same as accessing it via Node. When user asks Properties for a variable data in a Node, the process starts with finding the variable in the Properties data container and if it does not exist then get it from Node.
 * This means that the priority of data is with the one stored in Properties and then in Node.
 * @author Pooyan Dadvand
 */
class Properties : public IndexedObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Properties
    KRATOS_CLASS_POINTER_DEFINITION(Properties);

#ifdef  _WIN32 // work around for windows int64_t error
    using int64_t = __int64;
#endif
    using BaseType = IndexedObject;

    using ContainerType = DataValueContainer;

    using GeometryType = Geometry<Node> ;

    using IndexType = std::size_t;

    using TableType = Table<double>;

    using KeyType = IndexType;

    using AccessorPointerType = Accessor::UniquePointer;

    using AccessorsContainerType = std::unordered_map<KeyType, AccessorPointerType>;

    using TablesContainerType = std::unordered_map<std::size_t, TableType>; // This is a provisional implementation and should be changed to hash. Pooyan.

    /// Properties container. A vector set of properties with their Id's as key.
    using SubPropertiesContainerType = PointerVectorSet<Properties, IndexedObject>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit Properties(IndexType NewId = 0)
    : BaseType(NewId)
    , mData()
    , mTables()
    , mSubPropertiesList()
    , mAccessors() {}

    /// Default of properties with subproperties
    explicit Properties(IndexType NewId, const SubPropertiesContainerType& SubPropertiesList)
    : BaseType(NewId)
    , mData()
    , mTables()
    , mSubPropertiesList(SubPropertiesList)
    , mAccessors() {}

    /// Copy constructor.
    Properties(const Properties& rOther)
    : BaseType(rOther)
    , mData(rOther.mData)
    , mTables(rOther.mTables)
    , mSubPropertiesList(rOther.mSubPropertiesList)
    {
        for (auto& r_item : rOther.mAccessors) {
            const auto key = r_item.first;
            const auto& rp_accessor = r_item.second;
            mAccessors.emplace(key, rp_accessor->Clone());
        }
    }

    /// Destructor.
    ~Properties() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Properties& operator=(const Properties& rOther)
    {
        BaseType::operator=(rOther);
        mData = rOther.mData;
        mTables = rOther.mTables;
        mSubPropertiesList = rOther.mSubPropertiesList;
        for (auto& r_item : rOther.mAccessors) {
            const auto key = r_item.first;
            const auto& rp_accessor = r_item.second;
            mAccessors.emplace(key, rp_accessor->Clone());
        }
        return *this;
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rV)
    {
        return GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rV) const
    {
        return GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type& operator[](const TVariableType& rV)
    {
        return GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator[](const TVariableType& rV) const
    {
        return GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rV, Node& rThisNode)
    {
        return GetValue(rV, rThisNode);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rV, Node const& rThisNode) const
    {
        return GetValue(rV, rThisNode);
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rV, Node& rThisNode, IndexType SolutionStepIndex)
    {
        return GetValue(rV, rThisNode, SolutionStepIndex);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rV, Node const& rThisNode, IndexType SolutionStepIndex) const
    {
        return GetValue(rV, rThisNode, SolutionStepIndex);
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rV, Node& rThisNode, ProcessInfo const& rCurrentProcessInfo)
    {
        return GetValue(rV, rThisNode, rCurrentProcessInfo.GetSolutionStepIndex());
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rV, Node const& rThisNode, ProcessInfo const& rCurrentProcessInfo) const
    {
        return GetValue(rV, rThisNode, rCurrentProcessInfo.GetSolutionStepIndex());
    }

    ///@}
    ///@name Operations
    ///@{

    template<class TVariableType>
    void Erase(const TVariableType& rV)
    {
        mData.Erase(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rVariable)
    {
        return mData.GetValue(rVariable);
    }

    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rVariable) const
    {

        return mData.GetValue(rVariable);
    }

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rVariable, Node& rThisNode)
    {
        if (mData.Has(rVariable))
            return mData.GetValue(rVariable);
        return rThisNode.GetValue(rVariable);
    }

    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rVariable, Node const& rThisNode) const
    {
        if (mData.Has(rVariable))
            return mData.GetValue(rVariable);
        return rThisNode.GetValue(rVariable);
    }

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rVariable, Node& rThisNode, IndexType SolutionStepIndex)
    {
        if (mData.Has(rVariable))
            return mData.GetValue(rVariable);
        return rThisNode.GetValue(rVariable, SolutionStepIndex);
    }

    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rVariable, Node const& rThisNode, IndexType SolutionStepIndex) const
    {
        if (mData.Has(rVariable))
            return mData.GetValue(rVariable);
        return rThisNode.GetValue(rVariable, SolutionStepIndex);
    }

    /*
    Custom GetValue in which we check the Accessor
    */
    template<class TVariableType>
    typename TVariableType::Type GetValue(const TVariableType& rVariable, const GeometryType& rGeometry, const Vector& rShapeFunctionVector, const ProcessInfo& rProcessInfo) const
    {
        auto it_value = mAccessors.find(rVariable.Key());
        if (it_value != mAccessors.end()) {
            return (it_value->second)->GetValue(rVariable, *this, rGeometry, rShapeFunctionVector, rProcessInfo);
        } else {
            return mData.GetValue(rVariable);
        }
    }

    template<class TVariableType>
    void SetValue(TVariableType const& rV, typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rV, rValue);
    }

    bool HasVariables() const
    {
        return !mData.IsEmpty();
    }

    /**
     * @brief Set the Accessor object
     * This method sets a variable-accessor pair in current properties accessor container
     * @tparam TVariableType The variable type
     * @param rVariable Variable to which the accessor will refer to
     * @param pAccessor Pointer to the accessor instance
     */
    template <class TVariableType>
    void SetAccessor(const TVariableType& rVariable, AccessorPointerType pAccessor)
    {
        mAccessors.emplace(rVariable.Key(), std::move(pAccessor));
    }

    /**
     * @brief Get the Accessor object
     * If exists, this method returns a pointer to the requested variable accessor
     * If doesn't exist, the method throws an error
     * @tparam TVariableType The variable type
     * @param rVariable Variable to which the accessor refers to
     * @return AccessorPointerType& Pointer to the requested variable accessor
     */
    template <class TVariableType>
    Accessor& GetAccessor(const TVariableType& rVariable)
    {
        auto it_value = mAccessors.find(rVariable.Key());
        KRATOS_ERROR_IF(it_value == mAccessors.end())
            << "Trying to retrieve inexisting accessor for '" << rVariable.Name() << "' in properties " << Id() << "." << std::endl;
        return *(it_value->second);
    }

    /**
     * @brief Get the Accessor object
     * If exists, this method returns a pointer to the requested variable accessor
     * If doesn't exist, the method throws an error
     * @tparam TVariableType The variable type
     * @param rVariable Variable to which the accessor refers to
     * @return AccessorPointerType& Pointer to the requested variable accessor
     */
    template <class TVariableType>
    Accessor& GetAccessor(const TVariableType& rVariable) const
    {
        const auto it_value = mAccessors.find(rVariable.Key());
        KRATOS_ERROR_IF(it_value == mAccessors.end())
            << "Trying to retrieve inexisting accessor for '" << rVariable.Name() << "' in properties " << Id() << "." << std::endl;
        return *(it_value->second);
    }

    /**
     * @brief Get the Accessor object
     * If exists, this method returns a pointer to the requested variable accessor
     * If doesn't exist, the method throws an error
     * @tparam TVariableType The variable type
     * @param rVariable Variable to which the accessor refers to
     * @return AccessorPointerType& Pointer to the requested variable accessor
     */
    template <class TVariableType>
    AccessorPointerType& pGetAccessor(const TVariableType& rVariable)
    {
        const auto it_value = mAccessors.find(rVariable.Key());
        KRATOS_ERROR_IF(it_value == mAccessors.end())
            << "Trying to retrieve inexisting accessor for '" << rVariable.Name() << "' in properties " << Id() << "." << std::endl;
        return it_value->second;
    }

    /**
     * @brief Check if current properties have an accessor
     * This method checks if current properties have an accessor for the requested variable
     * @tparam TVariableType The variable type
     * @param rVariable Variable to which we are checking if an accessor exists
     * @return true If there is accessor for the requested variable
     * @return false If there is no accessor for the requested variable
     */
    template <class TVariableType>
    bool HasAccessor(const TVariableType& rVariable) const
    {
        return (mAccessors.find(rVariable.Key()) == mAccessors.end()) ? false : true;
    }

    template<class TXVariableType, class TYVariableType>
    TableType& GetTable(const TXVariableType& XVariable, const TYVariableType& YVariable)
    {
        return mTables[Key(XVariable.Key(), YVariable.Key())];
    }

    template<class TXVariableType, class TYVariableType>
    TableType const& GetTable(const TXVariableType& XVariable, const TYVariableType& YVariable) const
    {
        return mTables.at(Key(XVariable.Key(), YVariable.Key()));
    }

    template<class TXVariableType, class TYVariableType>
    void SetTable(const TXVariableType& XVariable, const TYVariableType& YVariable, TableType const& rThisTable)
    {
        mTables[Key(XVariable.Key(), YVariable.Key())] = rThisTable;
    }

    bool HasTables() const
    {
        return !mTables.empty();
    }

    bool IsEmpty() const
    {
        return !( HasVariables() || HasTables() );
    }

    int64_t Key(std::size_t XKey, std::size_t YKey) const
    {
        int64_t result_key = XKey;
        result_key = result_key << 32;
        result_key |= YKey; // I know that the key is less than 2^32 so I don't need zeroing the upper part
        return result_key;
    }

    /**
     * @brief This method returns the number of subproperties
     * @return The current number of subproperties
     */
    std::size_t NumberOfSubproperties() const
    {
        return mSubPropertiesList.size();
    }

    /**
     * @brief This method insert a new property into the list of subproperties
     * @param pNewSubProperty The new property to be added
     */
    void AddSubProperties(Properties::Pointer pNewSubProperty)
    {
        KRATOS_DEBUG_ERROR_IF(this->HasSubProperties(pNewSubProperty->Id())) << "SubProperties with ID: " << pNewSubProperty->Id() << " already defined" << std::endl;
        mSubPropertiesList.insert(mSubPropertiesList.begin(), pNewSubProperty);
    }

    /**
     * @brief This method checks if the subproperty exists from the index corresponding to the property id
     * @param SubPropertyIndex The index of the subproperty to be get
     * @return True if there is such subproperty, false otherwise
     */
    bool HasSubProperties(const IndexType SubPropertyIndex) const
    {
        return mSubPropertiesList.find(SubPropertyIndex) != mSubPropertiesList.end();
    }

    /**
     * @brief This method gets the subproperty from the index corresponding to the property id
     * @param SubPropertyIndex The index of the subproperty to be get
     * @return The pointer to the subproperty of interest
     */
    Properties::Pointer pGetSubProperties(const IndexType SubPropertyIndex)
    {
        // Looking into the database
        auto property_iterator = mSubPropertiesList.find(SubPropertyIndex);
        if (property_iterator != mSubPropertiesList.end()) {
            return *(property_iterator.base());
        } else {
            KRATOS_ERROR << "Subproperty ID: " << SubPropertyIndex << " is not defined on the current Properties ID: " << this->Id() << " creating a new one with ID: " << SubPropertyIndex << std::endl;
            return nullptr;
        }
    }

    /**
     * @brief This method gets the subproperty from the index corresponding to the property id (constant version)
     * @param SubPropertyIndex The index of the subproperty to be get
     * @return The pointer to the subproperty of interest
     */
    const Properties::Pointer pGetSubProperties(const IndexType SubPropertyIndex) const
    {
        // Looking into the database
        auto property_iterator = mSubPropertiesList.find(SubPropertyIndex);
        if (property_iterator != mSubPropertiesList.end()) {
            return *(property_iterator.base());
        } else {
            KRATOS_ERROR << "Subproperty ID: " << SubPropertyIndex << " is not defined on the current Properties ID: " << this->Id() << " creating a new one with ID: " << SubPropertyIndex << std::endl;
            return nullptr;
        }
    }

    /**
     * @brief This method gets the subproperty from the index corresponding to the property id
     * @param SubPropertyIndex The index of the subproperty to be get
     * @return The reference to the subproperty of interest
     */
    Properties& GetSubProperties(const IndexType SubPropertyIndex)
    {
        // Looking into the database
        auto property_iterator = mSubPropertiesList.find(SubPropertyIndex);
        if (property_iterator != mSubPropertiesList.end()) {
            return *(property_iterator);
        } else {
            KRATOS_ERROR << "Subproperty ID: " << SubPropertyIndex << " is not defined on the current Properties ID: " << this->Id() << " creating a new one with ID: " << SubPropertyIndex << std::endl;
            return *this;
        }
    }

    /**
     * @brief This method gets the subproperty from the index corresponding to the property id (constant version)
     * @param SubPropertyIndex The index of the subproperty to be get
     * @return The reference to the subproperty of interest
     */
    const Properties& GetSubProperties(const IndexType SubPropertyIndex) const
    {
        // Looking into the database
        if (mSubPropertiesList.find(SubPropertyIndex) != mSubPropertiesList.end()) {
            return *(mSubPropertiesList.find(SubPropertyIndex));
        } else {
            KRATOS_ERROR << "Subproperty ID: " << SubPropertyIndex << " is not defined on the current Properties ID: " << this->Id() << std::endl;
        }
    }

    /**
     * @brief This method returns the whole list of subproperties
     * @return The whole lis of subproperties
     */
    SubPropertiesContainerType& GetSubProperties()
    {
        return mSubPropertiesList;
    }

    /**
     * @brief This method returns the whole list of subproperties
     * @return The whole lis of subproperties
     */
    SubPropertiesContainerType const& GetSubProperties() const
    {
        return mSubPropertiesList;
    }

    /**
     * @brief This method set the whole list of subproperties
     * @param rSubPropertiesList The list of subproperties
     */
    void SetSubProperties(SubPropertiesContainerType& rSubPropertiesList)
    {
        mSubPropertiesList = rSubPropertiesList;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the whole data container
     * @return Data container
     */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use 'GetData()' instead.")
    ContainerType& Data()
    {
        return mData;
    }

    /**
     * @brief This method returns the whole data container (constant)
     * @return Data container
     */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use 'GetData()' instead.")
    ContainerType const& Data() const
    {
        return mData;
    }

    const ContainerType& GetData() const
    {
        return mData;
    }

    ContainerType& GetData()
    {
        return mData;
    }

    /**
     * @brief This method returns the tables
     * @return The whole lis of tables
     */
    TablesContainerType& Tables()
    {
        return mTables;
    }

    /**
     * @brief This method returns the tables (constant)
     * @return The whole lis of tables
     */
    TablesContainerType const& Tables() const
    {
        return mTables;
    }

    ///@}
    ///@name Inquiry
    ///@{

    template<class TVariableType>
    bool Has(TVariableType const& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    template<class TXVariableType, class TYVariableType>
    bool HasTable(const TXVariableType& XVariable, const TYVariableType& YVariable) const
    {
        return (mTables.find(Key(XVariable.Key(), YVariable.Key())) != mTables.end());
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Properties";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream <<  "Properties";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        // Id
        rOStream << "Id : " << this->Id() << "\n";

        // Data
        mData.PrintData(rOStream);

        // Tables
        if (mTables.size() > 0) {
            // Print the tables
            rOStream << "This properties contains " << mTables.size() << " tables\n";
            for (auto& r_table : mTables) {
                rOStream << "Table key: " << r_table.first << "\n";
                StringUtilities::PrintDataWithIdentation(rOStream, r_table.second);
            }
        }

        // Subproperties
        if (mSubPropertiesList.size() > 0) {
            // Print the subproperties
            rOStream << "\nThis properties contains " << mSubPropertiesList.size() << " subproperties\n";
            for (auto& r_subprop : mSubPropertiesList) {
                StringUtilities::PrintDataWithIdentation(rOStream, r_subprop);
            }
        }

        // Accessors
        if (mAccessors.size() > 0) {
            // Print the accessors
            rOStream << "\nThis properties contains " << mAccessors.size() << " accessors\n";
            for (auto& r_entry : mAccessors) {
                rOStream << "Accessor for variable key: " << r_entry.first << "\n";
                StringUtilities::PrintDataWithIdentation(rOStream, *r_entry.second);
            }
        }
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ContainerType mData;                        /// The data stored on the properties

    TablesContainerType mTables;                /// The tables contained on the properties

    SubPropertiesContainerType mSubPropertiesList; /// The vector containing the list of subproperties

    AccessorsContainerType mAccessors = {}; /// The map containing the variable and corresponding accessor pairs

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
        rSerializer.save("Data", mData);
        rSerializer.save("Tables", mTables);
        rSerializer.save("SubPropertiesList", mSubPropertiesList);
        std::vector<std::pair<const KeyType, Accessor*>> aux_accessors_container;
        for (auto& r_item : mAccessors) {
            const auto key = r_item.first;
            const auto& rp_accessor = r_item.second;
            aux_accessors_container.push_back(std::make_pair(key, &(*rp_accessor)));
        }
        rSerializer.save("Accessors", aux_accessors_container);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
        rSerializer.load("Data", mData);
        rSerializer.load("Tables", mTables);
        rSerializer.load("SubPropertiesList", mSubPropertiesList);
        std::vector<std::pair<const KeyType, Accessor*>> aux_accessors_container;
        rSerializer.load("Accessors", aux_accessors_container);
        for (auto& r_item : aux_accessors_container) {
            const auto key = r_item.first;
            const auto& rp_accessor = r_item.second;
            mAccessors.emplace(key, rp_accessor->Clone());
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class Properties

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Properties& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Properties& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.