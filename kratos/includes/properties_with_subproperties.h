//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Fernando Rastellini
//

#if !defined(KRATOS_PROPERTIES_WITH_SUBPROPERTIES_H_INCLUDED )
#define  KRATOS_PROPERTIES_WITH_SUBPROPERTIES_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/properties.h"

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
 * @class PropertiesWithSubProperties
 * @ingroup KratosCore
 * @brief This is a property which contains at the same time other subproperties
 * @details This class derives from the base Properties and it has all the methods. This class overrides the methods related with Subproperties from the base class
 * @author Vicente Mataix Ferrandiz
 * @author Fernando Rastellini
 */
class PropertiesWithSubProperties
    : public Properties
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PropertiesWithSubProperties
    KRATOS_CLASS_POINTER_DEFINITION(PropertiesWithSubProperties);

    /// The Properties base class
    typedef Properties BaseType;

    /// The index type
    typedef BaseType::IndexType IndexType;

    /// The list of subproperties type
    typedef BaseType::SubPropertiesListType SubPropertiesListType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param NewId The Id of the new properties
     */
    explicit PropertiesWithSubProperties(IndexType NewId = 0)
        : BaseType(NewId),
          mSubPropetiesList(SubPropertiesListType(0))
        {}

    /**
     * @brief Constructor with pointer to properties
     * @param pProperties A pointer to properties
     */
    explicit PropertiesWithSubProperties(Properties::Pointer pProperties)
        : BaseType(*pProperties),
          mSubPropetiesList(SubPropertiesListType(0))
        {}

    /**
     * @brief Constructor with pointer and list of subproperties
     * @param pProperties A pointer to properties
     * @param rSubPropetiesList The map containing the subproperties
     */
    explicit PropertiesWithSubProperties(
        Properties::Pointer pProperties,
        SubPropertiesListType& rSubPropetiesList
        ) : BaseType(*pProperties),
            mSubPropetiesList(rSubPropetiesList)
        {}

    /**
     * @brief Constructor with reference to properties
     * @param rProperties A reference to properties
     */
    explicit PropertiesWithSubProperties(Properties& rProperties)
        : BaseType(rProperties),
          mSubPropetiesList(SubPropertiesListType(0))
        {}

    /**
     * @brief Constructor with reference to properties and list of subproperties
     * @param rProperties A reference to properties
     * @param rSubPropetiesList The map containing the subproperties
     */
    explicit PropertiesWithSubProperties(
        Properties& rProperties,
        SubPropertiesListType& rSubPropetiesList
        ) : BaseType(rProperties),
            mSubPropetiesList(rSubPropetiesList)
        {}

    /// Copy constructor.
    PropertiesWithSubProperties(const PropertiesWithSubProperties& rOther)
        : BaseType(rOther),
          mSubPropetiesList(rOther.mSubPropetiesList)
        {}

    /// Destructor.
    ~PropertiesWithSubProperties() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PropertiesWithSubProperties& operator=(const PropertiesWithSubProperties& rOther)
    {
        BaseType::operator=(rOther);
        mSubPropetiesList = rOther.mSubPropetiesList;
        return *this;
    }

    /**
     * @brief This method returns the number of subproperties
     * @return The current number of subproperties
     */
    std::size_t NumberOfSubproperties() override
    {
        return mSubPropetiesList.size();
    }

    /**
     * @brief This method insert a new property into the list of subproperties
     * @param pNewSubProperty The new property to be added
     */
    void AddSubProperty(Properties::Pointer pNewSubProperty) override
    {
        mSubPropetiesList.insert(std::pair<IndexType, Properties::Pointer>(pNewSubProperty->Id(), pNewSubProperty));
    }

    /**
     * @brief This method gets the subproperty from the index corresponding to the property id
     * @param SubPropertyIndex The index of the subproperty to be get
     * @return The pointer to the subproperty of interest
     */
    Properties::Pointer GetSubProperty(const IndexType SubPropertyIndex) override
    {
        return mSubPropetiesList[SubPropertyIndex];
    }

    /**
     * @brief This method gets the subproperty from the index corresponding to the property id
     * @param SubPropertyIndex The index of the subproperty to be get
     * @return The pointer to the subproperty of interest
     */
    Properties::Pointer const& GetSubProperty(const IndexType SubPropertyIndex) const override
    {
        return (mSubPropetiesList.find(SubPropertyIndex))->second;
    }

    /**
     * @brief This method returns the whole list of subproperties
     * @return The whole lis of subproperties
     */
    SubPropertiesListType& GetSubProperties() override
    {
        return SubPropetiesList();
    }

    /**
     * @brief This method returns the whole list of subproperties
     * @return The whole lis of subproperties
     */
    SubPropertiesListType const& GetSubProperties() const override
    {
        return SubPropetiesList();
    }

    /**
     * @brief This method set the whole list of subproperties
     * @param rSubPropetiesList The list of subproperties
     */
    void SetSubProperties(SubPropertiesListType& rSubPropetiesList) override
    {
        mSubPropetiesList = rSubPropetiesList;
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the whole list of subproperties
     * @return The whole lis of subproperties
     */
    SubPropertiesListType& SubPropetiesList() override
    {
        return mSubPropetiesList;
    }

    /**
     * @brief This method returns the whole list of subproperties (constant)
     * @return The whole lis of subproperties
     */
    SubPropertiesListType const& SubPropetiesList() const override
    {
        return mSubPropetiesList;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PropertiesWithSubProperties";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream <<  "PropertiesWithSubProperties";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        rOStream << "\nThis properties contains the following subproperties " << mSubPropetiesList.size() << " subproperties" << std::endl;
        for (auto& subprop : mSubPropetiesList) {
            (subprop.second)->PrintData(rOStream);
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

    SubPropertiesListType mSubPropetiesList; /// The vector containing the list of subproperties

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
        rSerializer.save("SubPropetiesList", mSubPropetiesList);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
        rSerializer.load("SubPropetiesList", mSubPropetiesList);
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

}; // Class PropertiesWithSubProperties

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PropertiesWithSubProperties& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PropertiesWithSubProperties& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PROPERTIES_WITH_SUBPROPERTIES_H_INCLUDED  defined


