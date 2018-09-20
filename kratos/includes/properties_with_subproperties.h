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
 */
class PropertiesWithSubProperties
    : public Properties
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PropertiesWithSubProperties
    KRATOS_CLASS_POINTER_DEFINITION(PropertiesWithSubProperties);

    typedef Properties BaseType;

    typedef BaseType::IndexType IndexType;

    typedef BaseType::SubPropertiesListType SubPropertiesListType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PropertiesWithSubProperties(IndexType NewId = 0)
        : BaseType(NewId), mSubPropetiesList(SubPropertiesListType(0)) {}

    PropertiesWithSubProperties(Properties::Pointer pProperties)
        : BaseType(*pProperties), mSubPropetiesList(SubPropertiesListType(0)) {}

    PropertiesWithSubProperties(Properties::Pointer pProperties, SubPropertiesListType& rSubPropetiesList)
        : BaseType(*pProperties), mSubPropetiesList(rSubPropetiesList) {}

    PropertiesWithSubProperties(Properties& rProperties)
        : BaseType(rProperties), mSubPropetiesList(SubPropertiesListType(0)) {}

    PropertiesWithSubProperties(Properties& rProperties, SubPropertiesListType& rSubPropetiesList)
        : BaseType(rProperties), mSubPropetiesList(rSubPropetiesList) {}

    /// Copy constructor.
    PropertiesWithSubProperties(const PropertiesWithSubProperties& rOther) : BaseType(rOther), mSubPropetiesList(rOther.mSubPropetiesList) {}

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

    std::size_t NumberOfSubproperties() override
    {
        return mSubPropetiesList.size();
    }

    void AddSubProperty(Properties::Pointer pNewSubProperty) override
    {
        mSubPropetiesList.insert(std::pair<IndexType, Properties::Pointer>(pNewSubProperty->Id(), pNewSubProperty));
    }

    Properties::Pointer GetSubProperty(const IndexType SubPropertyIndex) override
    {
        return mSubPropetiesList[SubPropertyIndex];
    }

    SubPropertiesListType& GetSubProperties() override
    {
        return SubPropetiesList();
    }

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

    SubPropertiesListType& SubPropetiesList() override
    {
        return mSubPropetiesList;
    }

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
        rOStream << "This properties contains the following subproperties" << mSubPropetiesList.size() << " subproperties";
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


