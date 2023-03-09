//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Riccardo Rossi
//
//

# pragma once

// System includes
#include "includes/properties.h"
#include "geometries/geometry.h"

// External includes

// Project includes

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
 * PropertyAccessor sets a proper way of returning a value from the properties
 * It is used to properly interpolate material properties according to temperature, tables, etc...
*/
class PropertyAccessor
{
public:
    ///@name Type Definitions
    ///@{
    typedef Geometry<Node<3>> GeometryType;

    typedef std::function<double(const std::string&, const Properties&, const GeometryType&)> AccessorType;

    /// Pointer definition of ProcessInfo
    KRATOS_CLASS_POINTER_DEFINITION(PropertyAccessor);

    ///@}
    ///@name Life Cycle
    ///@{


    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief default constructor
     */
    PropertyAccessor(const double Value)
    {
        mValue = Value;
    }

    /// Assignment operator.
    PropertyAccessor& operator=(const PropertyAccessor& rOther)
    {
        mValue = rOther.mValue;
        mListOfAccessors = rOther.mListOfAccessors;
        return *this;
    }

    /**
     * @brief This method add a new acessor to a certain variable name
     */
    void AddAccessor(const std::string& rVariableName, AccessorType* pAccessor)
    {
        mListOfAccessors[rVariableName] = pAccessor;
    }

    /**
     * @brief This method return the value of the required variable
     */
    double GetProperty(const std::string& rVariableName, const Properties& rProperties, const GeometryType& rGeometry)
    {
        // If it is in the list, give back the corresponding accessor, otherwise give back value
        auto it = mListOfAccessors.find(rVariableName);
        if (it != mListOfAccessors.end()) {
            const auto& r_function = *(*it).second;
            return r_function(rVariableName, rProperties, rGeometry);
        } else {
            return mValue;
        }
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Solution Step Data
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

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

    double mValue;
    std::unordered_map<std::string, AccessorType*> mListOfAccessors;

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

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("Value", mValue);
        rSerializer.save("ListOfAccessors", mListOfAccessors);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("Value", mValue);
        rSerializer.load("ListOfAccessors", mListOfAccessors);
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

}; // Class ProcessInfo

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.


