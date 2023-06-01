//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Riccardo Rossi
//                   Carlos Roig
//                   Ruben Zorrilla
//

# pragma once

// System includes

// External includes

// Project includes
#include "includes/accessor.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class TableAccessor
 * @ingroup Kratos Core
 * @brief This class defines the way a certain property is accessed according to a table.
 * @brief The tables are supposed to relate double <-> double type entities. 
 * @brief The input variable is suposed to be a nodally accesible one (either historical or not)
 * @author Alejandro Cornejo, Riccardo Rossi, Carlos Roig and Ruben Zorrilla
 */
class KRATOS_API(KRATOS_CORE) TableAccessor : public Accessor
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry type definition
    using GeometryType = Geometry<Node>;

    /// Variable type
    using VariableType = Variable<double>;

    /// BaseType
    using BaseType = Accessor;

    using SizeType = std::size_t;

    /// Pointer definition of TableAccessor
    KRATOS_CLASS_POINTER_DEFINITION(TableAccessor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    TableAccessor() = default;

    /// Custom constructor
    TableAccessor(VariableType& rInputVariable) : mInputVariable(rInputVariable)
    {}

    /// Destructor.
    ~TableAccessor() = default;

    /// Copy constructor
    TableAccessor(const TableAccessor& rOther) 
    : BaseType(rOther),
        mInputVariable(rOther.mInputVariable) 
    {}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Custom method to retrieve double type properties
     * @param rVariable The variable considered (double type properties)
     * @param rProperties The properties considered
     * @param rGeometry The geometry considered
     * @param rShapeFunctionVector The shape function of the GP considered
     * @param rProcessInfo The process info considered
     * @return The double type properties
     */
    double GetValue(
        const Variable<double>& rVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo
        ) const override;

    /**
     * @brief This computes a material property according to a certain
     * nodal Variable<double> table
     */
    double GetValueFromTable(
        const Variable<double> &rIndependentVariable,
        const Variable<double> &rDependentVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const;


    /**
     * @brief Returns the member input variable
     */
    VariableType& GetInputVariable()
    {
        return mInputVariable;
    }

    /**
     * @brief Sets the InputVariable to a VariableType
     */
    // void SetInputVariableType(const VariableType& rVariable)
    // {
    //     mInputVariable = rVariable;
    // }

    // Getting a pointer to the class
    Accessor::UniquePointer Clone() const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "TableAccessor" ;

        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const  {rOStream << "TableAccessor";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {rOStream << "TableAccessor class";}

    ///@}

private:

    ///@name Member Variables
    ///@{

    VariableType mInputVariable;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
        rSerializer.save("InputVariable", mInputVariable);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
        rSerializer.load("InputVariable", mInputVariable);
    }


}; // class
///@}

///@} addtogroup block

} // namespace Kratos