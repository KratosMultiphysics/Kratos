//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_VARIABLE_DIFFERENCE_NORM_CALCULATION_UTILITY_H_INCLUDED)
#define KRATOS_RANS_VARIABLE_DIFFERENCE_NORM_CALCULATION_UTILITY_H_INCLUDED

// System includes
#include <tuple>

// External includes

// Project includes
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief This class is used to calculate difference norms of a given variable
 *
 * @tparam TDataType        Data type of the variable
 */
template <typename TDataType>
class RansVariableDifferenceNormsCalculationUtility
{
public:
    ///@name Pointer Definitions
    /// Pointer definition of RansVariableDifferenceNormsCalculationUtility
    KRATOS_CLASS_POINTER_DEFINITION(RansVariableDifferenceNormsCalculationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    RansVariableDifferenceNormsCalculationUtility(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable)
    : mrModelPart(rModelPart), mrVariable(rVariable)
    {
    }

    /**
     * Destructor
     */
    ~RansVariableDifferenceNormsCalculationUtility()
    {
        mData.clear();
    }

    ///@}
    ///@name Operations
    ///@{

    void InitializeCalculation();

    std::tuple<double, double> CalculateDifferenceNorm();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "RansVariableDifferenceNormsCalculationUtility";
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RansVariableDifferenceNormsCalculationUtility";
    }

private:
    ///@name Member Variables
    ///@{
    const ModelPart& mrModelPart;
    const Variable<TDataType>& mrVariable;
    std::vector<TDataType> mData;

    ///@}

}; // Class RansVariableDifferenceNormsCalculationUtility

///@}
///@name Input and output
///@{

template <typename TDataType>
inline std::istream& operator>>(
    std::istream& rIStream,
    RansVariableDifferenceNormsCalculationUtility<TDataType>& rThis);

/// output stream function
template <typename TDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansVariableDifferenceNormsCalculationUtility<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_VARIABLE_DIFFERENCE_NORM_CALCULATION_UTILITY_H_INCLUDED defined
