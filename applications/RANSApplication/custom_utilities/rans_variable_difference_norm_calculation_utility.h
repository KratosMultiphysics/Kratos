//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_VARIABLE_DIFFERENCE_NORM_CALCULATION_UTILITYH_INCLUDED)
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

template <typename TDataType>
class RansVariableDifferenceNormsCalculationUtility
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of RansVariableDifferenceNormsCalculationUtility
    KRATOS_CLASS_POINTER_DEFINITION(RansVariableDifferenceNormsCalculationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    RansVariableDifferenceNormsCalculationUtility(const ModelPart& rModelPart,
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
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void InitializeCalculation();

    std::tuple<double, double> CalculateDifferenceNorm();

    ///@}

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
    const ModelPart& mrModelPart;
    const Variable<TDataType>& mrVariable;
    std::vector<double> mData;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

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

}; // Class RansVariableDifferenceNormsCalculationUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

template <typename TDataType>
inline std::istream& operator>>(std::istream& rIStream,
                                RansVariableDifferenceNormsCalculationUtility<TDataType>& rThis);

/// output stream function
template <typename TDataType>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansVariableDifferenceNormsCalculationUtility<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_VARIABLE_DIFFERENCE_NORM_CALCULATION_UTILITYDED defined
