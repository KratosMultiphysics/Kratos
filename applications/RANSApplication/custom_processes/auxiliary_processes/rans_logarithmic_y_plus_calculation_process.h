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

#if !defined(KRATOS_RANS_LOGARITHMIC_Y_PLUS_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_LOGARITHMIC_Y_PLUS_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

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
 * @brief Calculates y_plus value base on the logarithmic law
 *
 * This process calculates $y^+$ value based on the following formula:
 *
 * \[
 * 	u^+ = \frac{||\vvel||}{\vel_\tau} =
 * \begin{cases}
 *	\frac{1}{\kappa}ln\left(y^+\right) + \beta &\text{ for } y^+ > y^+_{limit}
 *\\ y^+ &\text{ for } y^+ \leq y^+_{limit} \end{cases}
 * \]
 * Where,
 * \[
 *	y^+ = \frac{\vel_\tau y}{\nu}
 * \]
 * \[
 *	y^+_{limit} = \frac{1}{\kappa}ln\left(y^+_{limit}\right) + \beta = \vel^+
 * \]
 *
 */

class KRATOS_API(RANS_APPLICATION) RansLogarithmicYPlusCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansLogarithmicYPlusCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansLogarithmicYPlusCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansLogarithmicYPlusCalculationProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~RansLogarithmicYPlusCalculationProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void Execute() override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    Model& mrModel;
    Parameters mrParameters;
    std::string mModelPartName;

    unsigned int mEchoLevel;
    int mStep;

    unsigned int mMaxIterations;
    double mTolerance;

    double mVonKarman;
    double mBeta;
    double mLimitYPlus;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateLogarithmicWallLawYplus(NodeType& rNode,
                                          const Variable<array_1d<double, 3>>& rVelocityVariable);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RansLogarithmicYPlusCalculationProcess& operator=(RansLogarithmicYPlusCalculationProcess const& rOther);

    /// Copy constructor.
    RansLogarithmicYPlusCalculationProcess(RansLogarithmicYPlusCalculationProcess const& rOther);

    ///@}

}; // namespace Kratos

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansLogarithmicYPlusCalculationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_LOGARITHMIC_Y_PLUS_CALCULATION_PROCESS_H_INCLUDED  defined
