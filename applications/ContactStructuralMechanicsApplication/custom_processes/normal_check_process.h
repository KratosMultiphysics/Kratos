// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_NORMAL_CHECK_PROCESS_H_INCLUDED )
#define  KRATOS_NORMAL_CHECK_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"

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
 * @class NormalCheckProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process checks the normal
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) NormalCheckProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    // Some geometrical definitions
    typedef Node<3>                                              NodeType;
    typedef Point                                               PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;

    /// Definition of geometries
    typedef Geometry<NodeType>                               GeometryType;
    typedef Geometry<PointType>                         GeometryPointType;

    /// The definition of zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of NormalCheckProcess
    KRATOS_CLASS_POINTER_DEFINITION( NormalCheckProcess );

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the normal check process
     * @param rModelPart The model part to be considered
     * @param ThisParameters The configuration parameters
     */
    NormalCheckProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters =  Parameters(R"({})")
        ) : mrModelPart(rModelPart),
            mParameters(ThisParameters)
    {
        const Parameters default_parameters = GetDefaultParameters();
        mParameters.ValidateAndAssignDefaults(default_parameters);
    }

    virtual ~NormalCheckProcess()= default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /************************************ GET INFO *************************************/
    /***********************************************************************************/

    std::string Info() const override
    {
        return "NormalCheckProcess";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
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

    ModelPart& mrModelPart;  /// The model part to be considered
    Parameters mParameters;  /// The configuration parameters

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
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

}; // Class NormalCheckProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

inline std::istream& operator >> (std::istream& rIStream,
                                  NormalCheckProcess& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

inline std::ostream& operator << (std::ostream& rOStream,
                                  const NormalCheckProcess& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_NORMAL_CHECK_PROCESS_H_INCLUDED  defined
