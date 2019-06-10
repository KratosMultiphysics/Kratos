// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(TRESCA_FRICTIONAL_LAW_H_DEFINED )
#define  TRESCA_FRICTIONAL_LAW_H_DEFINED

// System includes

// External includes

// Project includes
#include "custom_frictional_laws/frictional_law.h"

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
 * @class TrescaFrictionalLaw
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class defines the Tresca frictional laws
 * @details This class defines the most basic friction model
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) TrescaFrictionalLaw
    : public FrictionalLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// Base class definition
    typedef FrictionalLaw BaseType;

    /// Node definition
    typedef Node<3> NodeType;

    /// Index type definition
    typedef std::size_t IndexType;

    /// Size type definition
    typedef std::size_t SizeType;

    /// Zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Counted pointer of TrescaFrictionalLaw
    KRATOS_CLASS_POINTER_DEFINITION( TrescaFrictionalLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    TrescaFrictionalLaw()
    {
    }

    ///Copy constructor  (not really required)
    TrescaFrictionalLaw(const TrescaFrictionalLaw& rhs)
    {
    }

    /// Destructor.
    ~TrescaFrictionalLaw()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    /**
     * @brief This method computes the threshold value considered for computing friction
     * @param rNode The node where the threshold value is obtained
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     */
    double GetThresholdValue(
        const NodeType& rNode,
        const Condition& rCondition,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    // TODO: Add derivatives

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
    std::string Info() const override
    {
        return "TrescaFrictionalLaw";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info() << std::endl;
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

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
    }

    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class TrescaFrictionalLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // TRESCA_FRICTIONAL_LAW_H_DEFINED  defined
