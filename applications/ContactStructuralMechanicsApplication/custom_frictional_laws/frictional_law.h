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

#if !defined(FRICTIONAL_LAW_H_DEFINED )
#define  FRICTIONAL_LAW_H_DEFINED

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/mortar_classes.h"

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
 * @class FrictionalLaw
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class defines the base class for frictional laws
 * @details This class does nothing, define derived frictional laws in order to make use of it
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TFrictional If we are solving a frictional or frictionless problem
 * @tparam TNormalVariation If we are consider normal variation
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// Node definition
    typedef Node<3> NodeType;

    /// Index type definition
    typedef std::size_t IndexType;

    /// Size type definition
    typedef std::size_t SizeType;

    /// Definition of the derivative data
    typedef DerivativeDataFrictional<TDim, TNumNodes, TNormalVariation, TNumNodesMaster> DerivativeDataType;

    /// Zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Counted pointer of FrictionalLaw
    KRATOS_CLASS_POINTER_DEFINITION( FrictionalLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    FrictionalLaw()
    {
    }

    ///Copy constructor  (not really required)
    FrictionalLaw(const FrictionalLaw& rhs)
    {
    }

    /// Destructor.
    ~FrictionalLaw()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns the friction coefficient
     * @param rNode The node where the friction coefficient is obtained
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     */
    virtual double GetFrictionCoefficient(
        const NodeType& rNode,
        const Condition& rCondition,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This method computes the threshold value considered for computing friction
     * @param rNode The node where the threshold value is obtained
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     */
    virtual double GetThresholdValue(
        const NodeType& rNode,
        const Condition& rCondition,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This method computes the threshold derivative value considered for computing friction
     * @param rNode The node where the threshold derivative value is obtained
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     * @param pDerivativeDatabase The pointer derivative database
     * @param IndexDerivative The derivative index
     */
    virtual double GetDerivativeThresholdValue(
        const NodeType& rNode,
        const Condition& rCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const DerivativeDataType& pDerivativeDatabase,
        const IndexType IndexDerivative
        );

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
    virtual std::string Info() const
    {
        return "FrictionalLaw";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class FrictionalLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // FRICTIONAL_LAW_H_DEFINED  defined
