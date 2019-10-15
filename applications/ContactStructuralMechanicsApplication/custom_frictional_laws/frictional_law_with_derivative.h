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

#if !defined(FRICTIONAL_LAW_WITH_DERIVATIVE_H_DEFINED )
#define  FRICTIONAL_LAW_WITH_DERIVATIVE_H_DEFINED

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
 * @class FrictionalLawWithDerivative
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class defines the base class for frictional laws with derivative
 * @details This class does nothing, define derived frictional laws in order to make use of it
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TFrictional If we are solving a frictional or frictionless problem
 * @tparam TNormalVariation If we are consider normal variation
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative
    : public FrictionalLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// Define the base class
    typedef FrictionalLaw BaseType;

    /// Node definition
    typedef Node<3> NodeType;

    /// Index type definition
    typedef std::size_t IndexType;

    /// Size type definition
    typedef std::size_t SizeType;

    /// Definition of the derivative data
    typedef DerivativeDataFrictional<TDim, TNumNodes, TNormalVariation, TNumNodesMaster> DerivativeDataType;

    /// The definition of the mortar operators
    typedef MortarOperatorWithDerivatives<TDim, TNumNodes, true, TNormalVariation, TNumNodesMaster> MortarConditionMatrices;

    /// Zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Counted pointer of FrictionalLawWithDerivative
    KRATOS_CLASS_POINTER_DEFINITION( FrictionalLawWithDerivative );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    FrictionalLawWithDerivative()
    {
    }

    ///Copy constructor  (not really required)
    FrictionalLawWithDerivative(const FrictionalLawWithDerivative& rhs)
    {
    }

    /// Destructor.
    ~FrictionalLawWithDerivative()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the threshold derivative value considered for computing friction
     * @param rNode The node where the threshold derivative value is obtained
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     * @param rDerivativeData The reference to the derivative database
     * @param rMortarConditionMatrices The container of the mortar operators
     * @param IndexDerivative The derivative index
     * @param IndexNode The corresponding node index on the condition geometry
     */
    virtual double GetDerivativeThresholdValue(
        const NodeType& rNode,
        const PairedCondition& rCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const DerivativeDataType& rDerivativeData,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const IndexType IndexDerivative,
        const IndexType IndexNode
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
    std::string Info() const override
    {
        return "FrictionalLawWithDerivative";
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
}; // Class FrictionalLawWithDerivative

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // FRICTIONAL_LAW_WITH_DERIVATIVE_H_DEFINED  defined
