// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Anna Rehr
//  Co-author   :    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CONTACT_SPR_ERROR_PROCESS)
#define KRATOS_CONTACT_SPR_ERROR_PROCESS

// System includes

// External includes

// Project includes
#include "custom_processes/spr_error_process.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// Definition of the size type
    typedef std::size_t SizeType;

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
 * @class ContactSPRErrorProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class is can be used to compute the metrics of the model part with a superconvergent patch recovery (SPR) approach
 * @details The formulation employed in order to compute the super patch recovery is based on the work of O. C. Zienkiewicz
J. Z. Zhu, and extended for contact mechanics. In the papers:
 * - The superconvergent patch recovery and a posteriori error estimates. Part 1: The recovery technique https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620330702
 * - The superconvergent patch recovery and a posteriori error estimates. Part 2: Error estimates and adaptivity https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620330703
 * This is a general recovery technique is developed for determining the derivatives (stresses) of the finite element solutions at nodes. The implementation of the recovery technique is simple and cost effective.
 * @tparam TDim The dimension to be computed
 * @author Anna Rehr
 */
template<SizeType TDim>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ContactSPRErrorProcess
    : public SPRErrorProcess<TDim>
{
public:

    ///@name Type Definitions
    ///@{

    // Basetype definition
    typedef SPRErrorProcess<TDim>                                                   BaseType;

    /// Containers definition
    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;

    /// The definition of the node type
    typedef Node <3>                                                                NodeType;

    /// Definition of the iterators
    typedef WeakPointerVector< Element >::iterator                         WeakElementItType;
    typedef NodesArrayType::iterator                                              NodeItType;
    typedef ElementsArrayType::iterator                                        ElementItType;

    /// Definition of the indextype
    typedef std::size_t                                                            IndexType;

    /// Pointer definition of ContactSPRErrorProcess
    KRATOS_CLASS_POINTER_DEFINITION(ContactSPRErrorProcess);

    /// The Voigt notation size
    static constexpr SizeType SigmaSize = (TDim == 2) ? 3 : 6;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */

    ContactSPRErrorProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    virtual ~ContactSPRErrorProcess() {}

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        BaseType::Execute();
    }

    ///@}
    ///@name Operations
    ///@{

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
    virtual std::string Info() const override
    {
        return "ContactSPRErrorProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ContactSPRErrorProcess";
    }

    /// Print object"s data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

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

    /**
     * @brief Calculates the recovered stress. Checks whatever this is a contact case or a standard one
     * @param itNode the node for which the recovered stress should be calculated
     * @param itPatchNode the center node of the patch
     * @param NeighbourSize Number of neighbour elements
     * @param rSigmaRecovered The recovered stress
     */
    void CalculatePatch(
        NodeItType itNode,
        NodeItType itPatchNode,
        SizeType NeighbourSize,
        Vector& rSigmaRecovered
        ) override;

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
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    double mPenaltyNormal;   /// The normal penalty
    double mPenaltyTangent;  /// The tangent penalty

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
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ContactSPRErrorProcess& operator=(ContactSPRErrorProcess const& rOther)
    {
        return *this;
    };

    /// Copy constructor.
    //ContactSPRErrorProcess(ContactSPRErrorProcess const& rOther);

};// class ContactSPRErrorProcess

};// namespace Kratos.
#endif /* KRATOS_CONTACT_SPR_ERROR_PROCESS defined */
