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

#if !defined(KRATOS_SPR_ERROR_PROCESS)
#define KRATOS_SPR_ERROR_PROCESS

// System includes
#ifdef _OPENMP
#include <omp.h>
#endif

// External includes

// Project includes
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "structural_mechanics_application_variables.h"
#include "spaces/ublas_space.h"
#include "utilities/math_utils.h"

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
 * @class SPRErrorProcess
 * @ingroup StructuralMechanicsApplication
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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SPRErrorProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{

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

    /// Pointer definition of SPRErrorProcess
    KRATOS_CLASS_POINTER_DEFINITION(SPRErrorProcess);

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

    SPRErrorProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    virtual ~SPRErrorProcess() {}

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief We initialize the metrics of the MMG sol using the Hessian metric matrix approach
     */
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
    virtual std::string Info() const override
    {
        return "SPRErrorProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SPRErrorProcess";
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

    ModelPart& mThisModelPart;                               /// The model part to compute
    Variable<Vector> mStressVariable = CAUCHY_STRESS_VECTOR; /// The stress variable considered
    SizeType mEchoLevel;                                     /// The echo level

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the superconvergent stresses
     */
    void CalculateSuperconvergentStresses();

    /**
     * @brief This method estimates the error
     * @param rEnergyNormOverall The mean of the energy norm
     * @param rErrorOverall The mean of the error
     */
    void CalculateErrorEstimation(
        double& rEnergyNormOverall,
        double& rErrorOverall
        );

    /**
     * @brief Calculates the recovered stress. Checks whatever this is a contact case or a standard one
     * @param itNode the node for which the recovered stress should be calculated
     * @param itPatchNode the center node of the patch
     * @param NeighbourSize Number of neighbour elements
     * @param rSigmaRecovered The recovered stress
     */
    virtual void CalculatePatch(
        NodeItType itNode,
        NodeItType itPatchNode,
        SizeType NeighbourSize,
        Vector& rSigmaRecovered
        );

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method performs a neighnour search
     * @param rModelPart The model part where to search
     */
    static inline void FindNodalNeighbours(ModelPart& rModelPart);

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
    SPRErrorProcess& operator=(SPRErrorProcess const& rOther)
    {
        return *this;
    };

    /// Copy constructor.
    //SPRErrorProcess(SPRErrorProcess const& rOther);

};// class SPRErrorProcess

};// namespace Kratos.
#endif /* KRATOS_SPR_ERROR_PROCESS defined */
