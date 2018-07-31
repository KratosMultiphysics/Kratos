//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Anna Rehr
//  Co-author   :    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SPR_ERROR_METRICS_PROCESS)
#define KRATOS_SPR_ERROR_METRICS_PROCESS

// System includes
#include <omp.h>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "structural_mechanics_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"
#include "spaces/ublas_space.h"
#include "utilities/math_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
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

    /// Definition of the size type
    typedef std::size_t                                                             SizeType;
    
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
 * @class SPRMetricProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This class is can be used to compute the metrics of the model part with a superconvergent patch recovery approach
 * @details The formulation employed in order to compute the super patch recovery is based on the work of O. C. Zienkiewicz
J. Z. Zhu, and extended for contact mechanics. In the papers:
 * - The superconvergent patch recovery and a posteriori error estimates. Part 1: The recovery technique https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620330702
 * - The superconvergent patch recovery and a posteriori error estimates. Part 2: Error estimates and adaptivity https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620330703
 * This is a general recovery technique is developed for determining the derivatives (stresses) of the finite element solutions at nodes. The implementation of the recovery technique is simple and cost effective.
 * @tparam TDim The dimension to be computed
 * @author Anna Rehr
 */
template<SizeType TDim>
class SPRMetricProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of SPRMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(SPRMetricProcess);
    
    /// The Voigt notation size
    static constexpr SizeType SigmaSize = (TDim == 2) ? 3 : 6;

    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */
    
    SPRMetricProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );
    
    /// Destructor.
    virtual ~SPRMetricProcess() {}
    
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
     * We initialize the metrics of the MMG sol using the Hessian metric matrix approach
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
        return "SPRMetricProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SPRMetricProcess";
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
    
    ModelPart& mThisModelPart;                /// The model part to compute

    double mMinSize;                          /// The minimal size of the elements
    double mMaxSize;                          /// The maximal size of the elements

    double mPenaltyNormal;                    /// The normal penalty
    double mPenaltyTangent;                   /// The tangent penalty

    SizeType mEchoLevel;                      /// The echo level

    bool mSetElementNumber;                   /// Determines if a target number of elements for the new mesh is set
    SizeType mElementNumber;                  /// The target number of elements for the new mesh
    double mTargetError;                      /// The overall target error for the new mesh
    bool mAverageNodalH;                      /// Determines if the nodal h is averaged from the surrounding elements or if the lowest value is taken
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the superconvergent stresses
     */
    void CalculateSuperconvergentStresses();

    /**
     * @brief This method estimates the error and computes the new element size
     * @param rEnergyNormOverall The mean of the energy norm
     * @param rErrorOverall The mean of the error
     */
    void CalculateErrorEstimationAndElementSize(
        double& rEnergyNormOverall,
        double& rErrorOverall
        );

    /**
     * @brief In this final step the metric is computed
     */
    void CalculateMetric();

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
        );

    /**
     * @brief Calculates the recovered stress at a node in the case of a standard patch without contact BC
     * @param itNode the node for which the recovered stress should be calculated
     * @param itPatchNode the center node of the patch
     * @param NeighbourSize Number of neighbour elements
     * @param rSigmaRecovered The recovered stress
     */
    void CalculatePatchStandard(
        NodeItType itNode,
        NodeItType itPatchNode,
        SizeType NeighbourSize,
        Vector& rSigmaRecovered
        );

    /**
     * @brief It calculates the recovered stress at a node where contact BCs are regarded
     * @param itNode the node for which the recovered stress should be calculated
     * @param itPatchNode the center node of the patch
     * @param NeighbourSize Number of neighbour elements
     * @param rSigmaRecovered The recovered stress
     */
    void CalculatePatchContact(
        NodeItType itNode,
        NodeItType itPatchNode,
        SizeType NeighbourSize,
        Vector& rSigmaRecovered
        );

    /**
     * @brief This computes the element size depending of the geometry and it assigns to the ELEMENT_H variable
     * @param itElement The element iterator
     */
    void ComputeElementSize(ElementItType itElement);
    
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
    SPRMetricProcess& operator=(SPRMetricProcess const& rOther)
    {
        return *this;
    };

    /// Copy constructor.
    //SPRMetricProcess(SPRMetricProcess const& rOther);

};// class SPRMetricProcess

};// namespace Kratos.
#endif /* KRATOS_SPR_ERROR_METRICS_PROCESS defined */
