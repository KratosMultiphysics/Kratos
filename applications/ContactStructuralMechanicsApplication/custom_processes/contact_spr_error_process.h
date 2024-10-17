// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Anna Rehr
//  Co-author   :    Vicente Mataix Ferrandiz
//

#pragma once

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
    using SizeType = std::size_t;

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

    /// Basetype definition
    using BaseType = SPRErrorProcess<TDim>;

    /// Containers definition
    using NodesArrayType = ModelPart::NodesContainerType;
    using ElementsArrayType = ModelPart::ElementsContainerType;
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    /// Definition of the iterators
    using WeakElementItType = GlobalPointersVector<Element>::iterator;
    using NodeItType = NodesArrayType::iterator;
    using ElementItType = ElementsArrayType::iterator;

    /// Definition of the indextype
    using IndexType = std::size_t;

    /// Pointer definition of ContactSPRErrorProcess
    KRATOS_CLASS_POINTER_DEFINITION(ContactSPRErrorProcess);

    /// The Voigt notation size
    static constexpr SizeType SigmaSize = (TDim == 2) ? 3 : 6;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */
    ContactSPRErrorProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~ContactSPRErrorProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        this->Execute();
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

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
    std::string Info() const override
    {
        return "ContactSPRErrorProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ContactSPRErrorProcess";
    }

    /// Print object"s data.
    void PrintData(std::ostream& rOStream) const override
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
        const SizeType NeighbourSize,
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

    /**
     * @brief This method computes the tangents and normal matrices
     * @param rNk Normal matrix
     * @param rTk1 First tangent matrix
     * @param rTk2 Second tangent matrix
     * @param rNormal The normal vector
     */
    void ComputeNormalTangentMatrices(
        BoundedMatrix<double, 1, SigmaSize>& rNk,
        BoundedMatrix<double, 1, SigmaSize>& rTk1,
        BoundedMatrix<double, 1, SigmaSize>& rTk2,
        const array_1d<double, 3>& rNormal
        );

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
