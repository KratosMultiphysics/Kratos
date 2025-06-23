// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_frictional_laws/frictional_law_with_derivative.h"

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
 * @class CoulombFrictionalLaw
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class defines the Coulomb frictional laws
 * @details This class defines the most typical friction model
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TFrictional If we are solving a frictional or frictionless problem
 * @tparam TNormalVariation If we are consider normal variation
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) CoulombFrictionalLaw
    : public FrictionalLawWithDerivative<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>
{
public:

    ///@name Type Definitions
    ///@{

    /// Base class definition
    using BaseType = FrictionalLawWithDerivative<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>;

    /// Definition of the derivative data
    using DerivativeDataType = typename BaseType::DerivativeDataType;

    /// The definition of the mortar operators
    using MortarConditionMatrices = typename BaseType::MortarConditionMatrices;

    /// Index type definition
    using IndexType= std::size_t;

    /// Size type definition
    using SizeType = std::size_t;

    /// Zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Counted pointer of CoulombFrictionalLaw
    KRATOS_CLASS_POINTER_DEFINITION( CoulombFrictionalLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    CoulombFrictionalLaw()
    {
    }

    ///Copy constructor  (not really required)
    CoulombFrictionalLaw(const CoulombFrictionalLaw& rhs)
    {
    }

    /// Destructor.
    ~CoulombFrictionalLaw()
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
        const Node& rNode,
        const PairedCondition& rCondition,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

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
    double GetDerivativeThresholdValue(
        const Node& rNode,
        const PairedCondition& rCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const DerivativeDataType& rDerivativeData,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const IndexType IndexDerivative,
        const IndexType IndexNode
        ) override;

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
        return "CoulombFrictionalLaw";
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
}; // Class CoulombFrictionalLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
