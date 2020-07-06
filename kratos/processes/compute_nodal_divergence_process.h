//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Miguel Maso Sotomayor
//

#ifndef KRATOS_COMPUTE_NODAL_DIVERGENCE_PROCESS_INCLUDED
#define KRATOS_COMPUTE_NODAL_DIVERGENCE_PROCESS_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/compute_nodal_normal_divergence_process.h"

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
 * @class ComputeNodalDivergenceProcess
 * @ingroup KratosCore
 * @brief Compute Nodal Divergence process
 * @details This process computes the divergence of a vector stored in the nodes
 * @author Miguel Maso Sotomayor
 * @tparam THistorical If the variable is historical or not
*/
template<bool THistorical>
class KRATOS_API(KRATOS_CORE) ComputeNodalDivergenceProcess
    : public ComputeNodalNormalDivergenceProcess<THistorical>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeNodalDivergenceProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNodalDivergenceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor. (component)
    ComputeNodalDivergenceProcess(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3> >& rOriginVariable,
        const Variable<double>& rDivergenceVariable,
        const Variable<double>& rAreaVariable = NODAL_AREA,
        const bool NonHistoricalOriginVariable = false
        );

    /// Destructor.
    ~ComputeNodalDivergenceProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

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
        return "ComputeNodalDivergenceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeNodalDivergenceProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
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

    /**
     * @brief This computes and nolmalizes the divergence
     */
    void ComputeDivergence(
        const Element::GeometryType& rThisGeometry,
        const Matrix& rDN_DX,
        const Variable<array_1d<double,3>>& rVariable,
        double& rDivergence
    ) override;

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

    ///@}

}; // Class ComputeNodalDivergenceProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_COMPUTE_NODAL_DIVERGENCE_PROCESS_INCLUDED  defined
