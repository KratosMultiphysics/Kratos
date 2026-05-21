//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_FLUX_LIMITER_H_INCLUDED
#define KRATOS_FLUX_LIMITER_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/node.h"
#include "includes/kratos_parameters.h"
#include "includes/global_pointer_variables.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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
 * @brief This is a helper class to separate the physics from the flux corrected scheme
 * @see FluxCorrectedShallowWaterScheme
*/
template<class TLocalVectorType>
class FluxLimiter
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FluxLimiter);

    struct AllowableIncrements {
        double max;
        double min;
        AllowableIncrements(double NewMax, double NewMin) {max = NewMax; min = NewMin;}
    };

    typedef Node NodeType;

    typedef GlobalPointersVector<Node> GlobalPointersVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    FluxLimiter() {
        AddFluxLimiters(this->GetDefaultParameters());
    }

    /**
     * @brief Constructor with parameters
     */
    FluxLimiter(Parameters ThisParameters) {
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
        AddFluxLimiters(ThisParameters);
    }

    /**
     * @brief Copy constructor
     */
    FluxLimiter(FluxLimiter const& rOther) {
        mUnlimitedFlux = rOther.mUnlimitedFlux;
        mAllowableIncrements = rOther.mAllowableIncrements;
    }

    /**
     * @brief Destructor
     */
    virtual ~FluxLimiter() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    IndexType size() {return mUnlimitedFlux.size();}

    double ComputeUnlimitedFlux(const TLocalVectorType& rContributions, const NodeType& rNode, IndexType iNode, IndexType iVar) {
        return mUnlimitedFlux[iVar](rContributions, rNode, iNode);
    }

    AllowableIncrements ComputeAllowableIncrements(const NodeType& rNode, IndexType iVar) {
        return mAllowableIncrements[iVar](rNode);
    }

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
        std::stringstream buffer;
        buffer << "FluxLimiter" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "FluxLimiter";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<std::function<double(const TLocalVectorType&, const NodeType&, IndexType)>> mUnlimitedFlux;

    std::vector<std::function<AllowableIncrements(const NodeType&)>> mAllowableIncrements;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void AddFluxLimiters(Parameters ThisParameters) {
        for (auto var_name : ThisParameters["limiting_variables"].GetStringArray()) {
            if (KratosComponents<Variable<double>>::Has(var_name)) {
                const auto& r_var = KratosComponents<Variable<double>>::Get(var_name);
                this->SetDoubleLimiter(r_var);
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(var_name)) {
                const auto& r_var = KratosComponents<Variable<array_1d<double,3>>>::Get(var_name);
                this->SetArrayLimiter(r_var);
            }
        }
    }

    virtual void SetDoubleLimiter(const Variable<double>& rVariable)
    {
        if (rVariable == HEIGHT) {
            mUnlimitedFlux.push_back(HeightUnlimitedFlux);
            mAllowableIncrements.push_back(HeightAllowableIncrements);
        } else if (rVariable == FREE_SURFACE_ELEVATION) {
            mUnlimitedFlux.push_back(FreeSurfaceUnlimitedFlux);
            mAllowableIncrements.push_back(FreeSurfaceAllowableIncrements);
        } else {
            KRATOS_ERROR << "FluxLimiter: The limiter for the variable " << rVariable << " is undefined." << std::endl;
        }
    }

    virtual void SetArrayLimiter(const Variable<array_1d<double,3>>& rVariable)
    {
        if (rVariable == MOMENTUM) {
            mUnlimitedFlux.push_back(FlowRateUnlimitedFlux);
            mAllowableIncrements.push_back(FlowRateAllowableIncrements);
        } else if (rVariable == VELOCITY) {
            mUnlimitedFlux.push_back(VelocityUnlimitedFlux);
            mAllowableIncrements.push_back(VelocityAllowableIncrements);
        } else {
            KRATOS_ERROR << "FluxLimiter: The limiter for the variable " << rVariable << " is undefined." << std::endl;
        }
    }

    static double HeightUnlimitedFlux(const TLocalVectorType& rContributions, const NodeType& rNode, IndexType iNode)
    {
        return rContributions(3*iNode + 2);
    }

    static double FreeSurfaceUnlimitedFlux(const TLocalVectorType& rContributions, const NodeType& rNode, IndexType iNode)
    {
        return rContributions(3*iNode + 2);
    }

    static double FlowRateUnlimitedFlux(const TLocalVectorType& rContributions, const NodeType& rNode, IndexType iNode)
    {
        array_1d<double,3> flux;
        flux[0] = rContributions(3*iNode);
        flux[1] = rContributions(3*iNode + 1);
        flux[2] = 0.0;
        const array_1d<double,3>& v = rNode.FastGetSolutionStepValue(VELOCITY);
        return inner_prod(flux, v); // projection onto velocity
    }

    static double VelocityUnlimitedFlux(const TLocalVectorType& rContributions, const NodeType& rNode, IndexType iNode)
    {
        array_1d<double,3> flux_q, flux_v;
        flux_q[0] = rContributions(3*iNode);
        flux_q[1] = rContributions(3*iNode + 1);
        flux_q[2] = 0.0;
        double flux_h = rContributions(3*iNode + 2);
        const double h = rNode.FastGetSolutionStepValue(HEIGHT);
        const array_1d<double,3>& v = rNode.FastGetSolutionStepValue(VELOCITY);
        const double next_h = h + flux_h;
        flux_v = flux_q / next_h - v * flux_h / next_h;
        return inner_prod(flux_v, v); // projection onto velocity
    }

    static AllowableIncrements HeightAllowableIncrements(const NodeType& rNode)
    {
        // Getting the maximum increments for the height variable
        AllowableIncrements increments(0.0, 0.0);
        const auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        const double h_i = rNode.FastGetSolutionStepValue(HEIGHT);
        for (IndexType j = 0; j < neigh_nodes.size(); ++j) {
            double delta_ij = neigh_nodes[j].FastGetSolutionStepValue(HEIGHT) - h_i;
            increments.max = std::max(increments.max, delta_ij);
            increments.min = std::min(increments.min, delta_ij);
        }
        return increments;
    }

    static AllowableIncrements FreeSurfaceAllowableIncrements(const NodeType& rNode)
    {
        // Getting the maximum increments for the free surface variable
        AllowableIncrements increments(0.0, 0.0);
        const auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        const double h_i = rNode.FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
        for (IndexType j = 0; j < neigh_nodes.size(); ++j) {
            double delta_ij = neigh_nodes[j].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - h_i;
            increments.max = std::max(increments.max, delta_ij);
            increments.min = std::min(increments.min, delta_ij);
        }
        return increments;
    }

    static AllowableIncrements FlowRateAllowableIncrements(const NodeType& rNode)
    {
        // Getting the maximum increments for the flow rate variable
        // The increments are projected to convert it to a scalar
        AllowableIncrements increments(0.0, 0.0);
        const auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        const array_1d<double,3>& q_i = rNode.FastGetSolutionStepValue(MOMENTUM);
        const array_1d<double,3>& v_i = rNode.FastGetSolutionStepValue(VELOCITY);
        for (IndexType j = 0; j < neigh_nodes.size(); ++j) {
            const array_1d<double,3>& delta_ij = neigh_nodes[j].FastGetSolutionStepValue(MOMENTUM) - q_i;
            double proj_ij = inner_prod(v_i, delta_ij);
            increments.max = std::max(increments.max, proj_ij);
            increments.min = std::min(increments.min, proj_ij);
        }
        return increments;
    }

    static AllowableIncrements VelocityAllowableIncrements(const NodeType& rNode)
    {
        // Getting the maximum increments for the velocity variable
        // The increments are projected to convert it to a scalar
        AllowableIncrements increments(0.0, 0.0);
        const auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        const array_1d<double,3>& v_i = rNode.FastGetSolutionStepValue(VELOCITY);
        for (IndexType j = 0; j < neigh_nodes.size(); ++j) {
            const array_1d<double,3>& delta_ij = neigh_nodes[j].FastGetSolutionStepValue(VELOCITY) - v_i;
            double proj_ij = inner_prod(v_i, delta_ij);
            increments.max = std::max(increments.max, proj_ij);
            increments.min = std::min(increments.min, proj_ij);
        }
        return increments;
    }

    /**
     * @brief This method provides the defaults parameters
     * @return The default parameters
     */
    virtual Parameters GetDefaultParameters() const
    {
        Parameters default_parameters = Parameters(R"({
            "limiting_variables" : ["HEIGHT"]
        })");
        return default_parameters;
    }

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
    FluxLimiter& operator=(FluxLimiter const& rOther){}

    ///@}

}; // Class FluxLimiter

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FLUX_LIMITER_H_INCLUDED  defined
