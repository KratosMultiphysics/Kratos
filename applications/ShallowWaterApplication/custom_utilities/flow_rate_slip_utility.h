//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_FLOW_RATE_SLIP_UTILITY_H_INCLUDED
#define KRATOS_FLOW_RATE_SLIP_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/coordinate_transformation_utilities.h"

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
 * @brief Tools to apply slip conditions
 * @detail A utility to rotate the local contributions of certain nodes to the system matrix, which is required to apply slip conditions in arbitrary directions.
*/
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
class FlowRateSlipUtility
    : public CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,TValueType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FlowRateSlipUtility
    KRATOS_CLASS_POINTER_DEFINITION(FlowRateSlipUtility);

    typedef CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,TValueType> BaseType;

    typedef std::size_t SizeType;

	typedef Node<3> NodeType;

	typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FlowRateSlipUtility() : BaseType(2,3,SLIP) {}

    /// Destructor.
    virtual ~FlowRateSlipUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

	/** 
     * @brief Apply slip boundary conditions to the rotated local contributions.
     * @detail This function takes the local system contributions rotated so each
	 * node's velocities are expressed using a base oriented with its normal
	 * and imposes that the normal velocity is equal to the mesh velocity in
	 * the normal direction.
     * @param rLocalMatrix A reference to the LHS local matrix
     * @param rLocalVector A reference to the RHS local vector
     * @param rGeometry A reference to the geometry of the element or condition
	 */
	virtual void ApplySlipCondition(
        TLocalMatrixType& rLocalMatrix,
		TLocalVectorType& rLocalVector,
		GeometryType& rGeometry) const override
    {
        const SizeType LocalSize = rLocalVector.size(); // We expect this to work both with elements and conditions

        if (LocalSize > 0)
        {
            for (SizeType it_node = 0; it_node < rGeometry.PointsNumber(); ++it_node)
            {
                if (this->IsSlip(rGeometry[it_node]))
                {
                    // We fix the first dof (normal velocity) for each rotated block
                    SizeType j = it_node * BaseType::GetBlockSize();

                    array_1d<double,3> vel = rGeometry[it_node].FastGetSolutionStepValue(MOMENTUM);
                    array_1d<double,3> n = rGeometry[it_node].FastGetSolutionStepValue(NORMAL);
                    this->Normalize(n);

                    for (SizeType i = 0; i < j; ++i) // Skip term (i,i)
                    {
                        rLocalMatrix(i,j) = 0.0;
                        rLocalMatrix(j,i) = 0.0;
                    }
                    for (SizeType i = j+1; i < LocalSize; ++i)
                    {
                        rLocalMatrix(i,j) = 0.0;
                        rLocalMatrix(j,i) = 0.0;
                    }

                    rLocalVector(j) = - inner_prod(n, vel);
                    rLocalMatrix(j,j) = 1.0;
                }
            }
        }
    }

	/** 
     * @brief RHS only version of ApplySlipCondition
     * @param rLocalVector A reference to the RHS local vector
     * @param rGeometry A reference to the geometry of the element or condition
	 */
	virtual void ApplySlipCondition(
        TLocalVectorType& rLocalVector,
		GeometryType& rGeometry) const override
    {
        if (rLocalVector.size() > 0)
        {
            for (SizeType it_node = 0; it_node < rGeometry.PointsNumber(); ++it_node)
            {
                if (this->IsSlip(rGeometry[it_node]))
                {
                    // We fix the first dof (normal velocity) for each rotated block
                    SizeType j = it_node * BaseType::GetBlockSize();

                    array_1d<double,3> vel = rGeometry[it_node].FastGetSolutionStepValue(MOMENTUM);
                    array_1d<double,3> n = rGeometry[it_node].FastGetSolutionStepValue(NORMAL);
                    this->Normalize(n);

                    rLocalVector[j] = inner_prod(n, vel);
                }
            }
        }
    }

	/**
     * @brief Transform nodal velocities to the rotated coordinates (aligned with each node's normal)
     * @param rModelPart A reference to the model part
     * @see RecoverVelocities
     */
	virtual void RotateVelocities(ModelPart& rModelPart) const override
    {
        struct TLS {
            TLocalVectorType vel;
            TLocalVectorType tmp;
        };
        TLS tls;
        tls.vel.resize(BaseType::GetDomainSize());
        tls.tmp.resize(BaseType::GetDomainSize());

        block_for_each(rModelPart.Nodes(), tls, [&](NodeType& rNode, TLS& rTLS){
            if (this->IsSlip(rNode))
            {
                // For shallow water problems, domain size is always 2
                BoundedMatrix<double,2,2> rot;
                BaseType::LocalRotationOperatorPure(rot, rNode);

                array_1d<double,3>& r_velocity = rNode.FastGetSolutionStepValue(MOMENTUM);
                for(SizeType i = 0; i < 2; i++) {
                    rTLS.vel[i] = r_velocity[i];
                }
                noalias(rTLS.tmp) = prod(rot, rTLS.vel);
                for(SizeType i = 0; i < 2; i++) {
                    r_velocity[i] = rTLS.tmp[i];
                }
            }
        });
    }

	/**
     * Transform nodal velocities from the rotated system to the original one
     * @param rModelPart A reference to the model part
     * @see RotateVelocities
     */
	virtual void RecoverVelocities(ModelPart& rModelPart) const override
    {
        struct TLS {
            TLocalVectorType vel;
            TLocalVectorType tmp;
        };
        TLS tls;
        tls.vel.resize(BaseType::GetDomainSize());
        tls.tmp.resize(BaseType::GetDomainSize());

        block_for_each(rModelPart.Nodes(), tls, [&](NodeType& rNode, TLS& rTLS){
            if( this->IsSlip(rNode) )
            {
                // For shallow water problems, domain size is always 2
                BoundedMatrix<double,2,2> rot;
                BaseType::LocalRotationOperatorPure(rot, rNode);

                array_1d<double,3>& r_velocity = rNode.FastGetSolutionStepValue(MOMENTUM);
                for(SizeType i = 0; i < 2; i++) {
                    rTLS.vel[i] = r_velocity[i];
                }
                noalias(rTLS.tmp) = prod(trans(rot), rTLS.vel);
                for(SizeType i = 0; i < 2; i++) {
                    r_velocity[i] = rTLS.tmp[i];
                }
            }
        });
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

	/**
     * Turn back information as a string.
     */
	virtual std::string Info() const override
	{
		std::stringstream buffer;
		buffer << "FlowRateSlipUtility";
		return buffer.str();
	}

	/**
     * Print information about this object.
     */
	virtual void PrintInfo(std::ostream& rOStream) const override
	{
		rOStream << "FlowRateSlipUtility";
	}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    // FlowRateSlipUtility& operator=(FlowRateSlipUtility const& rOther) {}

    /// Copy constructor.
    // FlowRateSlipUtility(FlowRateSlipUtility const& rOther) {}

    ///@}

}; // Class FlowRateSlipUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
inline std::istream& operator >> (
    std::istream& rIStream,
    FlowRateSlipUtility<TLocalMatrixType, TLocalVectorType,TValueType>& rThis)
{
	return rIStream;
}

/// output stream function
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const FlowRateSlipUtility<TLocalMatrixType, TLocalVectorType,TValueType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FLOW_RATE_SLIP_UTILITY_H_INCLUDED  defined
