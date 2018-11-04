//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#ifndef KRATOS_COMPRESSIBLE_ELEMENT_ROTATION_UTILITY
#define KRATOS_COMPRESSIBLE_ELEMENT_ROTATION_UTILITY

// system includes

// external includes

// kratos includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "utilities/coordinate_transformation_utilities.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
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

/// A utility to rotate the local contributions of certain nodes to the system matrix, which is required to apply slip conditions in arbitrary directions.
template<class TLocalMatrixType, class TLocalVectorType>
class CompressibleElementRotationUtility: public CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double> {
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of CompressibleElementRotationUtility
	KRATOS_CLASS_POINTER_DEFINITION(CompressibleElementRotationUtility);

	typedef Node<3> NodeType;

	typedef Geometry< Node<3> > GeometryType;

	///@}
	///@name Life Cycle
	///@{

	/// Constructor.
	/** @param DomainSize Number of space dimensions (2 or 3)
	 * @param NumRowsPerNode Number of matrix or vector rows associated to each node. Velocity DOFs are assumed to be the first mDomainSize rows in each block of rows.
	 * @param rVariable Kratos variable used to flag nodes where local system contributions will be rotated. All nodes with rVariable != Zero will be rotated.
	 * @param Zero The zero value for the variable.
	 */
	CompressibleElementRotationUtility(
        const unsigned int DomainSize,
		const Variable<double>& rVariable):
    CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double>(DomainSize,DomainSize+2,rVariable,0.0)
	{}

	/// Destructor.
	~CompressibleElementRotationUtility() override {}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	/// Rotate the local system contributions so that they are oriented with each node's normal.
	/**
	 @param rLocalMatrix Local system matrix
	 @param rLocalVector Local RHS vector
	 @param rGeometry A reference to the element's (or condition's) geometry
	 */
	void Rotate(
        TLocalMatrixType& rLocalMatrix,
		TLocalVectorType& rLocalVector,
		GeometryType& rGeometry) const override
	{
        if (this->GetDomainSize() == 2) this->template RotateAux<2,4,1>(rLocalMatrix,rLocalVector,rGeometry);
        else if (this->GetDomainSize() == 3) this->template RotateAux<3,5,1>(rLocalMatrix,rLocalVector,rGeometry);
	}

	/// RHS only version of Rotate
	void Rotate(
        TLocalVectorType& rLocalVector,
		GeometryType& rGeometry) const override
	{
        TLocalMatrixType dummy = ZeroMatrix(rLocalVector.size(),rLocalVector.size());
        if (this->GetDomainSize() == 2) this->template RotateAux<2,4,1>(dummy,rLocalVector,rGeometry);
        else if (this->GetDomainSize() == 3) this->template RotateAux<3,5,1>(dummy,rLocalVector,rGeometry);
	}

	/// Apply slip boundary conditions to the rotated local contributions.
	/** This function takes the local system contributions rotated so each
	 node's velocities are expressed using a base oriented with its normal
	 and imposes that the normal velocity is equal to the mesh velocity in
	 the normal direction.
	 */
	void ApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const override
	{
		const unsigned int LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

		if (LocalSize > 0)
		{
			for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
			{
				if(this->IsSlip(rGeometry[itNode]) )
				{
					// We fix the first momentum dof (normal component) for each rotated block
					unsigned int j = itNode * this->GetBlockSize() + 1; // +1 assumes DOF ordering as density-momentum-energy

					for( unsigned int i = 0; i < j; ++i)// Skip term (i,i)
					{
						rLocalMatrix(i,j) = 0.0;
						rLocalMatrix(j,i) = 0.0;
					}
					for( unsigned int i = j+1; i < LocalSize; ++i)
					{
						rLocalMatrix(i,j) = 0.0;
						rLocalMatrix(j,i) = 0.0;
					}

					rLocalVector(j) = 0.0;
					rLocalMatrix(j,j) = 1.0;
				}
			}
		}
	}

	/// RHS only version of ApplySlipCondition
	void ApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const override
	{
		if (rLocalVector.size() > 0)
		{
			for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
			{
				if( this->IsSlip(rGeometry[itNode]) )
				{
					// We fix the first momentum dof (normal component) for each rotated block
					unsigned int j = itNode * this->GetBlockSize() + 1; // +1 assumes DOF ordering as density-momentum-energy
					rLocalVector[j] = 0.0;
				}
			}
		}
	}

	/// Transform nodal velocities to the rotated coordinates (aligned with each node's normal)
	void RotateVelocities(ModelPart& rModelPart) const override
	{
		TLocalVectorType momentum(this->GetDomainSize());
		TLocalVectorType Tmp(this->GetDomainSize());

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
#pragma omp parallel for firstprivate(momentum,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			if( this->IsSlip(*itNode) )
			{
				//this->RotationOperator<TLocalMatrixType>(Rotation,);
				if(this->GetDomainSize() == 3)
				{
					BoundedMatrix<double,3,3> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rMomentum = itNode->FastGetSolutionStepValue(MOMENTUM);
					for(unsigned int i = 0; i < 3; i++) momentum[i] = rMomentum[i];
					noalias(Tmp) = prod(rRot,momentum);
					for(unsigned int i = 0; i < 3; i++) rMomentum[i] = Tmp[i];
				}
				else
				{
					BoundedMatrix<double,2,2> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rMomentum = itNode->FastGetSolutionStepValue(MOMENTUM);
					for(unsigned int i = 0; i < 2; i++) momentum[i] = rMomentum[i];
					noalias(Tmp) = prod(rRot,momentum);
					for(unsigned int i = 0; i < 2; i++) rMomentum[i] = Tmp[i];
				}
			}
		}
	}

	/// Transform nodal velocities from the rotated system to the original one
	void RecoverVelocities(ModelPart& rModelPart) const override
	{
		TLocalVectorType momentum(this->GetDomainSize());
		TLocalVectorType Tmp(this->GetDomainSize());

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
#pragma omp parallel for firstprivate(momentum,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			if( this->IsSlip(*itNode) )
			{
				if(this->GetDomainSize() == 3)
				{
					BoundedMatrix<double,3,3> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rMomentum = itNode->FastGetSolutionStepValue(MOMENTUM);
					for(unsigned int i = 0; i < 3; i++) momentum[i] = rMomentum[i];
					noalias(Tmp) = prod(trans(rRot),momentum);
					for(unsigned int i = 0; i < 3; i++) rMomentum[i] = Tmp[i];
				}
				else
				{
					BoundedMatrix<double,2,2> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rMomentum = itNode->FastGetSolutionStepValue(MOMENTUM);
					for(unsigned int i = 0; i < 2; i++) momentum[i] = rMomentum[i];
					noalias(Tmp) = prod(trans(rRot),momentum);
					for(unsigned int i = 0; i < 2; i++) rMomentum[i] = Tmp[i];
				}
			}
		}
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
	std::string Info() const override
	{
		std::stringstream buffer;
		buffer << "CompressibleElementRotationUtility";
		return buffer.str();
	}

	/// Print information about this object.
	void PrintInfo(std::ostream& rOStream) const override
	{
		rOStream << "CompressibleElementRotationUtility";
	}

	/// Print object's data.
	void PrintData(std::ostream& rOStream) const override {}

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
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	CompressibleElementRotationUtility& operator=(CompressibleElementRotationUtility const& rOther) {}

	/// Copy constructor.
	CompressibleElementRotationUtility(CompressibleElementRotationUtility const& rOther) {}

	///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TLocalMatrixType, class TLocalVectorType>
inline std::istream& operator >>(std::istream& rIStream,
		CompressibleElementRotationUtility<TLocalMatrixType, TLocalVectorType>& rThis) {
	return rIStream;
}

/// output stream function
template<class TLocalMatrixType, class TLocalVectorType>
inline std::ostream& operator <<(std::ostream& rOStream,
		const CompressibleElementRotationUtility<TLocalMatrixType, TLocalVectorType>& rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

///@}

///@} addtogroup block

}

#endif // KRATOS_COMPRESSIBLE_ELEMENT_ROTATION_UTILITY
