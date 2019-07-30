//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_NONCONFORMING_BOUNDARY_ROTATION_UTILITY
#define KRATOS_MPM_NONCONFORMING_BOUNDARY_ROTATION_UTILITY

// system includes

// external includes

// kratos includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/condition.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos {

///@addtogroup ParticleMechanicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

#if !defined(ROTATE_RECOVER)
#define ROTATE_RECOVER
    enum RotationConfiguration {RotateSystem = 0, RecoverSystem = 1};
#endif

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/* A utility to rotate the local contributions of NONCONFORMING conditions to the system matrix,
which is required to apply slip conditions (roller-type support) in arbitrary directions to the NONCONFORMING boundary conditions.*/
template<class TLocalMatrixType, class TLocalVectorType>
class MPMNonconformingBoundaryRotationUtility
{
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of MPMNonconformingBoundaryRotationUtility
	KRATOS_CLASS_POINTER_DEFINITION(MPMNonconformingBoundaryRotationUtility);

	typedef Node<3> NodeType;

	typedef Geometry< Node<3> > GeometryType;

	typedef Condition::Pointer ConditionPointerType;

	///@}
	///@name Life Cycle
	///@{

	/// Constructor.
	/** @param DomainSize Number of space dimensions (2 or 3)
	 * @param BlockSize Number of matrix or vector rows associated to each node. The BlockSize is used to adapt rotation well for various nonconforming boundary conditions
	 * 		  which utilize additional variables.
	 * @param rVariable Kratos variable used to flag nodes where local system contributions will be rotated. All nodes with rVariable != Zero will be rotated.
	 */
	MPMNonconformingBoundaryRotationUtility(
        const unsigned int DomainSize,
		const unsigned int BlockSize,
		const Kratos::Flags& rSelectionFlag = SLIP,
		const Variable<int>& rBoundaryType = MPC_BOUNDARY_CONDITION_TYPE):
	mDomainSize(DomainSize),
	mBlockSize(BlockSize),
	mrFlag(rSelectionFlag),
	mrFlagVariable(rBoundaryType)
	{}

	/// Destructor.
	virtual ~MPMNonconformingBoundaryRotationUtility() {}

	/// Assignment operator.
	MPMNonconformingBoundaryRotationUtility& operator=(MPMNonconformingBoundaryRotationUtility const& rOther) {}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	/// Rotate the local system contributions so that they are oriented with the condition's normal.
	/**
	 @param rLocalMatrix Local system matrix
	 @param rLocalVector Local RHS vector
	 @param pCondition A Nonconforming Boundary Condition
	 */
	virtual void Rotate(
        TLocalMatrixType& rLocalMatrix,
		TLocalVectorType& rLocalVector,
		ConditionPointerType pCondition) const
	{
		if (this->GetBlockSize() == this->GetDomainSize()) // irreducible case
		{
			if (this->GetDomainSize() == 2) RotateAndRecoverAuxPure<2>(rLocalMatrix,rLocalVector,pCondition);
			else if (this->GetDomainSize() == 3) RotateAndRecoverAuxPure<3>(rLocalMatrix,rLocalVector,pCondition);
		}
		else // cases with larger blocksize
		{
			KRATOS_ERROR << "Rotating Nonconforming conditions with BlockSize != DomainSize is not yet possible." << std::endl;
		}
	}

	/// RHS only version of Rotate
	virtual void Rotate(
        TLocalVectorType& rLocalVector,
		ConditionPointerType pCondition) const
	{
		if(this->GetDomainSize() == 2) RotateAndRecoverAuxPureRHS<2>(rLocalVector,pCondition);
		else if(this->GetDomainSize() == 3) RotateAndRecoverAuxPureRHS<3>(rLocalVector,pCondition);
	}


	/// Recover a rotated local system contributions to the original position.
	/**
	 @param rLocalMatrix Rotated local system matrix
	 @param rLocalVector Rotated local RHS vector
	 @param pCondition A Nonconforming Boundary Condition
	 */
	virtual void Recover(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			ConditionPointerType pCondition) const
	{
		if(this->GetBlockSize() == this->GetDomainSize())
		{
			if(this->GetDomainSize() == 2) RotateAndRecoverAuxPure<2>(rLocalMatrix,rLocalVector,pCondition,RecoverSystem);
			else if(this->GetDomainSize() == 3) RotateAndRecoverAuxPure<3>(rLocalMatrix,rLocalVector,pCondition,RecoverSystem);
		}
		else // cases with larger blocksize
		{
			KRATOS_ERROR << "Rotating Nonconforming conditions with BlockSize != DomainSize is not yet possible." << std::endl;
		}
	}

	/// RHS only version of Recover
	virtual void Recover(TLocalVectorType& rLocalVector,
			ConditionPointerType pCondition) const
	{
		if(this->GetDomainSize() == 2) RotateAndRecoverAuxPureRHS<2>(rLocalVector,pCondition,RecoverSystem);
		else if(this->GetDomainSize() == 3) RotateAndRecoverAuxPureRHS<3>(rLocalVector,pCondition,RecoverSystem);
	}

	/// Apply roler type boundary conditions to the rotated local contributions.
	/** This function takes the rotated local system contributions and impose a slip condition
		according to the different boundary types. Boundary type available are:
		1. Penalty Condition
	 */
	void ApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			ConditionPointerType pCondition) const
	{
		// 1. Penalty Condition
		if(this->GetBoundaryType(pCondition) == 1){
			this->ApplySlipPenaltyCondition(rLocalMatrix, rLocalVector, pCondition);
		}
		else // Boundary Type is not present
		{
			KRATOS_ERROR << "BoundaryType code number: " << this->GetBoundaryType(pCondition) << " is not present nor determined." << std::endl;
		}
	}

	/// RHS only version of ApplySlipCondition
	void ApplySlipCondition(TLocalVectorType& rLocalVector,
			ConditionPointerType pCondition) const
	{
		// 1. Penalty Condition
		if(this->GetBoundaryType(pCondition) == 1){
			this->ApplySlipPenaltyCondition(rLocalVector, pCondition);
		}
		else // Boundary Type is not present
		{
			KRATOS_ERROR << "BoundaryType code number: " << this->GetBoundaryType(pCondition) << " is not present nor determined." << std::endl;
		}
	}

	// A function to impose slip to penalty boundary condition
	void ApplySlipPenaltyCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			ConditionPointerType pCondition) const
	{
		const unsigned int local_size = rLocalVector.size();
		auto r_geometry = pCondition->GetGeometry();
		if (local_size > 0)
		{
			const unsigned int block_size = this->GetBlockSize();
			TLocalMatrixType temp_matrix = ZeroMatrix(rLocalMatrix.size1(),rLocalMatrix.size2());
			for(unsigned int itNode = 0; itNode < r_geometry.PointsNumber(); ++itNode)
			{
				if(this->IsSlip(pCondition))
				{
					// We fix the first displacement dof (normal component) for each rotated block
					unsigned int j = itNode * block_size;

					// Copy all normal value in LHS to the temp_matrix
					for (unsigned int i = j; i < rLocalMatrix.size1(); i+= block_size)
					{
						temp_matrix(i,j) = rLocalMatrix(i,j);
						temp_matrix(j,i) = rLocalMatrix(j,i);
					}

					// Remove all other value in RHS than the normal component
					for(unsigned int i = j; i < (j + block_size); ++i)
					{
						if (i!=j) rLocalVector[i] = 0.0;
					}
				}
			}

			rLocalMatrix = temp_matrix;
		}
	}

	/// RHS only version of ApplySlipPenaltyCondition
	void ApplySlipPenaltyCondition(TLocalVectorType& rLocalVector,
			ConditionPointerType pCondition) const
	{
		const unsigned int local_size = rLocalVector.size();
		auto r_geometry = pCondition->GetGeometry();
		if (local_size > 0)
		{
			const unsigned int block_size = this->GetBlockSize();
			for(unsigned int itNode = 0; itNode < r_geometry.PointsNumber(); ++itNode)
			{
				if( this->IsSlip(pCondition) )
				{
					// We fix the first momentum dof (normal component) for each rotated block
					unsigned int j = itNode * block_size;

					// Remove all other value than the normal component
					for(unsigned int i = j; i < (j + block_size); ++i)
					{
						if (i!=j) rLocalVector[i] = 0.0;
					}
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
	std::string Info() const
	{
		std::stringstream buffer;
		buffer << "MPMNonconformingBoundaryRotationUtility";
		return buffer.str();
	}

	/// Print information about this object.
	void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "MPMNonconformingBoundaryRotationUtility";
	}

	/// Print object's data.
	void PrintData(std::ostream& rOStream) const {}

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

	// To be used when there is only displacement in the DOF (no additional variable block)
	template<unsigned int TDim>
	void RotateAndRecoverAuxPure(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			ConditionPointerType pCondition,
			const RotationConfiguration ThisConfiguration = RotateSystem) const
	{
		const unsigned int local_size = rLocalVector.size();
		const unsigned int num_blocks = local_size / GetBlockSize();

		BoundedMatrix<double,TDim,TDim> rRot;
		LocalRotationOperatorPure<TDim>(rRot, pCondition);

		if (ThisConfiguration == RecoverSystem)
			rRot = trans(rRot);

		if(this->IsSlip(pCondition))
		{
			BoundedMatrix<double,TDim,TDim> mat_block, tmp;
			array_1d<double,TDim> aux, aux1;

			for(unsigned int i=0; i<num_blocks; i++)
			{
				for(unsigned int j=0; j<num_blocks; j++)
				{
					ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*GetBlockSize(), j*GetBlockSize());
					noalias(tmp) = prod(mat_block,trans(rRot));
					noalias(mat_block) = prod(rRot,tmp);
					WriteBlockMatrix<TDim>(mat_block, rLocalMatrix, i*GetBlockSize(), j*GetBlockSize());
				}

				for(unsigned int k=0; k<TDim; k++)
				aux[k] = rLocalVector[i*GetBlockSize()+k];

				noalias(aux1) = prod(rRot,aux);

				for(unsigned int k=0; k<TDim; k++)
				rLocalVector[i*GetBlockSize()+k] = aux1[k];
			}
		}
	}

	template<unsigned int TDim>
	void RotateAndRecoverAuxPureRHS(
			TLocalVectorType& rLocalVector,
			ConditionPointerType pCondition,
			const RotationConfiguration ThisConfiguration = RotateSystem) const
	{
		if (rLocalVector.size() > 0)
		{
			if(this->GetBlockSize() == this->GetDomainSize())
			{
				for(unsigned int j = 0; j < pCondition->GetGeometry().PointsNumber(); ++j)
				{
					if( this->IsSlip(pCondition) )
					{
						array_1d<double,TDim> aux,aux1;
						BoundedMatrix<double,TDim,TDim> rRot;
						LocalRotationOperatorPure<TDim>(rRot, pCondition);

						if (ThisConfiguration == RecoverSystem)
							rRot = trans(rRot);

						for(unsigned int k=0; k<TDim; k++)
						aux[k] = rLocalVector[j*GetBlockSize()+k];

						noalias(aux1) = prod(rRot,aux);

						for(unsigned int k=0; k<TDim; k++)
						rLocalVector[j*GetBlockSize()+k] = aux1[k];
					}
				}
			}
		}
	}

	template <unsigned int TDim>
	void LocalRotationOperatorPure(BoundedMatrix<double,TDim,TDim>& rRot,
			ConditionPointerType pCondition) const
	{
		array_1d<double,3> normal = pCondition->GetValue(MPC_NORMAL);
		ParticleMechanicsMathUtilities<double>::GetRotationMatrix<TDim>(rRot, normal);
	}

	///@}
	///@name Protected  Access
	///@{

	bool IsSlip(const ConditionPointerType pCurrentCondition) const
	{
		return pCurrentCondition->Is(mrFlag);
	}

	unsigned int GetBoundaryType(const ConditionPointerType pCurrentCondition) const
	{
		return pCurrentCondition->GetValue(mrFlagVariable);
	}

	unsigned int GetDomainSize() const
	{
		return mDomainSize;
	}

	unsigned int GetBlockSize() const
	{
		return mBlockSize;
	}

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

	const unsigned int mDomainSize;

	const unsigned int mBlockSize;

	const Kratos::Flags& mrFlag;

	const Variable<int>& mrFlagVariable;

	///@}
	///@name Member Variables
	///@{

	///@}
	///@name Private Operators
	///@{

	// Auxiliary functions

	template< unsigned int TBlockSize >
	void ReadBlockMatrix( BoundedMatrix<double,TBlockSize, TBlockSize>& rBlockMatrix,
						  const Matrix& rOriginMatrix, const unsigned int Ibegin, const unsigned int Jbegin) const
	{
		for(unsigned int i=0; i<TBlockSize; i++)
		{
			for(unsigned int j=0; j<TBlockSize; j++)
			{
				rBlockMatrix(i,j) = rOriginMatrix(Ibegin+i, Jbegin+j);
			}
		}
	}

	template< unsigned int TBlockSize >
	void WriteBlockMatrix( const BoundedMatrix<double,TBlockSize, TBlockSize>& rBlockMatrix,
						   Matrix& rDestinationMatrix, const unsigned int Ibegin, const unsigned int Jbegin) const
	{
		for(unsigned int i=0; i<TBlockSize; i++)
		{
			for(unsigned int j=0; j<TBlockSize; j++)
			{
				rDestinationMatrix(Ibegin+i, Jbegin+j) = rBlockMatrix(i,j);
			}
		}
	}

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
		MPMNonconformingBoundaryRotationUtility<TLocalMatrixType, TLocalVectorType>& rThis) {
	return rIStream;
}

/// output stream function
template<class TLocalMatrixType, class TLocalVectorType>
inline std::ostream& operator <<(std::ostream& rOStream,
		const MPMNonconformingBoundaryRotationUtility<TLocalMatrixType, TLocalVectorType>& rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

///@}

///@} addtogroup block

}

#endif // KRATOS_MPM_NONCONFORMING_BOUNDARY_ROTATION_UTILITY
