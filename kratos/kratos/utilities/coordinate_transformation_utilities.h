#ifndef KRATOS_COORDINATE_TRANSFORMATION_UTILITIES_H
#define KRATOS_COORDINATE_TRANSFORMATION_UTILITIES_H

// system includes

// external includes
#include "boost/numeric/ublas/matrix_proxy.hpp"

// kratos includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/variable.h"
#include "geometries/geometry.h"

namespace Kratos {

///@addtogroup KratosCore
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
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
class CoordinateTransformationUtils {
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of CoordinateTransformationUtils
	KRATOS_CLASS_POINTER_DEFINITION(CoordinateTransformationUtils);

	typedef Node<3> NodeType;

	typedef Geometry< Node<3> > GeometryType;

//     typedef boost::numeric::ublas::matrix_row<TLocalMatrixType>  LocalRowType;
//
//     typedef boost::numeric::ublas::matrix_range<TLocalMatrixType> MatrixBlockType;

	///@}
	///@name Life Cycle
	///@{

	/// Constructor.
	/** @param DomainSize Number of space dimensions (2 or 3)
	 * @param NumRowsPerNode Number of matrix or vector rows associated to each node. Velocity DOFs are assumed to be the first mDomainSize rows in each block of rows.
	 * @param rVariable Kratos variable used to flag nodes where local system contributions will be rotated. All nodes with rVariable != Zero will be rotated.
	 * @param Zero The zero value for the variable.
	 */
	CoordinateTransformationUtils(const unsigned int DomainSize,
			const unsigned int NumRowsPerNode,
			const Variable<TValueType>& rVariable,
			const TValueType Zero):
	mDomainSize(DomainSize),
	mBlockSize(NumRowsPerNode),
	mrFlagVariable(rVariable),
	mZero(Zero)
	{}

	/// Constructor.
	/** @param DomainSize Number of space dimensions (2 or 3)
	 * @param NumRowsPerNode Number of matrix or vector rows associated to each node. Velocity DOFs are assumed to be the first mDomainSize rows in each block of rows.
	 * @param rVariable Kratos variable used to flag nodes where local system contributions will be rotated. All nodes with rVariable != Zero will be rotated.
	 */
	CoordinateTransformationUtils(const unsigned int DomainSize,
			const unsigned int NumRowsPerNode,
			const Variable<TValueType> &rVariable):
	mDomainSize(DomainSize),
	mBlockSize(NumRowsPerNode),
	mrFlagVariable(rVariable),
	mZero(TValueType())
	{}

	/// Destructor.
	virtual ~CoordinateTransformationUtils() {}

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
	void Rotate(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		if(mBlockSize != mDomainSize) //Monolithic case
		{
			if(mDomainSize == 2) RotateAux<2>(rLocalMatrix,rLocalVector,rGeometry);
			if(mDomainSize == 3) RotateAux<3>(rLocalMatrix,rLocalVector,rGeometry);
		}
		else //fractional step case
		{
			if(mDomainSize == 2) RotateAuxPure<2>(rLocalMatrix,rLocalVector,rGeometry);
			if(mDomainSize == 3) RotateAuxPure<3>(rLocalMatrix,rLocalVector,rGeometry);
		}

	}

	/// RHS only version of Rotate
	void Rotate(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		//const unsigned int LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

		unsigned int Index = 0;

		if (rLocalVector.size() > 0)
		{
			if(mBlockSize != mDomainSize) //Monolithic case
			{
				for(unsigned int j = 0; j < rGeometry.PointsNumber(); ++j)
				{
					const NodeType& rNode = rGeometry[j];
					// The const reference is needed to ensure that we don't initialize the value for mrFlagVariable when it is not defined (which produces a race condtion in openmp).
					if( rNode.GetValue(mrFlagVariable) != mZero )
					{
						if(mDomainSize == 3)
						{
							array_1d<double,4> aux,aux1;
							boost::numeric::ublas::bounded_matrix<double,4,4> rRot;
							LocalRotationOperator(rRot,rGeometry[j]);

							for(unsigned int k=0; k<4; k++)
							aux[k] = rLocalVector[j*mBlockSize+k];

							noalias(aux1) = prod(rRot,aux);

							for(unsigned int k=0; k<4; k++)
							rLocalVector[j*mBlockSize+k] = aux1[k];
						}
						else
						{
							array_1d<double,3> aux,aux1;
							boost::numeric::ublas::bounded_matrix<double,3,3> rRot;
							LocalRotationOperator(rRot,rGeometry[j]);

							for(unsigned int k=0; k<3; k++)
							{
								aux[k] = rLocalVector[j*mBlockSize+k];
							}

							noalias(aux1) = prod(rRot,aux);

							for(unsigned int k=0; k<3; k++)
							rLocalVector[j*mBlockSize+k] = aux1[k];
						}
					}
					Index += mBlockSize;
				}

			}
			else //fractional step case
			{
				for(unsigned int j = 0; j < rGeometry.PointsNumber(); ++j)
				{
					const NodeType& rNode = rGeometry[j];
					// The const reference is needed to ensure that we don't initialize the value for mrFlagVariable when it is not defined (which produces a race condtion in openmp).
					if( rNode.GetValue(mrFlagVariable) != mZero )
					{
						if(mDomainSize == 3)
						{
							array_1d<double,3> aux,aux1;
							boost::numeric::ublas::bounded_matrix<double,3,3> rRot;
							LocalRotationOperatorPure(rRot,rGeometry[j]);

							for(unsigned int k=0; k<3; k++)
							aux[k] = rLocalVector[j*mBlockSize+k];

							noalias(aux1) = prod(rRot,aux);

							for(unsigned int k=0; k<3; k++)
							rLocalVector[j*mBlockSize+k] = aux1[k];
						}
						else
						{
							array_1d<double,2> aux,aux1;
							boost::numeric::ublas::bounded_matrix<double,2,2> rRot;
							LocalRotationOperatorPure(rRot,rGeometry[j]);

							for(unsigned int k=0; k<2; k++)
								aux[k] = rLocalVector[j*mBlockSize+k];

							noalias(aux1) = prod(rRot,aux);

							for(unsigned int k=0; k<2; k++)
							rLocalVector[j*mBlockSize+k] = aux1[k];
						}
					}
					Index += mBlockSize;
				}

			}

		}

	}

	/// Apply slip boundary conditions to the rotated local contributions.
	/** This function takes the local system contributions rotated so each
	 node's velocities are expressed using a base oriented with its normal
	 and imposes that the normal velocity is equal to the mesh velocity in
	 the normal direction.
	 */
	void ApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		const unsigned int LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

		if (LocalSize > 0)
		{
			for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
			{
				const NodeType& rNode = rGeometry[itNode]; // The const reference is needed to ensure that we don't initialize the value for mrFlagVariable when it is not defined.
				if( rNode.GetValue(mrFlagVariable) != mZero )
				{
					// We fix the first dof (normal velocity) for each rotated block
					unsigned int j = itNode * mBlockSize;
					//const double k = rLocalMatrix(j,j)+rLocalMatrix(j,j+1)+rLocalMatrix(j,j+2);

					// If the mesh is moving, we must impose v_normal = vmesh_normal
					array_1d<double,3> VMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
					VMesh -= rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
					array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
					this->Normalize(rN);

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

					rLocalVector(j) = inner_prod(rN,VMesh);
					rLocalMatrix(j,j) = 1.0;
				}
			}
		}
	}

	/// RHS only version of ApplySlipCondition
	void ApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		if (rLocalVector.size() > 0)
		{
			for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
			{
				const NodeType& rNode = rGeometry[itNode]; // The const reference is needed to ensure that we don't initialize the value for mrFlagVariable when it is not defined.
				if( rNode.GetValue(mrFlagVariable) != mZero )
				{
					// We fix the first dof (normal velocity) for each rotated block
					unsigned int j = itNode * mBlockSize;

					// If the mesh is moving, we must impose v_normal = vmesh_normal
					array_1d<double,3> VMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
					VMesh -= rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
					array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
					this->Normalize(rN);

					rLocalVector[j] = inner_prod(rN,VMesh);
				}
			}
		}
	}

	/// Transform nodal velocities to the rotated coordinates (aligned with each node's normal)
	void RotateVelocities(ModelPart& rModelPart) const
	{
		TLocalVectorType Vel(mDomainSize);
		TLocalVectorType Tmp(mDomainSize);

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
#pragma omp parallel for firstprivate(Vel,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			const NodeType& rNode = *itNode; // The const reference is needed to ensure that we don't initialize the value for mrFlagVariable when it is not defined.
			if( rNode.GetValue(mrFlagVariable) != mZero )
			{
				//this->RotationOperator<TLocalMatrixType>(Rotation,);
				if(mDomainSize == 3)
				{
					boost::numeric::ublas::bounded_matrix<double,3,3> rRot;
					LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
					for(unsigned int i = 0; i < 3; i++) Vel[i] = rVelocity[i];
					noalias(Tmp) = boost::numeric::ublas::prod(rRot,Vel);
					for(unsigned int i = 0; i < 3; i++) rVelocity[i] = Tmp[i];
				}
				else
				{
					boost::numeric::ublas::bounded_matrix<double,2,2> rRot;
					LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
					for(unsigned int i = 0; i < 2; i++) Vel[i] = rVelocity[i];
					noalias(Tmp) = boost::numeric::ublas::prod(rRot,Vel);
					for(unsigned int i = 0; i < 2; i++) rVelocity[i] = Tmp[i];
				}
			}
		}
	}

	/// Transform nodal velocities from the rotated system to the original one
	void RecoverVelocities(ModelPart& rModelPart) const
	{
		TLocalVectorType Vel(mDomainSize);
		TLocalVectorType Tmp(mDomainSize);

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
#pragma omp parallel for firstprivate(Vel,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			const NodeType& rNode = *itNode; // The const reference is needed to ensure that we don't initialize the value for mrFlagVariable when it is not defined.
			if( rNode.GetValue(mrFlagVariable) != mZero )
			{
				if(mDomainSize == 3)
				{
					boost::numeric::ublas::bounded_matrix<double,3,3> rRot;
					LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
					for(unsigned int i = 0; i < 3; i++) Vel[i] = rVelocity[i];
					noalias(Tmp) = boost::numeric::ublas::prod(boost::numeric::ublas::trans(rRot),Vel);
					for(unsigned int i = 0; i < 3; i++) rVelocity[i] = Tmp[i];
				}
				else
				{
					boost::numeric::ublas::bounded_matrix<double,2,2> rRot;
					LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
					for(unsigned int i = 0; i < 2; i++) Vel[i] = rVelocity[i];
					noalias(Tmp) = boost::numeric::ublas::prod(boost::numeric::ublas::trans(rRot),Vel);
					for(unsigned int i = 0; i < 2; i++) rVelocity[i] = Tmp[i];
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
	virtual std::string Info() const
	{
		std::stringstream buffer;
		buffer << "CoordinateTransformationUtils";
		return buffer.str();
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "CoordinateTransformationUtils";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const {}

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

	/// Number of spatial dimensions
	const unsigned int mDomainSize;

	/// Number of matrix or vector rows associated to each node.
	/** @note Velocity Dofs are assumed to be the first mDomainSize rows.
	 */
	const unsigned int mBlockSize;

	const Variable<TValueType>& mrFlagVariable;

	const TValueType mZero;

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	/// Normalize a vector.
	/**
	 * @param rThis the vector
	 * @return Original norm of the input vector
	 */
	template< class TVectorType >
	double Normalize(TVectorType& rThis) const
	{
		double Norm = 0;
		for(typename TVectorType::iterator iComponent = rThis.begin(); iComponent < rThis.end(); ++iComponent)
		Norm += (*iComponent)*(*iComponent);
		Norm = sqrt(Norm);
		for(typename TVectorType::iterator iComponent = rThis.begin(); iComponent < rThis.end(); ++iComponent)
		*iComponent /= Norm;
		return Norm;
	}

//     /// Compute a rotation matrix to transform values from the cartesian base to one oriented with the node's normal
//     /**
//      * The normal is read from solution step data NORMAL. Use NormalCalculationUtils::CalculateOnSimplex to
//      * obtain and store the nodal normal from the normals of the model's conditons.
//      * @param rRot The rotation matrix (output)
//      * @param rThisPoint The point used to orient the new coordinate system.
//      * @see NormalCalculationUtils
//      */
//     template<class TMatrixType>
//     void RotationOperator(TMatrixType& rRot,
//                           GeometryType::PointType& rThisPoint) const
//     {
//         typedef boost::numeric::ublas::matrix_row<TMatrixType> ThisRowType;
//         // Get the normal evaluated at the node
//         const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);
//
//         if(mDomainSize == 3)
//         {
//             // Define the new coordinate system, where the first vector is aligned with the normal
//             ThisRowType rN(rRot,0);
//             for( unsigned int i = 0; i < 3; ++i)
//                 rN[i] = rNormal[i];
//             this->Normalize(rN);
//
//             // To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
//             ThisRowType rT1(rRot,1);
//             rT1(0) = 1.0;
//             rT1(1) = 0.0;
//             rT1(2) = 0.0;
//
//             double dot = this->Dot(rN,rT1);
//
//             // It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
//             // If this is the case, repeat the procedure using (0,1,0)
//             if ( fabs(dot) > 0.99 )
//             {
//                 rT1(0) = 0.0;
//                 rT1(1) = 1.0;
//                 rT1(2) = 0.0;
//
//                 dot = this->Dot(rN,rT1);
//             }
//
//             // calculate projection and normalize
//             rT1 -= dot * rN;
//             this->Normalize(rT1);
//
//             // The third base component is choosen as N x T1, which is normalized by construction
//             ThisRowType rT2(rRot,2);
//             rT2(0) = rN(1)*rT1(2) - rN(2)*rT1(1);
//             rT2(1) = rN(2)*rT1(0) - rN(0)*rT1(2);
//             rT2(2) = rN(0)*rT1(1) - rN(1)*rT1(0);
//         }
//         else //if(mDomainSize == 2)
//         {
//             /* The basis for the new coordinate system is (normal,tangent)
//                Tangent vector is chosen (-normal_y, normal_x) so that the resulting base
//                is right-handed.
//              */
//             ThisRowType rN(rRot,0);
//             ThisRowType rT(rRot,1);
//
//             rN[0] = rNormal[0];
//             rN[1] = rNormal[1];
//             this->Normalize(rN);
//             rT[0] = -rN[1];
//             rT[1] = rN[0];
//         }
//
//     }

	template< class TVectorType >
	double Dot(const TVectorType& rV1,
			const TVectorType& rV2) const
	{
		double dot = 0.0;
		for( typename TVectorType::const_iterator iV1 = rV1.begin(),iV2 = rV2.begin(); iV1 != rV1.end(); ++iV1, ++iV2)
		{
			dot += (*iV1) * (*iV2);
		}
		return dot;
	}

	/// Transform a local contribution from cartesian coordinates to rotated ones
//     void ApplyRotation(TLocalMatrixType& rMatrix,
//                        const TLocalMatrixType& rRotation) const
//     {
//         // compute B = R*A*transpose(R)
//         const unsigned int LocalSize = rMatrix.size1();
//         const unsigned int NumBlocks = LocalSize / mBlockSize;
//         //TLocalMatrixType Tmp = ZeroMatrix(LocalSize,LocalSize);
// /*
//         for (unsigned int iBlock = 0; iBlock < NumBlocks; iBlock++)
//         {
//             for (unsigned int jBlock = 0; jBlock < NumBlocks; jBlock++)
//             {
//                 for (unsigned int i = iBlock*mBlockSize; i < (iBlock+1)*mBlockSize; i++)
//                 {
//                     for(unsigned int j = jBlock*mBlockSize; j < (jBlock+1)*mBlockSize; j++)
//                     {
//                         double& tij = Tmp(i,j);
//                         for(unsigned int k = iBlock*mBlockSize; k < (iBlock+1)*mBlockSize; k++)
//                         {
//                             for(unsigned int l = jBlock*mBlockSize; l < (jBlock+1)*mBlockSize; l++)
//                             {
//                                 tij += rRotation(i,k)*rMatrix(k,l)*rRotation(j,l);
//                             }
//                         }
//                     }
//                 }
//             }
//         }*/
//
// 	Matrix Tmp = prod(rMatrix,trans(rRotation));
// 	noalias(rMatrix) = prod(rRotation,Tmp);
//
// //         noalias(rMatrix) = Tmp;
//     }

	//auxiliary functions
	template< unsigned int TBlockSize >
	void ReadBlockMatrix( boost::numeric::ublas::bounded_matrix<double,TBlockSize, TBlockSize>& block, const Matrix& origin, const unsigned int Ibegin, const unsigned int Jbegin) const
	{
		for(unsigned int i=0; i<TBlockSize; i++)
		{
			for(unsigned int j=0; j<TBlockSize; j++)
			{
				block(i,j) = origin(Ibegin+i, Jbegin+j);
			}
		}
	}

	template< unsigned int TBlockSize >
	void WriteBlockMatrix( const boost::numeric::ublas::bounded_matrix<double,TBlockSize, TBlockSize>& block, Matrix& destination, const unsigned int Ibegin, const unsigned int Jbegin) const
	{
		for(unsigned int i=0; i<TBlockSize; i++)
		{
			for(unsigned int j=0; j<TBlockSize; j++)
			{
				destination(Ibegin+i, Jbegin+j) = block(i,j);
			}
		}
	}

	void LocalRotationOperator(boost::numeric::ublas::bounded_matrix<double,4,4>& rRot,
			GeometryType::PointType& rThisPoint) const
	{
		noalias(rRot) = IdentityMatrix(4,4);

		// Get the normal evaluated at the node
		const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

		double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1] + rNormal[2]*rNormal[2];
		aux = sqrt(aux);
		rRot(0,0) = rNormal[0]/aux;
		rRot(0,1) = rNormal[1]/aux;
		rRot(0,2) = rNormal[2]/aux;
		// Define the new coordinate system, where the first vector is aligned with the normal

		// To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
		array_1d<double,3> rT1;
		rT1(0) = 1.0;
		rT1(1) = 0.0;
		rT1(2) = 0.0;
		double dot = rRot(0,0);//this->Dot(rN,rT1);

		// It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
		// If this is the case, repeat the procedure using (0,1,0)
		if ( fabs(dot) > 0.99 )
		{
			rT1(0) = 0.0;
			rT1(1) = 1.0;
			rT1(2) = 0.0;

			dot = rRot(0,1); //this->Dot(rN,rT1);
		}

		// calculate projection and normalize
		rT1[0] -= dot*rRot(0,0);
		rT1[1] -= dot*rRot(0,1);
		rT1[2] -= dot*rRot(0,2);
		this->Normalize(rT1);
		rRot(1,0) = rT1[0];
		rRot(1,1) = rT1[1];
		rRot(1,2) = rT1[2];

		// The third base component is choosen as N x T1, which is normalized by construction
		rRot(2,0) = rRot(0,1)*rT1[2] - rRot(0,2)*rT1[1];
		rRot(2,1) = rRot(0,2)*rT1[0] - rRot(0,0)*rT1[2];
		rRot(2,2) = rRot(0,0)*rT1[1] - rRot(0,1)*rT1[0];
	}

	void LocalRotationOperator(boost::numeric::ublas::bounded_matrix<double,3,3>& rRot,
			GeometryType::PointType& rThisPoint) const
	{
		noalias(rRot) = IdentityMatrix(3,3);

		// Get the normal evaluated at the node
		const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

		double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1];
		aux = sqrt(aux);

		rRot(0,0) = rNormal[0]/aux;
		rRot(0,1) = rNormal[1]/aux;
		rRot(1,0) = -rNormal[1]/aux;
		rRot(1,1) = rNormal[0]/aux;

	}

	void LocalRotationOperatorPure(boost::numeric::ublas::bounded_matrix<double,3,3>& rRot,
			GeometryType::PointType& rThisPoint) const
	{

		// Get the normal evaluated at the node
		const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

		double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1] + rNormal[2]*rNormal[2];
		aux = sqrt(aux);
		rRot(0,0) = rNormal[0]/aux;
		rRot(0,1) = rNormal[1]/aux;
		rRot(0,2) = rNormal[2]/aux;
		// Define the new coordinate system, where the first vector is aligned with the normal

		// To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
		array_1d<double,3> rT1;
		rT1(0) = 1.0;
		rT1(1) = 0.0;
		rT1(2) = 0.0;
		double dot = rRot(0,0);//this->Dot(rN,rT1);

		// It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
		// If this is the case, repeat the procedure using (0,1,0)
		if ( fabs(dot) > 0.99 )
		{
			rT1(0) = 0.0;
			rT1(1) = 1.0;
			rT1(2) = 0.0;

			dot = rRot(0,1); //this->Dot(rN,rT1);
		}

		// calculate projection and normalize
		rT1[0] -= dot*rRot(0,0);
		rT1[1] -= dot*rRot(0,1);
		rT1[2] -= dot*rRot(0,2);
		this->Normalize(rT1);
		rRot(1,0) = rT1[0];
		rRot(1,1) = rT1[1];
		rRot(1,2) = rT1[2];

		// The third base component is choosen as N x T1, which is normalized by construction
		rRot(2,0) = rRot(0,1)*rT1[2] - rRot(0,2)*rT1[1];
		rRot(2,1) = rRot(0,2)*rT1[0] - rRot(0,0)*rT1[2];
		rRot(2,2) = rRot(0,0)*rT1[1] - rRot(0,1)*rT1[0];
	}

	void LocalRotationOperatorPure(boost::numeric::ublas::bounded_matrix<double,2,2>& rRot,
			GeometryType::PointType& rThisPoint) const
	{
		// Get the normal evaluated at the node
		const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

		double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1];
		aux = sqrt(aux);

		rRot(0,0) = rNormal[0]/aux;
		rRot(0,1) = rNormal[1]/aux;
		rRot(1,0) = -rNormal[1]/aux;
		rRot(1,1) = rNormal[0]/aux;

	}

	template<unsigned int TDim>
	void RotateAux(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		const unsigned int LocalSize = rLocalVector.size();

		unsigned int Index = 0;
		int rotations_needed = 0;
		const unsigned int NumBlocks = LocalSize / mBlockSize;
		boost::numeric::ublas::vector<bool> NeedRotation( NumBlocks, false);

		std::vector< boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1> > rRot(NumBlocks);
		for(unsigned int j = 0; j < NumBlocks; ++j)
		{
			const NodeType& rNode = rGeometry[j];
			// The const reference is needed to ensure that we don't initialize the value for mrFlagVariable when it is not defined (which produces a race condtion in openmp).
			if( rNode.GetValue(mrFlagVariable) != mZero )
			{
				NeedRotation[j] = true;
				rotations_needed++;

				LocalRotationOperator(rRot[j],rGeometry[j]);
			}

			Index += mBlockSize;
		}

		if(rotations_needed > 0)
		{
			boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1> mat_block, tmp;
			array_1d<double,TDim+1> aux, aux1;

			for(unsigned int i=0; i<NumBlocks; i++)
			{
				if(NeedRotation[i] == true)
				{
					for(unsigned int j=0; j<NumBlocks; j++)
					{
						if(NeedRotation[j] == true)
						{
							ReadBlockMatrix<TDim+1>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(mat_block,trans(rRot[j]));
							noalias(mat_block) = prod(rRot[i],tmp);
							WriteBlockMatrix<TDim+1>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
						else
						{
							ReadBlockMatrix<TDim+1>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(rRot[i],mat_block);
							WriteBlockMatrix<TDim+1>(tmp, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
					}

					for(unsigned int k=0; k<TDim+1; k++)
					aux[k] = rLocalVector[i*mBlockSize+k];

					noalias(aux1) = prod(rRot[i],aux);

					for(unsigned int k=0; k<TDim+1; k++)
					rLocalVector[i*mBlockSize+k] = aux1[k];

				}
				else
				{
					for(unsigned int j=0; j<NumBlocks; j++)
					{
						if(NeedRotation[j] == true)
						{
							ReadBlockMatrix<TDim+1>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(mat_block,trans(rRot[j]));
							WriteBlockMatrix<TDim+1>(tmp, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
					}
				}

			}
		}
	}

	//to be used when there is only velocity (no additional pressure or other var block)
	template<unsigned int TDim>
	void RotateAuxPure(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		const unsigned int LocalSize = rLocalVector.size();

		unsigned int Index = 0;
		int rotations_needed = 0;
		const unsigned int NumBlocks = LocalSize / mBlockSize;
		boost::numeric::ublas::vector<bool> NeedRotation( NumBlocks, false);

		std::vector< boost::numeric::ublas::bounded_matrix<double,TDim,TDim> > rRot(NumBlocks);
		for(unsigned int j = 0; j < NumBlocks; ++j)
		{
			const NodeType& rNode = rGeometry[j];
			// The const reference is needed to ensure that we don't initialize the value for mrFlagVariable when it is not defined (which produces a race condtion in openmp).
			if( rNode.GetValue(mrFlagVariable) != mZero )
			{
				NeedRotation[j] = true;
				rotations_needed++;

				LocalRotationOperatorPure(rRot[j],rGeometry[j]);
			}

			Index += mBlockSize;
		}

		if(rotations_needed > 0)
		{
			boost::numeric::ublas::bounded_matrix<double,TDim,TDim> mat_block, tmp;
			array_1d<double,TDim> aux, aux1;

			for(unsigned int i=0; i<NumBlocks; i++)
			{
				if(NeedRotation[i] == true)
				{
					for(unsigned int j=0; j<NumBlocks; j++)
					{
						if(NeedRotation[j] == true)
						{
							ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(mat_block,trans(rRot[j]));
							noalias(mat_block) = prod(rRot[i],tmp);
							WriteBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
						else
						{
							ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(rRot[i],mat_block);
							WriteBlockMatrix<TDim>(tmp, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
					}

					for(unsigned int k=0; k<TDim; k++)
					aux[k] = rLocalVector[i*mBlockSize+k];

					noalias(aux1) = prod(rRot[i],aux);

					for(unsigned int k=0; k<TDim; k++)
					rLocalVector[i*mBlockSize+k] = aux1[k];

				}
				else
				{
					for(unsigned int j=0; j<NumBlocks; j++)
					{
						if(NeedRotation[j] == true)
						{
							ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(mat_block,trans(rRot[j]));
							WriteBlockMatrix<TDim>(tmp, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
					}
				}

			}
		}
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
	CoordinateTransformationUtils& operator=(CoordinateTransformationUtils const& rOther) {}

	/// Copy constructor.
	CoordinateTransformationUtils(CoordinateTransformationUtils const& rOther) {}

	///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
inline std::istream& operator >>(std::istream& rIStream,
		CoordinateTransformationUtils<TLocalMatrixType, TLocalVectorType,
				TValueType>& rThis) {
	return rIStream;
}

/// output stream function
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
inline std::ostream& operator <<(std::ostream& rOStream,
		const CoordinateTransformationUtils<TLocalMatrixType, TLocalVectorType,
				TValueType>& rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

///@}

///@} addtogroup block

}

#endif // KRATOS_COORDINATE_TRANSFORMATION_UTILITIES_H
