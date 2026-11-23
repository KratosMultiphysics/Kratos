---
type: C++ Function
title: LocalRotationOperatorPure
resource: utilities/coordinate_transformation_utilities.h#L57-L192
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/coordinate_transformation_utilities/Normalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`class CoordinateTransformationUtils { public: ///@name Type Definitions ///@{ /// Pointer definition of CoordinateTransformationUtils KRATOS_CLASS_POINTER_DEFINITION(CoordinateTransformationUtils); using GeometryType = Geometry<Node>; ///@} ///@name Life Cycle ///@{ /// Constructor. /** @param DomainSize Number of space dimensions (2 or 3) * @param NumRowsPerNode Number of matrix or vector rows associated to each node. Velocity DOFs are assumed to be the first mDomainSize rows in each block of rows. * @param rSelectionFlag All nodes where the flag given by this argument is set to true will be transformed to a rotated coordinate system. */ CoordinateTransformationUtils( const unsigned int DomainSize, const unsigned int NumRowsPerNode, const Kratos::Flags& rSelectionFlag = SLIP) : mDomainSize(DomainSize) , mBlockSize(NumRowsPerNode) , mRotationDofPosition(0) , mrFlag(rSelectionFlag) { } CoordinateTransformationUtils( const unsigned int DomainSize, const unsigned int NumRowsPerNode, const unsigned int RotationDofPosition, const Kratos::Flags& rSelectionFlag = SLIP) : mDomainSize(DomainSize) , mBlockSize(NumRowsPerNode) , mRotationDofPosition(RotationDofPosition) , mrFlag(rSelectionFlag) { } /// Destructor. virtual ~CoordinateTransformationUtils() {} ///@} ///@name Operators ///@{ ///@} ///@name Operations ///@{ /** * @brief Calculates rotation operator for given point * * This metod calculates rotation matrix for a given point. Nodal normal variable should be * assigned properly since rotation is calculated based on it. * * @param rRotationMatrix Output rotation matrix * @param rThisPoint Current node */ virtual void CalculateRotationOperatorPure( TLocalMatrixType& rRotationMatrix, const GeometryType::PointType& rThisPoint) const { KRATOS_TRY if (mDomainSize == 2) { BoundedMatrix<double, 2, 2> local_matrix; this->LocalRotationOperatorPure(local_matrix, rThisPoint); if (rRotationMatrix.size1() != 2 || rRotationMatrix.size2() != 2) { rRotationMatrix.resize(2, 2, false); } noalias(rRotationMatrix) = local_matrix; } else if (mDomainSize == 3) { BoundedMatrix<double, 3, 3> local_matrix; this->LocalRotationOperatorPure(local_matrix, rThisPoint); if (rRotationMatrix.size1() != 3 || rRotationMatrix.size2() != 3) { rRotationMatrix.resize(3, 3, false); } noalias(rRotationMatrix) = local_matrix; } else { KRATOS_ERROR << "Unsupported domain size [ mDomainSize = " << mDomainSize << " ].\n"; } KRATOS_CATCH(""); } void LocalRotationOperatorPure( BoundedMatrix<double, 3, 3>& rRot, const GeometryType::PointType& rThisPoint) const`

# Calls

- [Normalize](../../../functions/utilities/coordinate_transformation_utilities/Normalize.md)