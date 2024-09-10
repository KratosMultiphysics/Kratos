// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//
#pragma once

#include "includes/condition.h"
//#include "condition.h"

namespace Kratos::Testing {

class GeoCustomCondition : public Condition {
public:



	GeoCustomCondition()
		: Condition()
	{
	};
	/// Constructor using Geometry
	GeoCustomCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{
	};

	GeoCustomCondition(IndexType NewId, const NodesArrayType& ThisNodes)
		: Condition(NewId, GeometryType::Pointer(new GeometryType(ThisNodes)))
	{
	}

	///Copy constructor
	GeoCustomCondition(GeoCustomCondition const& rOther);

	// Destructor
	~GeoCustomCondition() override
	{};

	///@}
	///@name Operators
	///@{

	/// Assignment operator.
	GeoCustomCondition& operator=(GeoCustomCondition const& rOther);


	Pointer Create(
		IndexType NewId,
		NodesArrayType const& ThisNodes,
		PropertiesType::Pointer pProperties
	) const override
	{
		return Kratos::make_intrusive<GeoCustomCondition>(NewId, GetGeometry().Create(ThisNodes));
	}


	Pointer Create(
		IndexType NewId,
		GeometryType::Pointer pGeom,
		PropertiesType::Pointer pProperties
	) const override
	{
		return Kratos::make_intrusive<GeoCustomCondition>(NewId, pGeom);
	}


	void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateMassMatrix(MatrixType& rMassMatrix, 
		const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateDampingMatrix(MatrixType& rDampingMatrix, 
		const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(VectorType& rRightHandSideVector,
		const ProcessInfo& rCurrentProcessInfo) override;

	void EquationIdVector(EquationIdVectorType& rResult, 
		const ProcessInfo& rCurrentProcessInfo) const override;

	void GetDofList(Condition::DofsVectorType& rDofList,
		const ProcessInfo& rCurrentProcessInfo) const override;
private:

	MatrixType mMassMatrix;
	MatrixType mDampingMatrix;
	MatrixType mStiffnessMatrix;
	VectorType mRhs; 

	bool mStiffnessMatrixSet = false;
	bool mMassMatrixSet = false;
	bool mDampingMatrixSet = false;
	bool mRhsSet = false;

};

} // namespace Kratos::Testing
