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

#include "includes/element.h"

namespace Kratos::Testing {

//class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoCustomElement : public Element {
class GeoCustomElement : public Element {
public:

	//KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoCustomElement);

	//using Element::Element;
	/**
    * Constructor using Geometry
    */
	GeoCustomElement()
		: Element()
	{
	};
 /// Constructor using Geometry
	GeoCustomElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
	};

	GeoCustomElement(IndexType NewId, const NodesArrayType& ThisNodes)
		: Element(NewId, GeometryType::Pointer(new GeometryType(ThisNodes)))
	{
	}

	///Copy constructor
	GeoCustomElement(GeoCustomElement const& rOther);

	// Destructor
	~GeoCustomElement() override
	{};

	///@}
	///@name Operators
	///@{

	/// Assignment operator.
	GeoCustomElement& operator=(GeoCustomElement const& rOther);


	Pointer GeoCustomElement::Create(
		IndexType NewId,
		NodesArrayType const& ThisNodes,
		PropertiesType::Pointer pProperties
	) const override
	{
		return Kratos::make_intrusive<GeoCustomElement>(NewId, GetGeometry().Create(ThisNodes));
	}


	Pointer Create(
		IndexType NewId,
		GeometryType::Pointer pGeom,
		PropertiesType::Pointer pProperties
	) const override
	{
		return Kratos::make_intrusive<GeoCustomElement>(NewId, pGeom);
	}


	void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateMassMatrix(MatrixType& rMassMatrix, 
		const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateDampingMatrix(MatrixType& rDampingMatrix, 
		const ProcessInfo& rCurrentProcessInfo) override;

	void EquationIdVector(EquationIdVectorType& rResult, 
		const ProcessInfo& rCurrentProcessInfo) const override;

	void CalculateRightHandSide(VectorType& rRightHandSideVector,
		const ProcessInfo& rCurrentProcessInfo) override;

	void GetDofList(Element::DofsVectorType& rDofList,
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
