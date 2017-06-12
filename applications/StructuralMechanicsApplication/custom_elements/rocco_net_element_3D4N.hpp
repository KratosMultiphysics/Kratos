// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Massimo Petracca
//


#if !defined(ROCCO_NET_ELEMENT_3D4N_H_INCLUDED )
#define  ROCCO_NET_ELEMENT_3D4N_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "utilities/quaternion.h"
#include "includes/define.h"
#include "includes/variables.h"

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

/** \brief ShellThickElement3D4N
 *
 * This element represents a 4-node bilinear Shell element
 * based on the Enhanced Assumed Strain Method (E.A.S.) for the membrane part
 * and on the Mixed Interpolation of Tensorial Components (M.I.T.C.)
 * for the trasverse shear part.
 * This element is formulated for small strains,
 * but can be used in Geometrically nonlinear problems
 * involving large displacements and rotations
 * using a Corotational Coordinate Transformation.
 * Material nonlinearity is handled by means of the cross section object.
 */
class RoccoNetElement3D4N : public Element
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(RoccoNetElement3D4N);

	typedef Element BaseType;
	typedef BaseType::GeometryType GeometryType;
	typedef BaseType::NodesArrayType NodesArrayType;
	typedef BaseType::PropertiesType PropertiesType;
	typedef BaseType::IndexType IndexType;
	typedef BaseType::SizeType SizeType;
	typedef BaseType::MatrixType MatrixType;
	typedef BaseType::VectorType VectorType;
	typedef BaseType::EquationIdVectorType EquationIdVectorType;
	typedef BaseType::DofsVectorType DofsVectorType;

    RoccoNetElement3D4N(IndexType NewId,
                          GeometryType::Pointer pGeometry);

    RoccoNetElement3D4N(IndexType NewId,
                          GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties);


    virtual ~RoccoNetElement3D4N();


    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
						 PropertiesType::Pointer pProperties) const;

    void Initialize();

    void EquationIdVector(EquationIdVectorType& rResult,
						ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList,
						ProcessInfo& CurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo);

    void GetValuesVector(Vector& rValues, int Step = 0) override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix,
							ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
							ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override;
	void CalculateLeftHandSide(
		MatrixType& rLeftHandSideMatrix,
		ProcessInfo& rCurrentProcessInfo) override;

	double CalculateActualCableLength();
	double CalculateActualDiagonalLength(const int node1, const int node2);

	double CalculateDiagonalLengthRefNodes(const int node1, const int node2);
	double CalculateCableLengthRefNodes();

	void CreateTransformationMatrix(Matrix& rRotationMatrix,
		const int node1, const int node2);

	void AddDiagonalStiffnessContribution(Matrix& rLeftHandSideMatrix);
	void AddCableStiffnessContribution(Matrix& rLeftHandSideMatrix);

	void AddDiagonalRHSForces(Vector& rInternalForces);
	void AddCableRHSForces(Vector& rInternalForces);

	double CalculateDiagonalStiffnessKb(const int node1, const int node2);

	void AddExplicitContribution(const VectorType& rRHSVector,
		const Variable<VectorType>& rRHSVariable,
		Variable<array_1d<double, 3> >& rDestinationVariable,
		const ProcessInfo& rCurrentProcessInfo) override;

private:


	double mDiameter, mCircumference, mThicknessWire, mNumberWindings;
	double mLengthDiagonalReference; //not equal diameter 
	double mKb, mKt, mdN1, mdN2;
	bool mNormalResistanceReached = false;
	bool mCableResistanceReached = false;
	Matrix mLHS;
	double mNodalMass_custom;

	RoccoNetElement3D4N() : Element() {}

    friend class Serializer;

};

}
#endif // ROCCO_NET_ELEMENT_3D4N_H_INCLUDED
