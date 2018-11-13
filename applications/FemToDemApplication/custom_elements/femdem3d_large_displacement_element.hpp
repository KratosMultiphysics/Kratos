//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
//

#if !defined(KRATOS_FEMDEM3D_LARGE_DISAPLCEMENT_ELEMENT_H_INCLUDED)
#define KRATOS_FEMDEM3D_LARGE_DISAPLCEMENT_ELEMENT_H_INCLUDED

#include "custom_elements/femdem3d_element.hpp"

namespace Kratos
{
class FemDem3DLargeDisplacementElement : public FemDem3DElement 
{

  public:
	/// Default constructors
	FemDem3DLargeDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry);

	FemDem3DLargeDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

	///Copy constructor
	FemDem3DLargeDisplacementElement(FemDem3DLargeDisplacementElement const &rOther);

	/// Destructor.
	virtual ~FemDem3DLargeDisplacementElement();

	/// Assignment operator.
	FemDem3DLargeDisplacementElement &operator=(FemDem3DLargeDisplacementElement const &rOther);

	Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;

	Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

	FemDem3DLargeDisplacementElement()
	{
	}

    void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo);
    void FinalizeNonLinearIteration(ProcessInfo &CurrentProcessInfo);

	void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                              VectorType &rRightHandSideVector,
							  ProcessInfo &rCurrentProcessInfo);

    int GetStrainSize(){return 6;}
    double CalculateDerivativesOnReferenceConfiguration(Matrix& rJ0,
                                                        Matrix& rInvJ0,
                                                        Matrix& rDN_DX,
                                                        const IndexType PointNumber,
                                                        IntegrationMethod ThisIntegrationMethod);

    void CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX);

    void CalculateAndAddMaterialK(MatrixType& rLeftHandSideMatrix,
                                    const Matrix& B,
                                    const Matrix& D,
                                    const double IntegrationWeight);

    void CalculateGeometricK(MatrixType& rLeftHandSideMatrix,
                            const Matrix& DN_DX,
                            const Vector& StressVector,
                            const double IntegrationWeight);

    void CalculateGreenLagrangeStrainVector(Vector& rStrainVector, const Matrix& F);
    
    void CalculateStressVectorPredictor(Vector& rStressVector, 
                                        const Matrix& rConstitutiveMAtrix, 
                                        const Vector& rStrainVector);

    void CalculateAndAddInternalForcesVector(Vector& rRightHandSideVector, 
                                            const Matrix& rB, 
                                            const Vector& rStressVector, 
                                            const double IntegrationWeight);

}; // Class
} // namespace Kratos
#endif // KRATOS_FEMDEM3D_LARGE_DISAPLCEMENT_ELEMENT_H_INCLUDED  defined