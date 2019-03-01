//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

#if !defined(KRATOS_LAPLACIAN_MESHMOVING_ELEMENT_INCLUDED)
#define KRATOS_LAPLACIAN_MESHMOVING_ELEMENT_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos {
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
/// This class implements a laplacian mesh-updating scheme
/**
 * This mesh-updating scheme solves the Laplace equation in order to update the
 * mesh.
 * to distribute the motion of the structure into the fluid flow domain.
 * It can be used for a) quadrilateral, b) triangular, c) tetrahedral and d)
 * hexaedral elements,
 * working best for b) and c),
*/
class LaplacianMeshMovingElement : public Element {
public:
  ///@name Type Definitions
  ///@{
  /// Counted pointer of LaplacianMeshMovingElement
  KRATOS_CLASS_POINTER_DEFINITION(LaplacianMeshMovingElement);

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

  typedef GeometryData::IntegrationMethod IntegrationMethod;
  ///@}

  ///@name Life Cycle
  /// Default constructor.
  LaplacianMeshMovingElement(IndexType NewId, GeometryType::Pointer pGeometry);

  /// Default constructor.
  LaplacianMeshMovingElement(IndexType NewId, GeometryType::Pointer pGeometry,
                             PropertiesType::Pointer pProperties);

  /// Destructor.
  virtual ~LaplacianMeshMovingElement() {}
  ///@{

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations

  BaseType::Pointer Create(IndexType NewId, NodesArrayType const &rThisNodes,
                           PropertiesType::Pointer pProperties) const override;

  BaseType::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override;

  void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                            VectorType &rRightHandSideVector,
                            ProcessInfo &rCurrentProcessInfo) override;

  void EquationIdVector(EquationIdVectorType &rResult,
                        ProcessInfo &rCurrentProcessInfo) override;

  void GetDofList(DofsVectorType &rElementalDofList,
                  ProcessInfo &rCurrentProcessInfo) override;

  void CalculateRightHandSide(VectorType &rRightHandSideVector,
                              ProcessInfo &rCurrentProcessInfo) override;

  int Check(const ProcessInfo& rCurrentProcessInfo) override;

  ///@{

  ///@}
  ///@name Access
  ///@{

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{
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

  ///@name Serialization
  ///@{
  friend class Serializer;
  LaplacianMeshMovingElement() {}
  ///@}

  ///@name Private Operators
  ///@{
  ///@}

  ///@name Private Operations
  ///@{

  void CalculateDeltaPosition(VectorType &IntermediateDisplacements,
                              ProcessInfo &rCurrentProcessInfo);

  void CheckElementMatrixDimension(MatrixType &rLeftHandSideMatrix,
                                   VectorType &rRightHandSideVector);

  MatrixType CalculateDerivatives(const int &rdimension,
                                  const double &rPointNumber);
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

}; // Class LaplacianMeshMovingElement

///@}

///@name Type Definitions
///@{
///@}

///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_LAPLACIAN_MESHMOVING_ELEMENT_INCLUDED  defined
