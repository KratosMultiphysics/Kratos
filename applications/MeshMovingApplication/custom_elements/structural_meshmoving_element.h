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

#if !defined(KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED)
#define KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED

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
/// This class implements a structural similarity mesh-updating scheme
/**
 * This mesh-updating scheme treats the mesh as a solid and therefore
 * solves the equations of solid mechanics using a modified linear elastic
 * consitutive law. The implementation is based on:
 * K. Stein, et al. Mesh moving techniques for fluid-structure interactions
 * with large displacements, ASME J. Appl. Mech. 70 (2003) 58-63.
 */

class StructuralMeshMovingElement : public Element {
public:
  ///@name Type Definitions
  ///@{
  /// Pointer definition of StructuralMeshMovingElement
  KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StructuralMeshMovingElement);

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
  ///@{

  StructuralMeshMovingElement(IndexType NewId, GeometryType::Pointer pGeometry);

  StructuralMeshMovingElement(IndexType NewId, GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties);

  virtual ~StructuralMeshMovingElement() {}

  ///@}
  ///@name Operators
  ///@{
  /// Assignment operator.
  ///@}

  ///@name Operations
  ///@{
  /**
  * Returns the currently selected integration method
  * @return current integration method selected
  */
  /**
   * creates a new total lagrangian updated element pointer
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  BaseType::Pointer Create(IndexType NewId, NodesArrayType const &rThisNodes,
                           PropertiesType::Pointer pProperties) const override;

  BaseType::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override;

  void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                            VectorType &rRightHandSideVector,
                            const ProcessInfo &rCurrentProcessInfo) override;

  /**
  * Sets on rResult the ID's of the element degrees of freedom
  */
  void EquationIdVector(EquationIdVectorType &rResult,
                        const ProcessInfo &rCurrentProcessInfo) const override;

  /**
  * Sets on rElementalDofList the degrees of freedom of the considered element
  * geometry
  */

  void GetDofList(DofsVectorType &rElementalDofList,
                  const ProcessInfo &rCurrentProcessInfo) const override;

  void CalculateRightHandSide(VectorType &rRightHandSideVector,
                              const ProcessInfo &rCurrentProcessInfo) override;

  void GetValuesVector(VectorType &rValues, int Step = 0) const override;

  int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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

  /**
   * Gets displacement values at nodes
   * @param rValues: reference to vector of nodal displacements
   */

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

  // A private default constructor necessary for serialization
  StructuralMeshMovingElement() {}

  MatrixType SetAndModifyConstitutiveLaw(const int Dimension,
                                         const double PointNumber) const;

  MatrixType CalculateBMatrix(const int Dimension, const double PointNumber) const;

  void CheckElementMatrixDimension(MatrixType &rLeftHandSideMatrix,
                                   VectorType &rRightHandSideVector) const;
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

  ///@name Serialization
  ///@{

  friend class Serializer;

  void save(Serializer& rSerializer) const override
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
  }

  void load(Serializer& rSerializer) override
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
  }

  ///@}

}; // Class StructuralMeshMovingElement

///@}

///@name Type Definitions
///@{
///@}
}
// namespace Kratos.

#endif // KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED
