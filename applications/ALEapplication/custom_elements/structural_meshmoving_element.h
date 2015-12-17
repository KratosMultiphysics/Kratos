/*
 ==============================================================================
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 Version 1.0 (Released on march 05, 2007).

 Copyright 2007
 Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
 pooyan@cimne.upc.edu
 rrossi@cimne.upc.edu
 janosch.stascheit@rub.de
 nagel@sd.rub.de
 - CIMNE (International Center for Numerical Methods in Engineering),
 Gran Capita' s/n, 08034 Barcelona, Spain
 - Ruhr-University Bochum, Institute for Structural Mechanics, Germany


 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 ==============================================================================
 */

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: Andreas.Mini@tum.de $
 *  Date:                $Date: November 2015 $
 *  Revision:            $Revision: 1.4 $
 * ***************************************************************************/

#if !defined( KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED )
#define  KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

/// This class implements a structural similarity mesh-updating scheme
/**
 * This mesh-updating scheme treats the mesh as a solid and therefore
 * solves the equations of solid mechanics using linear kinematics
 * and a linear elastic consitutive law. The stiffness of the elements
 * depends on their size and can be controlled by the Jacobian Determinant
 * weightened by an exponent.
 */

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

//template<unsigned int TDim>
class StructuralMeshMovingElement : public Element {
 public:
  ///@name Type Definitions
  ///@{
  /// Pointer definition of StructuralMeshMovingElement
  KRATOS_CLASS_POINTER_DEFINITION(StructuralMeshMovingElement);

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

  /// Default constructor.
  StructuralMeshMovingElement(IndexType NewId, GeometryType::Pointer pGeometry);

  /// Default constructor.
  StructuralMeshMovingElement(IndexType NewId, GeometryType::Pointer pGeometry,
      PropertiesType::Pointer pProperties);

  /// Destructor.
  virtual ~StructuralMeshMovingElement()
  {}

  ///@}
  ///@name Operators
  ///@{
  ///@}

  ///@name Operations
  ///@{

  /**
   * creates a new total lagrangian updated element pointer
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  BaseType::Pointer Create(IndexType NewId,
      NodesArrayType const& rThisNodes,
      PropertiesType::Pointer pProperties) const
  {
    const GeometryType& rGeom = this->GetGeometry();
    return BaseType::Pointer(new StructuralMeshMovingElement(NewId, rGeom.Create(rThisNodes), pProperties));
  }

  ///Bulid up system matrices
  void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo);

  void EquationIdVector(EquationIdVectorType& rResult,
      ProcessInfo& rCurrentProcessInfo);

  void GetDofList(DofsVectorType& rElementalDofList,
      ProcessInfo& rCurrentProcessInfo);

  ///Initialize the element (must be called before calculation is done)
  void Initialize();

  MatrixType SetAndModifyConstitutiveLaw(const int &dimension);

  MatrixType CalculateBMatrix(const int &dimension);

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
  void GetDisplacementValues(VectorType& rValues,
      const int Step = 0);

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
  IntegrationMethod mThisIntegrationMethod;
  GeometryType::JacobiansType mJ0;
  GeometryType::JacobiansType mInvJ0;
  VectorType mDetJ0;
  ///@}
  ///@name Member Variables
  ///@{
  ///@}

  ///@name Serialization
  ///@{

  friend class Serializer;

  StructuralMeshMovingElement() {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

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
  ///@}

};  // Class StructuralMeshMovingElement

///@}

///@name Type Definitions
///@{
///@}

}
  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED

