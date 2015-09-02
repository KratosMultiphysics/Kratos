/*
==============================================================================
KratosAdjointFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
(Released on march 05, 2007).

Copyright 2015
Mate Pentek, Michael Andre
mate.pentek@tum.de
michael.andre@tum.de
- Lehrstuhl fuer Statik, Technische Universitaet Muenchen, Arcisstrasse
21 80333 Munich, Germany

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
//
//   Project Name:        KratosAdjointFluidApplication $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:             March 2015 $
//   Revision:            $Revision:                0.0 $
//
//

#if !defined(KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "adjoint_fluid_application_variables.h"

namespace Kratos {

///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{


/**
 * @brief A discrete adjoint element for fluid shape sensitivity.
 *
 * The adjoint element computes the sensitivity of the drag acting on a structure
 * to the variation in the boundary's shape. The structure is marked by the
 * STRUCTURE flag and the boundary by the BOUNDARY flag. The element requires the
 * solution of the unperturbed fluid problem to be stored in nodal variables
 * PRIMAL_VELOCITY and PRIMAL_PRESSURE. ADJOINT_VELOCITY and ADJOINT_PRESSURE are
 * used to store the solution of the adjoint problem. The discrete adjoint is
 * based on the monolithic fluid element. Sensitivities are stored in the
 * variable SHAPE_SENSITIVITY and are calculated w.r.t one drag direction which
 * is given by the value of the ProcessInfo variable DRAG_FORCE_TYPE
 * (0: x-direction (default), 1: y-direction, 2: z-direction).
 *
 * @see VMS monolithic fluid element for the primal solution
 * @see AdjointFluidStrategy solution strategy for the adjoint fluid problem
 */
template<unsigned int TDim, unsigned int TNumNodes = TDim + 1>
class VMSAdjointElement : public Element {
public:

  ///@name Type Definitions
  ///@{

  /// Pointer definition of VMSAdjointElement
  KRATOS_CLASS_POINTER_DEFINITION(VMSAdjointElement);

  typedef Element::IndexType IndexType;

  typedef Element::SizeType SizeType;

  typedef Element::GeometryType GeometryType;

  typedef Element::PropertiesType PropertiesType;

  typedef Element::NodesArrayType NodesArrayType;

  typedef Element::VectorType VectorType;

  typedef Element::MatrixType MatrixType;

  typedef Element::DofsVectorType DofsVectorType;

  typedef Element::EquationIdVectorType EquationIdVectorType;

  typedef boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>
  ShapeFunctionDerivativesType;

  ///@}
  ///@name Life Cycle
  ///@{

  VMSAdjointElement(IndexType NewId = 0) :
      Element(NewId), mDragDirection(0)
  {}

  VMSAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry) :
      Element(NewId, pGeometry), mDragDirection(0)
  {}

  VMSAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry,
                    PropertiesType::Pointer pProperties) :
      Element(NewId, pGeometry, pProperties), mDragDirection(0)
  {}

  virtual ~VMSAdjointElement()
  {}

  ///@}
  ///@name Operations
  ///@{

  /**
   * @brief Creates a new element of this type.
   *
   * @return pointer to the newly created element
   */
  virtual Element::Pointer Create(IndexType NewId,
                                  NodesArrayType const& ThisNodes,
                                  PropertiesType::Pointer pProperties) const
  {
    return Element::Pointer(
        new VMSAdjointElement<TDim>(NewId, this->GetGeometry().Create(ThisNodes),
                                    pProperties) );
  }

  /**
   * @brief Checks for proper element geometry, nodal variables and dofs.
   *
   * @return 0 after successful completion.
   */
  virtual int Check(const ProcessInfo &/*rCurrentProcessInfo*/)
  {
    KRATOS_TRY;

    // Check the element id and geometry.
    ProcessInfo UnusedProcessInfo;
    int ReturnValue = Element::Check(UnusedProcessInfo);

    // Check if adjoint and fluid variables are defined.
    if (ADJOINT_VELOCITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument, "ADJOINT_VELOCITY Key is 0. Check if "
                   "the application was correctly registered.","");
    if (ADJOINT_PRESSURE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument, "ADJOINT_PRESSURE Key is 0. Check if "
                   "the application was correctly registered.","");
    if (PRIMAL_VELOCITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument, "PRIMAL_VELOCITY Key is 0. Check if the "
                   "application was correctly registered.","");
    if (PRIMAL_PRESSURE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument, "PRIMAL_PRESSURE Key is 0. Check if the "
                   "application was correctly registered.","");

    // Check if the nodes have adjoint and fluid variables and adjoint dofs.
    for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
    {
      if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_VELOCITY) == false)
        KRATOS_THROW_ERROR(std::invalid_argument,
                     "missing ADJOINT_VELOCITY variable on solution step data for node ",
                     this->GetGeometry()[iNode].Id());
      if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_PRESSURE) == false)
        KRATOS_THROW_ERROR(std::invalid_argument,
                     "missing ADJOINT_PRESSURE variable on solution step data for node ",
                     this->GetGeometry()[iNode].Id());
      if (this->GetGeometry()[iNode].SolutionStepsDataHas(PRIMAL_VELOCITY) == false)
        KRATOS_THROW_ERROR(std::invalid_argument,
                     "missing PRIMAL_VELOCITY variable on solution step data for node ",
                     this->GetGeometry()[iNode].Id());
      if (this->GetGeometry()[iNode].SolutionStepsDataHas(PRIMAL_PRESSURE) == false)
        KRATOS_THROW_ERROR(std::invalid_argument,
                     "missing PRIMAL_PRESSURE variable on solution step data for node ",
                     this->GetGeometry()[iNode].Id());
      if (this->GetGeometry()[iNode].HasDofFor(ADJOINT_VELOCITY_X) == false
          || this->GetGeometry()[iNode].HasDofFor(ADJOINT_VELOCITY_Y) == false
          || this->GetGeometry()[iNode].HasDofFor(ADJOINT_VELOCITY_Z) == false)
        KRATOS_THROW_ERROR(std::invalid_argument,
                     "missing ADJOINT_VELOCITY component degree of freedom on node ",
                     this->GetGeometry()[iNode].Id());
      if (this->GetGeometry()[iNode].HasDofFor(ADJOINT_PRESSURE) == false)
        KRATOS_THROW_ERROR(std::invalid_argument,
                     "missing ADJOINT_PRESSURE component degree of freedom on node ",
                     this->GetGeometry()[iNode].Id());
    }

    return ReturnValue;

    KRATOS_CATCH("");
  }

  /**
   * @brief Sets the direction for drag sensitivities.
   *
   * @param rCurrentProcessInfo process info containing the drag direction
   */
  virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

    mDragDirection = rCurrentProcessInfo[DRAG_FORCE_TYPE];

    if (mDragDirection >= TDim)
      KRATOS_THROW_ERROR(std::invalid_argument,
                   "invalid drag direction found in element ", this->Id());

    KRATOS_CATCH("");
  }

  /**
   * @brief Evaluates the elemental contribution to the adjoint fluid problem.
   *
   * The adjoint fluid problem is:
   *
   * -(\partial R / \partial W)^T \Lambda = (\partial Drag / \partial W)^T
   *
   * with R the fluid residual vector and W the vector of fluid variables
   * (i.e. PRIMAL_VELOCITY and PRIMAL_PRESSURE) on the nodes.
   *
   * @param rLeftHandSideMatrix system matrix of the adjoint fluid problem
   * @param rRightHandSideVector derivatives of drag w.r.t fluid variables
   */
  virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                    VectorType& rRightHandSideVector,
                                    ProcessInfo& /*rCurrentProcessInfo*/)
  {
    const SizeType BlockSize = TDim + 1;
    const SizeType LocalSize = BlockSize * TNumNodes;
    VectorType DragFlagVector = ZeroVector(LocalSize);

    // The partial derivative of the drag w.r.t a fluid variable (e.g.
    // PRIMAL_VELOCITY_X) on a node is calculated from the row of the transpose
    // of the Jacobian matrix (i.e. the adjoint system matrix) whose index
    // coincides with that fluid variable. Its value is computed by summing the
    // row's coefficients which coincide with the velocity component acting in
    // the drag direction and restricted to the structure's nodes.
    //
    // This sum is computed (for all elemental fluid variables) as
    // -prod(rLeftHandSideMatrix, DragFlagVector)
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
      if (this->GetGeometry()[iNode].Is(STRUCTURE))
        DragFlagVector(LocalIndex + mDragDirection) = 1.0;

      LocalIndex += BlockSize;
    }

    this->CalculateAdjointSystemMatrix(rLeftHandSideMatrix);

    // For solving the linear system in residual based form.
    // A*dx = b - A*x0, dx=x-x0
    VectorType AdjointValues = ZeroVector(LocalSize);
    this->GetAdjointValuesVector(AdjointValues);

    rRightHandSideVector = ZeroVector(LocalSize);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,
                                          DragFlagVector + AdjointValues);
  }

  virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                      ProcessInfo& /*rCurrentProcessInfo*/)
  {
    MatrixType TmpLHS;
    ProcessInfo UnusedProcessInfo;
    this->CalculateLocalSystem(TmpLHS,rRightHandSideVector,UnusedProcessInfo);
  }

  /**
   * @brief Calculates the shape sensitivity.
   *
   * The shape sensitivity is assembled and stored in the nodal variable
   * SHAPE_SENSITIVITY after the adjoint problem is solved.
   *
   * @see AdjointFluidStrategy
   */
  virtual void Calculate(const Variable<array_1d<double,3> >& rVariable,
                         array_1d<double,3>& /*rOutput*/,
                         const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

    if (rVariable == SHAPE_SENSITIVITY)
    {
      const SizeType BlockSize = TDim + 1;
      const SizeType FluidLocalSize = BlockSize * TNumNodes;
      const SizeType CoordLocalSize = TDim * TNumNodes;

      // If this element has no nodes on the boundary, this function can be
      // skipped.
      bool SkipThisElement = true;
      for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        if (this->GetGeometry()[iNode].Is(BOUNDARY))
        {
          SkipThisElement = false;
          break;
        }
      
      if (SkipThisElement)
        return;
      
      VectorType DragFlagVector = ZeroVector(FluidLocalSize);
      IndexType FluidLocalIndex = 0;
      for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
      {
        if (this->GetGeometry()[iNode].Is(STRUCTURE))
          DragFlagVector(FluidLocalIndex + mDragDirection) = 1.0;
        
        FluidLocalIndex += BlockSize;
      }
      
      MatrixType ShapeDerivatives(CoordLocalSize,FluidLocalSize);
      this->CalculateShapeDerivativesMatrix(ShapeDerivatives);

      // Local contribution to (\partial Drag / \partial X)^T
      VectorType LocalSensitivity = ZeroVector(CoordLocalSize);
      noalias(LocalSensitivity) += prod(ShapeDerivatives,DragFlagVector);

      // Orthogonal subscales are neglected !!!
      //
      //if (rCurrentProcessInfo[OSS_SWITCH] == 1)
      //  this->AddOrthogonalSubscaleShapeDerivatives(LocalSensitivity,
      //                                              DragFlagVector);
      
      // Local contribution from adjoint variables
      VectorType AdjointValues = ZeroVector(FluidLocalSize);
      this->GetAdjointValuesVector(AdjointValues);
      noalias(LocalSensitivity) += prod(ShapeDerivatives,AdjointValues);
      
      // Carefully write results to nodal variables.
      IndexType CoordLocalIndex = 0;
      for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
      {
        if (this->GetGeometry()[iNode].Is(BOUNDARY))
        {
          this->GetGeometry()[iNode].SetLock();
          array_1d<double,3>& rSensitivity =
              this->GetGeometry()[iNode].FastGetSolutionStepValue(SHAPE_SENSITIVITY);
          for (IndexType d = 0; d < TDim; ++d)
            rSensitivity[d] += LocalSensitivity(CoordLocalIndex++);
          this->GetGeometry()[iNode].UnSetLock();
        }
        else
        {
          // Skip this node block.
          CoordLocalIndex += TDim;
        }
      }
    }
    else
      KRATOS_THROW_ERROR(std::invalid_argument,"unsupported variable: ",rVariable);
    
    KRATOS_CATCH("");
  }

  virtual void GetDofList(DofsVectorType& rElementalDofList,
                          ProcessInfo& /*rCurrentProcessInfo*/);

  virtual void EquationIdVector(EquationIdVectorType& rResult,
                                ProcessInfo& /*rCurrentProcessInfo*/);

  ///@}
  ///@name Input and output
  ///@{

  virtual std::string Info() const
  {
    std::stringstream buffer;
    buffer << "VMSAdjointElement" << this->GetGeometry().WorkingSpaceDimension()
           << "D #" << this->Id();
    return buffer.str();
  }

  virtual void PrintInfo(std::ostream& rOStream) const
  {
    rOStream << "VMSAdjointElement"
             << this->GetGeometry().WorkingSpaceDimension() << "D #"
             << this->Id() << std::endl;
    rOStream << "Number of Nodes: " << this->GetGeometry().PointsNumber()
             << std::endl;
  }

  virtual void PrintData(std::ostream& rOStream) const
  {
    this->PrintInfo(rOStream);
    rOStream << "Geometry Data: " << std::endl;
    this->GetGeometry().PrintData(rOStream);
  }

  ///@}

protected:

  ///@name Protected Operations
  ///@{

  /// Returns the adjoint values stored in this element's nodes.
  void GetAdjointValuesVector(VectorType& rValues);

  /// Returns the velocity and pressure values stored in this element's nodes.
  void GetFluidValuesVector(VectorType& rValues);

  /**
   * @brief Returns the velocity at this integration point.
   *
   * @param rVelocity velocity vector
   * @param rN array of shape function values at this integration point
   */
  void CalculateAdvectiveVelocity(VectorType& rVelocity,
                                  const array_1d< double, TNumNodes >& rN);

  /**
   * @brief Returns the gradient matrix of the velocity.
   *
   * The row index corresponds to the velocity component and the column index to
   * the derivative.
   *
   * @param rGradVel velocity gradient matrix
   * @param rDN_DX shape functions' gradients
   */
  void CalculateVelocityGradient(MatrixType& rGradVel,
                                 const ShapeFunctionDerivativesType& rDN_DX);


  /**
   * @brief Returns the pressure gradient.
   *
   * @param rGradP pressure gradient
   * @param rDN_DX shape functions' gradients
   */
  void CalculatePressureGradient(VectorType& rGradP,
                                 const ShapeFunctionDerivativesType& rDN_DX)
  {
    rGradP = ZeroVector(TDim);

    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
      for (IndexType d = 0; d < TDim; ++d)
      {
        rGradP[d] += rDN_DX(iNode,d)
            * this->GetGeometry()[iNode].FastGetSolutionStepValue(PRIMAL_PRESSURE);
      }
    }
  }

  /**
   * @brief Returns the external body force at this integration point.
   *
   * @param rBodyForce body force vector
   * @param rN array of shape function values at this integration point
   */
  void CalculateBodyForce(VectorType& rBodyForce,
                          const array_1d< double, TNumNodes >& rN);

  /**
   * @brief Returns the element's size.
   *
   * @param Volume the volume (area in 2D) of the element
   */
  double CalculateElementSize(const double Volume);

  /**
   * @brief Returns derivatives of determinant of Jacobian w.r.t coordinates.
   *
   * The derivative of the determinant of the Jacobian w.r.t the jth coordinate
   * of the ith node is stored at the index (i * TDim + j).
   *
   * This function is only valid when the determinant of the Jacobian is constant
   * over the element.
   *
   * @see Triangle2D3
   * @see Tetrahedra3D4
   */
  void CalculateDeterminantOfJacobianDerivatives(VectorType& rDetJDerivatives);

  /**
   * @brief Returns the VMS stabilization parameters.
   *
   * @param rTauOne momentum stabilization parameter
   * @param rTauTwo divergence stabilization parameter
   * @param VelNorm Euclidean norm of the velocity
   * @param ElemSize size of this element
   * @param Density density of the fluid
   * @param Viscosity dynamic viscosity of the fluid
   */
  void CalculateStabilizationParameters(double& rTauOne, double& rTauTwo,
                                        const double VelNorm,
                                        const double ElemSize,
                                        const double Density,
                                        const double Viscosity)
  {
    double InvTau = Density * (2.0 * VelNorm / ElemSize)
        + 4.0 * Viscosity / (ElemSize * ElemSize);
    rTauOne = 1.0 / InvTau;
    rTauTwo = Viscosity + 0.5 * Density * ElemSize * VelNorm;
  }

  /**
   * @brief Returns stabilization parameters derived w.r.t a node's coordinate.
   *
   * @param rTauOneDeriv derivative of momentum stabilization parameter
   * @param rTauTwoDeriv derivative of divergence stabilization parameter
   * @param TauOne momentum stabilization parameter
   * @param TauTwo divergence stabilization parameter
   * @param VelNorm Euclidean norm of the velocity
   * @param ElemSize size of this element
   * @param Density density of the fluid
   * @param Viscosity dynamic viscosity of the fluid
   * @param DetJDeriv derivative of the determinant of the Jacobian
   */
  void CalculateStabilizationParametersDerivative(
      double& rTauOneDeriv,
      double& rTauTwoDeriv,
      const double TauOne,
      const double TauTwo,
      const double VelNorm,
      const double ElemSize,
      const double Density,
      const double Viscosity,
      const double DetJDeriv);

  /**
   * @brief Returns a scalar variable at this integration point.
   *
   * @param rResult the value of the scalar variable at this integration point
   * @param rVariable the variable to be evaluated
   * @param rShapeFunc array of shape function values at this integration point
   */
  void EvaluateInPoint(double& rResult, const Variable< double >& rVariable,
                       const array_1d< double, TNumNodes >& rShapeFunc)
  {
    rResult = rShapeFunc[0]
        * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
    for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
      rResult += rShapeFunc[iNode]
          * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
  }

  /**
   * @brief Adds viscous contributions to adjoint system matrix.
   *
   * @param rResult matrix to add viscous contributions to
   * @param rDN_DX shape functions' gradients
   * @param Weight integration weight including dynamic viscosity
   */
  void AddViscousTerm(MatrixType& rResult,
                      const ShapeFunctionDerivativesType& rDN_DX,
                      const double Weight);

  /**
   * @brief Adds derivative of viscous term w.r.t a node's coordinate.
   *
   * @param rResult matrix to add viscous contributions to
   * @param rDN_DX shape functions' gradients
   * @param rDN_DX_Deriv shape functions' gradients derived w.r.t the coordinate
   * @param Weight integration weight including dynamic viscosity
   * @param WeightDeriv integration weight derived w.r.t the coordinate
   *
   * @see AddViscousTerm
   */
  void AddViscousTermDerivative(
      MatrixType& rResult,
      const ShapeFunctionDerivativesType& rDN_DX,
      const ShapeFunctionDerivativesType& rDN_DX_Deriv,
      const double Weight,
      const double WeightDeriv);

  /**
   * @brief Returns the elemental contribution to the adjoint system matrix.
   *
   * The adjoint system matrix is computed as -(\partial R / \partial W)^T with
   * R the fluid residual vector and W the vector of fluid variables (i.e.
   * PRIMAL_VELOCITY and PRIMAL_PRESSURE) on the nodes.
   */
  void CalculateAdjointSystemMatrix(MatrixType& rAdjointMatrix)
   {
    const SizeType BlockSize = TDim + 1;
    const SizeType LocalSize = BlockSize * TNumNodes;

    rAdjointMatrix = ZeroMatrix(LocalSize,LocalSize);

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType DN_DX;
    array_1d<double, TNumNodes> N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

    // Density
    double Density;
    this->EvaluateInPoint(Density,DENSITY,N);

    // Dynamic viscosity
    double Viscosity;
    this->EvaluateInPoint(Viscosity,VISCOSITY,N);
    Viscosity *= Density;

    // u
    VectorType Velocity(TDim);
    this->CalculateAdvectiveVelocity(Velocity,N);

    // u * Grad(N)
    VectorType DensityVelGradN(TNumNodes);
    noalias(DensityVelGradN) = Density * prod(DN_DX,Velocity);

    // Grad(u)
    MatrixType DensityGradVel(TDim,TDim);
    this->CalculateVelocityGradient(DensityGradVel,DN_DX);

    // Div(u)
    double DivVel = 0.0;
    for (IndexType d = 0; d < TDim; ++d)
      DivVel += DensityGradVel(d,d);

    DensityGradVel *= Density;

    // Grad(p)
    VectorType GradP(TDim);
    this->CalculatePressureGradient(GradP,DN_DX);

    // ( Grad(u) * Grad(N) )^T
    MatrixType DN_DX_DensityGradVel(TNumNodes,TDim);
    noalias(DN_DX_DensityGradVel) = prod(DN_DX,DensityGradVel);

    // ( u * Grad(u) * Grad(N) )^T
    VectorType DN_DX_DensityGradVel_Vel(TNumNodes);
    noalias(DN_DX_DensityGradVel_Vel) = prod(DN_DX_DensityGradVel,Velocity);

    // u * Grad(u)
    VectorType DensityGradVel_Vel(TDim);
    noalias(DensityGradVel_Vel) = prod(DensityGradVel,Velocity);

    // Grad(N)^T * Grad(p)
    VectorType DN_DX_GradP(TNumNodes);
    noalias(DN_DX_GradP) = prod(DN_DX,GradP);

    // Grad(N)^T * BodyForce
    VectorType BodyForce(TDim);
    VectorType DN_DX_BodyForce(TNumNodes);
    this->CalculateBodyForce(BodyForce,N);
    BodyForce *= Density;
    noalias(DN_DX_BodyForce) = prod(DN_DX,BodyForce);

    // Stabilization parameters TauOne, TauTwo
    double VelNorm = norm_2(Velocity);
    double ElemSize = this->CalculateElementSize(Volume);
    double TauOne, TauTwo;
    this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
                                           Density,Viscosity);

    // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
    // depend on the definitions of TauOne and TauTwo and should be consistent
    // with the fluid element used to solve for PRIMAL_VELOCITY and
    // PRIMAL_PRESSURE.
    MatrixType TauOneDeriv(TNumNodes,TDim);
    MatrixType TauTwoDeriv(TNumNodes,TDim);

    if (VelNorm > 0.0)
    {
      double CoefOne =-2.0 * Density * TauOne * TauOne / (ElemSize * VelNorm);
      double CoefTwo = 0.5 * Density * ElemSize / VelNorm;

      for (IndexType i = 0; i < TNumNodes; ++i)
      {
        for (IndexType d = 0; d < TDim; ++d)
        {
          TauOneDeriv(i,d) = CoefOne * N[i] * Velocity[d];
          TauTwoDeriv(i,d) = CoefTwo * N[i] * Velocity[d];
        }
      }
    }

    // Here, -(\partial R / \partial W) is calculated. This is the discrete
    // derivative of the fluid residual w.r.t the fluid variables and therefore
    // includes many of the terms defined in the fluid element. Neglecting the
    // transient terms of the fluid element, this matrix is identical to the
    // Jacobian of the fluid residual used for Newton-Raphson iterations. The
    // matrix is transposed at the end to get the adjoint system matrix.

    IndexType FirstRow(0), FirstCol(0);
    // Loop over nodes
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
      for (IndexType j = 0; j < TNumNodes; ++j)
      {
        double diag = 0.0;

        // Convective term, v * (u * Grad(u))
        diag += N[i] * DensityVelGradN[j];

        // Stabilization, lsq convection
        // (u * Grad(v)) * TauOne * (u * Grad(u))
        diag += DensityVelGradN[i] * TauOne * DensityVelGradN[j];

        for (IndexType m = 0; m < TDim; ++m)
        {
          for (IndexType n = 0; n < TDim; ++n)
          {
            double valmn = 0.0;

            // Convective term, v * (u * Grad(u))
            valmn += N[i] * N[j] * DensityGradVel(m,n);

            // Stabilization, lsq convection
            // (u * Grad(v)) * TauOne * (u * Grad(u))
            valmn += DensityVelGradN[i] * TauOne * N[j] * DensityGradVel(m,n);
            valmn += DensityVelGradN[i] * TauOneDeriv(j,n)
                * DensityGradVel_Vel[m];
            valmn += Density * N[j] * DN_DX(i,n) * TauOne
                * DensityGradVel_Vel[m];

            // Stabilization, lsq divergence
            // Div(v) * TauTwo * Div(u)
            valmn += DN_DX(i,m) * TauTwo * DN_DX(j,n);
            valmn += DN_DX(i,m) * TauTwoDeriv(j,n) * DivVel;

            // Stabilization, convection-pressure
            // (u * Grad(v)) * TauOne * Grad(p)
            valmn += TauOneDeriv(j,n) * DensityVelGradN[i] * GradP[m];
            valmn += Density * TauOne * N[j] * DN_DX(i,n) * GradP[m];

            // Stabilization, convection-BodyForce
            // (u * Grad(v)) * TauOne * f
            valmn -= N[j] * DN_DX(i,n) * TauOne * Density * BodyForce[m];
            valmn -= DensityVelGradN[i] * TauOneDeriv(j,n) * BodyForce[m];

            rAdjointMatrix(FirstRow+m,FirstCol+n) += Volume * valmn;
          }

          rAdjointMatrix(FirstRow+m,FirstCol+m) += Volume * diag;

          double valmp = 0.0;
          double valpn = 0.0;

          // Pressure term
          // Div(v) * p
          valmp -= DN_DX(i,m) * N[j];

          // Stabilization, convection-pressure
          // (u * Grad(v)) * TauOne * Grad(p)
          valmp += TauOne * DensityVelGradN[i] * DN_DX(j,m);

          // Divergence term
          // q * Div(u)
          valpn += N[i] * DN_DX(j,m);

          // Stabilization, lsq pressure
          // TauOne * Grad(q) * Grad(p)
          valpn += DN_DX_GradP[i] * TauOneDeriv(j,m);

          // Stabilization, pressure-convection
          // Grad(q) * TauOne * (u * Grad(u))
          valpn += DN_DX(i,m) * TauOne * DensityVelGradN[j];
          valpn += DN_DX_DensityGradVel(i,m) * TauOne * N[j];
          valpn += DN_DX_DensityGradVel_Vel[i] * TauOneDeriv(j,m);

          // Stabilization, pressure-BodyForce
          // Grad(q) * TauOne * f
          valpn -= DN_DX_BodyForce[i] * TauOneDeriv(j,m);

          rAdjointMatrix(FirstRow+m,FirstCol+TDim) += Volume * valmp;
          rAdjointMatrix(FirstRow+TDim,FirstCol+m) += Volume * valpn;
        }

        // Stabilization, lsq pressure
        // TauOne * Grad(q) * Grad(p)
        double valpp = 0.0;
        for (IndexType d = 0; d < TDim; ++d)
        {
          valpp += DN_DX(i,d) * DN_DX(j,d);
        }
        valpp *= TauOne;

        rAdjointMatrix(FirstRow+TDim,FirstCol+TDim) += Volume * valpp;

        FirstCol += BlockSize;
      }  // Node block columns

      FirstRow += BlockSize;
      FirstCol = 0;
    }  // Node block rows

    // Viscous term
    this->AddViscousTerm(rAdjointMatrix,DN_DX,Viscosity * Volume);

    // Transpose to get adjoint system matrix.
    rAdjointMatrix = trans(rAdjointMatrix);
  }

  /**
   * @brief Returns the partial derivatives of residual w.r.t shape parameters.
   *
   * The shape derivatives matrix is computed as (\partial R / \partial X)^T
   * with R the fluid residual vector and X the vector of shape parameters
   * (i.e. the coordinates of the element's nodes).
   *
   * This function is only valid when the determinant of the Jacobian is constant
   * over the element.
   */
  void CalculateShapeDerivativesMatrix(MatrixType& rShapeDerivativesMatrix)
  {
    const SizeType BlockSize = TDim + 1;
    const SizeType FluidLocalSize = BlockSize * TNumNodes;
    const SizeType CoordLocalSize = TDim * TNumNodes;

    if (rShapeDerivativesMatrix.size1() != CoordLocalSize
        || rShapeDerivativesMatrix.size2() != FluidLocalSize)
      rShapeDerivativesMatrix.resize(CoordLocalSize,FluidLocalSize);

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType DN_DX;
    array_1d<double, TNumNodes> N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

    // Density
    double Density;
    this->EvaluateInPoint(Density,DENSITY,N);

    // Dynamic viscosity
    double Viscosity;
    this->EvaluateInPoint(Viscosity,VISCOSITY,N);
    Viscosity *= Density;

    // u
    VectorType Velocity(TDim);
    this->CalculateAdvectiveVelocity(Velocity,N);

    // u * Grad(N)
    VectorType DensityVelGradN(TNumNodes);
    noalias(DensityVelGradN) = Density * prod(DN_DX,Velocity);

    // Det(J)
    const double InvDetJ = 1.0 / this->GetGeometry().DeterminantOfJacobian(0);
    VectorType DetJDerivatives(CoordLocalSize);
    this->CalculateDeterminantOfJacobianDerivatives(DetJDerivatives);

    // Stabilization parameters TauOne, TauTwo
    double VelNorm = norm_2(Velocity);
    double ElemSize = this->CalculateElementSize(Volume);
    double TauOne, TauTwo;
    this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
                                           Density,Viscosity);

    // External body force
    VectorType BodyForce(TDim);
    this->CalculateBodyForce(BodyForce,N);
    BodyForce *= Density;

    // We compute the derivative of the residual w.r.t each coordinate of each
    // node and assign it to the corresponding row of the shape derivatives
    // matrix.
    for (IndexType iCoord = 0; iCoord < CoordLocalSize; ++iCoord)
    {
      // Det(J)'
      double DetJDeriv = DetJDerivatives[iCoord];

      // DN_DX'
      MatrixType DN_DX_Deriv(TNumNodes,TDim);
      for (IndexType i = 0; i < TNumNodes; ++i)
        for (IndexType d = 0; d < TDim; ++d)
          DN_DX_Deriv(i,d) = -DN_DX(iCoord / TDim,d) * DN_DX(i,iCoord % TDim);

      // Volume'
      double VolumeDeriv = Volume * InvDetJ * DetJDeriv;

      // u * Grad(N)'
      VectorType DensityVelGradNDeriv(TNumNodes);
      noalias(DensityVelGradNDeriv) = Density * prod(DN_DX_Deriv,Velocity);

      // TauOne', TauTwo'
      double TauOneDeriv, TauTwoDeriv;
      this->CalculateStabilizationParametersDerivative(
          TauOneDeriv,TauTwoDeriv,TauOne,TauTwo,VelNorm,ElemSize,Density,
          Viscosity,DetJDeriv);

      MatrixType LHS = ZeroMatrix(FluidLocalSize,FluidLocalSize);
      VectorType RHS = ZeroVector(FluidLocalSize);
      for (IndexType i = 0; i < TNumNodes; ++i)
      {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
          // Left-hand side matrix
          double diag = 0.0;
          double ddiag = 0.0;

          // Convective term, v * (u * Grad(u))
          diag += N[i] * DensityVelGradN[j];
          ddiag += N[i] * DensityVelGradNDeriv[j];

          // Stabilization, lsq convection
          // (u * Grad(v)) * TauOne * (u * Grad(u))
          diag += DensityVelGradN[i] * TauOne * DensityVelGradN[j];
          ddiag += DensityVelGradNDeriv[i] * TauOne * DensityVelGradN[j]
              + DensityVelGradN[i] * TauOneDeriv * DensityVelGradN[j]
              + DensityVelGradN[i] * TauOne * DensityVelGradNDeriv[j];

          for (IndexType m = 0; m < TDim; ++m)
          {
            for (IndexType n = 0; n < TDim; ++n)
            {
              // Stabilization, lsq divergence
              // Div(v) * TauTwo * Div(u)
              double valmn = DN_DX(i,m) * TauTwo * DN_DX(j,n);
              double dvalmn = DN_DX_Deriv(i,m) * TauTwo * DN_DX(j,n)
                  + DN_DX(i,m) * TauTwoDeriv * DN_DX(j,n)
                  + DN_DX(i,m) * TauTwo * DN_DX_Deriv(j,n);

              LHS(i*BlockSize+m,j*BlockSize+n) += VolumeDeriv*valmn
                  + Volume*dvalmn;
            }
            LHS(i*BlockSize+m,j*BlockSize+m) += VolumeDeriv * diag
                + Volume * ddiag;

            double valmp = 0.0;
            double dvalmp = 0.0;
            // Pressure term
            // Div(v) * p
            valmp -= DN_DX(i,m) * N[j];
            dvalmp -= DN_DX_Deriv(i,m) * N[j];

            // Stabilization, convection-pressure
            // (u * Grad(v)) * TauOne * Grad(p)
            valmp += TauOne * DensityVelGradN[i] * DN_DX(j,m);
            dvalmp += TauOneDeriv * DensityVelGradN[i] * DN_DX(j,m)
                + TauOne * DensityVelGradNDeriv[i] * DN_DX(j,m)
                + TauOne * DensityVelGradN[i] * DN_DX_Deriv(j,m);

            double valpn = 0.0;
            double dvalpn = 0.0;
            // Divergence term
            // q * Div(u)
            valpn += N[i] * DN_DX(j,m);
            dvalpn += N[i] * DN_DX_Deriv(j,m);

            // Stabilization, pressure-convection
            // Grad(q) * TauOne * (u * Grad(u))
            valpn += TauOne * DensityVelGradN[j] * DN_DX(i,m);
            dvalpn += TauOneDeriv * DensityVelGradN[j] * DN_DX(i,m)
                + TauOne * DensityVelGradNDeriv[j] * DN_DX(i,m)
                + TauOne * DensityVelGradN[j] * DN_DX_Deriv(i,m);

            LHS(i*BlockSize+m,j*BlockSize+TDim) += VolumeDeriv * valmp
                + Volume * dvalmp;
            LHS(i*BlockSize+TDim,j*BlockSize+m) += VolumeDeriv * valpn
                + Volume * dvalpn;
          }

          double valpp = 0.0;
          double dvalpp = 0.0;
          // Stabilization, lsq pressure
          // TauOne * Grad(q) * Grad(p)
          for (IndexType d = 0; d < TDim; ++d)
          {
            valpp += DN_DX(i,d) * DN_DX(j,d) * TauOne;
            dvalpp += DN_DX_Deriv(i,d) * DN_DX(j,d) * TauOne
                + DN_DX(i,d) * DN_DX_Deriv(j,d) * TauOne
                + DN_DX(i,d) * DN_DX(j,d) * TauOneDeriv;
          }

          LHS(i*BlockSize+TDim,j*BlockSize+TDim) += VolumeDeriv * valpp
              + Volume * dvalpp;
        } // Node block columns

        // Right-hand side vector
        double DN_DX_BodyForce = 0.0;
        double DN_DX_BodyForceDeriv = 0.0;
        for (IndexType d = 0; d < TDim; ++d)
        {
          DN_DX_BodyForce += DN_DX(i,d) * BodyForce[d];
          DN_DX_BodyForceDeriv += DN_DX_Deriv(i,d) * BodyForce[d];
        }

        for (IndexType m = 0; m < TDim; ++m)
        {
          double valm = 0.0;
          double dvalm = 0.0;

          // External body force
          valm += N[i] * BodyForce[m];

          // Stabilization, convection-BodyForce
          // (u * Grad(v)) * TauOne * f
          valm += TauOne * DensityVelGradN[i] * BodyForce[m];
          dvalm += TauOneDeriv * DensityVelGradN[i] * BodyForce[m]
              + TauOne * DensityVelGradNDeriv[i] * BodyForce[m];

          RHS[i*BlockSize+m] += VolumeDeriv * valm + Volume * dvalm;
        }

        double valp = TauOne * DN_DX_BodyForce;
        double dvalp = TauOneDeriv * DN_DX_BodyForce
            + TauOne * DN_DX_BodyForceDeriv;

        RHS[i*BlockSize+TDim] += VolumeDeriv * valp + Volume * dvalp;
      } // Node block rows

      this->AddViscousTermDerivative(LHS,DN_DX,DN_DX_Deriv,Viscosity * Volume,
                                     Viscosity * VolumeDeriv);

      // Assign the derivative of the residual w.r.t this coordinate to the
      // shape derivatives matrix.
      VectorType ResidualDerivative(FluidLocalSize);
      VectorType FluidValues(FluidLocalSize);
      this->GetFluidValuesVector(FluidValues);
      noalias(ResidualDerivative) = RHS - prod(LHS,FluidValues);
      for (IndexType k = 0; k < FluidLocalSize; ++k)
        rShapeDerivativesMatrix(iCoord,k) = ResidualDerivative[k];
    }
  }

  ///@}

private:

  ///@name Member Variables
  ///@{

  /**
   * direction of drag sensitivities
   */
  unsigned int mDragDirection;

  ///@}
  ///@name Serialization
  ///@{

  friend class Serializer;

  virtual void save(Serializer& rSerializer) const
  {
    KRATOS_TRY;

    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    rSerializer.save("DragDirection", mDragDirection);

    KRATOS_CATCH("");
  }

  virtual void load(Serializer& rSerializer)
  {
    KRATOS_TRY;

    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
    rSerializer.load("DragDirection", mDragDirection);

    KRATOS_CATCH("");
  }

  ///@}
  ///@name Unaccessible methods
  ///@{

  VMSAdjointElement& operator=(VMSAdjointElement const& rOther);

  VMSAdjointElement(VMSAdjointElement const& rOther);

  ///@}

};  // class VMSAdjointElement

///@} // Kratos classes

///@name Input and output
///@{

/// Defines an input stream operator that does nothing.
template<unsigned int TDim>
inline std::istream& operator >>(std::istream& rIStream,
                                 VMSAdjointElement<TDim>& rThis)
{
  return rIStream;
}

/// Defines an output stream operator that prints element info.
template<unsigned int TDim>
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const VMSAdjointElement<TDim>& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}

///@} // Adjoint Fluid Application group

}// namespace Kratos

#endif // KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED defined
