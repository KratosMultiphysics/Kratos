// ==============================================================================
/*
 KratosALEApllication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Pooyan Dadvand, Riccardo Rossi, Andreas Winterstein
                     pooyan@cimne.upc.edu
                     rrossi@cimne.upc.edu
                     a.winterstein@tum.de
- CIMNE (International Center for Numerical Methods in Engineering),
  Gran Capita' s/n, 08034 Barcelona, Spain
- Chair of Structural Analysis, Technical University of Munich
  Arcisstrasse 21 80333 Munich, Germany

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
*/
//==============================================================================

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: June 2016 $
 *  Revision:            $Revision: 1.5 $
 * ***************************************************************************/

#if !defined(KRATOS_LAPLACIAN_MESHMOVING_ELEMENT_INCLUDED)
#define KRATOS_LAPLACIAN_MESHMOVING_ELEMENT_INCLUDED

// System includes

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

/// This class implements a laplacian mesh-updating scheme
/**
 * This mesh-updating scheme solves the Laplace equation in order to update the
 * mesh.
 * to distribute the motion of the structure into the fluid flow domain.
 * It can be used for a) quadrilateral, b) triangular, c) tetrahedral and d)
 * hexaedral elements,
 * working best for b) and c),
*/

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

// template< unsigned int TDim >
class LaplacianMeshMovingElement : public Element
{
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
    LaplacianMeshMovingElement(IndexType NewId,
                               GeometryType::Pointer pGeometry,
                               PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~LaplacianMeshMovingElement()
    {
    }
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations

    BaseType::Pointer Create(IndexType NewId,
                             NodesArrayType const& rThisNodes,
                             PropertiesType::Pointer pProperties) const;

    BaseType::Pointer Create(IndexType NewId,
                             GeometryType::Pointer pGeom,
                             PropertiesType::Pointer pProperties) const;

    void Initialize();

    MatrixType CalculateDerivatives(const int& rdimension, const double& rPointNumber);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

    void CalculateDeltaPosition(VectorType& IntermediateDisplacements,
                                ProcessInfo& rCurrentProcessInfo);

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
    void GetDisplacementValues(VectorType& rValues, const int Step = 0);
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
    IntegrationMethod mThisIntegrationMethod;
    GeometryType::JacobiansType mJ0;
    GeometryType::JacobiansType mInvJ0;
    VectorType mDetJ0;
    double mTotalDomainInitialSize;
    ///@}

    ///@name Serialization
    ///@{
    friend class Serializer;
    LaplacianMeshMovingElement()
    {
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
