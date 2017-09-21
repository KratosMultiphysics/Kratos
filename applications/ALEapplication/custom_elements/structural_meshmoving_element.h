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
 *  Revision:            $Revision: 1.4 $
 * ***************************************************************************/

#if !defined(KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED)
#define KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED

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
 * solves the equations of solid mechanics using a modified linear elastic
 * consitutive law. The implementation is based on:
 * K. Stein, et al. Mesh moving techniques for fluid-structure interactions
 * with large displacements, ASME J. Appl. Mech. 70 (2003) 58-63.
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

// template<unsigned int TDim>
class StructuralMeshMovingElement : public Element
{
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

    StructuralMeshMovingElement(IndexType NewId, GeometryType::Pointer pGeometry);

    StructuralMeshMovingElement(IndexType NewId,
                                GeometryType::Pointer pGeometry,
                                PropertiesType::Pointer pProperties);

    virtual ~StructuralMeshMovingElement()
    {
    }

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
    BaseType::Pointer Create(IndexType NewId,
                             NodesArrayType const& rThisNodes,
                             PropertiesType::Pointer pProperties) const;

    BaseType::Pointer Create(IndexType NewId,
                             GeometryType::Pointer pGeom,
                             PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo);

    /**
    * Sets on rResult the ID's of the element degrees of freedom
    */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /**
    * Sets on rElementalDofList the degrees of freedom of the considered element
    * geometry
    */

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

    /// Initialize initial values of the element (must be called before
    /// calculation is done)
    void Initialize();

    MatrixType SetAndModifyConstitutiveLaw(const int& dimension, const double& rPointNumber);

    MatrixType CalculateBMatrix(const int& dimension, const double& rPointNumber);

    void CheckElementMatrixDimension(MatrixType& rLeftHandSideMatrix,
                                     VectorType& rRightHandSideVector);

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo);

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
    IntegrationMethod mThisIntegrationMethod;
    GeometryType::JacobiansType mJ0;
    GeometryType::JacobiansType mInvJ0;
    VectorType mDetJ0;
    double mTotalDomainInitialSize;
    SizeType mLocalSize;
    ///@}
    ///@name Member Variables
    ///@{
    ///@}

    StructuralMeshMovingElement()
    {
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

    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    ///@}

}; // Class StructuralMeshMovingElement

///@}

///@name Type Definitions
///@{
///@}
}
// namespace Kratos.

#endif // KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED
