/*
==============================================================================
KratosALEApplication
A library based on:
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

/* *********************************************************
*
*   Last Modified by:    $Author: AMini $
*   Date:                $Date: Mai 2015 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/

#if !defined( KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_NONLINEAR_INCLUDED )
#define  KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_NONLINEAR_INCLUDED



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
* solves the equations of solid mechanics using non-linear kinematics
* and a linear elastic consitutive law. The stiffness of the elements
* depends on their size and can be controlled by by the Jacobian Determinant
* weightened by an exponent.
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


template< unsigned int TDim >
class StructuralMeshMovingElementNonlinear
        : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of StructuralMeshMovingElem2D
    KRATOS_CLASS_POINTER_DEFINITION(StructuralMeshMovingElementNonlinear);

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
    typedef GeometryType::JacobiansType JacobiansType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructuralMeshMovingElementNonlinear(IndexType NewId,
                                         GeometryType::Pointer pGeometry):
        Element(NewId, pGeometry)
    {}

    /// Default constructor.
    StructuralMeshMovingElementNonlinear(IndexType NewId,
                                         GeometryType::Pointer pGeometry,
                                         PropertiesType::Pointer pProperties):
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~StructuralMeshMovingElementNonlinear()
    {}

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
        return BaseType::Pointer( new StructuralMeshMovingElementNonlinear(NewId, rGeom.Create(rThisNodes), pProperties));
    }


    ///Initialize function to get values of initial configuration
    void Initialize();

    ///Bulid up system matrices
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo);

    ///Equation id Vector
    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo);

    ///Get the degrees of freedom
    void GetDofList(DofsVectorType& rElementalDofList,
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

    /// Turn back information as a string.(Deactivated)
    //      virtual String Info() const;

    /// Print information about this object. (Deactivated)
    //      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.(Deactivated)
    //      virtual void PrintData(std::ostream& rOStream) const;

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

    IntegrationMethod mThisIntegrationMethod;
    JacobiansType mInvJ0;
    double mArea0;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
    * Gets displacement values at nodes
    * @param rValues: reference to vector of nodal displacements
    */
    void GetDisplacementValues(VectorType& rValues, const int Step = 0);


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
    StructuralMeshMovingElementNonlinear() {}

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

}; // Class StructuralMeshMovingElem2DNonlin

///@}

///@name Type Definitions
///@{


/// input stream function (Deactivated)
/*  inline std::istream& operator >> (std::istream& rIStream,
                    LaplacianMeshMovingElement& rThis);
*/
/// output stream function (Deactivated)
/*  inline std::ostream& operator << (std::ostream& rOStream,
                    const LaplacianMeshMovingElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/

///@}

}  // namespace Kratos.

#endif //  KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_NONLINEAR_INCLUDED defined


