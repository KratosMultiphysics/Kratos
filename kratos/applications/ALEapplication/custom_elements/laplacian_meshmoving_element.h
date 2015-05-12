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


#if !defined( KRATOS_LAPLACIAN_MESHMOVING_ELEMENT_INCLUDED )
#define  KRATOS_LAPLACIAN_MESHMOVING_ELEMENT_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

/// This class implements a laplacian mesh-updating scheme
/**
 * This mesh-updating scheme solves the Laplace equation in order to update the mesh.
 * It uses the L2 norm of the linear strain tensor
 * to distribute the motion of the structure into the fluid flow domain.
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
class LaplacianMeshMovingElement
        : public Element
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of LaplacianMeshMovingElement
    KRATOS_CLASS_POINTER_DEFINITION(LaplacianMeshMovingElement);
    ///@}

    ///@name Life Cycle
    /// Default constructor.
    LaplacianMeshMovingElement(IndexType NewId,
                               GeometryType::Pointer pGeometry)
    {}

    /// Default constructor.
    LaplacianMeshMovingElement(IndexType NewId,
                               GeometryType::Pointer pGeometry,
                               PropertiesType::Pointer pProperties):
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~LaplacianMeshMovingElement()
    {}
    ///@{

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations

    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    ///Bulid up system matrices
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    ///Get equation Id Vector
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    ///Get the degrees of freedom
    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
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

    /// Turn back information as a string.(Deactivated)
    //      virtual String Info() const;

    /// Print information about this object.(Deactivated)
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

#endif // KRATOS_LAPLACIAN_MESHMOVING_ELEMENT_INCLUDED  defined


