/*
==============================================================================
KratosStructuralApplication
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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: vladislav $
//   Date:                $Date: 2009-01-14 09:30:38 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_STRUCTURAL_APP_TIMOSHENKO_BEAM_ELEMENT_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APP_TIMOSHENKO_BEAM_ELEMENT_H_INCLUDED


// System includes


// External includes
#include "boost/smart_ptr.hpp"
#include "boost/math/constants/constants.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

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

/// Short class definition.
/** Detail class definition.
 */

class TimoshenkoBeamElement : public Element
{

public:
    ///@name Type Definitions
    ///@{
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    KRATOS_CLASS_POINTER_DEFINITION( TimoshenkoBeamElement );

    ///@}
    ///@name Life Cycle
    ///@{

    // A private default constructor necessary for serialization
    TimoshenkoBeamElement()
    {
    }

    /// Default constructor.
    TimoshenkoBeamElement( IndexType NewId, GeometryType::Pointer pGeometry );
    TimoshenkoBeamElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~TimoshenkoBeamElement();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const;

    virtual Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                           PropertiesType::Pointer pProperties) const;

    void Initialize();
    
    virtual void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    virtual void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);
    
    virtual void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
    
    virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
    
    void CalculateExternalLoadVector(Matrix& Rotation,  Vector& LocalBody, Vector& GlobalBody);
    
    void CalculateTransformationMatrix(Matrix& Rotation);
    
    void CalculateLocalMatrix(Matrix& LocalMatrix);
                           
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);    
    
    virtual int Check(const ProcessInfo& rCurrentProcessInfo);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
          rOStream << "TimoshenkoBeamElement #" << Id();
      }

    /// Print object's data.
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

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  Element );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer,  Element );
    }

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

//    static const double PI = boost::math::constants::pi<double>();

    double mArea;
    double mArea_y;
    double mArea_z;
    double mInertia_x;
    double mInertia_y;
    double mInertia_z;
    double mLength;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateSectionProperties();

    void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX, const Matrix& Ncontainer );

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //TimoshenkoBeamElement& operator=(const TimoshenkoBeamElement& rOther);

    /// Copy constructor.
    //TimoshenkoBeamElement(const TimoshenkoBeamElement& rOther);


    ///@}

}; // Class TimoshenkoBeamElement

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_APP_TIMOSHENKO_BEAM_ELEMENT_H_INCLUDED defined 
