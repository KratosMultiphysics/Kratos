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
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2007-10-18 16:23:41 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_PARTICLE_CONTACT_ELEMENT_H_INCLUDED )
#define  KRATOS_PARTICLE_CONTACT_ELEMENT_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
//#include "includes/constitutive_law.h"

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

/// Total Lagrangian element for 2D and 3D geometries.
/**
 * Implements a total Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 2D and 3D
 */

class Particle_Contact_Element
    : public Element
{

public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of Particle_Contact_Element
    KRATOS_CLASS_POINTER_DEFINITION( Particle_Contact_Element );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Particle_Contact_Element( IndexType NewId, GeometryType::Pointer pGeometry );
    Particle_Contact_Element( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~Particle_Contact_Element();

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
   

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

    void Initialize();

    

    void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

    

    

    void InitializeSolutionStep( ProcessInfo& CurrentProcessInfo );

    void FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo );

    void GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo);
  
    void Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo );
    
    void PrepareForPrinting();

    ///@}
    ///@name Access
    array_1d<double,3> mLocalContactForce;
    double mContactSigma;
    double mContactTau;
    double mContactFailure;
    double mFailureCriterionState;
    double mUnidimendionalDamage;

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
    //      virtual void PrintInfo(std::ostream& rOStream) const;

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
    
    //double mshearstress;
    //double mnormalstress;
    
    //std::array_1d<double,3> mlocalcontactforcehigh;
    //std::array_1d<double,3> mlocalcontactforcelow;
    //std::array_1d<double,3> mlocalcontactforcemean;
    
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
    /*  static Matrix msB;
    static Matrix msF;
    static Matrix msD;
    static Matrix msC;
    static Vector msStrainVector;
    static Vector msStressVector;
    static Matrix msDN_DX;
    */
    ///@}
    ///@name Member Variables
    ///@{
    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;
    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    double mTotalDomainInitialSize;
    std::vector< Matrix > mInvJ0;
    Vector mDetJ0;
    ///@}
    ///@name Private Operators
    ///@{

    void CalculateAndAddKm(
        MatrixType& K,
        Matrix& B,
        Matrix& D,
        double weight );

    /**
     * Calculation of the Geometric Stiffness Matrix. Kg = dB * S
     */
    void CalculateAndAddKg(
        MatrixType& K,
        Matrix& DN_DX,
        Vector& StressVector,
        double weight
    );

    void CalculateBodyForces(
        Vector& BodyForce,
        const ProcessInfo& CurrentProcessInfo
    );

    void InitializeVariables();

    

    double CalculateIntegrationWeight
    ( double GaussPointWeight,
      double DetJ0 );

    void CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        Vector& BodyForce,
        VectorType& mResidualVector,
        double weight
    );

    void CalculateStrain( const Matrix& C,
                          Vector& StrainVector );

    void CalculateB( Matrix& B,
                     Matrix& DN_DX,
                     unsigned int StrainSize );

    void  Comprobate_State_Vector( Vector& Result );

    std::string Info() const;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    Particle_Contact_Element() {}

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
    }

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //Particle_Contact_Element& operator=(const Particle_Contact_Element& rOther);
    /// Copy constructor.
    //Particle_Contact_Element(const Particle_Contact_Element& rOther);
    ///@}

}; // Class Particle_Contact_Element

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
Particle_Contact_Element& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
const Particle_Contact_Element& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}*/
///@}

}  // namespace Kratos.
#endif // KRATOS_PARTICLE_CONTACT_ELEMENT_H_INCLUDED  defined 
