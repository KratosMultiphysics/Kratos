//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
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
    
    void CalculateMeanContactArea(const bool has_mpi);

    ///@}
    ///@name Access
    array_1d<double,3> mLocalContactForce;
    double mContactSigma;
    double mContactTau;
    double mContactFailure;
    double mFailureCriterionState;
    double mUnidimendionalDamage;
    double mMeanContactArea;
    double mLocalContactAreaLow;
    double mLocalContactAreaHigh;


protected:


private:
  
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
