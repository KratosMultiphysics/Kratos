//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_PARTICLE_CONTACT_ELEMENT_H_INCLUDED )
#define  KRATOS_PARTICLE_CONTACT_ELEMENT_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
//#include "includes/constitutive_law.h"

namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) ParticleContactElement: public Element
{

public:

    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of ParticleContactElement
    KRATOS_CLASS_POINTER_DEFINITION( ParticleContactElement );

    /// Default constructor.
    ParticleContactElement( IndexType NewId, GeometryType::Pointer pGeometry );
    ParticleContactElement( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~ParticleContactElement();

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const override;

    void Initialize() override;

    void InitializeSolutionStep(ProcessInfo& r_process_info ) override;

    void FinalizeSolutionStep(ProcessInfo& r_process_info ) override;

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rOutput, const ProcessInfo& r_process_info) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& r_process_info) override;

    void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info ) override;

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

protected:


private:


    std::string Info() const override;

    friend class Serializer;

    // A private default constructor necessary for serialization
    ParticleContactElement() {}

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
    }

}; // Class ParticleContactElement

}  // namespace Kratos.
#endif // KRATOS_PARTICLE_CONTACT_ELEMENT_H_INCLUDED  defined
