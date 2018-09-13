#ifndef KRATOS_DEM_COMPOUND_CONSTITUTIVE_LAW_H
#define KRATOS_DEM_COMPOUND_CONSTITUTIVE_LAW_H

#include "includes/define.h"
#include "includes/serializer.h"
//#include "DEM_D_JKR_cohesive_law.h"
//#include "DEM_D_Hertz_viscous_Coulomb_CL.h"

namespace Kratos {

template <class MainCL, class CohesionCL>

class KRATOS_API(DEM_APPLICATION) DEM_compound_constitutive_law : public MainCL {

public:

    // Pointer types for DEM_compound_constitutive_law
    KRATOS_CLASS_POINTER_DEFINITION(DEM_compound_constitutive_law);

    DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_compound_constitutive_law<MainCL, CohesionCL>(*this));
        return p_clone;
    }

    void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_linear_viscous_Coulomb to Properties " << pProp->Id() << std::endl; // Print this correctly!!!
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    /// Destructor

    virtual ~DEM_compound_constitutive_law() {}

    double CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override {

        return mCCL.CalculateCohesiveNormalForce(element1, element2, indentation);

    };
    
    double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation) override {
    
        return mCCL.CalculateCohesiveNormalForceWithFEM(element, wall, indentation);
    };

private:

    CohesionCL mCCL;

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MainCL);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MainCL);
    }

}; // Class DEM_compound_constitutive_law : public MainCL

} // Namespace Kratos

#endif // KRATOS_DEM_COMPOUND_CONSTITUTIVE_LAW_H
