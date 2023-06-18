#ifndef KRATOS_DEM_COMPOUND_CONSTITUTIVE_LAW_FOR_PBM_H
#define KRATOS_DEM_COMPOUND_CONSTITUTIVE_LAW_FOR_PBM_H

#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos{

template <class BondCL, class UnbondCL>

class KRATOS_API(DEM_APPLICATION) DEM_compound_constitutive_law_for_PBM : public BondCL {

public:

    KRATOS_CLASS_POINTER_DEFINITION(DEM_compound_constitutive_law_for_PBM);

    DEMContinuumConstitutiveLaw::Pointer Clone() const override {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_compound_constitutive_law_for_PBM<BondCL, UnbondCL>(*this));
        return p_clone;
    }

    std::unique_ptr<DEMContinuumConstitutiveLaw> CloneUnique() {
        return std::unique_ptr<DEMContinuumConstitutiveLaw>{new DEM_compound_constitutive_law_for_PBM<BondCL, UnbondCL>(*this)};
    }

    virtual ~DEM_compound_constitutive_law_for_PBM() {}

    //add
    double ComputeNormalUnbondedForce(double indentation) override {
        return mCCL.CalculateNormalForce(indentation);
    }

    void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override {
        return mCCL.InitializeContact(element1, element2, indentation);
    }

    void CalculateUnbondedViscoDampingForce(double LocalRelVel[3],
                                                double UnbondedViscoDampingLocalContactForce[3],
                                                SphericParticle* const element1,
                                                SphericParticle* const element2) override{
        return mCCL.CalculateViscoDampingForce(LocalRelVel, UnbondedViscoDampingLocalContactForce, element1, element2);                                            
    }

private:

    UnbondCL mCCL;

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BondCL);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BondCL);
    }
};

} //namespace Kratos

#endif //KRATOS_DEM_COMPOUND_CONSTITUTIVE_LAW_FOR_PBM_H