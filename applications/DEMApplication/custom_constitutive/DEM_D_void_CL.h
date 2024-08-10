/////////////////////////////////////////////////
// Author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu, chengshun.shang1996@gmail.com
// Date: June 2023
/////////////////////////////////////////////////

#if !defined(DEM_D_VOID_CL_H_INCLUDED)
#define DEM_D_VOID_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_void : public DEMDiscontinuumConstitutiveLaw {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_void);

        DEM_D_void() {}

        ~DEM_D_void() {}

        std::string GetTypeOfLaw() override;

        void Check(Properties::Pointer pProp) const override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        std::unique_ptr<DEMDiscontinuumConstitutiveLaw> CloneUnique() override;

        void CalculateForces(const ProcessInfo& r_process_info,
                            const double OldLocalElasticContactForce[3],
                            double LocalElasticContactForce[3],
                            double LocalDeltDisp[3],
                            double LocalRelVel[3],
                            double indentation,
                            double previous_indentation,
                            double ViscoDampingLocalContactForce[3],
                            double& cohesive_force,
                            SphericParticle* element1,
                            SphericParticle* element2,
                            bool& sliding, double LocalCoordSystem[3][3]) override;

        void CalculateForcesWithFEM(const ProcessInfo& r_process_info,
                                    const double OldLocalElasticContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double LocalDeltDisp[3],
                                    double LocalRelVel[3],
                                    double indentation,
                                    double previous_indentation,
                                    double ViscoDampingLocalContactForce[3],
                                    double& cohesive_force,
                                    SphericParticle* const element,
                                    Condition* const wall,
                                    bool& sliding) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }

    }; //class DEM_D_void

} /* namespace Kratos.*/
#endif /* DEM_D_VOID_CL_H_INCLUDED  defined */
