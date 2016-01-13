//Authors: M.A. Celigueta and G. Casas (CIMNE)
//   Date: January 2016

#if !defined(DEM_D_BENTONITE_COLLOID_CL_H_INCLUDED)
#define  DEM_D_BENTONITE_COLLOID_CL_H_INCLUDED

#include <string>
#include <iostream>

#include "DEM_application.h"
#include "DEM_discontinuum_constitutive_law.h"
#include "../custom_elements/spheric_particle.h"

namespace Kratos {

    class DEM_Bentonite_Colloid : public DEMDiscontinuumConstitutiveLaw {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_Bentonite_Colloid);

        DEM_Bentonite_Colloid();

        void Initialize(const ProcessInfo& rCurrentProcessInfo);

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;
        
        std::string GetTypeOfLaw();

        ~DEM_Bentonite_Colloid() {
        }

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation);  

        void InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation, const double ini_delta = 0.0);

        void CalculateForces(ProcessInfo& rCurrentProcessInfo,
                             const double OldLocalContactForce[3],
                            double LocalElasticContactForce[3],
                            double LocalDeltDisp[3],
                            double LocalRelVel[3],            
                            double indentation,
                            double previous_indentation,
                            double ViscoDampingLocalContactForce[3],
                            double& cohesive_force,
                            SphericParticle* element1,
                            SphericParticle* element2,
                            bool& sliding);
        
        void CalculateForcesWithFEM(ProcessInfo& rCurrentProcessInfo,
                                    const double OldLocalContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double LocalDeltDisp[3],
                                    double LocalRelVel[3],            
                                    double indentation,
                                    double previous_indentation,
                                    double ViscoDampingLocalContactForce[3],
                                    double& cohesive_force,
                                    SphericParticle* const element,
                                    DEMWall* const wall,
                                    bool& sliding);


        double CalculateNormalForce(const double indentation);
        double CalculateVanDerWaalsForce(const double);
        double CalculateDiffuseDoubleLayerForce(const double);
        
        double CalculateCohesiveNormalForce(SphericParticle* const element1,
                                            SphericParticle* const element2,
                                            const double indentation);

        double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element,
                                            DEMWall* const wall,
                                            const double indentation);

        template <class NeighbourClassType>
        void CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                            const double OldLocalContactForce[3],
                                            double LocalElasticContactForce[3],
                                            double ViscoDampingLocalContactForce[3],
                                            const double LocalDeltDisp[3],            
                                            bool& sliding,
                                            SphericParticle* const element,
                                                   NeighbourClassType* const neighbour,
                                            double indentation,
                                            double previous_indentation);

        void CalculateViscoDampingForce(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        SphericParticle * const element1,
                                        SphericParticle* const element2);

        void CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        SphericParticle* const element,
                                        DEMWall* const wall);
        
        
    private:

        double mA_H;
        double mA_p;
        double mThickness;
        double mDDLCoefficient;
        double mDebyeLengthInv;
        
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} /* namespace Kratos.*/
#endif /* DEM_D_BENTONITE_COLLOID_CL_H_INCLUDED  defined */
