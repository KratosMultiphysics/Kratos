// Authors: Miguel Angel Celigueta, maceli@cimne.upc.edu
//          Guillermo Casas, gcasas@cimne.upc.edu

#if !defined(DEM_D_BENTONITE_COLLOID_CL_H_INCLUDED)
#define  DEM_D_BENTONITE_COLLOID_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {
    
    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Bentonite_Colloid : public DEMDiscontinuumConstitutiveLaw {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Bentonite_Colloid);

        DEM_D_Bentonite_Colloid();

        void Initialize(const ProcessInfo& r_process_info) override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;
        
        std::string GetTypeOfLaw() override;

        ~DEM_D_Bentonite_Colloid() {
        }

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;  

        void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0) override;

        void CalculateForces(const ProcessInfo& r_process_info,
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
                            bool& sliding, double LocalCoordSystem[3][3]) override;
        
        void CalculateForcesWithFEM(ProcessInfo& r_process_info,
                                    const double OldLocalContactForce[3],
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


        using  DEMDiscontinuumConstitutiveLaw::CalculateNormalForce; //To avoid Clang Warning
        double CalculateNormalForce(const double distance, const double cation_concentration);
        double CalculateVanDerWaalsForce(const double);
        double CalculateDiffuseDoubleLayerForce(const double distance, const double cation_concentration);
        
        double CalculateCohesiveNormalForce(SphericParticle* const element1,
                                            SphericParticle* const element2,
                                            const double indentation) override;

        double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element,
                                            Condition* const wall,
                                            const double indentation) override;

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
                                        Condition* const wall);
        
        double mA_H;
        double mA_p;
        double mD_p;
        double mThickness;
        double mDDLCoefficient;
        double mEquivRadius;
        
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
    };

} /* namespace Kratos.*/
#endif /* DEM_D_BENTONITE_COLLOID_CL_H_INCLUDED  defined */
