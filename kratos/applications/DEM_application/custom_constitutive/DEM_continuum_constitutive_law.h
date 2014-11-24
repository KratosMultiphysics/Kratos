
#if !defined(DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_utilities/properties_proxies.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "containers/flags.h"

#include "custom_utilities/GeometryFunctions.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"
namespace Kratos {

        /**
     * Base class of constitutive laws.
     */
    class /*__declspec( dllexport )*/ DEMContinuumConstitutiveLaw : public Flags {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMContinuumConstitutiveLaw);

        DEMContinuumConstitutiveLaw();

        DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw);

        virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        virtual ~DEMContinuumConstitutiveLaw();

        virtual DEMContinuumConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateContactForces(double mRadius,
                double mSqrtOfRealMass,
                double other_radius,
                double otherSqrtMass,
                double distance,
                double initial_delta,
                int& neighbour_failure_id,
                ProcessInfo& rCurrentProcessInfo,
                PropertiesProxy *myProperties,
                PropertiesProxy *neighbourProperties,
                int mapping_new_ini,
                int mapping_new_cont,
                unsigned int i_neighbour_count,
                double LocalElasticContactForce[3],
                double ViscoDampingLocalContactForce[3],
                double LocalDeltDisp[3],
                Vector mcont_ini_neigh_area,
                array_1d<double, 4 > &mHistory_mapping_new_cont,
                double mDempack_damping,
                int mDampType,
                int mIniNeighbourFailureId_mapping_new_ini,
                double LocalCoordSystem[3][3],
                double RelVel[3]); //FF


        virtual void CalculateViscoDamping(double LocalRelVel[3],
                double ViscoDampingLocalContactForce[3],
                double indentation,
                double equiv_visco_damp_coeff_normal,
                double equiv_visco_damp_coeff_tangential,
                bool sliding,
                int mDampType);




    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

    KRATOS_DEFINE_VARIABLE(DEMContinuumConstitutiveLaw::Pointer, DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER)

    inline void EvaluateFailureCriteria(
            const double contact_sigma,
            const double contact_tau,
            double& failure_criterion_state,
            bool& sliding,
            const int FailureCriterionOption,
            const double TauZero,
            const double TanContactInternalFriccion,
            const double SinContactInternalFriccion,
            const double CosContactInternalFriccion,
            int& NeighbourFailureId_i_neighbour_count, 
            int& IniNeighbourFailureId_mapping_new_ini, 
            const double TensionLimit) {

        //(1) MOHR-COULOMB FAILURE: (we don't consider rotational spring!!!!! here) need to be thought.
        if (FailureCriterionOption == 1) {//MOHR-COULOMB                    
            double sigma_max, sigma_min;

            if (contact_sigma >= 0) {
                sigma_max = contact_sigma;
                sigma_min = 0;
            } else {
                sigma_max = 0;
                sigma_min = contact_sigma;
            }

            //change into principal stresses
            double centre = 0.5 * (sigma_max + sigma_min);
            double radius = sqrt((sigma_max - centre)*(sigma_max - centre) + contact_tau * contact_tau);

            double sigma_I = centre + radius;
            double sigma_II = centre - radius;

            // Check:
            double distance_to_failure = (TauZero / (TanContactInternalFriccion) + centre) * SinContactInternalFriccion;
            failure_criterion_state = radius / distance_to_failure;

            if (sigma_I - sigma_II >= 2 * TauZero * CosContactInternalFriccion + (sigma_I + sigma_II) * SinContactInternalFriccion) {
                //breaks
                NeighbourFailureId_i_neighbour_count = 5; //mohr coulomb   
                IniNeighbourFailureId_mapping_new_ini = 5;
                failure_criterion_state = 1.0;
                sliding = true;
            }
        } //MOHR-COULOMB

        ///(2) UNCOUPLED FRACTURE
        if (FailureCriterionOption == 2) {//UNCOUPLED FRACTURE and DEMPACK        
            //double mTauZero = 0.5*sqrt(mCompressionLimit*mTensionLimit); 

            if (contact_sigma >= 0) {
                double tau_strength = TauZero + TanContactInternalFriccion*contact_sigma;
                failure_criterion_state = contact_tau / tau_strength;

                if (contact_tau > tau_strength) {
                    NeighbourFailureId_i_neighbour_count = 2; //shear in compression
                    IniNeighbourFailureId_mapping_new_ini = 2;
                    failure_criterion_state = 1.0;
                    sliding = true;
                }
            }//positive sigmas

            else { //negative sigmas            
                double tau_strength = TauZero;
                failure_criterion_state = GeometryFunctions::max(contact_tau / tau_strength, -contact_sigma / TensionLimit);

                if (contact_tau > tau_strength) {
                    NeighbourFailureId_i_neighbour_count = 3; //shear failure tension
                    IniNeighbourFailureId_mapping_new_ini = 3;
                    sliding = true;
                    failure_criterion_state = 1.0;                                     
                }
                
            } //negative values of sigma              

        } //UNCOUPLED FRACTURE
        if (failure_criterion_state > 1.0) failure_criterion_state = 1.0;

    }



} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

