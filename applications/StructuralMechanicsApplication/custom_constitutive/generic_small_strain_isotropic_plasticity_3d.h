// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo


#if !defined (KRATOS_GENERIC_SMALL_STRAIN_ISOTROPIC_PLASTICITY_3D_H_INCLUDED)
#define  KRATOS_GENERIC_SMALL_STRAIN_ISOTROPIC_PLASTICITY_3D_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"

#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

    template <class YieldSurfaceType, class ConstLawIntegratorType>
    class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericSmallStrainIsotropicPlasticity3D
     : public ConstitutiveLaw
    {

        public:
            KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainIsotropicPlasticity3D);
            /**  
            * Default constructor.
            */
            GenericSmallStrainIsotropicPlasticity3D()
            {
                //mpYieldSurface = YieldSurfaceType().Clone();
                //mpConstLawIntegrator = ConstLawIntegratorType().Clone();
            }

            /**  
            * Clone.
            */
            ConstitutiveLaw::Pointer Clone() const override
            {
                GenericSmallStrainIsotropicPlasticity3D<YieldSurfaceType, ConstLawIntegratorType>::Pointer p_clone
                    (new GenericSmallStrainIsotropicPlasticity3D< YieldSurfaceType,  ConstLawIntegratorType>(*this));
                return p_clone;
            }

            /**
            * Copy constructor.
            */
            GenericSmallStrainIsotropicPlasticity3D (const GenericSmallStrainIsotropicPlasticity3D& rOther)
            : ConstitutiveLaw(rOther)
            {
            }
            /**
            * Destructor.
            */
            ~GenericSmallStrainIsotropicPlasticity3D() override
            {
            }

            int GetVoigtSize(){return 6;}
            int GetWorkingSpaceDimension() {return 3;} 

            double GetThreshold() {return mThreshold;}
            double GetPlasticDissipation() {return mPlasticDissipation;}
            Vector GetPlasticStrain(){return mPlasticStrain;}

            double GetNonConvThreshold() {return mNonConvThreshold;}
            double GetNonConvPlasticDissipation() {return mNonConvPlasticDissipation;}
            Vector GetNonConvPlasticStrain(){return mNonConvPlasticStrain;}

            void SetThreshold(const double& toThreshold) {mThreshold = toThreshold;}
            void SetPlasticDissipation(const double& toCapap) {mPlasticDissipation = toCapap;}
            void SetPlasticStrain(const Vector& Ep){mPlasticStrain = Ep;}

            void SetNonConvThreshold(const double& toThreshold) {mNonConvThreshold = toThreshold;}
            void SetNonConvPlasticDissipation(const double& toCapap) {mNonConvPlasticDissipation = toCapap;}
            void SetNonConvPlasticStrain(const Vector& Ep){mNonConvPlasticStrain = Ep;}

        private:
            // Converged values
            double mPlasticDissipation = 0.0;
            double mThreshold = 0.0;
            Vector mPlasticStrain = ZeroVector(this->GetVoigtSize());

            // Non Converged values
            double mNonConvPlasticDissipation = 0.0;
            double mNonConvThreshold = 0.0;
            Vector mNonConvPlasticStrain = ZeroVector(this->GetVoigtSize());

            friend class Serializer;

            void save(Serializer& rSerializer) const override
            {
                KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
            }

            void load(Serializer& rSerializer) override
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
            }

        protected:

            //typename YieldSurfaceType::Pointer mpYieldSurface;
			//typename ConstLawIntegratorType::Pointer mpConstLawIntegrator;

            //YieldSurfaceType mYieldSurface;

            //YieldSurfaceType::Calculatebabbaa


    }; // class

} // namespace kratos
#endif