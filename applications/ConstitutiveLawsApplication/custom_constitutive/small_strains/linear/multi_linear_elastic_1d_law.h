// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#pragma once

// Project includes
#include "custom_constitutive/truss_constitutive_law.h"

namespace Kratos
{

    /**
     * @namespace MultiLinearElastic1DLaw
     *
     * @brief This constitutive law represents a multi linear elastic 1d claw
     *
     * parameters to define the strain-energy function are
     * MULTI_LINEAR_ELASTICITY_MODULI
     * MULTI_LINEAR_ELASTICITY_STRAINS
     *
     * Johannes Linhard - Numerisch-mechanische Betrachtung des Entwurfsprozesses von Membrantragwerken
     * https://mediatum.ub.tum.de/doc/682189/document.pdf
     *
     * @author Klaus B Sautter
     */


    class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) MultiLinearElastic1DLaw : public TrussConstitutiveLaw
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION( MultiLinearElastic1DLaw );

        /**
         * Default constructor.
         */
        MultiLinearElastic1DLaw();

        ConstitutiveLaw::Pointer Clone() const override;

        /**
         * Copy constructor.
         */
        MultiLinearElastic1DLaw (const MultiLinearElastic1DLaw& rOther);


        /**
         * Destructor.
         */
        ~MultiLinearElastic1DLaw() override;

        /**
         * This function provides the place to perform checks on the completeness of the input.
         * It is designed to be called only once (or anyway, not often) typically at the beginning
         * of the calculations, so to verify that nothing is missing from the input
         * or that no common error is found.
         * @param rMaterialProperties: The properties of the material
         * @param rElementGeometry: The geometry of the element
         * @param rCurrentProcessInfo: The current process info instance
         */
        int Check(
            const Properties& rMaterialProperties,
            const GeometryType& rElementGeometry,
            const ProcessInfo& rCurrentProcessInfo
        ) const override;

        double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
            const Variable<double>& rThisVariable,double& rValue) override;

    private:

        friend class Serializer;

        void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, TrussConstitutiveLaw);
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, TrussConstitutiveLaw);
        }


    }; // Class MultiLinearElastic1DLaw

}  // namespace Kratos.
