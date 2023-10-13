// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#pragma once


// Project includes
#include "custom_constitutive/linear_plane_stress.h"

namespace Kratos
{

    /**
     * @namespace MultiLinearIsotropicPlaneStress2D
     *
     * @brief This constitutive law represents a multi linear elastic 2d claw for plane stress
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
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) MultiLinearIsotropicPlaneStress2D
    : public LinearPlaneStress
{
public:

    /// Counted pointer of MultiLinearIsotropicPlaneStress2D
    KRATOS_CLASS_POINTER_DEFINITION( MultiLinearIsotropicPlaneStress2D );


    /**
     * Default constructor.
     */
    MultiLinearIsotropicPlaneStress2D();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    MultiLinearIsotropicPlaneStress2D (const MultiLinearIsotropicPlaneStress2D& rOther);


    /**
     * Destructor.
     */
    ~MultiLinearIsotropicPlaneStress2D() override;

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



protected:

    /**
    * It calculates the constitutive matrix C
    * @param C: The constitutive matrix
    * @param rValues Parameters of the constitutive law
    */
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues) override;

    /**
    * It calculates the stress vector
    * @param rStrainVector The strain vector in Voigt notation
    * @param rStressVector The stress vector in Voigt notation
    * @param rValues Parameters of the constitutive law
    */
    void CalculatePK2Stress(
        const Vector& rStrainVector,
        Vector& rStressVector,
        ConstitutiveLaw::Parameters& rValues
        ) override;


private:

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearPlaneStress)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearPlaneStress)
    }


}; // Class MultiLinearIsotropicPlaneStress2D
}  // namespace Kratos.
