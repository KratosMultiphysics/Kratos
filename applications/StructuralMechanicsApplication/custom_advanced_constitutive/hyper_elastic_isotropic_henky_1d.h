// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#if !defined (KRATOS_HYPER_ELASTIC_ISOTROPIC_HENKY_1D_LAW_H_INCLUDED)
#define  KRATOS_HYPER_ELASTIC_ISOTROPIC_HENKY_1D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{

/**
 * @namespace HyperElasticIsotropicHenky1D
 *
 * @brief This constitutive law represents the hyper-elastic henky 1D law
 *
 * @author Klaus B Sautter
 */


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) HyperElasticIsotropicHenky1D : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticIsotropicHenky1D
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticIsotropicHenky1D );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElasticIsotropicHenky1D();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    HyperElasticIsotropicHenky1D (const HyperElasticIsotropicHenky1D& rOther);


    /**
     * Destructor.
     */
    ~HyperElasticIsotropicHenky1D() override;

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 1;
    }

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
    ) override;

    array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rThisVariable,
        array_1d<double, 3 > & rValue) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,double& rValue) override;

    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue) override;

    array_1d<double, 3 > & CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue) override;

    void CalculateMaterialResponsePK2(Parameters& rValues) override;


    void FinalizeMaterialResponsePK2(Parameters& rValues) override
    {
        // plasticity law needs this function, so it is called in the truss element
    };

    //this functions calculates the current stress based on an element given (set)
    //strain
    double CalculateStressElastic(ConstitutiveLaw::Parameters& rParameterValues) const;

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    ///@}


    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw);
    }


}; // Class HyperElasticIsotropicHenky1D
}  // namespace Kratos.
#endif // KRATOS_HYPER_ELASTIC_ISOTROPIC_HENKY_1D_LAW_H_INCLUDED  defined
