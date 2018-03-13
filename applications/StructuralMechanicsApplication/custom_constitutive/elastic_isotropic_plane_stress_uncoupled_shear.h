// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philippe Bussetta
//

#if !defined (KRATOS_ELASTIC_ISOTROPIC_PLANE_STRESS_UNCOUPLED_SHEAR_LAW_H_INCLUDED)
#define  KRATOS_ELASTIC_ISOTROPIC_PLANE_STRESS_UNCOUPLED_SHEAR_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_plane_stress.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ElasticIsotropicPlaneStressUncoupledShear : public LinearPlaneStress
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    /**
     * Counted pointer of ElasticIsotropicPlaneStressUncoupledShear
     */
    KRATOS_CLASS_POINTER_DEFINITION( ElasticIsotropicPlaneStressUncoupledShear );

    ///@name Life Cycle
    ///@{

    /**
     * Default constructor.
     */
    ElasticIsotropicPlaneStressUncoupledShear();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    ElasticIsotropicPlaneStressUncoupledShear(const ElasticIsotropicPlaneStressUncoupledShear& rOther);

    /**
     * Destructor.
     */
    ~ElasticIsotropicPlaneStressUncoupledShear() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

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

    /**
    * It calculates the constitutive matrix C
    * @param C: The constitutive matrix
    * @param E: The Young Modulus
    * @param NU: The poisson coefficient
    */
    void CalculateElasticMatrix(Matrix& C, Parameters& rValues) override;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }


}; // Class ElasticIsotropicPlaneStressUncoupledShear
}  // namespace Kratos.
#endif // KRATOS_ELASTIC_ISOTROPIC_PLANE_STRESS_UNCOUPLED_SHEAR_LAW_H_INCLUDED  defined 
