//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Massimiliano Zecchetto
//  Collaborators:
//
//-------------------------------------------------------------
//

#if !defined(KRATOS_HYPOELASTIC_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED)
#define KRATOS_HYPOELASTIC_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/solid_laws/hypoelastic_2D_law.h"

namespace Kratos {
/**
 * Defines an Hypoelastic constitutive law for 2D
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO
 */
class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) HypoelasticTemperatureDependent2DLaw : public Hypoelastic2DLaw {
   public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    /**
     * Counted pointer of HypoelasticTemperatureDependent2DLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(HypoelasticTemperatureDependent2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HypoelasticTemperatureDependent2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    HypoelasticTemperatureDependent2DLaw(const HypoelasticTemperatureDependent2DLaw& rOther);

    /**
     * Destructor.
     */
    ~HypoelasticTemperatureDependent2DLaw() override;

    /**
     * Operations needed by the base class:
     */

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override;

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

    /// Get the effective density for the solid.
    double GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const override;

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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
    ///@}

};  // Class HypoelasticTemperatureDependent2DLaw

}  // namespace Kratos.

#endif  // KRATOS_HYPOELASTIC_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED  defined
