//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Ruben Zorilla
//  Collaborators:  Massimiliano Zecchetto
//
//-------------------------------------------------------------
//

#if !defined(KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED)
#define KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "newtonian_2D_law.h"

namespace Kratos {
/**
 * Defines a temperature dependent Newtonian constitutive law for 2D
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY
 */
class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) NewtonianTemperatureDependent2DLaw : public Newtonian2DLaw {
   public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    /**
     * Counted pointer of Newtonian3DLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(NewtonianTemperatureDependent2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    NewtonianTemperatureDependent2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    NewtonianTemperatureDependent2DLaw(const NewtonianTemperatureDependent2DLaw& rOther);

    /**
     * Destructor.
     */
    ~NewtonianTemperatureDependent2DLaw() override;

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
              const ProcessInfo& rCurrentProcessInfo) override;

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

    /// Get the effective viscosity (in dynamic units -- Pa s) for the fluid.
    double GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const override;

    /// Get the effective density for the fluid.
    double GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const override;

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

};  // Class NewtonianTemperatureDependent2DLaw

}  // namespace Kratos.

#endif  // KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED  defined
