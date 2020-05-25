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

#if !defined(KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED)
#define KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "fluid_constitutive_law.h"

namespace Kratos {
/**
 * Defines a Newtonian constitutive law for 2D
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY
 */
class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) Newtonian2DLaw : public PfemFluidConstitutiveLaw {
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
    KRATOS_CLASS_POINTER_DEFINITION(Newtonian2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    Newtonian2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    Newtonian2DLaw(const Newtonian2DLaw& rOther);

    /**
     * Destructor.
     */
    ~Newtonian2DLaw() override;

    /**
     * Operations needed by the base class:
     */

    /**
     * @return Working space dimension constitutive law
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @return Size of the strain vector (in Voigt notation) for the constitutive law
     */
    SizeType GetStrainSize() override;

    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

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

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
    ///@}

};  // Class Newtonian2DLaw

}  // namespace Kratos.

#endif  // KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED  defined
