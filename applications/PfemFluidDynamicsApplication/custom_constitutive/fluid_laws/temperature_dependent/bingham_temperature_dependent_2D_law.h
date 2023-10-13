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

#if !defined(KRATOS_BINGHAM_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED)
#define KRATOS_BINGHAM_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/fluid_laws/bingham_2D_law.h"

namespace Kratos
{
    /**
     * Defines a 2D Bingham non-Newtonian constitutive law
     * This material law is defined by the parameters:
     * 1) DYNAMIC_VISCOSITY
     * 2) YIELD_SHEAR
     * 3) ADAPTIVE_EXPONENT
     */
    class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) BinghamTemperatureDependent2DLaw : public Bingham2DLaw
    {
    public:
        /**
         * Type Definitions
         */
        typedef ProcessInfo ProcessInfoType;
        typedef ConstitutiveLaw BaseType;
        typedef std::size_t SizeType;

        /**
         * Counted pointer of BinghamTemperatureDependent2DLaw
         */
        KRATOS_CLASS_POINTER_DEFINITION(BinghamTemperatureDependent2DLaw);

        /**
         * Life Cycle
         */

        /**
         * Default constructor.
         */
        BinghamTemperatureDependent2DLaw();

        /**
         * Clone function (has to be implemented by any derived class)
         * @return a pointer to a new instance of this constitutive law
         */
        ConstitutiveLaw::Pointer Clone() const override;

        /**
         * Copy constructor.
         */
        BinghamTemperatureDependent2DLaw(const BinghamTemperatureDependent2DLaw &rOther);

        /**
         * Destructor.
         */
        ~BinghamTemperatureDependent2DLaw() override;

        /**
         * Operators
         */

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
        int Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                  const ProcessInfo &rCurrentProcessInfo) const override;

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

        ///@}
        ///@name Private  Access
        ///@{
        ///@}

        ///@}
        ///@name Serialization
        ///@{
        friend class Serializer;

        void save(Serializer &rSerializer) const override;

        void load(Serializer &rSerializer) override;
        ///@}

    }; // Class BinghamTemperatureDependent2DLaw

} // namespace Kratos.

#endif // KRATOS_BINGHAM_TEMPERATURE_DEPENDENT_2D_LAW_H_INCLUDED  defined
