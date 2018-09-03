//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined (KRATOS_EULER_LAW_2D_H_INCLUDED)
#define  KRATOS_EULER_LAW_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "fluid_constitutive_law.h"

namespace Kratos
{
    /**
     * Defines an Euler fluid constitutive law
     * This material law represents a null shear stress contribution.
     */
    class KRATOS_API(FLUID_DYNAMICS_APPLICATION) Euler2DLaw : public FluidConstitutiveLaw
    {
    public:
        /**
         * Type Definitions
         */
        typedef ProcessInfo      ProcessInfoType;
        typedef ConstitutiveLaw         BaseType;
        typedef std::size_t             SizeType;

        /**
         * Counted pointer of Euler2DLaw
         */
        KRATOS_CLASS_POINTER_DEFINITION(Euler2DLaw);

        /**
         * Life Cycle
         */

        /**
         * Default constructor.
         */
        Euler2DLaw();

        /**
         * Clone function (has to be implemented by any derived class)
         * @return a pointer to a new instance of this constitutive law
         */
        ConstitutiveLaw::Pointer Clone() const override;

        /**
         * Copy constructor.
         */
        Euler2DLaw (const Euler2DLaw& rOther);

        /**
         * Destructor.
         */
        ~Euler2DLaw() override;

        /**
         * Operators
         */

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

        void CalculateMaterialResponseCauchy (Parameters& rValues) override;

        int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

        /**
         * Input and output
         */

        /**
         * Turn back information as a string.
         */
        std::string Info() const override {
            std::stringstream buffer;
            buffer << "Euler2DLaw";
            return buffer.str();
        };

        /**
         * Print information about this object.
         */
        void PrintInfo(std::ostream& rOStream) const override {
            rOStream << Info();
        };

        /**
         * Print object's data.
         */
        void PrintData(std::ostream& rOStream) const override {
            rOStream << Info();
        };

    protected:

        ///@name Protected Operations
        ///@{

        /// Get the effective viscosity (in dynamic units -- Pa s) for the fluid.
        double GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const override;

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

    }; // Class Euler2DLaw
}  // namespace Kratos.
#endif // KRATOS_EULER_LAW_2D_H_INCLUDED  defined 
