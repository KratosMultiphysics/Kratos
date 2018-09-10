//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined (KRATOS_BINGHAM_LAW_3D_H_INCLUDED)
#define  KRATOS_BINGHAM_LAW_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "fluid_constitutive_law.h"

namespace Kratos
{
/**
 * Defines a 3D Bingham non-Newtonian constitutive law
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY
 * 2) YIELD_STRESS
 * 3) REGULARIZATION_COEFFICIENT
 */

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) Bingham3DLaw : public FluidConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of Bingham3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(Bingham3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    Bingham3DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    Bingham3DLaw (const Bingham3DLaw& rOther);


    /**
     * Destructor.
     */
    ~Bingham3DLaw() override;

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

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

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


}; // Class Bingham3DLaw
}  // namespace Kratos.
#endif // KRATOS_BINGHAM_LAW_3D_H_INCLUDED  defined 
