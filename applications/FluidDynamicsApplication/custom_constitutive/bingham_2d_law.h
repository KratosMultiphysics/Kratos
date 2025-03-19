//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//
#pragma once
// System includes

// External includes

// Project includes
#include "fluid_constitutive_law.h"

namespace Kratos
{

/**
 * @class Bingham2DLaw
 * @ingroup FluidDynamicsApplication
 * @brief This law defines a Bingham model for non-Newtonian fluids in two dimensions.
 * @details The Bingham model implementation has been done considering Papanastasiou regularization.
 * @author Uxue Chasco
*/

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) Bingham2DLaw : public FluidConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = 2;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 3;

    /**
     * Counted pointer of Bingham2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(Bingham2DLaw);

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    Bingham2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief Copy constructor.
     */
    Bingham2DLaw (const Bingham2DLaw& rOther);

    /**
     * @brief Destructor.
     */
    ~Bingham2DLaw() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
    return Dimension;
    };

    /**
     * @brief Size of the strain vector (in Voigt notation) for the constitutive law
     */
    SizeType GetStrainSize() const override
    {
    return VoigtSize;
    };

    /**
     * @brief Computes the material response: Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        return "Bingham2DLaw";
    };


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

    /**
     * @brief Get the effective viscosity (in dynamic units -- Pa s) for the fluid.
     */
    double GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const override
    {
        return rParameters.GetConstitutiveMatrix()(3,3);
    };

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

    void save(Serializer& rSerializer) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
    }

}; // Class Bingham2DLaw
}  // namespace Kratos.
