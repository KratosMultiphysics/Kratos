// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// System includes
#include "includes/define.h"

// Project includes
#include "custom_constitutive/small_strain_udsm_3D_law.hpp"

namespace Kratos
{
///@addtogroup ConstitutiveModelsApplication
///@{

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

/// Short class definition.
/** Detail class definition.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSM3DInterfaceLaw : public SmallStrainUDSM3DLaw
{
public:
    // The base class ConstitutiveLaw type definition
    using BaseType = ConstitutiveLaw;

    /// The size type definition
    using SizeType = std::size_t;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = N_DIM_3D;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = VOIGT_SIZE_3D_INTERFACE;

    /// Pointer definition of SmallStrainUDSM3DInterfaceLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUDSM3DInterfaceLaw);

    //@}
    //@name Life Cycle
    //@{

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override { return Dimension; }

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override { return VoigtSize; }

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default Green-Lagrange)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override { return StrainMeasure_Infinitesimal; }

    /**
     * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    StressMeasure GetStressMeasure() override { return StressMeasure_Cauchy; }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override { return "SmallStrainUDSM3DInterfaceLaw"; }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "SmallStrainUDSM3DInterfaceLaw Data";
    }

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
    ///@name Protected  Access
    ///@{
    void UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues) override;
    void SetExternalStressVector(Vector& rStressVector) override;
    void SetInternalStressVector(const Vector& rStressVector) override;
    void SetInternalStrainVector(const Vector& rStrainVector) override;
    void CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix) override;

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    indexStress3D getIndex3D(indexStress3DInterface index3D) const;

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SmallStrainUDSM3DLaw

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos