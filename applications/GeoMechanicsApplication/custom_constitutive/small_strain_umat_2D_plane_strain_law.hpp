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
#include "small_strain_umat_law.hpp"

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

// Currently, these UMAT constititutive laws are based on the 3D version of the SmallStrainUMAT law
// (using VOIGT_SIZE_3D). This seems counter-intuitive, for 2D laws, but currently this is needed
// because our UMATs are not implemented for 2D plane strain and interface conditions (but expect
// matrix/vector sizes to be consistent with a 3D model). Be careful with changing this, as it may
// lead to UMATs writing to out-of-bounds memory locations. Locally, the static definition of
// VoigtSize is used to ensure copying/using only the necessary data
class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUMAT2DPlaneStrainLaw
    : public SmallStrainUMATLaw<VOIGT_SIZE_3D>
{
public:
    // The base class ConstitutiveLaw type definition
    using BaseType = ConstitutiveLaw;

    /// The size type definition
    using SizeType = std::size_t;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = N_DIM_2D;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = VOIGT_SIZE_2D_PLANE_STRAIN;

    /// Pointer definition of SmallStrainUMAT2DPlaneStrainLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUMAT2DPlaneStrainLaw);

    explicit SmallStrainUMAT2DPlaneStrainLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

    //@}
    //@name Life Cycle
    //@{

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    using SmallStrainUMATLaw::GetValue;

    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using SmallStrainUMATLaw::SetValue;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    SmallStrainUMAT2DPlaneStrainLaw() = default;

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SmallStrainUMAT2DPlaneStrainLaw

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos