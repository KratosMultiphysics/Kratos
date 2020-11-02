// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Shahed Rezaei
//

#if !defined (KRATOS_LINEAR_PLANE_STRAIN_FRACTURE_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_PLANE_STRAIN_FRACTURE_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_fracture_3d.h"

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
    
/**
 * @class LinearPlaneStrainFracture
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a small deformation linear elastic constitutive model with phase-field-based fracture modelling for plane strain cases
 * @details This class derives from the linear elastic with fracture case on 3D
 * @author Reza Najian Asl
 * @author Shahed Rezaei
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearPlaneStrainFracture 
    : public ElasticIsotropicFracture3D
{
public:
    ///@name Type Definitions
    ///@{

    /// The process info definition
    typedef ProcessInfo      ProcessInfoType;
    
    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw       CLBaseType;
    
    /// The base class ElasticIsotropicFracture3D type definition
    typedef ElasticIsotropicFracture3D      BaseType;
    
    /// The size type definition
    typedef std::size_t             SizeType;
    
    /// Static definition of the dimension
    static constexpr SizeType Dimension = 2;
    
    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 3;
    
    /// Counted pointer of LinearPlaneStrainFracture
    KRATOS_CLASS_POINTER_DEFINITION( LinearPlaneStrainFracture );

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    LinearPlaneStrainFracture();

    /**
     * @brief The clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    LinearPlaneStrainFracture (const LinearPlaneStrainFracture& rOther);


    /**
     * @brief Destructor.
     */
    ~LinearPlaneStrainFracture() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     * @return The dimension were the law is working 
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    SizeType GetStrainSize() override
    {
        return VoigtSize;
    }

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

    /**
     * @brief It calculates the constitutive matrix C
     * @param C The constitutive matrix
     * @param rValues Parameters of the constitutive law
     */
    void CalculateElasticMatrix(
        Matrix& C, 
        ConstitutiveLaw::Parameters& rValues
        ) override;   

    /**
     * @brief It calculates the strain vector
     * @param rValues The internal values of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    void CalculateCauchyGreenStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector
        ) override;

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ElasticIsotropicFracture3D)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ElasticIsotropicFracture3D)
    }


}; // Class LinearPlaneStrainFracture
}  // namespace Kratos.
#endif // KRATOS_LINEAR_PLANE_STRAIN_FRACTURE_LAW_H_INCLUDED  defined
