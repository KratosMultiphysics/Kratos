// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#if !defined (KRATOS_THERMAL_DISPERSION_2D_LAW_GEO_H_INCLUDED)
#define  KRATOS_THERMAL_DISPERSION_2D_LAW_GEO_H_INCLUDED

// System includes

// External includes

// Project includes
//#include "custom_constitutive/linear_elastic_plane_strain_K0_law.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "geo_mechanics_application_variables.h"

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
 * @class GeoThermalDispersion2DLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the thermal dispersion for heat cases
 * @details This class derives from the linear elastic case on 3D
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalDispersion2DLaw : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw       CLBaseType;

    /// The base class ElasticIsotropicK03DLaw type definition
    //typedef LinearPlaneStrainK0Law      BaseType;

    /// The size type definition
    typedef std::size_t             SizeType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = N_DIM_2D;

    /// Static definition of the VoigtSize
    // for the time being
    static constexpr SizeType VoigtSize = VOIGT_SIZE_2D_PLANE_STRAIN;

    /// Counted pointer of LinearPlaneStrainK0Law
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalDispersion2DLaw);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    GeoThermalDispersion2DLaw();

    /**
     * @brief The clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    GeoThermalDispersion2DLaw(const GeoThermalDispersion2DLaw& rOther);


    /**
     * @brief Destructor.
     */
    ~GeoThermalDispersion2DLaw() override;

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
    //void GetLawFeatures(Features & rFeatures);

    /**
     * @brief Dimension of the law:
     * @return The dimension were the law is working
     */
    SizeType WorkingSpaceDimension()
    {
        return Dimension;
    }

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    SizeType GetStrainSize() const
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

    /**
     * @brief  Itreturns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue);

    /**
     * @brief It calculates the constitutive matrix C
     * @param C The constitutive matrix
     */
    void CalculateThermalDispersionMatrix(Matrix& C);

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



    // /**
    //  * @brief It calculates the strain vector
    //  * @param rValues The internal values of the law
    //  * @param rStrainVector The strain vector in Voigt notation
    //  */
    // void CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector) override;

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

    // stress vector indices
    // const int VOIGT_INDEX_XX = 0;
    // const int VOIGT_INDEX_YY = 1;

}; // Class LinearPlaneStrainK0Law
}  // namespace Kratos.
#endif // KRATOS_THERMAL_DISPERSION_2D_LAW_GEO_H_INCLUDED  defined
