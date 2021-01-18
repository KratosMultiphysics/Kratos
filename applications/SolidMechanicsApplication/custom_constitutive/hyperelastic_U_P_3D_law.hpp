//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_U_P_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_U_P_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic isotropic constitutive law in 3D Neohookean Model
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO
 * As there are no further parameters the functionality is limited
 * to large displacements elasticity.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) HyperElasticUP3DLaw : public HyperElastic3DLaw
{
protected:

    /**
     * Parameters to be used in the volumetric and deviatoric split
     */
    struct VectorSplit
    {
        Vector  Isochoric;
        Vector  Volumetric;
    };

    struct MatrixSplit
    {
        Matrix  Isochoric;
        Matrix  Volumetric;
    };


public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticUP3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticUP3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElasticUP3DLaw();


    /**
     * Copy constructor.
     */
    HyperElasticUP3DLaw (const HyperElasticUP3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticUP3DLaw& operator=(const HyperElasticUP3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~HyperElasticUP3DLaw() override;

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    //int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //String Info() const override;
    /**
     * Print information about this object.
     */
    //void PrintInfo(std::ostream& rOStream) const override;
    /**
     * Print object's data.
     */
    //void PrintData(std::ostream& rOStream) const override;

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
     * Calculates the Pressure of the domain (element)
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rPressure the calculated pressure to be returned
     */
    double& CalculateVolumetricPressure (const MaterialResponseVariables & rElasticVariables,
					 double & rPressure) override;


    /**
     * Calculates the Volumetric part factors
     * @param rElasticResponseVariables the material variables
     * @param rFactors Volumetric stress factors
     */
    Vector& CalculateVolumetricPressureFactors (const MaterialResponseVariables & rElasticVariables,
						Vector & rFactors) override;



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


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElastic3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElastic3DLaw )
    }



}; // Class HyperElasticUP3DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_U_P_3D_LAW_H_INCLUDED  defined
