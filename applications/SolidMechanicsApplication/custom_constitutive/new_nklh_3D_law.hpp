//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:               DAbadias $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_NEW_NKLH_3D_LAW_H_INCLUDED)
#define  KRATOS_NEW_NKLH_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_3D_law.hpp"

namespace Kratos
{
/**
 * Defines a linear isotropic constitutive law in 2D (Plane Strain)
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO 
 * As there are no further parameters the functionality is valid
 * for small and large displacements elasticity.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) NewNKLH3DLaw : public LinearElastic3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of NewNKLH3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(NewNKLH3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    NewNKLH3DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    NewNKLH3DLaw (const NewNKLH3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //LinearElastic3DLaw& operator=(const LinearElastic3DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~NewNKLH3DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& rMaterialProperties,
                             const GeometryType& rElementGeometry,
                             const Vector& rShapeFunctionsValues ) override;

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


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
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    Vector mStressPrevious;

    Vector mStrainPrevious;

    Vector mHistoricalVariablesPrevious;

    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}




    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */


    void CalculateLinearElasticStiffnessMatrix( Matrix& rConstitutiveMatrix,
            const double &rYoungModulus,
            const double &rPoissonCoefficient );

    /**
      * This function is designed to be called when before the material response
      * to check if all needed parameters for the constitutive are initialized
      * @param Parameters
      * @return
      */
    bool CheckParameters(Parameters& rValues) override;


    /**
      */
    Matrix & ConvertConstitutiveMatrixAppropiateSize(Matrix & rOutput, const Matrix & rInput, const int voigtSize);

    Matrix & ConvertStrainTensorTo3D( Matrix & StrainTensor);

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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }


}; // Class LinearElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined 
