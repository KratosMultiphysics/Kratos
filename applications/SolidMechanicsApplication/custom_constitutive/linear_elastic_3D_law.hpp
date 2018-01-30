//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"

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

class KRATOS_API(SOLID_MECHANICS_APPLICATION) LinearElastic3DLaw : public HyperElastic3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of LinearElastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    LinearElastic3DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    LinearElastic3DLaw (const LinearElastic3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //LinearElastic3DLaw& operator=(const LinearElastic3DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~LinearElastic3DLaw();

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
    void GetLawFeatures(Features& rFeatures);


    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues);


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);

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
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}


    /**
     * Calculates the stresses for given strain state
     * @param rStrainVector
     * @param rConstitutiveMatrix
     * @param rStressVector the stress vector corresponding to the deformation
     */
    virtual void CalculateStress( const Vector &rStrainVector,
                                  const Matrix &rConstitutiveMatrix,
                                  Vector& rStressVector);


    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */


    virtual void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
            const double &rYoungModulus,
            const double &rPoissonCoefficient );

    /**
      * This function is designed to be called when before the material response
      * to check if all needed parameters for the constitutive are initialized
      * @param Parameters
      * @return
      */
    bool CheckParameters(Parameters& rValues);
    
    virtual double& GetValue( const Variable<double>& rThisVariable, double& rValue );


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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElastic3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElastic3DLaw )
    }


}; // Class LinearElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined 
