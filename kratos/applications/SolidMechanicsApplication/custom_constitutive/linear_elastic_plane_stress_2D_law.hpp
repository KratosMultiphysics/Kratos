//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_LINEAR_ELASTIC_PLANE_STRESS_2D_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_PLANE_STRESS_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.hpp"

namespace Kratos
{
/**
 * Defines a linear isotropic constitutive law in 2D (Plane Stress)
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO
 * As there are no further parameters the functionality is valid
 * for small and large displacements elasticity.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) LinearElasticPlaneStress2DLaw : public LinearElasticPlaneStrain2DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of LinearElasticPlaneStress2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( LinearElasticPlaneStress2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    LinearElasticPlaneStress2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    LinearElasticPlaneStress2DLaw (const LinearElasticPlaneStress2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //LinearElasticPlaneStress2DLaw& operator=(const LinearElasticPlaneStress2DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~LinearElasticPlaneStress2DLaw();

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
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */


    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
                                       const double &rYoungModulus,
                                       const double &rPoissonCoefficient );



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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }


}; // Class LinearElasticPlaneStress2DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_PLANE_STRESS_2D_LAW_H_INCLUDED  defined 
