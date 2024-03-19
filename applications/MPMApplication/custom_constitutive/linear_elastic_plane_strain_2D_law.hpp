//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/linear_elastic_plane_strain_2D_law.hpp


#if !defined (KRATOS_LINEAR_ELASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define       KRATOS_LINEAR_ELASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED

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

class KRATOS_API(MPM_APPLICATION) LinearElasticPlaneStrain2DLaw : public LinearElastic3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of LinearElasticPlaneStrain2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( LinearElasticPlaneStrain2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    LinearElasticPlaneStrain2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    LinearElasticPlaneStrain2DLaw (const LinearElasticPlaneStrain2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //LinearElasticPlaneStrain2DLaw& operator=(const LinearElasticPlaneStrain2DLaw& rOther);


    /**
     * Destructor.
     */
    ~LinearElasticPlaneStrain2DLaw() override;

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 2;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        return 3;
    };

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


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
    ///@}

    /**
     * Calculates the GreenLagrange strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
            Vector& rStrainVector ) override;


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                         Vector& rStrainVector ) override;


    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */


    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
            const double &rYoungModulus,
            const double &rPoissonCoefficient ) override;



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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }


}; // Class LinearElasticPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined
