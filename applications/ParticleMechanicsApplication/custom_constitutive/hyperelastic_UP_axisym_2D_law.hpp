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
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/hyperelastic_U_P_axisym_2D_law.hpp


#if !defined (KRATOS_HYPERELASTIC_UP_AXISYM_2D_LAW_H_INCLUDED)
#define       KRATOS_HYPERELASTIC_UP_AXISYM_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_UP_3D_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic isotropic constitutive law in 2D Neohookean Model (Axisymmetric)
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO
 * As there are no further parameters the functionality is limited
 * to large displacements elasticity.
 */

class HyperElasticUPAxisym2DLaw : public HyperElasticUP3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticUPAxisym2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticUPAxisym2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElasticUPAxisym2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    HyperElasticUPAxisym2DLaw (const HyperElasticUPAxisym2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticUPAxisym2DLaw& operator=(const HyperElasticUPAxisym2DLaw& rOther);


    /**
     * Destructor.
     */
    ~HyperElasticUPAxisym2DLaw() override;

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
    SizeType GetStrainSize() override
    {
        return 4;
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
     * Calculates the isochoric constitutive matrix
     * @param rElasticVariables
     * @param rIsoStressVector the isochoric stress vector
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    void CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
            const Matrix & rIsoStressMatrix,
            Matrix& rConstitutiveMatrix) override;


    /**
     * Calculates the volumetric constitutive matrix
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
							Matrix& rConstitutiveMatrix) override;




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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticUP3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticUP3DLaw )
    }



}; // Class HyperElasticUPAxisym2DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_UP_AXISYM_2D_LAW_H_INCLUDED  defined
