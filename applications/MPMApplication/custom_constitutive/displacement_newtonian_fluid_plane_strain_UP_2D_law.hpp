//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Contri Alessandro
//
//  References:      This class is adapted from applications/ParticleMechanicsApplication/custom_constitutive/hyperelastic_UP_plane_strain_2D_law.cpp


#if !defined (KRATOS_DISPLACEMENT_NEWTONIAN_PLANE_STRAIN_UP_2D_LAW_H_INCLUDED)
#define       KRATOS_DISPLACEMENT_NEWTONIAN_PLANE_STRAIN_UP_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/displacement_newtonian_fluid_UP_3D_law.hpp"


namespace Kratos
{

/**
 * Defines a displacement-based newtonian fluid constitutive law
 * This material law is defined by the parameters:
 * 1) DYNAMIC VISCOSITY
 * 2) PRESSURE COEFFICIENT (Cole equation: p=BULK_MODULUS*((rho/rho_0)^(PRESSURE_COEFFICIENT)-1))
 * 1) BULK MODULUS
 * As there are no further parameters the functionality is limited
 * to large displacements.
 */

class KRATOS_API(MPM_APPLICATION) DispNewtonianFluidPlaneStrainUP2DLaw : public DispNewtonianFluidUP3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of DispNewtonianFluidPlaneStrainUP2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( DispNewtonianFluidPlaneStrainUP2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    DispNewtonianFluidPlaneStrainUP2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    DispNewtonianFluidPlaneStrainUP2DLaw (const DispNewtonianFluidPlaneStrainUP2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //DispNewtonianFluidPlaneStrainUP2DLaw& operator=(const DispNewtonianFluidPlaneStrainUP2DLaw& rOther);


    /**
     * Destructor.
     */
    ~DispNewtonianFluidPlaneStrainUP2DLaw() override;

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


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                 Vector& rStrainVector ) override;




    /**
     * Calculates the isochoric constitutive matrix
     * @param rViscousVariables
     * @param rIsoStressVector the isochoric stress vector
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    void CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables& rViscousVariables,
            const Matrix & rIsoStressMatrix,
            Matrix& rConstitutiveMatrix) override;



    /**
     * Calculates the volumetric constitutive matrix
     * @param rViscousVariables
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rViscousVariables,
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DispNewtonianFluidUP3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DispNewtonianFluidUP3DLaw )
    }



}; // Class DispNewtonianFluidPlaneStrainUP2DLaw
}  // namespace Kratos.
#endif // KRATOS_DISPLACEMENT_NEWTONIAN_PLANE_STRAIN_UP_2D_LAW_H_INCLUDED  defined
