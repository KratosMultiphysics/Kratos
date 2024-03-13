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
//  References:      This class is adapted from applications/MPMApplication/custom_constitutive/hyperelastic_plane_strain_2D_law.hpp


#if !defined (KRATOS_DISPLACEMENT_NEWTONIAN_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define       KRATOS_DISPLACEMENT_NEWTONIAN_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/displacement_newtonian_fluid_3D_law.hpp"


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

class KRATOS_API(MPM_APPLICATION) DispNewtonianFluidPlaneStrain2DLaw
    : public DispNewtonianFluid3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of DispNewtonianFluidPlaneStrain2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( DispNewtonianFluidPlaneStrain2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    DispNewtonianFluidPlaneStrain2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    DispNewtonianFluidPlaneStrain2DLaw (const DispNewtonianFluidPlaneStrain2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //DispNewtonianFluidPlaneStrain2DLaw& operator=(const DispNewtonianFluidPlaneStrain2DLaw& rOther);


    /**
     * Destructor.
     */
    ~DispNewtonianFluidPlaneStrain2DLaw() override;

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
    virtual void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                         Vector& rStrainVector ) override;


    /**
     * Calculates the constitutive matrix
     * @param rViscousVariables
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    virtual void CalculateConstitutiveMatrix ( const MaterialResponseVariables& rViscousVariables,
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DispNewtonianFluid3DLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DispNewtonianFluid3DLaw)
    }



}; // Class DispNewtonianFluidPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_DISPLACEMENT_NEWTONIAN_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined
