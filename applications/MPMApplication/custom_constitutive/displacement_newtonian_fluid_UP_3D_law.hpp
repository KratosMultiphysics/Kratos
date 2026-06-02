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
//  References:      This class is adapted from applications/ParticleMechanicsApplication/custom_constitutive/hyperelastic_UP_3D_law.hpp


#if !defined (KRATOS_DISPLACEMENT_NEWTONIAN_UP_3D_LAW_H_INCLUDED)
#define       KRATOS_DISPLACEMENT_NEWTONIAN_UP_3D_LAW_H_INCLUDED

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

class KRATOS_API(MPM_APPLICATION) DispNewtonianFluidUP3DLaw : public DispNewtonianFluid3DLaw
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
     * Counted pointer of DispNewtonianFluidUP3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( DispNewtonianFluidUP3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    DispNewtonianFluidUP3DLaw();


    /**
     * Copy constructor.
     */
    DispNewtonianFluidUP3DLaw (const DispNewtonianFluidUP3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //DispNewtonianFluidUP3DLaw& operator=(const DispNewtonianFluidUP3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~DispNewtonianFluidUP3DLaw() override;

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
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters & rValues) override;


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
    double& CalculateVolumetricPressure (const MaterialResponseVariables & rViscousVariables,
					 double & rPressure) override;


    /**
     * Calculates the Volumetric part factors
     * @param rViscousVariables the material variables
     * @param rFactors Volumetric stress factors
     */
    Vector& CalculateVolumetricPressureFactors (const MaterialResponseVariables & rViscousVariables,
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DispNewtonianFluid3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DispNewtonianFluid3DLaw )
    }



}; // Class DispNewtonianFluidUP3DLaw
}  // namespace Kratos.
#endif // KRATOS_DISPLACEMENT_NEWTONIAN_UP_3D_LAW_H_INCLUDED  defined
