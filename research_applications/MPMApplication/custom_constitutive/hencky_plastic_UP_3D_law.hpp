//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//


#if !defined (KRATOS_HENCKY_PLASTIC_UP_3D_LAW_H_INCLUDED)
#define       KRATOS_HENCKY_PLASTIC_UP_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hencky_plastic_3D_law.hpp"
#include "includes/ublas_interface.h"

namespace Kratos
{

/**
 * Defines a hencky-plastic isotropic constitutive law in 3D
 * The functionality is limited to large and small displacements
 */


class KRATOS_API(MPM_APPLICATION) HenckyElasticPlasticUP3DLaw : public HenckyElasticPlastic3DLaw
{
//protected:

    //struct MatrixSplit
    //{
    //Matrix  EigenValues;
    //Matrix  EigenVectors;
    //};


public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef MPMFlowRule::Pointer                MPMFlowRulePointer;
    typedef MPMYieldCriterion::Pointer    YieldCriterionPointer;
    typedef MPMHardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of HenckyElasticPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(HenckyElasticPlasticUP3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HenckyElasticPlasticUP3DLaw();


    HenckyElasticPlasticUP3DLaw(MPMFlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /**
     * Copy constructor.
     */
    HenckyElasticPlasticUP3DLaw (const HenckyElasticPlasticUP3DLaw& rOther)  ;


    /**
     * Assignment operator.
     */

    //HenckyElasticPlastic3DLaw& operator=(const HenckyElasticPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~HenckyElasticPlasticUP3DLaw() override;

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
        return 3;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        return 6;
    };

    void GetLawFeatures(Features& rFeatures) override;

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



    void CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables& rElasticVariables) override;

    void GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables);

    void CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen,const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables, const Properties& rProperties) override;

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

    void CalculatePrincipalStressTrial(const MaterialResponseVariables & rElasticVariables,Parameters & rValues,
                                       const MPMFlowRule::RadialReturnVariables& rReturnMappingVariables,
                                       Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix) override;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HenckyElasticPlastic3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HenckyElasticPlastic3DLaw )
    }




}; // Class HenckyElasticPlasticUP3DLaw

} //namespace Kratos

#endif  //KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED

