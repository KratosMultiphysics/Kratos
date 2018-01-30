// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined (KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearPlaneStrain : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of LinearPlaneStrain
     */

    KRATOS_CLASS_POINTER_DEFINITION( LinearPlaneStrain );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    LinearPlaneStrain();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    LinearPlaneStrain (const LinearPlaneStrain& rOther);


    /**
     * Destructor.
     */
    ~LinearPlaneStrain() override;

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
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 3;
    }

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;

    /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1 (Parameters & rValues) override;

    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters & rValues) override;
    
    /**
     * Finalizes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponsePK2 (Parameters & rValues) override;

    /**
     * Finalizes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponseKirchhoff (Parameters & rValues) override;

    /**
     * Finalizes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponsePK1 (Parameters & rValues) override;

    /**
     * Finalizes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */ 
    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override;
    
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
    ) override;


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

private:

    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    void CalculateElasticMatrix(Matrix& C, const double E, const double NU)
    {
        const double c0 = E / ((1.00 + NU)*(1-2*NU));
        const double c1 = (1.00 - NU)*c0;
        const double c2 = c0 * NU;
        const double c3 = (0.5 - NU)*c0;

        C(0,0) = c1;
        C(0,1) = c2;
        C(0,2) = 0.0;
        C(1,0) = c2;
        C(1,1) = c1;
        C(1,2) = 0.0;
        C(2,0) = 0.0;
        C(2,1) = 0.0;
        C(2,2) = c3;
    }

    void CalculatePK2Stress(const Vector& StrainVector, Vector& StressVector, const double E, const double NU)
    {
        const double c0 = E / ((1.00 + NU)*(1-2*NU));
        const double c1 = (1.00 - NU)*c0;
        const double c2 = c0 * NU;
        const double c3 = (0.5 - NU)*c0;


        StressVector[0] = c1*StrainVector[0] + c2 * (StrainVector[1])	;
        StressVector[1] = c1*StrainVector[1] + c2 * (StrainVector[0])	;
        StressVector[2] = c3*StrainVector[2];

    }

    void CalculateCauchyGreenStrain(
        Parameters& rValues,
        Vector& StrainVector
    )
    {
        //1.-Compute total deformation gradient
        const Matrix& F = rValues.GetDeformationGradientF();

        Matrix Etensor = prod(trans(F),F);
        Etensor -= IdentityMatrix(2,2);
        Etensor *= 0.5;

        noalias(StrainVector) = MathUtils<double>::StrainTensorToVector(Etensor);
    }


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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }


}; // Class LinearPlaneStrain
}  // namespace Kratos.
#endif // KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED  defined 
