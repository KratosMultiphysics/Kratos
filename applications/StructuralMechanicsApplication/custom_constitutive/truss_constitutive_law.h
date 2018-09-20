// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#if !defined (KRATOS_TRUSS_CONSTITUTIVE_LAW_H_INCLUDED)
#define  KRATOS_TRUSS_CONSTITUTIVE_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{

/**
 * @namespace TrussConstitutiveLaw
 *
 * @brief This constitutive law represents a linear elastic 1D law
 *
 * @author Klaus B Sautter
 */


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TrussConstitutiveLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of TrussConstitutiveLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( TrussConstitutiveLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    TrussConstitutiveLaw();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    TrussConstitutiveLaw (const TrussConstitutiveLaw& rOther);


    /**
     * Destructor.
     */
    ~TrussConstitutiveLaw() override;

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
        return 1;
    }

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

    array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rThisVariable,
        array_1d<double, 3 > & rValue) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,double& rValue) override;

    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue) override;

    array_1d<double, 3 > & CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue) override;

    void CalculateMaterialResponse(
        const Vector& rStrainVector,const Matrix& rDeformationGradient,
        Vector& rStressVector,Matrix& rAlgorithmicTangent,
        const ProcessInfo& rCurrentProcessInfo,const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues,
        bool CalculateStresses,int CalculateTangent,bool SaveInternalVariables) override;


    //empty because called in the element and this base class throws an error
    //if this is not overriden
    void FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                    const GeometryType& rElementGeometry,
                    const Vector& rShapeFunctionsValues,
                    const ProcessInfo& rCurrentProcessInfo) override {} ;

    //this functions calculates the current stress based on an element given (set)
    //strain
    double CalculateStressElastic(ConstitutiveLaw::Parameters& rParameterValues) const;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw);
    }


}; // Class TrussConstitutiveLaw
}  // namespace Kratos.
#endif // KRATOS_DUMMY_TRUSS_LAW_H_INCLUDED  defined
