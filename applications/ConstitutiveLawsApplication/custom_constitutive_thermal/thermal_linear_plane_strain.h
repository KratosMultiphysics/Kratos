// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined (KRATOS_THERMAL_KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ThermalLinearPlaneStrain
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a Thermo dependant CL, including the addition of thermal expansion strains
 * @details This class derives from the linear elastic case on 3D
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ThermalLinearPlaneStrain
    : public LinearPlaneStrain
{
public:

    ///@name Type Definitions
    ///@{

    typedef LinearPlaneStrain BaseType;

    /// Counted pointer of LinearPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearPlaneStrain);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ThermalLinearPlaneStrain() 
    {
        mGetMaterialValueFunction = [](const Variable<double>& rVariable,ConstitutiveLaw::Parameters& rParameters) -> double {
        const Properties& r_properties = rParameters.GetMaterialProperties();
        if (r_properties.HasTable(TEMPERATURE, rVariable))
            return AdvancedConstitutiveLawUtilities<3>::GetValueFromTable(TEMPERATURE, rVariable, rParameters);
        else
            return r_properties[rVariable];
        };
    }

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ThermalLinearPlaneStrain>(*this);
    }


    /**
     * @brief Destructor.
     */
    ~ThermalLinearPlaneStrain() override {}

    /**
    * Copy constructor.
    */
    ThermalLinearPlaneStrain(const ThermalLinearPlaneStrain &rOther)
        : BaseType(rOther),
          mReferenceTemperature(rOther.mReferenceTemperature)
    {
    }

    double& GetReferenceTemperature()
    {
        return mReferenceTemperature;
    }
    void SetReferenceTemperature(const double ToRefTemperature)
    {
        mReferenceTemperature = ToRefTemperature;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the material response:
     * @details PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief It calculates the value of a specified variable (Vector case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief It calculates the value of a specified variable (double case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(const Properties &rMaterialProperties,
                            const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

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

    double mReferenceTemperature = 0.0;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ElasticIsotropic3D)
        rSerializer.save("ReferenceTemperature", mReferenceTemperature);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ElasticIsotropic3D)
        rSerializer.load("ReferenceTemperature", mReferenceTemperature);
    }


}; // Class LinearPlaneStrain
}  // namespace Kratos.
#endif // KRATOS_THERMAL_KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED  defined
