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

#if !defined (KRATOS_THERMAL_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"
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
 * @class ThermalElasticIsotropic3D
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a Thermo dependant CL, including the addition of thermal expansion strains
 * @details This class derives from the linear elastic case on 3D
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ThermalElasticIsotropic3D
    : public ElasticIsotropic3D
{
public:

    ///@name Type Definitions
    ///@{

    /// Counted pointer of LinearPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION(ThermalElasticIsotropic3D);

    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
    * Copy constructor.
    */
    ThermalElasticIsotropic3D(const ThermalElasticIsotropic3D &rOther)
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

    /**
     * @brief This method retrieves a material property, in thermal CL
     * this method calculate the value of the property according to TEMPERATURE
     * @param rVariable The property of the material to be retrieved
     * @param rValues The constitutive parameters
     */
    double GetMaterialProperty(
        const Variable<double> &rVariable,
        ConstitutiveLaw::Parameters &rParameters
        ) override;
    {
        // return rValues.GetMaterialProperties()[rVariable];
        const Properties& r_properties = rParameters.GetMaterialProperties();
        double variable_value;

        if (r_properties.HasTable(TEMPERATURE, rVariable))
            variable_value = this->GetValueFromTable(TEMPERATURE, rVariable, rParameters);
        else
            variable_value = r_properties[rVariable];
        return variable_value;
    }

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
#endif // KRATOS_THERMAL_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED  defined
