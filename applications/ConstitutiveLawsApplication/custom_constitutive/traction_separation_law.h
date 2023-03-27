// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alireza Taherzadeh-Fard
//                   Alejandro Cornejo Velazquez
//                   Sergio Jimenez Reyes
//                   Lucia Gratiela Barbu
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/composites/rule_of_mixtures_law.h"


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
 * @class TractionSeparationLaw3D
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a parallel rule of mixture (classic law of mixture)
 * @details The constitutive law show have defined a subproperties in order to work properly
 * This law combines parallel CL considering the following principles:
 *  - All layer have the same strain
 *  - The total stress is the addition of the strain in each layer
 *  - The constitutive tensor is the addition of the constitutive tensor of each layer
 * @author Vicente Mataix Ferrandiz
 * @author Fernando Rastellini
 * @author Alejandro Cornejo Velazquez
 */
template<unsigned int TDim>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) TractionSeparationLaw3D
    : public ParallelRuleOfMixturesLaw<TDim>
{
public:

    ///@name Type Definitions
    ///@{
        typedef ParallelRuleOfMixturesLaw<TDim> BaseType;

        typedef std::size_t SizeType;

        /// The define the working dimension size, already defined in the integrator
        static constexpr SizeType VoigtSize = (TDim == 3) ? 6 : 3;

        /// The define the Voigt size, already defined in the  integrator
        static constexpr SizeType Dimension = TDim;

        /// Definition of the machine precision tolerance
        static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

        /// Pointer definition of TractionSeparationLaw3D
        KRATOS_CLASS_POINTER_DEFINITION(TractionSeparationLaw3D);

        ///@name Lyfe Cycle
        ///@{

        /**
         * @brief Default constructor.
         */
        TractionSeparationLaw3D();

        /**
         * @brief Constructor with values
         * @param rCombinationFactors The list of subproperties combination factors
         */
        TractionSeparationLaw3D(const std::vector<double> &rCombinationFactors);

        /**
         * @brief Copy constructor.
         */
        TractionSeparationLaw3D(const TractionSeparationLaw3D &rOther);

        /**
         * @brief Destructor.
         */
        ~TractionSeparationLaw3D() override;

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Clone operation
         */
        ConstitutiveLaw::Pointer Clone() const override;

        /**
         * @brief Creates a new constitutive law pointer
         * @param NewParameters The configuration parameters of the new constitutive law
         * @return a Pointer to the new constitutive law
         */
        ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

        

        

        /**
         * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
         */
        bool RequiresInitializeMaterialResponse() override
        {
            return true;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
      * @brief Is called to check whether the provided material parameters in the Properties match the requirements of current constitutive model.
      * @param rMaterialProperties the current Properties to be validated against.
      * @return true, if parameters are correct; false, if parameters are insufficient / faulty
      * @note  this has to be implemented by each constitutive model. Returns false in base class since no valid implementation is contained here.
      */
    bool ValidateInput(const Properties& rMaterialProperties) override;

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default linear strains)
     * @return the expected strain measure
     */
    ConstitutiveLaw::StrainMeasure GetStrainMeasure() override;

    /**
     * @brief Returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    ConstitutiveLaw::StressMeasure GetStressMeasure() override;

    /**
     * @brief Returns whether this constitutive model is formulated in incremental strains/stresses
     * @note By default, all constitutive models should be formulated in total strains
     * @return true, if formulated in incremental strains/stresses, false otherwise
     */
    bool IsIncremental() override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const ConstitutiveLaw::GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure& rStressMeasure
    );

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

    // std::vector<ConstitutiveLaw::Pointer> mConstitutiveLaws; /// The vector containing the constitutive laws (must be cloned, the ones contained on the properties can conflict between them)
    // std::vector<double> mCombinationFactors;                 /// The vector containing the combination factors of the different layers of the material
    Vector mdelamination_damage_mode_one;
    Vector mdelamination_damage_mode_two;
    Vector mthreshold_mode_one;
    Vector mthreshold_mode_two;
    
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

    // void save(Serializer& rSerializer) const override
    // {
    //     KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    //     rSerializer.save("ConstitutiveLaws", mConstitutiveLaws);
    //     rSerializer.save("CombinationFactors", mCombinationFactors);
    // }

    // void load(Serializer& rSerializer) override
    // {
    //     KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    //     rSerializer.load("ConstitutiveLaws", mConstitutiveLaws);
    //     rSerializer.load("CombinationFactors", mCombinationFactors);
    // }
    
    double MacaullyBrackets(const double Number)
    {
        return (Number > machine_tolerance) ? Number : 0.0;
    }


}; // Class TractionSeparationLaw3D
}  // namespace Kratos.
