// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:        BSD License
//                  license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"

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
 * @class JointedCohesiveFrictionalConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a CL for continuum rocks with an embedded joint
 * @details This CL is detailed in Modelling jointed rock mass as a continuum with an embedded cohesive-frictional model
 * Engineering Geology 228 (2017) 107–120, DOI:10.1016/j.enggeo.2017.07.011
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) JointedCohesiveFrictionalConstitutiveLaw
    : public LinearPlaneStrain
{
public:

    ///@name Type Definitions
    ///@{
    /// The define the working dimension size, already defined in the integrator
    /// The definition of the size type
    typedef std::size_t SizeType;

    static constexpr SizeType Dimension = 2;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = 3;

    /// The definition of the process info
    typedef ProcessInfo ProcessInfoType;

    /// The definition of the CL base  class
    using BaseType = LinearPlaneStrain;

    /// Definition of the machine precision tolerance
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

    static constexpr double tolerance = 1.0e-8;
    static constexpr double ratio_tolerance = 1.0e-4;

    /// The node definition
    using NodeType = Node;

    /// The geometry definition
    using GeometryType = Geometry<NodeType>;

    /// The definition of the Voigt array type
    using BoundedVectorType = array_1d<double, VoigtSize>;

    /// The definition of the bounded matrix type
    using BoundedMatrixType = BoundedMatrix<double, VoigtSize, VoigtSize>;

    /// Pointer definition of JointedCohesiveFrictionalConstitutiveLaw
    KRATOS_CLASS_POINTER_DEFINITION(JointedCohesiveFrictionalConstitutiveLaw);

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    JointedCohesiveFrictionalConstitutiveLaw()
    {}

    /**
     * @brief Destructor.
     */
    ~JointedCohesiveFrictionalConstitutiveLaw() override
    {}

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
    * Copy constructor.
    */
    JointedCohesiveFrictionalConstitutiveLaw(const JointedCohesiveFrictionalConstitutiveLaw &rOther)
        : BaseType(rOther)
        //   mPlasticDissipation(rOther.mPlasticDissipation),
        //   mDamageDissipation(rOther.mDamageDissipation),
        //   mThreshold(rOther.mThreshold),
        //   mPlasticStrain(rOther.mPlasticStrain),
        //   mOldStrain(rOther.mOldStrain),
        //   mComplianceMatrix(rOther.mComplianceMatrix),
        //   mComplianceMatrixCompression(rOther.mComplianceMatrixCompression)
    {
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

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
     * @brief Sets the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
     void SetValue(
        const Variable<double>& rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector >& rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;


    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;


    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(
        ConstitutiveLaw::ConstitutiveLaw::Parameters& rValues);


    /**
     * @brief This method add somehting if the increment is positive
     */
    void AddIfPositive(
        double&  rToBeAdded,
        const double Increment
        )
    {
        if (Increment > machine_tolerance)
            rToBeAdded += Increment;
    }

    /**
     * @brief This method evaluates the Macaulay brackets
     */
    double MacaullyBrackets(const double Number)
    {
        return (Number > machine_tolerance) ? Number : 0.0;
    }

    double Heaviside(const double Number)
    {
        return (Number >= 0) ? 1.0 : 0.0;
    }

    /**
     * @brief This method returns the value of the yield surface
     * defined in Eq. (12) in Lihn et al.
     */
    double YieldSurfaceValue(
        const double ts,
        const double D,
        const double mu0,
        const double mu,
        const double tn,
        const double ft,
        const double m,
        const double fc) const
    {
        const double one_minus_D = 1.0 - D;
        return std::pow(ts, 2) - (one_minus_D * std::pow(mu0, 2) + D * std::pow(mu, 2)) * std::pow((tn - one_minus_D * ft), 2) +
               m * fc * one_minus_D * (tn - one_minus_D * ft);
    }

    /**
     * @brief This method returns the normal n and rotation matrix R
     * Estimates the orientation of the crack by checking the maximum
     * yield stress with respect to theta angle
     */
    void CalculateJointOrientation(
        Matrix& rR,
        Matrix& rn,
        const Vector& rStressVectorTrial,
        const double D,
        const double mu0,
        const double mu,
        const double ft,
        const double m,
        const double fc)
    {
        const double pi_over_180 = Globals::Pi / 180.0;
        double theta = 0.0, y_max = -std::numeric_limits<double>::infinity(), theta_incr = 0.1;
        Matrix n_trial(VoigtSize, Dimension), R_trial(Dimension, Dimension);
        Vector tc_trial = ZeroVector(Dimension);
        n_trial.clear();
        R_trial.clear();

        if (rR.size1() != Dimension || rR.size1() != Dimension)
            rR.resize(Dimension, Dimension, false);

        if (rn.size1() != VoigtSize || rn.size1() != Dimension)
            rn.resize(VoigtSize, Dimension, false);

        while (theta <= 180.0) {
            const double angle_radians = theta * pi_over_180;
            const double n1 = std::cos(angle_radians);
            const double n2 = std::sin(angle_radians);
            const double m1 = -n2;
            const double m2 = n1;

            R_trial(0, 0) = n1;
            R_trial(0, 1) = n2;
            R_trial(1, 0) = m1;
            R_trial(1, 1) = m2;

            n_trial(0, 0) = n1;
            n_trial(1, 1) = n2;
            n_trial(2, 0) = n2;
            n_trial(2, 1) = n1;

            noalias(tc_trial) = prod(R_trial, Vector(prod(trans(n_trial), rStressVectorTrial)));
            const double y_trial = YieldSurfaceValue(tc_trial[1], D, mu0, mu, tc_trial[0], ft, m, fc);

            if (y_trial > y_max) {
                y_max = y_trial;
                noalias(rn) = n_trial;
                noalias(rR) = R_trial;
            }
            theta += theta_incr;
        }
    }

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

    double mDamage = 0.0;                            // Damage Variable "D"
    Vector mLocalTraction = ZeroVector(Dimension);   // The traction vector in local axes of the joint "tc"
    Vector mUc  = ZeroVector(Dimension);             // The local displacement jump "[[uc]]"
    Vector mUcp = ZeroVector(Dimension);             // The local plastic displacement jump "[[ucp]]" 
    Vector mOldStrainVector = ZeroVector(VoigtSize); // The previous converged strain vector
    Vector mOldStressVector = ZeroVector(VoigtSize); // The previous converged stress vector
    double mUp = 0.0;                                // 1d version of the plastic displacements
    bool mDoubleScale = false;                       // Indicates if the crack path has been opened

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

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        // rSerializer.save("PlasticDissipation", mPlasticDissipation);
        // rSerializer.save("DamageDissipation", mDamageDissipation);
        // rSerializer.save("Threshold", mThreshold);
        // rSerializer.save("PlasticStrain", mPlasticStrain);
        // rSerializer.save("OldStrain", mOldStrain);
        // rSerializer.save("ComplianceMatrix", mComplianceMatrix);
        // rSerializer.save("ComplianceMatrixCompression", mComplianceMatrixCompression);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        // rSerializer.load("PlasticDissipation", mPlasticDissipation);
        // rSerializer.load("DamageDissipation", mDamageDissipation);
        // rSerializer.load("Threshold", mThreshold);
        // rSerializer.load("PlasticStrain", mPlasticStrain);
        // rSerializer.load("OldStrain", mOldStrain);
        // rSerializer.load("ComplianceMatrix", mComplianceMatrix);
        // rSerializer.load("ComplianceMatrixCompression", mComplianceMatrixCompression);
    }


}; // Class JointedCohesiveFrictionalConstitutiveLaw
}  // namespace Kratos.
