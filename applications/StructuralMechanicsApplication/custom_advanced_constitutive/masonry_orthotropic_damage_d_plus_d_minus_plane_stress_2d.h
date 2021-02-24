// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philip Kalkbrenner
//                   Massimo Petracca
//                   Alejandro Cornejo
//  
//

#if !defined(KRATOS_MASONRY_ORTHOTROPIC_DAMAGE_D_PLUS_D_MINUS_PLANE_STRESS_2D_H_INCLUDED )
#define  KRATOS_MASONRY_ORTHOTROPIC_DAMAGE_D_PLUS_D_MINUS_PLANE_STRESS_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
    /**
    * @class MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw
    * @ingroup StructuralMechanicsApplication
    */
    class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw
        : public ConstitutiveLaw
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw);

        ///@name Type Definitions
        ///@{

        /// We define the working dimension size, already defined in the integrator
        static constexpr SizeType Dimension = 2;

        /// We define the Voigt size, already defined in the  integrator
        static constexpr SizeType VoigtSize = 3;

        /// Definition of the machine precision tolerance
        static constexpr double tolerance = std::numeric_limits<double>::epsilon();

        ///@}

        /// Default Constructor.
        MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw()
        {}

        /// Destructor.
        ~MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw() override
        {}

        /// Clone
        ConstitutiveLaw::Pointer Clone() const override
        {
            return ConstitutiveLaw::Pointer(new MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw());
        }

        /// @brief returns the working space dimension of the current constitutive law
        SizeType WorkingSpaceDimension() override
        {
            return Dimension;
        };

        /// @brief returns the size of the strain vector of the current constitutive law
        SizeType GetStrainSize() override
        {
            return VoigtSize;
        };

        struct DirectionalMaterialProperties {
            double E;  // Youngs modulus; 1->11; 2->22
            double nu; // Poissons ratio; 1->12; 2->21
            double G;  // Shear modulus;  1->12; 2->21

            double YieldStressTension;
            double FractureEnergyTension;

            double ElasticLimitStressCompression;
            double YieldStressCompression;
            double YieldStrainCompression;
            double ResidualStressCompression;
            double FractureEnergyCompression;

            double BezierControllerC1;
            double BezierControllerC2;
            double BezierControllerC3;

            double YieldStressShear;
            double YieldStressShearMinus;

            double CharacteristicLength;

            DirectionalMaterialProperties(
                double CharacteristicLength,
                double E, double nu, double G,
                double YieldStressShear,
                double YieldStressShearMinus,
                double YieldStressTension,
                double FractureEnergyTension,
                double ElasticLimitStressCompression,
                double YieldStressCompression,
                double YieldStrainCompression,
                double ResidualStressCompression,
                double FractureEnergyCompression,
                double BezierControllerC1 = 0.65,
                double BezierControllerC2 = 0.5,
                double BezierControllerC3 = 1.5)
                : CharacteristicLength(CharacteristicLength),
                E(E), nu(nu), G(G),
                YieldStressShear(YieldStressShear),
                YieldStressShearMinus(YieldStressShearMinus),
                YieldStressTension(YieldStressTension),
                FractureEnergyTension(FractureEnergyTension),
                ElasticLimitStressCompression(ElasticLimitStressCompression),
                YieldStressCompression(YieldStressCompression),
                YieldStrainCompression(YieldStrainCompression),
                ResidualStressCompression(ResidualStressCompression),
                FractureEnergyCompression(FractureEnergyCompression),
                BezierControllerC1(BezierControllerC1),
                BezierControllerC2(BezierControllerC2),
                BezierControllerC3(BezierControllerC3)
            {
            }

            DirectionalMaterialProperties(
                const DirectionalMaterialProperties& rOther)
                : CharacteristicLength(rOther.CharacteristicLength),
                E(rOther.E), nu(rOther.nu), G(rOther.G),
                YieldStressShear(rOther.YieldStressShear),
                YieldStressShearMinus(rOther.YieldStressShearMinus),
                YieldStressTension(rOther.YieldStressTension),
                FractureEnergyTension(rOther.FractureEnergyTension),
                ElasticLimitStressCompression(rOther.ElasticLimitStressCompression),
                YieldStressCompression(rOther.YieldStressCompression),
                YieldStrainCompression(rOther.YieldStrainCompression),
                ResidualStressCompression(rOther.ResidualStressCompression),
                FractureEnergyCompression(rOther.FractureEnergyCompression),
                BezierControllerC1(rOther.BezierControllerC1),
                BezierControllerC2(rOther.BezierControllerC2),
                BezierControllerC3(rOther.BezierControllerC3)
            {
            }
        };

        struct CalculationData {
            DirectionalMaterialProperties MaterialProperties1;
            DirectionalMaterialProperties MaterialProperties2;

            double BiaxialCompressionMultiplier;
            double ShearCompressionReductor;

            // Effective Stress Data
            array_1d<double, 3> EffectiveStressVectorOrthotropic;
            array_1d<double, 3> EffectiveStressVector;
            array_1d<double, 2> PrincipalStressVector;
            array_1d<double, 3> EffectiveTensionStressVector;
            array_1d<double, 3> EffectiveCompressionStressVector;

            Matrix ProjectionTensorTension;
            Matrix ProjectionTensorCompression;

            // Misc
            double DeltaTime;
            int TensionYieldModel;

            CalculationData(const DirectionalMaterialProperties& rMaterialProperties1,
                const DirectionalMaterialProperties& rMaterialProperties2)
                : MaterialProperties1(rMaterialProperties1),
                MaterialProperties2(rMaterialProperties2)
            {
                ProjectionTensorTension = ZeroMatrix(3, 3);
                ProjectionTensorCompression = ZeroMatrix(3, 3);
            }
        };

        struct TransformationMatrices {
            // Compute transformation from real anisotropic stress space to mapped isotropic stress space
            Matrix TransformationMatrixTension = ZeroMatrix(3, 3);
            Matrix TransformationMatrixCompression = ZeroMatrix(3, 3);
            // Compute transformation from mapped isotropic stress space to real anisotropic stress space
            Matrix InverseTransformationMatrixTension = ZeroMatrix(3, 3);
            Matrix InverseTransformationMatrixCompression = ZeroMatrix(3, 3);
        };

        /// Has variable <double>
        bool Has(const Variable<double> & rThisVariable) override;

        /// GetValue variable <double>
        double& GetValue(
            const Variable<double>& rThisVariable,
            double& rValue) override;

        /// SetValue variable <double>
        void SetValue(
            const Variable<double> & rVariable,
            const double& rValue,
            const ProcessInfo & rCurrentProcessInfo) override;

        /**
        * Is called to check whether the provided material parameters in the Properties
        * match the requirements of current constitutive model.
        * @param rMaterialProperties the current Properties to be validated against.
        * @return true, if parameters are correct; false, if parameters are insufficient / faulty
        * NOTE: this has to be implemented by each constitutive model. Returns false in base class since
        * no valid implementation is contained here.
        */
        bool ValidateInput(
            const Properties & rProperties) override;

        /**
        * returns the expected strain measure of this constitutive law (by default linear strains)
        * @return the expected strain measure
        */
        StrainMeasure GetStrainMeasure() override {
            return ConstitutiveLaw::StrainMeasure_Infinitesimal;
        }

        /**
        * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
        * @return the expected stress measure
        */
        StressMeasure GetStressMeasure() override {
            return ConstitutiveLaw::StressMeasure_Cauchy;
        };

        /**
        * returns whether this constitutive model is formulated in incremental strains/stresses
        * NOTE: by default, all constitutive models should be formulated in total strains
        * @return true, if formulated in incremental strains/stresses, false otherwise
        */
        bool IsIncremental() override {
            return false;
        }

        /// To check compatibility with element
        void GetLawFeatures(Features& rFeatures) override;


        /**
        * This is to be called at the very beginning of the calculation
        * (e.g. from InitializeElement) in order to initialize all relevant
        * attributes of the constitutive law
        * @param rMaterialProperties the Properties instance of the current element
        * @param rElementGeometry the geometry of the current element
        * @param rShapeFunctionsValues the shape functions values in the current integration point
        */
        void InitializeMaterial(
            const Properties & rProperties,
            const GeometryType & rGeometry,
            const Vector & rShapeFunctionsValues) override;

        /// Resets all internal variables as tresholds and damage parameters.
        void ResetMaterial(
            const Properties& rMaterialProperties,
            const GeometryType& rElementGeometry,
            const Vector& rShapeFunctionsValues) override;

        /**
        * to be called at the end of each solution step
        * (e.g. from Element::FinalizeSolutionStep)
        * @param rMaterialProperties the Properties instance of the current element
        * @param rElementGeometry the geometry of the current element
        * @param rShapeFunctionsValues the shape functions values in the current integration point
        * @param the current ProcessInfo instance
        */
        void FinalizeSolutionStep(
            const Properties & rProperties,
            const GeometryType & rGeometry,
            const Vector & rShapeFunctionsValues,
            const ProcessInfo & rCurrentProcessInfo) override;

        /// Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
        void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters & rValues) override;

        /// material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
        void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters & rValues) override;

        /// material response in terms of Kirchhoff stresses and constitutive tensor
        void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters & rValues) override;

        /// material response in terms of Cauchy stresses and constitutive tensor
        void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters & rValues) override;

        /**
        * This function is designed to be called once to perform all the checks needed
        * on the input provided. Checks can be "expensive" as the function is designed
        * to catch user's errors.
        * @param rMaterialProperties
        * @param rElementGeometry
        * @param rCurrentProcessInfo
        * @return
        */
        int Check(
            const Properties & rMaterialProperties,
            const GeometryType & rElementGeometry,
            const ProcessInfo & rCurrentProcessInfo) override;

    protected:

        ///@name Protected member Variables
        ///@{

        // Initialization
        bool InitializeDamageLaw = false;

        // Tension & Compression Thresholds
        double mCurrentThresholdTension = 0.0;
        double mCurrentThresholdCompression = 0.0; // at time step n:
        double mThresholdTension = 0.0;
        double mThresholdCompression = 0.0; // at time step n + 1:

        // Damage Parameters & Uniaxial Stresses
        double mDamageTension = 0.0;
        double mDamageCompression = 0.0;
        double mUniaxialStressTension = 0.0;
        double mUniaxialStressCompression = 0.0;

        // Misc
        double mInitialCharacteristicLength = 0.0;

        // for IMPLEX_Integration:
        double PreviousThresholdTension = 0.0;
        double PreviousThresholdCompression = 0.0;    // at time step n - 1:
        // end for IMPLEX_Integration

        // for IMPLEX Integration
        double CurrentDeltaTime = 0.0;    // at time step n + 1
        double PreviousDeltaTime = 0.0;   // at time step n

        // temporary internal variables saved for the implicit step in the finalize solution step of the IMPLEX 
        double TemporaryImplicitThresholdTension = 0.0;
        double TemporaryImplicitThresholdTCompression = 0.0;
        // end for IMPLEX Integration

        ///@}


        ///@name Protected Operators
        ///@{

        /**
        * @brief Initializes the CalculationData at the beginning of each SolutionStep
        * @param Properties& props         The Material Properties of the constitutive law from rValues
        *          GeometryType& geom     The Element Geometry from rValues
        *         ProcessInfo& pinfo        The ProcessInfo from rValues
        *         CalculationData
        */
        CalculationData GetCalculationData(
            const Properties& rProperties,
            const GeometryType& rGeometry,
            const ProcessInfo& rProcessInfo);

        /**
        * @brief Constructs the Linear Elasticity Matrix and stores it in the CalculationData
        * @param CalculationData
        */
        void CalculateOrthotropicElasticityMatrix(
            Matrix& rC,
            const DirectionalMaterialProperties& rMaterialProperties1,
            const DirectionalMaterialProperties& rMaterialProperties2);

        /**
        * @brief Assanbles the Transformation Matrix for Stress Mapping from Orthotropic to isotropic states and stores it in the CalculationData
        * @param CalculationData
        */
        void AssembleTransformationMatrix(
            const DirectionalMaterialProperties& rMaterialProperties1,
            const DirectionalMaterialProperties& rMaterialProperties2,
            TransformationMatrices& rTransformationMatrices);

        void CalculateProjectedFractureEnergyTension(
            const DirectionalMaterialProperties& rMaterialProperties1,
            const DirectionalMaterialProperties& rMaterialProperties2,
            double AngleToDamage,
            DirectionalMaterialProperties& rProjectedProperties);

        /**
        * @brief Splits the Effective Stress Vector into a positive (tension) and negative (compression) part
        * @param CalculationData
        */
        void TensionCompressionSplit(
            const array_1d<double, 3>& rEffectiveStressVector,
            array_1d<double, 2>& rPrincipalStressVector,
            array_1d<double, 3>& rEffectiveTensionStressVector,
            array_1d<double, 3>& rEffectiveCompressionStressVector);

        /**
        * @brief This method computes the positive (tension) and negative (compression) parts of the ProjectionMatrix
        * @param CalculationData
        */
        void ConstructProjectionTensors(
            const array_1d<double, 3>& rEffectiveStressVector,
            Matrix& rProjectionTensorTension,
            Matrix& rProjectionTensorCompression);

        /**
        * @brief This method computes the equivalent stress in Tension
        * @param CalculationData
        *          UniaxialStressTension The variable to be filled with the resulting value
        */
        void CalculateEquivalentStressTension(
            const CalculationData& data,
            const DirectionalMaterialProperties& rMaterialProperties,
            const array_1d<double, 3> rEffectiveStressVector,
            const array_1d<double, 2> rPrincipalStressVector,
            double& UniaxialStressTension) const;

        void CalculateProjectedFractureEnergyCompression(
            const DirectionalMaterialProperties& rMaterialProperties1,
            const DirectionalMaterialProperties& rMaterialProperties2,
            double AngleToDamage,
            DirectionalMaterialProperties& rProjectedProperties);

        /**
        * @brief This method computes the equivalent stress in Compression
        * @param CalculationData
        *          UniaxialStressCompression The variable to be filled with the resulting value
        */
        void CalculateEquivalentStressCompression(
            const CalculationData& data,
            const DirectionalMaterialProperties& rMaterialProperties,
            const array_1d<double, 3> rEffectiveStressVector,
            const array_1d<double, 2> rPrincipalStressVector,
            double& UniaxialStressCompression) const;

        /**
        * @brief This method computes the damage variable d+ of the tension law
        *        by considering an exponential softening behavior
        * @param CalculationData
        *        internal_variable   The internal variable that considers the considers the state of damage
        *        rDamage             The final damage variable to be filled by this method
        */
        void CalculateDamageTension(
            const DirectionalMaterialProperties& rMaterialProperties,
            const double TresholdTension,
            double& rDamageTension) const;

        /**
        *  BRIEF DOCUMENTATION OF THE USED UNIAXIAL SOFTENING BEHAVIOR IN COMPRESSION
        *  Entire documentation can be found in the the Phd Thesis of Massimo Petracca
        *  << Computational Multiscale Analysis of Masonry Structures>>
        *
        *  UNIAXIAL BEZIER COMPRESSION DAMAGE
        *  {I}   Linear Elastic
        *  {II}  Hardening Quadratic Bezier Curve
        *          Control nodes:  0=(e_0,s_0); I=(e_i,s_p); P=(e_p,s_p)
        *  {III} Softening Quadratic Bezier Curve
        *          Control nodes:  P=(e_p,s_p); J=(e_j,s_j); K=(e_k,s_k)
        *  {IV}  Softening Quadratic Bezier Curve
        *          Control nodes:  K=(e_k,s_k); R=(e_r,s_r); U=(e_u,s_u)
        *  {V}   Residual Strength
        *
        *    STRESS
        *       ^
        *      /|\
        *       |                     (P)
        * s_p = |------------(I)+----#####--+(J)
        * s_i = |               ' ###  ' ####
        * s_j   |              ###     '    ###
        *       |            ###'      '    ' ##
        * s_k   |-----------##--+------+----+--## (K)
        * s_0   |---------##(0) '      '    '   ###
        *       |        ## '   '      '    '    '###
        *       |       ##  '   '      '    '    '   ####
        *       |      ##   '   '      '    '    '      #####
        *       |     ##    '   '      '    '    '          #####
        *       |    ##     '   '      '    '    '    (R)       ######## (U)
        * s_r = |---##------+---+------'----+----+-----+-----------------######################
        * s_u   |  ##       '   '      '    '    '     '                 '
        *       |_##________+___+______+____+____+_____+_________________+______________________________\
        *                  e_0 e_i    e_p  e_j  e_k   e_r               e_u                             / STRAIN
        *        '          '          '         '                       '
        *        '   {I}    '   {II}   '  {III}  '        {IV}           '          {V}
        *        '          '          '         '                       '
        *
        */

        /**
        * @brief This method computes the damage variable d- of the compression law
        *        by considering the above explained bezier curved uniaxial behavior
        * @param CalculationData     Calculation Data for the CL
        *          internal_variable    The internal variable that considers the considers the state of damage
        *        rDamage             The final damage variable to be filled by this method
        */
        void CalculateDamageCompression(
            const DirectionalMaterialProperties& rMaterialProperties,
            const double TresholdCompression,
            double& rDamage) const;

        /**
        * 7@brief This method computes the energy of the uniaxial damage law before regularization
        * 7@param rBezierEnergy
        *         rBezierEnergy1
        *          double s_p, double s_k, double s_r, double e_p, double e_j, double e_k, double e_r, double e_u
        *                        As inputs for the energy calculation
        */

        void ComputeBezierEnergy(
            double& rBezierEnergy,
            double& rBezierEnergy1,
            double s_p, double s_k, double s_r,
            double e_p, double e_j, double e_k,
            double e_r, double e_u) const;


        /**
        * @brief This method evaluates the area below the bezier curves
        * @param x1, x2, x3 x-coordinates of the Bezier control points
        *          y1, y2, y3 y-coordinates of the Bezier control points
        */

        inline double EvaluateBezierArea(
            double x1, double x2, double x3,
            double y1, double y2, double y3) const;

        /**
        * @brief This method applies the stretcher to the strains, to regularize the fracture energy
        * @param stretcher            The stretch factor
        *          e_p                The reference strain from which on the stretching is applied
        *         e_j, e_k, e_r, e_u    The strains that have to be modified by the stretcher
        */

        inline void ApplyBezierStretcherToStrains(
            double stretcher, double e_p, double& e_j,
            double& e_k, double& e_r, double& e_u) const;

        /**
        * @brief This method evaluates the damage parameter by considering the Bezier law explained above
        * @param rDamageParameter    The parameter to obtain the damage d-
        *          xi                    The strain like counterpart
        *         x1, x2, x3         x-coordinates of the Bezier control points
        *          y1, y2, y3         y-coordinates of the Bezier control points
        */
        inline void EvaluateBezierCurve(
            double& rDamageParameter, double xi,
            double x1, double x2, double x3,
            double y1, double y2, double y3) const;

        /**
        * @brief This method computes the Charcteristic element length
        * @param GeometryType& geom        The element geometry data from rValues
        *           rCharacteristicLength    The characteristic Length
        */
        double ComputeCharacteristicLength(
            const GeometryType & geom,
            int DirectionIndex);


        /* Calculates the angle between the local coordinate and
        *  the principal stress, which is considered as the damage angle.
        */
        double CalculateDamageAngle(
            const array_1d<double, 3>& rStressVector);

        /**
        * @brief This method computes the internal material response strain to stress by applying cl
        * @param strain_vector        The strain vector
        *           stress_vector        The stress vector
        *         CalculationData    Calculation Data for the CL
        */
        void CalculateMaterialResponseInternal(
            Vector& rPredictiveStressVector,
            CalculationData& data,
            int IntegrationImplex);

        /**
         * @brief This method checks whether we are in loading/unloading/damage state
         * @param is_damaging_tension    The damage tension bool
         *          is_damaging_tension    The damage compression bool
         */
        void CheckDamageLoadingUnloading(
            bool& is_damaging_tension,
            bool& is_damaging_compression);

        /**
        * @brief This method computes the tangent tensor
        * @param rValues                 The constitutive law parameters and flags
        *        rStrainVector           The strain vector
        *        PredictiveStressVector  The stress vector
        *        CalculationData         Calculation Data for the CL
        *        IntegrationImplex       Using Integration implex if 1, elso not.
        */
        void CalculateTangentTensor(
            Parameters& rValues,
            const Matrix& rElasticityMatrix,
            const array_1d<double, 3>& rStrainVector,
            const array_1d<double, 3>& PredictiveStressVector,
            CalculationData& data,
            int IntegrationImplex);

        /**
         * @brief This method computes the secant tensor
         * @param rValues            The constitutive law parameters and flags
         *        CalculationData    Calculation Data for the CL
         */
        void CalculateSecantTensor(
            Parameters& rValues,
            const Matrix& rElasticityMatrix,
            const CalculationData & data);

        ///@}

    private:

        ///@name Serialization
        ///@{

        friend class Serializer;

        void save(Serializer & rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
        }

        void load(Serializer & rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
        }

        ///@}

    }; // Class MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw 

} // namespace Kratos
#endif // KRATOS_MASONRY_ORTHOTROPIC_DAMAGE_D_PLUS_D_MINUS_PLANE_STRESS_2D_H_INCLUDED  defined 
