//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#pragma once

// System includes

// External include

// Project includes

#include "custom_elements/generic_total_lagrangian_femdem_element.h"
#include "custom_utilities/constitutive_law_utilities.h"

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
 * @class GenericTotalLagrangianMixturesFemDemElement
 * @ingroup FemToDemApplication
 * @brief Total Lagrangian element taking into account a Classical Rule Of Mixtures
 * and j2-plasticity for the fiber. See Computational Methods for Plasticity - EA de Souza Neto, D Peric and DRJ Owen
 * @author Alejandro Cornejo
 */
template<unsigned int TDim, unsigned int TyieldSurf>
class KRATOS_API(FEM_TO_DEM_APPLICATION) GenericTotalLagrangianMixturesFemDemElement
    : public GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    ///definition of element type
    typedef Element ElementType;

    ///base type: an GeometricalObject that automatically has a unique number
    typedef GenericTotalLagrangianFemDemElement<TDim,TyieldSurf> BaseType;

    ///definition of node type (default is: Node)
    typedef Node NodeType;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    ///definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef GeometryData GeometryDataType;

    /// We define the dimension
    static constexpr SizeType VoigtSize = (TDim == 3) ? 6 : 3;

    /// We define the number of edges
    static constexpr SizeType NumberOfEdges = (TDim == 3) ? 6 : 3;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    /// The zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// The definition of the bounded vector type
    typedef array_1d<double, VoigtSize> BoundedVectorType;

    /// Counted pointer of GenericTotalLagrangianMixturesFemDemElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GenericTotalLagrangianMixturesFemDemElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    GenericTotalLagrangianMixturesFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry);
    GenericTotalLagrangianMixturesFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    // GenericTotalLagrangianMixturesFemDemElement(GenericTotalLagrangianMixturesFemDemElement const &rOther);

    /// Destructor.
    virtual ~GenericTotalLagrangianMixturesFemDemElement();

    /// Assignment operator.
    GenericTotalLagrangianMixturesFemDemElement(GenericTotalLagrangianMixturesFemDemElement const& rOther)
        : BaseType(rOther)
    {};

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    GenericTotalLagrangianMixturesFemDemElement()
    {
    }

protected:

    /**
     * this performs the smooting and integrates the CL and returns the integrated Stress
     */
    Vector IntegrateSmoothedConstitutiveLaw(const std::string &rYieldSurface, ConstitutiveLaw::Parameters &rValues,
                                            const BaseSolidElement::ConstitutiveVariables &rThisConstVars, const BaseSolidElement::KinematicVariables &rKinVariables,
                                            Vector &rStrainVector, double &rDamageElement, bool &rIsDamaging, const double CharacteristicLength,
                                            const bool SaveIntVars) override;

    /**
     * this integrates the stress according to plasticity
     */
    Vector IntegrateStressPlasticity(ConstitutiveLaw::Parameters &rValues,
                                     const Vector &rStrainVector,
                                     Vector& rPlasticStrainVector,
                                     double &rAcumulatedPlasticStrain,
                                     double &rThreshold,
                                     double& rUniaxialStress,
                                     bool &rIsPlastifying);
    /**
     * this method computes the plastic multiplier \dot{\lambda}
     * used for computing the plastic strain vector
     */
    void ComputePlasticMultiplier(const double UniaxialStress,
                                  const double Threshold,
                                  double &rPlasticMultiplier,
                                  ConstitutiveLaw::Parameters &rValues);
    /**
     * this method computes the plastic threshold  S = S0 + H*E^p
     */
    void ComputePlasticThreshold(const double AcumulatedPlasticStrain,
                                 double &rThreshold,
                                 ConstitutiveLaw::Parameters &rValues);

    /**
     * this integrates the perturbed strain
     */
    void IntegratePerturbedStrain(Vector &rPerturbedStressVector,
                                  const Vector &rPerturbedStrainVector,
                                  const Matrix &rElasticMatrix,
                                  ConstitutiveLaw::Parameters &rValues) override;

    void CalculateOnIntegrationPoints(
        const Variable<double> &rVariable,
        std::vector<double> &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    Vector CalculateAveragePlasticStrain()
    {
        Vector average_plastic_strain(VoigtSize);
        noalias(average_plastic_strain) = ZeroVector(VoigtSize);
        for (IndexType i = 0; i < NumberOfEdges; ++i)
            average_plastic_strain += mPlasticStrains[i];
        return average_plastic_strain / NumberOfEdges;
    }

    double CalculateAverageAcumulatedPlasticStrain()
    {
        double acumulated_strain = 0.0;
        for (IndexType i = 0; i < NumberOfEdges; ++i)
            acumulated_strain += mAcumulatedPlasticStrains[i];
        return acumulated_strain / NumberOfEdges;
    }


    void CheckIfEraseElement(
        const ProcessInfo &rCurrentProcessInfo,
        const Properties& rProperties
        ) override
    {
        if (this->mDamage >= 0.98 && this->CalculateAverageAcumulatedPlasticStrain() >= rProperties[MAX_PLASTIC_STRAIN]) {
            this->Set(ACTIVE, false);
            this->mDamage = 0.98;
            // We set a "flag" to generate the DEM
            // rCurrentProcessInfo[GENERATE_DEM] = true;
            this->SetValue(GENERATE_DEM, true);
        }
    }

    void CalculateAll(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag) override;

    Vector mAcumulatedPlasticStrains;
    Vector mPlasticityThresholds;
    std::vector<Vector> mPlasticStrains;
    double mUniaxialStress = 0.0;

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
}; // Class GenericTotalLagrangianMixturesFemDemElement

template<unsigned int TDim, unsigned int TyieldSurf> constexpr SizeType GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::VoigtSize;
template<unsigned int TDim, unsigned int TyieldSurf> constexpr SizeType GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::NumberOfEdges;

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.