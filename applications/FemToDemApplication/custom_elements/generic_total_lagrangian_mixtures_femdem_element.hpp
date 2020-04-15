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

#if !defined(KRATOS_GENERIC_TOTAL_LAGRANGIAN_MIXTURES_FEMDEM_ELEMENT_H_INCLUDED)
#define KRATOS_GENERIC_TOTAL_LAGRANGIAN_MIXTURES_FEMDEM_ELEMENT_H_INCLUDED

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
 * @brief Small Displacement element for the 2D and 3D cases
 * @author Alejandro Cornejo
 */
template<unsigned int TDim, unsigned int TyieldSurf>
class GenericTotalLagrangianMixturesFemDemElement 
    : public GenericTotalLagrangianFemDemElement<TDim,TyieldSurf> // Derived Element from SolidMechanics
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

    ///definition of node type (default is: Node<3>)
    typedef Node < 3 > NodeType;

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
                                             const ConstitutiveVariables &rThisConstVars, const KinematicVariables &rKinVariables, 
                                             Vector &rStrainVector, double& rDamageElement,  bool& rIsDamaging, const double CharacteristicLength,
                                             const bool SaveIntVars);
        
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

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif