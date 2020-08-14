//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mayu Sakuma
//

#if !defined(KRATOS_KINEMATIC_SIMULATION_SPECTRAL_CONSTANT_U_H_INCLUDED)
#define KRATOS_KINEMATIC_SIMULATION_SPECTRAL_CONSTANT_U_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/kinematic_simulation/kinematic_simulation_element.h"
#include "includes/checks.h"
#include "includes/define.h"

// Application includes
#include "includes/variables.h"
#include "rans_application_variables.h"

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

template <unsigned int TNumNodes>
class KinematicSimulationSpectralConstantUElement : public KinematicSimulationElement<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = KinematicSimulationElement<TNumNodes>;

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    using PropertiesType = Properties;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = GeometryType::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = Element::VectorType;

    /// Matrix type for local contributions to the linear system
    using MatrixType = Element::MatrixType;

    using IndexType = std::size_t;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of LaplaceElement
    KRATOS_CLASS_POINTER_DEFINITION(KinematicSimulationSpectralConstantUElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit KinematicSimulationSpectralConstantUElement(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    KinematicSimulationSpectralConstantUElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    KinematicSimulationSpectralConstantUElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    KinematicSimulationSpectralConstantUElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    KinematicSimulationSpectralConstantUElement(KinematicSimulationSpectralConstantUElement const& rOther)
        : BaseType(rOther)
    {
    }

    /**
     * Destructor
     */
    ~KinematicSimulationSpectralConstantUElement() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive<KinematicSimulationSpectralConstantUElement<TNumNodes>>(
            NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive<KinematicSimulationSpectralConstantUElement<TNumNodes>>(
            NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive<KinematicSimulationSpectralConstantUElement<TNumNodes>>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    const Variable<double>& GetVariable() const override
    {
        return SPECTRAL_CONSTANT_U;
    }


    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */

    /**
     * this is called during the CalculateLeftHandSide in order
     * to calculate integration w.r.t. wave number
     * @param TotalWaveNumber the total number of discretization
     * @param TurbulentEnergyDissipationRate
     * @param TurbulentKineticEnergy
     * @param EffectiveWaveNumber
     */
    double CalculateWaveNumberIntegration(const int TotalWaveNumber,
                                        const double TurbulentEnergyDissipationRate,
                                        const double TurbulentKineticEnergy,
                                        const double EffectiveWaveNumber,
                                        const double KinematicViscosity,
                                        const double SpectralConstantU,
                                        const double SpectralConstantV,
                                        const double SpectralConstantW,
                                        const double TurbulentKineticEnergyU,
                                        const double TurbulentKineticEnergyV,
                                        const double TurbulentKineticEnergyW) const override
    {
        double output = 0.0;
        double k1 = 2*M_PI*TurbulentEnergyDissipationRate*std::pow(TurbulentKineticEnergy, -1.5);
        double kN = std::pow(TurbulentEnergyDissipationRate, 0.25)*std::pow(KinematicViscosity, -0.75);
        double dk = (std::log(kN) - std::log(k1))/(TotalWaveNumber-1);
        double kn = k1;
        double kn_pre = 0.0; 
        for (int i = 0; i < TotalWaveNumber; ++i)
        {
            output += (2/EffectiveWaveNumber) * (kn-kn_pre) * std::exp(-2*std::pow((kn/kN), 2)) / (std::pow(1+std::pow((kn/EffectiveWaveNumber), 2), 0.8333333333));
            kn_pre = kn;
            kn = std::exp(std::log(k1)+(i+1)*dk);
        }
        return output;
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "KinematicSimulationSpectralConstantUElement #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KinematicSimulationSpectralConstantUElement #" << this->Id();
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    ///@}

}; // Class LaplaceElement

///@}

///@name Input and output
///@{

template <unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                KinematicSimulationSpectralConstantUElement<TNumNodes>& rThis);

/// output stream function
template <unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const KinematicSimulationSpectralConstantUElement<TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_LAPLACE_ELEMENT_H_INCLUDED defined
