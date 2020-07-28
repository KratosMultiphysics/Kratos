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

#if !defined(KRATOS_KINEMATIC_SIMULATION_EFFECTIVE_WAVE_NUMBER_H_INCLUDED)
#define KRATOS_KINEMATIC_SIMULATION_EFFECTIVE_WAVE_NUMBER_H_INCLUDED

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
class KinematicSimulationEffectiveWaveNumberElement : public KinematicSimulationElement<TNumNodes>
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

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of LaplaceElement
    KRATOS_CLASS_POINTER_DEFINITION(KinematicSimulationEffectiveWaveNumberElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit KinematicSimulationEffectiveWaveNumberElement(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    KinematicSimulationEffectiveWaveNumberElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    KinematicSimulationEffectiveWaveNumberElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    KinematicSimulationEffectiveWaveNumberElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    KinematicSimulationEffectiveWaveNumberElement(KinematicSimulationEffectiveWaveNumberElement const& rOther)
        : BaseType(rOther)
    {
    }

    /**
     * Destructor
     */
    ~KinematicSimulationEffectiveWaveNumberElement() override = default;

    ///@}
    ///@name Operators
    ///@{

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
        return Kratos::make_intrusive<KinematicSimulationEffectiveWaveNumberElement<TNumNodes>>(
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
        return Kratos::make_intrusive<KinematicSimulationEffectiveWaveNumberElement<TNumNodes>>(
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
        return Kratos::make_intrusive<KinematicSimulationEffectiveWaveNumberElement<TNumNodes>>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    const Variable<double>& GetVariable() const override
    {
        return EFFECTIVE_WAVE_NUMBER;
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side matrix only
     * @param rRHS the right hans side vector
     */
    void CalculateRightHandSideContribution(Vector& rRHS) const override
    {
        if (rRHS.size() != TNumNodes)
        {
            rRHS.resize(TNumNodes, false);
        }
        noalias(rRHS) = ZeroVector(TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        const GeometryType& r_geometry = this->GetGeometry();

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {

            const Vector& gauss_shape_functions = row(shape_functions, g);
            const double turbulent_energy_dissipation_rate = RansCalculationUtilities::EvaluateInPoint(r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);

            noalias(rRHS) += gauss_shape_functions * (gauss_weights[g] * turbulent_energy_dissipation_rate);
        }
    }

    /**
     * this is called during the CalculateLeftHandSide in order
     * to calculate integration w.r.t. wave number
     * @param TotalWaveNumber the total number of discretization
     * @param TurbulentEnergyDissipationRate at the node
     * @param TurbulentKineticEnergy at the node
     * @param EffectiveWaveNumber at the node
     * @param KinematicViscosity at the node
     * @param SpectralConstantU at the node (Au)
     * @param SpectralConstantV at the node (Av)
     * @param SpectralConstantW at the node (Aw)
     * @param TurbulentKineticEnergyU at the node
     * @param TurbulentKineticEnergyV at the node
     * @param TurbulentKineticEnergyW at the node
     */
    double CalculateWaveNumberIntegration(
                                        const int TotalWaveNumber,
                                        const double TurbulentEnergyDissipationRate,
                                        const double TurbulentKineticEnergy,
                                        const double EffectiveWaveNumber,
                                        const double KinematicViscosity,
                                        const double SpectralConstantU,
                                        const double SpectralConstantV,
                                        const double SpectralConstantW,
                                        const double TurbulentKineticEnergyU,
                                        const double TurbulentKineticEnergyV,
                                        const double TurbulentKineticEnergyW
                                        ) const override
    {
        double output = 0.0;
        double k1 = 2*M_PI*TurbulentEnergyDissipationRate*std::pow(TurbulentKineticEnergy, -1.5);
        double kN = std::pow(TurbulentEnergyDissipationRate, 0.25)*std::pow(KinematicViscosity, -0.75);
        double dk = (std::log(kN) - std::log(k1))/(TotalWaveNumber-1);
        double kn = k1;
        double kn_pre = 0.0; 
        for (int i = 0; i < TotalWaveNumber; ++i)
        {
            output += (2*std::pow(kn,2)/(std::pow(EffectiveWaveNumber,2))) * (kn-kn_pre) * std::exp(-2*std::pow((kn/kN), 2))
             * (SpectralConstantU*TurbulentKineticEnergyU
             + (std::pow((kn/EffectiveWaveNumber),2)
             *(SpectralConstantV*TurbulentKineticEnergyV+SpectralConstantW*TurbulentKineticEnergyW) 
             / (1+std::pow((kn/EffectiveWaveNumber), 2))))
             / (std::pow(1+std::pow((kn/EffectiveWaveNumber), 2), 0.83333333333333));
            kn_pre = kn;
            kn = std::exp(std::log(k1)+(i+1)*dk);
        }

        return output;
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                              const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

        //const double total_wave_number = rCurrentProcessInfo[TOTAL_WAVE_NUMBER];
        const int total_wave_number = rCurrentProcessInfo[ACTIVATION_LEVEL];

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        const GeometryType& r_geometry = this->GetGeometry();

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Vector& gauss_shape_functions = row(shape_functions, g);
            //const Matrix& r_shape_values = shape_functions[g];

            const double turbulent_kinetic_energy = RansCalculationUtilities::EvaluateInPoint(r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
            const double turbulent_energy_dissipation_rate = RansCalculationUtilities::EvaluateInPoint(r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
            const double effective_wave_number = RansCalculationUtilities::EvaluateInPoint(r_geometry, EFFECTIVE_WAVE_NUMBER, gauss_shape_functions);
            const double spectral_constant_u = RansCalculationUtilities::EvaluateInPoint(r_geometry, SPECTRAL_CONSTANT_U, gauss_shape_functions);
            const double spectral_constant_v = RansCalculationUtilities::EvaluateInPoint(r_geometry, SPECTRAL_CONSTANT_V, gauss_shape_functions);
            const double spectral_constant_w = RansCalculationUtilities::EvaluateInPoint(r_geometry, SPECTRAL_CONSTANT_W, gauss_shape_functions);
            const double turbulent_kinetic_energy_u = RansCalculationUtilities::EvaluateInPoint(r_geometry, TURBULENT_KINETIC_ENERGY_U, gauss_shape_functions);
            const double turbulent_kinetic_energy_v = RansCalculationUtilities::EvaluateInPoint(r_geometry, TURBULENT_KINETIC_ENERGY_V, gauss_shape_functions);
            const double turbulent_kinetic_energy_w = RansCalculationUtilities::EvaluateInPoint(r_geometry, TURBULENT_KINETIC_ENERGY_W, gauss_shape_functions);
            const double kinematic_viscosity = RansCalculationUtilities::EvaluateInPoint(r_geometry, KINEMATIC_VISCOSITY, gauss_shape_functions);
            const double r_wavenumber_integration=    this-> CalculateWaveNumberIntegration(total_wave_number,turbulent_energy_dissipation_rate,
                                                  turbulent_kinetic_energy,effective_wave_number, kinematic_viscosity,
                                                  spectral_constant_u, spectral_constant_v, spectral_constant_w,
                                                  turbulent_kinetic_energy_u, turbulent_kinetic_energy_v, turbulent_kinetic_energy_w);

            KRATOS_WATCH(turbulent_kinetic_energy_u)
            KRATOS_WATCH(turbulent_kinetic_energy_v)
            KRATOS_WATCH(turbulent_kinetic_energy_w)
            KRATOS_WATCH(turbulent_energy_dissipation_rate)
            KRATOS_WATCH(spectral_constant_u)
            KRATOS_WATCH(spectral_constant_v)
            KRATOS_WATCH(spectral_constant_w)
            KRATOS_WATCH(r_wavenumber_integration)

            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                for (IndexType b = 0; b < TNumNodes; ++b)
                {
                    rLeftHandSideMatrix(a, b) += gauss_weights[g] * gauss_shape_functions[a]*gauss_shape_functions[b] * 2 * kinematic_viscosity * r_wavenumber_integration;
                }
            }
        }

        KRATOS_CATCH("");
    }



    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int check = BaseType::Check(rCurrentProcessInfo);

        //const Variable<double>& r_variable = this->GetVariable();

        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            const NodeType& r_node = this->GetGeometry()[iNode];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(EFFECTIVE_WAVE_NUMBER, r_node);
        }

        return check;

        KRATOS_CATCH("");
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "KinematicSimulationEffectiveWaveNumberElement #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KinematicSimulationEffectiveWaveNumberElement #" << this->Id();
    }

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
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class LaplaceElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

template <unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                KinematicSimulationEffectiveWaveNumberElement<TNumNodes>& rThis);

/// output stream function
template <unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const KinematicSimulationEffectiveWaveNumberElement<TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_LAPLACE_ELEMENT_H_INCLUDED defined
