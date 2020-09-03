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

#if !defined(KRATOS_KINEMATIC_SIMULATION_ELEMENT_H_INCLUDED)
#define KRATOS_KINEMATIC_SIMULATION_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
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
class KinematicSimulationElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Element;

    /// Node type (default is: Node<3>)
    using NodeType = BaseType::NodeType;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    using PropertiesType = Properties;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = GeometryType::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = BaseType::VectorType;

    /// Matrix type for local contributions to the linear system
    using MatrixType = BaseType::MatrixType;

    using EquationIdVectorType = BaseType::EquationIdVectorType;

    using DofsVectorType = BaseType::DofsVectorType;

    using IndexType = std::size_t;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of LaplaceElement
    KRATOS_CLASS_POINTER_DEFINITION(KinematicSimulationElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit KinematicSimulationElement(IndexType NewId = 0) : Element(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    KinematicSimulationElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    KinematicSimulationElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    KinematicSimulationElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    KinematicSimulationElement(KinematicSimulationElement const& rOther) : Element(rOther)
    {
    }

    /**
     * Destructor
     */
    ~KinematicSimulationElement() override = default;

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
        KRATOS_ERROR << "Attempting to Create base "
                        "Kinematic Simulation instances."
                     << std::endl;
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
        KRATOS_ERROR << "Attempting to Create base "
                        "KinematicSimulationElement instances."
                     << std::endl;
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
        KRATOS_ERROR << "Attempting to Clone base "
                        "KinematicSimulationElement instances."
                     << std::endl;
        KRATOS_CATCH("");
    }

    virtual const Variable<double>& GetVariable() const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "KinematicSimulationElement "
                        "GetVariable method. Please implement it in the "
                        "derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const override
    {
        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        const Variable<double>& r_variable = this->GetVariable();

        for (unsigned int i = 0; i < TNumNodes; ++i)
            rResult[i] = Element::GetGeometry()[i].GetDof(r_variable).EquationId();
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const override
    {
        if (rElementalDofList.size() != TNumNodes)
            rElementalDofList.resize(TNumNodes);

        const Variable<double>& r_variable = this->GetVariable();

        for (unsigned int i = 0; i < TNumNodes; ++i)
            rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(r_variable);
    }

    void GetValuesVector(VectorType& rValues, int Step) override
    {
        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);

        const Variable<double>& r_variable = this->GetVariable();

        GeometryType& rGeom = this->GetGeometry();
        IndexType LocalIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        {
            rValues[LocalIndex++] =
                rGeom[iNode].FastGetSolutionStepValue(r_variable, Step);
        }
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        return GeometryData::GI_GAUSS_1;
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {

        // Calculate components of turbulent kinetic energy

        // Calculate RHS
        this->CalculateRightHandSideContribution(rRightHandSideVector);

        // Calculate LHS
        this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

        VectorType values;
        this->GetValuesVector(values, 0);

        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, values);
        // KRATOS_WATCH(rRightHandSideVector);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        this->CalculateRightHandSideContribution(rRightHandSideVector);

        Matrix lhs;
        this->CalculateLeftHandSide(lhs, rCurrentProcessInfo);

        Vector values;
        this->GetValuesVector(values, 0);

        noalias(rRightHandSideVector) -= prod(lhs, values);
        // KRATOS_WATCH(rRightHandSideVector)
        // KRATOS_WATCH(lhs)
        // KRATOS_WATCH(values)

        KRATOS_CATCH("");
    }


    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side matrix only
     * @param rRHS the right hans side vector
     */
    virtual void CalculateRightHandSideContribution(Vector& rRHS) const
    {
        if (rRHS.size() != TNumNodes)
        {
            rRHS.resize(TNumNodes, false);
        }
        noalias(rRHS) = ZeroVector(TNumNodes);

        // Get Shape function data
        // vector of gauss weights for each gauss point
        Vector gauss_weights;
        // matrix of gauss shape functions evaluated at each gauss point
        // rows of this matrix contains shape function values evaluated at that respective gauss point
        Matrix shape_functions;
        // vector of matrices of shape function derivatives, a matrix for each shape function evaluated at gauss point
        // matrix row index corresponds to shape function, column index to dimension
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        //const GeometryType& r_geometry = this->GetGeometry();

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Vector& gauss_shape_functions = row(shape_functions, g);
            noalias(rRHS) += gauss_shape_functions * gauss_weights[g];
        }
    }

    /**
     * this is called during the CalculateLeftHandSide in order
     * to calculate integration w.r.t. wave number
     * @param TotalWaveNumber the total number of discretization
     * @param TurbulentEnergyDissipationRate
     * @param TurbulentKineticEnergy
     * @param EffectiveWaveNumber
     */
    virtual double CalculateWaveNumberIntegration(
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
                                        const double TurbulentKineticEnergyW) const
    {
        double output = 0.0;
        double k1 = 2*M_PI*TurbulentEnergyDissipationRate*std::pow(TurbulentKineticEnergy, -1.5);
        double kN = std::pow(TurbulentEnergyDissipationRate, 0.25)*std::pow(KinematicViscosity, -0.75);
        double dk = (std::log(kN) - std::log(k1))/(TotalWaveNumber-1);
        double kn = k1;
        double kn_pre = 0.0;

        //KRATOS_WATCH("constant v,w");
        //KRATOS_WATCH(TurbulentKineticEnergy);
        //KRATOS_WATCH(TurbulentEnergyDissipationRate);
        //KRATOS_WATCH(EffectiveWaveNumber);
        //KRATOS_WATCH(TotalWaveNumber);

        for (int i = 0; i < TotalWaveNumber; ++i)
        {
            output += (2/EffectiveWaveNumber) * (kn-kn_pre) * std::pow((kn/EffectiveWaveNumber),2) * std::exp(-2*std::pow((kn/kN), 2)) / (std::pow(1+std::pow((kn/EffectiveWaveNumber), 2), 1.833333333));
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
    virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        const int total_wave_number = rCurrentProcessInfo[TOTAL_WAVE_NUMBER];

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        const GeometryType& r_geometry = this->GetGeometry();

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Vector gauss_shape_functions = row(shape_functions, g);

            const double turbulent_kinetic_energy = RansCalculationUtilities::EvaluateInPoint(r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
            const double turbulent_energy_dissipation_rate = RansCalculationUtilities::EvaluateInPoint(r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
            const double effective_wave_number = RansCalculationUtilities::EvaluateInPoint(r_geometry, EFFECTIVE_WAVE_NUMBER, gauss_shape_functions);
            const double kinematic_viscosity = RansCalculationUtilities::EvaluateInPoint(r_geometry, KINEMATIC_VISCOSITY, gauss_shape_functions);
            const double r_wavenumber_integration= this-> CalculateWaveNumberIntegration(total_wave_number,turbulent_energy_dissipation_rate,
                                                  turbulent_kinetic_energy,effective_wave_number, kinematic_viscosity,
                                                  0, 0, 0,
                                                  0, 0, 0);

            //KRATOS_WATCH("CalculateLHS_baseelement")
            //KRATOS_WATCH(turbulent_kinetic_energy)
            //KRATOS_WATCH(turbulent_energy_dissipation_rate)
            //KRATOS_WATCH(effective_wave_number)
            //KRATOS_WATCH(kinematic_viscosity)
            //KRATOS_WATCH(r_wavenumber_integration)

            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                for (IndexType b = 0; b < TNumNodes; ++b)
                {
                    rLeftHandSideMatrix(a, b) += gauss_weights[g] * gauss_shape_functions[a] * gauss_shape_functions[b] * r_wavenumber_integration;
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

        const Variable<double>& r_variable = this->GetVariable();

        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            const NodeType& r_node = this->GetGeometry()[iNode];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(EFFECTIVE_WAVE_NUMBER, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_U, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_V, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_W, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_variable, r_node);
            KRATOS_CHECK_DOF_IN_NODE(r_variable, r_node);
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
        buffer << "KinematicSimulationElement #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KinematicSimulationElement #" << this->Id();
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
    /**
     * @brief Calculates shape function data for this element
     *
     * @param rGaussWeights Gauss point weights list
     * @param rNContainer   Shape function values. Each row contains shape functions for respective gauss point
     * @param rDN_DX        List of matrices containing shape function derivatives for each gauss point
     */
    void CalculateGeometryData(Vector& rGaussWeights,
                                       Matrix& rNContainer,
                                       ShapeFunctionDerivativesArrayType& rDN_DX) const
    {
        const GeometryType& r_geometry = this->GetGeometry();

        RansCalculationUtilities::CalculateGeometryData(
            r_geometry, this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
    }

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
                                KinematicSimulationElement<TNumNodes>& rThis);

/// output stream function
template <unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const KinematicSimulationElement<TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_LAPLACE_ELEMENT_H_INCLUDED defined
