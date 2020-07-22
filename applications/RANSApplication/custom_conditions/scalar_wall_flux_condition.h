//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_SCALAR_WALL_FLUX_CONDITION_H_INCLUDED)
#define KRATOS_RANS_SCALAR_WALL_FLUX_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/condition.h"
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
class ScalarWallFluxCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Condition;

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = Vector;

    /// Matrix type for local contributions to the linear system
    using MatrixType = Matrix;

    using IndexType = std::size_t;

    using EquationIdVectorType = std::vector<IndexType>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    using ScalarWallFluxConditionDataType = TScalarWallFluxConditionData;

    using CurrentConditionType =
        ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ScalarWallFluxCondition
    KRATOS_CLASS_POINTER_DEFINITION(ScalarWallFluxCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit ScalarWallFluxCondition(
        IndexType NewId = 0)
    : Condition(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    ScalarWallFluxCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : Condition(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    ScalarWallFluxCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    ScalarWallFluxCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    ScalarWallFluxCondition(
        ScalarWallFluxCondition const& rOther)
    : Condition(rOther)
    {
    }

    /**
     * Destructor
     */
    ~ScalarWallFluxCondition() override = default;

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
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentConditionType>(
            NewId, Condition::GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentConditionType>(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentConditionType>(
            NewId, Condition::GetGeometry().Create(ThisNodes), Condition::pGetProperties());
        KRATOS_CATCH("");
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo) const override
    {
        if (rResult.size() != TNumNodes) {
            rResult.resize(TNumNodes, false);
        }

        const Variable<double>& r_variable =
            TScalarWallFluxConditionData::GetScalarVariable();

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rResult[i] = Condition::GetGeometry()[i].GetDof(r_variable).EquationId();
        }
    }

    /**
     * determines the elemental list of DOFs
     * @param ConditionalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(
        DofsVectorType& rConditionalDofList,
        const ProcessInfo& CurrentProcessInfo) const override
    {
        if (rConditionalDofList.size() != TNumNodes) {
            rConditionalDofList.resize(TNumNodes);
        }

        const Variable<double>& r_variable =
            TScalarWallFluxConditionData::GetScalarVariable();

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rConditionalDofList[i] = Condition::GetGeometry()[i].pGetDof(r_variable);
        }
    }

    void GetValuesVector(
        Vector& rValues,
        int Step = 0) const override
    {
        this->GetFirstDerivativesVector(rValues, Step);
    }

    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) const override
    {
        if (rValues.size() != TNumNodes) {
            rValues.resize(TNumNodes, false);
        }

        const auto& r_geometry = this->GetGeometry();
        const Variable<double>& r_variable =
            TScalarWallFluxConditionData::GetScalarVariable();

        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            rValues[local_index++] =
                r_geometry[i_node].FastGetSolutionStepValue(r_variable, Step);
        }
    }

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) const override
    {
        if (rValues.size() != TNumNodes) {
            rValues.resize(TNumNodes, false);
        }

        const auto& r_geometry = this->GetGeometry();
        const Variable<double>& r_variable =
            TScalarWallFluxConditionData::GetScalarRateVariable();

        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            rValues[local_index++] =
                r_geometry[i_node].FastGetSolutionStepValue(r_variable, Step);
        }
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
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        // Check sizes and initialize matrix
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

        // Calculate RHS
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rRightHandSideVector.size() != TNumNodes) {
            rRightHandSideVector.resize(TNumNodes, false);
        }

        noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

        if (RansCalculationUtilities::IsWallFunctionActive(*this)) {
            const auto& r_geometry = this->GetGeometry();
            // Get Shape function data
            Vector gauss_weights;
            Matrix shape_functions;
            RansCalculationUtilities::CalculateConditionGeometryData(
                r_geometry, TScalarWallFluxConditionData::GetIntegrationMethod(),
                gauss_weights, shape_functions);
            const IndexType num_gauss_points = gauss_weights.size();

            TScalarWallFluxConditionData r_current_data(r_geometry);

            r_current_data.CalculateConstants(rCurrentProcessInfo);

            if (r_current_data.IsWallFluxComputable()) {
                for (IndexType g = 0; g < num_gauss_points; ++g) {
                    const Vector& gauss_shape_functions = row(shape_functions, g);

                    const double flux =
                        r_current_data.CalculateWallFlux(gauss_shape_functions);

                    noalias(rRightHandSideVector) +=
                        gauss_shape_functions * (gauss_weights[g] * flux);
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
        TScalarWallFluxConditionData::Check(this->GetGeometry(), rCurrentProcessInfo);

        return check;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ScalarWallFluxCondition #" << Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SWF" << TScalarWallFluxConditionData::GetName();
    }

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);

        KRATOS_CATCH("");
    }

    ///@}

}; // Class ScalarWallFluxCondition

///@}
///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
inline std::istream& operator>>(
    std::istream& rIStream,
    ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_SCALAR_WALL_FLUX_CONDITION_H_INCLUDED defined
