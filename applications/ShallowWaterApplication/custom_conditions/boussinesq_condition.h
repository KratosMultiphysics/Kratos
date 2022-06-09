//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_BOUSSINESQ_CONDITION_H_INCLUDED
#define KRATOS_BOUSSINESQ_CONDITION_H_INCLUDED

// System includes


// External includes


// Project includes
#include "wave_condition.h"

namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @ingroup ShallowWaterApplication
 * @class BoussinesqCondition
 * @brief Implementation of a condition for shallow water waves problems
 * @author Miguel Maso Sotomayor
 */
template<std::size_t TNumNodes>
class BoussinesqCondition : public WaveCondition<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef WaveCondition<TNumNodes> WaveConditionType;

    typedef typename WaveConditionType::NodesArrayType NodesArrayType;

    typedef typename WaveConditionType::PropertiesType PropertiesType;

    typedef typename WaveConditionType::EquationIdVectorType EquationIdVectorType;

    typedef typename WaveConditionType::DofsVectorType DofsVectorType;

    typedef typename WaveConditionType::ConditionData ConditionData;

    typedef typename WaveConditionType::LocalVectorType LocalVectorType;

    typedef typename WaveConditionType::LocalMatrixType LocalMatrixType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BoussinesqCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    BoussinesqCondition() : WaveConditionType(){}

    /**
     * @brief Constructor using an array of nodes
     */
    BoussinesqCondition(IndexType NewId, const NodesArrayType& ThisNodes) : WaveConditionType(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    BoussinesqCondition(IndexType NewId, GeometryType::Pointer pGeometry) : WaveConditionType(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    BoussinesqCondition(IndexType NewId, GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties) : WaveConditionType(NewId, pGeometry, pProperties){}

    /**
     * @brief Destructor
     */
    ~ BoussinesqCondition() override {};

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<BoussinesqCondition<TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * @brief Create a new condition pointer
     * @param NewId: the ID of the new condition
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<BoussinesqCondition<TNumNodes>>(NewId, pGeom, pProperties);
    }

    /**
     * @brief Create a new condition pointer and clone the previous condition data
     * @param NewId the ID of the new condition
     * @param rThisNodes the nodes of the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        Condition::Pointer p_new_cond = Create(NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
        p_new_cond->SetData(this->GetData());
        p_new_cond->Set(Flags(*this));
        return p_new_cond;
    }

    /**
     * @brief Calculate the velocity laplacian projection
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containing the conditions
     */
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate the conditional contribution to the problem
     * @param rLeftHandSideMatrix Conditional left hand side matrix
     * @param rRightHandSideVector Conditional right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containing the condition
     */
    void CalculateLocalSystem(Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method provides the specifications/requirements of the element
     * @return specifications The required specifications/requirements
     */
    const Parameters GetSpecifications() const override;

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        return "BoussinesqCondition";
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    static constexpr IndexType mLocalSize = WaveConditionType::mLocalSize;

    ///@}
    ///@name Protected Operations
    ///@{

    const Variable<double>& GetUnknownComponent(int Index) const override;

    LocalVectorType GetUnknownVector(ConditionData& rData) override;

    void CalculateGaussPointData(
        ConditionData& rData,
        const IndexType PointIndex,
        const array_1d<double,TNumNodes>& rN) override;

    void AddLaplacianBoundary(
        LocalVectorType& rVelocityLaplacian,
        LocalVectorType& rMomentumLaplacian,
        const GeometryType& rParentGeometry,
        const ConditionData& rData,
        const array_1d<double,TNumNodes>& rN,
        const Matrix& rDN_DX,
        const double Weight = 1.0);

    void AddMomentumDispersionTerms(
        LocalVectorType& rLaplacianBoundary,
        const GeometryType& rParentGeometry,
        const ConditionData& rData,
        const array_1d<double,TNumNodes>& rN,
        const Matrix& rDN_DX,
        const double Weight = 1.0);

    void CalculateShapeFunctionDerivaties(
        Matrix& rDN_DX,
        const GeometryType& rParentGeometry,
        const Point& rPoint);

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class BoussinesqCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOUSSINESQ_CONDITION_H_INCLUDED  defined
