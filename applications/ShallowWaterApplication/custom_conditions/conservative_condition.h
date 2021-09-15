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

#ifndef KRATOS_CONSERVATIVE_CONDITION_H_INCLUDED
#define KRATOS_CONSERVATIVE_CONDITION_H_INCLUDED

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
 * @class ConservativeCondition
 * @brief Implementation of a condition for shallow water waves problems
 * @author Miguel Maso Sotomayor
 */
template<std::size_t TNumNodes>
class ConservativeCondition : public WaveCondition<TNumNodes>
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

    typedef array_1d<double, 3*TNumNodes> LocalVectorType;

    typedef BoundedMatrix<double, 3*TNumNodes, 3*TNumNodes> LocalMatrixType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConservativeCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    ConservativeCondition() : WaveConditionType(){}

    /**
     * @brief Constructor using an array of nodes
     */
    ConservativeCondition(IndexType NewId, const NodesArrayType& ThisNodes) : WaveConditionType(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    ConservativeCondition(IndexType NewId, GeometryType::Pointer pGeometry) : WaveConditionType(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    ConservativeCondition(IndexType NewId, GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties) : WaveConditionType(NewId, pGeometry, pProperties){}

    /**
     * @brief Destructor
     */
    virtual ~ ConservativeCondition(){}

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
        return Kratos::make_intrusive<ConservativeCondition<TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
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
        return Kratos::make_intrusive<ConservativeCondition<TNumNodes>>(NewId, pGeom, pProperties);
    }

    /**
     * @brief Create a new condition pointer and clone the previous condition data
     * @param NewId the ID of the new condition
     * @param rThisNodes the nodes of the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        Condition::Pointer p_new_elem = Create(NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
        p_new_elem->SetData(this->GetData());
        p_new_elem->Set(Flags(*this));
        return p_new_elem;
    }

    /**
     * @brief Check that all required data containers are properly initialized and registered in Kratos
     * @return 0 if no errors are detected.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief Fill given vector with the linear system row index for the condition's degrees of freedom
     * @param rResult
     * @param rCurrentProcessInfo
     */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief Fill given array with containing the condition's degrees of freedom
     * @param rConditionalDofList
     * @param rCurrentProcessInfo
     */
    void GetDofList(DofsVectorType& rConditionalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief Get the variable which defines the degrees of freedom
     */
    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        return "ConservativeCondition";
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

    void AddWaveTerms(
        LocalMatrixType& rMatrix,
        LocalVectorType& rVector,
        const ConditionData& rData,
        const array_1d<double,TNumNodes>& rN,
        const double Weight = 1.0) override;

    void AddFluxTerms(
        LocalVectorType& rVector,
        const ConditionData& rData,
        const array_1d<double,TNumNodes>& rN,
        const double Weight = 1.0) override;

    void AddMassTerms(
        LocalMatrixType& rMatrix,
        const ConditionData& rData,
        const array_1d<double,TNumNodes>& rN,
        const double Weight) override {}

    double StabilizationParameter(const ConditionData& rData) const override;

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

}; // Class ConservativeCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONSERVATIVE_CONDITION_H_INCLUDED  defined
