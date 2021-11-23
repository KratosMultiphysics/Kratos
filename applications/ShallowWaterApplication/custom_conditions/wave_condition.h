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

#ifndef KRATOS_WAVE_CONDITION_H_INCLUDED
#define KRATOS_WAVE_CONDITION_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"

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
 * @class WaveCondition
 * @brief Implementation of a condition for shallow water waves problems
 * @author Miguel Maso Sotomayor
 */
template<std::size_t TNumNodes>
class WaveCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    typedef Node<3> NodeType;

    typedef array_1d<double, 3*TNumNodes> LocalVectorType;

    typedef BoundedMatrix<double, 3*TNumNodes, 3*TNumNodes> LocalMatrixType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /// Pointer definition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(WaveCondition);

    ///@}
    ///@name Pointer definition
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    WaveCondition() : Condition(){}

    /**
     * @brief Constructor using an array of nodes
     */
    WaveCondition(IndexType NewId, const NodesArrayType& ThisNodes) : Condition(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    WaveCondition(IndexType NewId, GeometryType::Pointer pGeometry) : Condition(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    WaveCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Condition(NewId, pGeometry, pProperties){}

    /**
     * @brief Destructor
     */
    ~ WaveCondition() override {};

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
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<WaveCondition<TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * @brief Create a new condition pointer
     * @param NewId: the ID of the new condition
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<WaveCondition<TNumNodes>>(NewId, pGeom, pProperties);
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

    /**
     * @brief Get the time derivative of variable which defines the degrees of freedom
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    /**
     * @brief Get the second time derivative of variable which defines the degrees of freedom
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    /**
     * @brief Calculate the conditional contribution to the problem
     * @param rRightHandSideVector Conditional right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containing the condition
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate the conditional contribution to the problem
     * @param rLeftHandSideMatrix Conditional left hand side matrix
     * @param rRightHandSideVector Conditional right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containing the condition
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate the condition mass matrix
     * @param rMassMatrix the condition mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        return "WaveCondition";
    }

    /**
     * @brief Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << " : " << Id();
    }

    /**
     * @brief Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << GetGeometry();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    static constexpr IndexType mLocalSize = 3 * TNumNodes;

    ///@}
    ///@name Protected Classes
    ///@{

    struct ConditionData
    {
        bool integrate_by_parts;
        double gravity;
        double height;
        double length;
        double stab_factor;
        double relative_dry_height;
        array_1d<double,3> velocity;
        array_1d<double,3> normal;

        BoundedMatrix<double,3,3> A1;
        BoundedMatrix<double,3,3> A2;
        array_1d<double,3> b1;
        array_1d<double,3> b2;

        array_1d<double,TNumNodes> nodal_h;
        array_1d<double,TNumNodes> nodal_z;
        array_1d<array_1d<double,3>,TNumNodes> nodal_v;
        array_1d<array_1d<double,3>,TNumNodes> nodal_q;
    };
 
    ///@}
    ///@name Protected Operations
    ///@{

    virtual const Variable<double>& GetUnknownComponent(int Index) const;

    virtual LocalVectorType GetUnknownVector(ConditionData& rData);

    void CalculateGeometryData(
        Vector &rGaussWeights,
        Matrix &rNContainer) const;

    void InitializeData(
        ConditionData& rData,
        const ProcessInfo& rProcessInfo);

    virtual void CalculateGaussPointData(
        ConditionData& rData,
        const IndexType PointIndex,
        const array_1d<double,TNumNodes>& rN);

    void AddWaveTerms(
        LocalMatrixType& rMatrix,
        LocalVectorType& rVector,
        const ConditionData& rData,
        const array_1d<double,TNumNodes>& rN,
        const double Weight = 1.0);

    void AddFluxTerms(
        LocalVectorType& rVector,
        const ConditionData& rData,
        const array_1d<double,TNumNodes>& rN,
        const double Weight = 1.0);

    void AddMassTerms(
        LocalMatrixType& rMatrix,
        const ConditionData& rData,
        const array_1d<double,TNumNodes>& rN,
        const double Weight);

    const array_1d<double,3> VectorProduct(const array_1d<array_1d<double,3>,TNumNodes>& rV, const array_1d<double,TNumNodes>& rN) const;

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

}; // Class WaveCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_WAVE_CONDITION_H_INCLUDED  defined
