#if !defined(KRATOS_CURVE_BASE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_CURVE_BASE_DISCRETE_CONDITION_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"

// External includes

// Project includes
#include "custom_conditions/base_discrete_condition.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** CurveBaseDiscreteCondition deals as base class for curve structured condition formulations.
*/
class  CurveBaseDiscreteCondition
    : public BaseDiscreteCondition
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of CurveBaseDiscreteCondition
    KRATOS_CLASS_POINTER_DEFINITION(CurveBaseDiscreteCondition);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
     // Constructor using an array of nodes
    CurveBaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteCondition(NewId, pGeometry)
    {};
     // Constructor using an array of nodes with properties
    CurveBaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : BaseDiscreteCondition(NewId, pGeometry, pProperties)
    {};

    CurveBaseDiscreteCondition() : BaseDiscreteCondition()
    {};

    /// Destructor.
    virtual ~CurveBaseDiscreteCondition() override
    {};

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_ERROR << "Trying to create a \"CurveBaseDiscreteCondition\"" << std::endl;
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * Called to initialize the condition.
    * Must be called before any calculation is done
    */
    void Initialize() override;



    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"CurveBaseDiscreteCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"CurveBaseDiscreteCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
protected:
    ///@name Static Member Variables
    ///@{
    Vector mBaseVector0;

    ///@}
    ///@name Operations
    ///@{

    void CalculateHessianSurface(
        Matrix& Hessian, 
        const Matrix& DDN_DDe, 
        const int rDimension);

    void CalculateBaseVector(
        Vector& rBaseVector, 
        const Matrix& rDN_De);

    void GetBaseVectorsSurface(
        const Matrix& DN_De,
        Vector& g1,
        Vector& g2,
        Vector& g3);

    /**
    * GetBoundaryEdgeBaseVector computes t3 of the boundary edge
    * @param DN_De derivatives of shape functions.
    * @param Tangents in Parameter space
    * @see rBaseVector t3 of the edge
    */
    void GetBoundaryEdgeBaseVector(const Matrix& DN_De,
        const array_1d<double, 2>& Tangents,
        Vector& rBaseVector);

    void CalculateNormalVector(
        Vector& rNormalVector,
        const Matrix& rDN_De
    );

private:
    ///@name Operations
    ///@{

    ///@}

    ///@name Static Member Variables

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseDiscreteCondition)
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseDiscreteCondition)
    }
    ///@}
};     // Class CurveBaseDiscreteCondition
///@}
}  // namespace Kratos.

#endif // KRATOS_CURVE_BASE_DISCRETE_CONDITION_H_INCLUDED  defined