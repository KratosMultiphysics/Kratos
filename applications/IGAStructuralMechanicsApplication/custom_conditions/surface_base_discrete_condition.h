#if !defined(KRATOS_SURFACE_BASE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_SURFACE_BASE_DISCRETE_CONDITION_H_INCLUDED


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
/** Discrete surface element deals as base class for thin walled structures.
KRATOS_API(IGA_STRUCTURAL_MECHANICS_APPLICATION)
*/
class  SurfaceBaseDiscreteCondition
    : public BaseDiscreteCondition
{
protected:
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SurfaceBaseDiscreteCondition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SurfaceBaseDiscreteCondition #" << Id();
    }

public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of SurfaceBaseDiscreteCondition
    KRATOS_CLASS_POINTER_DEFINITION(SurfaceBaseDiscreteCondition);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
     // Constructor using an array of nodes
    SurfaceBaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteCondition(NewId, pGeometry)
    {};
     // Constructor using an array of nodes with properties
    SurfaceBaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : BaseDiscreteCondition(NewId, pGeometry, pProperties)
    {};

    SurfaceBaseDiscreteCondition() : BaseDiscreteCondition()
    {};

    /// Destructor.
    virtual ~SurfaceBaseDiscreteCondition() override
    {};

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_ERROR << "Trying to create a \"SurfaceBaseDiscreteCondition\"" << std::endl;
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * Called to initialize the element.
    * Must be called before any calculation is done
    */
    void Initialize() override;

    ///@}
protected:

    ///@name Operations
    ///@{

    /**
    * This functions calculates both the RHS and the LHS
    * @param rLeftHandSideMatrix: The LHS
    * @param rRightHandSideVector: The RHS
    * @param rCurrentProcessInfo: The current process info instance
    * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
    * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
    */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    /**
    * Hessian calculates the Hessian for the given system with the
    shape function derivatives DN_De.
    *
    * @param[in] DN_De derivatives of shape functions.
    * @param[out] Hessian calculated Hessian. Is always of size 3x3.
    */
    void CalculateHessian(
        Matrix& Hessian,
        const Matrix& DDN_DDe,
        const int rDimension = 3);


    void CalculateBaseVector(
        Vector& rBaseVector,
        const Matrix& rDN_De);
    ///@}

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
};     // Class SurfaceBaseDiscreteCondition
///@}
}  // namespace Kratos.

#endif // KRATOS_SURFACE_BASE_DISCRETE_CONDITION_H_INCLUDED  defined