#pragma once 

#include "custom_elements/lagrangian_particle.h"
#include "includes/define.h"

namespace Kratos
{
template<class TKernelType>
class SmallDisplacementParticle2D : LagrangianParticle<TKernelType>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementParticle2D);

    using BaseType = SmallDisplacementParticle2D;

    // Constructor void 
    SmallDisplacementParticle2D();

    // Constructor using an array of nodes 
    SmallDisplacementParticle2D(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry);

    // Constructor using an array of nodes with properties 
    SmallDisplacementParticle2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties);

    // Copy constructor
    SmallDisplacementParticle2D(SmallDisplacementParticle2D const& rOther)
        : BaseType(rOther), mThisConstitutiveLaw(rOther.mThisConstitutiveLaw){};

    /// Destructor.
    ~SmallDisplacementParticle2D() override;

protected:

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS matrix
     * @param rRightHandSideVector The RHS vector
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        );

};

}