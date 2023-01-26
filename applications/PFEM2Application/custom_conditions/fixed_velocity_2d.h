#if !defined(KRATOS_FIXED_VELOCITY_CONDITION_H_INCLUDED )
#define  KRATOS_FIXED_VELOCITY_CONDITION_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"


namespace Kratos
{
class FixedVelocity2D : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of FixedVelocity2D
    KRATOS_CLASS_POINTER_DEFINITION(FixedVelocity2D);

    /// Default constructor.
    FixedVelocity2D(IndexType NewId, GeometryType::Pointer pGeometry);

    FixedVelocity2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~FixedVelocity2D() override;


    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ConditionalDofList,const ProcessInfo& CurrentProcessInfo) const override;

protected:


private:

    friend class Serializer;

    // A private default constructor necessary for serialization
    FixedVelocity2D() : Condition()
    {
    }


}; // Class FixedVelocity2D

} //namespace kratos
#endif
