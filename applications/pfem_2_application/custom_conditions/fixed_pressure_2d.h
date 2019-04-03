#if !defined(KRATOS_FIXED_PRESSURE_CONDITION_H_INCLUDED )
#define  KRATOS_FIXED_PRESSURE_CONDITION_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{
class FixedPressure2D : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of FixedPressure2D
    KRATOS_CLASS_POINTER_DEFINITION(FixedPressure2D);

    /// Default constructor.
    FixedPressure2D(IndexType NewId, GeometryType::Pointer pGeometry);

    FixedPressure2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~FixedPressure2D() override;


    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo) override;

protected:


private:

    friend class Serializer;

    // A private default constructor necessary for serialization
    FixedPressure2D() : Condition()
    {
    }


  }; // Class FixedPressure2D

} //namespace kratos
#endif
