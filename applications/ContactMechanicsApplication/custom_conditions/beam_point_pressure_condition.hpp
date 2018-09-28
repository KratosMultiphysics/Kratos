//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2013 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_BEAM_PRESSURE_CONDITION_H_INCLUDED )
#define  KRATOS_BEAM_PRESSURE_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_conditions/beam_point_rigid_contact_condition.hpp"

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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) BeamPointPressureCondition
: public Condition
{
protected:


public:



   ///@name Type Definitions

    ///Tensor order 1 definition
    typedef Vector                              VectorType;
    typedef Node<3>::Pointer              PointPointerType;

    std::vector<PointPointerType> mNodesList;
    array_1d<double, 3 > mForce;

    ///@{
    // Counted pointer of BeamPointRigidContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( BeamPointPressureCondition );
    ///@}

    /// Default constructor.
    BeamPointPressureCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    BeamPointPressureCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    BeamPointPressureCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);


    /// Copy constructor
    BeamPointPressureCondition( BeamPointPressureCondition const& rOther);


    /// Destructor.
    virtual ~BeamPointPressureCondition();


    Condition::Pointer Create(IndexType NewId, NodesArrayType const&
                              ThisNodes,  PropertiesType::Pointer pProperties) const;
    virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);
    virtual void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);
    void         GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo);
    void         EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
    void AddExplicitContribution(const VectorType& rRHSVector,
						 const Variable<VectorType>& rRHSVariable,
						 Variable<array_1d<double,3> >& rDestinationVariable,
						 const ProcessInfo& rCurrentProcessInfo);
    void InitializeSystemMatrices( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector);
    void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );
    void CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo );
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );
    void GetValuesVector(Vector& rValues, int Step);
    void GetFirstDerivativesVector( Vector& rValues, int Step );
    void GetSecondDerivativesVector( Vector& rValues, int Step );


protected:

    BeamPointPressureCondition() {};

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition )
    }


}; // Class BeamPointPressureCondition


}  // namespace Kratos.

#endif // KRATOS_BEAM_PRESSURE_CONDITION_H_INCLUDED  defined
