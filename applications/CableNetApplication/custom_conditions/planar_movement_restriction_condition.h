//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           license: CABLE_NET_APPLICATION/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

#if !defined(PLANAR_MOVEMENT_RESTRICTION_CONDITION_H_INCLUDED )
#define  PLANAR_MOVEMENT_RESTRICTION_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "cable_net_application_variables.h"

namespace Kratos
{


/**
 * @class PlanarMovementRestrictionCondition3D1N
 */
class KRATOS_API(CABLE_NET_APPLICATION)  PlanarMovementRestrictionCondition3D1N
    : public Condition
{
public:


    /// We define the base class Condition
    typedef Condition BaseType;

    // Counted pointer of PlanarMovementRestrictionCondition3D1N
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( PlanarMovementRestrictionCondition3D1N );

    // Constructor void
    PlanarMovementRestrictionCondition3D1N()
    {};

    // Constructor using an array of nodes
    PlanarMovementRestrictionCondition3D1N( IndexType NewId, GeometryType::Pointer pGeometry ):Condition(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    PlanarMovementRestrictionCondition3D1N( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):Condition(NewId,pGeometry,pProperties)
    {};

    ///Copy constructor
    PlanarMovementRestrictionCondition3D1N(PlanarMovementRestrictionCondition3D1N const& rOther);

    // Destructor
    ~PlanarMovementRestrictionCondition3D1N() override
    {};


    /// Assignment operator.
    PlanarMovementRestrictionCondition3D1N& operator=(PlanarMovementRestrictionCondition3D1N const& rOther);

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    Condition::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& ThisNodes
        ) const override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    
    void GetDofList(
        DofsVectorType& ElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    void GetValuesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo
        ) override;


    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;


    void AddExplicitContribution(const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;


    double CalculateNormalDistance(const Vector& rCurrentDisplacements) const;

    int Check( const ProcessInfo& rCurrentProcessInfo ) const override;


protected:

    static constexpr int msNumberOfNodes = 1;
    static constexpr int msDimension = 3;
    static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;

private:
  
    friend class Serializer;

    void save( Serializer& rSerializer ) const override;

    void load( Serializer& rSerializer ) override;

}; // class PlanarMovementRestrictionCondition3D1N.


} // namespace Kratos.

#endif // PLANAR_MOVEMENT_RESTRICTION_CONDITION_H_INCLUDED  defined
