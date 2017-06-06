// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_SURFACE_LOAD_CONDITION_3D_H_INCLUDED )
#define  KRATOS_SURFACE_LOAD_CONDITION_3D_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_conditions/base_load_condition.h"
#include "includes/ublas_interface.h"

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
    
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  SurfaceLoadCondition3D
    : public BaseLoadCondition
{
public:

    ///@name Type Definitions
    ///@{
    
    // Counted pointer of SurfaceLoadCondition3D
    KRATOS_CLASS_POINTER_DEFINITION( SurfaceLoadCondition3D );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    SurfaceLoadCondition3D();

    // Constructor using an array of nodes
    SurfaceLoadCondition3D( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    SurfaceLoadCondition3D( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~SurfaceLoadCondition3D();
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // Name Operations
    Condition::Pointer Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const;
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo );

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo );

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    
private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag );

    void CalculateAndSubKp(
        Matrix& K,
        array_1d<double, 3>& ge,
        array_1d<double, 3>& gn,
        const Matrix& DN_De,
        const Vector& N,
        double pressure,
        double weight );

    void MakeCrossMatrix(
        boost::numeric::ublas::bounded_matrix<double, 3, 3>& M,
        array_1d<double, 3>& U );

    void CalculateAndAdd_PressureForce(
        VectorType& residualvector,
        const Vector& N,
        const array_1d<double, 3>& v3,
        double pressure,
        double weight,
        const ProcessInfo& rCurrentProcessInfo );
    
    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseLoadCondition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseLoadCondition );
    }


}; // class SurfaceLoadCondition3D.

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_CONDITION_3D_H_INCLUDED  defined 
