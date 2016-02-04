//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                 $
//   Revision:            $Revision:               $
//

#if !defined(KRATOS_POINT_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_LOAD_CONDITION_H_INCLUDED

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

class PointLoadCondition : public Condition
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION( PointLoadCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    PointLoadCondition();
    
    // Constructor 1
    PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    //Constructor 2
    PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    //Destructor
    virtual ~PointLoadCondition();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const;
    
    void GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );

    void EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, 
                                bool CalculateLHSMatrixFlag, bool CalculateResidualVectorFlag);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }
    
}; // class PointLoadCondition.

} // namespace Kratos.

#endif // KRATOS_POINT_LOAD_CONDITION_H_INCLUDED defined 
