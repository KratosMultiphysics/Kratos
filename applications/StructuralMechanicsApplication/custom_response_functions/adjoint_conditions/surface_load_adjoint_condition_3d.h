// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder 
//

#if !defined(KRATOS_SURFACE_LOAD_ADJOINT_CONDITION_3D_H_INCLUDED )
#define  KRATOS_SURFACE_LOAD_ADJOINT_CONDITION_3D_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_conditions/surface_load_condition_3d.h"
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
    
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  SurfaceLoadAdjointCondition3D
    : public SurfaceLoadCondition3D
{
public:

    ///@name Type Definitions
    ///@{
    
    // Counted pointer of SurfaceLoadAdjointCondition3D
    KRATOS_CLASS_POINTER_DEFINITION( SurfaceLoadAdjointCondition3D );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    SurfaceLoadAdjointCondition3D();

    // Constructor using an array of nodes
    SurfaceLoadAdjointCondition3D( 
        IndexType NewId, 
        GeometryType::Pointer pGeometry 
        );

    // Constructor using an array of nodes with properties
    SurfaceLoadAdjointCondition3D( 
        IndexType NewId,
        GeometryType::Pointer pGeometry, 
        PropertiesType::Pointer pProperties 
        );

    // Destructor
    ~SurfaceLoadAdjointCondition3D() override;
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // Name Operations
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;
    
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties 
        ) const override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo ) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo ) override;

    void GetValuesVector(Vector& rValues, int Step = 0 ) override;   

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override; 

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;    

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    int Check( const ProcessInfo& rCurrentProcessInfo ) override; 
    

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

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SurfaceLoadCondition3D );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SurfaceLoadCondition3D);
    }


}; // class SurfaceLoadAdjointCondition3D.

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_ADJOINT_CONDITION_3D_H_INCLUDED  defined 
