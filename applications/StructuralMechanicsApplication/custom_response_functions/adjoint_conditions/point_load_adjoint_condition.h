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

// System includes
#if !defined(KRATOS_POINT_LOAD_ADJOINT_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_LOAD_ADJOINT_CONDITION_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_conditions/point_load_condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  PointLoadAdjointCondition
    : public PointLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PointLoadAdjointCondition
    KRATOS_CLASS_POINTER_DEFINITION( PointLoadAdjointCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointLoadAdjointCondition( 
        IndexType NewId, 
        GeometryType::Pointer pGeometry 
        );
    
    PointLoadAdjointCondition( 
        IndexType NewId, 
        GeometryType::Pointer pGeometry,  
        PropertiesType::Pointer pProperties 
        );

    /// Destructor.
    ~PointLoadAdjointCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;
    
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const override;

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

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

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
    
    // A protected default constructor necessary for serialization
    PointLoadAdjointCondition(): PointLoadCondition(){};

    

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointLoadCondition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointLoadCondition );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //PointLoadAdjointCondition& operator=(const PointLoadAdjointCondition& rOther);

    /// Copy constructor.
    //PointLoadAdjointCondition(const PointLoadAdjointCondition& rOther);


    ///@}

}; // Class PointLoadAdjointCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        PointLoadAdjointCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const PointLoadAdjointCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_LOAD_ADJOINT_CONDITION_H_INCLUDED  defined 


