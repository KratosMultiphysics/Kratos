// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder
//

// System includes
#if !defined(KRATOS_POINT_LOAD_CONDITION_FOR_SA_H_INCLUDED )
#define  KRATOS_POINT_LOAD_CONDITION_FOR_SA_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_conditions/base_load_condition.h"
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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  PointLoadConditionForSA
    : public BaseLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PointLoadConditionForSA
    KRATOS_CLASS_POINTER_DEFINITION( PointLoadConditionForSA );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointLoadConditionForSA( 
        IndexType NewId, 
        GeometryType::Pointer pGeometry 
        );
    
    PointLoadConditionForSA( 
        IndexType NewId, 
        GeometryType::Pointer pGeometry,  
        PropertiesType::Pointer pProperties 
        );

    /// Destructor.
    ~PointLoadConditionForSA() override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @param rResult: The vector containing the equation id
     * @param rCurrentProcessInfo: The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo 
        ) override;
    
    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList: The vector containing the dof of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Sets on rValues the nodal displacements
     * @param rValues: The values of displacements
     * @param Step: The step to be computed
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0 
        ) override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
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

    /**
     * this is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				       ProcessInfo& rCurrentProcessInfo) override; 
    //fusseder TODO: maybe remove this function to base_load_condition                   

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

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;    

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;                                  

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
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        ) override;


        
    /**
     * It calcules the integration load for the point load 
     */
    virtual double GetPointLoadIntegrationWeight();
        
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
    PointLoadConditionForSA() {};

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseLoadCondition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseLoadCondition );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //PointLoadConditionForSA& operator=(const PointLoadConditionForSA& rOther);

    /// Copy constructor.
    //PointLoadConditionForSA(const PointLoadConditionForSA& rOther);


    ///@}

}; // Class PointLoadConditionForSA

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        PointLoadConditionForSA& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const PointLoadConditionForSA& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_LOAD_CONDITION_FOR_SA_H_INCLUDED  defined 


