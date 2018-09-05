//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LINE_HEAT_FLUX_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_HEAT_FLUX_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"


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

class LineHeatFluxCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of LineHeatFluxCondition
    KRATOS_CLASS_POINTER_DEFINITION( LineHeatFluxCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LineHeatFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    LineHeatFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    ~LineHeatFluxCondition() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const override;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo ) override;

    void GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
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

    void CalculateElementalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                   ProcessInfo& rCurrentProcessInfo,
                                   bool CalculateStiffnessMatrixFlag,
                                   bool CalculateResidualVectorFlag);


    void CalculateAndSubKheatflux(Matrix& rK,
                                  const Matrix& rDN_De,
                                  const Vector& rN,
                                  double rFlux,
                                  double rIntegrationWeight);


    void CalculateAndAddFaceHeatFlux (Vector& rF,
                                      const Vector& rN,
                                      Vector& rNormal,
                                      double rFlux,
                                      double rIntegrationWeight );

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

    // A private default constructor necessary for serialization
    LineHeatFluxCondition() {};

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //LineHeatFluxCondition& operator=(const LineHeatFluxCondition& rOther);

    /// Copy constructor.
    //LineHeatFluxCondition(const LineHeatFluxCondition& rOther);


    ///@}

}; // Class LineHeatFluxCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  // namespace Kratos.

#endif // KRATOS_LINE_HEAT_FLUX_CONDITION_H_INCLUDED  defined


