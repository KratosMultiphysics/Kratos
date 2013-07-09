//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_FACE_FORCE_2D_CONDITION_H_INCLUDED )
#define  KRATOS_FACE_FORCE_2D_CONDITION_H_INCLUDED



// System includes


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

class FaceForce2DCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of FaceForce2DCondition
    KRATOS_CLASS_POINTER_DEFINITION( FaceForce2DCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FaceForce2DCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    FaceForce2DCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~FaceForce2DCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo );

    void GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo );

    void MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo );

    void DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo );

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );
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
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Face Force 2D Condition #" << Id();
    }


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
    void CalculateConditionalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                   ProcessInfo& rCurrentProcessInfo,
                                   bool CalculateStiffnessMatrixFlag,
                                   bool CalculateResidualVectorFlag );

    void CalculateAndSubKp(Matrix& rK,
                           const Matrix& rDN_De,
                           const Vector& rN,
                           double rPressure,
                           double rIntegrationWeight
                          );


    void CalculateAndAddFacePressure (Vector& rF,
                                      const Vector& rN,
                                      Vector& rNormal,
                                      double rPressure,
                                      double rIntegrationWeight );

    void CalculateAndAddFaceForce(Vector& rF,
                                  const Vector& rN,
                                  Vector& rForce,
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
    FaceForce2DCondition() {};

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //FaceForce2DCondition& operator=(const FaceForce2DCondition& rOther);

    /// Copy constructor.
    //FaceForce2DCondition(const FaceForce2DCondition& rOther);


    ///@}

}; // Class FaceForce2DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        FaceForce2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const FaceForce2DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_FACE_FORCE_2D_CONDITION_H_INCLUDED  defined 


