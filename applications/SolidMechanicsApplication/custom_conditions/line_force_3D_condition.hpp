//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LINE_FORCE_3D_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_FORCE_3D_CONDITION_H_INCLUDED



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

class LineForce3DCondition
    : public Condition
{
public:

    // Counted pointer of LineForce3DCondition
    KRATOS_CLASS_POINTER_DEFINITION( LineForce3DCondition );


    // Constructor void
    LineForce3DCondition();

    // Constructor using an array of nodes
    LineForce3DCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    LineForce3DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~LineForce3DCondition();


    // Name Operations

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo );

    void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo );

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo );

    void MassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo );

    void DampMatrix(
        MatrixType& rDampMatrix,
        ProcessInfo& rCurrentProcessInfo );

//    void GetValuesVector(
//        Vector& values,
//        int Step = 0 );

//    void GetFirstDerivativesVector(
//        Vector& values,
//        int Step = 0 );

//    void GetSecondDerivativesVector(
//        Vector& values,
//        int Step = 0 );

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );


protected:


private:
    ///@name Static Member Variables

    /// privat variables


    // privat name Operations

//    void CalculateAll(
//        MatrixType& rLeftHandSideMatrix,
//        VectorType& rRightHandSideVector,
//        const ProcessInfo& rCurrentProcessInfo,
//        bool CalculateStiffnessMatrixFlag,
//        bool CalculateResidualVectorFlag );

//    void CalculateAndSubKp(
//        Matrix& K,
//        array_1d<double, 3>& ge,
//        array_1d<double, 3>& gn,
//        const Matrix& DN_De,
//        const Vector& N,
//        double pressure,
//        double weight );

//    void MakeCrossMatrix(
//        boost::numeric::ublas::bounded_matrix<double, 3, 3>& M,
//        array_1d<double, 3>& U );

//    void CrossProduct(
//        array_1d<double, 3>& cross,
//        array_1d<double, 3>& a,
//        array_1d<double, 3>& b );

//    void SubtractMatrix(
//        MatrixType& Destination,
//        boost::numeric::ublas::bounded_matrix<double, 3, 3>& InputMatrix,
//        int InitialRow,
//        int InitialCol );

//    void ExpandReducedMatrix(
//        Matrix& Destination,
//        Matrix& ReducedMatrix );

//    void CalculateAndAdd_PressureForce(
//        VectorType& residualvector,
//        const Vector& N,
//        const array_1d<double, 3>& v3,
//        double pressure,
//        double weight,
//        const ProcessInfo& rCurrentProcessInfo );

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

}; // class LineForce3DCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_FORCE_3D_H_INCLUDED  defined
