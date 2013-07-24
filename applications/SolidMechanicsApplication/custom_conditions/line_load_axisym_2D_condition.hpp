//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LINE_LOAD_AXISYM_2D_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_LOAD_AXISYM_2D_CONDITION_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_conditions/line_load_2D_condition.hpp"

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

class LineLoadAxisym2DCondition
    : public LineLoad2DCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of LineLoadAxisym2DCondition
    KRATOS_CLASS_POINTER_DEFINITION( LineLoadAxisym2DCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LineLoadAxisym2DCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    LineLoadAxisym2DCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Copy constructor
    LineLoadAxisym2DCondition( LineLoadAxisym2DCondition const& rOther);


    /// Destructor.
    virtual ~LineLoadAxisym2DCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    //virtual int Check( const ProcessInfo& rCurrentProcessInfo );
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
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Line Load Axisym 2D Condition #" << Id();
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

    void CalculateConditionalSystem( MatrixType& rLeftHandSideMatrix, 
				     VectorType& rRightHandSideVector,
				     ProcessInfo& rCurrentProcessInfo,
				     bool CalculateStiffnessMatrixFlag,
				     bool CalculateResidualVectorFlag );



    void CalculateRadius(double & rCurrentRadius,
			 double & rReferenceRadius,
			 const Vector& rN);
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

    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;

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

    // A private default constructor necessary for serialization
    LineLoadAxisym2DCondition() {};

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
    //LineLoadAxisym2DCondition& operator=(const LineLoadAxisym2DCondition& rOther);

    /// Copy constructor.
    //LineLoadAxisym2DCondition(const LineLoadAxisym2DCondition& rOther);


    ///@}

}; // Class LineLoadAxisym2DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        LineLoadAxisym2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const LineLoadAxisym2DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_LINE_LOAD_2D_CONDITION_H_INCLUDED  defined 


