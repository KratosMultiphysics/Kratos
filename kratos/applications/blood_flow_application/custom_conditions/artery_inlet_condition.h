//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_ARTERY_INLET_CONDITION_H_INCLUDED )
#define  KRATOS_ARTERY_INLET_CONDITION_H_INCLUDED



// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "includes/serializer.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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
class ArteryInletCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ArteryInletCondition
    KRATOS_CLASS_POINTER_DEFINITION(ArteryInletCondition);

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructor.
    ArteryInletCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    ArteryInletCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~ArteryInletCondition();


    ///@}
    ///@name Operators
    ///@{
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

    int Check(const ProcessInfo& rCurrentProcessInfo);
	
	void Initialize();


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{
	
	     ///@}
     ///@name Serialization
     ///@{	
	
		friend class Serializer;     

        virtual void save(Serializer& rSerializer) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ArteryInletCondition #";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << Id();
    }



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

    array_1d<double,1> mInitialArea;
	double mH0;
        double mBeta;


    ///@}
    ///@name Private Operators
    ///@{
    ArteryInletCondition() : Condition(){}


    ///@}
    ///@name Private Operations
    ///@{

    double UpdateArea(double Beta, double Density, double Flow);


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
 //   ArteryInletCondition& operator=(ArteryInletCondition const& rOther) {};

    /// Copy constructor.
//    ArteryInletCondition(ArteryInletCondition const& rOther) {};


    ///@}

}; // Class ArteryInletCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ArteryInletCondition& rThis)
{
    return rIStream;
};

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ArteryInletCondition& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ARTERY_INLET_CONDITION_H_INCLUDED  defined


