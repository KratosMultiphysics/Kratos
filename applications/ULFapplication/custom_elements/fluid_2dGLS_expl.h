//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Pavel Ryzhakov and Julio Marti 
//



#if !defined(KRATOS_FLUID_2D_GLS_H_INCLUDED )
#define  KRATOS_FLUID_2D_GLS_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"

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
class Fluid2DGLS_expl
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Fluid2DGLS_expl
    KRATOS_CLASS_POINTER_DEFINITION(Fluid2DGLS_expl);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Fluid2DGLS_expl(IndexType NewId, GeometryType::Pointer pGeometry);
    Fluid2DGLS_expl(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~Fluid2DGLS_expl();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void Calculate(const Variable<double >& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

    void Calculate(const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo);


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
    virtual std::string Info() const
    {
        return "Fluid2DGLS_expl #" ;
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << Id();
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
    IntegrationMethod mThisIntegrationMethod;
    std::vector< Matrix > mInvJ0;
    Vector mDetJ0;
    ///@name Static Member Variables
    ///@{
//		static boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors;
//		static boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
//  		static array_1d<double,3> msN; //dimension = number of nodes
    //static Matrix msDN_DX;
    //static Matrix msMassFactors;
//		static array_1d<double,2> ms_vel_gauss; //dimesion coincides with space dimension
//		static array_1d<double,3> ms_temp_vec_np; //dimension = number of nodes
//		static array_1d<double,3> ms_u_DN;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    Fluid2DGLS_expl() : Element()
    {
    }

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Operators
    ///@{
    void CalculateLumpedMass();

    void CalculateGalerkinMomentumResidual(VectorType& Galerkin_RHS);
    //calculates stabilized RHS and writes it directly to the nodes
    void CalculateRHSVector(VectorType& GalerkinRHS, double& dt);

    //void RungeKuttaTimeIntegration(VectorType& RHS, double& d_t);
    //SECOND STEP OF FRAC STEP IS DONE WITHIN CALCULATE LOCAL SYSTEM
    //void SecondStepOfFractionalStep(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, double& dt);
    void FinalFractionalStep(const ProcessInfo&);
    void ComputeTimeStep(double CFLNumber);

    //void Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex);
    //void Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);


    //inline void CalculateGeometryData(Matrix& msDN_DX, Vector& N, double& Area)
//	  inline void CalculateGeometryData(boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, array_1d<double,3>& N, double& Area);

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //Fluid2DGLS_expl& operator=(const Fluid2DGLS_expl& rOther);

    /// Copy constructor.
    //Fluid2DGLS_expl(const Fluid2DGLS_expl& rOther);


    ///@}

}; // Class Fluid2DGLS_expl

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    Fluid2DGLS_expl& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const Fluid2DGLS_expl& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_FLUID_2D_GLS_H_INCLUDED  defined 


