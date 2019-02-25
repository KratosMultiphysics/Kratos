//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//


//#define FIRST

// System includes

// External includes


#if !defined(KRATOS_FLUID_2Dp_GLS_H_INCLUDED )
#define  KRATOS_FLUID_2Dp_GLS_H_INCLUDED



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

  /// This element implements a Multi-stage element (2D case) to be used in conjuntion with @see FractionalStepStrategy

  /** The element is designed for the solution of the Navier-Stokes equations. Velocity components are considered to be uncoupled, and
   * laplacian formulation is used for the viscous term.
   * OSS (Orthogonal Sub-grid Scale) stabilization is used for both the incompressibility constraint and for the convective term.
   * smagorinsky turbulence model is optionally implemented and controlled by the value of the C_SMAGORINSKY constant, which is passed thorugh the
   * Processinfo.
   * The computation of the "tau" used in the stabilization allows the user to take in account a term depending on 1/Dt
   * this option is controlled by the variable ProcessInfo[DYNAMIC_TAU]. Setting it to 0.0 implies NOT considering a dependence
   * of tau on Dt.
   * The class is organized mainly in 3 stages
   * Stage1 - computes the velocity (designed for non-linear iteration)
   * Stage2 - computes the pressure
   * Stage3 - corrects the velocity taking in account the pressure variation computed in the second step
   */
  class QFluid2D
    : public Element
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of QFluid2D
    KRATOS_CLASS_POINTER_DEFINITION(QFluid2D);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    QFluid2D(IndexType NewId, GeometryType::Pointer pGeometry);
    QFluid2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~QFluid2D() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    //void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
    //void CalculateViscousMatrix(MatrixType& K, const BoundedMatrix<double,3,2>& DN_DX, const double& nu);
    void CalculateViscousMatrix(MatrixType& K  , const BoundedMatrix<double,3,2>& DN_DX, const double& nu, const double& deltat, const double& bulk); //, const double& k);

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

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

    virtual std::string Info() const override
    {
      return "QFluid2D #";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const override
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
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    IntegrationMethod mThisIntegrationMethod;
    std::vector< Matrix > mInvJ0;
    Vector mDetJ0;


    ///@}
    ///@name Private Operators
    ///@{
    void Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);//, unsigned int ComponentIndex);
    void Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    // inline double CalculateTau(BoundedMatrix<double, 3,2 > & DN_DX, array_1d<double, 2 > & vel_gauss, const double h, const double nu, const double norm_u, const ProcessInfo& CurrentProcessInfo);
    double ComputeSmagorinskyViscosity(const BoundedMatrix<double, 3, 2 > & DN_DX,
                                       const double& h,
                                       const double& C,
                                       const double nu
				       );

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

    QFluid2D() : Element()
    {
    }

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //QFluid2D& operator=(const QFluid2D& rOther);

    /// Copy constructor.
    //QFluid2D(const QFluid2D& rOther);


    ///@}

  }; // Class QFluid2D

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  /*  inline std::istream& operator >> (std::istream& rIStream,
      QFluid2D& rThis);
  */
  /// output stream function
  /*  inline std::ostream& operator << (std::ostream& rOStream,
      const QFluid2D& rThis)
      {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
      }*/
  ///@}

} // namespace Kratos.

#endif // KRATOS_FLUID_2D_H_INCLUDED  defined
