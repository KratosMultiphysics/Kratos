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


#if !defined(KRATOS_FLUID_3Dp_GLS_H_INCLUDED )
#define  KRATOS_FLUID_3Dp_GLS_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
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
  class QFluid3D
    : public Element
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of QFluid3DCoupled
    KRATOS_CLASS_POINTER_DEFINITION(QFluid3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    QFluid3D(IndexType NewId, GeometryType::Pointer pGeometry);
    QFluid3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~QFluid3D() override;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo) override;

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

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


    ///@}

  private:
    ///@name Static Member Variables
    ///@{
    static boost::numeric::ublas::bounded_matrix<double,4,4> msaux_matrix;
    static boost::numeric::ublas::bounded_matrix<double,4,4> msMassFactors;
    static boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX; //cartesian coords
    static array_1d<double,4> msN; //dimension = number of nodes
    static array_1d<double,3> ms_aux; //dimesion coincides with space dimension
    static array_1d<double,3> ms_vel_gauss; //dimesion coincides with space dimension
    static array_1d<double,4> ms_temp_vec_np; //dimension = number of nodes
    static array_1d<double,4> ms_u_DN;

    static void InitializeAuxiliaries()
    {
      msMassFactors = ZeroMatrix(4,4);
      msMassFactors(0,0) = 0.25;
      msMassFactors(1,1) = 0.25;
      msMassFactors(2,2) = 0.25;
      msMassFactors(3,3) = 0.25;
    }

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    QFluid3D() {}


    ///@}
    ///@name Private Operators
    ///@{
    void Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    void Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    inline double CalculateH(double Volume);

    void CalculateViscousMatrix(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX, const double& nu,const double& dt );

    inline void  ExpandAndAddReducedMatrix(MatrixType& Destination,	boost::numeric::ublas::bounded_matrix<double,4,4>& ReducedMatrix, const unsigned int dimension);
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
    //QFluid3DCoupled& operator=(const QFluid3DCoupled& rOther);

    /// Copy constructor.
    //QFluid3DCoupled(const QFluid3DCoupled& rOther);


    ///@}

  }; // Class QFluid3DCoupled

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  /*  inline std::istream& operator >> (std::istream& rIStream,
      QFluid3DCoupled& rThis);
  */
  /// output stream function
  /*  inline std::ostream& operator << (std::ostream& rOStream,
      const QFluid3DCoupled& rThis)
      {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
      }*/
  ///@}

}  // namespace Kratos.

#endif // KRATOS_FLUID_3D_COUPLED_H_INCLUDED  defined
