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
 

#if !defined(KRATOS_CONVDIFF_ELEM_3D_H_INCLUDED )
#define  KRATOS_CONVDIFF_ELEM_3D_H_INCLUDED



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
  class Rad3D
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of ConvDiff3D
      typedef GeometryData::IntegrationMethod IntegrationMethod;
      KRATOS_CLASS_POINTER_DEFINITION(Rad3D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Rad3D(IndexType NewId, GeometryType::Pointer pGeometry);
      Rad3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~Rad3D();
      
      IntegrationMethod GetIntegrationMethod1();
      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

      void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

      inline double CalculateH(boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, double Volume);
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
// 		static boost::numeric::ublas::bounded_matrix<double,4,4> msMassFactors;
// 		static boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX;
//   		static array_1d<double,4> msN; //dimension = number of nodes
// 		static array_1d<double,3> ms_vel_gauss; //dimesion coincides with space dimension
//   		static array_1d<double,4> ms_temp_vec_np; //dimension = number of nodes
// 		static array_1d<double,4> ms_u_DN;
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
		
        IntegrationMethod mThisIntegrationMethod;
        std::vector< Matrix > mInvJ0;
        Vector mDetJ0;
        
    friend class Serializer;

    // A private default constructor necessary for serialization
    Rad3D() : Element()
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
    
    template<class T>                                     
    bool InvertMatrix(const T& input, T& inverse)  ;    

    double CalculateTriangleArea3D(	array_1d<double,3>& Point1, array_1d<double,3>& Point2, array_1d<double,3>& Point3	);
   
    double Length(array_1d<double,3>& Point1, array_1d<double,3>& Point2);
    
    void qi( const unsigned int ndivisionsp, std::vector< Matrix > edges_tauxp, std::vector< Matrix > nodes_auxp , std::vector< Matrix > rGradientauxp,  array_1d<double,6> conductivitiesp, Matrix& Kaux1p,Matrix& Lenrichaux1p);

    void Heat_Source(VectorType& rRightHandSideVector,const int ndivisionsp, array_1d<double,6>& volumesp,array_1d<double,6>& conductivitiesp,boost::numeric::ublas::bounded_matrix<double,6, 8 >& Ngaussnewp, const double heatp);   
 

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
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      //ConvDiff3D& operator=(const ConvDiff3D& rOther);

      /// Copy constructor.
      //ConvDiff3D(const ConvDiff3D& rOther);

        
      ///@}    
        
    }; // Class ConvDiff3D 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    ConvDiff3D& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ConvDiff3D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_CONVDIFF_ELEM_3D_H_INCLUDED  defined 


