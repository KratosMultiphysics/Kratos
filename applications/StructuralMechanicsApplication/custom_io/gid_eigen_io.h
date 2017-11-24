//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher 
//

#if !defined(KRATOS_GID_EIGEN_IO_H_INCLUDED )
#define  KRATOS_GID_EIGEN_IO_H_INCLUDED


// System includes
#include <string>
#include <iostream> 

// External includes 

// Project includes
#include "includes/gid_io.h"


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
//   template<class TGaussPointContainer = GidGaussPointsContainer, class TMeshContainer = GidMeshContainer>
  class GidEigenIO : public GidIO<>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of GidEigenIO
      KRATOS_CLASS_POINTER_DEFINITION(GidEigenIO);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      GidEigenIO( const std::string& rDatafilename,
                  GiD_PostMode Mode,
                  MultiFileFlag use_multiple_files_flag,
                  WriteDeformedMeshFlag write_deformed_flag,
                  WriteConditionsFlag write_conditions_flag) : 
        GidIO<>(rDatafilename, Mode, use_multiple_files_flag, write_deformed_flag, write_conditions_flag)
        {}

      /// Destructor.
      virtual ~GidEigenIO(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
        
        /**
     * writes nodal results for variables of type array_1d<double, 3>
     * (e.g. DISPLACEMENT)
     */
    void WriteEigenResults( Variable<array_1d<double, 3> > const& rVariable,
                            NodesContainerType& rNodes,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {

        Timer::Start("Writing Eigen Results");

        GiD_fBeginResult(mResultFile,(char*)(rVariable.Name().c_str()), "Kratos_Eigen",
                         SolutionTag, GiD_Vector,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable,
                                        SolutionStepNumber );
            GiD_fWriteVector( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Eigen Results");

    }
      
      
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
	std::stringstream buffer;
        buffer << "GidEigenIO" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "GidEigenIO";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}
      
            
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
    //   GidEigenIO& operator=(GidEigenIO const& rOther){}

      /// Copy constructor.
    //   GidEigenIO(GidEigenIO const& rOther){}

        
      ///@}    
        
    }; // Class GidEigenIO 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    GidEigenIO& rThis){}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const GidEigenIO& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_GID_EIGEN_IO_H_INCLUDED  defined 


