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

// External includes 

// Project includes
#include "includes/gid_io.h"


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
  
  /// GidIO specialized for writting Eigenvalue Results
  
  /** The main functionality of this class is to write the custom format 
  * that the postprocessing of the Eigenvalues creates. 
  * Also the result labels contain the Eigenvalue or -frequency
  */
class GidEigenIO : public GidIO<>
{
  public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of GidEigenIO
    KRATOS_CLASS_POINTER_DEFINITION(GidEigenIO);

    typedef std::size_t SizeType;    

    typedef std::vector<std::vector<double>> EigenResults;
    /* Structure of this EigenResults (results for one Eigenvalue)
    Nodes
        Coordinates
    */

    ///@}
    ///@name Life Cycle 
    ///@{ 
    
    /// Default constructor.
    GidEigenIO( const std::string& rDatafilename,
                GiD_PostMode Mode,
                MultiFileFlag use_multiple_files_flag,
                WriteDeformedMeshFlag write_deformed_flag,
                WriteConditionsFlag write_conditions_flag) : 
            GidIO<>(rDatafilename, 
                    Mode, 
                    use_multiple_files_flag, 
                    write_deformed_flag, 
                    write_conditions_flag) { }

    /// Destructor.
    ~GidEigenIO() = default;

    // Explicitly delete the other constructors
    GidEigenIO(const GidEigenIO&) = delete;
    GidEigenIO& operator=(const GidEigenIO&) = delete;

    ///@}
    ///@name Operators 
    ///@{
    
    
    ///@}
    ///@name Operations
    ///@{
    
    /**
    * Write the post-processed eigensolver-results
    * The label is the Eigenvalue or -frequency
    */
    template< class TVarType>
    inline void WriteEigenResults( ModelPart& rModelPart, 
                                   const TVarType& rVariable,
                                   std::string Label, 
                                   const SizeType AnimationStepNumber );
      
      
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

    
    ///@}    
    
}; // Class GidEigenIO 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 

template<>
inline void GidEigenIO::WriteEigenResults<Variable<double>>(ModelPart& rModelPart, 
                                                  const Variable<double>& rVariable,
                                                  std::string Label, 
                                                  const SizeType AnimationStepNumber) 
{
    Label += "_" + rVariable.Name();  
    GiD_fBeginResult( mResultFile, (char*)Label.c_str() , "EigenVector_Animation",
                        AnimationStepNumber, GiD_Scalar,
                        GiD_OnNodes, NULL, NULL, 0, NULL );

    for (const auto& node : rModelPart.Nodes())
    {
        const double& nodal_result = node.FastGetSolutionStepValue(rVariable);
        GiD_fWriteScalar( mResultFile, node.Id(), nodal_result );
    }

    GiD_fEndResult(mResultFile);
}

template<>
inline void GidEigenIO::WriteEigenResults<Variable<array_1d<double, 3>>>(ModelPart& rModelPart, 
                                                               const Variable<array_1d<double, 3>>& rVariable,
                                                               std::string Label, 
                                                               const SizeType AnimationStepNumber) 
{
    Label += "_" + rVariable.Name();
    GiD_fBeginResult( mResultFile, (char*)Label.c_str() , "EigenVector_Animation",
                        AnimationStepNumber, GiD_Vector,
                        GiD_OnNodes, NULL, NULL, 0, NULL );

    for (auto& node : rModelPart.Nodes())
    {
        const array_1d<double, 3>& nodal_result = node.FastGetSolutionStepValue(rVariable);
        GiD_fWriteVector(mResultFile, node.Id(), nodal_result[0], nodal_result[1], nodal_result[2]);
    }

    GiD_fEndResult(mResultFile);
}
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_GID_EIGEN_IO_H_INCLUDED  defined 


