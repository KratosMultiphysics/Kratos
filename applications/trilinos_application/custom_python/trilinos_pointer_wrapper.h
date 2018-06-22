//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					     Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_TRILINOS_POINTER_WRAPPER_H_INCLUDED)
#define KRATOS_TRILINOS_POINTER_WRAPPER_H_INCLUDED

// System includes

// External includes

//Trilinos includes

// Project includes
#include "trilinos_application.h"
#include "trilinos_space.h"

namespace Kratos
{

///@addtogroup TrilinosApplication
///@{

///@name Kratos Globals
///@{ 

///@} 
///@name Type Definitions
///@{ 

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;

///@} 
///@name  Enum's
///@{
    
///@}
///@name  Functions 
///@{
    
///@}
///@name Kratos Classes
///@{

/// Auxiliarty Trilinos matrix pointer wrapper class
/** This class is intended to handle the Trilinos space matrices exportation to Python.
 * Since the matrix pointer cannot be exported to Python by means of PyBind, it is stored
 * in this auxilary class, which is the object that is exported to python. The exportation
 * of the matrix wrapper is done in add_trilinos_space_to_python.cpp. Then, the matrix is
 * retrieved in Python by calling the GetReference() method of the AuxiliaryMatrixWrapper.
 * Alternatively, one can define an auxiliary function that does this operation before
 * the Python export.
*/
class AuxiliaryMatrixWrapper
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename TrilinosSparseSpaceType::MatrixType TrilinosMatrixType;
    typedef typename TrilinosSparseSpaceType::MatrixPointerType TrilinosMatrixPointerType;
  
    ///@}
    ///@name Life Cycle 
    ///@{ 
      
    /// Default constructor.
    AuxiliaryMatrixWrapper(TrilinosMatrixPointerType p) : mp(p){};

    /// Destructor.
    virtual ~AuxiliaryMatrixWrapper(){}
      
    ///@}
    ///@name Operators 
    ///@{
        
    ///@}
    ///@name Operations
    ///@{
     
    /**
     * Get the matrix pointer
     * @return the member pointer to the matrix
     */
    TrilinosMatrixPointerType& GetPointer() { return mp; }

    /**
     * Get a reference to the matrix in the wrapper
     * @return the refernce to the wrapper matrix
     */
    TrilinosMatrixType& GetReference() { return *mp; }

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
    virtual std::string Info() const {
	    std::stringstream buffer;
        buffer << "AuxiliaryMatrixWrapper" ;
        return buffer.str();
    }
      
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "AuxiliaryMatrixWrapper";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}
      
private:
    ///@}      
    ///@name Friends
    ///@{
       
    ///@}
    ///@name Static Member Variables 
    ///@{ 
          
    ///@} 
    ///@name Member Variables 
    ///@{ 

    TrilinosMatrixPointerType mp;
        
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
    AuxiliaryMatrixWrapper& operator=(const AuxiliaryMatrixWrapper &rOther) = delete;

    ///@}        
}; // Class AuxiliaryMatrixWrapper

/// Auxiliarty Trilinos vector pointer wrapper class
/** This class is intended to handle the Trilinos space matrices exportation to Python.
 * Since the vector pointer cannot be exported to Python by means of PyBind, it is stored
 * in this auxilary class, which is the object that is exported to python. The exportation
 * of the vector wrapper is done in add_trilinos_space_to_python.cpp. Then, the vector is
 * retrieved in Python by calling the GetReference() method of the AuxiliaryVectorWrapper.
 * Alternatively, one can define an auxiliary function that does this operation before
 * the Python export.
*/
class AuxiliaryVectorWrapper
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename TrilinosSparseSpaceType::VectorType TrilinosVectorType;
    typedef typename TrilinosSparseSpaceType::VectorPointerType TrilinosVectorPointerType;
      
    ///@}
    ///@name Life Cycle 
    ///@{ 
      
    /// Default constructor.
    AuxiliaryVectorWrapper(TrilinosVectorPointerType p) : mp(p){};

    /// Destructor.
    virtual ~AuxiliaryVectorWrapper(){}
      
    ///@}
    ///@name Operators 
    ///@{
        
    ///@}
    ///@name Operations
    ///@{
     
    /**
     * Get the matrix pointer
     * @return the member pointer to the matrix
     */
    TrilinosVectorPointerType& GetPointer() { return mp; }

    /**
     * Get a reference to the matrix in the wrapper
     * @return the refernce to the wrapper matrix
     */
    TrilinosVectorType& GetReference() { return *mp; }

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
    virtual std::string Info() const {
	    std::stringstream buffer;
        buffer << "AuxiliaryVectorWrapper" ;
        return buffer.str();
    }
      
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "AuxiliaryVectorWrapper";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}
      
private:
    ///@}      
    ///@name Friends
    ///@{
       
    ///@}
    ///@name Static Member Variables 
    ///@{ 
          
    ///@} 
    ///@name Member Variables 
    ///@{ 

    TrilinosVectorPointerType mp;
        
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
    AuxiliaryVectorWrapper& operator=(AuxiliaryVectorWrapper const& rOther) = delete;

    ///@}        
}; // Class AuxiliaryVectorWrapper

} // namespace Kratos.

#endif // KRATOS_TRILINOS_POINTER_WRAPPER_H_INCLUDED  defined