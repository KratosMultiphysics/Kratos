/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 

#if !defined(KRATOS_ILU_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_ILU_PRECONDITIONER_H_INCLUDED




// System includes 



// External includes 
#include "boost/smart_ptr.hpp"



// Project includes
#include "includes/define.h"



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
  
  ///@name  Preconditioners 
  ///@{ 
  
  /// ILUPreconditioner class. 
  /**   */
  template<class TSparseSpaceType, class TDenseSpaceType>
    class ILUPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of ILUPreconditioner
      typedef boost::shared_ptr<ILUPreconditioner> Pointer;


      typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;


      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ILUPreconditioner()
      { 
        L = NULL;
        iL = NULL;
        jL = NULL;
        U = NULL;
        iU = NULL;
        jU = NULL;
      }


      /// Copy constructor.
      ILUPreconditioner(const ILUPreconditioner& Other){}


      /// Destructor.
      virtual ~ILUPreconditioner()
      { 
        if ( L!=NULL) delete[]  L;
        if (iL!=NULL) delete[] iL;
        if (jL!=NULL) delete[] jL;
        if ( U!=NULL) delete[]  U;
        if (iU!=NULL) delete[] iU;
        if (jU!=NULL) delete[] jU;
        
        L = NULL;
	iL = NULL;
	jL = NULL;
	U = NULL;
	iU = NULL;
	jU = NULL;
      }
      


      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      ILUPreconditioner& operator=(const ILUPreconditioner& Other)
        {
          mILUSize = Other.mILUSize;
          unsigned int size = Other.iL[mILUSize];
          L = new double[size];
          U = new double[size];
          iL = new int[mILUSize+1];
          jL = new int[size];
          iU = new int[mILUSize+1];
          jU = new int[size];


          std::copy(Other.L, Other.L+size, L); 
          std::copy(Other.U, Other.U+size, U); 
          std::copy(Other.iL, Other.iL+mILUSize+1, iL); 
          std::copy(Other.jL, Other.jL+size, jL); 
          std::copy(Other.iU, Other.iU+mILUSize+1, iU); 
          std::copy(Other.jU, Other.jU+size, jU); 


          return *this;
        }


      
      ///@}
      ///@name Operations
      ///@{
      


      virtual void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
      {
        VectorType z = rX;
        TSparseSpaceType::Mult(rA,z, rY);
        ApplyLeft(rY);
      }
      
      virtual void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
      {
        VectorType z = rX;
        ApplyTransposeLeft(z);
        TSparseSpaceType::TransposeMult(rA,z, rY);
      }
      
      /** multiply first rX by L^-1 and store result in temp
          then multiply temp by U^-1 and store result in rX 
          @param rX  Unknows of preconditioner suystem
      */
      virtual VectorType& ApplyLeft(VectorType& rX)
      {
        const int size = TSparseSpaceType::Size(rX);
        VectorType temp(size);
        double sum;
        int i, indexj;
        for (i=0; i<size; i++) {
                sum=rX[i];
                for (indexj=iL[i]; indexj<iL[i+1]; indexj++) {
                        sum=sum-L[indexj]*temp[jL[indexj]];
                }
                temp[i]=sum;
        }
        for (i=size-1; i>=0; i--) {
                sum=temp[i];
                for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++) {
                        sum=sum-U[indexj]*rX[jU[indexj]];
                }
                rX[i]=sum/U[iU[i]];
        }
        return rX;
      }
      
      /** Multiply first rX by U^-T and store result in temp
          then multiply temp by L^-T and store result in rX
          @param rX  Unknows of preconditioner suystem
      */    
      virtual VectorType& ApplyTransposeLeft(VectorType& rX)
      {
        const int size = TSparseSpaceType::Size(rX);
        VectorType temp(size);
        int i, indexj;
        double tempi, rxi;
        for (i=0; i<size; i++) temp[i]=rX[i];
        for (i=0; i<size; i++) {
                temp[i]=temp[i]/U[iU[i]];
                tempi=temp[i];
                for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++) {
                        temp[jU[indexj]]=temp[jU[indexj]]-tempi*U[indexj];
                }
        }
        for (i=0; i<size; i++) rX[i]=temp[i];
        for (i=size-1; i>=0; i--) {
                rxi=rX[i];
                for (indexj=iL[i]; indexj<iL[i+1]; indexj++) {
                        rX[jL[indexj]]=rX[jL[indexj]]-rxi*L[indexj];
                }
        }
        return rX;
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
      
      /// Return information about this object.
	  virtual std::string Info() const
        {
          return "ILUPreconditioner";
        }


      /// Print information about this object.
      virtual void  PrintInfo(std::ostream& OStream) const
        {
          OStream << "ILUPreconditioner";
        }


      virtual void PrintData(std::ostream& OStream) const
      {
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
        
      unsigned int mILUSize;
      int *iL, *jL, *iU, *jU;
      double *L, *U;
        
        
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
        
    }; // Class ILUPreconditioner 
  
  ///@} 


  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
  inline std::istream& operator >> (std::istream& IStream, 
                                    ILUPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
	   return IStream;
   }


  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
  inline std::ostream& operator << (std::ostream& OStream, 
                                    const ILUPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);


      return OStream;
   }
  ///@} 
  
  
}  // namespace Kratos.


#endif // KRATOS_ILU_PRECONDITIONER_H_INCLUDED  defined 

