//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_FLAGS_H_INCLUDED )
#define  KRATOS_FLAGS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


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
  
  /// Short class definition.
  /** Detail class definition.
  */
  class Flags
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Flags
      KRATOS_CLASS_POINTER_DEFINITION(Flags);
      
      typedef int64_t BlockType;
      
      typedef std::size_t IndexType;
      
      enum FlagList
      {
	Flag0 = BlockType(1),
	Flag1 = BlockType(1) << 1,
	Flag2 = BlockType(1) << 2,
	Flag3 = BlockType(1) << 3,
	Flag4 = BlockType(1) << 4,
	Flag5 = BlockType(1) << 5,
	Flag6 = BlockType(1) << 6,
	Flag7 = BlockType(1) << 7,
	Flag8 = BlockType(1) << 8,
	Flag9 = BlockType(1) << 9,
	
 	Flag10 = BlockType(1) << 10,
	Flag11 = BlockType(1) << 11,
	Flag12 = BlockType(1) << 12,
	Flag13 = BlockType(1) << 13,
	Flag14 = BlockType(1) << 14,
	Flag15 = BlockType(1) << 15,
	Flag16 = BlockType(1) << 16,
	Flag17 = BlockType(1) << 17,
	Flag18 = BlockType(1) << 18,
	Flag19 = BlockType(1) << 19,

	Flag20 = BlockType(1) << 20,
	Flag21 = BlockType(1) << 21,
	Flag22 = BlockType(1) << 22,
	Flag23 = BlockType(1) << 23,
	Flag24 = BlockType(1) << 24,
	Flag25 = BlockType(1) << 25,
	Flag26 = BlockType(1) << 26,
	Flag27 = BlockType(1) << 27,
	Flag28 = BlockType(1) << 28,
	Flag29 = BlockType(1) << 29,

	Flag30 = BlockType(1) << 30,
	Flag31 = BlockType(1) << 31
      };
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Flags() : mFlags(BlockType()){}

      /// Copy constructor.
      Flags(Flags const& rOther) : mFlags(rOther.mFlags)
      {
      }
      
      template<class TFlagsType>
      Flags(TFlagsType const& rOther)
      {
	mFlags=static_cast<BlockType>(rOther);
      }

      /// Destructor.
      virtual ~Flags(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      Flags& operator=(Flags const& rOther)
      {
	mFlags = rOther.mFlags;
	return *this;
      }
      
      
      operator bool()
      {
	return mFlags;
      }
      
      Flags operator~() const
      {
	return  ~mFlags;
      }
      
      bool operator!() const
      {
	return  !mFlags;
      }
      
      ///@}
      ///@name Operations
      ///@{
      
	void Set(IndexType Position, bool Value=true )
	{
	  if(Value)
		mFlags |= (1 << Position);
	  else
		mFlags &= ~(1 << Position);
	}

	bool Get(IndexType Position) const
	{
		return (mFlags & (1 << Position));
	}
      
	void Clear(IndexType Position)
	{
		Set(Position,false);
	}
	
	void Clear()
	{
	  mFlags = BlockType();
	}
	
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
	bool Is(BlockType const& Flag)
	{
	  return (mFlags & Flag);
	}
      
	template<class TFlagsType>
	bool Is(TFlagsType const& Flag)
	{
	  return (mFlags & Flag);
	}
      
	bool Is(Flags const& rOther)
	{
	  return (mFlags & rOther.mFlags);
	}

	bool IsNot(BlockType const& Flag)
	{
	  return !(mFlags & Flag);
	}
      
	template<class TFlagsType>
	bool IsNot(TFlagsType const& Flag )
	{
	  return !(mFlags & Flag);
	}
      
	bool IsNot(Flags const& rOther)
	{
	  return !(mFlags & rOther.mFlags);
	}
      
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	std::stringstream buffer;
        buffer << "Flags" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const 
      {
	rOStream << "Flags";
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const 
      {
	
	for(std::size_t i = 0 ; i < sizeof(BlockType)*8 ; i++)
	  rOStream << bool(mFlags & (BlockType(1) << i));
      }
      
            
      ///@}      
      ///@name Friends
      ///@{
	
	
	friend bool operator==(const Flags& Left, const Flags& Right )
	{
	  return (Left.mFlags == Right.mFlags);
	}
	
	friend bool operator!=(const Flags& Left, const Flags& Right )
	{
	  return (Left.mFlags != Right.mFlags);
	}
	
	friend Flags operator|(const Flags& Left, const Flags& Right )
	{
	  return (Left.mFlags | Right.mFlags);
	}
	
	friend Flags operator&(const Flags& Left, const Flags& Right )
	{
	  return (Left.mFlags & Right.mFlags);
	}
	
	const Flags& operator|=(const Flags& Other )
	{
	  mFlags |= Other.mFlags;
	  return *this;
	}
	
	const Flags operator&=(const Flags& Other )
	{
	  mFlags &= Other.mFlags;
	  return *this;
	}
	
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
      
	BlockType mFlags;
        
        
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
        
    }; // Class Flags 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Flags& rThis)
  {
    return rIStream;
  }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Flags& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << " : ";
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FLAGS_H_INCLUDED  defined 


