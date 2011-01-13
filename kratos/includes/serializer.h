//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SERIALIZER_H_INCLUDED )
#define  KRATOS_SERIALIZER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <map>
#include <set>


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/buffer.h"
#include "containers/weak_pointer_vector.h"
// #include "containers/variable.h"


#define KRATOS_SERIALIZATION_DIRECT_LOAD(type)                           \
	void load(std::string const & rTag, type& rValue)                \
	{                                                                \
	  mBuffer >> rValue;                                             \
	}  								 \
	   								 \
	void load_base(std::string const & rTag, type& rValue)           \
	{                                                                \
	  mBuffer >> rValue;                                             \
	}  

#define KRATOS_SERIALIZATION_DIRECT_SAVE(type)                           \
	void save(std::string const & rTag, type const & rValue)         \
	{                                                                \
	  mBuffer << rValue;                                             \
	}    								 \
	   								 \
	void save_base(std::string const & rTag, type const & rValue)    \
	{                                                                \
	  mBuffer << rValue;                                             \
	}  

#define KRATOS_SERIALIZATION_DIRECT_CREATE(type)                         \
	void* create(std::string const & rTag, type* prototype)          \
	{                                                                \
	  type* p_new = new type;                                        \
	  load(rTag, *p_new);                                            \
	  return p_new;                                                  \
	}    								 
	
namespace Kratos
{

  class ModelPart;
  class VariableData;
  template <class TDataType> class Variable;
  template <class TDataType> class KratosComponents;


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
  class Serializer
    {
    public:
      ///@name  Enum's
      ///@{
      
      enum PointerType{SP_INVALID_POINTER, SP_BASE_CLASS_POINTER, SP_DERIVED_CLASS_POINTER};

      ///@}
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Serializer
      KRATOS_CLASS_POINTER_DEFINITION(Serializer);
  
      typedef std::size_t SizeType;
      
      typedef void* (*ObjectFactoryType)();
      
      typedef std::map<void*, void*> LoadedPointersContainerType;
            
      typedef std::map<std::string, ObjectFactoryType> RegisteredObjectsContainerType;

      typedef std::set<const void*> SavedPointersContainerType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Serializer(SizeType NewBufferSize = 1024) : mBuffer(NewBufferSize)
      {
      }
      

      /// Destructor.
      virtual ~Serializer()
      {
      }
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
	template<class TDataType>
	static void* Create()
	{
	  return new TDataType;
	}
      
	template<class TDataType>
	static void Register(std::string const & rName, TDataType const& pPrototype)
	{
	  msRegisteredObjects.insert(RegisteredObjectsContainerType::value_type(rName,Create<TDataType>));
// 	  msRegisteredObjects.insert(RegisteredObjectsContainerType::value_type(rName,&pPrototype));
	}
      
	template<class TDataType>
	void load(std::string const & rTag, TDataType& rObject)
	{
	  rObject.load(*this);
	}
      
      template<class TDataType>
      void load(std::string const & rTag, boost::shared_ptr<TDataType>& pValue)
      {
	PointerType pointer_type = SP_INVALID_POINTER;
	void* p_pointer;
	mBuffer >> pointer_type;
	
	if(pointer_type != SP_INVALID_POINTER)
	{
	  mBuffer >> p_pointer;
	  LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
	  if(i_pointer == mLoadedPointers.end())
	  {
	    if(pointer_type == SP_BASE_CLASS_POINTER)
	    {
	      if(!pValue)
		pValue = boost::shared_ptr<TDataType>(new TDataType);
	      
	      load(rTag, *pValue);
	    }
	    else if(pointer_type == SP_DERIVED_CLASS_POINTER)
	    {
	      std::string object_name;
	      mBuffer >> object_name;
	      typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);
	      
	      if(i_prototype == msRegisteredObjects.end())
		KRATOS_ERROR(std::runtime_error, "There is no object registered in Kratos with name : ", object_name)
	      
	      if(!pValue)
		pValue = boost::shared_ptr<TDataType>(static_cast<TDataType*>((i_prototype->second)()));

	      pValue->load(*this);
	      
	    }
	    mLoadedPointers[p_pointer]=&pValue;
	  }
	  else
	    pValue = *static_cast<boost::shared_ptr<TDataType>*>((i_pointer->second));
	}
      }
      
      template<class TDataType>
      void load(std::string const & rTag, TDataType*& pValue)
      {
	PointerType pointer_type = SP_INVALID_POINTER;
	void* p_pointer;
	mBuffer >> pointer_type;
	
	if(pointer_type != SP_INVALID_POINTER)
	{
	  mBuffer >> p_pointer;
	  LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
	  if(i_pointer == mLoadedPointers.end())
	  {
	    if(pointer_type == SP_BASE_CLASS_POINTER)
	    {
	      if(!pValue)
		pValue = new TDataType;
	      
	      load(rTag, *pValue);
	    }
	    else if(pointer_type == SP_DERIVED_CLASS_POINTER)
	    {
	      if(!pValue)
		pValue = new TDataType;
	      
	      load(rTag, *pValue);
	      
	    }
	    mLoadedPointers[p_pointer]=&pValue;
	  }
	  else 
	  {
	    pValue = *static_cast<TDataType**>((i_pointer->second));
	  }
	}
      }
      
       template<class TDataType>
	void load(std::string const & rTag, boost::weak_ptr<TDataType>& pValue)
	{
	  // This is for testing. I have to change it. Pooyan.
	  KRATOS_ERROR(std::logic_error, "The serialization for weak_ptrs is not implemented yet", "")
// 	  mBuffer >> *pValue;
	}
	
	template<class TDataType>
	void load(std::string const & rTag, WeakPointerVector<TDataType>& pValue)
	{
	  // This is for testing. I have to change it. Pooyan.
	  KRATOS_ERROR(std::logic_error, "The serialization for weak_ptrs is not implemented yet", "")
// 	  mBuffer >> *pValue;
	}

	template<class TDataType>
	void load(std::string const & rTag, const Variable<TDataType>* pVariable);
	
 	template<class TDataType>
	void load(std::string const & rTag, std::vector<TDataType>& rObject)
	{
	  mBuffer >> rObject;
	}
	
 	template<class TDataType>
	void load(std::string const & rTag, boost::numeric::ublas::vector<TDataType>& rObject)
	{
	  mBuffer >> rObject;
	}
	
 	template<class TDataType, std::size_t TDimension>
	void load(std::string const & rTag, array_1d<TDataType, TDimension>& rObject)
	{
	  mBuffer >> rObject;
	}
     
        KRATOS_SERIALIZATION_DIRECT_LOAD(bool)
        KRATOS_SERIALIZATION_DIRECT_LOAD(int)
        KRATOS_SERIALIZATION_DIRECT_LOAD(long)
        KRATOS_SERIALIZATION_DIRECT_LOAD(double)
        KRATOS_SERIALIZATION_DIRECT_LOAD(std::size_t)
        KRATOS_SERIALIZATION_DIRECT_LOAD(unsigned int)
        KRATOS_SERIALIZATION_DIRECT_LOAD(std::string)
        KRATOS_SERIALIZATION_DIRECT_LOAD(Matrix)

	

	template<class TDataType>
	void save(std::string const & rTag, std::vector<TDataType> const& rObject)
	{
	  mBuffer << rObject;
	}
      
	template<class TDataType>
	void save(std::string const & rTag, boost::numeric::ublas::vector<TDataType> const& rObject)
	{
	  mBuffer << rObject;
	}
      
	template<class TDataType, std::size_t TDimension>
	void save(std::string const & rTag, array_1d<TDataType, TDimension> const& rObject)
	{
	  mBuffer << rObject;
	}
      
	template<class TDataType>
	void save(std::string const & rTag, TDataType const& rObject)
	{
	  rObject.save(*this);
	}
      
	template<class TDataType>
	void save(std::string const & rTag, const Variable<TDataType>* pVariable)
	{
	  mBuffer << pVariable->Name();
	}
	
     
	template<class TDataType>
	void save(std::string const & rTag, boost::shared_ptr<TDataType> pValue)
	{
	  save(rTag, pValue.get());
	}
     
	template<class TDataType>
	void save(std::string const & rTag, const TDataType * pValue)
	{
	  if(pValue)
	  {
	    if(IsDerived(pValue))
	      mBuffer << SP_DERIVED_CLASS_POINTER;
	    else
	      mBuffer << SP_BASE_CLASS_POINTER;
	    
	    SavePointer(rTag,pValue);	      
	  }
	  else
	  {
	    mBuffer << SP_INVALID_POINTER;
	  }
	}
	
	template<class TDataType>
	bool IsDerived(TDataType * pValue)
	{
	  bool is_derived = (typeid(TDataType) != typeid(*pValue));
// 	  std::cout << "for TDataType : " << typeid(TDataType).name() << " and *pValue type : " << typeid(*pValue).name() << " is derived : " << is_derived << std::endl;
	  return is_derived;
	}
	
     
	template<class TDataType>
	void save(std::string const & rTag, TDataType * pValue)
	{
	  if(pValue)
	  {
	    if(IsDerived(pValue))
	    {
	      mBuffer << SP_DERIVED_CLASS_POINTER;
	    }
	    else
	    {
	      mBuffer << SP_BASE_CLASS_POINTER;
	    }
	    
	    SavePointer(rTag,pValue);	      
	  }
	  else
	  {
	    mBuffer << SP_INVALID_POINTER;
	  }
	}
     
	template<class TDataType>
	void save(std::string const & rTag, boost::weak_ptr<TDataType> pValue)
	{
	  // This is for testing. I have to implement it. Pooyan.
	  KRATOS_ERROR(std::logic_error, "The serialization for weak_ptrs is not implemented yet", "")
// 	  mBuffer << *pValue;
	}
     
	template<class TDataType>
	void save(std::string const & rTag, Kratos::WeakPointerVector<TDataType> pValue)
	{
	  // This is for testing. I have to implement it. Pooyan.
	  KRATOS_ERROR(std::logic_error, "The serialization for weak_ptrs is not implemented yet", "")
// 	  mBuffer << *pValue;
	}
     
	template<class TDataType>
	void save(std::string const & rTag, boost::shared_ptr<const TDataType> pValue)
	{
	  // This is for testing. I have to change it. Pooyan.
	  mBuffer << *pValue;
	}

	void save(std::string const & rTag, const char * pValue)
	{
	  mBuffer << std::string(pValue);
	}

	KRATOS_SERIALIZATION_DIRECT_SAVE(bool)
	KRATOS_SERIALIZATION_DIRECT_SAVE(int)
        KRATOS_SERIALIZATION_DIRECT_SAVE(long)
        KRATOS_SERIALIZATION_DIRECT_SAVE(double)
        KRATOS_SERIALIZATION_DIRECT_SAVE(std::size_t)
        KRATOS_SERIALIZATION_DIRECT_SAVE(unsigned int)
        KRATOS_SERIALIZATION_DIRECT_SAVE(std::string)
        KRATOS_SERIALIZATION_DIRECT_SAVE(Vector)
        KRATOS_SERIALIZATION_DIRECT_SAVE(Matrix)
       
      
	template<class TDataType>
	void load_base(std::string const & rTag, TDataType& rObject)
	{
	  rObject.TDataType::load(*this);
	}
 
	
 	template<class TDataType>
	void load_base(std::string const & rTag, std::vector<TDataType>& rObject)
	{
	  mBuffer >> rObject;
	}
	
 	template<class TDataType>
	void load_base(std::string const & rTag, boost::numeric::ublas::vector<TDataType>& rObject)
	{
	  mBuffer >> rObject;
	}
	
 	template<class TDataType, std::size_t TDimension>
	void load_base(std::string const & rTag, array_1d<TDataType, TDimension>& rObject)
	{
	  mBuffer >> rObject;
	}

	template<class TDataType>
	void save_base(std::string const & rTag, std::vector<TDataType> const& rObject)
	{
	  mBuffer << rObject;
	}
      
	template<class TDataType>
	void save_base(std::string const & rTag, boost::numeric::ublas::vector<TDataType> const& rObject)
	{
	  mBuffer << rObject;
	}
      
	template<class TDataType, std::size_t TDimension>
	void save_base(std::string const & rTag, array_1d<TDataType, TDimension> const& rObject)
	{
	  mBuffer << rObject;
	}
      
	template<class TDataType>
	void save_base(std::string const & rTag, TDataType const& rObject)
	{
	  rObject.TDataType::save(*this);
	}
	
	KRATOS_SERIALIZATION_DIRECT_CREATE(bool)
	KRATOS_SERIALIZATION_DIRECT_CREATE(int)
        KRATOS_SERIALIZATION_DIRECT_CREATE(long)
        KRATOS_SERIALIZATION_DIRECT_CREATE(double)
        KRATOS_SERIALIZATION_DIRECT_CREATE(std::size_t)
        KRATOS_SERIALIZATION_DIRECT_CREATE(unsigned int)
        KRATOS_SERIALIZATION_DIRECT_CREATE(std::string)
        KRATOS_SERIALIZATION_DIRECT_CREATE(Vector)
        KRATOS_SERIALIZATION_DIRECT_CREATE(Matrix)

    
        
      ///@}
      ///@name Access
      ///@{ 
      
	Buffer& GetBuffer()
	{
	  return mBuffer;
	}
	
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	return "Serializer";
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {}
      
            
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
      
	static RegisteredObjectsContainerType msRegisteredObjects;
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
      
	Buffer mBuffer;
	
	SavedPointersContainerType mSavedPointers;
	LoadedPointersContainerType mLoadedPointers;
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
     
	template<class TDataType>
	void SavePointer(std::string const & rTag, const TDataType * pValue)
	{
	    mBuffer << pValue;
	    if(mSavedPointers.find(pValue) == mSavedPointers.end())
	    {
	      save(rTag,*pValue);
// 	      pValue->save(*this);
	      mSavedPointers.insert(pValue);
	    }
	}
        
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
      Serializer& operator=(Serializer const& rOther);

      /// Copy constructor.
      Serializer(Serializer const& rOther);

        
      ///@}    
        
    }; // Class Serializer 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
  
//   template<class TDataType>
//   inline Serializer& operator >> (Serializer& rThis, TDataType& rObject)
//   {
//     rThis.load(rObject);
    
//     return rThis;
//   }

  
//   template<class TDataType>
//   inline Serializer& operator << (Serializer& rThis, TDataType& rObject)
//   {
//     rThis.save(rObject, KRATOS_VERSION);
    
//     return rThis;
//   }
  /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
// 				    Serializer& rThis);

  /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream, 
// 				    const Serializer& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
  ///@} 
  
  
}  // namespace Kratos.

#undef KRATOS_SERIALIZATION_DIRECT_LOAD
#undef KRATOS_SERIALIZATION_DIRECT_SAVE

#endif // KRATOS_SERIALIZER_H_INCLUDED  defined 


