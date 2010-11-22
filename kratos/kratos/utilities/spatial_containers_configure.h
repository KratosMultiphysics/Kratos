//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine  $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SPATIAL_CONTAINERS_CONFIGURE_INCLUDED )
#define  KRATOS_SPATIAL_CONTAINERS_CONFIGURE_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>



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
  

template <class T> 
class ContactPair 
{
   public:
   T value[2];
   
  ContactPair(){}
   
   ContactPair(const T& First,const T& Second)
       {
	 value[0]  = First;
	 value[1]  = Second;
       }
       
    ~ContactPair(){}
     T& operator[](std::size_t index)
             { 
		return value[index];
	     }

    
    T const& operator[](std::size_t index) const 
    {
         return value[index];
    }
    
    
    ContactPair& operator = (ContactPair& Pair)
             { 
		value[0] = Pair[0];
		value[1] = Pair[1];
		return *this; 
	     }
   
    
    inline bool operator == (const ContactPair& Pair)
             { 
	       return  (value[0] == Pair[0]) && (value[1] == Pair[1]) ;
	     }
	     
    inline bool  operator == (const ContactPair& Pair) const
             { 
                return ( (value[0] == Pair[0]) && (value[1] == Pair[1]) ) ;
	     }
	     
    inline std::ostream& operator << (const ContactPair& Pair)
         {
	    std::ostream OStream;
	    OStream << "Object 1 =  "<<  Pair[0] << "  " << "Object 2 =" << Pair[1] << std::endl;  
	    return OStream;  
	 }
	 
    inline std::ostream& operator << (const ContactPair& Pair) const
         {
	    std::ostream OStream;
	    OStream << "Object 1 =  "<<  Pair[0] << "  " << "Object 2 =" << Pair[1] << std::endl; 
	    return OStream;
	 }	 
    
    
};  
 
template <std::size_t TDimension>
class SpatialContainersConfigure
    {
    public:
      ///@name Type Definitions
      ///@{
         
      enum { Dimension = TDimension };
      typedef Point<Dimension, double>                        PointType;
      typedef std::vector<double>::iterator                   DistanceIteratorType;
      typedef ModelPart::ElementsContainerType::ContainerType ContainerType;
      typedef ContainerType::value_type                       PointerType;
      typedef ContainerType::iterator                         IteratorType; 
      typedef ModelPart::ElementsContainerType::ContainerType ResultContainerType;    
      typedef ContainerType::iterator                         ResultIteratorType; 
      
      /// Contact Pairs
//       typedef std::pair<PointerType, PointerType>            ContactPairType;
//       typedef std::vector<ContactPairType>                   ContainerContactType; 
//       typedef std::vector<ContactPairType>::iterator         IteraratorContactType; 
      
      /// Contact Pairs
      typedef ContactPair<PointerType>                       ContactPairType;
      //typedef array_1d<PointerType, 2>                       ContactPairType;
      typedef std::vector<ContactPairType>                   ContainerContactType; 
      typedef std::vector<ContactPairType>::iterator         IteratorContactType; 
      
      
        
      /// Pointer definition of SpatialContainersConfigure
      KRATOS_CLASS_POINTER_DEFINITION(SpatialContainersConfigure);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      SpatialContainersConfigure(){}

      /// Destructor.
      virtual ~SpatialContainersConfigure(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
         
///******************************************************************************************************************
///******************************************************************************************************************  
	
   static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    { 
    rHighPoint = rObject->GetGeometry().GetPoint(0);
    rLowPoint  = rObject->GetGeometry().GetPoint(0);        
     for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)      
       {
 	for(std::size_t i = 0; i<TDimension; i++)
	  {
 	      rLowPoint[i]  =  (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i]; 
 	      rHighPoint[i] =  (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
 	  }
        }           
    }
    
///******************************************************************************************************************
///******************************************************************************************************************   

     static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
      { 
      	      Element::GeometryType& geom_1 = rObj_1->GetGeometry();
	      Element::GeometryType& geom_2 = rObj_2->GetGeometry();
	      return  geom_1.HasIntersection(geom_2); 
      
      }


///******************************************************************************************************************
///******************************************************************************************************************   
    
      static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
      { 
	 return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint); 
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
      virtual std::string Info() const {return " Spatial Containers Configure"; }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {}

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
      SpatialContainersConfigure& operator=(SpatialContainersConfigure const& rOther);

      /// Copy constructor.
      SpatialContainersConfigure(SpatialContainersConfigure const& rOther);

        
      ///@}    
        
    }; // Class SpatialContainersConfigure 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template <std::size_t TDimension>
  inline std::istream& operator >> (std::istream& rIStream, 
				    SpatialContainersConfigure<TDimension> & rThis){return rIStream; }

  /// output stream function
  template <std::size_t TDimension>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const SpatialContainersConfigure<TDimension>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_SPATIAL_CONTAINERS_CONFIGURE_INCLUDED  defined 


