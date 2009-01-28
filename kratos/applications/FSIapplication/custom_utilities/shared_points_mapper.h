//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:42 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SHARED_POINTS_MAPPER_H_INCLUDED )
#define  KRATOS_SHARED_POINTS_MAPPER_H_INCLUDED



// System includes
#include <set> 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "containers/pointer_vector.h"
#include "includes/node.h"


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
		class SharedPointsMapper 
    {
    public:
		///@name Type Definitions
		///@{
		
		/// Counted pointer of SharedPointsMapper
		KRATOS_CLASS_POINTER_DEFINITION(SharedPointsMapper);
		
		///@}
		///@name Life Cycle 
		///@{ 
		
		/// Constructor with given array of Nodes.
		//**************************************************************************************************
		//**************************************************************************************************
		SharedPointsMapper(
			const ModelPart::NodesContainerType& OriginNodes, 
			const ModelPart::NodesContainerType& DestinationNodes,
            double tol = 1e-9)
		{
		KRATOS_TRY

		if( OriginNodes.size()!=0 && DestinationNodes.size()!=0)
		{	
			if( OriginNodes.size() != DestinationNodes.size() )
				KRATOS_ERROR(std::logic_error,"wrong number of nodes","");

			mOriginNodes.reserve( OriginNodes.size() );
			mDestinationNodes.reserve( DestinationNodes.size() );
			
			for(ModelPart::NodesContainerType::const_iterator origin_it = OriginNodes.begin(); origin_it != OriginNodes.end(); origin_it++)
			{
				for(ModelPart::NodesContainerType::const_iterator destination_it = DestinationNodes.begin(); destination_it != DestinationNodes.end(); destination_it++)
				{
						if	( 
							fabs(origin_it->X() - destination_it->X() ) < tol  &&
							fabs(origin_it->Y() - destination_it->Y() ) < tol  &&
							fabs(origin_it->Z() - destination_it->Z() ) < tol 				 
							)
						{							
							mOriginNodes.push_back(  *(origin_it.base() ) );
							mDestinationNodes.push_back(  *(destination_it.base() ) );
						}		              
				}
	            
				}
		}

		KRATOS_CATCH("")
         }
		
	
		/// Destructor.
		virtual ~SharedPointsMapper(){}
		
		
		///@}
		///@name Operators 
		///@{
		
		
		///@}
		///@name Operations
		///@{
		//**************************************************************************************************
		//**************************************************************************************************
		void ScalarMap( const Variable<double>& rOriginVariable, 
				const Variable<double>& rDestinationVariable)		
		{
			KRATOS_TRY
					
			PointerVector< Node<3> >::iterator it_origin = mOriginNodes.begin();
   			PointerVector< Node<3> >::iterator it_destination = mDestinationNodes.begin();
			
			for(unsigned int i = 0 ; i < mOriginNodes.size() ; ++i)
			{
			   (it_destination++ )->FastGetSolutionStepValue(rDestinationVariable) = 
			   			   (it_origin++ )->FastGetSolutionStepValue(rOriginVariable);
			   
			}

			KRATOS_CATCH("")			
		}

		//**************************************************************************************************
		//**************************************************************************************************
		void VectorMap( const Variable<array_1d<double,3> >& rOriginVariable, 
				const Variable<array_1d<double,3> >& rDestinationVariable)		
		{
			KRATOS_TRY
					
			PointerVector< Node<3> >::iterator it_origin = mOriginNodes.begin();
   			PointerVector< Node<3> >::iterator it_destination = mDestinationNodes.begin();
			
			for(unsigned int i = 0 ; i < mOriginNodes.size() ; ++i)
			{
			   noalias( (it_destination++ )->FastGetSolutionStepValue(rDestinationVariable) ) = 
			   			   (it_origin++ )->FastGetSolutionStepValue(rOriginVariable);
			   
			}

			KRATOS_CATCH("")			
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
		
		/// Print information about this object.
		virtual void PrintInfo(std::ostream& OStream) const
		{
		}
		
		
		/// Print object's data.
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
		PointerVector< Node<3> > mOriginNodes;
        
		PointerVector< Node<3> > mDestinationNodes;
		        
        
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
		SharedPointsMapper& operator=(const SharedPointsMapper& Other);
		
		/// Copy constructor.
		//SharedPointsMapper(const SharedPointsMapper& Other);
		
        
		///@}    
        
    }; // Class SharedPointsMapper 
	
	///@} 
	
	///@name Type Definitions       
	///@{ 
	
	
	///@} 
	///@name Input and output 
	///@{ 
	
	///@} 
	
	
}  // namespace Kratos.

#endif // KRATOS_SHARED_POINTS_MAPPER_H_INCLUDED  defined 


