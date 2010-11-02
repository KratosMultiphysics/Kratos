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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-07-22 11:16:11 $
//   Revision:            $Revision: 1.17 $
//
//


#if !defined(KRATOS_NODE_H_INCLUDED )
#define  KRATOS_NODE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/point.h"
#include "includes/dof.h"
//#include "includes/neighbours.h"
#include "containers/buffer.h"
#include "containers/vector_map.h"
#include "containers/pointer_vector_set.h"
//#include "containers/fix_data_value_container.h"
#include "containers/variables_list_data_value_container.h"
#include "utilities/indexed_object.h"

#include "containers/weak_pointer_vector.h"

#ifdef _OPENMP
 #include "omp.h"
#endif


namespace Kratos
{

  class Element;

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
  template<std::size_t TDimension, class TDofType = Dof<double> >
    class Node : public Point<TDimension>,  public IndexedObject
    {
      class GetDofKey : public std::unary_function<TDofType, VariableData::KeyType> {
      public:
	VariableData::KeyType operator()(TDofType  const & This){return This.GetVariable().Key();}
      };

    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Node
      KRATOS_CLASS_POINTER_DEFINITION(Node);
  
      typedef Node<TDimension, TDofType> NodeType;

      typedef Point<TDimension> BaseType;

      typedef Point<TDimension> PointType;
      
      typedef TDofType DofType;

      typedef std::size_t IndexType;
      
      typedef typename std::size_t SizeType;

      //typedef Neighbours<NodeType, Element> NeighboursType;

      /** Element array type for neighbor elements. */
      //typedef typename NeighboursType::NeighbourElementsArrayType NeighbourElementsArrayType;

      /** Element array type for neighbor elements. */
      //typedef typename NeighboursType::NeighbourNodesArrayType NeighbourNodesArrayType;

      typedef PointerVectorSet<TDofType, GetDofKey> DofsContainerType;
      
      //typedef Buffer<std::vector<FixDataValueContainer> > SolutionStepsNodalDataContainerType;
      
//        typedef Buffer<std::vector<VariablesListDataValueContainer> > SolutionStepsNodalDataContainerType;
       typedef VariablesListDataValueContainer SolutionStepsNodalDataContainerType; 

      typedef VariablesListDataValueContainer::BlockType BlockType;

      //typedef VectorMap<IndexType, NeighboursType > SolutionStepsNeighboursContainerType;
      
      typedef Variable<double> DoubleVariableType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Node() 
	: BaseType()
	, IndexedObject(0)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition()
	{
		CreateSolutionStepData();

                #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
	}
	
      Node(IndexType NewId ) 
	: BaseType()
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition()
	{
		KRATOS_ERROR(std::logic_error, "calling the default constructor for the node ... illegal operation!!","");
	  CreateSolutionStepData();

          #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
	}

      /// 1d constructor.
      Node(IndexType NewId, double const& NewX) 
	: BaseType(NewX)
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition(NewX)
      {
	  CreateSolutionStepData();

          #ifdef _OPENMP
              omp_init_lock(&mnode_lock);
          #endif
      }
      
      /// 2d constructor.
      Node(IndexType NewId, double const& NewX, double const& NewY) 
	: BaseType(NewX, NewY)
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition(NewX, NewY)
      {
	  CreateSolutionStepData();

            #ifdef _OPENMP
               omp_init_lock(&mnode_lock);
            #endif
      }
      
      /// 3d constructor.
      Node(IndexType NewId, double const& NewX, double const& NewY, double const& NewZ) 
	: BaseType(NewX, NewY, NewZ)
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition(NewX, NewY, NewZ)
      {
	  CreateSolutionStepData();

          #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
      }
      
      /// Point constructor.
      Node(IndexType NewId, PointType const& rThisPoint) 
	: BaseType(rThisPoint)
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition(rThisPoint)
      {
	  CreateSolutionStepData();

          #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
      }
      
      /** Copy constructor. Initialize this node with given node.*/
      Node(Node const& rOtherNode)
	: BaseType(rOtherNode)
	, IndexedObject(rOtherNode)
	, mDofs(rOtherNode.mDofs)
        , mData(rOtherNode.mData)
	, mSolutionStepsNodalData(rOtherNode.mSolutionStepsNodalData)
	, mInitialPosition(rOtherNode.mInitialPosition)
	{
		//TODO ... this copy constructor should be removed sometimes as it is often source of error
		//KRATOS_ERROR(std::logic_error, "copying Nodes is not allowed", "");

          #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
	}
      
      /** Copy constructor from a node with different dimension.*/
 template<SizeType TOtherDimension>
      Node(Node<TOtherDimension> const& rOtherNode) 
	: BaseType(rOtherNode)
	, IndexedObject(rOtherNode)
	, mDofs(rOtherNode.mDofs)
	, mData(rOtherNode.mData)
	, mSolutionStepsNodalData(rOtherNode.mSolutionStepsNodalData)
	, mInitialPosition(rOtherNode.mInitialPosition)
      {

     #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
      }
      
      /** Copy constructor from a point with different dimension.*/
      template<SizeType TOtherDimension>
      Node(IndexType NewId, Point<TOtherDimension> const& rThisPoint) 
	: BaseType(rThisPoint)
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition(rThisPoint)
      {
	  CreateSolutionStepData();

          #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
      }
      
        /** Constructor using coordinates stored in given array. Initialize
        this point with the coordinates in the array. */
      template<class TVectorType>
      Node(IndexType NewId, vector_expression<TVectorType> const&  rOtherCoordinates)
	: BaseType(rOtherCoordinates)
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition(rOtherCoordinates)
      {
	  CreateSolutionStepData();

          #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
      }
      

        
        /** Constructor using coordinates stored in given std::vector. Initialize
        this point with the coordinates in the array. */
      Node(IndexType NewId, std::vector<double> const&  rOtherCoordinates)
	: BaseType(rOtherCoordinates)
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData()
	, mInitialPosition()
      {
	  CreateSolutionStepData();

          #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
      }

       /// 3d with variables list and data constructor.
      Node(IndexType NewId, double const& NewX, double const& NewY, double const& NewZ, VariablesList*  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1) 
	: BaseType(NewX, NewY, NewZ)
	, IndexedObject(NewId)
	, mDofs()
	, mData()
	, mSolutionStepsNodalData(pVariablesList,ThisData,NewQueueSize)
	, mInitialPosition(NewX, NewY, NewZ)
      {
          #ifdef _OPENMP
                    omp_init_lock(&mnode_lock);
                #endif
      }
      
       
      /// Destructor.
      virtual ~Node(){
          #ifdef _OPENMP
            omp_destroy_lock(&mnode_lock);
          #endif
	}

      void SetId(IndexType NewId)
      {
            IndexedObject::SetId(NewId);
            Node<3>::DofsContainerType& my_dofs = (this)->GetDofs();
            for(Node<3>::DofsContainerType::iterator iii = my_dofs.begin();    iii != my_dofs.end(); iii++)
            {
                iii->SetId(NewId);
            }
      }

      #ifdef _OPENMP
      omp_lock_t& GetLock(){return mnode_lock;}
      #endif

      inline void SetLock()
      {
           //does nothing if openMP is not present
           #ifdef _OPENMP
          omp_set_lock(&mnode_lock);
          #endif
      }

      inline void UnSetLock()
      {
           //does nothing if openMP is not present
           #ifdef _OPENMP
          omp_unset_lock(&mnode_lock);
          #endif
      }

      ///@}
      ///@name Operators 
      ///@{

      Node& operator=(const Node& rOther)
	{
	  BaseType::operator=(rOther);
	  
	  mDofs = rOther.mDofs;
	  mData = rOther.mData;
	  mSolutionStepsNodalData = rOther.mSolutionStepsNodalData;
	  mInitialPosition = rOther.mInitialPosition;

	  return *this;
	}

      /// Assignment operator.
      template<SizeType TOtherDimension>
	Node& operator=(const Node<TOtherDimension>& rOther)
	{
	  BaseType::operator=(rOther);
	  IndexedObject::operator=(rOther);
	  mDofs = rOther.mDofs;
	  mData = rOther.mData;
	  mSolutionStepsNodalData = rOther.mSolutionStepsNodalData;
	  mInitialPosition = rOther.mInitialPosition;

	  return *this;
	
      }
      
      template<class TVariableType> typename TVariableType::Type& operator()(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
	{
	  return GetSolutionStepValue(rThisVariable, SolutionStepIndex);
	}
      
      template<class TVariableType> typename TVariableType::Type& operator()(const TVariableType& rThisVariable)
	{
	  return GetSolutionStepValue(rThisVariable);
	}
      
      template<class TDataType> TDataType& operator[](const Variable<TDataType>& rThisVariable)
	{
	  return GetValue(rThisVariable);
	}
      
      template<class TDataType> const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
	{
	  return GetValue(rThisVariable);
	}
      
      template<class TAdaptorType> typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable)
	{
	  return GetValue(rThisVariable);
	}
      
      template<class TAdaptorType> const typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return GetValue(rThisVariable);
	}
       
      double& operator[](IndexType ThisIndex)
	{
	  return BaseType::operator[](ThisIndex);
	}
      
      double operator[](IndexType ThisIndex) const
	{
	  return BaseType::operator[](ThisIndex);
	}
      
       
      
      ///@}
      ///@name Nodal Data
      ///@{
      
//       void CreateSolutionStepData(IndexType SolutionStepIndex = IndexType())
// 	{
// 	  if(mSolutionStepsNodalData->find(SolutionStepIndex) == mSolutionStepsNodalData->end())
// 	    mSolutionStepsNodalData->insert(SolutionStepIndex, DataValueContainer());

// 	}

//       void CloneSolutionStepNodalData(IndexType SolutionStepIndex, IndexType SourceSolutionStepIndex)
// 	{
// 	  SolutionStepsNodalDataContainerType::iterator i;
// 	  if(mSolutionStepsNodalData->find(SolutionStepIndex) == mSolutionStepsNodalData->end())
// 	    {
// 	      if((i = mSolutionStepsNodalData->find(SourceSolutionStepIndex)) == mSolutionStepsNodalData->end())
// 		mSolutionStepsNodalData->insert(SolutionStepIndex, DataValueContainer());
// 	      else
// 		mSolutionStepsNodalData->insert(SolutionStepIndex, *i);
// 	    }
// 	}

//       void CloneSolutionStepNodalData(IndexType SolutionStepIndex, DataValueContainer const &  SourceSolutionStepData)
// 	{
// 	  if(mSolutionStepsNodalData->find(SolutionStepIndex) == mSolutionStepsNodalData->end())
// 		mSolutionStepsNodalData->insert(SolutionStepIndex, SourceSolutionStepData);
// 	}

      void CreateSolutionStepData()
	{
	    mSolutionStepsNodalData.PushFront();
// 	  VariablesListDataValueContainer temp(&Globals::DefaultVariablesList);
// 	  if(!mSolutionStepsNodalData.empty())
// 	      mSolutionStepsNodalData.push_front(temp);
// 	  else
// 	      mSolutionStepsNodalData.resize(1,temp);
// 	  mSolutionStepsNodalData->push_front(DataValueContainer());

	}

      void CloneSolutionStepData()
	{
	    mSolutionStepsNodalData.CloneFront();
// 	  if(!mSolutionStepsNodalData.empty())
// 	      mSolutionStepsNodalData.push_front(mSolutionStepsNodalData[SourceSolutionStepIndex]);
	}

      void OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
	{
	    mSolutionStepsNodalData.AssignData(mSolutionStepsNodalData.Data(SourceSolutionStepIndex), DestinationSourceSolutionStepIndex);
	}

//       void CloneSolutionStepNodalData(VariablesListDataValueContainer const &  SourceSolutionStepData)
// 	{
// 	  if(!mSolutionStepsNodalData.empty())
//   	      mSolutionStepsNodalData.push_front(SourceSolutionStepData);
// 	}

//       void OverwriteSolutionStepNodalData(IndexType SolutionStepIndex, IndexType SourceSolutionStepIndex, IndexType OldSolutionStepIndex)
// 	{
// 	  SolutionStepsNodalDataContainerType::iterator i;
// 	  SolutionStepsNodalDataContainerType::iterator i_old;
// 	  if(mSolutionStepsNodalData.find(SolutionStepIndex) == mSolutionStepsNodalData.end())
// 	    {
// 	      if((i_old = mSolutionStepsNodalData.find(OldSolutionStepIndex)) == mSolutionStepsNodalData.end())
// 		{
// 		  if((i = mSolutionStepsNodalData.find(SourceSolutionStepIndex)) == mSolutionStepsNodalData.end())
// 		    mSolutionStepsNodalData.insert(SolutionStepIndex, DataValueContainer());
// 		  else
// 		    mSolutionStepsNodalData.insert(SolutionStepIndex, *i);
// 		}
// 	      else
// 		{
// 		  i_old.base()->first = SolutionStepIndex;
// 		  i_old.base()->second = *i;
// 		}
// 	    }
// 	}

//       void OverwriteSolutionStepNodalData(IndexType SolutionStepIndex, DataValueContainer const &  SourceSolutionStepData)
// 	{
// 	  if(mSolutionStepsNodalData.find(SolutionStepIndex) == mSolutionStepsNodalData.end())
// 	      if((i_old = mSolutionStepsNodalData.find(OldSolutionStepIndex)) == mSolutionStepsNodalData.end())
// 		mSolutionStepsNodalData.insert(SolutionStepIndex, SourceSolutionStepData);
// 	      else
// 		{
// 		  i_old.base()->first = SolutionStepIndex;
// 		  i_old.base()->second = SourceSolutionStepData;
// 		}
// 	}

//       void RemoveSolutionStepNodalData(IndexType SolutionStepIndex)
// 	{
// 	    mSolutionStepsNodalData.erase(SolutionStepIndex);
// 	}

      void ClearSolutionStepsData()
	{
	    mSolutionStepsNodalData.Clear();
	}

      void SetSolutionStepVariablesList(VariablesList* pVariablesList)
	{
	    mSolutionStepsNodalData.SetVariablesList(pVariablesList);
	}

      
//       DataValueContainer::Pointer pSolutionStepNodalData(IndexType SolutionStepIndex)
// 	{
// 	  SolutionStepsNodalDataContainerType::iterator i;
// 	  if((i = mSolutionStepsNodalData.find(SolutionStepIndex)) == mSolutionStepsNodalData.end())
// 	    KRATOS_ERROR(std::invalid_argument, "Solution step index out of range.", *this);

// 	  return (i.base()->second);
// 	}
      
//       template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable)
// 	{
// 	  if(mSolutionStepsNodalData.empty())
// 	     CreateSolutionStepData(0);

// 	  return mSolutionStepsNodalData.back().GetValue(rThisVariable);
// 	}

//       template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
// 	{
// 	  SolutionStepsNodalDataContainerType::iterator i;
// 	  if((i = mSolutionStepsNodalData.find(SolutionStepIndex)) == mSolutionStepsNodalData.end())
// 	    KRATOS_ERROR(std::invalid_argument, "Solution step index out of range.", *this);

// 	  return i->GetValue(rThisVariable);
// 	}

//       VariablesListDataValueContainer& SolutionStepNodalData()
// 	{
// 	  return mSolutionStepsNodalData;
// 	}
      
      VariablesListDataValueContainer& SolutionStepData()
	{
	  return mSolutionStepsNodalData;
	}
      
      DataValueContainer& Data()
	{
	  return mData;
	}

      template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable)
	{
	  return mSolutionStepsNodalData.GetValue(rThisVariable);
	}

      template<class TVariableType> typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable) const
	{
	  return mSolutionStepsNodalData.GetValue(rThisVariable);
	}

      template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
	{
	  return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);
	}
 
	template<class TVariableType> typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
	{
	  return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);
	}


      template<class TDataType> bool SolutionStepsDataHas(const Variable<TDataType>& rThisVariable) const
	{
	  return mSolutionStepsNodalData.Has(rThisVariable);
	}
      template<class TAdaptorType> bool SolutionStepsDataHas(const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return mSolutionStepsNodalData.Has(rThisVariable);
	}

	//*******************************************************************************************
	//By Riccardo
	//very similar to the one before BUT throws an error if the variable does not exist
	template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable)
	{ return mSolutionStepsNodalData.FastGetValue(rThisVariable);	}
      
	template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
	{  return mSolutionStepsNodalData.FastGetValue(rThisVariable, SolutionStepIndex);	}

	template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex, IndexType ThisPosition)
	{  return mSolutionStepsNodalData.FastGetValue(rThisVariable, SolutionStepIndex, ThisPosition);	}

	template<class TVariableType> typename TVariableType::Type& FastGetCurrentSolutionStepValue(const TVariableType& rThisVariable, IndexType ThisPosition)
	{  return mSolutionStepsNodalData.FastGetCurrentValue(rThisVariable, ThisPosition);	}
//*******************************************************************************************

      template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
	{
	  return mData.GetValue(rThisVariable);
	}

      template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable) const
	{
	  return mData.GetValue(rThisVariable);
	}

      template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
	{
	  if(!mData.Has(rThisVariable))
		return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);
	  
	  return mData.GetValue(rThisVariable);
	}

      template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
	{
	  if(!mData.Has(rThisVariable))
		return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);
	  
	  return mData.GetValue(rThisVariable);
	}

	template<class TVariableType> 
		void SetValue(const TVariableType& rThisVariable, typename TVariableType::Type const& rValue)
	{
	  mData.SetValue(rThisVariable, rValue);
	}

      template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
	{
	  return mData.Has(rThisVariable);
	}
      template<class TAdaptorType> bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return mData.Has(rThisVariable);
	}
      
//       template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
// 	{
// 	  if(!mData.Has(rThisVariable))
// 	    {
// 	      SolutionStepsNodalDataContainerType::iterator i;
// 	      if((i = mSolutionStepsNodalData.find(SolutionStepIndex)) != mSolutionStepsNodalData.end())
// 		return i->GetValue(rThisVariable);
// 	    }
	  
// 	  return mData.GetValue(rThisVariable);
// 	}

      ///@}
      ///@name Neighbours
      ///@{
      
/*      void CreateSolutionStepNeighbours(IndexType SolutionStepIndex)
	{
	  if(mSolutionStepsNeighbours.find(SolutionStepIndex) == mSolutionStepsNeighbours.end())
	    mSolutionStepsNeighbours.insert(SolutionStepIndex, typename NeighboursType::Pointer(new NeighboursType()));
	}

      void CreateSolutionStepNeighbours(IndexType SolutionStepIndex, IndexType SourceSolutionStepIndex)
	{
	  typename SolutionStepsNeighboursContainerType::iterator i;
	  if(mSolutionStepsNeighbours.find(SolutionStepIndex) == mSolutionStepsNeighbours.end())
	    {
	      if((i = mSolutionStepsNeighbours.find(SourceSolutionStepIndex)) == mSolutionStepsNeighbours.end())
		mSolutionStepsNeighbours.insert(SolutionStepIndex, typename NeighboursType::Pointer(new NeighboursType));
	      else
		mSolutionStepsNeighbours.insert(SolutionStepIndex, *i);
	    }
	}

      void CreateSolutionStepNeighbours(IndexType SolutionStepIndex, NeighboursType const &  SourceSolutionStepData)
	{
	  if(mSolutionStepsNeighbours.find(SolutionStepIndex) == mSolutionStepsNeighbours.end())
		mSolutionStepsNeighbours.insert(SolutionStepIndex, SourceSolutionStepData);
	}

      void RemoveSolutionStepNeighbours(IndexType SolutionStepIndex)
	{
	    mSolutionStepsNeighbours.erase(SolutionStepIndex);
	}

      void ClearSolutionStepsNeighbours()
	{
	    mSolutionStepsNeighbours.clear();
	}

      
      typename NeighboursType::Pointer pSolutionStepNeighbours(IndexType SolutionStepIndex)
	{
	  typename SolutionStepsNeighboursContainerType::iterator i;
	  if((i = mSolutionStepsNeighbours.find(SolutionStepIndex)) == mSolutionStepsNeighbours.end())
	    KRATOS_ERROR(std::invalid_argument, "Solution step index out of range.", *this);

	  return *(i.base());
	}
      
      NeighboursType& SolutionStepNeighbours(IndexType SolutionStepIndex)
	{
	  return *pSolutionStepNeighbours(SolutionStepIndex);
	}
      
      NeighbourElementsArrayType& NeighbourElements(IndexType SolutionStepIndex)
	{return SolutionStepNeighbours(SolutionStepIndex).NeighbourElements();}

      NeighbourNodesArrayType& NeighbourNodes(IndexType SolutionStepIndex)
	{return SolutionStepNeighbours(SolutionStepIndex).NeighbourNodes();}
  */    
      
      ///@}
      ///@name Operations
      ///@{
      
      template<class TVariableType>
      inline void Fix(const TVariableType& rDofVariable)
	{
	  pAddDof(rDofVariable)->FixDof();
	}
      
      template<class TVariableType>
      inline void Free(const TVariableType& rDofVariable)
	{
	  pAddDof(rDofVariable)->FreeDof();
	}
      
      IndexType GetBufferSize()
      {
	return mSolutionStepsNodalData.QueueSize();
      }

      void SetBufferSize(IndexType NewBufferSize)
      {
	  mSolutionStepsNodalData.Resize(NewBufferSize);
      }

      ///@}
      ///@name Access
      ///@{ 
      
      const PointType& GetInitialPosition() const
	{return mInitialPosition;}
      PointType& GetInitialPosition()
	{return mInitialPosition;}
      
      double& X0()  
	{return mInitialPosition.X();}
      double& Y0()  
	{return mInitialPosition.Y();}
      double& Z0()  
	{return mInitialPosition.Z();}
      
      double X0() const
	{return mInitialPosition.X();}
      double Y0() const 
	{return mInitialPosition.Y();}
      double Z0() const  
	{return mInitialPosition.Z();}
      
      void SetInitialPosition(const PointType& NewInitialPosition) 
	{
	  mInitialPosition.X() = NewInitialPosition.X();
	  mInitialPosition.Y() = NewInitialPosition.Y();
	  mInitialPosition.Z() = NewInitialPosition.Z();
	}
      
      void SetInitialPosition(double X,double Y, double Z)
	{
	  mInitialPosition.X() = X;
	  mInitialPosition.Y() = Y;
	  mInitialPosition.Z() = Z;
	}

	VariablesList * pGetVariablesList()
	{
		return mSolutionStepsNodalData.pGetVariablesList();
	}


        /// TODO: remove this function when removing data_file_io object.
//        IndexType& DepricatedIdAccess()
//        {
//            return BaseType::DepricatedIdAccess();
//        }

      ///@}
      ///@name Dofs
      ///@{

	//advanced functions by Riccardo
      template<class TVariableType>
	inline unsigned int GetDofPosition(TVariableType const& rDofVariable)
	{
	  typename DofsContainerType::iterator it=mDofs.find(rDofVariable.Key());
      return it - mDofs.begin();
	}

	template<class TVariableType>
	inline DofType& GetDof(TVariableType const& rDofVariable, int pos)
	{
	  typename DofsContainerType::iterator it_begin=mDofs.begin();	  
	  typename DofsContainerType::iterator it_end=mDofs.end();	
	  typename DofsContainerType::iterator it;	
	  //if the guess is exact return the guess
	  if(pos < it_end-it_begin)
	  {
		  it = it_begin + pos;
		  if( (it)->GetVariable() == rDofVariable)
			  return *it;
	  }

	  //otherwise do a find 
	  it = mDofs.find(rDofVariable.Key());
	  return *(it);
	}

      /** retuns the Dof asociated with variable  */
      template<class TVariableType>
      inline DofType& GetDof(TVariableType const& rDofVariable)
	{
	  typename DofsContainerType::iterator it=mDofs.find(rDofVariable.Key());
	  if ( it!= mDofs.end() )
	      return *it;

	  std::stringstream buffer;
	  buffer << "Not existant DOF in node #" << Id() << " for variable : " << rDofVariable.Name();
	  KRATOS_ERROR(std::invalid_argument, buffer.str(), "");

	}
	
	/** retuns all of the Dofs  */
	DofsContainerType& GetDofs(){return mDofs;}
      
      /** retuns a counted pointer to the Dof asociated with variable  */
      template<class TVariableType>
      inline typename DofType::Pointer pGetDof(TVariableType const& rDofVariable)
	{
	  typename DofsContainerType::iterator it=mDofs.find(rDofVariable.Key());
	  if ( it!= mDofs.end() )
	    return *(it.base());

	  std::stringstream buffer;
	  buffer << "Not existant DOF in node #" << Id() << " for variable : " << rDofVariable.Name();
	  KRATOS_ERROR(std::invalid_argument, buffer.str(), "");

	}
      
	/** adds a Dof to the node and return new added dof or existed one. */
	template<class TVariableType>
			inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable)
	{
		KRATOS_TRY_LEVEL_3
	    
				typename DofsContainerType::iterator i_dof = mDofs.find(rDofVariable);
		if(i_dof != mDofs.end())
			return *(i_dof.base());
		return *(mDofs.insert(mDofs.begin(), DofType(Id(), &mSolutionStepsNodalData, rDofVariable)).base());

		KRATOS_CATCH_LEVEL_3(*this);
	  
	}
      
	/** adds a Dof to the node and return new added dof or existed one. */
			inline typename DofType::Pointer pAddDof(DofType const& SourceDof)
	{
		KRATOS_TRY_LEVEL_3
	    
		typename DofsContainerType::iterator i_dof = mDofs.find(SourceDof.GetVariable());
		if(i_dof != mDofs.end())
			return *(i_dof.base());
		
		typename DofType::Pointer p_new_dof(  *(mDofs.insert(mDofs.begin(), SourceDof).base())  );

                p_new_dof->SetId(Id());
//		p_new_dof->Id() = Id();
		
		p_new_dof->SetSolutionStepsData(&mSolutionStepsNodalData);
		
		return p_new_dof;

		KRATOS_CATCH_LEVEL_3(*this);
	  
	}
      
      /** adds a Dof to the node and return new added dof or existed one. */
      template<class TVariableType, class TReactionType>
      inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable, TReactionType const& rDofReaction)
	{
	  KRATOS_TRY_LEVEL_3
	    
	  typename DofsContainerType::iterator i_dof = mDofs.find(rDofVariable);
	  if(i_dof != mDofs.end())
	  {
	      i_dof->SetReaction(rDofReaction);
	      return *(i_dof.base());
	  }
	  return *(mDofs.insert(mDofs.begin(), DofType(Id(), &mSolutionStepsNodalData, rDofVariable, rDofReaction)).base());

	  KRATOS_CATCH_LEVEL_3(*this);
	  
	}
      
      /** adds a Dof to the node and return new added dof or existed one. */
      template<class TVariableType>
      inline DofType& AddDof(TVariableType const& rDofVariable)
	{
	  KRATOS_TRY_LEVEL_3
	  
	  typename DofsContainerType::iterator i_dof = mDofs.find(rDofVariable);
	  if(i_dof != mDofs.end())
	    return *i_dof;
	  return *(mDofs.insert(mDofs.begin(), DofType(Id(), &mSolutionStepsNodalData, rDofVariable)));

	  KRATOS_CATCH_LEVEL_3(*this);
	  
	}
      
      /** adds a Dof to the node and return new added dof or existed one. */
      template<class TVariableType, class TReactionType>
      inline DofType& AddDof(TVariableType const& rDofVariable, TReactionType const& rDofReaction)
	{
	  KRATOS_TRY_LEVEL_3
	    
	  typename DofsContainerType::iterator i_dof = mDofs.find(rDofVariable);
	  if(i_dof != mDofs.end())
	    return *i_dof;
	  return *(mDofs.insert(mDofs.begin(), DofType(Id(), &mSolutionStepsNodalData, rDofVariable, rDofReaction)));

	  KRATOS_CATCH_LEVEL_3(*this);
	  
	}

      ///@}
      ///@name Inquiry
      ///@{
      
      /** Return true if the dof of freedom is present on the node */
      inline bool HasDofFor(const VariableData& rDofVariable) const
	{
	  return (mDofs.find(rDofVariable) != mDofs.end());
	}
      
      inline  bool IsFixed(const VariableData& rDofVariable) const
	{
	  typename DofsContainerType::const_iterator i;
	  return (((i= mDofs.find(rDofVariable)) == mDofs.end()) ? false : i->IsFixed());
	}
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  std::stringstream buffer;
	  buffer << "Node #" << Id();
	  return buffer.str();
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info();
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  BaseType::PrintData(rOStream);
	  if(!mDofs.empty())
	      rOStream << std::endl << "    Dofs :" << std::endl;
	  for(typename DofsContainerType::const_iterator i = mDofs.begin() ; i != mDofs.end() ; i++)
	      rOStream << "        " << i->Info() << std::endl;
// 	  rOStream << "        " << "solution steps  : " << *mSolutionStepsNodalData;
// 	  rOStream << "        " << "solution steps capacity : " << mSolutionStepsNodalData.GetContainer().capacity();
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

     #ifdef _OPENMP
        omp_lock_t mnode_lock;
     #endif
        
        
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
        
      /** storage for the dof of the node */
      DofsContainerType  mDofs;

      /** A pointer to data related to this node. */
      DataValueContainer mData;

      SolutionStepsNodalDataContainerType mSolutionStepsNodalData;

      //SolutionStepsNeighboursContainerType mSolutionStepsNeighbours;

      ///Initial Position of the node
      PointType mInitialPosition;
      
        
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
        
    }; // Class Node 


  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<std::size_t TDimension, class TDofType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Node<TDimension, TDofType>& rThis);

  /// output stream function
  template<std::size_t TDimension, class TDofType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Node<TDimension, TDofType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << " : ";
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  


	//*********************************************************************************
	//*********************************************************************************
	//*********************************************************************************
	//definition of the NEIGHBOUR_NODES variable
	//*********************************************************************************
	//*********************************************************************************
	//*********************************************************************************
	KRATOS_DEFINE_VARIABLE(WeakPointerVector<Node<3> >, NEIGHBOUR_NODES)

//     namespace Globals
//     {
// 	extern Node<3> DefaultNode3;
//     }


  
}  // namespace Kratos.

#endif // KRATOS_NODE_H_INCLUDED  defined 


