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
 
////   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-06-20 17:08:15 $
//   Revision:            $Revision: 1.6 $
//
//



// System includes
#include <string>
#include <iostream> 
#include <fstream>


// External includes 
// External includes 
#include <boost/version.hpp>

#if BOOST_VERSION >= 103800
   #include <boost/spirit/include/classic_core.hpp>
   #include <boost/spirit/include/classic_symbols.hpp>
   #include <boost/spirit/include/classic_assign_key_actor.hpp>
   
   #include <boost/spirit/include/classic_assign_key_actor.hpp>
  #include <boost/spirit/include/classic_confix.hpp>
  #include <boost/spirit/include/classic_lists.hpp>
  #include <boost/spirit/include/classic_clear_actor.hpp>
  #include <boost/spirit/include/classic_increment_actor.hpp>
#else
   #include <boost/spirit/core.hpp> 
   #include <boost/spirit/symbols/symbols.hpp>
   #include <boost/spirit/actor/assign_key_actor.hpp>
   
   #include <boost/spirit/actor/assign_key_actor.hpp>
  #include <boost/spirit/utility/confix.hpp>
  #include <boost/spirit/utility/lists.hpp>
  #include <boost/spirit/actor/clear_actor.hpp>
  #include <boost/spirit/actor/increment_actor.hpp>
#endif
//#include <boost/spirit/core.hpp>
// #include <boost/spirit/symbols/symbols.hpp>
// #include <boost/spirit/actor/assign_key_actor.hpp>
// #include <boost/spirit/utility/confix.hpp>
// #include <boost/spirit/utility/lists.hpp>
// #include <boost/spirit/actor/clear_actor.hpp>
// #include <boost/spirit/actor/increment_actor.hpp>

// Project includes
#include "includes/define.h"
#include "includes/datafile_io.h"
#include "includes/io.h"
#include "includes/kratos_components.h"
#include "includes/mesh.h"
#include "includes/model_part.h"
#include "utilities/timer.h"

#ifdef KRATOS_INDEX_PARSER
#undef KRATOS_INDEX_PARSER
#endif
#define KRATOS_INDEX_PARSER(index) \
	( \
	'[' \
	>> uint_p[assign_a(index)] \
	>> ']' \
	)

#ifdef KRATOS_COORDINATES_PARSER
#undef KRATOS_COORDINATES_PARSER
#endif
#define KRATOS_COORDINATES_PARSER(point) \
	( \
	real_p[assign_a(point.X())] \
	>> ',' \
	>> real_p[assign_a(point.Y())] \
	>> ','  \
	>> real_p[assign_a(point.Z())] \
	)

#ifdef KRATOS_NODE_COORDINATES_PARSER
#undef KRATOS_NODE_COORDINATES_PARSER
#endif
#define KRATOS_NODE_COORDINATES_PARSER(point) \
	( \
	real_p[assign_a(point.X())][assign_a(point.X0())] \
	>> ',' \
	>> real_p[assign_a(point.Y())][assign_a(point.Y0())] \
	>> ','  \
	>> real_p[assign_a(point.Z())][assign_a(point.Z0())] \
	)

#ifdef KRATOS_NODE_PARSER
#undef KRATOS_NODE_PARSER
#endif
#define KRATOS_NODE_PARSER(node) \
	( \
	str_p("NODES") \
	>> KRATOS_INDEX_PARSER(node.DepricatedIdAccess()) \
	>> '='  \
	>> str_p("Node")  \
	>> '(' \
	>> 	KRATOS_NODE_COORDINATES_PARSER(node) \
	>> ')' \
	)

#ifdef KRATOS_NODE_DATA_PARSER
#undef KRATOS_NODE_DATA_PARSER
#endif
#define KRATOS_NODE_DATA_PARSER(node) \
	( \
	ch_p('[') \
	>> uint_p[assign_a(node.DepricatedIdAccess())] \
	>> ',' \
	>> 	KRATOS_NODE_COORDINATES_PARSER(node) \
	>> ']' \
	)

#ifdef KRATOS_PROPERTIES_LHS_PARSER
#undef KRATOS_PROPERTIES_LHS_PARSER
#endif
#define KRATOS_PROPERTIES_LHS_PARSER(index, type) \
	( \
	str_p("PROPERTIES") \
	  >> KRATOS_INDEX_PARSER(index.first) \
	  >> '[' \
	  >> ComponentParser<type >()[assign_a(index.second)] \
          >> ']' \
	)

#ifdef KRATOS_PROPERTIES_TEMPORARY_VARIABLES
#undef KRATOS_PROPERTIES_TEMPORARY_VARIABLES
#endif
#define KRATOS_PROPERTIES_TEMPORARY_VARIABLES(variable, type, index, value) \
    std::pair<unsigned int,boost::reference_wrapper<variable  const> > index(1,boost::cref(variable::StaticObject())); \
    type value;

#ifdef KRATOS_ELEMENT_DATA_PARSER
#undef KRATOS_ELEMENT_DATA_PARSER
#endif
#define KRATOS_ELEMENT_DATA_PARSER(index, element_nodes, element_properties, nodes, properties) \
	( \
	ch_p('[') \
	>> uint_p[assign_a(index)] \
	>> ',' \
	>> '[' \
	>> uint_p[push_back_pointer_from_a(element_nodes, nodes)] \
	>> *(',' \
	     >> uint_p[push_back_pointer_from_a(element_nodes, nodes)]) \
	>> ']' \
	>> ',' \
        >> uint_p[set_pointer_from_a(element_properties, properties)] \
	>> ']' \
	)

#ifdef KRATOS_CONDITIONS_FIX_PARSER
#undef KRATOS_CONDITIONS_FIX_PARSER
#endif
#define KRATOS_CONDITIONS_FIX_PARSER(nodes, index, type, variable) \
	( \
	str_p("NODES") \
	>> KRATOS_INDEX_PARSER(index) \
	>> '.' \
	>> str_p("Fix") \
	>> '(' \
	>> ComponentParser<type >()[assign_a(variable)] \
	>> ')' \
	)[node_fix_a(nodes, variable, index)] \

#ifdef KRATOS_CONDITIONS_TEMPORARY_VARIABLES
#undef KRATOS_CONDITIONS_TEMPORARY_VARIABLES
#endif
#define KRATOS_CONDITIONS_TEMPORARY_VARIABLES(type, variable) \
    boost::reference_wrapper<type const> variable(boost::cref(type::StaticObject()));

#ifdef KRATOS_ARRAY_1D_3_PARSER
#undef KRATOS_ARRAY_1D_3_PARSER
#endif
#define KRATOS_ARRAY_1D_3_PARSER(value) \
	( \
        ch_p('[') \
	>> real_p[assign_a(value[0])] \
	>> ',' \
	>> real_p[assign_a(value[1])] \
	>> ',' \
	>> real_p[assign_a(value[2])] \
	>> ']' \
	)

#ifdef KRATOS_VECTOR_PARSER
#undef KRATOS_VECTOR_PARSER
#endif
#define KRATOS_VECTOR_PARSER(vector) \
	( \
        ch_p('[') \
	>> real_p[push_back_a(vector)] \
	>> *(',' \
	     >> real_p[push_back_a(vector)]) \
	>> ']' \
	)
	
#ifdef KRATOS_MATRIX_PARSER
#undef KRATOS_MATRIX_PARSER
#endif
#define KRATOS_MATRIX_PARSER(row_vector, matrix) \
	( \
	ch_p('[') \
	>> KRATOS_VECTOR_PARSER(row_vector)[add_row_a(matrix, row_vector)][clear_a(row_vector)] \
	>> *(',' >> KRATOS_VECTOR_PARSER(row_vector)[add_row_a(matrix, row_vector)][clear_a(row_vector)]) \
	>> ']' \
	)
	
#ifdef KRATOS_NODE_INITIALIZE_LHS_PARSER
#undef KRATOS_NODE_INITIALIZE_LHS_PARSER
#endif
#define KRATOS_NODE_INITIALIZE_LHS_PARSER(index, type, node, nodes) \
	( \
 	str_p("NODES") \
	>> '[' \
	>> uint_p[set_pointer_from_a(node, nodes)] \
	>> ']' \
	>> '(' \
	>> ComponentParser<type >()[assign_a(index.second)] \
	>> ',' \
	>> uint_p[assign_a(index.first)] \
	>> ')' \
	)

namespace Kratos
{

    template<
        typename T,
        typename ValueT,
        typename ActionT
    >
    class ref_ref_value_actor : public ActionT
    {
    private:
        T& ref;
        ValueT& value_ref;
    public:
        ref_ref_value_actor(
            T& ref_,
            ValueT& value_ref_
            )
        :
            ref(ref_),
            value_ref(value_ref_)
        {}


        template<typename T2>
        void operator()(T2 const& val_) const
        {
            this->act(ref,value_ref,val_); // defined in ActionT
        }


        template<typename IteratorT>
            void operator()(
            IteratorT const& first_,
            IteratorT const& last_
            ) const
        {
            this->act(ref,value_ref,first_,last_); // defined in ActionT
        }
    };

      template<class TComponentType>
      class  ComponentParser : public KRATOS_BOOST_SPIRIT::symbols<boost::reference_wrapper<TComponentType const> >
      {
	  public:
    	ComponentParser()
	  {
	    typedef typename KratosComponents<TComponentType>::ComponentsContainerType::const_iterator iterator_type;
	    for(iterator_type i_component = KratosComponents<TComponentType>::GetComponents().begin() ; 
		i_component != KratosComponents<TComponentType>::GetComponents().end() ; 
		i_component++)
	      add(i_component->first.c_str(), i_component->second);
    	}
	
      };


    struct set_properties_action
    {
        template<
            typename T,
            typename ValueT,
            typename KeyT
        >
        void act(T& ref_, ValueT const& value_, KeyT const& key_) const
        {
            ref_[ key_.first ][ key_.second.get() ] = value_;
        }

    };

    struct set_pointer_from_action
    {
        template<
            typename T,
            typename ValueT,
            typename KeyT
        >
        void act(T& ref_, ValueT & value_, KeyT const& key_) const
        {
            ref_ = value_( key_);
        }

    };

    struct node_init_action
    {
        template<
            typename T,
            typename ValueT,
            typename KeyT
        >
        void act(T& ref_, ValueT const& value_, KeyT const& key_) const
        {
            ref_->GetSolutionStepValue(key_.second.get(), key_.first) = value_;
        }

    };

    struct push_back_pointer_from_action
    {
        template<
            typename T,
            typename ValueT,
            typename KeyT
        >
        void act(T& ref_, ValueT & value_, KeyT const& key_) const
        {
            ref_.push_back(value_(key_));
        }

    };

    struct create_element_action
    {
        template<
            typename T,
            typename ValueT,
            typename KeyT
        >
        void act(T& ref_, ValueT & value_, KeyT const& key_) const
        {
             ref_.push_back(key_.second.get().Create(key_.first, value_.first, value_.second));
        }

    };

    struct node_fix_action
    {
        template<
            typename T,
            typename ValueT,
            typename KeyT
        >
        void act(T& ref_, ValueT & value_, KeyT const& key_) const
        {
             ref_[key_].Fix(value_.get());
        }

    };

    struct add_row_action
    {
        template<
            typename T,
            typename ValueT
        >
        void act(T& ref_, ValueT const& row_) const
        {
			std::size_t size1 = ref_.size1();
			if(ref_.size1() == 0)
			{
				std::size_t size2 = row_.size();
				ref_.resize(size1+1,size2);
				for(std::size_t i = 0 ; i < size2 ; i++)
					ref_(size1,i) = row_[i];
			}
			else
			{
				std::size_t size2 = ref_.size2();
				ref_.resize(size1+1,size2);
				for(std::size_t i = 0 ; i < size2 ; i++)
					ref_(size1,i) = row_[i];
			}
        }

    };


    template<
        typename T,
        typename ValueT,
        typename KeyT
    >
    inline KRATOS_BOOST_SPIRIT::ref_const_ref_const_ref_actor<
        T,
        ValueT,
        KeyT,
        set_properties_action
    > 
        set_properties_a(
            T& ref_,
            ValueT const& value_,
            KeyT const& key_
    )
    {
        return KRATOS_BOOST_SPIRIT::ref_const_ref_const_ref_actor<
            T,
            ValueT,
            KeyT,
            set_properties_action
        >(
            ref_,
            value_,
            key_
            );
    }

    template<
        typename T,
        typename ValueT
    >
    inline KRATOS_BOOST_SPIRIT::ref_const_ref_actor<
        T,
        ValueT,
        add_row_action
    > 
        add_row_a(
            T& ref_,
            ValueT const& row_
    )
    {
        return KRATOS_BOOST_SPIRIT::ref_const_ref_actor<
            T,
            ValueT,
            add_row_action
        >(
            ref_,
            row_
         );
    }

    template<
        typename T,
        typename ValueT
    >
    inline ref_ref_value_actor<
        T,
        ValueT,
        set_pointer_from_action
    > 
        set_pointer_from_a(
            T& ref_,
            ValueT & value_
    )
    {
        return ref_ref_value_actor<
            T,
	    ValueT,
            set_pointer_from_action
        >(
            ref_,
            value_
            );
    }

    template<
        typename T,
        typename ValueT
    >
    inline ref_ref_value_actor<
        T,
        ValueT,
        push_back_pointer_from_action
    > 
        push_back_pointer_from_a(
            T& ref_,
            ValueT & value_
    )
    {
        return ref_ref_value_actor<
            T,
            ValueT,
            push_back_pointer_from_action
        >(
            ref_,
            value_
            );
    }

    template<
        typename T,
        typename ValueT,
        typename KeyT
    >
    inline KRATOS_BOOST_SPIRIT::ref_const_ref_const_ref_actor<
        T,
        ValueT,
        KeyT,
        create_element_action
    > 
        create_element_a(
            T& ref_,
            ValueT const& value_,
            KeyT const& key_
    )
    {
        return KRATOS_BOOST_SPIRIT::ref_const_ref_const_ref_actor<
            T,
            ValueT,
            KeyT,
            create_element_action
        >(
            ref_,
            value_,
            key_
            );
    }

    template<
        typename T,
        typename ValueT,
        typename KeyT
    >
    inline KRATOS_BOOST_SPIRIT::ref_const_ref_const_ref_actor<
        T,
        ValueT,
        KeyT,
        node_fix_action
    > 
        node_fix_a(
            T& ref_,
            ValueT const& value_,
            KeyT const& key_
    )
    {
        return KRATOS_BOOST_SPIRIT::ref_const_ref_const_ref_actor<
            T,
            ValueT,
            KeyT,
            node_fix_action
        >(
            ref_,
            value_,
            key_
            );
    }


    template<
        typename T,
        typename ValueT,
        typename KeyT
    >
    inline KRATOS_BOOST_SPIRIT::ref_const_ref_const_ref_actor<
        T,
        ValueT,
        KeyT,
        node_init_action
    > 
        node_init_a(
            T& ref_,
            ValueT const& value_,
            KeyT const& key_
    )
    {
        return KRATOS_BOOST_SPIRIT::ref_const_ref_const_ref_actor<
            T,
            ValueT,
            KeyT,
            node_init_action
        >(
            ref_,
            value_,
            key_
            );
    }
  

      /// IO constructor.
      DatafileIO::DatafileIO(const std::string& rNodeDatafile, 
	     const std::string& rPropertiesDatafile,
	     const std::string& rElementDatafile, 
	     const std::string& rConditionDatafile, 
	     const std::string& rInitialValueDatafile)
	: mNodesRead(0)
	, mPropertiesRead(0)
	, mElementsRead(0)
	, mConditionsRead(0)
	, mInitialValuesRead(0)
	, mNodeDatafileIterator(rNodeDatafile)
	, mPropertiesDatafileIterator(rPropertiesDatafile)
	, mElementDatafileIterator(rElementDatafile)
	, mConditionDatafileIterator(rConditionDatafile)
	, mInitialValueDatafileIterator(rInitialValueDatafile)
	, mNodeDatafileStream(rNodeDatafile.c_str(), std::ios::out | std::ios::app)
	, mPropertiesDatafileStream(rPropertiesDatafile.c_str(), std::ios::out | std::ios::app)
	, mElementDatafileStream(rElementDatafile.c_str(), std::ios::out | std::ios::app)
	, mConditionDatafileStream(rConditionDatafile.c_str(), std::ios::out | std::ios::app)
	, mInitialValueDatafileStream(rInitialValueDatafile.c_str(), std::ios::out | std::ios::app)
	{
	}

      /// Single stream IO constructor.
      DatafileIO::DatafileIO(const std::string& rDatafile)
	: mNodesRead(0)
	, mPropertiesRead(0)
	, mElementsRead(0)
	, mConditionsRead(0)
	, mInitialValuesRead(0)
	//, mNodeDatafileIterator(rDatafile + ".node")
	//, mPropertiesDatafileIterator(rDatafile + ".prop")
	//, mElementDatafileIterator(rDatafile + ".elem")
	//, mConditionDatafileIterator(rDatafile + ".cond")
	//, mInitialValueDatafileIterator(rDatafile + ".init")
	//, mNodeDatafileStream((rDatafile + ".node").c_str(), std::ios::out | std::ios::app)
	//, mPropertiesDatafileStream((rDatafile + ".prop").c_str(), std::ios::out | std::ios::app)
	//, mElementDatafileStream((rDatafile + ".elem").c_str(), std::ios::out | std::ios::app)
	//, mConditionDatafileStream((rDatafile + ".cond").c_str(), std::ios::out | std::ios::app)
	//, mInitialValueDatafileStream((rDatafile + ".init").c_str(), std::ios::out | std::ios::app)
	{
		std::string node_filename = rDatafile + ".node";
		std::string properties_filename = rDatafile + ".prop";
		std::string element_filename = rDatafile + ".elem";
		std::string condition_filename = rDatafile + ".cond";
		std::string init_filename = rDatafile + ".init";
	 mNodeDatafileIterator = FileIterator(node_filename);
	mPropertiesDatafileIterator = FileIterator(properties_filename);
    mElementDatafileIterator = FileIterator(element_filename);
	 mConditionDatafileIterator = FileIterator(condition_filename);
	 mInitialValueDatafileIterator = FileIterator(init_filename);

	  Timer::SetOuputFile(rDatafile + ".time");

		std::cout << rDatafile << " opened for io" << std::endl;
		std::cout << node_filename << " opened for io" << std::endl;
		std::cout << properties_filename << " opened for io" << std::endl;
		std::cout << element_filename << " opened for io" << std::endl;
		std::cout << condition_filename << " opened for io" << std::endl;
		std::cout << init_filename << " opened for io" << std::endl;

	mNodeDatafileStream.open(node_filename.c_str(), std::ios::out | std::ios::app);
	mPropertiesDatafileStream.open(properties_filename.c_str(), std::ios::out | std::ios::app);
	mElementDatafileStream.open(element_filename.c_str(), std::ios::out | std::ios::app);
	mConditionDatafileStream.open(condition_filename.c_str(), std::ios::out | std::ios::app);
	mInitialValueDatafileStream.open(init_filename.c_str(), std::ios::out | std::ios::app);

	}


      /// Destructor.
      DatafileIO::~DatafileIO(){}
       

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      bool DatafileIO::ReadNode(NodeType& rThisNode)
      {
		using namespace KRATOS_BOOST_SPIRIT;
		FileIterator end_of_file = mNodeDatafileIterator.make_end();

		return parse(mNodeDatafileIterator, end_of_file,

            //  Begin grammar
			(
				KRATOS_NODE_PARSER(rThisNode)[increment_a(mNodesRead)]
			)
            ,
            //  End grammar

            (space_p | comment_p("//") | comment_p("/*", "*/"))).full;
	
      }
      
      bool DatafileIO::ReadNodes(NodesContainerType& rThisNodes)
	{
		KRATOS_TRY
	  using namespace KRATOS_BOOST_SPIRIT;

	  ParseNodes(rThisNodes, str_p("NODES"));

		return true;
		KRATOS_CATCH("")
	}
      
      void DatafileIO::WriteNodes(NodesContainerType const& rThisNodes)
      {
		KRATOS_TRY
		  mNodeDatafileStream << "NODES = NodesList([" << std::endl;
		  for(NodesContainerType::const_iterator i_node = rThisNodes.begin() ; i_node != rThisNodes.end()-1 ; i_node++)
			  mNodeDatafileStream << "[ " << i_node->Id() 
					      << " , " << i_node->X()  
					      << " , " << i_node->Y() 
					      << " , " << i_node->Z() 
					      << "]," << std::endl;

		  mNodeDatafileStream << "[ " << rThisNodes.back().Id() 
				      << " , " << rThisNodes.back().X()  
				      << " , " << rThisNodes.back().Y() 
				      << " , " << rThisNodes.back().Z() 
				      << ']' << std::endl;

		  mNodeDatafileStream << "])" << std::endl;
		KRATOS_CATCH("")
      }

      void DatafileIO::ReadProperties(PropertiesContainerType& rThisProperties)
      {
		KRATOS_TRY
	  using namespace KRATOS_BOOST_SPIRIT;

	  typedef array_1d<double, 3> array_1d_3_type;
	  KRATOS_PROPERTIES_TEMPORARY_VARIABLES(Variable<int>, int, int_index, int_value);
	  KRATOS_PROPERTIES_TEMPORARY_VARIABLES(Variable<double>, double, double_index, double_value);
	  KRATOS_PROPERTIES_TEMPORARY_VARIABLES(Variable<array_1d_3_type>, array_1d_3_type, array3_double_index, array3_double_value);
	  KRATOS_PROPERTIES_TEMPORARY_VARIABLES(Variable<Matrix>, Matrix, matrix_double_index,  matrix_double_value);

	  std::vector<double> first_row_vector;

	  std::vector<double> temp_vector;

	  FileIterator end_of_file = mPropertiesDatafileIterator.make_end();

	  parse(mPropertiesDatafileIterator, end_of_file,

		  //  Begin grammar
		(
			*( 
			   (
			    KRATOS_PROPERTIES_LHS_PARSER(int_index, Variable<int>)
			    >> '='
			    >> int_p[assign_a(int_value)]
			    )[set_properties_a(rThisProperties, int_value, int_index)][increment_a(mPropertiesRead)]
			   | 
			   (
			    KRATOS_PROPERTIES_LHS_PARSER(double_index, Variable<double>)
			    >> '='
			    >> real_p[assign_a(double_value)]
			    )[set_properties_a(rThisProperties, double_value, double_index)][increment_a(mPropertiesRead)]
			   | 
			   (
			    KRATOS_PROPERTIES_LHS_PARSER(array3_double_index, Variable<array_1d_3_type>)
			    >> '='
			    >> KRATOS_ARRAY_1D_3_PARSER(array3_double_value)
			    )[set_properties_a(rThisProperties, array3_double_value, array3_double_index)][increment_a(mPropertiesRead)]
			   | 
			   (
			    KRATOS_PROPERTIES_LHS_PARSER(matrix_double_index, Variable<Matrix>)
			    >> '='
				>> KRATOS_MATRIX_PARSER(temp_vector, matrix_double_value)
			    )[set_properties_a(rThisProperties, matrix_double_value, matrix_double_index)][increment_a(mPropertiesRead)]
			   | 
			   anychar_p
			)
			
		) 
		//  End grammar
		,
        (space_p | comment_p("//") | comment_p("/*", "*/")));

 		KRATOS_CATCH("")
     }
      
      void DatafileIO::ReadElements(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements)
      {
		KRATOS_TRY
	  using namespace KRATOS_BOOST_SPIRIT;

	  Element dummy(0, Element::GeometryType( Element::GeometryType::PointsArrayType()));

	  std::pair<std::size_t, boost::reference_wrapper<Element const> >  index(0,boost::cref(dummy));

	  Element::Pointer p_element;

	  std::pair<Element::NodesArrayType, Properties::Pointer> element_data;

	  Element::NodesArrayType nodes;

	  Properties::Pointer  p_properties;

	  FileIterator end_of_file = mElementDatafileIterator.make_end();

	  parse(mElementDatafileIterator, end_of_file,

		//  Begin grammar
		(
			*( 
			   (
			    str_p("ELEMENTS")
			    >> '='
			    >> str_p("ElementsList")
			    >> '('
			    >> ComponentParser<Element>()[assign_a(index.second)]
			    >> ','
			    >> '['
			    >> KRATOS_ELEMENT_DATA_PARSER(index.first, element_data.first, element_data.second, rThisNodes, rThisProperties)
			    [create_element_a(rThisElements, element_data, index)][clear_a(element_data.first)][increment_a(mElementsRead)]
			    >> *(','
				 >> KRATOS_ELEMENT_DATA_PARSER(index.first, element_data.first, element_data.second, rThisNodes, rThisProperties))
			    [create_element_a(rThisElements, element_data, index)][clear_a(element_data.first)][increment_a(mElementsRead)]
			    >> ']'
			    >> ')'
			    )
			   | 
			   (
			    str_p("ELEMENTS")
			    >> KRATOS_INDEX_PARSER(index.first)
			    >> '='
			    >> ComponentParser<Element>()[assign_a(index.second)]
			    >> '('
			    >> '['
			    >> uint_p[push_back_pointer_from_a(element_data.first, rThisNodes)]
			    >> *(','
				 >> uint_p[push_back_pointer_from_a(element_data.first, rThisNodes)])
			    >> ']'
			    >> ','
			    >> uint_p[set_pointer_from_a(element_data.second,rThisProperties)]
			    >> ')'
			    )[create_element_a(rThisElements, element_data, index)][clear_a(element_data.first)][increment_a(mElementsRead)]
			   | 
			   anychar_p
			)
			
		) 
		//  End grammar
		,
        (space_p | comment_p("//") | comment_p("/*", "*/")));
 		KRATOS_CATCH("")
     }

       std::size_t DatafileIO::ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
      {
	KRATOS_TRY
	  using namespace KRATOS_BOOST_SPIRIT;

	std::size_t index;

	std::vector<std::size_t> node_indices;

	std::size_t number_of_elements = 0;

	  FileIterator end_of_file = mElementDatafileIterator.make_end();

	  parse(mElementDatafileIterator, end_of_file,

		//  Begin grammar
		(
			*( 
			   (
			    str_p("ELEMENTS")
			    >> '='
			    >> str_p("ElementsList")
			    >> '('
			    >> ComponentParser<Element>()
			    >> ','
			    >> '['
			    >> 	(
				 ch_p('[') 
				 >> uint_p
				 >> ',' 
				 >> '[' 
				 >> uint_p[push_back_a(node_indices)] 
				 >> *(',' 
				      >> uint_p[push_back_a(node_indices)]) 
				 >> ']' 
				 >> ',' 
				 >> uint_p 
				 >> ']' 
				 )[push_back_a(rElementsConnectivities, node_indices)][clear_a(node_indices)][increment_a(number_of_elements)]
	
			    >> *(','
				 >> (
				     ch_p('[') 
				     >> uint_p
				     >> ',' 
				     >> '[' 
				     >> uint_p[push_back_a(node_indices)] 
				     >> *(',' 
					  >> uint_p[push_back_a(node_indices)]) 
				     >> ']' 
				     >> ',' 
				     >> uint_p 
				     >> ']' 
				     ) [push_back_a(rElementsConnectivities, node_indices)][clear_a(node_indices)][increment_a(number_of_elements)]
				 )
			    
			    >> ']'
			    >> ')'
			    )
			   | 
			   (
			    str_p("ELEMENTS")
			    >> KRATOS_INDEX_PARSER(index)
			    >> '='
			    >> ComponentParser<Element>()
			    >> '('
			    >> '['
			    >> uint_p[push_back_a(node_indices)]
			    >> *(','
				 >> uint_p[push_back_a(node_indices)]
				 )
			    >> ']'
			    >> ','
			    >> uint_p
			    >> ')'
			    )[push_back_a(rElementsConnectivities, node_indices)][clear_a(node_indices)][increment_a(number_of_elements)]
			   | 
			   anychar_p
			)
		
		) 
		//  End grammar
		,
        (space_p | comment_p("//") | comment_p("/*", "*/")));
	  
	  return number_of_elements;

 		KRATOS_CATCH("")
     }

     void DatafileIO::WriteElements(ElementsContainerType const& rThisElements)
      {
		KRATOS_TRY
		  mElementDatafileStream << "ELEMENTS = ElementsList([" << std::endl;
		  for(ElementsContainerType::const_iterator i_element = rThisElements.begin() ; i_element != rThisElements.end()-1 ; i_element++)
		    {
		      mElementDatafileStream << "[ "  << i_element->Id() << " , [ ";
		      if(i_element->GetGeometry().size() != 0)
			mElementDatafileStream << i_element->GetGeometry()[0].Id();
		      for(std::size_t i = 1 ; i < i_element->GetGeometry().size() ; i++)
			mElementDatafileStream << " , " << i_element->GetGeometry()[i].Id();
			 mElementDatafileStream  << " ] , " << i_element->GetProperties().Id() 
						 << "]," << std::endl;
		    }

		  mElementDatafileStream << "[ " << rThisElements.back().Id() << " , [ ";
		  if(rThisElements.back().GetGeometry().size() != 0)
		    mElementDatafileStream << rThisElements.back().GetGeometry()[0].Id();
		  for(std::size_t i = 1 ; i < rThisElements.back().GetGeometry().size() ; i++)
		    mElementDatafileStream << " , " << rThisElements.back().GetGeometry()[i].Id();
		  mElementDatafileStream << " ] , " << rThisElements.back().GetProperties().Id()  
					 << ']' << std::endl;

		  mElementDatafileStream << "])" << std::endl;
		KRATOS_CATCH("")
      }

      void DatafileIO::ReadConditions(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions)
      {
		KRATOS_TRY
	  using namespace KRATOS_BOOST_SPIRIT;

	  Condition dummy(0, Condition::GeometryType( Condition::GeometryType::PointsArrayType()));

	  std::pair<std::size_t, boost::reference_wrapper<Condition const> >  index(0,boost::cref(dummy));

	  Condition::Pointer p_condition;

	  std::pair<Condition::NodesArrayType, Properties::Pointer> condition_data;

	  Condition::NodesArrayType nodes;

	  Properties::Pointer  p_properties;

	  typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

	  KRATOS_CONDITIONS_TEMPORARY_VARIABLES(Variable<double>, double_variable);

	  KRATOS_CONDITIONS_TEMPORARY_VARIABLES(array_1d_component_type, array_1d_component_variable);

	  FileIterator end_of_file = mConditionDatafileIterator.make_end();

	  parse(mConditionDatafileIterator, end_of_file,

		//  Begin grammar
		(
			*( 
			    KRATOS_CONDITIONS_FIX_PARSER(rThisNodes, index.first, Variable<double >, double_variable)[increment_a(mConditionsRead)]
			   | 
			    KRATOS_CONDITIONS_FIX_PARSER(rThisNodes, index.first, array_1d_component_type, array_1d_component_variable)[increment_a(mConditionsRead)]
			   | 
			   (
			    str_p("CONDITIONS")
			    >> '='
			    >> str_p("ConditionsList")
			    >> '('
			    >> ComponentParser<Condition>()[assign_a(index.second)]
			    >> ','
			    >> '['
			    >> KRATOS_ELEMENT_DATA_PARSER(index.first, condition_data.first, condition_data.second, rThisNodes, rThisProperties)
			    [create_element_a(rThisConditions, condition_data, index)][clear_a(condition_data.first)][increment_a(mConditionsRead)]
			    >> *(','
				 >> KRATOS_ELEMENT_DATA_PARSER(index.first, condition_data.first, condition_data.second, rThisNodes, rThisProperties))
			    [create_element_a(rThisConditions, condition_data, index)][clear_a(condition_data.first)][increment_a(mConditionsRead)]
			    >> ']'
			    >> ')'
			    )
			   | 
			   (
			    str_p("CONDITIONS")
			    >> KRATOS_INDEX_PARSER(index.first)
			    >> '='
			    >> ComponentParser<Condition>()[assign_a(index.second)]
			    >> '('
			    >> '['
			    >> uint_p[push_back_pointer_from_a(condition_data.first, rThisNodes)]
			    >> *(','
				 >> uint_p[push_back_pointer_from_a(condition_data.first, rThisNodes)])
			    >> ']'
			    >> ','
			    >> uint_p[set_pointer_from_a(condition_data.second,rThisProperties)]
			    >> ')'
			    )[create_element_a(rThisConditions, condition_data, index)][clear_a(condition_data.first)][increment_a(mConditionsRead)]
			   | 
			   anychar_p
			)
			
		) 
		//  End grammar
		,
        (space_p | comment_p("//") | comment_p("/*", "*/")));

 		KRATOS_CATCH("")
     }


      std::size_t  DatafileIO::ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
      {
	KRATOS_TRY
	  using namespace KRATOS_BOOST_SPIRIT;

	std::size_t index;

	std::vector<std::size_t> node_indices;

	std::size_t number_of_conditions = 0;

	  FileIterator end_of_file = mElementDatafileIterator.make_end();

	  parse(mElementDatafileIterator, end_of_file,

		//  Begin grammar
		(
			*( 
			   (
			    str_p("CONDITIONS")
			    >> '='
			    >> str_p("ConditionsList")
			    >> '('
			    >> ComponentParser<Condition>()
			    >> ','
			    >> '['
			    >> 	(
				 ch_p('[') 
				 >> uint_p
				 >> ',' 
				 >> '[' 
				 >> uint_p[push_back_a(node_indices)] 
				 >> *(',' 
				      >> uint_p[push_back_a(node_indices)]) 
				 >> ']' 
				 >> ',' 
				 >> uint_p 
				 >> ']' 
				 )[push_back_a(rConditionsConnectivities, node_indices)][clear_a(node_indices)][increment_a(number_of_conditions)]
	
			    >> *(','
				 >> (
				     ch_p('[') 
				     >> uint_p
				     >> ',' 
				     >> '[' 
				     >> uint_p[push_back_a(node_indices)] 
				     >> *(',' 
					  >> uint_p[push_back_a(node_indices)]) 
				     >> ']' 
				     >> ',' 
				     >> uint_p 
				     >> ']' 
				     )[push_back_a(rConditionsConnectivities, node_indices)][clear_a(node_indices)][increment_a(number_of_conditions)] 
				 )
			    
			    >> ']'
			    >> ')'
			    )
			   | 
			   (
			    str_p("CONDITIONS")
			    >> KRATOS_INDEX_PARSER(index)
			    >> '='
			    >> ComponentParser<Condition>()
			    >> '('
			    >> '['
			    >> uint_p[push_back_a(node_indices)]
			    >> *(','
				 >> uint_p[push_back_a(node_indices)]
				 )
			    >> ']'
			    >> ','
			    >> uint_p
			    >> ')'
			    )[push_back_a(rConditionsConnectivities, node_indices)][clear_a(node_indices)][increment_a(number_of_conditions)]
			   | 
			   anychar_p
			)
		
		) 
		//  End grammar
		,
        (space_p | comment_p("//") | comment_p("/*", "*/")));

	  return number_of_conditions;

 		KRATOS_CATCH("")
     }

      void DatafileIO::ReadInitialValues(NodesContainerType& rThisNodes, ElementsContainerType& rThisElements, ConditionsContainerType& rThisConditions)
      {
		KRATOS_TRY
	//NODES[1].GetSolutionStepValue(TEMPERATURE, 0) = 10.00;
 	  using namespace KRATOS_BOOST_SPIRIT;

 	  //std::size_t index;
		std::cout << " reading initiali values " << std::endl;
	  NodeType::Pointer p_node;

	  typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_variable_componet_type;

	  KRATOS_PROPERTIES_TEMPORARY_VARIABLES(Variable<double>, double, double_index, double_value);
	  KRATOS_PROPERTIES_TEMPORARY_VARIABLES(array_1d_variable_componet_type, double, array_1d_variable_componet_index, array_1d_variable_componet_value);

 	  FileIterator end_of_file = mInitialValueDatafileIterator.make_end();

 	  parse(mInitialValueDatafileIterator, end_of_file,

 		//  Begin grammar
 		(
 			*( 
 			   (
			    KRATOS_NODE_INITIALIZE_LHS_PARSER(double_index, Variable<double>, p_node, rThisNodes)
			    >> '='
			    >> real_p[assign_a(double_value)]
 			    )[node_init_a(p_node, double_value, double_index)][increment_a(mInitialValuesRead)]
 			   | 
 			   (
			    KRATOS_NODE_INITIALIZE_LHS_PARSER(array_1d_variable_componet_index, array_1d_variable_componet_type, p_node, rThisNodes)
			    >> '='
			    >> real_p[assign_a(array_1d_variable_componet_value)]
 			    )[node_init_a(p_node, array_1d_variable_componet_value, array_1d_variable_componet_index)][increment_a(mInitialValuesRead)]
 			   | 
 			   anychar_p
 			)
			
 		) 
 		//  End grammar
 		,
         (space_p | comment_p("//") | comment_p("/*", "*/")));
std::cout << " finished reading initial values " << std::endl;
 		KRATOS_CATCH("")
     }

//       void DatafileIO::ReadGeometries(NodesContainerType& rThisNodes, GeometriesContainerType& rResults)
// 	{
// 	  std::size_t id, property_id , id1, id2, id3;

// 	  while(true)
// 	    {
// 	      char buffer[256];
// 	      node_file.getline(buffer, 256);
// 	      sscanf(buffer, "ELEMENTS[%d] = TriangularDiffusionConvectionElement([%d,%d,%d],%d);", &id, &id1, &id2, &id3, &property_id);
// 	      if(strstr(buffer, "ELEMENTS") == NULL)
// 		break;
// 	      rResults.push_back(Triangle2D<Node<3> >(rThisNodes(id1),rThisNodes(id2),rThisNodes(id3)));
// 	    }
// 	}

      void DatafileIO::ReadMesh(MeshType & rThisMesh)
	{
		KRATOS_TRY
	  std::cout << "reading nodes" << std::endl;
	  ReadNodes(rThisMesh.Nodes());
	  std::cout << "Nodes succesfully read" << std::endl;
	  
	  std::cout << "reading Properties" << std::endl;
	  ReadProperties(rThisMesh.Properties());
	  std::cout << "Properties succesfully read" << std::endl;
	  
	  std::cout << "reading Elements" << std::endl;
	  ReadElements(rThisMesh.Nodes(), rThisMesh.Properties(), rThisMesh.Elements());
	  std::cout << "Properties succesfully read" << std::endl;
	  
	  std::cout << "reading Conditions" << std::endl;
	  ReadConditions(rThisMesh.Nodes(), rThisMesh.Properties(), rThisMesh.Conditions());
	  std::cout << "Properties succesfully read" << std::endl;
	  
	  std::cout << "reading Initial Values" << std::endl;
	  ReadInitialValues(rThisMesh.Nodes(), rThisMesh.Elements(), rThisMesh.Conditions());
	  std::cout << "Initial Values succesfully read" << std::endl;

	  std::cout << "Reading completed succesfully" << std::endl;

		KRATOS_CATCH("")
	}

      void DatafileIO::ReadModelPart(ModelPart & rThisModelPart)
	{
 	 // using namespace KRATOS_BOOST_SPIRIT;
	  //ParseNodes(rThisModelPart.Nodes(), str_p("MODEL") >> '[' >> str_p(rThisModelPart.Name().c_str()) >> ']' >> '.' >>  str_p("Nodes"));
		KRATOS_TRY

	  MeshType& r_mesh = rThisModelPart.GetMesh();
	  std::cout << "reading nodes" << std::endl;
	  ReadNodes(r_mesh.Nodes());
	  rThisModelPart.SetNodalSolutionStepVariablesList();
	  std::cout << "Nodes succesfully read" << std::endl;
	  
	  std::cout << "reading Properties" << std::endl;
	  ReadProperties(r_mesh.Properties());
	  std::cout << "Properties succesfully read" << std::endl;
	  
	  std::cout << "reading Elements" << std::endl;
	  ReadElements(r_mesh.Nodes(), r_mesh.Properties(), r_mesh.Elements());
	  std::cout << "Properties succesfully read" << std::endl;
	  
	  std::cout << "reading Conditions" << std::endl;
	  ReadConditions(r_mesh.Nodes(), r_mesh.Properties(), r_mesh.Conditions());
	  std::cout << "Properties succesfully read" << std::endl;
	  
	  std::cout << "reading Initial Values" << std::endl;
	  ReadInitialValues(r_mesh.Nodes(), r_mesh.Elements(), r_mesh.Conditions());
	  std::cout << "Initial Values succesfully read" << std::endl;

	  std::cout << "Reading completed succesfully" << std::endl;

		KRATOS_CATCH("")

	}

      /// Turn back information as a string.
      std::string DatafileIO::Info() const
	{
	  return "datafile io";
	}
      
      /// Print information about this object.
      void DatafileIO::PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info();
	}

      /// Print object's data.
      void DatafileIO::PrintData(std::ostream& rOStream) const
	{
	  rOStream << "    Number of nodes statement read         : " << mNodesRead << std::endl;
 	  rOStream << "    Number of properties statement read    : " << mPropertiesRead << std::endl;
	  rOStream << "    Number of elements statement read      : " << mElementsRead << std::endl;
	  rOStream << "    Number of conditions statement read    : " << mConditionsRead << std::endl;
	  rOStream << "    Number of inital values statement read : " << mInitialValuesRead << std::endl;
	}
      
            
      
      template<class TParserType>
      void DatafileIO::ParseNodes(NodesContainerType& rThisNodes, TParserType const& rThisParser)
	{
		KRATOS_TRY
#if BOOST_VERSION >= 103800
   using namespace KRATOS_BOOST_SPIRIT;
#else
   using namespace KRATOS_BOOST_SPIRIT;
#endif
	  
	  FileIterator end_of_file = mNodeDatafileIterator.make_end();
	  
	  typedef rule<scanner<FileIterator> > rule_type;
	  Node<3> temp(0,0.0,0.0,0.0);
	  parse(mNodeDatafileIterator, end_of_file,

		  //  Begin grammar
		(
			*( 
			  (rThisParser 
			   >>  KRATOS_INDEX_PARSER(temp.DepricatedIdAccess())
			   >> '=' 
			   >> str_p("Node")  
			   >> '(' 
			   >> 	KRATOS_NODE_COORDINATES_PARSER(temp) 
			   >> ')' 
			   )[push_back_a(rThisNodes, temp)][increment_a(mNodesRead)]
			|  (rThisParser 
					>>  '=' 
					>> str_p("NodesList") 
					>>'(' 
					>> '['
					>> KRATOS_NODE_DATA_PARSER(temp)[push_back_a(rThisNodes, temp)][increment_a(mNodesRead)]
					>> *(ch_p(',') 
						>> KRATOS_NODE_DATA_PARSER(temp))[push_back_a(rThisNodes, temp)][increment_a(mNodesRead)]
					>> ']' 
					>> ')'
					)
			| anychar_p
			)
			
		) 
		//  End grammar
		,
        (space_p | comment_p("//") | comment_p("/*", "*/")));
					
		KRATOS_WATCH(temp)
		KRATOS_CATCH("") 
	}
      
      ///@} 

    DatafileIO& operator >> (DatafileIO& rInput, IO::NodeType& rNode)
      {
	rInput.ReadNode(rNode);

	return rInput;
      }

    DatafileIO& operator >> (DatafileIO& rInput, IO::NodesContainerType& rNodes)
      {
	rInput.ReadNodes(rNodes);

	return rInput;
      }

    DatafileIO& operator >> (DatafileIO& rInput, IO::PropertiesContainerType& rProperties)
      {
	rInput.ReadProperties(rProperties);

	return rInput;
      }
 
    DatafileIO& operator >> (DatafileIO& rInput, IO::MeshType& rMesh)
      {
	rInput.ReadMesh(rMesh);

	return rInput;
      }
 
    DatafileIO& operator >> (DatafileIO& rInput, ModelPart& rModelPart)
      {
	rInput.ReadModelPart(rModelPart);

	return rInput;
      }
 
    DatafileIO& operator << (DatafileIO& rOutput, IO::NodesContainerType& rNodes)
      {
	rOutput.WriteNodes(rNodes);

	return rOutput;
      }
 
    DatafileIO& operator << (DatafileIO& rOutput, IO::ElementsContainerType& rElements)
      {
	rOutput.WriteElements(rElements);

	return rOutput;
      }
 
  /// input stream function
/*   inline std::istream& operator >> (std::istream& rIStream,  */
/* 				    DataIO& rThis); */

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const DatafileIO& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#undef KRATOS_INDEX_PARSER
#undef KRATOS_COORDINATES_PARSER
#undef KRATOS_NODE_COORDINATES_PARSER
#undef KRATOS_NODE_PARSER
#undef KRATOS_NODE_DATA_PARSER
#undef KRATOS_PROPERTIES_LHS_PARSER
#undef KRATOS_PROPERTIES_TEMPORARY_VARIABLES
#undef KRATOS_ELEMENT_DATA_PARSER
#undef KRATOS_ARRAY_1D_3_PARSER
#undef KRATOS_VECTOR_PARSER
#undef KRATOS_MATRIX_PARSER
#undef KRATOS_CONDITIONS_TEMPORARY_VARIABLES
#undef KRATOS_CONDITIONS_FIX_PARSER
#undef KRATOS_NODE_INITIALIZE_LHS_PARSER



