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
//   Date:                $Date: 2008-01-21 17:06:57 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_DATAFILE_IO_BASE_H_INCLUDED)
#define  KRATOS_DATAFILE_IO_BASE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <fstream>

#include <boost/version.hpp>

// External includes 
#if BOOST_VERSION >= 103800
   #include <boost/spirit/include/classic_file_iterator.hpp>
   #define KRATOS_BOOST_SPIRIT boost::spirit::classic
#else
   #include <boost/spirit/iterator/file_iterator.hpp>
   #define KRATOS_BOOST_SPIRIT boost::spirit
#endif

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/kratos_components.h"
#include "includes/mesh.h"
#include "includes/model_part.h"

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
    class DatafileIO : public IO
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of DatafileIO
      KRATOS_CLASS_POINTER_DEFINITION(DatafileIO);

#if BOOST_VERSION >= 103800
   typedef boost::spirit::classic::file_iterator<> FileIterator;
#else
   typedef boost::spirit::file_iterator<> FileIterator;
#endif
      
  
      typedef IO::ConnectivitiesContainerType ConnectivitiesContainerType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// IO constructor.
      DatafileIO(const std::string& rNodeDatafile, 
	     const std::string& rPropertiesDatafile,
	     const std::string& rElementDatafile, 
	     const std::string& rConditionDatafile, 
		 const std::string& rInitialValueDatafile);

      /// Single stream IO constructor.
      DatafileIO(const std::string& rDatafile);



      /// Destructor.
      virtual ~DatafileIO();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      bool ReadNode(NodeType& rThisNode);
      
      bool ReadNodes(NodesContainerType& rThisNodes);
      
      void WriteNodes(NodesContainerType const& rThisNodes);

      void ReadProperties(PropertiesContainerType& rThisProperties);
      
      void ReadElements(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements);

      std::size_t ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities);

      void WriteElements(ElementsContainerType const& rThisElements);

      void ReadConditions(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions);

      std::size_t ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities);

      void ReadInitialValues(NodesContainerType& rThisNodes, ElementsContainerType& rThisElements, ConditionsContainerType& rThisConditions);

//       void ReadGeometries(NodesContainerType& rThisNodes, GeometriesContainerType& rResults);

      void ReadMesh(MeshType & rThisMesh);

      void ReadModelPart(ModelPart & rThisModelPart);

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
      virtual std::string Info() const;
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;
      
            
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
        
      std::size_t mNodesRead;
      std::size_t mPropertiesRead;
      std::size_t mElementsRead;
      std::size_t mConditionsRead;
      std::size_t mInitialValuesRead;
        
      FileIterator mNodeDatafileIterator;
      FileIterator mPropertiesDatafileIterator;
      FileIterator mElementDatafileIterator;
      FileIterator mConditionDatafileIterator;
      FileIterator mInitialValueDatafileIterator;
        
      std::ofstream mNodeDatafileStream;
      std::ofstream mPropertiesDatafileStream;
      std::ofstream mElementDatafileStream;
      std::ofstream mConditionDatafileStream;
      std::ofstream mInitialValueDatafileStream;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{
      
      template<class TParserType>
	void ParseNodes(NodesContainerType& rThisNodes, TParserType const& rThisParser);
	
      
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
      DatafileIO& operator=(DatafileIO const& rOther);

      /// Copy constructor.
	  DatafileIO(DatafileIO const& rOther) : mNodesRead(rOther.mNodesRead),
		  mPropertiesRead(rOther.mPropertiesRead),
		  mElementsRead(rOther.mElementsRead),
		  mConditionsRead(rOther.mConditionsRead),
		  mInitialValuesRead(rOther.mInitialValuesRead),
		  mNodeDatafileIterator(rOther.mConditionDatafileIterator),
		  mPropertiesDatafileIterator(rOther.mPropertiesDatafileIterator),
		  mElementDatafileIterator(rOther.mElementDatafileIterator),
		  mConditionDatafileIterator(rOther.mConditionDatafileIterator),
		  mInitialValueDatafileIterator(rOther.mInitialValueDatafileIterator)
		  //mNodeDatafileStream(rOther.mNodeDatafileStream),
		  //mPropertiesDatafileStream(rOther.mPropertiesDatafileStream),
		  //mElementDatafileStream(rOther.mElementDatafileStream),
		  //mConditionDatafileStream(rOther.mConditionDatafileStream),
		  //mInitialValueDatafileStream(rOther.mInitialValueDatafileStream)
	  {}
        
      ///@}    
        
    }; // Class DatafileIO 

  ///@}

  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 

  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_DATAFILE_IO_BASE_H_INCLUDED  defined 


