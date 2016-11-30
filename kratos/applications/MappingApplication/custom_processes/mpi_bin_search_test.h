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

#if !defined(KRATOS_CustomMpiSearchTest_PROCESS_H_INCLUDED )
#define  KRATOS_CustomMpiSearchTest_PROCESS_H_INCLUDED


// System includes

#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

// External includes
//@{KRATOS_EXTERNA_INCLUDES}
#include "includes/kratos_flags.h"
#include <boost/python.hpp>

// Project includes

#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

#include "mpi_node_configure.h"
#include "custom_utilities/bins_dynamic_objects_mpi.h"

// External includes
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif

#include <omp.h>
#include "utilities/timer.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
	typedef array_1d<double,3> array_3d;
	typedef Node < 3 > PointType;
	typedef Node < 3 > ::Pointer PointTypePointer;
	typedef std::vector<PointType::Pointer> PointVector;
	typedef std::vector<PointType::Pointer>::iterator PointIterator;
	typedef std::vector<double> DistanceVector;
	typedef std::vector<double>::iterator DistanceIterator;
	typedef ModelPart::ConditionsContainerType ConditionsArrayType;

	typedef Element BaseType;
	typedef BaseType::GeometryType GeometryType;

  typedef std::vector<PointTypePointer> neighborsVector;
  // typedef PointerVector<PointType> neighborsVector;

  typedef std::map<unsigned int, unsigned int> neighborMap;

  typedef ModelPart::NodesContainerType NodesContainerType;
//  typedef ModelPart::NodesContainerType NodePointer;

  typedef typename NodeConfigure::IteratorType            IteratorType;

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

class CustomMpiSearchTest : public Process {
public:

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of CustomMpiSearchTest
  KRATOS_CLASS_POINTER_DEFINITION(CustomMpiSearchTest);

  ///@}
  ///@name Life Cycle
  ///@{
  CustomMpiSearchTest(ModelPart& i_m_model_part):
	  m_model_part(i_m_model_part) {
      // BinsObjectDynamic<NodeConfigure> bin_structure;
      //BinsObjectDynamic<NodeConfigure> bin_structure = BinsObjectDynamic<NodeConfigure>(i_m_model_part.Nodes().ptr_begin(), i_m_model_part.Nodes().ptr_end());

      m_num_nodes = m_model_part.Nodes().size();

      m_initial_search_radius = 0.01;
      m_max_search_iterations = 15;

      std::cout << "MPI Bin Stuff initialized" << std::endl;

      int mpiRank;
      int mpiSize;

      MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
      MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

      constexpr std::size_t maxNumR = 10000;
      constexpr double      searchR = .1f;

      std::size_t numThreads = omp_get_max_threads();

      MpiNodeConfigure::ContainerType node_list;

      std::vector<double> rdsList(m_num_nodes, searchR);
      std::vector<std::vector<MpiNodeConfigure::PointerType>> resList(m_num_nodes, std::vector<MpiNodeConfigure::PointerType>(maxNumR));
      std::vector<std::vector<double>> resDist(m_num_nodes, std::vector<double>(maxNumR));
      std::vector<std::size_t> resNumb(m_num_nodes, 0);

      // for(std::size_t i = m_num_nodes * mpiRank; i < m_num_nodes * (mpiRank + 1); i++) {
      //   node_list.push_back(Point<3, double>::Pointer(new Point<3, double>(i,i,i)));
      // }

      // for (auto &node_itr : m_model_part.Nodes()) {
      //     node_list.push_back(node_itr);
      // }

      for(auto node_itr = m_model_part.Nodes().ptr_begin(); node_itr != m_model_part.Nodes().ptr_end(); node_itr++) {
          node_list.push_back(*node_itr);
          // break;
      }

      if(mpiRank == 0) {
        std::cout << "NumThreads per mpi process: " << numThreads << std::endl;
        std::cout << "NodesList Size per process: " << node_list.size() << std::endl;
      }


      auto startGen = std::chrono::steady_clock::now();


      BinsObjectDynamicMpi<MpiNodeConfigure> testBins(m_model_part.Nodes().ptr_begin(), m_model_part.Nodes().ptr_end());

      std::unordered_set<std::size_t> partitionList;
      if(mpiRank == 0) {
        // void SearchPartition(SearchStructureType & Box, const PointerType& i_object, std::unordered_set<std::size_t> & partitionList)
        for (auto & itr : node_list) {
            testBins.SearchPartition(itr, partitionList);
            std::cout << "" << std::endl;

            for(auto a : partitionList ) {
                std::cout << mpiRank << " partition " << a << std::endl;
            }
            std::cout << "" << std::endl;
            partitionList.clear();
        }
        std::cout << "" << std::endl;

        for(auto a : partitionList ) {
          std::cout << mpiRank << " partition " << a << std::endl;
        }
        std::cout << "" << std::endl;

      }

      auto endGen = std::chrono::steady_clock::now();
      auto difGen = endGen - startGen;

      if(mpiRank == 0) {
        std::cout << std::chrono::duration<double, std::milli>(difGen).count() << " ms" << std::endl;
      }


      if(mpiRank == 0) {
        auto cellLocalNumb = testBins.GetDivisions();
        auto cellLocalSize = testBins.GetCellSize();

        auto cellPartNumb = testBins.GetDivisions();
        auto cellPartSize = testBins.GetCellSize();

        // std::cout << "cellPartNumb: " << cellPartNumb[0] << " " << cellPartNumb[1] << " " << cellPartNumb[2] << std::endl;
        // std::cout << "cellPartSize: " << cellPartSize[0] << " " << cellPartSize[1] << " " << cellPartSize[2] << std::endl;

        // auto cells = testBins.GetCellContainer();
        // for(std::size_t i = 0 ; i < cells.size(); i++) {
        //   auto myCell = (*cells[i].GetObject(0))();
        //   std::cout << "cells(" << i << ") " << myCell.size() << std::endl;
        //   for(auto p : myCell) {
        //     std::cout << "\t" << p << std::endl;
        //   }
        // }

        std::cout << "Searching" << std::endl;
      }


      auto startSearch = std::chrono::steady_clock::now();

      testBins.SearchObjectsInRadius(
        node_list.begin(), node_list.end(), node_list.size(),
        rdsList, resList, resDist, resNumb, maxNumR, &m_model_part.GetCommunicator()
      );


      auto endSearch = std::chrono::steady_clock::now();
      auto difSearch = endSearch - startSearch;

      if(mpiRank == 0) {
        std::cout << std::chrono::duration<double, std::milli>(difSearch).count() << " ms" << std::endl;
      }


      if(mpiRank == 0) {
        auto results = 0;

        for(int i = 0; i < m_num_nodes; i++) {
          results += resNumb[i];
        }

        std::cout << "results:" << results << std::endl;
      }


  }

  /// Destructor.
  virtual ~CustomMpiSearchTest() {
  }

  ///@}
  ///@name Operators
  ///@{

  void operator()() {
    Execute();
  }

  ///@}
  ///@name Operations
  ///@{

  virtual void Execute() {
  }

  virtual void Clear() {
  }

  /// Turn back information as a string.
  virtual std::string Info() const {
      return "CustomMpiSearchTest";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << "CustomMpiSearchTest";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const {
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

  int m_num_nodes;

  double m_initial_search_radius;
  int m_max_search_iterations;

  ModelPart& m_model_part; 	// The model part of the interface origin

  /// Assignment operator.
  CustomMpiSearchTest& operator=(CustomMpiSearchTest const& rOther);

  /// Copy constructor.
  //CustomMpiSearchTest(CustomMpiSearchTest const& rOther);

  ///@}

}; // Class CustomMpiSearchTest

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


inline std::istream & operator >>(std::istream& rIStream,
		CustomMpiSearchTest& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const CustomMpiSearchTest& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_CUSTOM_NEAREST_NEIGHBOR_MAPPER_PROCESS_H_INCLUDED  defined
