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

#if !defined(KRATOS_CUSTOM_MORTAR_MAPPER_PROCESS_H_INCLUDED )
#define  KRATOS_CUSTOM_MORTAR_MAPPER_PROCESS_H_INCLUDED


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

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"
#include "custom_base_mapper_process.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
	typedef typename SparseSpaceType::MatrixType SparseMatrixType;
	typedef typename SparseSpaceType::VectorType VectorType;

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

class CustomMortarMapperProcess : public Process {
public:

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of CustomMortarMapperProcess
  KRATOS_CLASS_POINTER_DEFINITION(CustomMortarMapperProcess);

  ///@}
  ///@name Life Cycle
  ///@{
  CustomMortarMapperProcess(ModelPart& i_mr_master_model_part, ModelPart& i_mr_slave_model_part):
	  mr_master_model_part(i_mr_master_model_part), mr_slave_model_part(i_mr_slave_model_part)  {

	  // Initialize required matrices
	        numNodesOnMaster = mr_master_model_part.Nodes().size();
	        numNodesOnSlave = mr_slave_model_part.Nodes().size();

	        numElementsOnMaster = mr_master_model_part.Elements().size();
	        numElementsOnSlave = mr_slave_model_part.Elements().size();
  }
/*
  CustomRbfMapperProcess(ModelPart& i_mr_master_model_part, ModelPart& i_mr_slave_model_part, double i_radius, int i_orderOfBasis=1):
                    mr_master_model_part(i_mr_master_model_part), mr_slave_model_part(i_mr_slave_model_part)  {
      this -> orderOfBasis 	= i_orderOfBasis;
      this -> radius 		= i_radius;
      // Initialize requred matrices
      numNodesOnMaster = mr_master_model_part.Nodes().size();
      numNodesOnSlave = mr_slave_model_part.Nodes().size();

      this->m_Maa.resize(this->numNodesOnMaster, this->numNodesOnMaster);
      this->m_Mba.resize(this->numNodesOnSlave, this->numNodesOnMaster);
  }
*/
  /// Destructor.
  virtual ~CustomMortarMapperProcess() {
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

  /* This methods computes the mapping matrices*/
  bool Initialize() {

	  // ***************************************************
/*
	  1. Find Neighbors(Elements/Nodes) (What about the normals, do we need to flip elements?)

	  2. Put many GPs (only on Slave Side)
	  	  2.1 Put GPs in Parametric Domain
	  	  2.2 Map them into Physical Domain

	  3. Project GPs from Slave to Master using an orthogonal projection
	  	  3.1 Project them onto all the Neighbors found in 1. and select the closest
	  	  3.2 In Case some GPs can not be projected use Nearest Neighbor
	  	  If no GP is projected for one element then we do nearest neighbor for that (i.e. zero row in mapping matrix)
	  	  	  3.3.1 Exclude these Nodes from the Mortar Matrix, i.e. modify the linear system accordingly

	  4. Perform integration: Compute C_BB & C_BA (OpenMP parallel)

	  5. Optional: Enforce Consistency


	  actually steps 2-4 should be able to be done in parallel
*/
	  PointType::CoordinatesArrayType cords;
	  PointType::CoordinatesArrayType cordsRes;
	  /*
      for(ModelPart::ElementIterator it = (mr_master_model_part).ElementsBegin();it != (mr_master_model_part).ElementsEnd(); ++it){
		  std::cout << "Element ID = " << it->Id() << std::endl;

		  cords=( ZeroVector( 3 ) );
		  cordsRes=( ZeroVector( 3 ) );
			std::cout << it->GetGeometry().IsInside(cords,cordsRes) << std::endl;
			std::cout << cordsRes << std::endl;

	  }



	  GeometryType::Pointer rGeom = mr_master_model_part.Elements()[0].pGetGeometry();


	  if (rGeom==NULL){
		  std::cout << "NULL!!!" << std::endl;

	  }


	  cords[0]=0;
	  cords[1]=0;
	  cords[2]=0;



	  cords=( ZeroVector( 3 ) );
	  cordsRes=( ZeroVector( 3 ) );
	  cordsRes[0]=0;
	  cordsRes[1]=0;
	  cordsRes[2]=0;*/

	  //std::string INFO = mr_master_model_part.Elements()[0].Info();

	  //std::cout << rGeom->Info() << std::endl;

	  //rGeom.IsInside(cords,cordsRes);

	  //std::cout << cords[0] << std::endl;
	  //std::cout << cordsRes[0] << std::endl;

	  //std::cout << rGeom << std::endl;

	  //rGeom->

	  //std::cout << rGeom << std::endl;


	  // ***************************************************

	  std::cout << "Mortar Initialized" << std::endl;

      return true;
  }

  /* This function maps a variable from Master to Slave */
  //bool MapFromMasterToSlave(const Variable< array_1d<double,3> > & rMasterVar, Variable< array_1d<double,3> > & rSlaveVar){
  bool MapFromMasterToSlave(){ // always consistent

	  // Template it such that it can take scalar and vector values
	  // same for MapFromSlaveToMaster
	  /*
	  std::cout << "Mortar MapFromMasterToSlave" << std::endl;

	  int numNodes = mr_master_model_part.NumberOfNodes();
	  std::cout << "numNodes = " << numNodes << std::endl;

	  //mr_master_model_part.CreateNewNode(1, 1.00,0.00,0.00);
	  std::cout << "numNodes = " << mr_master_model_part.NumberOfNodes() << std::endl;
	  */
	  int mypid = mr_master_model_part.GetCommunicator().MyPID();
	  int noProc = mr_master_model_part.GetCommunicator().TotalProcesses();

	  std::cout << "MASTER; MyPID = " << mypid << " / " << noProc << std::endl;

      return true;
  }

  /* This function maps a variable from Slave to Master*/
  //bool MapFromSlaveToMaster(const Variable< array_1d<double,3> > & rSlaveVar, Variable< array_1d<double,3> > & rMasterVar){
  bool MapFromSlaveToMaster(){ // always conservative

	  //std::cout << "Mortar MapFromSlaveToMaster" << std::endl;

	  int mypid = mr_slave_model_part.GetCommunicator().MyPID();
	  int noProc = mr_slave_model_part.GetCommunicator().TotalProcesses();

	  std::cout << "SLAVE; MyPID = " << mypid << " / " << noProc << std::endl;

	  return true;
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
  virtual std::string Info() const {
      return "CustomMortarMapperProcess";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << "CustomMortarMapperProcess";
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

  unsigned int numNodesOnMaster;
  unsigned int numNodesOnSlave;

  unsigned int numElementsOnMaster;
  unsigned int numElementsOnSlave;

  ModelPart& mr_master_model_part; 	// The model part of the interface master
  ModelPart& mr_slave_model_part;   // The model part of the interface slave

  // These are the actual mapping matrices
  SparseMatrixType m_Maa;
  SparseMatrixType m_Mba;

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
  CustomMortarMapperProcess& operator=(CustomMortarMapperProcess const& rOther);

  /// Copy constructor.
  //CustomMortarMapperProcess(CustomMortarMapperProcess const& rOther);

  ///@}

}; // Class CustomMortarMapperProcess

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


inline std::istream & operator >>(std::istream& rIStream,
                                  CustomMortarMapperProcess& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const CustomMortarMapperProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_CUSTOM_MORTAR_MAPPER_PROCESS_H_INCLUDED  defined
