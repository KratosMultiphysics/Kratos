//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

#if !defined(KRATOS_CUSTOM_RBF_MAPPER_PROCESS_H_INCLUDED )
#define  KRATOS_CUSTOM_RBF_MAPPER_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include "math.h"
#include <vector>

// External includes
//@{KRATOS_EXTERNA_INCLUDES}
#include "includes/kratos_flags.h"
#include "boost/smart_ptr.hpp"
#include "boost/current_function.hpp"
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// Project includes

#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"
#include "../../kratos/spatial_containers/spatial_containers.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
// #include "../../ExternalSolversApplication/external_includes/amgcl_solver.h"
// #include "../../ExternalSolversApplication/external_includes/amgcl_ns_solver.h"
// #include "../../ExternalSolversApplication/external_includes/amgcl_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"

namespace Kratos {

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
/*
This class implements a RBF mapping process with is taken from Numerical Recipes text book.
One can find the theoretical description of this method on page number 163 of the text book Numerical Recipes.
The varient implemented here in this process is of Shepard interpolation.
*/
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class CustomRbfMapperProcess : public Process{
public:

  ///@name Type Definitions
  ///@{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
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
  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of CustomRbfMapperProcess
  KRATOS_CLASS_POINTER_DEFINITION(CustomRbfMapperProcess);

  ///@}
  ///@name Life Cycle
  ///@{

  /* This constructor is specific to the RBF based mapping*/
  CustomRbfMapperProcess(typename TLinearSolver::Pointer pNewLinearSolver, ModelPart& i_mr_master_model_part, ModelPart& i_mr_slave_model_part, double i_radius, int i_orderOfBasis=1, double i_tol=1e-3):
                    mr_master_model_part(i_mr_master_model_part), mr_slave_model_part(i_mr_slave_model_part)  {
      this -> orderOfBasis 	= i_orderOfBasis;
      this -> radius 		= i_radius;
      this -> tol = i_tol;
      this->formulated = false;
      // Initialize requred matrices
      numNodesOnMaster = mr_master_model_part.Nodes().size();
      numNodesOnSlave = mr_slave_model_part.Nodes().size();

      this->m_Maa.resize(this->numNodesOnMaster + 1, this->numNodesOnMaster + 1);
      this->m_Mba.resize(this->numNodesOnSlave, this->numNodesOnMaster);

      this->mLinearEquationsSolver = pNewLinearSolver;



      Initialize();
  }


  /// Destructor.
  virtual ~CustomRbfMapperProcess() {
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
void Initialize() {

	  KRATOS_TRY

	  // Building coefficient matrix

	  for(long int i=0; i<this->numNodesOnMaster; i++){
		  for(long int j=0; j<this->numNodesOnMaster; j++){
			  double r = 0.0; //ComputeSquaredDistance(this->mr_master_model_part.Nodes()[i], this->mr_master_model_part.Nodes()[j]);
			  this->m_Maa(i,j) = ComputeMultiquadricRBF(r);
		  }
	  }

	  long int j = 0;
      for(ModelPart::NodeIterator i = (mr_master_model_part).NodesBegin();i != (mr_master_model_part).NodesEnd(); ++i){
		  double x = i->Coordinates()[0];
		  this->m_Maa(numNodesOnMaster,j) = x;
		  j = j + 1;
	  }

      j = 0;
      for(ModelPart::NodeIterator i = (mr_master_model_part).NodesBegin();i != (mr_master_model_part).NodesEnd(); ++i){
		  double x = i->Coordinates()[0];
		  this->m_Maa(j,numNodesOnMaster) = x;
		  j = j + 1;
	  }/* Commented because it throws a compiler warning, Philipp, 26.09.2016
      for(int i = 0; i<this->numNodesOnMaster+1; i++){
    	  std::cout<<std::endl;
    	  for(int j= 0; j<this->numNodesOnMaster+1; j++)
    		  std::cout<<" "<<m_Maa(i,j);
    	  std::cout<<std::endl;
      }*/

	  KRATOS_CATCH("Initializing of the mapper failed !!")

  }

  /* This function maps a variable from Master to Slave */
  bool MapFromMasterToSlave(const Variable< array_1d<double,3> > & rMasterVar,
		Variable< array_1d<double,3> > & rSlaveVar){

	if(!this->formulated){
		// Building the RHS vector
		this->rhsVector.resize(this->numNodesOnMaster,true);
		this->mappingConstants.resize(this->numNodesOnMaster,true);
		long int j = 0;
		for(ModelPart::NodeIterator i = (mr_master_model_part).NodesBegin();i != (mr_master_model_part).NodesEnd(); ++i){
			rhsVector(3*j+0) = i->FastGetSolutionStepValue(rMasterVar)[0];
			rhsVector(3*j+1) = i->FastGetSolutionStepValue(rMasterVar)[1];
			rhsVector(3*j+2) = i->FastGetSolutionStepValue(rMasterVar)[2];

			std::cout<<" "<<rhsVector(3*j+0)<<" "<<rhsVector(3*j+1)<<" "<<rhsVector(3*j+2)<<std::endl;

			j = j + 1;
		}
		// Solving for the mapping constants
		this->mLinearEquationsSolver->Solve(this->m_Maa, this->mappingConstants, rhsVector);

		this->formulated = true;
		rhsVector.clear();
		std::cout<<"Number of elements in rhs vector :: "<<rhsVector.size()<<" and in mappingConstants :: "<<this->mappingConstants.size()<<std::endl;

		return true;
	}


	// Solving for

	return false;
  }

  /* This function maps a variable from Slave to Master*/
  bool MapFromSlaveToMaster(const Variable< array_1d<double,3> > & rSlaveVar,
            Variable< array_1d<double,3> > & rMasterVar){


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
/*
  /// Turn back information as a string.
  virtual std::string Info() const {
      return "CustomRbfMapperProcess";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << "CustomRbfMapperProcess";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const {
  }
*/
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

  double radius;     // The radius of the basis function -> Here we use a spherical shape function.
  int orderOfBasis;  // The order of radial basis function to be used -> The more the better.
  double tol; 		 // The tolerance to be considered for spacial neighbourhod
  bool formulated;
  unsigned int numNodesOnMaster;
  unsigned int numNodesOnSlave;
  ModelPart& mr_master_model_part; 	// The model part of the interface master
  ModelPart& mr_slave_model_part;   // The model part of the interface slave

  // These are the actual mapping matrices
  SparseMatrixType m_Maa;
  SparseMatrixType m_Mba;

  // Mapping Constants
  typename SparseSpaceType::VectorType mappingConstants;
  typename SparseSpaceType::VectorType rhsVector;

  // Linear equations solver
  LinearSolverType::Pointer mLinearEquationsSolver;

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  /* This method calculates the distance between two given nodes
   */
  double ComputeSquaredDistance(ModelPart::NodeType& node1, ModelPart::NodeType& node2){
	  double r_squared = 0.0;
	  array_3d coord1 =  node1.Coordinates();
	  array_3d coord2 =  node2.Coordinates();
	  r_squared += (coord1[0]-coord2[0])*(coord1[0]-coord2[0]);
	  r_squared += (coord1[1]-coord2[1])*(coord1[1]-coord2[1]);
	  r_squared += (coord1[2]-coord2[2])*(coord1[2]-coord2[2]);

	  return r_squared;
  }



  /* This method calculates the radial basis function value for a give set of points
   * In this function implements a "Thin Plate Spline" function .
   * 			RBF(r) = r² * log(r)
   */
  double ComputeThinPlateSplineRBF(double r_squared){

      double rbf = r_squared * log(r_squared)/2;
      return rbf;
  }


  /* This method calculates the radial basis function value for a give set of points
   * In this function implements a "Multiquadric" function.
   * 			RBF(r) = (r² + (r_0)²)^0.5
   *  where r_0 is the scalling factor
   */
  double ComputeMultiquadricRBF(double r_squared){

        double rbf = r_squared + this->radius*this->radius;
        return sqrt(rbf);
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
  CustomRbfMapperProcess& operator=(CustomRbfMapperProcess const& rOther);

  /// Copy constructor.
  //CustomRbfMapperProcess(CustomRbfMapperProcess const& rOther);

  ///@}

}; // Class CustomRbfMapperProcess

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function

/*
inline std::istream & operator >>(std::istream& rIStream,
                                  CustomRbfMapperProcess& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const CustomRbfMapperProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}
*/
///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_CUSTOM_RBF_MAPPING_PROCESS_H_INCLUDED  defined
