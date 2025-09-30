//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2014-03-08 08:56:42 $
//
//

#if !defined(KRATOS_DERIVATIVE_RECOVERY)
#define  KRATOS_DERIVATIVE_RECOVERY

// System includes
#include <string>
#include <iostream>
#include <cstdlib>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "mollification/density_function_polynomial.h"
#include "custom_functions.h"

//Database includes
#include "spaces/ublas_space.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/search/point_point_search.h"
#include "custom_utilities/fields/velocity_field.h"

#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

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

/// This class constructs nodal approximations of the derivatives of a field given a Lagrangian linear FEM approximation of it.
/** @author  Guillermo Casas Gonzalez <gcasas@cimne.upc.edu>
*/

template <std::size_t TDim>
class KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery
{
public:
///@name Type Definitions
///@{
typedef ModelPart::ElementsContainerType                      ElementsArrayType;
typedef ElementsArrayType::ContainerType                      ResultElementsContainerType;
typedef std::vector<ResultElementsContainerType>              VectorResultElementsContainerType;

typedef ModelPart::NodesContainerType                         NodesArrayType;
typedef NodesArrayType::ContainerType                         ResultNodesContainerType;
typedef std::vector<ResultNodesContainerType>                 VectorResultNodesContainerType;
typedef std::vector<Node::Pointer>                         NodalPointersContainerType;
typedef ModelPart::NodesContainerType::iterator               NodeIteratorType;
typedef ModelPart::ElementsContainerType::iterator            ElementIteratorType;

typedef std::size_t                                           ListIndexType;
typedef SpatialSearch::DistanceType                           DistanceType;
typedef SpatialSearch::VectorDistanceType                     VectorDistanceType;

// For the full projection L2 recover
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef UblasSpace<double, DenseMatrix<double>, DenseVector<double>> DenseSpaceType;
// typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
// typedef typename TSparseSpace::MatrixType TSystemMatrixType;
// typedef typename TSparseSpace::VectorType TSystemVectorType;
// typedef typename TSparseSpace::MatrixType MatrixType;
// typedef typename TSparseSpace::VectorType VectorType;

/// Pointer definition of DerivativeRecovery
typedef DerivativeRecovery<TDim> DerivativeRecovery_TDim;
KRATOS_CLASS_POINTER_DEFINITION(DerivativeRecovery_TDim);

///@}
///@name Life Cycle
///@{

/// Default constructor.

DerivativeRecovery(ModelPart& r_model_part, Parameters& r_parameters):
    mModelPart(r_model_part),
    mMyCustomFunctions(),
    mFirstGradientRecovery(true),
    mFirstLaplacianRecovery(true),
    mSomeCloudsDontWork(false),
    mCalculatingTheGradient(false),
    mCalculatingTheLaplacian(false),
    mCalculatingGradientAndLaplacian(false),
    mFirstTimeAppending(true)
{
    mStoreFullGradient = r_parameters.GetValue("store_full_gradient_option").GetBool();
    mComputeExactL2 = r_parameters.GetValue("compute_exact_L2").GetBool();
    mMassMatrixAlreadyComputed = false;
}

/// Destructor.
virtual ~DerivativeRecovery(){}

///@}
///@name Operators
///@{

///@}
///@name Operations
///@{

//***************************************************************************************************************
//***************************************************************************************************************
void AddTimeDerivative(ModelPart& r_model_part, Variable<array_1d<double, 3> >& material_derivative_container);

void AddTimeDerivativeComponent(ModelPart& r_model_part, Variable<array_1d<double, 3> >& material_derivative_container, const int i_component);

void RecoverGradientOfAScalar(const VariableData& origin_variable, const VariableData& destination_variable);

template <class TScalarVariable>
void RecoverSuperconvergentGradient(ModelPart& r_model_part,  TScalarVariable& scalar_container, Variable<array_1d<double, 3> >& gradient_container, unsigned int& ord);

template <class TScalarVariable>
void RecoverSuperconvergentGradientAlt(ModelPart& r_model_part,  TScalarVariable& scalar_container, Variable<array_1d<double, 3> >& gradient_container, unsigned int& ord);

void RecoverSuperconvergentMatDeriv(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& mat_deriv_container, unsigned int& ord);

void RecoverSuperconvergentLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container, unsigned int& ord);

void RecoverSuperconvergentVelocityLaplacianFromGradient(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container, unsigned int& ord);

void RecoverSuperconvergentMatDerivAndLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& mat_deriv_container, Variable<array_1d<double, 3> >& laplacian_container, unsigned int& ord);

void CalculateLocalMassMatrix(const unsigned&, const ModelPart::ElementsContainerType::iterator&, Matrix&);

void CalculateLocalRHS(const unsigned&, const ModelPart::ElementsContainerType::iterator&, Vector&);

void ComputeNonZeroMassMatrixIndex(ModelPart&, const std::vector<std::map<unsigned, unsigned>>&, std::vector<std::unordered_set<unsigned>>&);

void ConstructMassMatrixStructure(ModelPart&, const unsigned&, const std::vector<std::map<unsigned, unsigned>>&, SparseSpaceType::MatrixType&);

void AssembleMassMatrix(SparseSpaceType::MatrixType& global_matrix, const Matrix& local_lhs, std::map<unsigned, unsigned>&);

void CalculateFieldL2Projection(ModelPart&, VelocityField::Pointer flow_field, Variable<array_1d<double, 3> >&, Variable<array_1d<double, 3> >&, Variable<array_1d<double, 3> >&);

void CalculateVectorMaterialDerivativeExactL2Parallel(ModelPart&, Variable<array_1d<double, 3> >&, Variable<array_1d<double, 3> >&, Variable<array_1d<double, 3> >&);

void CalculateVectorMaterialDerivative(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container);

void CalculateVectorMaterialDerivativeExactL2(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container);

void RecoverLagrangianAcceleration(ModelPart& r_model_part);

void CalculateVectorMaterialDerivativeFromGradient(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_gradient_container_x, Variable<array_1d<double, 3> >& vector_gradient_container_y, Variable<array_1d<double, 3> >& vector_gradient_container_z, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3>  >& material_derivative_container);

void CalculateVectorMaterialDerivativeComponent(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_component_gradient_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3>  >& material_derivative_container);

void CalculateVorticityFromGradient(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_gradient_container_x, Variable<array_1d<double, 3> >& vector_gradient_container_y, Variable<array_1d<double, 3> >& vector_gradient_container_z, Variable<array_1d<double, 3>  >& vorticity_container);

void CalculateVorticityContributionOfTheGradientOfAComponent(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_component_gradient_container, Variable<array_1d<double, 3>  >& vorticity_container);

template <class TScalarVariable>
void CalculateGradient(ModelPart& r_model_part, TScalarVariable& scalar_container, Variable<array_1d<double, 3> >& gradient_container);

void SmoothVectorField(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_field, Variable<array_1d<double, 3> >& auxiliary_veriable);

void CalculateVectorLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container);

void CalculateVelocityLaplacianRate(ModelPart& r_model_part);

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************


///@}
///@name Access
///@{

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.
virtual std::string Info() const
{
    return "";
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const {}

///@}
///@name Friends
///@{
DenseVector<unsigned int>& GetElementPartition()
{
  return (mElementsPartition);
}

DenseVector<unsigned int>& GetNodePartition()
{
  return (mNodesPartition);
}

ElementsArrayType::iterator GetElementPartitionBegin(ModelPart& r_model_part, unsigned int k)
{
  ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
  return (pElements.ptr_begin() + this->GetElementPartition()[k]);
}

ElementsArrayType::iterator GetElementPartitionEnd(ModelPart& r_model_part, unsigned int k)
{
  ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
  return (pElements.ptr_begin() + this->GetElementPartition()[k + 1]);
}

NodesArrayType::iterator GetNodePartitionBegin(ModelPart& r_model_part, unsigned int k)
{
  NodesArrayType& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
  return (pNodes.ptr_begin() + this->GetNodePartition()[k]);
}

NodesArrayType::iterator GetNodePartitionEnd(ModelPart& r_model_part, unsigned int k)
{
  NodesArrayType& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
  return (pNodes.ptr_begin() + this->GetNodePartition()[k + 1]);
}
///@}

protected:

DenseVector<unsigned int> mElementsPartition;
DenseVector<unsigned int> mNodesPartition;

private:

ModelPart& mModelPart;
CustomFunctionsCalculator<TDim> mMyCustomFunctions;
bool mFirstGradientRecovery;
bool mFirstLaplacianRecovery;
bool mSomeCloudsDontWork;
bool mCalculatingTheGradient;
bool mCalculatingTheLaplacian;
bool mCalculatingGradientAndLaplacian;
bool mFirstTimeAppending;
bool mStoreFullGradient;
bool mComputeExactL2;
bool mMassMatrixAlreadyComputed;
double mLastMeasurementTime;
double mLastPressureVariation;
double mTotalVolume;
std::vector<double> mPressures;
std::vector<DenseVector<double> > mFirstRowsOfB;
bool mMustCalculateMaxNodalArea;
double mMinFluidFraction;
double mMaxNodalAreaInv;
int mCouplingType;
int mViscosityModificationType;
int mParticlesPerDepthDistance;
VariablesList mDEMCouplingVariables;
VariablesList mFluidCouplingVariables;
PointPointSearch::Pointer mpPointPointSearch;

SparseSpaceType::MatrixType m_global_mass_matrix;  // To be declared as member variable

// neighbour lists (for mCouplingType = 3)
std::vector<double>  mSearchRadii; // list of nodal search radii (filter radii). It is a vector since spatial search is designed for varying radius
VectorResultNodesContainerType mVectorsOfNeighNodes; // list of arrays of pointers to the particle's nodal neighbours
VectorDistanceType mVectorsOfDistances; // list of arrays of distances to the particle's neighbours
VectorDistanceType mVectorsOfRadii;

struct IsCloser{
    bool operator()(std::pair<unsigned int, double> const& first_pair, std::pair<unsigned int, double> const& second_pair)
    {
        return(first_pair.second < second_pair.second || (first_pair.second == second_pair.second && first_pair.first < second_pair.first));
    }
};

inline int Factorial(const unsigned int n);
bool IsEdgeNode(Geometry<Node>::GeometriesArrayType&, Node::Pointer&);
void SetEdgeNodesAndWeights(ModelPart&);
void SetNeighboursAndWeights(ModelPart& r_model_part, const unsigned int&);
void SetNeighboursAndWeightsForTheLaplacian(ModelPart& r_model_part, const unsigned int&);
void OrderByDistance(Node::Pointer &p_node, GlobalPointersVector<Node >& neigh_nodes);
bool SetInitialNeighboursAndWeights(ModelPart& r_model_part, Node::Pointer &p_node, const unsigned int&);
bool SetNeighboursAndWeights(ModelPart& r_model_part, Node::Pointer& p_node, const unsigned int&);
double SecondDegreeTestPolynomial(const array_1d <double, 3>& coordinates);
double ThirdDegreeTestPolynomial(const array_1d <double, 3>& coordinates);
double SecondDegreeGenericPolynomial(DenseMatrix<double> C, const array_1d <double, 3>& coordinates);
double ThirdDegreeGenericPolynomial(DenseMatrix<double> C, const array_1d <double, 3>& coordinates);
bool SetWeightsAndRunLeastSquaresTest(ModelPart& r_model_part, Node::Pointer& p_node, const unsigned int&);
unsigned int GetNumberOfUniqueNeighbours(const int my_id, const GlobalPointersVector<Element>& my_neighbour_elements);
double CalculateTheMaximumDistanceToNeighbours(Node::Pointer& p_node);
double CalculateTheMaximumEdgeLength(ModelPart& r_model_part);
double CalculateTheMinumumEdgeLength(ModelPart& r_model_part);
void ComputeCoefficientsMatrix(Node::Pointer& p_node, DenseMatrix<double>& CoeffsMatrix, const unsigned int& ord, bool& is_matrix_successfully_computed);
void ClassifyEdgeNodes(ModelPart& r_model_part);
void ComputeDerivativeMonomialsVector(unsigned int& ord, DenseMatrix<double>& result);
void ComputeDerivativeMonomialsVector(array_1d<double, 3>& position, Node& r_node, unsigned int& ord, DenseMatrix<double>& result);
void ComputeGradientForVertexNode(Node& inode, unsigned int& ord, Vector& gradient, Variable<double>& scalar_container);
void ComputeGradientForEdgeNode(Node& inode, unsigned int& ord, Vector& gradient, Variable<double>& scalar_container);
void CheckNeighbours(ModelPart& r_model_part);
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************



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
DerivativeRecovery& operator=(DerivativeRecovery const& rOther);

///@}

}; // Class DerivativeRecovery

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DerivativeRecovery<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_DERIVATIVE_RECOVERY  defined
