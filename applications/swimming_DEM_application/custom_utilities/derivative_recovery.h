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
#include <stdlib.h>

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
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/search/point_point_search.h"

#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"
//#include "../../DEM_application/DEM_application.h"

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
typedef std::vector<Node<3>::Pointer>                         NodalPointersContainerType;
typedef ModelPart::NodesContainerType::iterator               NodeIteratorType;

typedef std::size_t                                           ListIndexType;
typedef SpatialSearch::DistanceType                           DistanceType;
typedef SpatialSearch::VectorDistanceType                     VectorDistanceType;

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
void RecoverSuperconvergentGradient(ModelPart& r_model_part,  TScalarVariable& scalar_container, Variable<array_1d<double, 3> >& gradient_container);

void RecoverSuperconvergentMatDeriv(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& mat_deriv_container);

void RecoverSuperconvergentLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container);

void RecoverSuperconvergentVelocityLaplacianFromGradient(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container);

void RecoverSuperconvergentMatDerivAndLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& mat_deriv_container, Variable<array_1d<double, 3> >& laplacian_container);

void CalculateVectorMaterialDerivative(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container);

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

void SetNeighboursAndWeights(ModelPart& r_model_part);
void SetNeighboursAndWeightsForTheLaplacian(ModelPart& r_model_part);
void OrderByDistance(Node<3>::Pointer &p_node, WeakPointerVector<Node<3> >& neigh_nodes);
bool SetInitialNeighboursAndWeights(ModelPart& r_model_part, Node<3>::Pointer &p_node);
bool SetNeighboursAndWeights(ModelPart& r_model_part, Node<3>::Pointer& p_node);
double SecondDegreeTestPolynomial(const array_1d <double, 3>& coordinates);
double SecondDegreeGenericPolynomial(DenseMatrix<double> C, const array_1d <double, 3>& coordinates);
inline int Factorial(const unsigned int n);
bool SetWeightsAndRunLeastSquaresTest(ModelPart& r_model_part, Node<3>::Pointer& p_node);
unsigned int GetNumberOfUniqueNeighbours(const int my_id, const WeakPointerVector<Element>& my_neighbour_elements);
double CalculateTheMaximumDistanceToNeighbours(Node<3>::Pointer& p_node);
double CalculateTheMaximumEdgeLength(ModelPart& r_model_part);
double CalculateTheMinumumEdgeLength(ModelPart& r_model_part);

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
