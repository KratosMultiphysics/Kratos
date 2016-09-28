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
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "density_function_polynomial.h"
#include "custom_functions.h"

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/point_point_search.h"

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

class DerivativeRecovery
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

DerivativeRecovery(ModelPart& r_model_part): mModelPart(r_model_part), mMyCustomFunctions(){}

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

void RecoverGradientOfAScalar(const VariableData& origin_variable, const VariableData& destination_variable);

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateVectorMaterialDerivative(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container);

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void CalculateGradient(ModelPart& r_model_part, Variable<double>& scalar_container, Variable<array_1d<double, 3> >& gradient_container);
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void RecoverSuperconvergentLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container);

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void CalculateVectorLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container);


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
vector<unsigned int>& GetElementPartition()
{
  return (mElementsPartition);
}

vector<unsigned int>& GetNodePartition()
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

vector<unsigned int> mElementsPartition;
vector<unsigned int> mNodesPartition;

private:
private:

bool mPressuresFilled;
bool mFirstGradientRecovery;
bool mFirstLaplacianRecovery;
bool mSomeCloudsDontWork;
bool mCalculatingTheGradient;
bool mCalculatingTheLaplacian;
bool mFirstTimeAppending;
double mLastMeasurementTime;
double mLastPressureVariation;
double mTotalVolume;
std::vector<double> mPressures;
std::vector<vector<double> > mFirstRowsOfB;
bool mMustCalculateMaxNodalArea;
double mMinFluidFraction;
double mMaxNodalAreaInv;
int mCouplingType;
int mViscosityModificationType;
int mParticlesPerDepthDistance;
VariablesList mDEMCouplingVariables;
VariablesList mFluidCouplingVariables;
PointPointSearch::Pointer mpPointPointSearch;
ModelPart& mModelPart;
CustomFunctionsCalculator<TDim> mMyCustomFunctions;

// neighbour lists (for mCouplingType = 3)
std::vector<double>  mSearchRadii; // list of nodal search radii (filter radii). It is a vector since spatial search is designed for varying radius
VectorResultNodesContainerType mVectorsOfNeighNodes; // list of arrays of pointers to the particle's nodal neighbours
VectorDistanceType mVectorsOfDistances; // list of arrays of distances to the particle's neighbours
VectorDistanceType mVectorsOfRadii;

//***************************************************************************************************************
//***************************************************************************************************************

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SetNeighboursAndWeights(ModelPart& r_model_part)
{
    // Finding elements concurrent to each node. The nodes of these elements will form the initial cloud of points
    FindNodalNeighboursProcess neighbour_finder = FindNodalNeighboursProcess(r_model_part);
    neighbour_finder.Execute();
    const unsigned int n_max_iterations = 100;

    unsigned int i = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        bool the_cloud_of_neighbours_is_successful = SetInitialNeighboursAndWeights(r_model_part, *(inode.base()));
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);

        unsigned int iteration = 0;
        while (!the_cloud_of_neighbours_is_successful && iteration < n_max_iterations){
            the_cloud_of_neighbours_is_successful = SetNeighboursAndWeights(r_model_part, *(inode.base()));
            iteration++;
        }

        if (iteration >= n_max_iterations){ // giving up on this method, settling for the default method
            mSomeCloudsDontWork = true;
            neigh_nodes.clear();
            inode->FastGetSolutionStepValue(NODAL_WEIGHTS).clear();
            std::cout << "Warning!, for the node with id " << inode->Id() << " it has not been possible to form an adequate cloud of neighbours\n";
            std::cout << "for the gradient recovery. A lower accuracy method has been employed for this node.";
        }
        i++;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SetNeighboursAndWeightsForTheLaplacian(ModelPart& r_model_part)
{
    // Finding elements concurrent to each node. The nodes of these elements will form the initial cloud of points
    FindNodalNeighboursProcess neighbour_finder = FindNodalNeighboursProcess(r_model_part);
    neighbour_finder.Execute();
    const unsigned int n_max_iterations = 100;

    unsigned int i = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        bool the_cloud_of_neighbours_is_successful = SetInitialNeighboursAndWeights(r_model_part, *(inode.base()));
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);

        unsigned int iteration = 0;
        while (!the_cloud_of_neighbours_is_successful && iteration < n_max_iterations){
            the_cloud_of_neighbours_is_successful = SetNeighboursAndWeights(r_model_part, *(inode.base()));
            iteration++;
        }

        if (iteration >= n_max_iterations){ // giving up on this method, settling for the default method
            mSomeCloudsDontWork = true;
            neigh_nodes.clear();
            inode->FastGetSolutionStepValue(NODAL_WEIGHTS).clear();
            std::cout << "Warning!, for the node with id " << inode->Id() << " it has not been possible to form an adequate cloud of neighbours\n";
            std::cout << "for the gradient recovery. A lower accuracy method has been employed for this node.";
        }
        i++;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

struct IsCloser{
    bool operator()(std::pair<unsigned int, double> const& first_pair, std::pair<unsigned int, double> const& second_pair)
    {
        return(first_pair.second < second_pair.second || (first_pair.second == second_pair.second && first_pair.first < second_pair.first));
    }
};

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void OrderByDistance(Node<3>::Pointer &p_node, WeakPointerVector<Node<3> >& neigh_nodes)
{
    const unsigned int n_nodes = neigh_nodes.size();
    std::vector<double> distances_squared;
    distances_squared.resize(n_nodes);
    const array_1d <double, 3>& origin = p_node->Coordinates();

    for (unsigned int i = 0; i < n_nodes; ++i){
        const array_1d <double, 3> rel_coordinates = (neigh_nodes[i] - origin);
        distances_squared[i] = DEM_INNER_PRODUCT_3(rel_coordinates, rel_coordinates);
    }
    std::vector <std::pair<unsigned int, double> > ordering;
    ordering.resize(n_nodes);

    for (unsigned int i = 0; i < n_nodes; ++i){
        ordering[i] = std::make_pair(i, distances_squared[i]);
    }
    std::sort(ordering.begin(), ordering.end(), IsCloser());
    WeakPointerVector<Node<3> > ordered_neighbours;

    for (unsigned int i = 0; i < n_nodes; ++i){
        Node<3>::WeakPointer p_neigh = neigh_nodes(ordering[i].first);
        ordered_neighbours.push_back(p_neigh);
    }

    ordered_neighbours.swap(neigh_nodes);
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

bool SetInitialNeighboursAndWeights(ModelPart& r_model_part, Node<3>::Pointer &p_node)
{
    WeakPointerVector<Element>& neigh_elems = p_node->GetValue(NEIGHBOUR_ELEMENTS);
    WeakPointerVector<Node<3> >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    std::map<std::size_t, std::size_t> ids; // map to keep track of all different ids corresponding to already added neighbours to avoid repetition
    ids[p_node->Id()] = p_node->Id();

    unsigned int i = 0;

    for (unsigned int i_el = 0; i_el < neigh_elems.size(); ++i_el){
        Geometry<Node<3> >& geom = neigh_elems[i_el].GetGeometry();

        unsigned int jj = 0; // index of the node in geom corresponding to neighbour neigh_elems[i_el]
        if (geom[jj].Id() == p_node->Id()){ // skipping itself
            jj++;
        }

        for (unsigned int j = 0; j < TDim; ++j){
            Node<3>::Pointer p_neigh = geom(jj);

            if (ids.find(p_neigh->Id()) == ids.end()){
                neigh_nodes.push_back(p_neigh);
                ids[p_neigh->Id()] = p_neigh->Id();
            }
        }
        i += TDim;
    }

    OrderByDistance(p_node, neigh_nodes);

    if (neigh_nodes.size() < 10){ // Not worthwhile checking, since there are 10 independent coefficients to be determined
        return false;
    }

    else {
        return(SetWeightsAndRunLeastSquaresTest(r_model_part, p_node));
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

bool SetNeighboursAndWeights(ModelPart& r_model_part, Node<3>::Pointer& p_node)
{
    WeakPointerVector<Node<3> >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    const unsigned int node_increase_per_neighbour = 1;
    const unsigned int node_increase_overall = 1;
    std::map<std::size_t, std::size_t> ids;
    ids[p_node->Id()] = p_node->Id();

    for (unsigned int i = 0; i < (unsigned int)neigh_nodes.size(); ++i){
        Node<3>::Pointer p_neigh = neigh_nodes(i).lock();
        ids[p_neigh->Id()] = p_neigh->Id();
    }

    const unsigned int n_neigh = neigh_nodes.size();

    for (unsigned int i = 0; i < n_neigh; ++i){
        Node<3>::Pointer p_neigh = neigh_nodes(i).lock();
        WeakPointerVector<Node<3> >& neigh_neigh_nodes = p_neigh->GetValue(NEIGHBOUR_NODES);
        unsigned int n_new_nodes = 0;
        for (unsigned int j = 0; j < (unsigned int)neigh_neigh_nodes.size(); ++j){
            Node<3>::Pointer p_neigh_neigh = neigh_neigh_nodes(j).lock();
            if (ids.find(p_neigh_neigh->Id()) == ids.end()){
                neigh_nodes.push_back(p_neigh_neigh);
                ids[p_neigh_neigh->Id()] = p_neigh_neigh->Id();
                n_new_nodes++;
            }

            if (n_new_nodes >= node_increase_per_neighbour){
                break;
            }
        }
    }

    OrderByDistance(p_node, neigh_nodes);
    const unsigned int new_size = std::min(n_neigh + node_increase_overall, (unsigned int)neigh_nodes.size());
    neigh_nodes.resize(new_size); // keeping only nearest nodes

    if (neigh_nodes.size() < 10){ // it is not worthwhile checking, since there are 10 independent coefficients to be determined
        return false;
    }

    else {
        return(SetWeightsAndRunLeastSquaresTest(r_model_part, p_node));
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SecondDegreeTestPolynomial(const array_1d <double, 3>& coordinates)
{
    const double x = coordinates[0];
    const double y = coordinates[1];
    const double z = coordinates[2];
    return(1.0 + x + y + z + x * y + x * z + y * z + x * x + y * y + z * z);
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SecondDegreeGenericPolynomial(boost::numeric::ublas::matrix<double> C, const array_1d <double, 3>& coordinates)
{
    const double x = coordinates[0];
    const double y = coordinates[1];
    const double z = coordinates[2];
    return(C(0,0) + C(1,0) * x + C(2,0) * y + C(3,0) * z + C(4,0) * x * y + C(5,0) * x * z + C(6,0) * y * z + C(7,0) * x * x + C(8,0) * y * y + C(9,0) * z * z);
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

inline int Factorial(const unsigned int n){

    if (n == 0){
        return 1;
    }

    unsigned int k = n;

    for (unsigned int i = n - 1; i > 0; --i){
        k *= i;
    }

    return k;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

bool SetWeightsAndRunLeastSquaresTest(ModelPart& r_model_part, Node<3>::Pointer& p_node)
{
    using namespace boost::numeric::ublas;

    unsigned int n_poly_terms = Factorial(TDim + 2) / (2 * Factorial(TDim)); // 2 is the polynomial order

    if (TDim == 2){
        KRATOS_THROW_ERROR(std::runtime_error,"Gradient recovery not implemented yet in 2D!)","");
    }

    WeakPointerVector<Node<3> >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    unsigned int n_nodal_neighs = (unsigned int)neigh_nodes.size();
    const double h_inv = 1.0 / CalculateTheMaximumDistanceToNeighbours(p_node); // we use it as a scaling parameter to improve stability
    const array_1d <double, 3> origin = p_node->Coordinates();
    matrix<double> TestNodalValues(n_nodal_neighs, 1);
    matrix<double> A(n_nodal_neighs, n_poly_terms);

    for (unsigned int i = 0; i < n_nodal_neighs; ++i){
        A(i, 0) = 1.0;

        if (TDim == 3){
            Node<3>& neigh = neigh_nodes[i];
            const array_1d <double, 3> rel_coordinates = (neigh.Coordinates() - origin) * h_inv;
            TestNodalValues(i, 0) = SecondDegreeTestPolynomial(rel_coordinates);

            for (unsigned int d = 1; d < 10; ++d){
                if (d < 4){
                    A(i, d) = rel_coordinates[d - 1];
                }
                else if (d == 4){
                    A(i, d) = rel_coordinates[0] * rel_coordinates[1];
                }
                else if (d == 5){
                    A(i, d) = rel_coordinates[0] * rel_coordinates[2];
                }
                else if (d == 6){
                    A(i, d) = rel_coordinates[1] * rel_coordinates[2];
                }
                else {
                    A(i, d) = rel_coordinates[d - 7] * rel_coordinates[d - 7];
                }
            }
        }

        else {
            KRATOS_THROW_ERROR(std::runtime_error,"Gradient recovery not implemented yet in 2D!)","");
        }
    }

    matrix<double>AtransA(n_poly_terms, n_poly_terms);
    noalias(AtransA) = prod(trans(A), A);

    if (fabs(mMyCustomFunctions.template determinant< matrix<double> >(AtransA)) < 0.01){
        return false;
    }

    else {

        unsigned int n_relevant_terms = 0;

        if (mCalculatingTheGradient){
            n_relevant_terms = TDim;
        }

        else if (mCalculatingTheLaplacian){
            n_relevant_terms = n_poly_terms - (TDim + 1);
        }

        std::vector<unsigned int> relevant_terms;
        relevant_terms.resize(n_relevant_terms);
        double normalization =  1.0;

        if (mCalculatingTheGradient){
            normalization = h_inv;
            relevant_terms[0] = 1;
            relevant_terms[1] = 2;
            relevant_terms[2] = 3;
        }

        else if (mCalculatingTheLaplacian){
            normalization = h_inv * h_inv;
            relevant_terms[0] = 4;
            relevant_terms[1] = 5;
            relevant_terms[2] = 6;
            relevant_terms[3] = 7;
            relevant_terms[4] = 8;
            relevant_terms[5] = 9;
        }

        Vector& nodal_weights = p_node->FastGetSolutionStepValue(NODAL_WEIGHTS);
        nodal_weights.resize(n_relevant_terms * n_nodal_neighs);
        matrix<double>AtransAinv(n_poly_terms, n_poly_terms);
        noalias(AtransAinv) = mMyCustomFunctions.Inverse(AtransA);
//        for (unsigned i = 0; i < n_poly_terms; i++){
//            for (unsigned j = 0; j < n_poly_terms; j++){
//                if (abs(AtransAinv(i,j))>1e6){
//                    return false;
//                }
//            }
//        }
        matrix<double>AtransAinvAtrans(n_poly_terms, n_nodal_neighs);
        noalias(AtransAinvAtrans) = prod(AtransAinv, trans(A));

        for (unsigned int i = 0; i < n_nodal_neighs; ++i){
            for (unsigned int d = 0; d < n_relevant_terms; ++d){
                nodal_weights(n_relevant_terms * i + d) = AtransAinvAtrans(relevant_terms[d], i) * normalization;
            }
        }

        matrix<double> C(n_nodal_neighs, 1);
        C = prod(AtransAinvAtrans, TestNodalValues);

        double abs_difference = 0.0;

        for (unsigned int i = 0; i < n_nodal_neighs; ++i){
            const array_1d <double, 3>& rel_coordinates = (neigh_nodes[i].Coordinates() - origin) * h_inv;
            abs_difference += fabs(SecondDegreeGenericPolynomial(C, rel_coordinates) - SecondDegreeTestPolynomial(rel_coordinates));
        }
//        abs_difference = abs(C(7, 0) - 1.0) + abs(C(8, 0) - 1.0) + abs(C(8, 0) - 1.0);
        const double tolerance = 0.001; // recommended by E Ortega

        return (abs_difference > tolerance? false : true);
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

unsigned int GetNumberOfUniqueNeighbours(const int my_id, const WeakPointerVector<Element>& my_neighbour_elements)
{
    std::vector<int> ids;
    ids.push_back(my_id);

    for (unsigned int i_el = 0; i_el < my_neighbour_elements.size(); ++i_el){
        const Geometry<Node<3> >& geom = my_neighbour_elements[i_el].GetGeometry();
        for (unsigned int jj = 0; jj < TDim + 1; ++jj){
            int id = (int)geom[jj].Id();
            std::vector<int>::iterator it;
            it = find(ids.begin(), ids.end(), id);

            if (it >= ids.end()){
                ids.push_back(id);
            }
        }
    }

    return((int)ids.size());
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateTheMaximumDistanceToNeighbours(Node<3>::Pointer& p_node)
{
    double max_distance_yet = 0.0;
    const array_1d <double, 3>& coors = p_node->Coordinates();
    WeakPointerVector<Node<3> >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);

    for (unsigned int i = 0; i < (unsigned int)neigh_nodes.size(); ++i){
        array_1d <double, 3> delta = neigh_nodes[i].Coordinates() - coors;
        double distance_2 = DEM_INNER_PRODUCT_3(delta, delta);
        max_distance_yet = max_distance_yet > distance_2 ? max_distance_yet : distance_2;
    }

    return(std::sqrt(max_distance_yet));
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateTheMaximumEdgeLength(ModelPart& r_model_part)
{
    double max_distance_yet = 0.0;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
        Geometry<Node<3> >& geom = ielem->GetGeometry();
        unsigned int n_nodes = static_cast<unsigned int>(TDim + 1);

        for (unsigned int k = 1; k < n_nodes - 1; ++k){
            for (unsigned int i = k; i < n_nodes; ++i){
                array_1d <double, 3> delta_i = geom[k - 1] - geom[i];
                double distance_2 = DEM_INNER_PRODUCT_3(delta_i, delta_i);
                max_distance_yet = max_distance_yet > distance_2 ? max_distance_yet : distance_2;
            }
        }
    }

    return(std::sqrt(max_distance_yet));
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateTheMinumumEdgeLength(ModelPart& r_model_part)
{
    double min_distance_yet = 0.0;

    bool first_node = true;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
        Geometry<Node<3> >& geom = ielem->GetGeometry();

        if (first_node){ // assign the distance (squared) between any two nodes to min_distance_yet
            array_1d <double, 3> delta = geom[0] - geom[1];
            double distance_2 = DEM_INNER_PRODUCT_3(delta, delta);
            min_distance_yet = distance_2;
        }

        unsigned int n_nodes = static_cast<unsigned int>(TDim + 1);

        for (unsigned int k = 1; k < n_nodes - 1; ++k){
            for (unsigned int i = k; i < n_nodes; ++i){
                array_1d <double, 3> delta_i = geom[k - 1] - geom[i];
                double distance_2 = DEM_INNER_PRODUCT_3(delta_i, delta_i);

                min_distance_yet = min_distance_yet < distance_2 ? min_distance_yet : distance_2;
            }
        }
    }

    return(std::sqrt(min_distance_yet));
}

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
