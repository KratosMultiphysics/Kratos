// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Anna Rehr
//  Co-author   :    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/find_nodal_neighbours_process.h"
#include "custom_processes/spr_error_process.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
template<SizeType TDim>
SPRErrorProcess<TDim>::SPRErrorProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ): mThisModelPart(rThisModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "stress_vector_variable" : "CAUCHY_STRESS_VECTOR",
        "echo_level"             : 0
    })"
    );

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mStressVariable = KratosComponents<Variable<Vector>>::Get(ThisParameters["stress_vector_variable"].GetString());
    mEchoLevel = ThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SPRErrorProcess<TDim>::Execute()
{
    // Getting process info
    ProcessInfo::Pointer p_process_info = mThisModelPart.pGetProcessInfo();

    // Some initializations
    VariableUtils().SetNonHistoricalVariable(ELEMENT_ERROR, 0.0, mThisModelPart.Elements());
    VariableUtils().SetNonHistoricalVariable(ELEMENT_H, 0.0, mThisModelPart.Elements());

    /************************************************************************
    --1-- Calculate superconvergent stresses (at the nodes) --1--
    ************************************************************************/
    CalculateSuperconvergentStresses();

    /******************************************************************************
    --2-- Calculate error estimation --2--
    ******************************************************************************/
    double energy_norm_overall = 0.0;
    double error_overall = 0.0;
    CalculateErrorEstimation(energy_norm_overall, error_overall);

    // Final calculations
    const double tolerance = std::numeric_limits<double>::epsilon();
    const double denominator = std::sqrt(std::pow(error_overall, 2) + std::pow(energy_norm_overall , 2));
    const double coeff = denominator < tolerance ? 1.0 : 1.0/denominator;
    KRATOS_WARNING_IF("SPRErrorProcess", denominator < tolerance) << "Denominator of error estimate zero or almost zero " << denominator << std::endl;
    p_process_info->SetValue(ENERGY_NORM_OVERALL, energy_norm_overall);
    p_process_info->SetValue(ERROR_OVERALL, error_overall);
    p_process_info->SetValue(ERROR_RATIO, coeff * error_overall);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SPRErrorProcess<TDim>::CalculateSuperconvergentStresses()
{
    // We do a find of neighbours
    FindNodalNeighbours(mThisModelPart);

    // Iteration over all nodes -- construction of patches
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());

    #pragma omp parallel for
    for(int i_node = 0; i_node < num_nodes; ++i_node) {
        auto it_node = nodes_array.begin() + i_node;

        KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(NEIGHBOUR_ELEMENTS)) << "SPRErrorProcess:: Search didn't work with elements" << std::endl;
        const SizeType neighbour_size = it_node->GetValue(NEIGHBOUR_ELEMENTS).size();

        Vector sigma_recovered(SigmaSize, 0.0);

        if(neighbour_size > TDim) {
            CalculatePatch(it_node, it_node, neighbour_size,sigma_recovered);
            it_node->SetValue(RECOVERED_STRESS, sigma_recovered);

            KRATOS_INFO_IF("SPRErrorProcess", mEchoLevel > 2) << "Recovered sigma: " << sigma_recovered << std::endl;
        } else {
            KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(NEIGHBOUR_NODES)) << "SPRErrorProcess:: Search didn't work with nodes" << std::endl;
            auto& neigh_nodes = it_node->GetValue(NEIGHBOUR_NODES);
            for(auto it_neighbour_nodes = neigh_nodes.begin(); it_neighbour_nodes != neigh_nodes.end(); it_neighbour_nodes++) {

                Vector sigma_recovered_i(SigmaSize,0);

                IndexType count_i = 0;
                for(int i_node_loop = 0; i_node_loop < num_nodes; ++i_node_loop) { // FIXME: Avoid this double loop, extreamily expensive
                    auto it_node_loop = nodes_array.begin() + i_node_loop;
                    const SizeType size_elem_neigh = it_node_loop->GetValue(NEIGHBOUR_ELEMENTS).size();
                    if (it_node_loop->Id() == it_neighbour_nodes->Id() && size_elem_neigh > TDim){
                        CalculatePatch(it_node, it_node_loop, neighbour_size, sigma_recovered_i);
                        ++count_i;
                    }
                }

                // Average solution from different patches
                if(count_i != 0)
                    sigma_recovered = sigma_recovered*(count_i-1)/count_i + sigma_recovered_i/count_i;
            }

            it_node->SetValue(RECOVERED_STRESS,sigma_recovered);
            KRATOS_INFO_IF("SPRErrorProcess", mEchoLevel > 2) << "Recovered sigma: " << sigma_recovered << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SPRErrorProcess<TDim>::CalculateErrorEstimation(
    double& rEnergyNormOverall,
    double& rErrorOverall
    )
{
    // Loop over all elements:
    ElementsArrayType& elements_array = mThisModelPart.Elements();
    const int num_elem = static_cast<int>(elements_array.size());

    // Compute the error estimate per element
    double error_overall= 0.0;
    double energy_norm_overall = 0.0;
    #pragma omp parallel for reduction(+:error_overall, energy_norm_overall)
    for(int i_elem = 0; i_elem < num_elem; ++i_elem){
        auto it_elem = elements_array.begin() + i_elem;

        std::vector<double> error_integration_point;
        const auto& process_info = mThisModelPart.GetProcessInfo();
        it_elem->GetValueOnIntegrationPoints(ERROR_INTEGRATION_POINT, error_integration_point, process_info);

        KRATOS_INFO_IF("SPRErrorProcess", mEchoLevel > 2) << "Error GP:" << error_integration_point << std::endl;

        // We compute the error overall
        double error_energy_norm = 0.0;
        for(IndexType i = 0;i < error_integration_point.size();++i)
            error_energy_norm += error_integration_point[i];
        error_overall += error_energy_norm;
        error_energy_norm = std::sqrt(error_energy_norm);
        it_elem->SetValue(ELEMENT_ERROR, error_energy_norm);


        // We compute now the energy norm
        std::vector<double> strain_energy;
        it_elem->GetValueOnIntegrationPoints(STRAIN_ENERGY, strain_energy, process_info);

        double energy_norm = 0.0;
        for(IndexType i = 0;i < strain_energy.size(); ++i)
            energy_norm += 2.0 * strain_energy[i];
        energy_norm_overall += energy_norm;
        energy_norm= std::sqrt(energy_norm);

        KRATOS_INFO_IF("SPRErrorProcess", mEchoLevel > 2) << "Element Id:" << it_elem->Id() << ". Element error: " << error_energy_norm << ". Energy norm: " << energy_norm << std::endl;
    }

    rErrorOverall = std::sqrt(error_overall);
    rEnergyNormOverall = std::sqrt(energy_norm_overall);
    double error_percentage = rErrorOverall/std::sqrt((std::pow(rErrorOverall, 2) + std::pow(rEnergyNormOverall, 2)));

    KRATOS_INFO_IF("SPRErrorProcess", mEchoLevel > 1)
        << "Overall error norm: " << rErrorOverall << std::endl
        << "Overall energy norm: "<< rEnergyNormOverall << std::endl
        << "Error in percent: " << error_percentage << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SPRErrorProcess<TDim>::CalculatePatch(
    NodeItType itNode,
    NodeItType itPatchNode,
    SizeType NeighbourSize,
    Vector& rSigmaRecovered
    )
{
    // Triangle and tetrahedra have only one GP by default
    std::vector<Vector> stress_vector(1);
    std::vector<array_1d<double,3>> coordinates_vector(1);

    // Our interest is to assemble the system A and b to solve a local problem for the element and estimate the new element size
    BoundedMatrix<double, TDim + 1, TDim + 1> A = ZeroMatrix(TDim + 1,TDim + 1);
    BoundedMatrix<double, TDim + 1, SigmaSize> b = ZeroMatrix(TDim + 1,SigmaSize);
    BoundedMatrix<double, 1, TDim + 1> p_k;

    auto& neigh_elements = itPatchNode->GetValue(NEIGHBOUR_ELEMENTS);
    for( WeakElementItType it_elem = neigh_elements.begin(); it_elem != neigh_elements.end(); ++it_elem) {

        it_elem->GetValueOnIntegrationPoints(mStressVariable,stress_vector,mThisModelPart.GetProcessInfo());
        it_elem->GetValueOnIntegrationPoints(INTEGRATION_COORDINATES,coordinates_vector,mThisModelPart.GetProcessInfo());

        KRATOS_INFO_IF("SPRErrorProcess", mEchoLevel > 3)
        << "\tStress: " << stress_vector[0] << std::endl
        << "\tx: " << coordinates_vector[0][0] << "\ty: " << coordinates_vector[0][1] << "\tz_coordinate: " << coordinates_vector[0][2] << std::endl;

        BoundedMatrix<double, 1, SigmaSize> sigma;
        for(IndexType j = 0; j < SigmaSize; ++j)
            sigma(0,j) = stress_vector[0][j];
        p_k(0,0) = 1.0;
        p_k(0,1) = coordinates_vector[0][0] - itPatchNode->X();
        p_k(0,2) = coordinates_vector[0][1] - itPatchNode->Y();
        if(TDim == 3)
            p_k(0,3)=coordinates_vector[0][2] - itPatchNode->Z();

        // Finally we add the contributiosn to our local system (A, b)
        noalias(A) += prod(trans(p_k), p_k);
        noalias(b) += prod(trans(p_k), sigma);
    }

    double det;
    BoundedMatrix<double, TDim + 1, TDim + 1> invA = MathUtils<double>::InvertMatrix<TDim + 1>(A, det, -1.0); // We consider a negative tolerance in order to avoid error

    KRATOS_INFO_IF("SPRErrorProcess", mEchoLevel > 3) << A << std::endl << invA << std::endl << det<< std::endl;

    // We do a little correction trick in case of almost zero determinant
    const double tolerance = std::numeric_limits<double>::epsilon();
    if(det < tolerance){
        KRATOS_WARNING_IF("SPRErrorProcess", mEchoLevel == 2) << A << std::endl;
        for( IndexType i = 0; i < TDim + 1;i++){
            for( IndexType j = 0; j < TDim + 1; j++)
                A(i,j) += 0.001;
        }
        invA = MathUtils<double>::InvertMatrix<TDim + 1>(A,det);
        KRATOS_WARNING_IF("SPRErrorProcess", mEchoLevel > 0) << "det: " << det << std::endl;
    }

    const BoundedMatrix<double, TDim + 1, SigmaSize> coeff = prod(invA, b);

    if(NeighbourSize > TDim) {
        noalias(rSigmaRecovered) = row(coeff, 0);
    } else {
        p_k(0,1) = itNode->X() - itPatchNode->X();
        p_k(0,2) = itNode->Y() - itPatchNode->Y();
        if(TDim ==3)
            p_k(0,3) = itNode->Z() - itPatchNode->Z();
        const BoundedMatrix<double, 1, SigmaSize> sigma = prod(p_k,coeff);
        noalias(rSigmaRecovered) = row(sigma, 0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
inline void SPRErrorProcess<TDim>::FindNodalNeighbours(ModelPart& rModelPart)
{
    FindNodalNeighboursProcess find_neighbours(rModelPart);
    // Array of nodes
    NodesArrayType& nodes_array = rModelPart.Nodes();
    if (nodes_array.begin()->Has(NEIGHBOUR_ELEMENTS)) find_neighbours.ClearNeighbours();
    find_neighbours.Execute();
}

/***********************************************************************************/
/***********************************************************************************/

template class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SPRErrorProcess<2>;
template class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SPRErrorProcess<3>;

};// namespace Kratos.
