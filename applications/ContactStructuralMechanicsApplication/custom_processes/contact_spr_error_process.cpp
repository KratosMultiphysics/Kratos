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
#include "custom_processes/contact_spr_error_process.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
template<SizeType TDim>
ContactSPRErrorProcess<TDim>::ContactSPRErrorProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ): BaseType(rThisModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "stress_vector_variable"              : "CAUCHY_STRESS_VECTOR",
        "penalty_normal"                      : 1.0e4,
        "penalty_tangential"                  : 1.0e4,
        "echo_level"                          : 0
    })"
    );

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mPenaltyNormal = ThisParameters["penalty_normal"].GetDouble();
    mPenaltyTangent = ThisParameters["penalty_tangential"].GetDouble();

    BaseType::mStressVariable = KratosComponents<Variable<Vector>>::Get(ThisParameters["stress_vector_variable"].GetString());
    BaseType::mEchoLevel = ThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void ContactSPRErrorProcess<TDim>::CalculatePatch(
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
    CompressedMatrix A(SigmaSize * (TDim + 1), SigmaSize * (TDim + 1), 0.0);
    BoundedMatrix<double, SigmaSize * (TDim+1),1> b = ZeroMatrix(SigmaSize * (TDim+1), 1);
    BoundedMatrix<double, SigmaSize, SigmaSize * (TDim+1)> p_k;
    BoundedMatrix<double, 1, SigmaSize> N_k;
    BoundedMatrix<double, 1, SigmaSize> T_k1;
    BoundedMatrix<double, 1, SigmaSize> T_k2;  // in case of 3D: second tangential vector
    BoundedMatrix<double, SigmaSize, 1> sigma;

    /* Computation A and b */
    // PART 1: Contributions from the neighboring elements
    auto& neigh_elements = itPatchNode->GetValue(NEIGHBOUR_ELEMENTS);
    for( WeakElementItType it_elem = neigh_elements.begin(); it_elem != neigh_elements.end(); ++it_elem) {

        auto& process_info = BaseType::mThisModelPart.GetProcessInfo();
        it_elem->GetValueOnIntegrationPoints(BaseType::mStressVariable,stress_vector, process_info);
        it_elem->GetValueOnIntegrationPoints(INTEGRATION_COORDINATES,coordinates_vector, process_info);

        KRATOS_INFO_IF("ContactSPRErrorProcess", BaseType::mEchoLevel > 3)
        << "\tElement: " << it_elem->Id() << std::endl
        << "\tStress: " << stress_vector[0] << std::endl
        << "\tX: " << coordinates_vector[0][0] << "\tY: " << coordinates_vector[0][1] << "\tZ: " << coordinates_vector[0][2] << std::endl;

        for( IndexType j = 0; j < SigmaSize; ++j)
            sigma(j,0) = stress_vector[0][j];

        for ( IndexType j = 0; j < SigmaSize; ++j){
            p_k(j, j * (TDim+1) + 0) = 1.0;
            p_k(j, j * (TDim+1) + 1) = coordinates_vector[0][0] - itPatchNode->X();
            p_k(j, j * (TDim+1) + 2) = coordinates_vector[0][1] - itPatchNode->Y();
            if(TDim == 3)
                p_k(j, j * (TDim+1) + 3) = coordinates_vector[0][2] - itPatchNode->Z();
        }

        // Finally we add the contributiosn to our local system (A, b)
        noalias(A) += prod(trans(p_k), p_k);
        noalias(b) += prod(trans(p_k), sigma);
    }

    // Computing A and b
    BoundedMatrix<double, SigmaSize * (TDim + 1), 1> A1;
    BoundedMatrix<double, 1, SigmaSize * (TDim + 1)> A2;
    for (IndexType j = 0; j < SigmaSize; ++j){
        p_k(j,j * (TDim + 1) + 1)= itNode->X() - itPatchNode->X();
        p_k(j,j * (TDim + 1) + 2)= itNode->Y() - itPatchNode->Y();
        if(TDim == 3)
            p_k(j,j * (TDim + 1) + 3)= itNode->Z() - itPatchNode->Z();
    }

    // Set the normal and tangential vectors in Voigt Notation
    const array_1d<double, 3>& normal = itNode->GetValue(NORMAL);

    if(TDim == 2) {
        N_k(0,0) = normal[0] * normal[0];
        N_k(0,1) = normal[1] * normal[1];
        N_k(0,2) = 2.0 * normal[0] * normal[1];

        T_k1(0,0) =  normal[0] * normal[1];
        T_k1(0,1) = -normal[0] * normal[1];
        T_k1(0,2) =  normal[1] * normal[1] - normal[0] * normal[0];
    } else if (TDim == 3) {
        N_k(0,0) = normal[0] * normal[0];
        N_k(0,1) = normal[1] * normal[1];
        N_k(0,1) = normal[2] * normal[2];
        N_k(0,3) = 2.0 * normal[1] * normal[2];
        N_k(0,4) = 2.0 * normal[2] * normal[0];
        N_k(0,5) = 2.0 * normal[0] * normal[1];

        // Set tangential vectors
        array_1d<double, 3> t1(3, 0.0);
        array_1d<double, 3> t2(3, 0.0);

        const double tolerance = std::numeric_limits<double>::epsilon();
        if(std::abs(normal[0]) > tolerance || std::abs(normal[1]) > tolerance) {
            const double norm = std::sqrt((t1[0]*t1[0]+t1[1]*t1[1]));
            t1[0] = normal[1]/norm;
            t1[1] = -normal[0]/norm;

            t2[0] = -normal[0]*normal[2]/norm;
            t2[1] = -normal[1]*normal[2]/norm;
            t2[2] = normal[0]*normal[0]+normal[1]*normal[1]/norm;
        } else{
            t1[0] = 1.0;
            t2[1] = 1.0;
        }

        T_k1(0,0) = normal[0]*t1[0];
        T_k1(0,1) = normal[1]*t1[1];
        T_k1(0,2) = normal[2]*t1[2];
        T_k1(0,3) = normal[1]*t1[2]+normal[2]*t1[1];
        T_k1(0,4) = normal[2]*t1[0]+normal[0]*t1[2];
        T_k1(0,5) = normal[0]*t1[1]+normal[1]*t1[0];

        T_k2(0,0) = normal[0]*t2[0];
        T_k2(0,1) = normal[1]*t2[1];
        T_k2(0,2) = normal[2]*t2[2];
        T_k2(0,3) = normal[1]*t2[2]+normal[2]*t2[1];
        T_k2(0,4) = normal[2]*t2[0]+normal[0]*t2[2];
        T_k2(0,5) = normal[0]*t2[1]+normal[1]*t2[0];
    }

    noalias(A1) = prod(trans(p_k),trans(N_k));
    noalias(A2) = prod(N_k,p_k);
    noalias(A) += mPenaltyNormal * prod(A1, A2);

    noalias(A1) = prod(trans(p_k),trans(T_k1));
    noalias(A2) = prod(T_k1,p_k);

    // Finally we add the contributiosn to our local system (A, b)
    noalias(A) += mPenaltyTangent*prod(A1, A2);
    noalias(b) += mPenaltyNormal*prod(trans(p_k),trans(N_k)) * itNode->GetValue(CONTACT_PRESSURE);

    //PART 2: Contributions from contact nodes: regard all nodes from the patch which are in contact
    // Patch center node:
    if (itPatchNode->Has(CONTACT_PRESSURE)){
        const array_1d<double, 3>& normal_patch_node = itPatchNode->GetValue(NORMAL);
        p_k(0,1)=0.0;
        p_k(0,2)=0.0;
        p_k(1,4)=0.0;
        p_k(1,5)=0.0;
        p_k(2,7)=0.0;
        p_k(2,8)=0.0;
        N_k(0,0) = normal_patch_node[0]*normal_patch_node[0];
        N_k(0,1) = normal_patch_node[1]*normal_patch_node[1];
        N_k(0,2) = 2*normal_patch_node[0]*normal_patch_node[1];
        T_k1(0,0) = normal_patch_node[0]*normal_patch_node[1];
        T_k1(0,1) = -normal_patch_node[0]*normal_patch_node[1];
        T_k1(0,2) = normal_patch_node[1]*normal_patch_node[1]-normal_patch_node[0]*normal_patch_node[0];

        noalias(A1) = prod(trans(p_k),trans(N_k));
        noalias(A2) = prod(N_k,p_k);
        noalias(A) += mPenaltyNormal*prod(A1, A2);

        noalias(A1) = prod(trans(p_k),trans(T_k1));
        noalias(A2) = prod(T_k1,p_k);

        // Finally we add the contributiosn to our local system (A, b)
        noalias(A) += mPenaltyTangent*prod(A1, A2);
        // NOTE: Code commented. Could be of interest in the future
//         noalias(A) += mPenaltyNormal*prod(prod(trans(p_k),trans(N_k)),prod(N_k,p_k));
//         noalias(A) += mPenaltyTangent*prod(prod(prod(trans(p_k),trans(T_k1)),T_k1),p_k);
        noalias(b) -= mPenaltyNormal*prod(trans(p_k),trans(N_k))*itPatchNode->GetValue(CONTACT_PRESSURE);
    }

    // Neighboring nodes:
    for( auto& i_neighbour_node : itPatchNode->GetValue(NEIGHBOUR_NODES)) {
        if (i_neighbour_node.Has(CONTACT_PRESSURE)){
            const array_1d<double, 3>& normal_neigh_node = i_neighbour_node.GetValue(NORMAL);
            const double x_patch = itPatchNode->X();
            const double y_patch = itPatchNode->Y();
            const double x_neigh = i_neighbour_node.X();
            const double y_neigh = i_neighbour_node.Y();
            p_k(0,1)= x_neigh-x_patch;
            p_k(0,2)= y_neigh-y_patch;
            p_k(1,4)= x_neigh-x_patch;
            p_k(1,5)= y_neigh-y_patch;
            p_k(2,7)= x_neigh-x_patch;
            p_k(2,8)= y_neigh-y_patch;
            N_k(0,0) = normal_neigh_node[0]*normal_neigh_node[0];
            N_k(0,1) = normal_neigh_node[1]*normal_neigh_node[1];
            N_k(0,2) = 2*normal_neigh_node[0]*normal_neigh_node[1];
            T_k1(0,0) = normal_neigh_node[0]*normal_neigh_node[1];
            T_k1(0,1) = -normal_neigh_node[0]*normal_neigh_node[1];
            T_k1(0,2) = normal_neigh_node[1]*normal_neigh_node[1]-normal_neigh_node[0]*normal_neigh_node[0];

            noalias(A1) = prod(trans(p_k),trans(N_k));
            noalias(A2) = prod(N_k,p_k);
            noalias(A) += mPenaltyNormal*prod(A1, A2);

            noalias(A1) = prod(trans(p_k),trans(T_k1));
            noalias(A2) = prod(T_k1,p_k);

            // Finally we add the contributiosn to our local system (A, b)
            noalias(A) += mPenaltyTangent*prod(A1, A2);
            noalias(b) += mPenaltyNormal*prod(trans(p_k),trans(N_k))*i_neighbour_node.GetValue(CONTACT_PRESSURE);
        }
    }

    // Computing coefficients a: A*a=b
    KRATOS_INFO_IF("ContactSPRErrorProcess", BaseType::mEchoLevel > 3) << A << std::endl;

    Vector coeff(SigmaSize * (TDim+1));
    const Vector& b_vector = row(b, 0);
    MathUtils<double>::Solve(A, coeff, b_vector);

    for (IndexType j = 0; j < SigmaSize;++j){
        p_k(j,j*(TDim + 1) + 1) = itNode->X() - itPatchNode->X();
        p_k(j,j*(TDim + 1) + 2) = itNode->Y() - itPatchNode->Y();
        if (TDim == 3)
            p_k(j,j*(TDim + 1) + 3) = itNode->Z() - itPatchNode->Z();
    }

    BoundedMatrix<double, SigmaSize*(TDim + 1), 1> coeff_matrix;
    for (IndexType i=0; i<SigmaSize*(TDim + 1); ++i)
        coeff_matrix(i, 0) = coeff[i];

    // Asssign the stress recovered
    noalias(sigma) = prod(p_k, coeff_matrix);
    noalias(rSigmaRecovered) = row(sigma,0);

    KRATOS_INFO_IF("ContactSPRErrorProcess", BaseType::mEchoLevel > 1) <<" Recovered pressure: "<< prod(N_k,sigma) <<", LM: "<<itNode->GetValue(CONTACT_PRESSURE)<<std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template class ContactSPRErrorProcess<2>;
template class ContactSPRErrorProcess<3>;

};// namespace Kratos.
