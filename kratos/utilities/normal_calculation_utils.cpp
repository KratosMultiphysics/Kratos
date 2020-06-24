//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

// System includes

// External includes

// Project includes
#include "utilities/normal_calculation_utils.h"

namespace Kratos
{

void NormalCalculationUtils::CalculateOnSimplex(
    ConditionsArrayType& rConditions,
    const std::size_t Dimension
    )
{
    KRATOS_TRY

    // Calculating the normals and storing on the conditions
    array_1d<double,3> An;
    if(Dimension == 2)  {
        for(ConditionsArrayType::iterator it =  rConditions.begin(); it !=rConditions.end(); it++) {
            if (it->GetGeometry().PointsNumber() == 2)
                CalculateNormal2D(it,An);
        }
    } else if(Dimension == 3) {
        array_1d<double,3> v1, v2;
        for(ConditionsArrayType::iterator it =  rConditions.begin(); it !=rConditions.end(); it++) {
            // Calculate the normal on the given condition
            if (it->GetGeometry().PointsNumber() == 3)
                CalculateNormal3D(it,An,v1,v2);
        }
    }

    // Adding the normals to the nodes
    for(ConditionsArrayType::iterator it =  rConditions.begin(); it !=rConditions.end(); it++) {
        Geometry<Node<3> >& pGeometry = (it)->GetGeometry();
        double coeff = 1.00/pGeometry.size();
        const array_1d<double,3>& r_normal = it->GetValue(NORMAL);
        for(unsigned int i = 0; i<pGeometry.size(); i++) {
            noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * r_normal;
        }
    }


    KRATOS_CATCH("")

}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplex(
    ModelPart& rModelPart,
    const std::size_t Dimension
    )
{
    // Resetting the normals
    const array_1d<double,3> zero = ZeroVector(3);

    if (rModelPart.GetCommunicator().GetDataCommunicator().IsDistributed()) {
        // If Parallel make sure normals are reset in all partitions
        VariableUtils().SetFlag(VISITED, false, rModelPart.Nodes());

        for(auto& r_cond : rModelPart.Conditions()) {
            for(auto& r_node: r_cond.GetGeometry()) {
                r_node.Set(VISITED, true);
            }
        }

        rModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);

        for(auto& r_node: rModelPart.Nodes()) {
            if(r_node.Is(VISITED)) {
                r_node.FastGetSolutionStepValue(NORMAL) = zero;
            }
        }
    } else {
        // In serial iteratre normally over the condition nodes
        for(auto& r_cond: rModelPart.Conditions()) {
            for(auto& r_node: r_cond.GetGeometry()) {
                r_node.FastGetSolutionStepValue(NORMAL) = zero;
            }
        }
    }

    const auto& r_process_info = rModelPart.GetProcessInfo();
    const bool has_domain_size = r_process_info.Has(DOMAIN_SIZE);
    KRATOS_ERROR_IF(has_domain_size && Dimension == 0) << "Dimension not defined" << std::endl;
    const SizeType dimension_in_model_part = has_domain_size ? r_process_info.GetValue(DOMAIN_SIZE) : Dimension;
    KRATOS_WARNING_IF("NormalCalculationUtils", dimension_in_model_part != Dimension) << "Inconsistency between DOMAIN_SIZE and Dimension provided" << std::endl;
    this->CalculateOnSimplex(rModelPart.Conditions(), dimension_in_model_part);
    rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::SwapNormals(ModelPart& rModelPart)
{
    KRATOS_TRY

    for(auto& r_cond : rModelPart.Conditions()) {
        GeometryType& r_geometry = r_cond.GetGeometry();
        Node<3>::Pointer paux = r_geometry(0);
        r_geometry(0) = r_geometry(1);
        r_geometry(1) = paux;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormal2D(
    ConditionsArrayType::iterator it,
    array_1d<double,3>& An
    )
{
    Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

    An[0] =    pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;

    array_1d<double,3>& normal = (it)->GetValue(NORMAL);
    noalias(normal) = An;
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormal3D(
    ConditionsArrayType::iterator it,
    array_1d<double,3>& An,
    array_1d<double,3>& v1,
    array_1d<double,3>& v2
    )
{
    Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;

    array_1d<double,3>& normal = (it)->GetValue(NORMAL);
    noalias(normal) = An;
}

} // namespace Kratos
