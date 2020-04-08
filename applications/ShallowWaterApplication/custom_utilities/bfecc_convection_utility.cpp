//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "bfecc_convection_utility.h"
#include "geometries/tetrahedra_3d_4.h"


namespace Kratos
{

template<std::size_t TDim>
BFECCConvectionUtility<TDim>::BFECCConvectionUtility(ModelPart& rThisModelPart, Parameters ThisParameters)
 : mrModelPart(rThisModelPart)
 , mSearchStructure(mrModelPart)
{
    Parameters default_parameters(R"({
        "maximum_results"   : 10000
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    mMaxResults = ThisParameters["maximum_results"].GetDouble();
    
    mSearchStructure.UpdateSearchDatabase();
}


template<std::size_t TDim>
void BFECCConvectionUtility<TDim>::UpdateSearchDatabase()
{
    mSearchStructure.UpdateSearchDatabase();
}


template<std::size_t TDim>
template<class TVarType, class TType>
void BFECCConvectionUtility<TDim>::Convect(const TVarType& rVar, const Variable<array_1d<double,3>>& rConvVar)
{
    const double dt = mrModelPart.GetProcessInfo()[DELTA_TIME];

    Vector N(TDim + 1);
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(mMaxResults);

    const int nparticles = mrModelPart.Nodes().size();

    PointerVector< Element > elem_backward(mrModelPart.Nodes().size());
    std::vector< Vector > Ns(mrModelPart.Nodes().size());
    std::vector< bool > found(mrModelPart.Nodes().size());

    // First loop: estimate rVar(n+1)
    #pragma omp parallel for firstprivate(results, N)
    for (int i = 0; i < nparticles; i++)
    {
        typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

        ModelPart::NodesContainerType::iterator i_particle = mrModelPart.NodesBegin() + i;

        Element::Pointer p_element;

        array_1d<double,3> bck_pos = i_particle->Coordinates();
        const array_1d<double,3>& vel = i_particle->FastGetSolutionStepValue(rConvVar);
        bool is_found = RK2Convect(dt, bck_pos, vel, N, p_element, result_begin, -1, rConvVar);
        found[i] = is_found;

        if(is_found) {
            //save position backwards
            elem_backward(i) = p_element;
            Ns[i] = N;

            Geometry< Node < 3 > >& geom = p_element->GetGeometry();
            TType phi1 = N[0] * ( geom[0].FastGetSolutionStepValue(rVar,1));
            for (unsigned int k = 1; k < geom.size(); k++) {
                phi1 += N[k] * ( geom[k].FastGetSolutionStepValue(rVar,1) );
            }

            i_particle->FastGetSolutionStepValue(rVar) = phi1;
        }
    }

    // Second loop: obtain the value at time step N by taking it from N+1
    #pragma omp parallel for firstprivate(results, N)
    for (int i = 0; i < nparticles; i++)
    {
        typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

        ModelPart::NodesContainerType::iterator i_particle = mrModelPart.NodesBegin() + i;

        Element::Pointer p_element;

        array_1d<double,3> fwd_pos = i_particle->Coordinates();
        const array_1d<double,3>& vel = i_particle->FastGetSolutionStepValue(rConvVar,1);
        bool is_found = RK2Convect(dt, fwd_pos, vel, N, p_element, result_begin, 1, rConvVar);

        if(is_found) {
            Geometry< Node < 3 > >& geom = p_element->GetGeometry();
            TType phi_old = N[0] * ( geom[0].FastGetSolutionStepValue(rVar));

            for (unsigned int k = 1; k < geom.size(); k++) {
                phi_old  += N[k] * ( geom[k].FastGetSolutionStepValue(rVar) );
            }

            // Store the correction
            i_particle->SetValue(rVar, 1.5*i_particle->FastGetSolutionStepValue(rVar,1) - 0.5 * phi_old);
        }
        else
        {
            i_particle->SetValue(rVar, i_particle->FastGetSolutionStepValue(rVar,1));
        }
    }

    // Third loop: apply the correction
    #pragma omp parallel for
    for (int i = 0; i < nparticles; i++)
    {
        ModelPart::NodesContainerType::iterator i_particle = mrModelPart.NodesBegin() + i;
        bool is_found = found[i];
        if(is_found) {
            Vector N = Ns[i];
            Geometry< Node < 3 > >& geom = elem_backward[i].GetGeometry();
            TType phi1 = N[0] * ( geom[0].GetValue(rVar));
            for (unsigned int k = 1; k < geom.size(); k++) {
                phi1 += N[k] * ( geom[k].GetValue(rVar) );
            }

            i_particle->FastGetSolutionStepValue(rVar) = phi1;
        }
    }
}


template<std::size_t TDim>
bool BFECCConvectionUtility<TDim>::RK2Convect(
    const double Dt,
    array_1d<double,3>& rPosition,
    const array_1d<double,3>& rInitialVelocity,
    Vector& rN,
    Element::Pointer& pElement,
    typename BinBasedFastPointLocator<TDim>::ResultIteratorType& rResultBegin,
    const int VelocitySign,
    const Variable<array_1d<double,3> >& rConvVar)
{
    bool is_found = false;
    array_1d<double,3> pos_step1 = rPosition + 0.5 * Dt * VelocitySign * rInitialVelocity;
    is_found = mSearchStructure.FindPointOnMesh(pos_step1, rN, pElement, rResultBegin, mMaxResults);
    if (is_found)
    {
        Geometry<Node<3>>& geom = pElement->GetGeometry();
        array_1d<double,3> vel_step1 = 0.5 * rN[0] * (geom[0].FastGetSolutionStepValue(rConvVar) + geom[0].FastGetSolutionStepValue(rConvVar,1));
        for (std::size_t i = 1; i < geom.size(); ++i) {
            vel_step1 += 0.5 * rN[i] * (geom[i].FastGetSolutionStepValue(rConvVar) + geom[i].FastGetSolutionStepValue(rConvVar,1));
        }
        rPosition += Dt * VelocitySign * vel_step1;
        is_found = mSearchStructure.FindPointOnMesh(rPosition, rN, pElement, rResultBegin, mMaxResults);
    }
    return is_found;
}


template<std::size_t TDim>
template<class TVarType>
void BFECCConvectionUtility<TDim>::ResetBoundaryConditions(const TVarType& rVar)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        if (it_node->IsFixed(rVar))
        {
            it_node->FastGetSolutionStepValue(rVar) = it_node->FastGetSolutionStepValue(rVar,1);
        }
    }
}


template<std::size_t TDim>
template<class TVarType>
void BFECCConvectionUtility<TDim>::CopyVariableToPreviousTimeStep(const TVarType& rVar)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(rVar,1) = it_node->FastGetSolutionStepValue(rVar);
    }
}


template class BFECCConvectionUtility<2>;
template class BFECCConvectionUtility<3>;

template void BFECCConvectionUtility<2>::Convect<Variable<double>, double>(const Variable<double>&, const Variable<array_1d<double,3>>&);
template void BFECCConvectionUtility<3>::Convect<Variable<double>, double>(const Variable<double>&, const Variable<array_1d<double,3>>&);

template void BFECCConvectionUtility<2>::Convect<Variable<array_1d<double,3>>, array_1d<double,3>>(const Variable<array_1d<double,3>>&, const Variable<array_1d<double,3>>&);
template void BFECCConvectionUtility<3>::Convect<Variable<array_1d<double,3>>, array_1d<double,3>>(const Variable<array_1d<double,3>>&, const Variable<array_1d<double,3>>&);

template void BFECCConvectionUtility<2>::ResetBoundaryConditions(const Variable<double>&);
template void BFECCConvectionUtility<2>::ResetBoundaryConditions(const VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>&);

template void BFECCConvectionUtility<2>::CopyVariableToPreviousTimeStep(const Variable<double>&);
template void BFECCConvectionUtility<2>::CopyVariableToPreviousTimeStep(const Variable<array_1d<double,3>>&);

}  // namespace Kratos