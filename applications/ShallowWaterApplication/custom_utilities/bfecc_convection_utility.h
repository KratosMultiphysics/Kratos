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


#ifndef KRATOS_BFECC_CONVECTION_UTILITY_H_INCLUDED
#define KRATOS_BFECC_CONVECTION_UTILITY_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/variables.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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

/// Explicit convection utility
/** Convection of scalars and vectors for shallow water equations using BFECC correction
*/
template<std::size_t TDim> class BFECCConvectionUtility
{

public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    KRATOS_CLASS_POINTER_DEFINITION(BFECCConvectionUtility<TDim>);

    ///@}
    ///@name Life Cycle
    ///@{

    BFECCConvectionUtility(ModelPart& rThisModelPart, Parameters ThisParameters = Parameters()) :
        mrModelPart(rThisModelPart),
        mSearchStructure(mrModelPart)
    {
        Parameters default_parameters(R"({
            "maximum_results"   : 10000
        })");

        ThisParameters.ValidateAndAssignDefaults(default_parameters);
        mMaxResults = ThisParameters["maximum_results"].GetDouble();

        mSearchStructure.UpdateSearchDatabase();
    }


    ~BFECCConvectionUtility() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method sets at Tn+1 the material value from Tn, using the velocity variable to convect
     * @param rVar The variable to convect
     * @param rConvVar The velocity variable
     */
    template<class TVarType, class TType>
    void Convect(const TVarType& rVar, const Variable<array_1d<double,3>>& rConvVar)
    {
        const double dt = mrModelPart.GetProcessInfo()[DELTA_TIME];
        const int n_particles = mrModelPart.Nodes().size();

        PointerVector< Element > elem_backward(mrModelPart.Nodes().size());
        std::vector< Vector > Ns(mrModelPart.Nodes().size());
        std::vector< bool > found(mrModelPart.Nodes().size());

        struct TLS {
            Vector N;
            typename BinBasedFastPointLocator<TDim>::ResultContainerType results;
        };
        TLS tls;
        tls.N.resize(TDim + 1);
        tls.results.resize(mMaxResults);

        IndexPartition<int>(n_particles).for_each(tls, [&](int i, TLS& rTLS){
            auto i_particle = mrModelPart.NodesBegin() + i;
            auto result_begin = rTLS.results.begin();
            Element::Pointer p_element;

            array_1d<double,3> bck_pos = i_particle->Coordinates();
            const array_1d<double,3>& vel = i_particle->FastGetSolutionStepValue(rConvVar);
            bool is_found = RK2Convect(dt, bck_pos, vel, rTLS.N, p_element, result_begin, -1, rConvVar);
            found[i] = is_found;

            if(is_found) {
                //save position backwards
                elem_backward(i) = p_element;
                Ns[i] = rTLS.N;

                Geometry< Node < 3 > >& geom = p_element->GetGeometry();
                TType phi1 = rTLS.N[0] * ( geom[0].FastGetSolutionStepValue(rVar,1));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += rTLS.N[k] * ( geom[k].FastGetSolutionStepValue(rVar,1) );
                }

                i_particle->FastGetSolutionStepValue(rVar) = phi1;
            }
        });

        // Second loop: obtain the value at time step N by taking it from N+1
        IndexPartition<int>(n_particles).for_each(tls, [&](int i, TLS& rTLS){
            auto i_particle = mrModelPart.NodesBegin() + i;
            auto result_begin = rTLS.results.begin();
            Element::Pointer p_element;

            array_1d<double,3> fwd_pos = i_particle->Coordinates();
            const array_1d<double,3>& vel = i_particle->FastGetSolutionStepValue(rConvVar,1);
            bool is_found = RK2Convect(dt, fwd_pos, vel, rTLS.N, p_element, result_begin, 1, rConvVar);

            if(is_found) {
                Geometry< Node < 3 > >& geom = p_element->GetGeometry();
                TType phi_old = rTLS.N[0] * ( geom[0].FastGetSolutionStepValue(rVar));

                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi_old  += rTLS.N[k] * ( geom[k].FastGetSolutionStepValue(rVar) );
                }

                // Store the correction
                i_particle->SetValue(rVar, 1.5*i_particle->FastGetSolutionStepValue(rVar,1) - 0.5 * phi_old);
            }
            else
            {
                i_particle->SetValue(rVar, i_particle->FastGetSolutionStepValue(rVar,1));
            }
        });

        // Third loop: apply the correction
        IndexPartition<int>(n_particles).for_each([&](int i){
            auto i_particle = mrModelPart.NodesBegin() + i;
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
        });
    }

    /**
     * @brief This method updates the search structure if the mesh has been modified
     */
    void UpdateSearchDatabase()
    {
        mSearchStructure.UpdateSearchDatabase();
    }

    /**
     * @brief This method copy the value from Tn to Tn+1 if the variable is fixed
     * @param The variable to reset if is fixed
     */
    template<class TVarType>
    void ResetBoundaryConditions(const TVarType& rVar)
    {
        block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
            if (rNode.IsFixed(rVar))
                rNode.FastGetSolutionStepValue(rVar) = rNode.FastGetSolutionStepValue(rVar,1);
        });
    }

    /**
     * This method copies the variable from Tn+1 to Tn
     * @param rVar The variable to copy to the previous time step
     */
    template<class TVarType>
    void CopyVariableToPreviousTimeStep(const TVarType& rVar)
    {
        block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
            rNode.FastGetSolutionStepValue(rVar,1) = rNode.FastGetSolutionStepValue(rVar);
        });
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


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    BinBasedFastPointLocator<TDim> mSearchStructure;
    int mMaxResults;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    bool RK2Convect(
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

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

};

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_BFECC_CONVECTION_UTILITY_H_INCLUDED  defined