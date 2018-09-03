//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla
//

#if !defined( KRATOS_NODAL_UPDATE_UTILITIES )
#define  KRATOS_NODAL_UPDATE_UTILITIES


/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/mesh_moving_variables.h"
#include "includes/fsi_variables.h"
#include "containers/array_1d.h"
#include "includes/model_part.h"
#include "includes/communicator.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */

/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */
/**@name  Enum's */
/*@{ */

/*@} */
/**@name  Functions */
/*@{ */

/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.
Detail class definition.
*/
template <unsigned int TDim>
class NodalUpdateBaseClass
{

public:

    /** Type Definitions
    */

    /*@{ */

    //~ /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( NodalUpdateBaseClass );
    /*@} */

    /** Constructor.
     */

    /**
    * Empty constructor
    */
    NodalUpdateBaseClass() {}

    /*@} */

    /** Copy constructor.
    */

    /*@{ */
    NodalUpdateBaseClass(const NodalUpdateBaseClass& Other);
    /*@{ */

    /** Destructor.
     */

    /*@{ */
    virtual ~NodalUpdateBaseClass() {}

    /*@} */
    /**@name Public Operators*/
    /*@{ */

    /**
     * Computes the displacement time derivatives according to the computed displacement values.
     * @param rInterfaceModelPart: modelpart in where the nodal update is to be performed
     * @param timeStep: time step value
     */
    virtual void UpdateMeshTimeDerivatives(ModelPart& rModelPart,
                                           const double timeStep) {
        KRATOS_ERROR << "Calling the nodal update base class UpdateMeshTimeDerivatives() method. Call the proper time scheme derived one.";
    }

    /**
     * Sets the fluid interface time derivatives as the mesh displacement computed values.
     * @param rInterfaceModelPart: modelpart in where the nodal update is to be performed
     */
    virtual void SetMeshTimeDerivativesOnInterface(ModelPart& rInterfaceModelPart) {
        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;

            array_1d<double, 3>& v_node = it_node->FastGetSolutionStepValue(VELOCITY);  // Current step interface velocity
            noalias(v_node) = it_node->FastGetSolutionStepValue(MESH_VELOCITY);         // Set the current interface velocity as the mesh velocity;
        }

        rInterfaceModelPart.GetCommunicator().SynchronizeVariable(VELOCITY);
    }

    /*@} */

protected:

    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /**@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

private:

    /*@} */
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

}; /* Class NodalUpdateBaseClass */


/** Short class definition.
Detail class definition.
*/
template <unsigned int TDim>
class NodalUpdateNewmark : public NodalUpdateBaseClass<TDim>
{

public:

    /** Type Definitions
    */

    /*@{ */

    //~ /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( NodalUpdateNewmark );
    /*@} */

    /** Constructor.
     */

    /**
    * Empty constructor
    */
    NodalUpdateNewmark(const double BossakAlpha = -0.3) {
        const double bossak_f = 0.0;
        const double bossak_beta = 0.25;
        const double bossak_gamma = 0.5;

        mBossakBeta = std::pow((1.0 + bossak_f - BossakAlpha), 2) * bossak_beta;
        mBossakGamma = bossak_gamma + bossak_f - BossakAlpha;
    }


    /*@} */

    /** Copy constructor.
    */

    /*@{ */
    NodalUpdateNewmark(const NodalUpdateNewmark& Other);
    /*@{ */

    /** Destructor.
     */

    /*@{ */
    virtual ~NodalUpdateNewmark() {}

    /*@} */
    /**@name Public Operators*/
    /*@{ */

    /**
     * Computes the displacement time derivatives according to the computed displacement values using the Newmark (Bossak) scheme.
     * @param rModelPart: modelpart in where the nodal update is to be performed
     * @param timeStep: time step value
     */
    void UpdateMeshTimeDerivatives(ModelPart &rModelPart,
                                   const double timeStep)  override{

        auto& rLocalMesh = rModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k = 0; k < static_cast<int>(rLocalMesh.NumberOfNodes()); ++k) {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;

            const array_1d<double, 3>& umesh_n = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);  // Previous step mesh displacement
            const array_1d<double, 3>& vmesh_n = it_node->FastGetSolutionStepValue(MESH_VELOCITY, 1);      // Previous step mesh velocity
            const array_1d<double, 3>& amesh_n = it_node->FastGetSolutionStepValue(MESH_ACCELERATION, 1);  // Previous step mesh acceleration

            const array_1d<double, 3>& umesh_n1 = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);    // Current step mesh displacement
            array_1d<double, 3>& vmesh_n1 = it_node->FastGetSolutionStepValue(MESH_VELOCITY);              // Current step mesh velocity (to be updated)
            array_1d<double, 3>& amesh_n1 = it_node->FastGetSolutionStepValue(MESH_ACCELERATION);          // Current step mesh acceleration (to be updated)

            const double const_u = mBossakGamma / (timeStep * mBossakBeta);
            const double const_v = 1.0 - mBossakGamma / mBossakBeta;
            const double const_a = timeStep * (1.0 - mBossakGamma / (2.0 * mBossakBeta));

            for (unsigned int d=0; d<TDim; ++d) {
                vmesh_n1[d] = const_u * (umesh_n1[d] - umesh_n[d]) + const_v * vmesh_n[d] + const_a * amesh_n[d];
                amesh_n1[d] = (1.0 / (timeStep * mBossakGamma)) * (vmesh_n1[d] - vmesh_n[d]) - ((1 - mBossakGamma) / mBossakGamma) * amesh_n[d];
            }
        }

        rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
        rModelPart.GetCommunicator().SynchronizeVariable(MESH_ACCELERATION);

    }

    /*@} */

protected:

    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    double mBossakBeta;
    double mBossakGamma;

    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /**@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

private:

    /*@} */
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

}; /* Class NodalUpdateNewmark */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_PARTITIONED_FSI_UTILITIES  defined */
