//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla
//

#if !defined(KRATOS_PARTITIONED_FSI_UTILITIES )
#define  KRATOS_PARTITIONED_FSI_UTILITIES


/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "fsi_application.h"
#include "includes/define.h"
#include "containers/array_1d.h"
#include "includes/model_part.h"
#include "includes/communicator.h"
#include "utilities/openmp_utils.h"


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

template<class TSpace, unsigned int TDim>
class PartitionedFSIUtilities
{

public:

    /** Type Definitions
    */

    /*@{ */
    typedef typename TSpace::VectorType                             VectorType;
    // typedef typename TSpace::MatrixType                             MatrixType;

    typedef typename TSpace::VectorPointerType               VectorPointerType;
    // typedef typename TSpace::MatrixPointerType               MatrixPointerType;

    //~ /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( PartitionedFSIUtilities );
    /*@} */

    /** Constructor.
     */

    /*@{ */
    PartitionedFSIUtilities(ModelPart& rFluidInterfaceModelPart,
                            ModelPart& rStructureInterfaceModelPart):
                            mrFluidInterfaceModelPart(rFluidInterfaceModelPart),
                            mrStructureInterfaceModelPart(rStructureInterfaceModelPart)
    {
    }
    /*@} */

    /** Copy constructor.
    */

    /*@{ */
    PartitionedFSIUtilities(const PartitionedFSIUtilities& Other);
    /*@{ */

    /** Destructor.
     */

    /*@{ */
    ~PartitionedFSIUtilities()
    {
    }

    /*@} */
    /**@name Public Operators*/
    /*@{ */

    // Get the fluid interface number of nodes
    unsigned int GetFluidInterfaceProblemSize()
    {
        return mrFluidInterfaceModelPart.NumberOfNodes();
    }

    // Get the fluid interface number of nodes
    unsigned int GetStructureInterfaceProblemSize()
    {
        return mrStructureInterfaceModelPart.NumberOfNodes();
    }

    // Get the fluid interface number of nodes
    unsigned int GetFluidInterfaceResidualSize()
    {
        return (mrFluidInterfaceModelPart.NumberOfNodes())*TDim;
    }

    // Get the fluid interface number of nodes
    unsigned int GetStructureInterfaceResidualSize()
    {
        return (mrStructureInterfaceModelPart.NumberOfNodes())*TDim;
    }

    // Set the fluid interface velocity
    void SetFluidInterfaceVelocity(VectorType& rFluidInterfaceVelocity)
    {
        unsigned int i = 0;

        if (TDim == 2)
        {
            for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
            {
                array_1d<double,3> aux_velocity;
                aux_velocity[0] = rFluidInterfaceVelocity[i];
                aux_velocity[1] = rFluidInterfaceVelocity[i+1];
                aux_velocity[2] = 0.0;

                it_node->Fix(VELOCITY_X);
                it_node->Fix(VELOCITY_Y);
                it_node->Fix(VELOCITY_Z);

                it_node->FastGetSolutionStepValue(VELOCITY) = aux_velocity;
                i += 2;
            }
        }
        else
        {
            for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
            {
                array_1d<double,3> aux_velocity;
                aux_velocity[0] = rFluidInterfaceVelocity[i];
                aux_velocity[1] = rFluidInterfaceVelocity[i+1];
                aux_velocity[2] = rFluidInterfaceVelocity[i+2];

                it_node->Fix(VELOCITY_X);
                it_node->Fix(VELOCITY_Y);
                it_node->Fix(VELOCITY_Z);

                it_node->FastGetSolutionStepValue(VELOCITY) = aux_velocity;
                i += 3;
            }
        }
    };

    // Set the fluid interface nodal force
    void SetFluidInterfaceNodalFlux(VectorType& rFluidInterfaceNodalForce)
    {
        unsigned int i = 0;

        if (TDim == 2)
        {
            for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
            {
                array_1d<double,3> aux_force;
                aux_force[0] = rFluidInterfaceNodalForce[i];
                aux_force[1] = rFluidInterfaceNodalForce[i+1];
                aux_force[2] = 0.0;

                it_node->Fix(FORCE_X);
                it_node->Fix(FORCE_Y);
                it_node->Fix(FORCE_Z);

                it_node->FastGetSolutionStepValue(FORCE) = aux_force;
                i += 2;
            }
        }
        else
        {
            for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
            {
                array_1d<double,3> aux_force;
                aux_force[0] = rFluidInterfaceNodalForce[i];
                aux_force[1] = rFluidInterfaceNodalForce[i+1];
                aux_force[2] = rFluidInterfaceNodalForce[i+2];

                it_node->Fix(FORCE_X);
                it_node->Fix(FORCE_Y);
                it_node->Fix(FORCE_Z);

                it_node->FastGetSolutionStepValue(FORCE) = aux_force;
                i += 3;
            }
        }
    };

    // Compute fluid interface velocity residual
    VectorType ComputeFluidInterfaceVelocityResidual()
    {
        // Compute the residual
        unsigned int residual_size = this->GetFluidInterfaceResidualSize();
        double err_j = 0.0;
        double total_weight = 0.0;
        double fluid_interface_residual_norm = 0.0;
        VectorType fluid_interface_residual(residual_size);

        unsigned int i = 0;

        for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
        {
            double err_node = 0.0;
            array_1d<double,3> velocity_fluid = it_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double,3> velocity_fluid_projected = it_node->FastGetSolutionStepValue(VECTOR_PROJECTED);

            for (int j=0; j<static_cast<int>(TDim); j++)
            {
                err_j = velocity_fluid[j] - velocity_fluid_projected[j];
                fluid_interface_residual[i+j] = err_j;
                err_node += std::pow(err_j,2);
            }

            double weight = it_node->GetValue(NODAL_AREA);
            total_weight += weight;
            fluid_interface_residual_norm += weight*err_j;

            i += TDim;
        }

        mrFluidInterfaceModelPart.GetCommunicator().SumAll(total_weight);
        mrFluidInterfaceModelPart.GetCommunicator().SumAll(fluid_interface_residual_norm);

        // Store its weighted L2 norm in the fluid process info
        mrFluidInterfaceModelPart.GetProcessInfo().GetValue(FSI_INTERFACE_RESIDUAL_NORM) = fluid_interface_residual_norm/total_weight;

        return fluid_interface_residual;
    };

    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    ModelPart&      mrFluidInterfaceModelPart;
    ModelPart&      mrStructureInterfaceModelPart;

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
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

    //~ ModelPart& mr_model_part;

    //~ bool mMoveMeshFlag;


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

}; /* Class PartitionedFSIUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_PARTITIONED_FSI_UTILITIES  defined */
