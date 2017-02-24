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

template<class TSpace, unsigned int TDim>
class PartitionedFSIUtilities
{

public:

    /** Type Definitions
    */

    /*@{ */
    typedef typename TSpace::VectorType                             VectorType;
    typedef typename TSpace::MatrixType                             MatrixType;

    // typedef typename TSpace::VectorPointerType               VectorPointerType;
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
        if (TDim == 2)
        {
            #pragma omp parallel for
            for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
                unsigned int base_i = k*2;

                array_1d<double,3> aux_velocity;
                aux_velocity[0] = rFluidInterfaceVelocity[base_i];
                aux_velocity[1] = rFluidInterfaceVelocity[base_i+1];
                aux_velocity[2] = 0.0;

                it_node->Fix(VELOCITY_X);
                it_node->Fix(VELOCITY_Y);
                it_node->Fix(VELOCITY_Z);

                noalias(it_node->FastGetSolutionStepValue(VELOCITY)) = aux_velocity;
            }
        }
        else
        {
            # pragma omp parallel for
            for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
                unsigned int base_i = k*3;

                array_1d<double,3> aux_velocity;
                aux_velocity[0] = rFluidInterfaceVelocity[base_i];
                aux_velocity[1] = rFluidInterfaceVelocity[base_i+1];
                aux_velocity[2] = rFluidInterfaceVelocity[base_i+2];

                it_node->Fix(VELOCITY_X);
                it_node->Fix(VELOCITY_Y);
                it_node->Fix(VELOCITY_Z);

                noalias(it_node->FastGetSolutionStepValue(VELOCITY)) = aux_velocity;
            }
        }
    };

    // Set the fluid interface nodal force
    void SetFluidInterfaceNodalFlux(VectorType& rFluidInterfaceNodalForce)
    {
        if (TDim == 2)
        {
            #pragma omp parallel for
            for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
                unsigned int base_i = k*2;

                array_1d<double,3> aux_force;
                aux_force[0] = rFluidInterfaceNodalForce[base_i];
                aux_force[1] = rFluidInterfaceNodalForce[base_i+1];
                aux_force[2] = 0.0;

                it_node->Fix(FORCE_X);
                it_node->Fix(FORCE_Y);
                it_node->Fix(FORCE_Z);

                noalias(it_node->FastGetSolutionStepValue(FORCE)) = aux_force;
            }
        }
        else
        {
            #pragma omp parallel for
            for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
                unsigned int base_i = k*3;

                array_1d<double,3> aux_force;
                aux_force[0] = rFluidInterfaceNodalForce[base_i];
                aux_force[1] = rFluidInterfaceNodalForce[base_i+1];
                aux_force[2] = rFluidInterfaceNodalForce[base_i+2];

                it_node->Fix(FORCE_X);
                it_node->Fix(FORCE_Y);
                it_node->Fix(FORCE_Z);

                noalias(it_node->FastGetSolutionStepValue(FORCE)) = aux_force;
            }
        }
    };

    // Compute fluid interface velocity residual vector and store its norm in the ProcessInfo.
    // TODO: MPI parallelization
    VectorType ComputeFluidInterfaceVelocityResidual()
    {

        VectorType fluid_interface_residual = ZeroVector(this->GetFluidInterfaceResidualSize());

        // Initialize the interface residual
        VariableUtils().SetToZero_VectorVar(FSI_INTERFACE_RESIDUAL, mrFluidInterfaceModelPart.Nodes());

        // Compute the consistent residual TODO: Separate the condition consistent residual computation in a private method
        #pragma omp parallel for
        for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfConditions()); ++k)
        {
            ModelPart::ConditionIterator it_cond = mrFluidInterfaceModelPart.ConditionsBegin()+k;

            const Condition::GeometryType& rGeom = it_cond->GetGeometry();
            const unsigned int BlockSize = TDim;
            const unsigned int NumNodes = rGeom.PointsNumber();

            // Set the nodal values residual vector
            VectorType ResVect = ZeroVector(BlockSize*TDim);

            for (int j = 0; j < static_cast<int>(NumNodes); ++j)
            {
                array_1d<double, 3> velocity_fluid = rGeom[j].FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3> velocity_fluid_projected = rGeom[j].FastGetSolutionStepValue(VECTOR_PROJECTED);

                for (k = 0; k < static_cast<int>(TDim); ++k)
                {
                    ResVect[j*BlockSize+k] = velocity_fluid[k] - velocity_fluid_projected[k];
                }
            }

            // Compute the condition mass matrix
            MatrixType MassMat = ZeroMatrix(BlockSize*TDim, BlockSize*TDim);

            const Condition::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
            const unsigned int NumGauss = IntegrationPoints.size();
            const Matrix NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

            for (int g=0; g<static_cast<int>(NumGauss); ++g)
            {
                const Kratos::Vector& N = row(NContainer,g);
                const double DetJacGauss = rGeom.DeterminantOfJacobian(g, GeometryData::GI_GAUSS_2);
                const double GaussWeight = DetJacGauss * IntegrationPoints[g].Weight();

                unsigned int RowIndex = 0;
                unsigned int ColIndex = 0;

                for (unsigned int i=0; i<NumNodes; ++i)
                {
                    for (unsigned int j=0; j<NumNodes; ++j)
                    {
                        double Mij = GaussWeight * N[i] * N[j];

                        for (unsigned int d=0; d<TDim; d++)
                        {
                            MassMat(RowIndex+d,ColIndex+d) += Mij;
                        }

                        ColIndex += BlockSize;
                    }

                    RowIndex += BlockSize;
                    ColIndex = 0;
                }
            }

            // Accumulate the obtained consistent residual values
            VectorType ConsResVect(BlockSize*TDim);
            TSpace::Mult(MassMat, ResVect, ConsResVect);

            if (TDim == 2)
            {
                for (int i=0; i<static_cast<int>(NumNodes); ++i)
                {
                    array_1d<double, 3> aux_val;
                    aux_val[0] += ConsResVect[i*BlockSize];
                    aux_val[1] += ConsResVect[i*BlockSize+1];
                    it_cond->GetGeometry()[i].FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL) += aux_val;
                }
            }
            else
            {
                for (int i=0; i<static_cast<int>(NumNodes); ++i)
                {
                    array_1d<double, 3> aux_val;
                    aux_val[0] += ConsResVect[i*BlockSize];
                    aux_val[1] += ConsResVect[i*BlockSize+1];
                    aux_val[2] += ConsResVect[i*BlockSize+2];
                    it_cond->GetGeometry()[i].FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL) += aux_val;
                }
            }
        }

        // Assemble the final consistent residual values
        if (TDim == 2)
        {
            #pragma omp parallel for
            for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
                unsigned int base_i = k*TDim;

                fluid_interface_residual[base_i] = it_node->FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL_X);
                fluid_interface_residual[base_i+1] = it_node->FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL_Y);
            }
        }
        else
        {
            #pragma omp parallel for
            for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
                unsigned int base_i = k*TDim;

                fluid_interface_residual[base_i] = it_node->FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL_X);
                fluid_interface_residual[base_i+1] = it_node->FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL_Y);
                fluid_interface_residual[base_i+2] = it_node->FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL_Y);
            }
        }

        // Store the L2 norm of the error in the fluid process info
        mrFluidInterfaceModelPart.GetProcessInfo().GetValue(FSI_INTERFACE_RESIDUAL_NORM) = TSpace::TwoNorm(fluid_interface_residual);

        return fluid_interface_residual;
    };

    // Compute fluid interface mesh velocity residual norm (TODO: Compute the mesh residual in a consistent manner)
    double ComputeFluidInterfaceMeshVelocityResidualNorm()
    {

        double total_weight = 0.0;
        double fluid_mesh_interface_residual_norm = 0.0;

        for (ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node ++)
        {
            double err_j = 0.0;
            double err_node = 0.0;
            array_1d<double, 3> velocity_fluid = it_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3> mesh_velocity_fluid = it_node->FastGetSolutionStepValue(MESH_VELOCITY);

            for (int j=0; j<static_cast<int>(TDim); j++)
            {
                err_j = velocity_fluid[j] - mesh_velocity_fluid[j];
                err_node += std::pow(err_j, 2);
            }

            double weight = it_node->FastGetSolutionStepValue(NODAL_AREA);
            total_weight += weight;
            // fluid_mesh_interface_residual_norm += weight*err_node;
            fluid_mesh_interface_residual_norm += err_node;
        }

        mrFluidInterfaceModelPart.GetCommunicator().SumAll(total_weight);
        mrFluidInterfaceModelPart.GetCommunicator().SumAll(fluid_mesh_interface_residual_norm);

        // return std::sqrt(fluid_mesh_interface_residual_norm/total_weight);
        return std::sqrt(fluid_mesh_interface_residual_norm);

    }

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
