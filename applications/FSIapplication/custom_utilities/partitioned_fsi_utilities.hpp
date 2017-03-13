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
#include "utilities/math_utils.h"
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
    typedef typename TSpace::VectorType                     VectorType;
    typedef typename TSpace::MatrixType                     MatrixType;

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

    /**
     * This function computes the fluid interface residual size as the
     * number of fluid interface nodes times the problem domain size.
     * @return the fluid domain residual size
     */
    unsigned int GetFluidInterfaceResidualSize()
    {
        return (mrFluidInterfaceModelPart.NumberOfNodes())*TDim;
    }

    /**
     * This function computes the structure interface residual size as the
     * number of structure interface nodes times the problem domain size.
     * @return the structure domain residual size
     */
    unsigned int GetStructureInterfaceResidualSize()
    {
        return (mrStructureInterfaceModelPart.NumberOfNodes())*TDim;
    }

    /**
     * This function returns the fluid interface length in 2D or the fluid interface area in 3D.
     * @return the fluid interface length
     */
    double GetFluidInterfaceArea()
    {
        double FluidInterfaceArea = 0.0;

        #pragma omp parallel for reduction(+:FluidInterfaceArea)
        for(int k=0; k < static_cast<int>(mrFluidInterfaceModelPart.NumberOfConditions()); ++k)
        {
            ModelPart::ConditionIterator it_cond = mrFluidInterfaceModelPart.ConditionsBegin()+k;

            const Condition::GeometryType& rGeom = it_cond->GetGeometry();
            FluidInterfaceArea += rGeom.Length();
        }

        mrFluidInterfaceModelPart.GetCommunicator().SumAll(FluidInterfaceArea);

        return FluidInterfaceArea;
    }

    /**
     * This function returns the structure interface length in 2D or the structure interface area in 3D.
     * @return the structure interface length
     */
    double GetStructureInterfaceArea()
    {
        double StructureInterfaceArea = 0.0;

        #pragma omp parallel for reduction(+:StructureInterfaceArea)
        for(int k=0; k < static_cast<int>(mrStructureInterfaceModelPart.NumberOfConditions()); ++k)
        {
            ModelPart::ConditionIterator it_cond = mrStructureInterfaceModelPart.ConditionsBegin()+k;

            const Condition::GeometryType& rGeom = it_cond->GetGeometry();
            StructureInterfaceArea += rGeom.Length();
        }

        mrStructureInterfaceModelPart.GetCommunicator().SumAll(StructureInterfaceArea);

        return StructureInterfaceArea;
    }

    /**
     * This function sets the variable data contained in a vector over the the
     * fluid interface. The variable can be fixed or not.
     * @param rVariable: variable to be set
     * @param FixVariable: decide wether the variable is fixed or not
     * @param rFluidInterfaceDataVector: vector containing the data values to be set
     */
    void SetAndFixFluidInterfaceVectorVariable(const Variable<array_1d<double, 3 > >& rVariable,
                                               const bool FixVariable,
                                               const VectorType& rFluidInterfaceDataVector)
    {
        // Initialize the variable value
        VariableUtils().SetToZero_VectorVar(rVariable, mrFluidInterfaceModelPart.Nodes());

        #pragma omp parallel for
        for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
        {
            ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
            unsigned int base_i = k*TDim;

            array_1d<double,3>& value_to_set = it_node->FastGetSolutionStepValue(rVariable);
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                value_to_set[jj] = rFluidInterfaceDataVector[base_i+jj];
            }
        }

        // If needed, apply fixity to rVariable
        if (FixVariable)
        {
            // Get the variable components
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;

            const std::string variable_name = rVariable.Name();
            const component_type varx = KratosComponents< component_type >::Get(variable_name+std::string("_X"));
            const component_type vary = KratosComponents< component_type >::Get(variable_name+std::string("_Y"));
            const component_type varz = KratosComponents< component_type >::Get(variable_name+std::string("_Z"));

            // Fix the variable components
            VariableUtils().ApplyFixity(varx, true, mrFluidInterfaceModelPart.Nodes());
            VariableUtils().ApplyFixity(vary, true, mrFluidInterfaceModelPart.Nodes());
            VariableUtils().ApplyFixity(varz, true, mrFluidInterfaceModelPart.Nodes());
        }
    }

    /**
     * This function computes the velocity residual vector over the fluid interface.
     * The velocity residual is defined as the fluid velocity value minus the velocity
     * value mapped from the structure (stored in the VECTOR_PROJECTED variable).
     * The nodal values of the residual are stored in the FSI_INTERFACE_RESIDUAL variable.
     * Besides, the norm of the residual vector is stored in the ProcessInfo using
     * the FSI_INTERFACE_RESIDUAL_NORM variable.
     * @param fluid_interface_residual: reference to the residual vector
     */
    void ComputeFluidInterfaceVelocityResidual(VectorType& fluid_interface_residual) // TODO: MPI parallelization
    {
        fluid_interface_residual = ZeroVector(this->GetFluidInterfaceResidualSize());

        // Compute node-by-node residual
        // this->ComputeNodeByNodeResidual(VELOCITY, VECTOR_PROJECTED, FSI_INTERFACE_RESIDUAL);

        // Compute consitent residual
        this->ComputeConsistentResidual(VELOCITY, VECTOR_PROJECTED, FSI_INTERFACE_RESIDUAL);

        // Assemble the final consistent residual values
        #pragma omp parallel for
        for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
            const unsigned int base_i = k*TDim;

            const array_1d<double,3>& fsi_res = it_node->FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL);
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                fluid_interface_residual[base_i+jj] = fsi_res[jj];
            }
        }

        // Store the L2 norm of the error in the fluid process info
        mrFluidInterfaceModelPart.GetProcessInfo().GetValue(FSI_INTERFACE_RESIDUAL_NORM) = TSpace::TwoNorm(fluid_interface_residual);

    };

    /**
     * This function computes the mesh velocity residual over the fluid interface.
     * The mesh velocity residual is defined as the fluid velocity value minus the
     * mesh velocity value.
     * The nodal values of the residual are stored in the FSI_INTERFACE_MESH_RESIDUAL variable.
     * Besides, the norm of the mesh residual vector is stored in the ProcessInfo using
     * the FSI_INTERFACE_MESH_RESIDUAL_NORM variable.
     */
    void ComputeFluidInterfaceMeshVelocityResidualNorm() // TODO: MPI parallelization
    {

        VectorType fluid_interface_mesh_residual = ZeroVector(this->GetFluidInterfaceResidualSize());

        // Compute node-by-node residual
        // this->ComputeNodeByNodeResidual(VELOCITY, MESH_VELOCITY, FSI_INTERFACE_MESH_RESIDUAL);

        // Compute consitent residual
        this->ComputeConsistentResidual(VELOCITY, MESH_VELOCITY, FSI_INTERFACE_MESH_RESIDUAL);

        // Assemble the final consistent residual values
        #pragma omp parallel for
        for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
            const unsigned int base_i = k*TDim;

            const array_1d<double,3>& fsi_mesh_res = it_node->FastGetSolutionStepValue(FSI_INTERFACE_MESH_RESIDUAL);
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                fluid_interface_mesh_residual[base_i+jj]   = fsi_mesh_res[jj];
            }
        }

        // Store the L2 norm of the error in the fluid process info
        mrFluidInterfaceModelPart.GetProcessInfo().GetValue(FSI_INTERFACE_MESH_RESIDUAL_NORM) = TSpace::TwoNorm(fluid_interface_mesh_residual);

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


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /**
     * This function computes the nodal error of a vector magnitude. The error is defined
     * as OriginalVariable minus ModifiedVariable.
     * @param rOriginalVariable: variable with the reference value
     * @param rModifiedVariable: variable with the computed vvalue
     * @param rErrorStorageVariable: variable to store the error nodal value
     */
    void ComputeNodeByNodeResidual(const Variable<array_1d<double, 3 > >& rOriginalVariable,
                                   const Variable<array_1d<double, 3 > >& rModifiedVariable,
                                   const Variable<array_1d<double, 3 > >& rErrorStorageVariable)
    {
        // Initialize the residual storage variable
        VariableUtils().SetToZero_VectorVar(rErrorStorageVariable, mrFluidInterfaceModelPart.Nodes());

        #pragma omp parallel for
        for(int k=0; k<static_cast<int>(mrFluidInterfaceModelPart.NumberOfNodes()); ++k)
        {
            ModelPart::NodeIterator it_node = mrFluidInterfaceModelPart.NodesBegin()+k;
            array_1d<double, 3>& rErrorStorage = it_node->FastGetSolutionStepValue(rErrorStorageVariable);

            const array_1d<double, 3>& velocity_fluid = it_node->FastGetSolutionStepValue(rOriginalVariable);
            const array_1d<double, 3>& velocity_fluid_projected = it_node->FastGetSolutionStepValue(rModifiedVariable);

            rErrorStorage = velocity_fluid - velocity_fluid_projected;
        }
    }

    /**
     * This function computes the nodal error of a vector magnitude in a consistent manner.
     * The error is defined as the integral over the interface of a tests function times
     * the difference between rOriginalVariable and rModifiedVariable. 
     * @param rOriginalVariable: variable with the reference value
     * @param rModifiedVariable: variable with the computed vvalue
     * @param rErrorStorageVariable: variable to store the error nodal value
     */
    void ComputeConsistentResidual(const Variable<array_1d<double, 3 > >& rOriginalVariable,
                                   const Variable<array_1d<double, 3 > >& rModifiedVariable,
                                   const Variable<array_1d<double, 3 > >& rErrorStorageVariable)
    {
        // Initialize the interface residual variable
        VariableUtils().SetToZero_VectorVar(rErrorStorageVariable, mrFluidInterfaceModelPart.Nodes());

        #pragma omp parallel for
        for(int k=0; k < static_cast<int>(mrFluidInterfaceModelPart.NumberOfConditions()); ++k)
        {
            ModelPart::ConditionIterator it_cond = mrFluidInterfaceModelPart.ConditionsBegin()+k;

            const Condition::GeometryType& rGeom = it_cond->GetGeometry();
            const unsigned int BlockSize = TDim;
            const unsigned int NumNodes = rGeom.PointsNumber();

            // Set the residual vector nodal values
            VectorType ResVect = ZeroVector(BlockSize*NumNodes);
            for (int jj = 0; jj < static_cast<int>(NumNodes); ++jj)
            {
                const array_1d<double, 3>& velocity_fluid = rGeom[jj].FastGetSolutionStepValue(rOriginalVariable);
                const array_1d<double, 3>& velocity_fluid_projected = rGeom[jj].FastGetSolutionStepValue(rModifiedVariable);

                for (int kk = 0; kk < static_cast<int>(TDim); ++kk)
                {
                    ResVect[jj*BlockSize+kk] = velocity_fluid[kk] - velocity_fluid_projected[kk];
                }
            }

            // Compute the condition mass matrix
            MatrixType MassMat = ZeroMatrix(BlockSize*NumNodes, BlockSize*NumNodes);
            const Condition::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
            const unsigned int NumGauss = IntegrationPoints.size();
            const Matrix NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
            VectorType JacGauss;
            rGeom.DeterminantOfJacobian(JacGauss, GeometryData::GI_GAUSS_2);

            for (int g=0; g<static_cast<int>(NumGauss); ++g)
            {
                const Kratos::Vector& N = row(NContainer,g);
                const double GaussWeight = JacGauss[g] * IntegrationPoints[g].Weight();

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
            VectorType ConsResVect(BlockSize*NumNodes);
            TSpace::Mult(MassMat, ResVect, ConsResVect);

            for (int ii=0; ii<static_cast<int>(NumNodes); ++ii)
            {
                array_1d<double, 3> aux_val;
                aux_val[0] = 0.0;
                aux_val[1] = 0.0;
                aux_val[2] = 0.0;
                for (unsigned int jj=0; jj<TDim; ++jj)
                {
                    aux_val[jj] = ConsResVect[ii*BlockSize+jj];
                }
                it_cond->GetGeometry()[ii].SetLock(); // So it is safe to write in the condition node in OpenMP
                it_cond->GetGeometry()[ii].FastGetSolutionStepValue(rErrorStorageVariable) += aux_val;
                it_cond->GetGeometry()[ii].UnSetLock(); // Free the condition node for other threads
            }
        }
    }

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
