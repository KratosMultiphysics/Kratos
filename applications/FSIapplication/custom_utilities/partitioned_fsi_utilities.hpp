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
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ale_variables.h"
#include "includes/fsi_variables.h"
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

    typedef typename TSpace::VectorPointerType              VectorPointerType;
    typedef typename TSpace::MatrixPointerType              MatrixPointerType;

    //~ /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( PartitionedFSIUtilities );
    /*@} */

    /** Constructor.
     */

    /**
    * Empty constructor
    */
    PartitionedFSIUtilities()
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
    virtual ~PartitionedFSIUtilities()
    {
    }

    /*@} */
    /**@name Public Operators*/
    /*@{ */

    /**
     * This function computes the interface residual size as the
     * number of interface nodes times the problem domain size.
     * @return the model part residual size
     */
    int GetInterfaceResidualSize(ModelPart& rInterfaceModelPart)
    {
        int local_number_of_nodes = (rInterfaceModelPart.GetCommunicator().LocalMesh().NumberOfNodes())*TDim;
        rInterfaceModelPart.GetCommunicator().SumAll(local_number_of_nodes);

        return local_number_of_nodes;
    }

    /**
     * This function returns the interface length in 2D or the interface area in 3D.
     * @param rInterfaceModelPart: interface modelpart in where the are is computed
     * @return the given modelpart interface length
     */
    double GetInterfaceArea(ModelPart& rInterfaceModelPart)
    {
        double interface_area = 0.0;

        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::ConditionIterator local_mesh_conditions_begin = rLocalMesh.ConditionsBegin();
        #pragma omp parallel for firstprivate(local_mesh_conditions_begin) reduction(+:interface_area)
        for(int k=0; k < static_cast<int>(rLocalMesh.NumberOfConditions()); ++k)
        {
            const ModelPart::ConditionIterator it_cond = local_mesh_conditions_begin+k;
            const Condition::GeometryType& rGeom = it_cond->GetGeometry();
            interface_area += rGeom.Length();
        }

        rInterfaceModelPart.GetCommunicator().SumAll(interface_area);

        return interface_area;
    }

    /**
     * This function sets the variable data contained in a vector over the the
     * fluid interface.
     * @param rInterfaceModelPart: interface modelpart in where the vector variable is set
     * @param rVariable: variable to be set
     * @param rInterfaceDataVector: vector containing the data values to be set
     */
    // void SetInterfaceVectorVariable(ModelPart& rInterfaceModelPart,
    //                                 const Variable<array_1d<double, 3 > >& rVariable,
    //                                 const VectorType& rInterfaceDataVector)
    // {
    //     // Initialize the variable value
    //     VariableUtils().SetToZero_VectorVar(rVariable, rInterfaceModelPart.Nodes());
    //
    //     #pragma omp parallel for
    //     for(int k=0; k<static_cast<int>(rInterfaceModelPart.NumberOfNodes()); ++k)
    //     {
    //         ModelPart::NodeIterator it_node = rInterfaceModelPart.NodesBegin()+k;
    //         unsigned int base_i = k*TDim;
    //
    //         array_1d<double,3>& value_to_set = it_node->FastGetSolutionStepValue(rVariable);
    //         for (unsigned int jj=0; jj<TDim; ++jj)
    //         {
    //             value_to_set[jj] = rInterfaceDataVector[base_i+jj];
    //         }
    //     }
    // }

    /**
     * This function sets the variable data contained in a vector over the the
     * fluid interface. The variable can be fixed or not.
     * @param rInterfaceModelPart: interface modelpart in where the vector variable is set
     * @param rVariable: variable to be set
     * @param rFluidInterfaceDataVector: vector containing the data values to be set
     * @param FixVariable: if true, fixes the variable in the fluid interface model part
     */
    // void SetAndFixInterfaceVectorVariable(ModelPart& rInterfaceModelPart,
    //                                       const Variable<array_1d<double, 3 > >& rVariable,
    //                                       const VectorType& rFluidInterfaceDataVector,
    //                                       const bool FixVariable)
    // {
    //     this->SetInterfaceVectorVariable(rInterfaceModelPart, rVariable, rFluidInterfaceDataVector);
    //
    //     // If needed, apply fixity to rVariable
    //     if (FixVariable)
    //     {
    //         // Get the variable components
    //         typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
    //
    //         const std::string variable_name = rVariable.Name();
    //         const component_type varx = KratosComponents< component_type >::Get(variable_name+std::string("_X"));
    //         const component_type vary = KratosComponents< component_type >::Get(variable_name+std::string("_Y"));
    //         const component_type varz = KratosComponents< component_type >::Get(variable_name+std::string("_Z"));
    //
    //         // Fix the variable components
    //         VariableUtils().ApplyFixity(varx, true, rInterfaceModelPart.Nodes());
    //         VariableUtils().ApplyFixity(vary, true, rInterfaceModelPart.Nodes());
    //         VariableUtils().ApplyFixity(varz, true, rInterfaceModelPart.Nodes());
    //     }
    // }

    virtual void SetUpInterfaceVector(ModelPart& rInterfaceModelPart,
                                      VectorPointerType& pInterfaceVector)
    {
        unsigned int ResidualSize = this->GetInterfaceResidualSize(rInterfaceModelPart);
        if ( TSpace::Size(*pInterfaceVector) != ResidualSize )
        {
            TSpace::Resize(pInterfaceVector,ResidualSize);
        }

        TSpace::SetToZero(*pInterfaceVector);
    }

    /**
     * This function computes (and stores in a vector) the residual of a vector variable over the fluid interface.
     * The residual is defined as the OriginalVariable value minus the ModifiedVariable value.
     * The nodal values of the residual are stored in the FSI_INTERFACE_RESIDUAL variable.
     * Besides, the norm of the residual vector is stored in the ProcessInfo using
     * the FSI_INTERFACE_RESIDUAL_NORM variable.
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @param interface_residual: reference to the residual vector
     */
    virtual void ComputeInterfaceVectorResidual(ModelPart& rInterfaceModelPart,
                                                const Variable<array_1d<double, 3 > >& rOriginalVariable,
                                                const Variable<array_1d<double, 3 > >& rModifiedVariable,
                                                VectorType& interface_residual)
    {
        TSpace::SetToZero(interface_residual);

        // Compute node-by-node residual
        this->ComputeNodeByNodeResidual(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, FSI_INTERFACE_RESIDUAL);

        // Compute consitent residual
        // this->ComputeConsistentResidual(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, FSI_INTERFACE_RESIDUAL);

        // Assemble the final consistent residual values
        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            const unsigned int base_i = k*TDim;

            const array_1d<double,3>& fsi_res = it_node->FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL);
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                this->SetLocalValue(interface_residual, base_i+jj, fsi_res[jj]);
            }
        }

        // Store the L2 norm of the error in the fluid process info
        rInterfaceModelPart.GetProcessInfo().GetValue(FSI_INTERFACE_RESIDUAL_NORM) = TSpace::TwoNorm(interface_residual);

    }

    /**
     * This function computes the mesh velocity residual over the fluid interface.
     * The mesh velocity residual is defined as the fluid velocity value minus the
     * mesh velocity value.
     * The nodal values of the residual are stored in the FSI_INTERFACE_MESH_RESIDUAL variable.
     * Besides, the norm of the mesh residual vector is stored in the ProcessInfo using
     * the FSI_INTERFACE_MESH_RESIDUAL_NORM variable.
     * @param rFluidInterfaceModelPart: interface modelpart in where the residual is computed
     */
    void ComputeFluidInterfaceMeshVelocityResidualNorm(ModelPart& rFluidInterfaceModelPart)
    {
        VectorPointerType pFluidInterfaceMeshResidual = TSpace::CreateEmptyVectorPointer();
        this->SetUpInterfaceVector(rFluidInterfaceModelPart, pFluidInterfaceMeshResidual);

        // Compute node-by-node residual
        this->ComputeNodeByNodeResidual(rFluidInterfaceModelPart, VELOCITY, MESH_VELOCITY, FSI_INTERFACE_MESH_RESIDUAL);

        // Compute consitent residual
        // this->ComputeConsistentResidual(rFluidInterfaceModelPart, VELOCITY, MESH_VELOCITY, FSI_INTERFACE_MESH_RESIDUAL);

        // Assemble the final consistent residual values
        auto& rLocalMesh = rFluidInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            const unsigned int base_i = k*TDim;

            const array_1d<double,3>& fsi_mesh_res = it_node->FastGetSolutionStepValue(FSI_INTERFACE_MESH_RESIDUAL);
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                this->SetLocalValue(*pFluidInterfaceMeshResidual, base_i+jj, fsi_mesh_res[jj]);
            }
        }

        // Store the L2 norm of the error in the fluid process info
        rFluidInterfaceModelPart.GetProcessInfo().GetValue(FSI_INTERFACE_MESH_RESIDUAL_NORM) = TSpace::TwoNorm(*pFluidInterfaceMeshResidual);

    }

    /**
     * Sets the values in the corrected guess vector in the selected nodal variable.
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @param rSolutionVariable: variable in where the corrected solution is to be stored
     * @param rCorrectedGuess: vector containing the interface corrected values
     */
    virtual void UpdateInterfaceValues(ModelPart& rInterfaceModelPart,
                                       const Variable<array_1d<double, 3 > >& rSolutionVariable,
                                       VectorType& rCorrectedGuess)
    {
        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            const unsigned int base_i = k*TDim;

            array_1d<double,3>& updated_value = it_node->FastGetSolutionStepValue(rSolutionVariable);
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                updated_value[jj] = this->GetLocalValue( rCorrectedGuess, base_i+jj );
            }
        }

        rInterfaceModelPart.GetCommunicator().SynchronizeVariable(rSolutionVariable);
    }

    /**
     * Computes the displacement time derivatives according to the corrected interface values.
     * Note that this is done using the Bossak formulaes.
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @param alphaBossak: Bossak scheme alpha coefficient
     * @param timeStep: time step value
     * @param rCorrectedGuess: vector containing the interface corrected values
     */
    virtual void ComputeCorrectedInterfaceDisplacementDerivatives(ModelPart& rInterfaceModelPart,
                                                                  const double alphaBossak,
                                                                  const double timeStep,
                                                                  VectorType& rCorrectedGuess)
    {
        const double gamma = 0.5*(1-2*alphaBossak);
        const double beta = ((1-alphaBossak)*(1-alphaBossak))/4;

        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            const unsigned int base_i = k*TDim;

            const array_1d<double, 3>& u_n = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1); // Previous step mesh displacement
            const array_1d<double, 3>& v_n = it_node->FastGetSolutionStepValue(VELOCITY, 1);          // Previous step velocity
            const array_1d<double, 3>& a_n = it_node->FastGetSolutionStepValue(ACCELERATION, 1);      // Previous step acceleration

            array_1d<double, 3> u_n1 = ZeroVector(3);
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                u_n1[jj] = this->GetLocalValue( rCorrectedGuess, base_i+jj );  // Current step displacement (taken from the corrected interface value)
            }

            array_1d<double, 3>& a_n1 = it_node->FastGetSolutionStepValue(ACCELERATION);  // Current step acceleration (computed with Bossak scheme)
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                a_n1[jj] = (u_n1[jj] - u_n[jj] - timeStep*v_n[jj] - (timeStep*timeStep)*(0.5-beta+beta*alphaBossak)*a_n[jj])/((timeStep*timeStep)*beta*(1-alphaBossak));
            }

            array_1d<double, 3>& v_n1 = it_node->FastGetSolutionStepValue(VELOCITY);            // Current step velocity (computed with Bossak scheme)
            array_1d<double, 3>& vmesh_n1 = it_node->FastGetSolutionStepValue(MESH_VELOCITY);   // Current step mesh velocity (equal to current step velocity)
            for (unsigned int jj=0; jj<TDim; ++jj)
            {
                const double updated_velocity = v_n[jj] + timeStep*(1-gamma)*a_n[jj] + timeStep*gamma*(1-alphaBossak)*a_n1[jj] + timeStep*gamma*alphaBossak*a_n[jj];
                v_n1[jj] = updated_velocity;
                vmesh_n1[jj] = updated_velocity;
            }

        }

        rInterfaceModelPart.GetCommunicator().SynchronizeVariable(VELOCITY);
        rInterfaceModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
        rInterfaceModelPart.GetCommunicator().SynchronizeVariable(ACCELERATION);

    }

    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    // ModelPart&      mrFluidInterfaceModelPart;
    // ModelPart&      mrStructureInterfaceModelPart;

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /**
     * This function computes the nodal error of a vector magnitude. The error is defined
     * as OriginalVariable minus ModifiedVariable.
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @param rOriginalVariable: variable with the reference value
     * @param rModifiedVariable: variable with the computed vvalue
     * @param rErrorStorageVariable: variable to store the error nodal value
     */
    void ComputeNodeByNodeResidual(ModelPart& rInterfaceModelPart,
                                   const Variable<array_1d<double, 3 > >& rOriginalVariable,
                                   const Variable<array_1d<double, 3 > >& rModifiedVariable,
                                   const Variable<array_1d<double, 3 > >& rErrorStorageVariable)
    {
        // Initialize the residual storage variable
        VariableUtils().SetToZero_VectorVar(rErrorStorageVariable, rInterfaceModelPart.GetCommunicator().LocalMesh().Nodes());

        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;

            array_1d<double, 3>& rErrorStorage = it_node->FastGetSolutionStepValue(rErrorStorageVariable);
            const array_1d<double, 3>& value_fluid = it_node->FastGetSolutionStepValue(rOriginalVariable);
            const array_1d<double, 3>& value_fluid_projected = it_node->FastGetSolutionStepValue(rModifiedVariable);

            rErrorStorage = value_fluid - value_fluid_projected;
        }
    }

    /**
     * This function computes the nodal error of a vector magnitude in a consistent manner.
     * The error is defined as the integral over the interface of a tests function times
     * the difference between rOriginalVariable and rModifiedVariable.
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @param rOriginalVariable: variable with the reference value
     * @param rModifiedVariable: variable with the computed vvalue
     * @param rErrorStorageVariable: variable to store the error nodal value
     */
    void ComputeConsistentResidual(ModelPart& rInterfaceModelPart,
                                   const Variable<array_1d<double, 3 > >& rOriginalVariable,
                                   const Variable<array_1d<double, 3 > >& rModifiedVariable,
                                   const Variable<array_1d<double, 3 > >& rErrorStorageVariable)
    {
        // Initialize the interface residual variable
        VariableUtils().SetToZero_VectorVar(rErrorStorageVariable, rInterfaceModelPart.Nodes());

        #pragma omp parallel for
        for(int k=0; k < static_cast<int>(rInterfaceModelPart.NumberOfConditions()); ++k)
        {
            ModelPart::ConditionIterator it_cond = rInterfaceModelPart.ConditionsBegin()+k;

            const Condition::GeometryType& rGeom = it_cond->GetGeometry();
            const unsigned int BlockSize = TDim;
            const unsigned int NumNodes = rGeom.PointsNumber();

            // Set the residual vector nodal values
            VectorType ResVect = ZeroVector(BlockSize*NumNodes);
            for (int jj = 0; jj < static_cast<int>(NumNodes); ++jj)
            {
                const array_1d<double, 3>& value_fluid = rGeom[jj].FastGetSolutionStepValue(rOriginalVariable);
                const array_1d<double, 3>& value_fluid_projected = rGeom[jj].FastGetSolutionStepValue(rModifiedVariable);

                for (int kk = 0; kk < static_cast<int>(TDim); ++kk)
                {
                    ResVect[jj*BlockSize+kk] = value_fluid[kk] - value_fluid_projected[kk];
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
                array_1d<double, 3> aux_val = ZeroVector(3);
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

    virtual void SetLocalValue(VectorType& rVector, int LocalRow, double Value) const
    {
        TSpace::SetValue(rVector,LocalRow,Value);
    }

    virtual double GetLocalValue(VectorType& rVector, int LocalRow) const
    {
        return TSpace::GetValue(rVector,LocalRow);
    }

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

}; /* Class PartitionedFSIUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_PARTITIONED_FSI_UTILITIES  defined */
