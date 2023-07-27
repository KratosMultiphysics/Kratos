//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_BOUNDARY_ROTATION_UTILITY
#define KRATOS_MPM_BOUNDARY_ROTATION_UTILITY

// system includes

// external includes

// kratos includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "utilities/coordinate_transformation_utilities.h"

// Application includes
#include "particle_mechanics_application_variables.h"

namespace Kratos {

///@addtogroup ParticleMechanicsApplication
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

/* A utility to rotate the local contributions of certain nodes to the system matrix,
which is required to apply slip conditions (roller-type support) in arbitrary directions to the boundary nodes.*/
template<class TLocalMatrixType, class TLocalVectorType>
class MPMBoundaryRotationUtility: public CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double> {
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of MPMBoundaryRotationUtility
	KRATOS_CLASS_POINTER_DEFINITION(MPMBoundaryRotationUtility);

	using CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double>::Rotate;

	typedef Node NodeType;

	typedef Geometry< Node > GeometryType;

	///@}
	///@name Life Cycle
	///@{

	/// Constructor.
	/** @param DomainSize Number of space dimensions (2 or 3)
	 * @param NumRowsPerNode Number of matrix or vector rows associated to each node. Displacement DOFs are assumed to be the first mDomainSize rows in each block of rows.
	 * @param rVariable Kratos variable used to flag nodes where local system contributions will be rotated. All nodes with rVariable != Zero will be rotated.
	 */
	MPMBoundaryRotationUtility(
        const unsigned int DomainSize,
		const unsigned int BlockSize,
		const Variable<double>& rVariable):
    CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double>(DomainSize,BlockSize,SLIP), mrFlagVariable(rVariable)
	{}

	/// Destructor.
	~MPMBoundaryRotationUtility() override {}

	/// Assignment operator.
	MPMBoundaryRotationUtility& operator=(MPMBoundaryRotationUtility const& rOther) {}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	/// Rotate the local system contributions so that they are oriented with each node's normal.
	/**
	 @param rLocalMatrix Local system matrix
	 @param rLocalVector Local RHS vector
	 @param rGeometry A reference to the element's (or condition's) geometry
	 */
	void Rotate(
        TLocalMatrixType& rLocalMatrix,
		TLocalVectorType& rLocalVector,
		GeometryType& rGeometry) const override
	{
		if (this->GetBlockSize() == this->GetDomainSize()) // irreducible case
		{
			if (this->GetDomainSize() == 2) this->template RotateAuxPure<2>(rLocalMatrix,rLocalVector,rGeometry);
			else if (this->GetDomainSize() == 3) this->template RotateAuxPure<3>(rLocalMatrix,rLocalVector,rGeometry);
		}
		else // mixed formulation case
		{
			if (this->GetDomainSize() == 2) this->template RotateAux<2,3>(rLocalMatrix,rLocalVector,rGeometry);
			else if (this->GetDomainSize() == 3) this->template RotateAux<3,4>(rLocalMatrix,rLocalVector,rGeometry);
		}

	}

	/// RHS only version of Rotate
	void RotateRHS(
        TLocalVectorType& rLocalVector,
		GeometryType& rGeometry) const
	{
		this->Rotate(rLocalVector,rGeometry);
	}

    // Auxiliary function to clear friction-related flags --
    // MUST be called before (re-)building the RHS in a given non-linear iteration
    static void ClearFrictionFlag(const ModelPart &rModelPart) {
        KRATOS_TRY
        // Loop over the grid nodes performed to clear INLET flag for indicating that nodal friction has alr been set
        // and remove the accumulated normal forces
        for(NodeType &curr_node : rModelPart.Nodes()){
            curr_node.SetLock();
            curr_node.Reset(INLET);
            curr_node.Reset(OUTLET);
            curr_node.FastGetSolutionStepValue(NORMAL_REACTION, 0) = 0.0;
            curr_node.UnSetLock();
        }
        KRATOS_CATCH( "" )
    }

	/// Apply roler type boundary conditions to the rotated local contributions.
	/** This function takes the rotated local system contributions so each
	 node's displacement are expressed using a base oriented with its normal
	 and imposes that the normal displacement is equal to the mesh displacement in
	 the normal direction.
	 */
	void ApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const override
	{
		const unsigned int LocalSize = rLocalVector.size();

		if (LocalSize > 0)
		{
			for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
			{
				if(this->IsSlip(rGeometry[itNode]) )
				{
					// We fix the first displacement dof (normal component) for each rotated block
					unsigned int j = itNode * this->GetBlockSize();

					// Get the displacement of the boundary mesh, this does not assume that the mesh is moving.
					// If the mesh is moving, need to consider the displacement of the moving mesh into account.
					const array_1d<double,3> & displacement = rGeometry[itNode].FastGetSolutionStepValue(DISPLACEMENT);

					// Get Normal Vector of the boundary
					array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
					this->Normalize(rN);

                    // Zero out row/column corresponding to normal displacement DoF except diagonal term (set to 1)
                    // Applied IFF the local matrix passed is not empty [otherwise does nothing -- RHS only case]
                    if (rLocalMatrix.size1() != 0) {
                        for( unsigned int i = 0; i < LocalSize; ++i)
                        {
                            rLocalMatrix(i,j) = 0.0;
                            rLocalMatrix(j,i) = 0.0;
                        }
                        rLocalMatrix(j, j) = 1.0; // set diagonal term to 1.0
                    }

                    /// Computation of nodal reaction forces due to conforming SLIP BC -- use NORMAL_REACTION
                    /// as frame of reference is aligned with node normals [& not global frame of ref per REACTION]
                    // Accumulate RHS values along normal direction to FORCE_RESIDUAL [nodal reaction forces due to conforming SLIP]
                    // -- when converged, RHS ~= 0 -> FORCE_RESIDUAL = -RHS value [computed without adding reaction force]

                    rGeometry[itNode].SetLock();
                    rGeometry[itNode].FastGetSolutionStepValue(NORMAL_REACTION) -= rLocalVector[j];
                    rGeometry[itNode].UnSetLock();

                    // Set value of normal displacement at node directly to the normal displacement of the boundary mesh
					rLocalVector[j] = inner_prod(rN,displacement);

                    // TODO: refactor? [e.g. set INLET flag and local flag, unlock node and do another if based on local flag]
                    // Prescribe a constant force of FRICTION_FORCE Newtons in the 1st tangential direction to all nodes
                    // INLET flag (reset in MPMResidualBasedBossakScheme.InitializeSolutionStep) ensures that friction
                    // is applied exactly once per node
                    rGeometry[itNode].SetLock();
                    if(!rGeometry[itNode].Is(INLET)) {
                        // obtain nodal velocity and rotate it to the same frame of reference as the local geometry
                        array_1d<double, 3> nodal_velocity = ZeroVector(3);

                        nodal_velocity = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);

                        // TODO: make threshold configurable?
                        this->RotateAndNormalizeVector(nodal_velocity, rGeometry[itNode], VELOCITY_THRESHOLD);

                        // apply friction force in opposite direction of tangential velocity components
                        for (unsigned dim = 1; dim < this->GetDomainSize(); dim++) {
                            rLocalVector[j + dim] -= FRICTION_FORCE * nodal_velocity[dim];
                        }

                        rGeometry[itNode].Set(INLET);
                    }
                    rGeometry[itNode].UnSetLock();
				}
			}
		}
	}

	/// RHS only version of ApplySlipCondition
	void ApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const override
	{
        // creates an empty dummy matrix to pass into the 'full' ApplySlipCondition -- this dummy matrix is
        // ignored, effectively only updating the RHS
        TLocalMatrixType dummyMatrix;
        this->ApplySlipCondition(dummyMatrix, rLocalVector, rGeometry);
	}

	// An extra function to distinguish the application of slip in element considering penalty imposition
	void ElementApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		// If it is not a penalty element, do as standard
		// Otherwise, if it is a penalty element, don't do anything
		if (!this->IsPenalty(rGeometry))
		{
			this->ApplySlipCondition(rLocalMatrix, rLocalVector, rGeometry);
		}
	}

	// An extra function to distinguish the application of slip in element considering penalty imposition (RHS Version)
	void ElementApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		// If it is not a penalty element, do as standard
		// Otherwise, if it is a penalty element, don't do anything
		if (!this->IsPenalty(rGeometry))
		{
			this->ApplySlipCondition(rLocalVector, rGeometry);
		}
	}

	// An extra function to distinguish the application of slip in condition considering penalty imposition
	void ConditionApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry,
            const ProcessInfo& rCurrentProcessInfo) const
	{
		// If it is not a penalty condition, do as standard
		if (!this->IsPenalty(rGeometry))
		{
			this->ApplySlipCondition(rLocalMatrix, rLocalVector, rGeometry);
		}
		// Otherwise, do the following modification [ONLY applied to penalty conditions]
		else
		{
			const unsigned int LocalSize = rLocalVector.size();

			if (LocalSize > 0)
			{
				const unsigned int block_size = this->GetBlockSize();
                const unsigned int domain_size = this->GetDomainSize();
				TLocalMatrixType temp_matrix = ZeroMatrix(rLocalMatrix.size1(),rLocalMatrix.size2());
				for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
				{
                    // if FLAG_VARIABLE set to 1 in rCurrentProcessInfo, treat SLIP nodes as fixed
					if(this->IsSlip(rGeometry[itNode]) && rCurrentProcessInfo.GetValue(FLAG_VARIABLE) <= 0)
					{
						// First displacement dof (normal component) for each rotated block is always constrained
						unsigned int j = itNode * block_size;

                        double mu = rGeometry[itNode].GetValue(FRICTION_COEFFICIENT);

                        if (mu > 0) {
                            // Positive friction coefficient -- friction active

                            // ONLY need to modify the LHS/RHS for friction if node is in sliding state, i.e. the
                            // max tangential force has been exceeded for the node => friction_state == true
                            int friction_state = rGeometry[itNode].FastGetSolutionStepValue(FRICTION_STATE, 0);

                            if (friction_state == SLIDING) {
                                // check if dynamic friction has already been set
                                rGeometry[itNode].SetLock();
                                bool dyn_friction_set = rGeometry[itNode].Is(OUTLET);
                                rGeometry[itNode].Set(OUTLET);
                                rGeometry[itNode].UnSetLock();
//                                KRATOS_WATCH(dyn_friction_set);

                                // if dynamic friction has not been set, set it -- otherwise zero out current contribution
//                                for (unsigned int i = 1; i < domain_size; i++)
                                    rLocalVector[j + 1] = dyn_friction_set ? 0.0 : rGeometry[itNode].GetValue(MAX_TANGENT_FORCE);
//                                    KRATOS_WATCH(rLocalVector[j+1]);
//                                    KRATOS_WATCH(rGeometry[itNode]);
                                // allow slip along tangential direction -- zero out tangential penalty terms on LHS
                                for (unsigned int k = j + 1; k < j + block_size; ++k) {
                                    for (unsigned int i = 0; i < rLocalMatrix.size1(); ++i) {
                                        rLocalMatrix(i, k) = 0.0;
                                        rLocalMatrix(k, i) = 0.0;
                                    }
                                }
                            } // Else, do nothing -- max tangent force not exceeded, node is fixed
                        }
                        else {
                            // Friction off -- remove all other value in RHS & LHS than the normal component
                            for (unsigned int i = j; i < (j + block_size); ++i) {
                                if (i != j)
                                    rLocalVector[i] = 0.0;

                                for (unsigned int k = 0; k < rLocalMatrix.size1(); ++k) {
                                    rLocalMatrix(i, k) = 0.0;
                                    rLocalMatrix(k, i) = 0.0;
                                }
                            }
                        }
                    }
                }
			}
		}
	}

	// An extra function to distinguish the application of slip in condition considering penalty imposition (RHS Version)
	void ConditionApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo) const
	{
        // creates an empty dummy matrix to pass into the 'full' ConditionApplySlipCondition -- this dummy matrix is
        // ignored, effectively only updating the RHS
        TLocalMatrixType dummyMatrix;
        this->ConditionApplySlipCondition(dummyMatrix, rLocalVector, rGeometry, rCurrentProcessInfo);
	}

	// Checking whether it is normal element or penalty element
	bool IsPenalty(GeometryType& rGeometry) const
	{
		bool is_penalty = false;
		for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
		{
			if(this->IsSlip(rGeometry[itNode]) )
			{
				const double identifier = rGeometry[itNode].FastGetSolutionStepValue(mrFlagVariable);
				const double tolerance  = 1.e-6;
				if (identifier > 1.00 + tolerance)
				{
					is_penalty = true;
					break;
				}
			}
		}

		return is_penalty;
	}

	/// Same functionalities as RotateVelocities, just to have a clear function naming
	virtual	void RotateDisplacements(ModelPart& rModelPart) const
	{
		this->RotateVelocities(rModelPart);
	}

	/// Transform nodal displacement to the rotated coordinates (aligned with each node's normal)
	/// The name is kept to be Rotate Velocities, since it is currently a derived class of coordinate_transformation_utilities in the core
	void RotateVelocities(ModelPart& rModelPart) const override
	{
		TLocalVectorType displacement(this->GetDomainSize());
		TLocalVectorType Tmp(this->GetDomainSize());

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
		#pragma omp parallel for firstprivate(displacement,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			if( this->IsSlip(*itNode) )
			{
				//this->RotationOperator<TLocalMatrixType>(Rotation,);
				if(this->GetDomainSize() == 3)
				{
					BoundedMatrix<double,3,3> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
					for(unsigned int i = 0; i < 3; i++) displacement[i] = rDisplacement[i];
					noalias(Tmp) = prod(rRot,displacement);
					for(unsigned int i = 0; i < 3; i++) rDisplacement[i] = Tmp[i];
				}
				else
				{
					BoundedMatrix<double,2,2> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
					for(unsigned int i = 0; i < 2; i++) displacement[i] = rDisplacement[i];
					noalias(Tmp) = prod(rRot,displacement);
					for(unsigned int i = 0; i < 2; i++) rDisplacement[i] = Tmp[i];
				}
			}
		}
	}

	/// Same functionalities as RecoverVelocities, just to have a clear function naming
	virtual void RecoverDisplacements(ModelPart& rModelPart) const
	{
		this->RecoverVelocities(rModelPart);
	}

	/// Transform nodal displacement from the rotated system to the original configuration
	/// The name is kept to be Recover Velocities, since it is currently a derived class of coordinate_transformation_utilities in the core
	void RecoverVelocities(ModelPart& rModelPart) const override
	{
		TLocalVectorType displacement(this->GetDomainSize());
		TLocalVectorType Tmp(this->GetDomainSize());

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
		#pragma omp parallel for firstprivate(displacement,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			if( this->IsSlip(*itNode) )
			{
				if(this->GetDomainSize() == 3)
				{
					BoundedMatrix<double,3,3> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
					for(unsigned int i = 0; i < 3; i++) displacement[i] = rDisplacement[i];
					noalias(Tmp) = prod(trans(rRot),displacement);
					for(unsigned int i = 0; i < 3; i++) rDisplacement[i] = Tmp[i];
				}
				else
				{
					BoundedMatrix<double,2,2> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
					for(unsigned int i = 0; i < 2; i++) displacement[i] = rDisplacement[i];
					noalias(Tmp) = prod(trans(rRot),displacement);
					for(unsigned int i = 0; i < 2; i++) rDisplacement[i] = Tmp[i];
				}
			}
		}
	}

    // Constraints velocities and accelerations normal to SLIP condition to zero
    void ConstraintDerivatives(array_1d<double, 3>& rCurrentVelocity,
                               array_1d<double, 3>& rCurrentAcceleration,
                               const NodeType& rNode){
        // rotate to local normal-tangential frame
        this->RotateVector(rCurrentVelocity, rNode);
        this->RotateVector(rCurrentAcceleration, rNode);

        rCurrentVelocity[0] = 0;
        rCurrentAcceleration[0] = 0;

        // rotate back to global coordinates
        this->RotateVector(rCurrentVelocity, rNode, true);
        this->RotateVector(rCurrentAcceleration, rNode, true);
    }

    // Sets FRICTION_STATE for each friction node to indicate its stick/sliding state.
    // Also stores the maximum tangent force allowed in each direction.
    void AssignFrictionState(ModelPart& rModelPart){
        ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();

        bool isFirstTimestep = (rModelPart.GetProcessInfo()[FLAG_VARIABLE] == 1.0);

        for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
        {
            ModelPart::NodeIterator itNode = it_begin+iii;

            const double mu = itNode->GetValue(FRICTION_COEFFICIENT);
            const double nodal_mass = itNode->FastGetSolutionStepValue(NODAL_MASS, 0);

            // for all SLIP nodes with active friction that contain material
            if( this->IsSlip(*itNode) &&  mu > 0 && nodal_mass > std::numeric_limits<double>::epsilon())
            {
                // Rotate REACTION to normal-tangential frame of reference
                // [by-value to avoid modifying REACTION value stored at nodes]
                array_1d<double,3> reaction = itNode->FastGetSolutionStepValue(REACTION);
                this->RotateVector(reaction, *itNode);

                int& r_friction_state = itNode->FastGetSolutionStepValue(FRICTION_STATE, 0);

                // update normal force norm assoc. with current timestep
                // [note: no friction if REACTION[0] > 0 [i.e. contact lost]]
                itNode->FastGetSolutionStepValue(NORMAL_REACTION, 0) = fmax(-reaction[0], 0.0);

                // obtain normal (prev timestep) and tangent forces (current) assoc. with node
                const double prev_normal_force_norm = itNode->FastGetSolutionStepValue(NORMAL_REACTION, 1);
                const double tangent_force1 = reaction[1];
                const double tangent_force2 = reaction[2];

                double tangent_force_norm = sqrt(tangent_force1 * tangent_force1 + tangent_force2 * tangent_force2);

                double max_tangent_force_norm = prev_normal_force_norm * mu;

//                KRATOS_WATCH(*itNode);
//                KRATOS_WATCH(prev_normal_force_norm);
//                KRATOS_WATCH(tangent_force_norm);
//                KRATOS_WATCH(max_tangent_force_norm);

                r_friction_state = (tangent_force_norm >= max_tangent_force_norm) ? SLIDING : STICK;

                if (r_friction_state == SLIDING) {
                    // TODO: extend for 3d and add check for tangent_force_norm being close to zero
                    const double tangent_force_dir1 = tangent_force1 / tangent_force_norm;

                    itNode->GetValue(MAX_TANGENT_FORCE) = tangent_force_dir1 * max_tangent_force_norm;
                }

//                KRATOS_WATCH(r_friction_state);
            } else {
                // If friction no longer active OR node contains no material, set FRICTION_STATE to SLIDING
                itNode->FastGetSolutionStepValue(FRICTION_STATE) = SLIDING;
            }
        }
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

	/// Turn back information as a string.
	std::string Info() const override
	{
		std::stringstream buffer;
		buffer << "MPMBoundaryRotationUtility";
		return buffer.str();
	}

	/// Print information about this object.
	void PrintInfo(std::ostream& rOStream) const override
	{
		rOStream << "MPMBoundaryRotationUtility";
	}

	/// Print object's data.
	void PrintData(std::ostream& rOStream) const override {}

	///@}
	///@name Friends
	///@{

	///@}

protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:
	///@name Static Member Variables
	///@{

	const Variable<double>& mrFlagVariable;

    const double FRICTION_FORCE = 0000;
    const double VELOCITY_THRESHOLD = 1e-10;

    const int SLIDING = 0;
    const int STICK = 1;

	///@}
	///@name Member Variables
	///@{

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{


    /// Helper function to rotate a 3-vector to and from the coordinate system defined by the NORMAL defined at rNode
    /**
     @param rVector Vector to be rotated
     @param rNode A reference to the node associated with the vector
     @param toGlobalCoordinates If true, instead rotates the vector back to the global coordinates [default: false]
     */
    void RotateVector(array_1d<double, 3>& rVector,
                      const Node& rNode,
                      const bool toGlobalCoordinates = false) const {

        array_1d<double, 3> rotated_nodal_vector = ZeroVector(3);
        BoundedMatrix<double, 3, 3> rotation_matrix = ZeroMatrix(3);

        this->LocalRotationOperatorPure(rotation_matrix, rNode);
        noalias(rotated_nodal_vector) = prod(toGlobalCoordinates ? trans(rotation_matrix) : rotation_matrix, rVector);

        rVector = rotated_nodal_vector;
    }

    /// Additionally normalizes the rotated vector
    /**
     @param rVector Vector to be rotated
     @param rNode A reference to the node associated with the vector
     @param toGlobalCoordinates If true, instead rotates the vector back to the global coordinates [default: false]
     @param threshold Value below which the value of the component is considered 0 [ default value: machine epsilon ]
     */
    void RotateAndNormalizeVector( array_1d<double, 3> &rVector,
                      const Node &rNode,
                      const bool toGlobalCoordinates = false,
                      const double threshold = std::numeric_limits<double>::epsilon() ) const {
        RotateVector(rVector, rNode, toGlobalCoordinates);

        // Check if velocity is close to zero [ALL components below threshold]
        bool is_zero_vector = true;

        // Evaluates to false if ANY of the rVector component is above threshold
        for (auto component : rVector) {
            is_zero_vector = is_zero_vector && (abs(component) < threshold);
        }

        // If not close to zero, rotate the vector and obtain its norm
        // Otherwise do nothing
        if(!is_zero_vector){
            this->Normalize(rVector);
        }
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

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TLocalMatrixType, class TLocalVectorType>
inline std::istream& operator >>(std::istream& rIStream,
		MPMBoundaryRotationUtility<TLocalMatrixType, TLocalVectorType>& rThis) {
	return rIStream;
}

/// output stream function
template<class TLocalMatrixType, class TLocalVectorType>
inline std::ostream& operator <<(std::ostream& rOStream,
		const MPMBoundaryRotationUtility<TLocalMatrixType, TLocalVectorType>& rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

///@}

///@} addtogroup block

}

#endif // KRATOS_MPM_BOUNDARY_ROTATION_UTILITY
