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
    CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double>(DomainSize,BlockSize,SLIP), mrFlagVariable(rVariable), mDomainSize(DomainSize)
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

                    // Set value of normal displacement at node directly to the normal displacement of the boundary mesh
					rLocalVector[j] = inner_prod(rN,displacement);
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
		// Otherwise, do the following modification
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
					if(this->IsSlip(rGeometry[itNode]))
					{
                        // First displacement dof of a Node is the normal component
						unsigned int j = itNode * block_size;

                        double mu = rGeometry[itNode].GetValue(FRICTION_COEFFICIENT);

                        if (mu > 0) {
                            // Friction active

                            int friction_state = rGeometry[itNode].FastGetSolutionStepValue(FRICTION_STATE, 0);

                            if (friction_state == SLIDING) {
                                // check if dynamic friction has already been set
                                rGeometry[itNode].SetLock();
                                bool dyn_friction_set = rGeometry[itNode].Is(OUTLET);
                                rGeometry[itNode].Set(OUTLET);
                                rGeometry[itNode].UnSetLock();

                                // Zero out tangential penalty terms of RHS corresponding to SLIDING node on RHS, adding
                                // the dynamic friction computed in AssignFrictionState if not yet added
                                for (unsigned int i = 1; i < domain_size; i++)
                                    rLocalVector[j + i] = dyn_friction_set ? 0.0 : rGeometry[itNode].GetValue(FRICTION_CONTACT_FORCE)[i];

                                // Zero out tangential penalty terms of LHS corresponding to SLIDING node
                                for (unsigned int k = j + 1; k < j + block_size; ++k) {
                                    for (unsigned int i = 0; i < rLocalMatrix.size1(); ++i) {
                                        rLocalMatrix(i, k) = 0.0;
                                        rLocalMatrix(k, i) = 0.0;
                                    }
                                }
                            } // Else, do nothing -- max tangent force not exceeded, node is fixed
                        }
                        else {
                            // Friction inactive -- zero out tangential entries in RHS & LHS
                            for (unsigned int i = j + 1; i < (j + block_size); ++i) {
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

    // Auxiliary function to clear friction-related flags --
    // MUST be called before (re-)building the RHS in a given non-linear iteration
    static void ClearFrictionFlag(const ModelPart& rModelPart) {
        KRATOS_TRY
        // Loop over the grid nodes performed to clear OUTLET flag for indicating that nodal friction has alr been set
        // and remove the accumulated normal forces
        for(NodeType &curr_node : rModelPart.Nodes()){
            curr_node.Reset(OUTLET);
            curr_node.GetValue(FRICTION_CONTACT_FORCE).clear();
            curr_node.FastGetSolutionStepValue(FRICTION_CONTACT_FORCE, 0).clear();
        }
        KRATOS_CATCH( "" )
    }

    // Sets FRICTION_STATE for each friction node to indicate its stick/sliding state.
    // Also stores the maximum tangent force allowed in each direction.
    void AssignFrictionState(ModelPart& rModelPart){
        ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();

        bool isInitialLoop = (rModelPart.GetProcessInfo()[FLAG_VARIABLE] > 0);

        for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
        {
            ModelPart::NodeIterator itNode = it_begin+iii;

            const double mu = itNode->GetValue(FRICTION_COEFFICIENT);
            const double nodal_mass = itNode->FastGetSolutionStepValue(NODAL_MASS, 0);

            // for all SLIP nodes
            if( this->IsSlip(*itNode) )
            {
                // Rotate REACTION to normal-tangential frame of reference
                array_1d<double,3>& r_reaction = itNode->FastGetSolutionStepValue(REACTION);
                this->RotateVector(r_reaction, *itNode);

                int& r_friction_state = itNode->FastGetSolutionStepValue(FRICTION_STATE, 0);

                // Limit/zero out tangential forces for friction/slip case resp.
                if ( mu > 0  && nodal_mass > std::numeric_limits<double>::epsilon() ) {
                    // update nodal normal & tangent forces assoc. with current timestep
                    // [note: no friction if normal component < 0 [direction chosen to be consistent with contact algo]]
                    // [since normal point away from the boundary/interface, a normal component < 0 indicates a force
                    //  pulling towards the normal (i.e. tensile force)]
                    itNode->FastGetSolutionStepValue(FRICTION_CONTACT_FORCE_X, 0) = fmax(r_reaction[0], 0.0);
                    itNode->FastGetSolutionStepValue(FRICTION_CONTACT_FORCE_Y, 0) = r_reaction[1];
                    itNode->FastGetSolutionStepValue(FRICTION_CONTACT_FORCE_Z, 0) = r_reaction[2];

                    // obtain normal and tangent forces assoc. with node at the desired timestep
                    // [ currently: normal forces from prev timestep, tangent forces current timestep ]
                    const double normal_force_norm = itNode->FastGetSolutionStepValue(FRICTION_CONTACT_FORCE_X, 1);
                    const double tangent_force1 = itNode->FastGetSolutionStepValue(FRICTION_CONTACT_FORCE_Y, 0);
                    const double tangent_force2 = itNode->FastGetSolutionStepValue(FRICTION_CONTACT_FORCE_Z, 0);

                    const double tangent_force_norm = sqrt(tangent_force1 * tangent_force1 + tangent_force2 * tangent_force2);
                    const double max_tangent_force_norm = normal_force_norm * mu;

                    // special treatment for initial loop
                    if (isInitialLoop) {
                        if (itNode->Is(INLET)) { // used to mark nodes with non-zero momentum in the initial timestep
                            r_friction_state = SLIDING;
                        }
                        else {
                            r_friction_state = STICK;
                        }
                    }
                    else {
                        r_friction_state = (tangent_force_norm >= max_tangent_force_norm) ? SLIDING : STICK;
                    }

                    if (r_friction_state == SLIDING) {
                        double tangent_force_dir1 = 0.0;
                        double tangent_force_dir2 = 0.0;

                        if (tangent_force_norm > std::numeric_limits<double>::epsilon()) {
                            tangent_force_dir1 = tangent_force1 / tangent_force_norm;
                            tangent_force_dir2 = tangent_force2 / tangent_force_norm;
                        }

                        // SLIDING -> sets max tangent forces for use in ConditionApplyCondition and modify r_reaction
                        // [ note the use of SetValue instead of FastGetSolutionStepValue ]
                        itNode->SetValue(FRICTION_CONTACT_FORCE_Y, tangent_force_dir1 * max_tangent_force_norm);
                        itNode->SetValue(FRICTION_CONTACT_FORCE_Z, tangent_force_dir2 * max_tangent_force_norm);

                        r_reaction[1] = itNode->GetValue(FRICTION_CONTACT_FORCE_Y);
                        r_reaction[2] = itNode->GetValue(FRICTION_CONTACT_FORCE_Z);
                    }
                }
                else {
                    // If friction not active OR node contains no material, set FRICTION_STATE to SLIDING
                    // and treat as pure slip case
                    r_friction_state = SLIDING;
                    itNode->SetValue(FRICTION_CONTACT_FORCE_Y, 0.0);
                    itNode->SetValue(FRICTION_CONTACT_FORCE_Z, 0.0);

                    r_reaction[1] = itNode->GetValue(FRICTION_CONTACT_FORCE_Y);
                    r_reaction[2] = itNode->GetValue(FRICTION_CONTACT_FORCE_Z);
                }

                // Rotate REACTION back to global frame of reference
                this->RotateVector(r_reaction, *itNode, true);
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

    const int SLIDING = 0;
    const int STICK = 1;

    /// Number of spatial dimensions
    const unsigned int mDomainSize;

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

        if(mDomainSize == 3){
            // 3D case
            array_1d<double, 3> rotated_nodal_vector = ZeroVector(3);
            BoundedMatrix<double, 3, 3> rotation_matrix = ZeroMatrix(3);
            this->LocalRotationOperatorPure(rotation_matrix, rNode);

            noalias(rotated_nodal_vector) = prod(toGlobalCoordinates ? trans(rotation_matrix) : rotation_matrix, rVector);

            rVector = rotated_nodal_vector;
        } else {
            // 2D case
            array_1d<double, 2> rotated_nodal_vector = ZeroVector(2);
            BoundedMatrix<double, 2, 2> rotation_matrix = ZeroMatrix(2);
            this->LocalRotationOperatorPure(rotation_matrix, rNode);

            array_1d<double, 2> tmp = ZeroVector(2);
            tmp(0) = rVector(0);
            tmp(1) = rVector(1);

            noalias(rotated_nodal_vector) = prod(toGlobalCoordinates ? trans(rotation_matrix) : rotation_matrix, tmp);

            rVector(0) = rotated_nodal_vector(0);
            rVector(1) = rotated_nodal_vector(1);
        }


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
