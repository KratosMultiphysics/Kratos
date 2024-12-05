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
#include "utilities/parallel_utilities.h"
#include "mpm_application_variables.h"

namespace Kratos {

///@addtogroup MPMApplication
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
    using CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double>::RevertRotate;

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
		const unsigned int BlockSize):
    CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double>(DomainSize,BlockSize,SLIP)
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

    void RevertRotate(
        TLocalMatrixType& rLocalMatrix,
        TLocalVectorType& rLocalVector,
        GeometryType& rGeometry) const
    {
        if (this->GetBlockSize() == this->GetDomainSize()) // irreducible case
        {
            if (this->GetDomainSize() == 2) this->template RotateAuxPure<2,true>(rLocalMatrix,rLocalVector,rGeometry);
            else if (this->GetDomainSize() == 3) this->template RotateAuxPure<3,true>(rLocalMatrix,rLocalVector,rGeometry);
        }
        else // mixed formulation case
        {
            if (this->GetDomainSize() == 2) this->template RotateAux<2,3,true>(rLocalMatrix,rLocalVector,rGeometry);
            else if (this->GetDomainSize() == 3) this->template RotateAux<3,4,true>(rLocalMatrix,rLocalVector,rGeometry);
        }

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
                // Checks for conforming slip
                // [ non-conforming SLIP BCs are treated within the resp. condition itself ]
				if( this->IsConformingSlip(rGeometry[itNode]) )
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

                    // Modifications to introduce friction force
                    const double mu = rGeometry[itNode].GetValue(FRICTION_COEFFICIENT);
                    const double tangential_penalty_factor = rGeometry[itNode].GetValue(TANGENTIAL_PENALTY_FACTOR);

                    if (mu > 0 && tangential_penalty_factor > 0) { // Friction active
                        array_1d<double,3> & r_stick_force = rGeometry[itNode].FastGetSolutionStepValue(STICK_FORCE);

                        // accumulate normal forces (RHS vector) for subsequent re-determination of friction state
                        // and check if friction force has been computed for the current node
                        rGeometry[itNode].SetLock();
                        rGeometry[itNode].FastGetSolutionStepValue(STICK_FORCE)[0] -= rLocalVector[j];
                        const bool friction_assigned = rGeometry[itNode].GetValue(FRICTION_ASSIGNED);
                        rGeometry[itNode].SetValue(FRICTION_ASSIGNED, true);
                        rGeometry[itNode].UnSetLock();

                        if (!friction_assigned) {
                            // obtain displacement in normal-tangential coordinates
                            // [ use a copy and not reference to avoid rotating the displacement in SolutionStepValue ]
                            array_1d<double,3> displacement_copy = rGeometry[itNode].FastGetSolutionStepValue(DISPLACEMENT);
                            RotateVector(displacement_copy, rGeometry[itNode]);

                            // determine penalty-based sticking force in the tangential direction for the current node
                            // [ displacement in MPM is the displacement update -> penalize any and all displacement
                            //   updates in the tangential direction [this assumes a stationary background grid]
                            for( unsigned int i = 1; i < this->GetDomainSize(); ++i) {
                                r_stick_force[i] = -tangential_penalty_factor * displacement_copy[i];
                            }

                            const int friction_state = rGeometry[itNode].FastGetSolutionStepValue(FRICTION_STATE);

                            if(friction_state == STICK) {
                                if (rLocalMatrix.size1() != 0) {
                                    for (unsigned int i = 1; i < this->GetDomainSize(); ++i)
                                        rLocalMatrix(j + i, j + i) += tangential_penalty_factor;
                                }
                            } // else, sliding: nothing needed for LHS

                            // add friction to RHS
                            const array_1d<double,3> & r_friction_force = rGeometry[itNode].FastGetSolutionStepValue(REACTION);

                            for( unsigned int i = 1; i < this->GetDomainSize(); ++i)
                                rLocalVector[j + i] += r_friction_force[i];
                        }
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

    bool IsParticleBasedSlip(const NodeType& rNode) const
    {
        return rNode.GetValue(PARTICLE_BASED_SLIP);
    }

    // Checking whether it is normal element or penalty element
    bool IsParticleBasedSlip(const GeometryType& rGeometry) const
    {
        for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
        {
            if(IsParticleBasedSlip(rGeometry[itNode]))
                return true;
        }

        return false;
    }

    bool IsConformingSlip(const NodeType& rNode) const {
        return rNode.Is(SLIP) && !IsParticleBasedSlip(rNode);
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

            // rotate conforming slip ONLY -- particle-based slip is handled by the relevant condition
            if( this->IsConformingSlip(*itNode) )
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
            // rotate conforming slip ONLY -- particle-based slip is handled by the relevant condition
            if( this->IsConformingSlip(*itNode) )
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


    /// Helper function to rotate a 3-vector to and from the coordinate system defined by the NORMAL defined at rNode
    /**
     @param rVector Vector to be rotated
     @param rNode A reference to the node associated with the vector
     @param toGlobalCoordinates If true, instead rotates the vector back to the global coordinates [default: false]
     */
    inline void RotateVector(array_1d<double, 3>& rVector,
                      const Node& rNode,
                      const bool toGlobalCoordinates = false) const {

        TLocalVectorType rotated_nodal_vector(this->GetDomainSize());
        TLocalVectorType tmp(this->GetDomainSize());

        if(this->GetDomainSize() == 3){
            // 3D case
            BoundedMatrix<double, 3, 3> rotation_matrix;
            this->LocalRotationOperatorPure(rotation_matrix, rNode);

            noalias(rotated_nodal_vector) = prod(toGlobalCoordinates ? trans(rotation_matrix) : rotation_matrix, rVector);

            rVector = rotated_nodal_vector;
        } else {
            // 2D case
            BoundedMatrix<double, 2, 2> rotation_matrix;
            this->LocalRotationOperatorPure(rotation_matrix, rNode);

            tmp(0) = rVector(0);
            tmp(1) = rVector(1);

            noalias(rotated_nodal_vector) = prod(toGlobalCoordinates ? trans(rotation_matrix) : rotation_matrix, tmp);

            rVector(0) = rotated_nodal_vector(0);
            rVector(1) = rotated_nodal_vector(1);
        }


    }


    // Sets FRICTION_STATE for a SLIP node to indicate its stick/sliding state and stores the friction force in REACTION
    void ComputeFrictionAndResetFlags(ModelPart& rModelPart) const {
        const bool is_initial_loop = !rModelPart.GetProcessInfo()[INITIAL_LOOP_COMPLETE];

        block_for_each(rModelPart.Nodes(), [&](NodeType& rNode)
        {
            const Node& rConstNode = rNode; // const Node reference to avoid issues with previously unset GetValue()

            int& r_friction_state = rNode.FastGetSolutionStepValue(FRICTION_STATE, 0);
            const double mu = rConstNode.GetValue(FRICTION_COEFFICIENT);

            array_1d<double, 3>& r_friction_force = rNode.FastGetSolutionStepValue(REACTION);

            // Limit tangential forces for friction
            if (mu > 0) {
                // obtain normal and tangent forces assoc. with node at the desired timestep
                const double normal_force = rNode.FastGetSolutionStepValue(STICK_FORCE_X, 1);
                const double tangent_force1 = rNode.FastGetSolutionStepValue(STICK_FORCE_Y, 0);
                const double tangent_force2 = rNode.FastGetSolutionStepValue(STICK_FORCE_Z, 0);

                // forces unmodified in normal direction
                r_friction_force[0] = normal_force;

                // [note: no friction if normal component < 0 [direction chosen to be consistent with contact algo]]
                // [since normal point away from the boundary/interface, a normal component < 0 indicates a force
                //  pulling towards the normal (i.e. tensile force)]
                const double normal_force_norm = fmax(normal_force, 0.0);

                const double tangent_force_norm = sqrt(tangent_force1 * tangent_force1 + tangent_force2 * tangent_force2);
                const double max_tangent_force_norm = normal_force_norm * mu;

                // special treatment for initial loop
                if (is_initial_loop) {
                    if (rConstNode.GetValue(HAS_INITIAL_MOMENTUM)) {
                        r_friction_state = SLIDING;
                    }
                    else {
                        r_friction_state = STICK;
                    }
                } else {
                    r_friction_state = (tangent_force_norm >= max_tangent_force_norm) ? SLIDING : STICK;
                }

                if (r_friction_state == SLIDING) {
                    double tangent_force_dir1 = 0.0;
                    double tangent_force_dir2 = 0.0;

                    if (tangent_force_norm > std::numeric_limits<double>::epsilon()) {
                        tangent_force_dir1 = tangent_force1 / tangent_force_norm;
                        tangent_force_dir2 = tangent_force2 / tangent_force_norm;
                    }

                    r_friction_force[1] = tangent_force_dir1 * max_tangent_force_norm;
                    r_friction_force[2] = tangent_force_dir2 * max_tangent_force_norm;
                }
                else { // STICK
                    r_friction_force[1] = tangent_force1;
                    r_friction_force[2] = tangent_force2;
                }

            }
        });
    }

	///@}
	///@name Access
	///@{
	int GetSlidingState() const{
        return SLIDING;
    }

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
    constexpr static int SLIDING = 0;
    constexpr static int STICK = 1;

	///@}
	///@name Member Variables
	///@{

    ///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

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
