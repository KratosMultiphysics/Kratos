//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#ifndef BEHR_WALL_CONDITION_H
#define BEHR_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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

/// Implements a wall condition for the Navier-Stokes monolithic formulation.
/**
  It is intended to be used in combination with ASGS Navier-Stokes symbolic elements or their
  derived classes and the ResidualBasedIncrementalUpdateStaticSchemeSlip time scheme, which supports
  slip conditions.
  @see NavierStokes,EmbeddedNavierStokes,ResidualBasedIncrementalUpdateStaticSchemeSlip
 */
template< unsigned int TDim, unsigned int TNumNodes = TDim >
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) BehrWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NavierStokesWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(BehrWallCondition);

    struct ConditionDataStruct
    {
        double wGauss;                  // Gauss point weight
        double charVel;                 // Problem characteristic velocity (used in the outlet inflow prevention)
        double delta;                   // Non-dimensional positive sufficiently small constant (used in the outlet inflow prevention)
        array_1d<double, 3> Normal;     // Condition normal
        array_1d<double, TNumNodes> N;  // Gauss point shape functions values
    };

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef Geometry<NodeType>::PointType PointType;

	typedef Geometry<NodeType>::GeometriesArrayType GeometriesArrayType;

	typedef Element::WeakPointer ElementWeakPointerType;

	typedef Element::Pointer ElementPointerType;

	typedef std::size_t SizeType;
    /// type for shape function values
	typedef Vector ShapeFunctionsType;
	/// Type for a matrix containing the shape function gradients
	typedef Matrix ShapeFunctionDerivativesType;
	/// Type for an array of shape function gradient matrices
	typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    typedef array_1d<double, 3> Normal;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new         // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);condition
      */
    BehrWallCondition(IndexType NewId = 0):Condition(NewId)
    {
        std::cout << "BEHR" << std::endl;
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    BehrWallCondition(IndexType NewId, const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
        std::cout << "BEHR" << std::endl;
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    BehrWallCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
        std::cout << "BEHR" << std::endl;
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    BehrWallCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
        std::cout << "BEHR" << std::endl;
    }

    /// Copy constructor.
    BehrWallCondition(BehrWallCondition const& rOther):
        Condition(rOther)
    {
        std::cout << "BEHR" << std::endl;
    }

    /// Destructor.
    ~BehrWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    BehrWallCondition & operator=(BehrWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new NavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<BehrWallCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new NavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< BehrWallCondition >(NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    /// Calculates the LHS and RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_gauss;
        BoundedMatrix<double,MatrixSize, MatrixSize> lhs_gauss;

        // LHS and RHS contributions initialization
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Compute condition unit normal vector
        this->CalculateNormal(data.Normal); //this already contains the area
        const double A = norm_2(data.Normal);
        data.Normal /= A;

        // Store the outlet inflow prevention constants in the data structure
        data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
        const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // const to avoid race conditions on data_value_container access/initialization
        data.charVel = rProcessInfo[CHARACTERISTIC_VELOCITY];

        // Gauss point information
        GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        Vector GaussPtsJDet = ZeroVector(NumGauss);
        rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        // Loop on gauss points
        for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = GaussPtsJDet[igauss];
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            ComputeGaussPointRHSContribution(rhs_gauss, data);
            ComputeGaussPointLHSContribution(lhs_gauss, data);

            noalias(rLeftHandSideMatrix) += lhs_gauss;
            noalias(rRightHandSideVector) += rhs_gauss;
        }

        KRATOS_CATCH("")
    }

    /// Calculates the RHS condition contributions                                           ----- TO DO -----------
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        std::cout << "BEHR CONDITION was initilized" << std::endl;

        GeometryType& rGeom = this->GetGeometry();     // one edge or one face

        std::vector<IndexType> NodeIds(TNumNodes), ElementNodeIds;
        std::vector< Matrix > NodalNormals(TNumNodes);

        // Retrieving the nodal normal vectors and storing them in matrix format
		for (SizeType i=0; i < TNumNodes; i++){
			NodeIds[i] = rGeom[i].Id();
            NodalNormals[i].resize(3,1);
            for (int j = 0; j < 3; j++)
            {
                NodalNormals[i](j,0) = rGeom[i].GetValue(NORMAL)(j);
            }
		}

        // Computation of NodalMultMatrix = ( [I] - (na)(na) )
        std::vector<MatrixType> NodalMultMatrix(TNumNodes);

        const MatrixType auxIdentMatrix = identity_matrix<double>(3);
        MatrixType auxMatrix = zero_matrix<double>(3,3);
        
        for (SizeType i=0; i < TNumNodes; i++){
            auxMatrix = prod( NodalNormals[i], trans(NodalNormals[i]) );
            NodalMultMatrix[i] = auxIdentMatrix - auxMatrix;
        }

        // Finding the parent element
        // this->FindParentElement();


        const ShapeFunctionsType& N = row(this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_2),0);

        // this->GetGeometry().ShapeFunctionsLocalGradients

        // const ShapeFunctionDerivativesType& DN;
        
        KRATOS_CATCH("")
    }


    void Initialize() override
	{
		KRATOS_TRY;

        std::cout << "BEHR CONDITION was initilized" << std::endl;

		const array_1d<double,3>& rNormal = this->GetValue(NORMAL);
		if (norm_2(rNormal) == 0.0){
		    std::cout << "error on condition -> " << this->Id() << std::endl;
		    KRATOS_THROW_ERROR(std::logic_error, "NORMAL must be calculated before using this condition","");
		}

		if (mInitializeWasPerformed){
		        return;
		}

		mInitializeWasPerformed = true;

		double EdgeLength;
		array_1d<double,3> Edge;
		GeometryType& rGeom = this->GetGeometry();
		WeakPointerVector<Element> ElementCandidates;
		for (SizeType i = 0; i < TNumNodes; i++){

			WeakPointerVector<Element>& rNodeElementCandidates = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);

			for (SizeType j = 0; j < rNodeElementCandidates.size(); j++){
				ElementCandidates.push_back(rNodeElementCandidates(j));
			}
		}

		std::vector<IndexType> NodeIds(TNumNodes), ElementNodeIds;
		for (SizeType i=0; i < TNumNodes; i++){
			NodeIds[i] = rGeom[i].Id();
		}

		std::sort(NodeIds.begin(), NodeIds.end());

		for (SizeType i=0; i < ElementCandidates.size(); i++){
			GeometryType& rElemGeom = ElementCandidates[i].GetGeometry();
			ElementNodeIds.resize(rElemGeom.PointsNumber());

			for (SizeType j=0; j < rElemGeom.PointsNumber(); j++){
				ElementNodeIds[j] = rElemGeom[j].Id();
			}

			std::sort(ElementNodeIds.begin(), ElementNodeIds.end());

			if ( std::includes(ElementNodeIds.begin(), ElementNodeIds.end(), NodeIds.begin(), NodeIds.end()) )
			{
				mParentElement = ElementCandidates(i);

				Edge = rElemGeom[1].Coordinates() - rElemGeom[0].Coordinates();
				mMinEdgeLength = Edge[0]*Edge[0];
				for (SizeType d=1; d < TDim; d++){
					mMinEdgeLength += Edge[d]*Edge[d];
				}

				for (SizeType j=2; j < rElemGeom.PointsNumber(); j++){

					for(SizeType k=0; k < j; k++){
						Edge = rElemGeom[j].Coordinates() - rElemGeom[k].Coordinates();
						EdgeLength = Edge[0]*Edge[0];

						for (SizeType d = 1; d < TDim; d++){
							EdgeLength += Edge[d]*Edge[d];
						}

						mMinEdgeLength = (EdgeLength < mMinEdgeLength) ? EdgeLength : mMinEdgeLength;
					}
				}
				mMinEdgeLength = sqrt(mMinEdgeLength);
				return;
			}
		}

		std::cout << "error in condition -> " << this->Id() << std::endl;
		KRATOS_THROW_ERROR(std::logic_error, "Condition cannot find parent element","");
		KRATOS_CATCH("");
	}


    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rRightHandSideVector reference to the RHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        noalias(rRightHandSideVector) = ZeroVector(TDim+1);

        KRATOS_CATCH("")
    }


    /// Condition check
    /**
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            // Check that all required variables have been registered
            if(VELOCITY.Key() == 0)
                KRATOS_ERROR << "VELOCITY Key is 0. Check if the application was correctly registered.";
            if(MESH_VELOCITY.Key() == 0)
                KRATOS_ERROR << "MESH_VELOCITY Key is 0. Check if the application was correctly registered.";
            if(ACCELERATION.Key() == 0)
                KRATOS_ERROR << "ACCELERATION Key is 0. Check if the application was correctly registered.";
            if(PRESSURE.Key() == 0)
                KRATOS_ERROR << "PRESSURE Key is 0. Check if the application was correctly registered.";
            if(DENSITY.Key() == 0)
                KRATOS_ERROR << "DENSITY Key is 0. Check if the application was correctly registered.";
            if(DYNAMIC_VISCOSITY.Key() == 0)
                KRATOS_ERROR << "DYNAMIC_VISCOSITY Key is 0. Check if the application was correctly registered.";
            if(EXTERNAL_PRESSURE.Key() == 0)
                KRATOS_ERROR << "EXTERNAL_PRESSURE Key is 0. Check if the application was correctly registered.";

            // Checks on nodes
            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
            {
                if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                    KRATOS_ERROR << "missing VELOCITY variable on solution step data for node " << this->GetGeometry()[i].Id();
                if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                    KRATOS_ERROR << "missing PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();
                if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                    KRATOS_ERROR << "missing MESH_VELOCITY variable on solution step data for node " << this->GetGeometry()[i].Id();
                if(this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
                    KRATOS_ERROR << "missing ACCELERATION variable on solution step data for node " << this->GetGeometry()[i].Id();
                if(this->GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_PRESSURE) == false)
                    KRATOS_ERROR << "missing EXTERNAL_PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();
                if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                   this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                   this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                    KRATOS_ERROR << "missing VELOCITY component degree of freedom on node " << this->GetGeometry()[i].Id();
                if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                    KRATOS_ERROR << "missing PRESSURE component degree of freedom on node " << this->GetGeometry()[i].Id();
            }

            return Check;
        }

        KRATOS_CATCH("");
    }

    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo) override;

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
        buffer << "BehrWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BehrWallCondition";
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

    void CalculateNormal(array_1d<double,3>& An);

    void ComputeGaussPointLHSContribution(BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs, const ConditionDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);

    void ComputeRHSNeumannContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);
    void ComputeRHSOutletInflowContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);

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


    ///@}
    ///@name Member Variables
    ///@{

    bool mInitializeWasPerformed;
	double mMinEdgeLength;
	ElementWeakPointerType mParentElement;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

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

}; // Class NavierStokesWallCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, BehrWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const BehrWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_NAVIER_STOKES_WALL_CONDITION_H
