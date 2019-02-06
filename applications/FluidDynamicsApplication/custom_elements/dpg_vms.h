//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Kazem Kamran
//

#if !defined(KRATOS_TWO_FLUID_DPGVMS_H_INCLUDED )
#define  KRATOS_TWO_FLUID_DPGVMS_H_INCLUDED
// System includes
#include <string>
#include <iostream>
// External includes
// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "custom_elements/vms.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "utilities/split_tetrahedra.h"
#include "utilities/enrichment_utilities.h"
// Application includes
#include "fluid_dynamics_application_variables.h"
#include "vms.h"
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
template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
class DPGVMS : public VMS<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition of DPGVMS
    KRATOS_CLASS_POINTER_DEFINITION(DPGVMS);
    ///base type: an IndexedObject that automatically has a unique number
    typedef IndexedObject BaseType;
    ///Element from which it is derived
    typedef VMS<TDim, TNumNodes> ElementBaseType;
    ///definition of node type (default is: Node<3>)
    typedef Node < 3 > NodeType;
    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;
    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;
    ///definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef typename ElementBaseType::MatrixType MatrixType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef std::vector< Dof<double>::Pointer > DofsVectorType;
    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;
    ///@}
    ///@name Life Cycle
    ///@{
    //Constructors.
    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    DPGVMS(IndexType NewId = 0) :
        ElementBaseType(NewId)
    {
    }
    ///Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    DPGVMS(IndexType NewId, const NodesArrayType& ThisNodes) :
        ElementBaseType(NewId, ThisNodes)
    {
    }
    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    DPGVMS(IndexType NewId, GeometryType::Pointer pGeometry) :
        ElementBaseType(NewId, pGeometry)
    {
    }
    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    DPGVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        ElementBaseType(NewId, pGeometry, pProperties)
    {
    }
    /// Destructor.
    ~DPGVMS() override
    {
    }
    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /// Create a new element of this type
    /**
     * Returns a pointer to a new DPGVMS element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<DPGVMS>(NewId, (this->GetGeometry()).Create(ThisNodes), pProperties);
    }

    /// Create a new element of this type.
	/**
	 @param NewId Index of the new element
     @param pGeom A pointer to the geometry of the new element
	 @param pProperties Pointer to the element's properties
	 */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< DPGVMS >(NewId,pGeom,pProperties);
    }

    /// Call at teh begining of each step, ita decides if element is cutted or no!
    /**
      */
    void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo) override
    {
// 	for (unsigned int jj = 0; jj < 4; jj++){
// 	      this->GetGeometry()[jj].FastGetSolutionStepValue(WET_VOLUME ) = 0.0;
// 	      this->GetGeometry()[jj].FastGetSolutionStepValue(CUTTED_AREA ) = 0.0;
// 	}
    }

    /// Call at teh begining of each iteration, ita decides if element is cutted or no!
    /**
      */
    void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo) override
    {
	  // Calculate this element's geometric parameters
	  double Area;
	  array_1d<double, TNumNodes> N;
	  BoundedMatrix<double, TNumNodes, TDim> DN_DX;
	  GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
	  //get position of the cut surface
	  Vector distances(TNumNodes);
	  Matrix Nenriched(6, 1);
	  Vector volumes(6);
	  Matrix coords(TNumNodes, TDim);
	  Matrix Ngauss(6, TNumNodes);
	  Vector signs(6);
	  std::vector< Matrix > gauss_gradients(6);
	  //fill coordinates
	  for (unsigned int i = 0; i < TNumNodes; i++)
	  {
	      const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
	      volumes[i] = 0.0;
	      distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
	      for (unsigned int j = 0; j < TDim; j++)
		  coords(i, j) = xyz[j];
	  }
	  this->GetValue(AUX_INDEX) = 0.0;
	  for (unsigned int i = 0; i < 6; i++)
	      gauss_gradients[i].resize(1, TDim, false);

      array_1d<double,6> edge_areas;
	  unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);

	  if(ndivisions == 1)
	    this->is_cutted = 0;
	  else{
	    this->is_cutted = 1;
		this->GetValue(AUX_INDEX) = 1.0;
	  }
    }


    /// Provides local contributions from body forces and OSS projection terms
    /**
     * This is called during the assembly process and provides the terms of the
     * system that are either constant or computed explicitly (from the 'old'
     * iteration variables). In this case this means the body force terms and the
     * OSS projections, that are treated explicitly.
     * @param rLeftHandSideMatrix the elemental left hand side matrix. Not used here, required for compatibility purposes only.
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
//         this->IsCutted();
        unsigned int LocalSize = (TDim + 1) * TNumNodes;

	if( this->is_cutted == 1)
	{
	    LocalSize += 1;
	    // Check sizes and initialize matrix
	    if (rLeftHandSideMatrix.size1() != LocalSize)
		rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

	    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
            this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
	}
	else
	{
	    // Check sizes and initialize matrix
	    if (rLeftHandSideMatrix.size1() != LocalSize)
		rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

	    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
            ElementBaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

	}

    }


    /// Provides local contributions from body forces and projections to the RHS
    /**
     * This is called during the assembly process and provides the RHS terms of the
     * system that are either constant or computed explicitly (from the 'old'
     * iteration variables). In this case this means the body force terms and the
     * OSS projections, that are treated explicitly.
     * @param rRightHandSideVector Will be filled with the elemental right hand side
     * @param rCurrentProcessInfo ProcessInfo instance from the ModelPart. It is
     * expected to contain values for OSS_SWITCH, DYNAMIC_TAU and DELTA_TIME
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
      	if( this->is_cutted == 1)
	{
	  unsigned int LocalSize = (TDim + 1) * TNumNodes + 1;


	  // Check sizes and initialize
	  if (rRightHandSideVector.size() != LocalSize)
	      rRightHandSideVector.resize(LocalSize, false);
	  noalias(rRightHandSideVector) = ZeroVector(LocalSize);
	  // Calculate this element's geometric parameters
	  double Area;
	  array_1d<double, TNumNodes> N;
	  BoundedMatrix<double, TNumNodes, TDim> DN_DX;
	  GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
	  //get position of the cut surface
	  Vector distances(TNumNodes);
	  Matrix Nenriched(6, 1);
	  Vector volumes(6);
	  Matrix coords(TNumNodes, TDim);
	  Matrix Ngauss(6, TNumNodes);
	  Vector signs(6);
	  std::vector< Matrix > gauss_gradients(6);
	  //fill coordinates
	  for (unsigned int i = 0; i < TNumNodes; i++)
	  {
	      const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
	      volumes[i] = 0.0;
	      distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
	      for (unsigned int j = 0; j < TDim; j++)
		  coords(i, j) = xyz[j];
	  }
	  for (unsigned int i = 0; i < 6; i++)
	      gauss_gradients[i].resize(1, TDim, false);

      array_1d<double,6> edge_areas;
	  unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);
	  //do integration
	  for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
	  {
	      //assigning the gauss data
	      for (unsigned int k = 0; k < TNumNodes; k++)
		  N[k] = Ngauss(igauss, k);
	      double wGauss = volumes[igauss];
	      // Calculate this element's fluid properties
	      double Density;
	      this->EvaluateInPoint(Density, DENSITY, N);
	      // Calculate Momentum RHS contribution
	      this->AddMomentumRHS(rRightHandSideVector, Density, N, wGauss);
	      // For OSS: Add projection of residuals to RHS
  //             if (rCurrentProcessInfo[OSS_SWITCH] == 1)
  //             {
  //                 array_1d<double, 3 > AdvVel;
  //                 this->GetAdvectiveVel(AdvVel, N);
  //                 double KinViscosity;
  //                 this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
  //                 double Viscosity;
  //                 this->GetEffectiveViscosity(Density, KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);
  //                 // Calculate stabilization parameters
  //                 double TauOne, TauTwo;
  // //                    if (ndivisions == 1)
  //                 this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);
  // //                else
  // //                {
  // //                    TauOne = 0.0;
  // //                    TauTwo = 0.0;
  // //                }
  // //                    this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);
  //                 this->AddProjectionToRHS(rRightHandSideVector, AdvVel, Density, TauOne, TauTwo, N, DN_DX, wGauss, rCurrentProcessInfo[DELTA_TIME]);
  //             }
	  }
	}
	else
	 ElementBaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    }
    /// Computes local contributions to the mass matrix
    /**
     * Provides the local contributions to the mass matrix, which is defined here
     * as the matrix associated to velocity derivatives. Note that the mass
     * matrix implemented here is lumped.
     * @param rMassMatrix Will be filled with the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
//       this->IsCutted();
      if( this->is_cutted == 0)
	  ElementBaseType::CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
      else
       {
         const unsigned int LocalSize = (TDim + 1) * TNumNodes + 1;
        // Resize and set to zero
        if (rMassMatrix.size1() != LocalSize)
            rMassMatrix.resize(LocalSize, LocalSize, false);
        rMassMatrix = ZeroMatrix(LocalSize, LocalSize);
        // Get the element's geometric parameters
        double Area;
        array_1d<double, TNumNodes> N;
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
        //get position of the cut surface
        Vector distances(TNumNodes);
        Matrix Nenriched(6, 1);
        Vector volumes(6);
        Matrix coords(TNumNodes, TDim);
        Matrix Ngauss(6, TNumNodes);
        Vector signs(6);
        std::vector< Matrix > gauss_gradients(6);
        //fill coordinates
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
            volumes[i] = 0.0;
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            for (unsigned int j = 0; j < TDim; j++)
                coords(i, j) = xyz[j];
        }
        for (unsigned int i = 0; i < 6; i++)
            gauss_gradients[i] = ZeroMatrix(1,TDim);//.resize(1, TDim, false);

        array_1d<double,6> edge_areas;
        unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);
        //mass matrix
        for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
        {
            //assigning the gauss data
            for (unsigned int k = 0; k < TNumNodes; k++)
                N[k] = Ngauss(igauss, k);
            double wGauss = volumes[igauss];
            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);
            // Consisten Mass Matrix
            this->AddConsistentMassMatrixContribution(rMassMatrix, N, Density, wGauss);
        }
        this->LumpMassMatrix(rMassMatrix);

        //stabilization terms
        for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
        {
            //assigning the gauss data
            for (unsigned int k = 0; k < TNumNodes; k++)
                N[k] = Ngauss(igauss, k);
            double wGauss = volumes[igauss];
            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);
            /* For ASGS: add dynamic stabilization terms.
             * These terms are not used in OSS, as they belong to the finite element
             * space and cancel out with their projections.
             */
            if (rCurrentProcessInfo[OSS_SWITCH] != 1)
            {
                double ElemSize = this->ElementSize(Area);
                double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize, rCurrentProcessInfo);

                // Get Advective velocity
                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                // Calculate stabilization parameters
                double TauOne, TauTwo;
                this->CalculateTau(TauOne, TauTwo, AdvVel, ElemSize, Density, Viscosity, rCurrentProcessInfo);

                // Add dynamic stabilization terms ( all terms involving a delta(u) )
                this->AddMassStabTerms(rMassMatrix, Density, AdvVel, TauOne, N, DN_DX, wGauss,gauss_gradients[igauss]);
            }
        }
       }
    }
    /// Computes the local contribution associated to 'new' velocity and pressure values
    /**
     * Provides local contributions to the system associated to the velocity and
     * pressure terms (convection, diffusion, pressure gradient/velocity divergence
     * and stabilization).
     * @param rDampingMatrix Will be filled with the velocity-proportional "damping" matrix
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) override
    {
//       this->IsCutted();
      if( this->is_cutted == 0){
	    ElementBaseType::CalculateLocalVelocityContribution(rDampingMatrix, rRightHandSideVector, rCurrentProcessInfo);

	    //compute boundary term
	    int boundary_nodes = 0;
	    //unsigned int inside_index = -1;
	    for (unsigned int i = 0; i < TNumNodes; i++)
	    {
	      double nd_flag = this->GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE);
	      if (nd_flag == 5.0)
		boundary_nodes++;
	      //else
		//inside_index = i;
	    }


	/*  if(boundary_nodes == TDim)
	  {
	    // Get this element's geometric properties
	    double Volume;
	    array_1d<double, TNumNodes> N;
	    BoundedMatrix<double, TNumNodes, TDim> DN_DX;
	    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);


	    BoundedMatrix<double, 16, 16 > boundary_damp_matrix;
	    noalias(boundary_damp_matrix) = ZeroMatrix(16,16);
	    array_1d<double,3> face_normal;
	    face_normal[0] = -DN_DX(inside_index,0);
	    face_normal[1] = -DN_DX(inside_index,1);
	    face_normal[2] = -DN_DX(inside_index,2);
	    const double fn = norm_2(face_normal);
	    face_normal/= fn;
	    double face_area = 3.0*Volume*fn;

// 	    int LocalIndex = 0;
// 	    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
// 	    {
// 		const double& pnode = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
// 		for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
// 		{
// 		    if(iNode != inside_index)
// 		      rRightHandSideVector[LocalIndex] -= pnode*face_normal[d]*face_area/3.0  ;
// 		    ++LocalIndex;
// 		}
// 		++LocalIndex;
// 	    }

// 	    int LocalIndex = 0;
// 	    double pface = 0.0;
// 	    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
// 	    {
// 		const double& pnode = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
// 		if(iNode != inside_index) pface += pnode;
// 	    }
// 	    pface/=3.0;
//
// 	    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
// 	    {
// 		for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
// 		{
// 		    if(iNode != inside_index)
// 		      rRightHandSideVector[LocalIndex] -= pface*face_normal[d]*face_area/3.0  ;
// 		    ++LocalIndex;
// 		}
// 		++LocalIndex;
// 	    }
        AddBoundaryTerm(boundary_damp_matrix, DN_DX, N, face_normal, face_area, Volume, rCurrentProcessInfo);

	    VectorType U = ZeroVector(16);
	    int LocalIndex = 0;

	    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
	    {
		array_1d< double, 3 > & rVel = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
		for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
		{
		    U[LocalIndex] = rVel[d];
		    ++LocalIndex;
		}
		U[LocalIndex] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
		++LocalIndex;
	    }
	    noalias(rRightHandSideVector) -= prod(boundary_damp_matrix, U);
	    noalias(rDampingMatrix) += boundary_damp_matrix;
	  }	     */
      }
      else
       {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes + 1;
        // Resize and set to zero the matrix
        // Note that we don't clean the RHS because it will already contain body force (and stabilization) contributions
        if (rDampingMatrix.size1() != LocalSize)
            rDampingMatrix.resize(LocalSize, LocalSize, false);
        noalias(rDampingMatrix) = ZeroMatrix(LocalSize, LocalSize);
        // Get this element's geometric properties
        double Area;
        array_1d<double, TNumNodes> N;

        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
        //get position of the cut surface
        Vector distances(TNumNodes);
        Matrix Nenriched(6, 1);
        Vector volumes(6);
        Matrix coords(TNumNodes, TDim);
        Matrix Ngauss(6, TNumNodes);
        Vector signs(6);
        std::vector< Matrix > gauss_gradients(6);
        //fill coordinates
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
            volumes[i] = 0.0;
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            for (unsigned int j = 0; j < TDim; j++)
                coords(i, j) = xyz[j];
        }
        for (unsigned int i = 0; i < 6; i++)
            gauss_gradients[i] = ZeroMatrix(1,TDim);

        array_1d<double,6> edge_areas;
        unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);
//         Vector enrichment_terms_vertical = ZeroVector(LocalSize);
//         Vector enrichment_terms_horizontal = ZeroVector(LocalSize);
//         double enrichment_diagonal = 0.0;
//         double enriched_rhs = 0.0;
        array_1d<double,3> bf = ZeroVector(3);

        double positive_volume = 0.0;
        double negative_volume = 0.0;
        //do integration
        for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
        {
            //assigning the gauss data
            for (unsigned int k = 0; k < TNumNodes; k++)
                N[k] = Ngauss(igauss, k);
            double wGauss = volumes[igauss];
            if(signs[igauss] > 0)
                positive_volume += wGauss;
            else
                negative_volume += wGauss;
            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            double ElemSize = this->ElementSize(Area);
            double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize, rCurrentProcessInfo);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Calculate stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne, TauTwo, AdvVel, ElemSize, Density, Viscosity, rCurrentProcessInfo);

            this->AddIntegrationPointVelocityContribution(rDampingMatrix, rRightHandSideVector, Density, Viscosity, AdvVel, TauOne, TauTwo, N, DN_DX, wGauss,Nenriched(igauss, 0),gauss_gradients[igauss]);
//             if (ndivisions > 1)
//             {
//                 //compute enrichment terms contribution
//                 for (unsigned int inode = 0; inode < TNumNodes; inode++)
//                 {
//                     int base_index = (TDim + 1) * inode;
//                     array_1d<double,TNumNodes> AGradN = ZeroVector(TNumNodes);
//                     this->GetConvectionOperator(AGradN,AdvVel,DN_DX);
//                     //momentum term
//                     for (unsigned int k = 0; k < TDim; k++)
//                     {
//                         double ConvTerm = wGauss * TauOne * gauss_gradients[igauss](0,k)* Density * AGradN[inode];
//                         enrichment_terms_vertical[base_index + k] += ConvTerm - wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
//                         enrichment_terms_horizontal[base_index + k] += ConvTerm + wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
// //                             enrichment_terms_vertical[base_index + k] +=wGauss*N[inode]*gauss_gradients[igauss](0, k); //-= wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
// //                            enrichment_terms_horizontal[base_index + k] -=Density*wGauss*N[inode]*gauss_gradients[igauss](0, k); //   += Density*wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
//                     }
//                     //pressure term
//                     for (unsigned int k = 0; k < TDim; k++)
//                     {
//                         double temp =  wGauss * TauOne* DN_DX(inode, k) * gauss_gradients[igauss](0, k);
//                         enrichment_terms_vertical[base_index + TDim] += temp;
//                         enrichment_terms_horizontal[base_index + TDim] += temp;
//                     }
//                     //add acceleration enrichment term
// 					//const array_1d<double,3>& vnode = this->GetGeometry()[inode].FastGetSolutionStepValue(VELOCITY);
//                     const array_1d<double,3>& old_vnode = this->GetGeometry()[inode].FastGetSolutionStepValue(VELOCITY,1);
//                     for (unsigned int k = 0; k < TDim; k++)
//                     {
//                         double coeff = wGauss * TauOne *Density *  gauss_gradients[igauss](0,k)*N[inode] * 2.0/ Dt;
//                         enrichment_terms_horizontal[base_index + k] += coeff;
//                         enriched_rhs += coeff * (old_vnode[k]);
// //                             enrichment_terms_vertical[base_index + k] +=wGauss*N[inode]*gauss_gradients[igauss](0, k); //-= wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
// //                            enrichment_terms_horizontal[base_index + k] -=Density*wGauss*N[inode]*gauss_gradients[igauss](0, k); //   += Density*wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
//                     }
//                 }
//                 //compute diagonal enrichment term
// 				array_1d<double,3> OldAcceleration = N[0]*this->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION,1);
// 				for(unsigned int jjj=0; jjj<(this->GetGeometry()).size(); jjj++)
// 					OldAcceleration += N[jjj]*(this->GetGeometry())[jjj].FastGetSolutionStepValue(ACCELERATION,1);
//
//                 for (unsigned int k = 0; k < TDim; k++)
//                 {
//                     const Matrix& enriched_grad = gauss_gradients[igauss];
//                     enrichment_diagonal += wGauss  * TauOne * pow(enriched_grad(0, k), 2);
//                     enriched_rhs += wGauss * TauOne *Density * enriched_grad(0,k)*(bf[k]+OldAcceleration[k]);
//                 }
//             }
         }

//         if (ndivisions > 1)
//         {
//             //add to LHS enrichment contributions
//             double inverse_diag_term = 1.0 / ( enrichment_diagonal);
//
//
//
//             for (unsigned int i = 0; i < LocalSize; i++)
//                 for (unsigned int j = 0; j < LocalSize; j++)
//                     rDampingMatrix(i, j) -= inverse_diag_term * enrichment_terms_vertical[i] * enrichment_terms_horizontal[j];
//             rRightHandSideVector -= (inverse_diag_term*enriched_rhs )*enrichment_terms_vertical;
//
//         }

	  // Now calculate an additional contribution to the residual: r -= rDampingMatrix * (u,p)
	  VectorType U = ZeroVector(LocalSize);
	  int LocalIndex = 0;
	  for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
	  {
	      array_1d< double, 3 > & rVel = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
	      for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
	      {
		  U[LocalIndex] = rVel[d];
		  ++LocalIndex;
	      }
	      U[LocalIndex] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
	      ++LocalIndex;
	  }
	  const double enriched_pr = this->GetValue(PRESSUREAUX);
	  U[LocalIndex] = enriched_pr;
	  noalias(rRightHandSideVector) -= prod(rDampingMatrix, U);
// 	  KRATOS_WATCH(enriched_pr);
	}
    }

    /// Implementation of FinalizeNonLinearIteration to compute enriched pressure.
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
//       this->IsCutted();
      if( this->is_cutted == 0)
	 ElementBaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);
      else
      {
      //fill vector of solution
      const unsigned int LocalSize = (TDim + 1) * TNumNodes;
      VectorType DU = ZeroVector(LocalSize);
      int LocalIndex = 0;
      for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
      {
	  const array_1d< double, 3 >  rVel = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
	  const array_1d< double, 3 >  old_rVel = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY,1);
	  for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
	  {
	      DU[LocalIndex] = rVel[d]-old_rVel[d];
	      ++LocalIndex;
	  }
	  DU[LocalIndex] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE) - this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE,1); // Pressure Dof
	  ++LocalIndex;
      }

      /*call Residual vector for condensation, it is filled in the
      schme and consists of last row of LHS plus last member of RHS*/
       VectorType residual_enr_vector = ZeroVector(LocalSize + 2);
       residual_enr_vector = this->GetValue(GAPS);

       //compute dp_enr = D^-1 * ( f - C U)
       double C_Dup(0.0);
       for (unsigned int ii = 0; ii < LocalSize; ++ii)
	  C_Dup += residual_enr_vector[ii]*DU[ii];

       //take the value of enriched pressure from the last step and add teh increment
       double enr_p =  this->GetValue(AUX_INDEX);
       if( residual_enr_vector[LocalSize] == 0.0 )
	 KRATOS_THROW_ERROR(std::invalid_argument,"Diagonla member of enriched RES is zero !!!!!","")
       else
	enr_p += (residual_enr_vector[LocalSize + 1] - C_Dup)/residual_enr_vector[LocalSize] ;

       //current iteration value of the enriched pressure
//        enr_p = 0.0;
       this->SetValue(PRESSUREAUX,enr_p);
      }
    }

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */

    void GetFirstDerivativesVector(Vector& values, int Step) override
    {
// 	this->IsCutted();
	if( this->is_cutted == 0)
	   ElementBaseType::GetFirstDerivativesVector(values, Step);
	else
	 {
	    unsigned int MatSize = (TDim + 1) * TNumNodes + 1;
	    if (values.size() != MatSize) values.resize(MatSize, false);
	    for (unsigned int i = 0; i < TNumNodes; i++)
	    {
		unsigned int index = i * (TDim + 1);
		values[index] =  this->GetGeometry()[i].GetSolutionStepValue(VELOCITY_X, Step);
		values[index + 1] = this->GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y, Step);
		values[index + 2] = this->GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z, Step);
		values[index + 3] = this->GetGeometry()[i].GetSolutionStepValue(PRESSURE, Step);

	    }
	    //add teh enriched component
	    int last_index = (TDim + 1) * TNumNodes;
	    values[last_index] = this->GetValue(PRESSUREAUX);
	  }
    }
    /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */

      void GetSecondDerivativesVector(Vector& values, int Step) override
      {
// 	this->IsCutted();
	if( this->is_cutted == 0)
	   ElementBaseType::GetSecondDerivativesVector(values, Step);
	else
	{
	  unsigned int MatSize = (TDim + 1) * TNumNodes + 1;
	  if (values.size() != MatSize) values.resize(MatSize, false);
	  for (unsigned int i = 0; i < TNumNodes; i++)
	  {
	      unsigned int index = i * (TDim + 1);
	      values[index] = this->GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X, Step);
	      values[index + 1] = this->GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y, Step);
	      values[index + 2] = this->GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z, Step);
	      values[index + 3] = 0.0;
	  }
	   //add teh enriched component
	    int last_index = (TDim + 1) * TNumNodes;
	    values[last_index] = 0.0;
	}
      }


    /// Implementation of Calculate to compute the local OSS projections.
    /**
     * If rVariable == ADVPROJ, This function computes the OSS projection
     * terms using pressure and velocity values from the previous iteration. The
     * projections are then added to the nodal variables ADVPROJ (Momentum residual)
     * and DIVPROJ (Mass continuity residual). It is assumed that the scheme will
     * divide the result by the assembled NODAL_AREA, which is equivalent to a
     * nodal interpolation using a lumped mass matrix.
     * @param rVariable Use ADVPROJ
     * @param Output Will be overwritten with the elemental momentum error
     * @param rCurrentProcessInfo Process info instance (unused)
     */
    void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                           array_1d<double, 3 > & rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == ADVPROJ) // Compute residual projections for OSS
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
            array_1d< double, 3 > ElementalMomRes = ZeroVector(3);
            double ElementalMassRes(0);
            //get position of the cut surface
            Vector distances(TNumNodes);
            Matrix Nenriched(6, 1);
            Vector volumes(6);
            Matrix coords(TNumNodes, TDim);
            Matrix Ngauss(6, TNumNodes);
            Vector signs(6);
            std::vector< Matrix > gauss_gradients(6);
            //fill coordinates
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
                volumes[i] = 0.0;
                distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
                for (unsigned int j = 0; j < TDim; j++)
                    coords(i, j) = xyz[j];
            }
            for (unsigned int i = 0; i < 6; i++)
                gauss_gradients[i].resize(1, TDim, false);

            array_1d<double,6> edge_areas;
            unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);
            //do integration
            for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
            {
                //assigning the gauss data
                for (unsigned int k = 0; k < TNumNodes; k++)
                    N[k] = Ngauss(igauss, k);
                double wGauss = volumes[igauss];
                // Calculate this element's fluid properties
                double Density;;
                this->EvaluateInPoint(Density, DENSITY, N);

                // Get Advective velocity
                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);
                // Output containers
                ElementalMomRes = ZeroVector(3);
                ElementalMassRes = 0.0;
                this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, wGauss);
                if (rCurrentProcessInfo[OSS_SWITCH] == 1)
                {
                    // Carefully write results to nodal variables, to avoid parallelism problems
                    for (unsigned int i = 0; i < TNumNodes; ++i)
                    {
                        this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
                        array_1d< double, 3 > & rAdvProj = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
                        for (unsigned int d = 0; d < TDim; ++d)
                            rAdvProj[d] += N[i] * ElementalMomRes[d];
                        this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += N[i] * ElementalMassRes;
                        this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += wGauss * N[i];
                        this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
                    }
                }
            }
            /// Return output
            rOutput = ElementalMomRes;
        }
        else if (rVariable == SUBSCALE_VELOCITY)
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
            array_1d< double, 3 > ElementalMomRes = ZeroVector(3);
            double ElementalMassRes(0);
            //get position of the cut surface
            Vector distances(TNumNodes);
            Matrix Nenriched(6, 1);
            Vector volumes(6);
            Matrix coords(TNumNodes, TDim);
            Matrix Ngauss(6, TNumNodes);
            Vector signs(6);
            std::vector< Matrix > gauss_gradients(6);
            //fill coordinates
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
                volumes[i] = 0.0;
                distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
                for (unsigned int j = 0; j < TDim; j++)
                    coords(i, j) = xyz[j];
            }
            for (unsigned int i = 0; i < 6; i++)
                gauss_gradients[i].resize(1, TDim, false);

            array_1d<double,6> edge_areas;
            unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);
            //do integration
            for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
            {
                //assigning the gauss data
                for (unsigned int k = 0; k < TNumNodes; k++)
                    N[k] = Ngauss(igauss, k);
                double wGauss = volumes[igauss];
                // Calculate this element's fluid properties
                double Density;
                this->EvaluateInPoint(Density, DENSITY, N);

                // Get Advective velocity
                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                // Output containers
                ElementalMomRes = ZeroVector(3);
                ElementalMassRes = 0.0;
                this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, wGauss);
                if (rCurrentProcessInfo[OSS_SWITCH] == 1)
                {
                    /* Projections of the elemental residual are computed with
                     * Newton-Raphson iterations of type M(lumped) dx = ElemRes - M(consistent) * x
                     */
                    const double Weight = ElementBaseType::ConsistentMassCoef(wGauss); // Consistent mass matrix is Weigth * ( Ones(TNumNodes,TNumNodes) + Identity(TNumNodes,TNumNodes) )
                    // Carefully write results to nodal variables, to avoid parallelism problems
                    for (unsigned int i = 0; i < TNumNodes; ++i)
                    {
                        this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
                        // Add elemental residual to RHS
                        array_1d< double, 3 > & rMomRHS = this->GetGeometry()[i].GetValue(ADVPROJ);
                        double& rMassRHS = this->GetGeometry()[i].GetValue(DIVPROJ);
                        for (unsigned int d = 0; d < TDim; ++d)
                            rMomRHS[d] += N[i] * ElementalMomRes[d];
                        rMassRHS += N[i] * ElementalMassRes;
                        // Write nodal area
                        this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += wGauss * N[i];
                        // Substract M(consistent)*x(i-1) from RHS
                        for (unsigned int j = 0; j < TNumNodes; ++j) // RHS -= Weigth * Ones(TNumNodes,TNumNodes) * x(i-1)
                        {
                            for (unsigned int d = 0; d < TDim; ++d)
                                rMomRHS[d] -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(ADVPROJ)[d];
                            rMassRHS -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ);
                        }
                        for (unsigned int d = 0; d < TDim; ++d) // RHS -= Weigth * Identity(TNumNodes,TNumNodes) * x(i-1)
                            rMomRHS[d] -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ)[d];
                        rMassRHS -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
                        this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
                    }
                }
            }
            /// Return output
            rOutput = ElementalMomRes;
        }
    }

/**
 * @see DPGVMS::GetValueOnIntegrationPoints
 */

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override
      {

	  if (rVariable == PRESSUREAUX)
	  {
	      for (unsigned int PointNumber = 0;
		      PointNumber < 1; PointNumber++)
	      {
		  //	KRATOS_WATCH(this->GetValue(IS_WATER));
		  //	KRATOS_WATCH(this->Info());
		  rValues[PointNumber] = this->GetValue(PRESSUREAUX);;
	      }
	  }
	  else if(rVariable == AUX_INDEX)
	  {
            double Area;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);
            double ElemSize = this->ElementSize(Area);

            rValues.resize(1, false);

            rValues[0] = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);
	  }

      }
    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if (ErrorCode != 0) return ErrorCode;
        // Check that all required variables have been registered
        if (DISTANCE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "DISTANCE Key is 0. Check if the application was correctly registered.", "");
        if (VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY Key is 0. Check if the application was correctly registered.", "");
        if (MESH_VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "MESH_VELOCITY Key is 0. Check if the application was correctly registered.", "");
        if (ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "ACCELERATION Key is 0. Check if the application was correctly registered.", "");
        if (PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "PRESSURE Key is 0. Check if the application was correctly registered.", "");
        if (DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY Key is 0. Check if the application was correctly registered.", "");
        if (VISCOSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "VISCOSITY Key is 0. Check if the application was correctly registered.", "");
        if (OSS_SWITCH.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "OSS_SWITCH Key is 0. Check if the application was correctly registered.", "");
        if (DYNAMIC_TAU.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "DYNAMIC_TAU Key is 0. Check if the application was correctly registered.", "");
        if (DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "DELTA_TIME Key is 0. Check if the application was correctly registered.", "");
        if (ADVPROJ.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "ADVPROJ Key is 0. Check if the application was correctly registered.", "");
        if (DIVPROJ.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "DIVPROJ Key is 0. Check if the application was correctly registered.", "");
        if (NODAL_AREA.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "NODAL_AREA Key is 0. Check if the application was correctly registered.", "");
        if (C_SMAGORINSKY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "C_SMAGORINSKY Key is 0. Check if the application was correctly registered.", "");
        if (ERROR_RATIO.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "ERROR_RATIO Key is 0. Check if the application was correctly registered.", "");
        // Additional variables, only required to print results:
        // SUBSCALE_VELOCITY, SUBSCALE_PRESSURE, TAUONE, TAUTWO, MU, VORTICITY.
        // Checks on nodes
        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for (unsigned int i = 0; i<this->GetGeometry().size(); ++i)
        {
            if (this->GetGeometry()[i].SolutionStepsDataHas(DISTANCE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing DISTANCE variable on solution step data for node ", this->GetGeometry()[i].Id());
            if (this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY variable on solution step data for node ", this->GetGeometry()[i].Id());
            if (this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE variable on solution step data for node ", this->GetGeometry()[i].Id());
            if (this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing MESH_VELOCITY variable on solution step data for node ", this->GetGeometry()[i].Id());
            if (this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing ACCELERATION variable on solution step data for node ", this->GetGeometry()[i].Id());
            if (this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                    this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                    this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY component degree of freedom on node ", this->GetGeometry()[i].Id());
            if (this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE component degree of freedom on node ", this->GetGeometry()[i].Id());
        }
        // Not checking OSS related variables NODAL_AREA, ADVPROJ, DIVPROJ, which are only required as SolutionStepData if OSS_SWITCH == 1
        // If this is a 2D problem, check that nodes are in XY plane
        if (this->GetGeometry().WorkingSpaceDimension() == 2)
        {
            for (unsigned int i = 0; i<this->GetGeometry().size(); ++i)
            {
                if (this->GetGeometry()[i].Z() != 0.0)
                    KRATOS_THROW_ERROR(std::invalid_argument, "Node with non-zero Z coordinate found. Id: ", this->GetGeometry()[i].Id());
            }
        }
        return 0;
        KRATOS_CATCH("");
    }
    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Elemental Data
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
        buffer << "DPGVMS #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DPGVMS" << TDim << "D";
    }
    //        /// Print object's data.
    //        virtual void PrintData(std::ostream& rOStream) const;
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
    /// Add lumped mass matrix
    void LumpMassMatrix(MatrixType& rMassMatrix)
    {
        for (unsigned int i = 0; i < rMassMatrix.size1(); ++i)
        {
            double diag_factor = 0.0;
            for (unsigned int j = 0; j < rMassMatrix.size2(); ++j)
            {
                diag_factor += rMassMatrix(i, j);
                rMassMatrix(i, j) = 0.0;
            }
            rMassMatrix(i, i) = diag_factor;
        }
    }
    /// Add the weighted value of a variable at a point inside the element to a double
    /**
     * Evaluate a scalar variable in the point where the form functions take the
     * values given by rShapeFunc and add the result, weighted by Weight, to
     * rResult. This is an auxiliary function used to compute values in integration
     * points.
     * @param rResult The double where the value will be added to
     * @param rVariable The nodal variable to be read
     * @param rShapeFunc The values of the form functions in the point
     * @param Step The time Step (Defaults to 0 = Current)
     * @param Weight The variable will be weighted by this value before it is added to rResult
     */
    virtual void AddPointContribution(double& rResult,
                                      const Variable< double >& rVariable,
                                      const array_1d< double, TNumNodes >& rShapeFunc,
                                      const double Weight = 1.0)
    {
        double temp = 0.0;
        this->EvaluateInPoint(temp, rVariable, rShapeFunc);
        rResult += Weight*temp;
    }
    /// Write the value of a variable at a point inside the element to a double
    /**
     * Evaluate a scalar variable in the point where the form functions take the
     * values given by rShapeFunc and write the result to rResult.
     * This is an auxiliary function used to compute values in integration points.
     * @param rResult The double where the value will be added to
     * @param rVariable The nodal variable to be read
     * @param rShapeFunc The values of the form functions in the point
     * @param Step The time Step (Defaults to 0 = Current)
     */
    void EvaluateInPoint(double& rResult,
                                 const Variable< double >& rVariable,
                                 const array_1d< double, TNumNodes >& rShapeFunc) override
    {
        //compute sign of distance on gauss point
        double dist = 0.0;
        for (unsigned int i = 0; i < TNumNodes; i++)
            dist += rShapeFunc[i] * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);

        double navg = 0.0;
        double value = 0.0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            if ( (dist * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)) > 0.0)
            {
                navg += 1.0;
                value += this->GetGeometry()[i].FastGetSolutionStepValue(rVariable);
            }
        }
        if(navg != 0)
            value /= navg;
        else
            KRATOS_THROW_ERROR(std::invalid_argument,"IT IS A CUTTED ELEMENT navg MUST BE NON ZERO !!!!","");
        rResult = value;
//            if(rVariable==DENSITY)
//                KRATOS_WATCH(rResult)
    }
    /// Add the weighted value of a variable at a point inside the element to a vector
    /**
     * Evaluate a vector variable in the point where the form functions take the
     * values given by rShapeFunc and add the result, weighted by Weight, to
     * rResult. This is an auxiliary function used to compute values in integration
     * points.
     * @param rResult The vector where the value will be added to
     * @param rVariable The nodal variable to be read
     * @param rShapeFunc The values of the form functions in the point
     * @param Weight The variable will be weighted by this value before it is added to rResult
     */
    virtual void AddPointContribution(array_1d< double, 3 > & rResult,
                                      const Variable< array_1d< double, 3 > >& rVariable,
                                      const array_1d< double, TNumNodes>& rShapeFunc,
                                      const double Weight = 1.0)
    {
        array_1d<double, 3 > temp = ZeroVector(3);
        this->EvaluateInPoint(temp, rVariable, rShapeFunc);
        rResult += Weight*temp;
    }
    /// Write the value of a variable at a point inside the element to a double
    /**
     * Evaluate a scalar variable in the point where the form functions take the
     * values given by rShapeFunc and write the result to rResult.
     * This is an auxiliary function used to compute values in integration points.
     * @param rResult The double where the value will be added to
     * @param rVariable The nodal variable to be read
     * @param rShapeFunc The values of the form functions in the point
     */
    void EvaluateInPoint(array_1d< double, 3 > & rResult,
                                 const Variable< array_1d< double, 3 > >& rVariable,
                                 const array_1d< double, TNumNodes >& rShapeFunc) override
    {
        //compute sign of distance on gauss point
        double dist = 0.0;
        for (unsigned int i = 0; i < TNumNodes; i++)
            dist += rShapeFunc[i] * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        double navg = 0.0;
        array_1d< double, 3 > value = ZeroVector(3);
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            if (dist * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE) > 0.0)
            {
                navg += 1.0;
                value += this->GetGeometry()[i].FastGetSolutionStepValue(rVariable);
            }
        }
        if(navg != 0)
            value /= navg;
        else
            ElementBaseType::EvaluateInPoint(value,rVariable,rShapeFunc);
        rResult = value;
    }

//     virtual void IsCutted(){
//       int negative = 0;
//        this->is_cutted = 0;
//       for( int ii = 0; ii<TNumNodes; ++ii)
//         if( this->GetGeometry()[ii].FastGetSolutionStepValue(DISTANCE) < 0.0)
// 	  negative++;
//
// 	if( negative != TNumNodes && negative != 0)
// 	  this->is_cutted = 1;
//
//     }
     /// Add mass-like stabilization terms to LHS.
    /**
     * This function is only used in ASGS. For OSS, we avoid computing these
     * terms, as they shoud cancel out with the dynamic part of the projection
     * (which is not computed either)
     * @param rLHSMatrix Left hand side of the velocity-pressure system
     * @param Density Density on integration point
     * @param rAdvVel Advective velocity on integration point
     * @param TauOne Stabilization parameter for momentum equation
     * @param rShapeFunc Shape funcitions evaluated on integration point
     * @param rShapeDeriv Shape function derivatives evaluated on integration point
     * @param Weight Area (or volume) times integration point weight
     */
    void AddMassStabTerms(MatrixType& rLHSMatrix,
                          const double Density,
                          const array_1d<double, 3 > & rAdvVel,
                          const double TauOne,
                          const array_1d<double, TNumNodes>& rShapeFunc,
                          const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                          const double Weight,
			  const Matrix gauss_enriched_gradients)
    {
        const unsigned int BlockSize = TDim + 1;

        double Coef = Weight * TauOne;
        unsigned int FirstRow(0), FirstCol(0);
        double K; // Temporary results

        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        // Note: Dof order is (vx,vy,[vz,]p) for each node
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            // Loop over columns
            for (unsigned int j = 0; j < TNumNodes; ++j)
            {
                // Delta(u) * TauOne * [ AdvVel * Grad(v) ] in velocity block
                K = Coef * Density * AGradN[i] * rShapeFunc[j];

                for (unsigned int d = 0; d < TDim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                {
                    rLHSMatrix(FirstRow + d, FirstCol + d) += K;
                    // Delta(u) * TauOne * Grad(q) in q * Div(u) block
                    rLHSMatrix(FirstRow + TDim, FirstCol + d) += Coef * Density * rShapeDeriv(i, d) * rShapeFunc[j];
                }
                // Update column index
                FirstCol += BlockSize;
            }
            // Update matrix indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }

        //add Delta(u) * TauOne * Grad(q_ENR)
        // Loop over columns
        for (unsigned int j = 0; j < TNumNodes; ++j)
            {
                for (unsigned int d = 0; d < TDim; ++d)
                {
	          rLHSMatrix(FirstRow , FirstCol + d) +=  Coef * Density * gauss_enriched_gradients(0,d) * rShapeFunc[j];
		}
	      FirstCol += BlockSize;
	    }
    }
    /// Add a the contribution from a single integration point to the velocity contribution
    void AddIntegrationPointVelocityContribution(MatrixType& rDampingMatrix,
            VectorType& rDampRHS,
            const double Density,
            const double Viscosity,
            const array_1d< double, 3 > & rAdvVel,
            const double TauOne,
            const double TauTwo,
            const array_1d< double, TNumNodes >& rShapeFunc,
            const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
            const double Weight,
            const double gauss_N_en,
	    Matrix& gauss_enriched_gradients)
    {
        const unsigned int BlockSize = TDim + 1;

        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        // Build the local matrix and RHS
        unsigned int FirstRow(0), FirstCol(0); // position of the first term of the local matrix that corresponds to each node combination
        double K, G, PDivV, L, qF; // Temporary results

        // Note that we iterate first over columns, then over rows to read the Body Force only once per node
        for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over colums
        {
            // Get Body Force
            const array_1d<double, 3 > & rBodyForce = this->GetGeometry()[j].FastGetSolutionStepValue(BODY_FORCE);

            for (unsigned int i = 0; i < TNumNodes; ++i) // iterate over rows
            {
                // Calculate the part of the contributions that is constant for each node combination

                // Velocity block
                K = Density * rShapeFunc[i] * AGradN[j]; // Convective term: v * ( a * Grad(u) )
                K += TauOne * Density * AGradN[i] * Density * AGradN[j]; // Stabilization: (a * Grad(v)) * TauOne * (a * Grad(u))
                K *= Weight;

                // q-p stabilization block (reset result)
                L = 0;

                for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
                {
                    // Velocity block
//                        K += Weight * Viscosity * rShapeDeriv(i, m) * rShapeDeriv(j, m); // Diffusive term: Viscosity * Grad(v) * Grad(u)

                    // v * Grad(p) block
                    G = TauOne * Density * AGradN[i] * rShapeDeriv(j, m); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                    PDivV = rShapeDeriv(i, m) * rShapeFunc[j]; // Div(v) * p

                    // Write v * Grad(p) component
                    rDampingMatrix(FirstRow + m, FirstCol + TDim) += Weight * (G - PDivV);
                    // Use symmetry to write the q * Div(u) component
                     rDampingMatrix(FirstCol + TDim, FirstRow + m) += Weight * (G + PDivV);
//		    rDampingMatrix(FirstCol + TDim, FirstRow + m) += Weight * (G - rShapeDeriv(j, m) * rShapeFunc[i]);

                    // q-p stabilization block
                    L += rShapeDeriv(i, m) * rShapeDeriv(j, m); // Stabilization: Grad(q) * TauOne * Grad(p)

                    for (unsigned int n = 0; n < TDim; ++n) // iterate over u components (ux,uy[,uz])
                    {
                        // Velocity block
                        rDampingMatrix(FirstRow + m, FirstCol + n) += Weight * TauTwo * rShapeDeriv(i, m) * rShapeDeriv(j, n); // Stabilization: Div(v) * TauTwo * Div(u)
                    }

                }

                // Write remaining terms to velocity block
                for (unsigned int d = 0; d < TDim; ++d)
                    rDampingMatrix(FirstRow + d, FirstCol + d) += K;

                // Write q-p stabilization block
                rDampingMatrix(FirstRow + TDim, FirstCol + TDim) += Weight * TauOne * L;

                // Operate on RHS
                qF = 0.0;
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    rDampRHS[FirstRow + d] += Weight * TauOne * Density * AGradN[i] * rShapeFunc[j] * Density * rBodyForce[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
                    qF += rShapeDeriv(i, d) * rShapeFunc[j] * rBodyForce[d];
                }
                rDampRHS[FirstRow + TDim] += Weight * Density * TauOne * qF; // Grad(q) * TauOne * (Density * BodyForce)

                // Update reference row index for next iteration
                FirstRow += BlockSize;
            }

            // Update reference indices
            FirstRow = 0;
            FirstCol += BlockSize;
        }

//            this->AddBTransCB(rDampingMatrix,rShapeDeriv,Viscosity*Weight);
        this->AddViscousTerm(rDampingMatrix,rShapeDeriv,Viscosity*Weight);


	//add enrichment terms
	for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over colums
        {
            // Get Body Force
            const array_1d<double, 3 > & rBodyForce = this->GetGeometry()[j].FastGetSolutionStepValue(BODY_FORCE);
	    L = 0.0;
	    qF = 0.0;
	    for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
                {
                    // v * Grad(p) block
                    G = TauOne * Density * AGradN[j] * gauss_enriched_gradients(0,m); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                    PDivV = rShapeDeriv(j, m) * gauss_N_en; // Div(v) * p
                    double VGradP = -gauss_enriched_gradients(0,m)*rShapeFunc[j]; // Grad(p_star) * v

                    // Write v * Grad(p) component
                    rDampingMatrix(FirstRow + m, FirstCol) += Weight * (G - VGradP);
                    // Use symmetry to write the q * Div(u) component
                    rDampingMatrix(FirstCol , FirstRow + m) += Weight * (G + PDivV);//Weight * (G + PDivV);

                    // q-p stabilization block
                    L += rShapeDeriv(j, m) * gauss_enriched_gradients(0,m); // Stabilization: Grad(q) * TauOne * Grad(p)

                    // grad_q*body_force
                    qF += gauss_enriched_gradients(0,m) * rShapeFunc[j] * rBodyForce[m];
		}

                // Write q-p stabilization block
                rDampingMatrix(FirstRow + TDim, FirstCol ) += Weight * TauOne * L;
                rDampingMatrix(FirstCol, FirstRow + TDim ) += Weight * TauOne * L;

		// Operate on RHS
                rDampRHS[FirstCol] += Weight * Density * TauOne * qF; // Grad(q) * TauOne * (Density * BodyForce)

                // Update reference indices
		FirstRow += BlockSize;
       }
	//Write q_enr-p_enr stabilization
	for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
	      rDampingMatrix(FirstCol, FirstCol ) += Weight * TauOne * gauss_enriched_gradients(0,m) *  gauss_enriched_gradients(0,m);

    }

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
     unsigned int is_cutted;
    ///@}
    ///@name Member Variables
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ElementBaseType);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ElementBaseType);
    }
    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{
void AddBoundaryTerm(BoundedMatrix<double, 16, 16 >& rDampingMatrix,
                              const BoundedMatrix<double,4,3>& rShapeDeriv,
		              const array_1d<double, 4>& N_shape,
			      const array_1d<double,3>& nn,
                              const double& Weight,
                     const double ElementVolume,
			      ProcessInfo& rCurrentProcessInfo)
{

    const double OneThird = 1.0 / 3.0;

    double Density,Viscosity;
    ElementBaseType::EvaluateInPoint(Density, DENSITY, N_shape);
    // NODAL viscosity (no Smagorinsky/non-newtonian behaviour)
    ElementBaseType::EvaluateInPoint(Viscosity,VISCOSITY,N_shape);

    double viscos_weight = Weight *  OneThird * Viscosity * Density;

    unsigned int FirstRow(0),FirstCol(0);
    for (unsigned int j = 0; j < TNumNodes; ++j)
    {
      double jj_nd_flag = this->GetGeometry()[j].FastGetSolutionStepValue(FLAG_VARIABLE);
      if(jj_nd_flag == 5.0){
	//pressrue terms (three gauss points)
	array_1d<double,3> node_normal = this->GetGeometry()[j].FastGetSolutionStepValue(NORMAL);
	node_normal /= norm_2(node_normal);
// 	double dot_prod = MathUtils<double>::Dot(node_normal,nn);

// 	rDampingMatrix(FirstRow,FirstRow+3) = Weight *  OneThird * (nn[0]);// - dot_prod * node_normal[0]);//nn[0];
// 	rDampingMatrix(FirstRow+1,FirstRow+3) = Weight *  OneThird * (nn[1]);// - dot_prod * node_normal[1]);//nn[1];
// 	rDampingMatrix(FirstRow+2,FirstRow+3) = Weight *  OneThird * (nn[2]);// - dot_prod * node_normal[2]);//nn[2];

	  for (unsigned int i = 0; i < TNumNodes; ++i)
	  {
	    double ii_nd_flag = this->GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE);
	    if(ii_nd_flag == 5.0){
	      // nxdn/dx + nydn/dy + nz dn/dz
	      const double Diag =  nn[0]*rShapeDeriv(i,0) + nn[1]*rShapeDeriv(i,1) + nn[2]*rShapeDeriv(i,2);

	      // First Row
	      rDampingMatrix(FirstRow,FirstCol) = -(viscos_weight *  ( nn[0]*rShapeDeriv(i,0) + Diag ) );// * nn[0]*nn[0];
	      rDampingMatrix(FirstRow,FirstCol+1) = -(viscos_weight * nn[1]*rShapeDeriv(i,0) );//* nn[0]*nn[1] ;
	      rDampingMatrix(FirstRow,FirstCol+2) = -(viscos_weight * nn[2]*rShapeDeriv(i,0) );//* nn[0]*nn[2];;


	      // Second Row
	      rDampingMatrix(FirstRow+1,FirstCol) = -(viscos_weight * nn[0]*rShapeDeriv(i,1) );//* nn[1]*nn[0] ;
	      rDampingMatrix(FirstRow+1,FirstCol+1) = -(viscos_weight * ( nn[1]*rShapeDeriv(i,1) + Diag ) );//*nn[1]*nn[1] ;
	      rDampingMatrix(FirstRow+1,FirstCol+2) = -(viscos_weight * nn[2]*rShapeDeriv(i,1) );// * nn[1]*nn[2];


	      // Third Row
	      rDampingMatrix(FirstRow+2,FirstCol) = -(viscos_weight * nn[0]*rShapeDeriv(i,2) );// * nn[2]*nn[0];
	      rDampingMatrix(FirstRow+2,FirstCol+1) = -(viscos_weight * nn[1]*rShapeDeriv(i,2) );// * nn[2]*nn[1];
	      rDampingMatrix(FirstRow+2,FirstCol+2) = -(viscos_weight * ( nn[2]*rShapeDeriv(i,2) + Diag ) );//* nn[2]*nn[2] ;
	    }

	     // Update Counter
	      FirstCol += 4;
	  }
      }
	  FirstCol = 0;
	  FirstRow += 4;
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
    /// Assignment operator.
    DPGVMS & operator=(DPGVMS const& rOther);
    /// Copy constructor.
    DPGVMS(DPGVMS const& rOther);
    ///@}
}; // Class DPGVMS
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::istream & operator >>(std::istream& rIStream,
                                  DPGVMS<TDim, TNumNodes>& rThis)
{
    return rIStream;
}
/// output stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const DPGVMS<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}
///@} // Fluid Dynamics Application group
} // namespace Kratos.
#endif // KRATOS_TWO_FLUID_DPGVMSS_H_INCLUDED  defined
