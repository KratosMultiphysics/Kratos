//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Jordi Cotela
//

#if !defined(KRATOS_TWO_FLUID_VMS_H_INCLUDED )
#define  KRATOS_TWO_FLUID_VMS_H_INCLUDED
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
#include "utilities/math_utils.h"
#include "utilities/split_tetrahedra.h"
// #include "utilities/enrichment_utilities.h"
#include "utilities/enrichment_utilities_duplicate_dofs.h"
// Application includes
#include "fluid_dynamics_application_variables.h"
#include "vms.h"

namespace Kratos
{

/**The "TwoFluidVMS" element is an element based on the Variation Multiscale Stabilization technique (VMS)
 * which is designed for the solution of a two fluid problem.
 *
 * A distinctive feature of the element is the use of 4 LOCAL enrichment functions, which allows to model
 * a discontinuity in both the pressure field and in its gradient.
 * The enrichment functions are obtained by duplicating all of the degrees of freedom of the element.
 * Since the enrichment is performed elementwise, a purely local static condensation
 * step is performed.
 *
 * Since a jump in the pressure can be considered, the element shall be able to habdle moderate changes of the viscosity
 * between the two fluids to be considered
 *
 * WARNING: From the implementation point of view, the element hard codes a BDF2 scheme within the element
 * this is different from the VMS base element which supports the use of an arbitrary time integrator.
 * In the practice this implies that the element can ONLY be used in conjunction with the scheme implemented
 * in "residualbased_predictorcorrector_velocity_bdf_scheme_turbulent.h"
 *
 * a buffer size of dimension 3 is needed since the current step and two steps in the past need to be stored.
 *
 * TODO: so far only ASGS stabilization is implemented. OSS stabilization is possible but has not yet been implemented
 *
 *
 *
 */
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
class TwoFluidVMS : public VMS<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition of TwoFluidVMS
    KRATOS_CLASS_POINTER_DEFINITION(TwoFluidVMS);
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
    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;
    ///@}
    ///@name Life Cycle
    ///@{
    //Constructors.
    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    TwoFluidVMS(IndexType NewId = 0) :
        ElementBaseType(NewId)
    {
    }
    ///Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    TwoFluidVMS(IndexType NewId, const NodesArrayType& ThisNodes) :
        ElementBaseType(NewId, ThisNodes)
    {
    }
    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    TwoFluidVMS(IndexType NewId, GeometryType::Pointer pGeometry) :
        ElementBaseType(NewId, pGeometry)
    {
    }
    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    TwoFluidVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        ElementBaseType(NewId, pGeometry, pProperties)
    {
    }
    /// Destructor.
    ~TwoFluidVMS() override
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
     * Returns a pointer to a new TwoFluidVMS element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< TwoFluidVMS >(NewId, (this->GetGeometry()).Create(ThisNodes), pProperties);
    }
    Element::Pointer Create(IndexType NewId,
                           GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< TwoFluidVMS >(NewId, pGeom, pProperties);
    }

    /// Provides local contributions from body forces to the RHS
    /**
     * This is called during the assembly process and provides the RHS terms of the
     * system that are either constant or computed explicitly (from the 'old'
     * iteration variables). In this case this means the body force terms
     * @param rRightHandSideVector Will be filled with the elemental right hand side
     * @param rCurrentProcessInfo ProcessInfo instance from the ModelPart. It is
     * expected to contain values for DYNAMIC_TAU and DELTA_TIME
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int local_size = (TDim+1)*(TDim+1);
        Matrix tmp(local_size,local_size);
        CalculateLocalSystem(tmp,rRightHandSideVector,rCurrentProcessInfo);
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
        KRATOS_THROW_ERROR(std::logic_error,"MassMatrix function shall not be called when using this type of element","");
    }


    /// Calculate the element's local contribution to the system for the current step.
    /// this function is essentially identical to the one of the father element, to which it only
    /// adds a term in the momentum equation to allow imposing weakly the tangential component of the velocity
    /// on the cut elements
        void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

        const ProcessInfo& rConstProcessInfo = rCurrentProcessInfo; // Taking const reference for thread safety

        //****************************************************
        // Resize and set to zero the RHS
        if(rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize,false);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

        // Resize and set to zero the LHS
        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

       Matrix MassMatrix = ZeroMatrix(LocalSize, LocalSize);

       //****************************************************
        //Get Vector of BDF coefficients
        const Vector& BDFVector = rConstProcessInfo[BDF_COEFFICIENTS];

       //****************************************************
        // Get this element's geometric properties
        double Area;
        array_1d<double, TNumNodes> N;
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);


        //estimate a minimal h
        /*double h=0.0;
        if(TDim == 3) h = pow(6.0*Area, 1.0/3.0);
        else h = sqrt(2.0*Area);*/



        //input data for enrichment function
        Vector distances(TNumNodes);
        Matrix coords(TNumNodes, TDim);

        //output data for enrichment function
        Matrix Nenriched;
        Vector volumes;
        Matrix Ngauss;
        Vector signs(6);
        std::vector< Matrix > gauss_gradients;


        //fill coordinates
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
//             volumes[i] = 0.0;
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            for (unsigned int j = 0; j < TDim; j++)
                coords(i, j) = xyz[j];
        }

/*        for (unsigned int i = 0; i < 6; i++)
            gauss_gradients[i].resize(1, TDim, false);
*/
//         unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
         unsigned int ndivisions = EnrichmentUtilitiesDuplicateDofs::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
        const unsigned int nenrichments = Nenriched.size2();

        Matrix enrichment_terms_vertical   = ZeroMatrix(LocalSize,nenrichments);
        Matrix enrichment_terms_horizontal = ZeroMatrix(nenrichments,LocalSize);

        Matrix enrichment_diagonal = ScalarMatrix(nenrichments,nenrichments,0.0);
        Vector enriched_rhs(nenrichments,0.0);
        array_1d<double,3> bf = ZeroVector(3);

        double positive_volume = 0.0;
        double negative_volume = 0.0;

        if(ndivisions == 1) //compute gauss points for exact integration of a tetra element
        {
            const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints(GeometryData::GI_GAUSS_2);

            Ngauss = this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

            volumes.resize(IntegrationPoints.size(),false);

            for (unsigned int g = 0; g < this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
	    {
                volumes[g] = 6.0*Area * IntegrationPoints[g].Weight();
		        signs[g] = signs[0];
	    }

        }
//         else
//         {
//             Matrix aux = Ngauss;
//             Ngauss.resize(ndivisions,Ngauss.size2(),false);
//             for(unsigned int i=0; i<ndivisions; i++)
//                 for(unsigned int k=0; k<Ngauss.size2(); k++)
//                     Ngauss(i,k) = aux(i,k);
//         }
/*  KRATOS_WATCH(this->Id());
KRATOS_WATCH(volumes);
KRATOS_WATCH(Ngauss);  */

        for (unsigned int igauss = 0; igauss < Ngauss.size1(); igauss++)
        {
            double wGauss = volumes[igauss];


            if(signs[igauss] > 0) //check positive and negative volume
                positive_volume += wGauss;
            else
                negative_volume += wGauss;
        }


        const double min_area_ratio = 1e-6;
//         if(positive_volume/Area < min_area_ratio)
//         {
//             for (unsigned int igauss = 0; igauss < Ngauss.size1(); igauss++)
//             {
//                  if(signs[igauss] > 0) //check positive and negative volume
//                     volumes[igauss] =0.0;
//             }
//         }
//
//         if(negative_volume/Area < min_area_ratio)
//         {
//             for (unsigned int igauss = 0; igauss < Ngauss.size1(); igauss++)
//             {
//                  if(signs[igauss] < 0) //check positive and negative volume
//                     volumes[igauss] =0.0;
//             }
//         }

        // Porous media losses
        const Properties& r_properties = this->GetProperties();
        const double c1 = r_properties[LIN_DARCY_COEF];
        const double c2 = r_properties[NONLIN_DARCY_COEF];


        //****************************************************
        //compute LHS and RHS + first part of mass computation
        for (unsigned int igauss = 0; igauss < Ngauss.size1(); igauss++)
        {
            //assigning the gauss data
            for (unsigned int k = 0; k < TNumNodes; k++)
                N[k] = Ngauss(igauss, k);
            double wGauss = volumes[igauss];


//             if(signs[igauss] > 0) //check positive and negative volume
//                 positive_volume += wGauss;
//             else
//                 negative_volume += wGauss;

            //****************************************************
            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            double ElemSize = this->ElementSize(Area);
            double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

            //compute RHS contributions
            this->AddMomentumRHS(rRightHandSideVector, Density, N, wGauss);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);
            const double VelNorm = MathUtils<double>::Norm3(AdvVel);

            const double DarcyTerm = this->CalculateDarcyTerm(Density, Viscosity, c1, c2, N);
            // Calculate stabilization parameters
            double TauOne, TauTwo;

            //compute stabilization parameters
            this->CalculateStabilizationTau(TauOne, TauTwo, VelNorm, ElemSize, Density, Viscosity, DarcyTerm, rCurrentProcessInfo);

            this->AddIntegrationPointVelocityContribution(rLeftHandSideMatrix, rRightHandSideVector, Density, Viscosity, AdvVel, DarcyTerm, TauOne, TauTwo, N, DN_DX, wGauss);

            //compute mass matrix - terms related to real mass
            this->AddConsistentMassMatrixContribution(MassMatrix, N, Density, wGauss);


            //****************************************************
            //enrichment variables
            if (ndivisions > 1)
            {
//                 KRATOS_WATCH(Nenriched);
                this->EvaluateInPoint(bf, BODY_FORCE, N);

                //note that here we compute only a part of the acceleration term
                //this is done like this since the velocity*BDFVector[0] is treated implicitly
                array_1d<double,3> OldAcceleration = ZeroVector(3);
                for(unsigned int step=1; step<BDFVector.size(); step++)
                {
                    for(unsigned int jjj=0; jjj<(this->GetGeometry()).size(); jjj++)
                        OldAcceleration += N[jjj] * BDFVector[step] * (this->GetGeometry())[jjj].FastGetSolutionStepValue(VELOCITY,step);
                }

                for(unsigned int enriched_id = 0; enriched_id < nenrichments; enriched_id++)
                {
                    const Matrix& enriched_grad = gauss_gradients[igauss];
//                     KRATOS_WATCH(enriched_grad);


                    //compute enrichment terms contribution
                    for (unsigned int inode = 0; inode < TNumNodes; inode++)
                    {
                        int base_index = (TDim + 1) * inode;
                        array_1d<double,TNumNodes> AGradN = ZeroVector(TNumNodes);
                        this->GetConvectionOperator(AGradN,AdvVel,DN_DX);
                        //momentum term
                        for (unsigned int k = 0; k < TDim; k++)
                        {
                            double convection_stab = wGauss * TauOne * enriched_grad(enriched_id,k)* Density * AGradN[inode];
                            double darcy_stab = wGauss * TauOne * enriched_grad(enriched_id,k) * DarcyTerm * N[enriched_id];

                            //                      enrichment_terms_vertical[base_index + k] += velocity_stab + wGauss*N[inode]*enriched_grad(0, k);
                            enrichment_terms_vertical(base_index + k,enriched_id) += convection_stab - darcy_stab - wGauss * DN_DX(inode, k) * Nenriched(igauss, enriched_id);
                            enrichment_terms_horizontal(enriched_id,base_index + k) += convection_stab + darcy_stab + wGauss * DN_DX(inode, k) * Nenriched(igauss, enriched_id);
    //                             enrichment_terms_vertical[base_index + k] +=wGauss*N[inode]*enriched_grad(0, k); //-= wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
    //                            enrichment_terms_horizontal[base_index + k] -=Density*wGauss*N[inode]*enriched_grad(0, k); //   += Density*wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
                        }
                        //pressure term
                        for (unsigned int k = 0; k < TDim; k++)
                        {
                            double temp =  wGauss * TauOne* DN_DX(inode, k) * enriched_grad(enriched_id, k);
                            enrichment_terms_vertical(base_index + TDim,enriched_id) += temp;
                            enrichment_terms_horizontal(enriched_id,base_index + TDim) += temp;
                        }
                        //add acceleration enrichment term
                        for (unsigned int k = 0; k < TDim; k++)
                        {
                            double coeff = wGauss * TauOne *Density *  enriched_grad(enriched_id,k)*N[inode] * BDFVector[0];
                            enrichment_terms_horizontal(enriched_id,base_index + k) += coeff;
                            //i believe this shall not be here!! enriched_rhs += coeff * (old_vnode[k]);
    //                             enrichment_terms_vertical[base_index + k] +=wGauss*N[inode]*gauss_gradients[igauss](0, k); //-= wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
    //                            enrichment_terms_horizontal[base_index + k] -=Density*wGauss*N[inode]*gauss_gradients[igauss](0, k); //   += Density*wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
                        }
                    }
                    //compute diagonal enrichment terms

                    for (unsigned int k = 0; k < TDim; k++)
                    {
                        for(unsigned int lll=0; lll<nenrichments; lll++)
                        {
                            enrichment_diagonal(enriched_id,lll) += wGauss  * TauOne * enriched_grad(enriched_id, k) * enriched_grad(lll,k);
                        }


                        enriched_rhs[enriched_id] += wGauss * TauOne *Density * enriched_grad(enriched_id,k)*(bf[k]-OldAcceleration[k]); //changed the sign of the acc term
                    }
                }
            }
        }


//         if (ndivisions > 1)
//         {
//             KRATOS_WATCH(this->Id());
//             KRATOS_WATCH(ndivisions);
//         KRATOS_WATCH( (positive_volume+negative_volume )/Area);
//         }
//    KRATOS_WATCH("line 438");
        //lump mass matrix
        this->LumpMassMatrix(MassMatrix);

        //add mass matrix stabilization contributions
        for (unsigned int igauss = 0; igauss < Ngauss.size1(); igauss++)
        {
            //assigning the gauss data
            for (unsigned int k = 0; k < TNumNodes; k++)
                N[k] = Ngauss(igauss, k);
            double wGauss = volumes[igauss];
            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            double ElemSize = this->ElementSize(Area);
            double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);
            const double VelNorm = MathUtils<double>::Norm3(AdvVel);

            const double DarcyTerm = this->CalculateDarcyTerm(Density, Viscosity, c1, c2, N);

            double TauOne,TauTwo;
            this->CalculateStabilizationTau(TauOne, TauTwo, VelNorm, ElemSize, Density, Viscosity, DarcyTerm, rCurrentProcessInfo);

            // Add dynamic stabilization terms ( all terms involving a delta(u) )
            this->AddMassStabTerms(MassMatrix, Density, AdvVel, DarcyTerm, TauOne, N, DN_DX, wGauss);

        }

        //****************************************************
        //consider contributions of mass to LHS and RHS
        //add Mass Matrix to the LHS with the correct coefficient
        noalias(rLeftHandSideMatrix) += BDFVector[0]*MassMatrix;

        //do RHS -= MassMatrix*(BDFVector[1]*un + BDFVector[2]*u_(n-1))
        //note that the term related to BDFVector[0] is included in the LHS
        array_1d<double,LocalSize> aaa = ZeroVector(LocalSize);
        for(unsigned int k = 0; k<TNumNodes; k++)
        {
            unsigned int base=k*(TDim+1);
            for(unsigned int step=1; step<BDFVector.size(); step++)
            {
                const array_1d<double,3>& u = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY,step);
                const double& bdf_coeff = BDFVector[step];
                aaa[base]   += bdf_coeff*u[0];
                aaa[base+1] += bdf_coeff*u[1];
                aaa[base+2] += bdf_coeff*u[2];
            }
        }
        noalias(rRightHandSideVector) -= prod(MassMatrix,aaa);

        //****************************************************
        //finalize computation of the residual
        // Now calculate an additional contribution to the residual: r -= rLeftHandSideMatrix * (u,p)
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
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, U);
//                KRATOS_WATCH("line 517");

        //****************************************************
        //finalize computation of enrichment terms
        //(do static condensation) of enrichment terms
        //note that it each step we assume that the enrichment starts from 0
        if (ndivisions > 1)
        {
                //finalize the computation of the rhs
                noalias(enriched_rhs) -= prod(enrichment_terms_horizontal,U);

                double max_diag = 0.0;
                for(unsigned int k=0; k<TDim+1; k++)
                    if(fabs(enrichment_diagonal(k,k) ) > max_diag) max_diag = fabs(enrichment_diagonal(k,k) );
                if(max_diag == 0) max_diag = 1.0;



                if(positive_volume/Area < min_area_ratio)
                {
//                     KRATOS_WATCH(this->Id());
//                     KRATOS_WATCH(positive_volume/Area)
                    for(unsigned int i=0; i<TDim+1; i++)
                    {
                        if(distances[i] >= 0.0)
                        {
/*
                            for(unsigned int k=0; k<TDim+1; k++)
                            {
                                enrichment_diagonal(i,k) = 0.0;
                                enrichment_diagonal(k,i) = 0.0;
                            }*/
                            enrichment_diagonal(i,i) += 1000.0*max_diag;
//                             enriched_rhs[i] = 0.0;
//
//                             for(unsigned int k=0; k<enrichment_terms_horizontal.size2();k++)
//                             {
//                                 enrichment_terms_horizontal(i,k) = 0.0;
//                                 enrichment_terms_vertical(k,i) = 0.0;
//                             }
                        }
                    }
                }
                 if(negative_volume/Area < min_area_ratio)
                {
//                     KRATOS_WATCH(this->Id());
//                     KRATOS_WATCH(negative_volume/Area)
                    for(unsigned int i=0; i<TDim+1; i++)
                    {
                        if(distances[i] < 0.0)
                        {
                            enrichment_diagonal(i,i) += 1000.0*max_diag;
//                             for(unsigned int k=0; k<TDim+1; k++)
//                             {
//                                 enrichment_diagonal(i,k) = 0.0;
//                                 enrichment_diagonal(k,i) = 0.0;
//                             }
//                             enrichment_diagonal(i,i) = 1.0; //max_diag;
//                             enriched_rhs[i] = 0.0;
//
//
//                             for(unsigned int k=0; k<enrichment_terms_horizontal.size2();k++)
//                             {
//                                 enrichment_terms_horizontal(i,k) = 0.0;
//                                 enrichment_terms_vertical(k,i) = 0.0;
//                             }
                        }
                    }
                }


                //ensure the matrix is invertible
//                 for(unsigned int i=0; i<nenrichments; i++)
//                 {
//                     if(fabs(enrichment_diagonal(i,i)) < 1e-30)
//                     {
//                         enrichment_diagonal(i,i) = 1.0;
//                         enriched_rhs[i] = 0.0;
//                     }
//                 }

                //"weakly" impose continuity
//                 KRATOS_WATCH("line 541");
                for(unsigned int i=0; i<TDim; i++)
                {
                    const double di = fabs(distances[i]);

                    for(unsigned int j=i+1; j<TDim+1; j++)
                    {
                        const double dj =  fabs(distances[j]);

                        if( distances[i]*distances[j] < 0.0) //cut edge
                        {
                            double sum_d = di+dj;
                            double Ni = dj/sum_d;
                            double Nj = di/sum_d;

                            double penalty_coeff = max_diag*0.001; // h/BDFVector[0];
                            enrichment_diagonal(i,i) += penalty_coeff * Ni*Ni;
                            enrichment_diagonal(i,j) -= penalty_coeff * Ni*Nj;
                            enrichment_diagonal(j,i) -= penalty_coeff * Nj*Ni;
                            enrichment_diagonal(j,j) += penalty_coeff * Nj*Nj;

                        }
                    }
                }
//                 KRATOS_WATCH("line 565");


//                 KRATOS_WATCH(enrichment_diagonal);


                //add to LHS enrichment contributions
                Matrix inverse_diag(nenrichments, nenrichments);
                double det;
                MathUtils<double>::InvertMatrix(enrichment_diagonal,inverse_diag,det);

                  //  double inverse_diag_term = 1.0 / ( enrichment_diagonal);
//        KRATOS_WATCH(this->Id());
//        KRATOS_WATCH(enrichment_terms_horizontal);
//        KRATOS_WATCH(enrichment_terms_vertical);
//        KRATOS_WATCH(inverse_diag);
                Matrix tmp = prod(inverse_diag,enrichment_terms_horizontal);
                noalias(rLeftHandSideMatrix) -= prod(enrichment_terms_vertical,tmp);

                Vector tmp2 = prod(inverse_diag,enriched_rhs);
                noalias(rRightHandSideVector) -= prod(enrichment_terms_vertical,tmp2);

/*                for (unsigned int i = 0; i < LocalSize; i++)
                    for (unsigned int j = 0; j < LocalSize; j++)
                        rLeftHandSideMatrix(i, j) -= inverse_diag_term * enrichment_terms_vertical[i] * enrichment_terms_horizontal[j];
                noalias(rRightHandSideVector) -= (inverse_diag_term*enriched_rhs )*enrichment_terms_vertical;
      */      }

//       KRATOS_WATCH("finished elem")
    }

    /** does nothing for this element
     * @param rVariable
     * @param Output
     * @param rCurrentProcessInfo
     */
    void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                           array_1d<double, 3 > & rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override
    {

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
        if (PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "PRESSURE Key is 0. Check if the application was correctly registered.", "");
        if (DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY Key is 0. Check if the application was correctly registered.", "");
        if (VISCOSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "VISCOSITY Key is 0. Check if the application was correctly registered.", "");
        if (DYNAMIC_TAU.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "DYNAMIC_TAU Key is 0. Check if the application was correctly registered.", "");
        if (DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "DELTA_TIME Key is 0. Check if the application was correctly registered.", "");
//         if (ADVPROJ.Key() == 0)
//             KRATOS_THROW_ERROR(std::invalid_argument, "ADVPROJ Key is 0. Check if the application was correctly registered.", "");
//         if (DIVPROJ.Key() == 0)
//             KRATOS_THROW_ERROR(std::invalid_argument, "DIVPROJ Key is 0. Check if the application was correctly registered.", "");
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
            if (this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                    this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                    this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY component degree of freedom on node ", this->GetGeometry()[i].Id());
            if (this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE component degree of freedom on node ", this->GetGeometry()[i].Id());
        }

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
        buffer << "TwoFluidVMS #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TwoFluidVMS" << TDim << "D";
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
            ElementBaseType::EvaluateInPoint(value,rVariable,rShapeFunc);
        rResult = value;
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


        /// Add a the contribution from a single integration point to the velocity contribution
    void AddIntegrationPointVelocityContribution(MatrixType& rDampingMatrix,
            VectorType& rDampRHS,
            const double Density,
            const double Viscosity,
            const array_1d< double, 3 > & rAdvVel,
            const double ReactionTerm,
            const double TauOne,
            const double TauTwo,
            const array_1d< double, TNumNodes >& rShapeFunc,
            const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
            const double Weight)
    {
        const unsigned int BlockSize = TDim + 1;

        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        // Build the local matrix and RHS
        unsigned int FirstRow(0), FirstCol(0); // position of the first term of the local matrix that corresponds to each node combination
        double K, PDivV, L, qF; // Temporary results

        array_1d<double,3> BodyForce = ZeroVector(3);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,rShapeFunc);
        BodyForce *= Density;

        array_1d<double, TNumNodes> StabilizationOperator = Density*AGradN;
        noalias(StabilizationOperator) -= ReactionTerm * rShapeFunc;
        StabilizationOperator *= TauOne;

        for (unsigned int i = 0; i < TNumNodes; ++i) // iterate over rows
        {
            for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over columns
            {
                // Convection + Reaction: v *( a*Grad(u) + sigma*u )
                // For Darcy: sigma = A + B|u|
                K = rShapeFunc[i] * ( Density*AGradN[j] + ReactionTerm*rShapeFunc[j] );
                // Stabilization: (a * Grad(v) - sigma * N) * TauOne *( a*Grad(u) + sigma*u )
                K += StabilizationOperator[i] * ( Density*AGradN[j] + ReactionTerm*rShapeFunc[j] );
                K *= Weight;

                // q-p stabilization block (reset result)
                L = 0;

                const array_1d<double,3>& OldVel = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);

                for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
                {
                    // v-p block (pressure gradient)
                    double div_v_p = rShapeDeriv(i, m) * rShapeFunc[j];
                    double stab_grad_p = StabilizationOperator[i] * rShapeDeriv(j,m);
                    rDampingMatrix(FirstRow + m, FirstCol + TDim) += Weight * (stab_grad_p - div_v_p);

                    // q-u block (velocity divergence)
                    double q_div_u = rShapeFunc[i] * rShapeDeriv(j,m);
                    double stab_div_u = TauOne*rShapeDeriv(i,m)* ( Density*AGradN[j] + ReactionTerm * rShapeFunc[j] );
                    rDampingMatrix(FirstRow + TDim, FirstCol + m) += Weight * ( q_div_u + stab_div_u );

                    PDivV = rShapeDeriv(i, m) * rShapeFunc[j]; // Div(v) * p
                    rDampRHS[FirstCol + TDim] -=  Weight * PDivV*OldVel[m];

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

                // Update reference column index for next iteration
                FirstCol += BlockSize;
            }

            // Operate on RHS
            qF = 0.0;
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rDampRHS[FirstRow + d] += Weight * StabilizationOperator[i] * BodyForce[d]; // (a * Grad(v) - sigma * N) * TauOne * (Density * BodyForce)
                qF += rShapeDeriv(i, d) * BodyForce[d];
            }
            rDampRHS[FirstRow + TDim] += Weight * TauOne * qF; // Grad(q) * TauOne * (Density * BodyForce)

            // Update reference indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }

        this->AddViscousTerm(rDampingMatrix,rShapeDeriv,Viscosity*Weight);
    }

    void AddMassStabTerms(MatrixType& rLHSMatrix,
                          const double Density,
                          const array_1d<double, 3 > & rAdvVel,
                          const double ReactionTerm,
                          const double TauOne,
                          const array_1d<double, TNumNodes>& rShapeFunc,
                          const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                          const double Weight)
    {
        const unsigned int BlockSize = TDim + 1;

        unsigned int FirstRow(0), FirstCol(0);
        double K; // Temporary results

        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        array_1d<double, TNumNodes> StabilizationOperator = Density*AGradN;
        noalias(StabilizationOperator) -= ReactionTerm * rShapeFunc;
        StabilizationOperator *= TauOne;

        // Note: Dof order is (vx,vy,[vz,]p) for each node
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            // Loop over columns
            for (unsigned int j = 0; j < TNumNodes; ++j)
            {
                // Delta(u) * TauOne * [ AdvVel * Grad(v) - sigma * N ] in velocity block
                K = Weight * StabilizationOperator[i] * Density * rShapeFunc[j];

                for (unsigned int d = 0; d < TDim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                {
                    rLHSMatrix(FirstRow + d, FirstCol + d) += K;
                    // Delta(u) * TauOne * Grad(q) in q * Div(u) block
                    rLHSMatrix(FirstRow + TDim, FirstCol + d) += Weight * TauOne * rShapeDeriv(i, d) * Density * rShapeFunc[j];
                }
                // Update column index
                FirstCol += BlockSize;
            }
            // Update matrix indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }
    }

    void CalculateStabilizationTau(
        double& TauOne,
        double& TauTwo,
        const double VelNorm,
        const double ElemSize,
        const double Density,
        const double DynamicViscosity,
        const double ReactionTerm,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const double DynamicTerm = rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME];
        double InvTau = Density * ( DynamicTerm + 2.0*VelNorm / ElemSize ) + 4.0*DynamicViscosity/ (ElemSize * ElemSize) + ReactionTerm;
        TauOne = 1.0 / InvTau;

        TauTwo = DynamicViscosity + 0.5 * Density * ElemSize * VelNorm;
    }

    virtual double CalculateDarcyTerm(
        const double Density,
        const double DynamicViscosity,
        const double LinearCoefficient,
        const double NonlinearCoefficient,
        const array_1d<double, TNumNodes>& rShapefunctions) {

        array_1d<double,3> velocity;
        this->GetAdvectiveVel(velocity, rShapefunctions);
        const double velocity_norm = MathUtils<double>::Norm3(velocity);

        return DynamicViscosity * LinearCoefficient + Density * NonlinearCoefficient*velocity_norm;
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
    TwoFluidVMS & operator=(TwoFluidVMS const& rOther);
    /// Copy constructor.
    TwoFluidVMS(TwoFluidVMS const& rOther);
    ///@}
}; // Class TwoFluidVMS
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
                                  TwoFluidVMS<TDim, TNumNodes>& rThis)
{
    return rIStream;
}
/// output stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const TwoFluidVMS<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}
///@} // Fluid Dynamics Application group
} // namespace Kratos.
#endif // KRATOS_TWO_FLUID_VMS_H_INCLUDED  defined
