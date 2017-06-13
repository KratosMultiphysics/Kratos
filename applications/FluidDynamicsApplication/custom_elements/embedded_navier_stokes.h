//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_EMBEDDED_NAVIER_STOKES)
#define  KRATOS_EMBEDDED_NAVIER_STOKES

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"
#include "custom_elements/navier_stokes.h"
#include "utilities/split_tetrahedra_utilities.h"

namespace Kratos
{

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

/**
* This is a 2D and 3D Navier-Stokes embedded element, stabilized by employing an ASGS stabilization
* Both the formulation and the symbolic implementation can be found in the symbolic_generation
* folder of the FluidDynamicsApplication.
*/
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
class EmbeddedNavierStokes : public NavierStokes<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedNavierStokes);

    typedef NavierStokes<TDim, TNumNodes>                              BaseType;

    typedef typename BaseType::ElementDataStruct                ElementDataType;

    typedef typename BaseType::VectorType                            VectorType;

    typedef typename BaseType::MatrixType                            MatrixType;

    typedef typename BaseType::IndexType                              IndexType;

    typedef typename BaseType::GeometryType::Pointer        GeometryPointerType;

    typedef typename BaseType::NodesArrayType                    NodesArrayType;

    typedef typename BaseType::PropertiesType::Pointer    PropertiesPointerType;

    struct ElementSplittingDataStruct
    {

        MatrixType      N_container;                        // Container with the shape functions values in each partition Gauss points
        MatrixType      N_container_cut;                    // Container with the shape functions values evaluated at the intersection (interface) Gauss points
        MatrixType      N_container_int;                    // Containter with the interior (structure) shape functions values in each intersection (interface) Gauss point
        MatrixType      N_container_out;                    // Containter with the outside (fluid) shape functions values in each intersection (interface) Gauss point

        VectorType      gauss_volumes;                      // Container with the Gauss points volumes in each partition
        VectorType      partition_signs;                    // Indicates if the edge is cut or not
        VectorType      edge_areas;                         // Vector containing the edge intersection surface area associated to each intersection point in the cut edges
        VectorType      cut_edge_areas;                     // Similar to edge_areas but only containing the intersected edges areas

        std::vector<unsigned int>   int_vec_identifiers;    // Interior (structure) nodes identifiers
        std::vector<unsigned int>   out_vec_identifiers;    // Outside (fluid) nodes identifiers

        array_1d<double, TDim>      intersection_normal;    // Intersection plane unit normal vector

        unsigned int    ndivisions;                         // Number of element subdivisions
        unsigned int    ncutpoints;                         // Number of intersected edges
        unsigned int    n_pos;                              // Number of postivie distance nodes
        unsigned int    n_neg;                              // Number of negative distance nodes

    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    EmbeddedNavierStokes(IndexType NewId, GeometryPointerType pGeometry)
        : NavierStokes<TDim, TNumNodes>(NewId, pGeometry)
    {}

    EmbeddedNavierStokes(IndexType NewId, GeometryPointerType pGeometry, PropertiesPointerType pProperties)
        : NavierStokes<TDim, TNumNodes>(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~EmbeddedNavierStokes() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, Element::PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< EmbeddedNavierStokes < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }


    Element::Pointer Create(IndexType NewId, Element::GeometryType::Pointer pGeom, Element::PropertiesType::Pointer pProperties) const override
    {
        return boost::make_shared< EmbeddedNavierStokes < TDim, TNumNodes > >(NewId, pGeom, pProperties);
    }


    /**
     * Clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param rThisNodes: the nodes of the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Element::Pointer pNewElement = Create(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

        pNewElement->SetData(this->GetData());
        pNewElement->SetFlags(this->GetFlags());

        return pNewElement;
    }


    /**
     * Calculates both LHS and RHS contributions
     * @param rLeftHandSideMatrix: reference to the LHS matrix
     * @param rRightHandSideVector: reference to the RHS vector
     * @param rCurrentProcessInfo: reference to the ProcessInfo
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
        {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); // false says to not preserve existing storage!!
        }
        else if (rLeftHandSideMatrix.size2() != MatrixSize)
        {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); // false says to not preserve existing storage!!
        }

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false);            // false says to not preserve existing storage!!

        // Struct to pass around the data
        ElementDataType data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Getting the nodal distance values
        array_1d<double, TNumNodes> distances;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Number of positive and negative distance function values
        unsigned int npos=0, nneg=0;
        for (unsigned int i = 0; i<TNumNodes; i++)
        {
            if(distances[i] > 0.0)
                npos++;
            else
                nneg++;
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize, MatrixSize);   // LHS initialization
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);              // RHS initialization

        // Element LHS and RHS contributions computation
        if(npos == TNumNodes) // All nodes belong to fluid domain
        {
            ComputeElementAsFluid<MatrixSize>(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo);
        }
        else if(nneg == TNumNodes) // All nodes belong to structure domain
        {
            rLeftHandSideMatrix.clear();
            rRightHandSideVector.clear();
        }
        else // Element intersects both fluid and structure domains
        {
            ComputeElementAsMixed<MatrixSize>(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo, distances);
        }

        KRATOS_CATCH("Error in embedded Navier-Stokes element")
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
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        // Base element check
        int ErrorCode = NavierStokes<TDim, TNumNodes>::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        // Specific embedded element check
        if(DISTANCE.Key() == 0)
            KRATOS_ERROR << "DISTANCE Key is 0. Check if the application was correctly registered.";

        return 0;

        KRATOS_CATCH("");
    }


    /**
     * Calculates both LHS and RHS elemental contributions for those cases in where
     * all the nodes belong to the fluid domain.
     * @param lhs_local: reference to the Gauss pt. LHS matrix
     * @param rhs_local: reference to the Gauss pt. RHS vector
     * @param rLeftHandSideMatrix: reference to the LHS matrix
     * @param rRightHandSideVector: reference to the RHS vector
     * @param rData: reference to element data structure
     * @param rCurrentProcessInfo: reference to the ProcessInfo
     */
    template<unsigned int MatrixSize>
    void ComputeElementAsFluid(MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               ElementDataType& rData,
                               ProcessInfo& rCurrentProcessInfo)
    {
        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_local;

        // Shape functions Gauss points values
        bounded_matrix<double, TNumNodes, TNumNodes> Ncontainer; // Container with the evaluation of the 4 shape functions in the 4 Gauss pts.
        BaseType::GetShapeFunctionsOnGauss(Ncontainer);

        // Loop on gauss point
        for(unsigned int igauss = 0; igauss<Ncontainer.size2(); igauss++)
        {
            noalias(rData.N) = row(Ncontainer, igauss);

            BaseType::ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            BaseType::ComputeGaussPointLHSContribution(lhs_local, rData);
            BaseType::ComputeGaussPointRHSContribution(rhs_local, rData);

            // All the Gauss pts. have the same weight so the accumulated contributions can be multiplied by volume/n_nodes at the end
            noalias(rLeftHandSideMatrix) += lhs_local;
            noalias(rRightHandSideVector) += rhs_local;
        }

        rLeftHandSideMatrix *= rData.volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= rData.volume/static_cast<double>(TNumNodes);;
    }


    /**
    * Calculates both LHS and RHS elemental contributions for those cases in where
    * the element has both fluid and structure nodes.
    * @param lhs_local: reference to the Gauss pt. LHS matrix
    * @param rhs_local: reference to the Gauss pt. RHS vector
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rCurrentProcessInfo: reference to the ProcessInfo
    * @param distances: reference to an array containing the nodal distance
    */
    template<unsigned int MatrixSize>
    void ComputeElementAsMixed(MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               ElementDataType& rData,
                               ProcessInfo& rCurrentProcessInfo,
                               array_1d<double, TNumNodes>& rDistances)
    {
        constexpr unsigned int nEdges = (TDim-1)*3;     // Edges per element

        MatrixType Ncontainer;
        VectorType gauss_volumes;
        VectorType signs(nEdges); //ATTENTION: this shall be initialized of size 6 in 3D and 3 in 2D
        VectorType edge_areas(nEdges);

        // Creation of the data structure to store the intersection information
        ElementSplittingDataStruct SplittingData;

        // Splitting to determine the new Gauss pts.
        this->ComputeSplitting(rData, rDistances, SplittingData);

        if(SplittingData.ndivisions == 1)
        {
            // Cases exist when the element is like not subdivided due to the characteristics of the provided distance.
            // Use the distance value at the central Gauss pt. to determine if the element is wether fluid or structure.
            array_1d<double, TNumNodes> Ncenter;
            for(unsigned int i=0; i<TNumNodes; i++) Ncenter[i]=0.25;
            const double dgauss = inner_prod(rDistances, Ncenter);

            // Gauss pt. is FLUID (the element is close to be full of fluid)
            if(dgauss > 0.0)
            {
                ComputeElementAsFluid<MatrixSize>(rLeftHandSideMatrix, rRightHandSideVector, rData, rCurrentProcessInfo);
            }
        }
        else
        {
            // Allocate memory needed
            array_1d<double, MatrixSize> rhs_local;
            bounded_matrix<double,MatrixSize, MatrixSize> lhs_local;

            // Loop over subdivisions
            for(unsigned int division = 0; division<SplittingData.ndivisions; division++)
            {
                // Loop on gauss points
                for(unsigned int igauss = 0; igauss<4; igauss++)
                {
                    noalias(rData.N) = row(SplittingData.N_container, division*4+igauss); // Take the new Gauss pts. shape functions values

                    const double dgauss = inner_prod(rDistances, rData.N);  // Compute the distance on the gauss point

                    // Gauss pt. is FLUID
                    if(dgauss > 0.0)
                    {
                        BaseType::ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

                        BaseType::ComputeGaussPointLHSContribution(lhs_local, rData);
                        BaseType::ComputeGaussPointRHSContribution(rhs_local, rData);

                        noalias(rLeftHandSideMatrix) += SplittingData.gauss_volumes[division*4+igauss]*lhs_local;
                        noalias(rRightHandSideVector) += SplittingData.gauss_volumes[division*4+igauss]*rhs_local;
                    }
                }
            }
        }

        // Add level set boundary terms, penalty and modified Nitche contributions
        AddBoundaryConditionElementContribution(rLeftHandSideMatrix, rRightHandSideVector, rDistances, rData, SplittingData);
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

    virtual std::string Info() const override
    {
        return "EmbeddedNavierStokes3D #";
    }

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

    EmbeddedNavierStokes() : NavierStokes<TDim, TNumNodes>()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

    /**
    * Computes the element splitting according to the distance values.
    * @param rData: reference to element data structure
    * @param rDistances: reference to an array containing the nodal distance values
    * @param rSplittingData: reference to the intersection data structure
    */
    void ComputeSplitting(const ElementDataType& rData,
                          const array_1d<double,TNumNodes>& rDistances,
                          ElementSplittingDataStruct& rSplittingData)
    {

        rSplittingData.edge_areas = ZeroVector((TDim-1)*3);         // Should be initialized to the edge number
        rSplittingData.partition_signs = ZeroVector((TDim-1)*3);    // Should be initialized to the edge number

        VectorType NodalDistances = rDistances;
        MatrixType NoSplitGradients = rData.DN_DX;

        // Fill nodal coordinates
        MatrixType NodalCoords(TNumNodes, TDim);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
            for (unsigned int j = 0; j < TDim; j++)
            {
                NodalCoords(i, j) = xyz[j];
            }
        }

        rSplittingData.ndivisions = SplitTetrahedraUtilities::CalculateSplitTetrahedraShapeFuncions(NodalCoords, NoSplitGradients, NodalDistances,
                                                                                                    rSplittingData.gauss_volumes,
                                                                                                    rSplittingData.N_container,
                                                                                                    rSplittingData.partition_signs,
                                                                                                    rSplittingData.edge_areas);

        if (rSplittingData.ndivisions > 1)
        {
            this->GetIntersectionGeometryData(rDistances, rData, rSplittingData);
        }
    }

    /**
    * If the element is split, this function computes and stores in a data structure
    * the intersection geometry information.
    * @param rDistances: reference to an array containing the nodal distance values
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void GetIntersectionGeometryData(const array_1d<double,TNumNodes>& rDistances,
                                     const ElementDataType& rData,
                                     ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int nDOFs = TDim+1;          // DOFs per node
        constexpr unsigned int nEdges = (TDim-1)*3;     // Edges per element

        rSplittingData.n_pos = 0;
        rSplittingData.n_neg = 0;

        (rSplittingData.int_vec_identifiers).clear();
        (rSplittingData.out_vec_identifiers).clear();

        for (unsigned int inode = 0; inode<TNumNodes; inode++)
        {
            if(rDistances[inode] > 0)
            {
                (rSplittingData.int_vec_identifiers).push_back(inode);
                rSplittingData.n_pos++; // Number of positive distance value nodes (fluid)
            }
            else
            {
                (rSplittingData.out_vec_identifiers).push_back(inode);
                rSplittingData.n_neg++; // Number of negative distance value nodes (structure)
            }
        }

        // Identify the cut edges
        rSplittingData.ncutpoints = 0;                  // Number of cut edges (or cut points)
        array_1d<unsigned int, nEdges> i_edges;         // i-node of the edge
        array_1d<unsigned int, nEdges> j_edges;         // j-node of the edge
        array_1d<unsigned int, nEdges> cut_edges;       // 0: non cut edge 1: cut edge

        unsigned int aux_count = 0;
        for (unsigned int i = 0; i<TDim; i++)
        {
            for(unsigned int j = i+1; j<nDOFs; j++)
            {
                i_edges(aux_count) = i;     // Construct the edges i-nodes vector
                j_edges(aux_count) = j;     // Construct the edges j-nodes vector
                cut_edges(aux_count) = 0;

                if((rDistances[i]*rDistances[j])<0)
                {
                    rSplittingData.ncutpoints++;
                    cut_edges(aux_count) = 1; // Flag 1 says that this edge is cut
                }

                aux_count++;
            }
        }

        // Store the non-zero edge areas in an auxiliar vector cut_areas
        (rSplittingData.cut_edge_areas).resize(rSplittingData.ncutpoints, false);

        // Compute the shape function values at the intersection points
        rSplittingData.N_container_cut = ZeroMatrix(rSplittingData.ncutpoints, TNumNodes);                     // Shape function values at the interface (intersection)
        rSplittingData.N_container_int = ZeroMatrix(rSplittingData.ncutpoints, rSplittingData.n_pos);          // Interior nodes shape function values at the interface
        rSplittingData.N_container_out = ZeroMatrix(rSplittingData.ncutpoints, rSplittingData.n_neg);          // Exterior nodes shape function values at the interface

        // rNcontainer_cut matrix fill
        unsigned int icut = 0;
        for (unsigned int i = 0; i<nEdges; i++)
        {
            if (cut_edges(i) == 1)
            {
                unsigned int i_node = i_edges(i);
                unsigned int j_node = j_edges(i);

                // We take advantage of the fact that the sh. funct. that do not belong to the edge vanish along it (linear elements)
                double aux = fabs(rDistances[i_node])/(fabs(rDistances[i_node])+fabs(rDistances[j_node])+1e-30);

                // The columns associated to the nodes that conform the cut edge are filled. The other two columns remain as zero
                rSplittingData.N_container_cut(icut,i_node) = 1.0-aux;
                rSplittingData.N_container_cut(icut,j_node) = aux;

                rSplittingData.cut_edge_areas(icut) = rSplittingData.edge_areas(i);

                icut++;
            }
        }

        // Ncontainer_int and Ncontainer_out matrices fill
        for (unsigned int i = 0; i<rSplittingData.ncutpoints; i++)
        {
            for (unsigned int j = 0; j<rSplittingData.n_pos; j++)
            {
                rSplittingData.N_container_int(i,j) = rSplittingData.N_container_cut(i, rSplittingData.int_vec_identifiers[j]);
            }

            for (unsigned int j = 0; j<rSplittingData.n_neg; j++)
            {
                rSplittingData.N_container_out(i,j) = rSplittingData.N_container_cut(i, rSplittingData.out_vec_identifiers[j]);
            }
        }

        // Compute intersection normal (gradient of the distance function)
        noalias(rSplittingData.intersection_normal) = -prod(trans(rData.DN_DX), rDistances);
        rSplittingData.intersection_normal /= (norm_2(rSplittingData.intersection_normal)+1e-15);

    }


    /**
    * This functions adds the contribution of the boundary terms in the level set cut
    * These terms, which do not vanish at the level set since the test function is not zero
    * at the intersection points, come from the integration by parts of the stress term.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddIntersectionBoundaryTermsContribution(MatrixType& rLeftHandSideMatrix,
                                                  VectorType& rRightHandSideVector,
                                                  const ElementDataType& rData,
                                                  const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
            prev_sol(i*BlockSize+TDim) = rData.p(i);
        }

        // Declare auxiliar arrays
        bounded_matrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Get the normal projection matrix in Voigt notation
        bounded_matrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix = ZeroMatrix(TDim, (TDim-1)*3);
        SetVoigtNormalProjectionMatrix(rSplittingData, voigt_normal_projection_matrix);

        // Get the strain matrix (constant since linear elements are used)
        bounded_matrix<double, (TDim-1)*3, TNumNodes*TDim> B_matrix = ZeroMatrix((TDim-1)*3, TNumNodes*TDim);
        SetStrainMatrix(rData, B_matrix);

        // Expand the B matrix to set 0 in the pressure rows.
        bounded_matrix<double, (TDim-1)*3, MatrixSize> B_matrix_exp = ZeroMatrix((TDim-1)*3, MatrixSize);
        for (unsigned int i=0; i<TNumNodes; ++i)
        {
            for (unsigned int j=0; j<TDim; ++j)
            {
                for (unsigned int k=0; k<(TDim-1)*3; ++k)
                {
                    B_matrix_exp(k, i*BlockSize+j) = B_matrix(k, i*TDim+j);
                }
            }
        }

        // Compute some element constant matrices
        const bounded_matrix<double, TDim, (TDim-1)*3> aux_matrix_AC = prod(voigt_normal_projection_matrix, rData.C);
        const bounded_matrix<double, (TDim-1)*3, MatrixSize> aux_matrix_ACB = prod(aux_matrix_AC, B_matrix_exp);

        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

            // Fill the pressure to Voigt notation operator matrix
            bounded_matrix<double, (TDim-1)*3, MatrixSize> pres_to_voigt_matrix_op = ZeroMatrix((TDim-1)*3, MatrixSize);
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    pres_to_voigt_matrix_op(comp, i*BlockSize+TDim) = aux_cut(i);
                }
            }

            // Set the shape functions auxiliar transpose matrix
            bounded_matrix<double, MatrixSize, TDim> N_aux_trans = ZeroMatrix(MatrixSize, TDim);
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    N_aux_trans(i*BlockSize+comp, comp) = aux_cut(i);
                }
            }

            // Contribution coming fron the shear stress operator
            auxLeftHandSideMatrix += weight*prod(N_aux_trans, aux_matrix_ACB);

            // Contribution coming from the pressure terms
            const bounded_matrix<double, MatrixSize, (TDim-1)*3> N_voigt_proj_matrix= prod(N_aux_trans, voigt_normal_projection_matrix);
            auxLeftHandSideMatrix -= weight*prod(N_voigt_proj_matrix, pres_to_voigt_matrix_op);

        }

        // LHS assembly
        rLeftHandSideMatrix -= auxLeftHandSideMatrix;

        // RHS assembly
        rRightHandSideVector += prod(auxLeftHandSideMatrix, prev_sol);

    }


    /**
    * This function computes the penalty coefficient for the level set BC imposition
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    double ComputePenaltyCoefficient(MatrixType& rLeftHandSideMatrix,
                                     const ElementDataType& rData,
                                     const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Get the intersection area
        double intersection_area = 0.0;
        for (unsigned int icut = 0; icut<rSplittingData.ncutpoints; icut++)
        {
            intersection_area += rSplittingData.cut_edge_areas(icut);
        }

        // Compute the penalty coefficient as K*max(LHS(i,i))*IntArea (we integrate P_gamma over the intersection area)
        double diag_max = 0.0;
        for (unsigned int i=0; i<MatrixSize; i++)
        {
            if ((fabs(rLeftHandSideMatrix(i,i)) > diag_max) && (i%BlockSize != 0.0))
            {
                diag_max = fabs(rLeftHandSideMatrix(i,i)); // Maximum diagonal value (associated to velocity)
            }
        }

        const double K = 100.0;
        const double denominator = std::max(0.0001*rData.h*rData.h, intersection_area);
        const double pen_coef = K*diag_max/denominator;

        return pen_coef;
    }


    /**
    * This functions adds the penalty extra term level set contribution.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddBoundaryConditionPenaltyContribution(MatrixType& rLeftHandSideMatrix,
                                                 VectorType& rRightHandSideVector,
                                                 const ElementDataType& rData,
                                                 const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
        }

        MatrixType P_gamma(TNumNodes, TNumNodes);    // Penalty matrix
        noalias(P_gamma) = ZeroMatrix(TNumNodes, TNumNodes);

        for (unsigned int icut = 0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

            P_gamma += weight*outer_prod(aux_cut,aux_cut);
        }

        // Multiply the penalty matrix by the penalty coefficient
        double pen_coef = ComputePenaltyCoefficient(rLeftHandSideMatrix, rData, rSplittingData);
        P_gamma *= pen_coef;

        VectorType auxRightHandSideVector = ZeroVector(MatrixSize);
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // LHS penalty contribution assembly (symmetric mass matrix)
        for (unsigned int i = 0; i<TNumNodes; i++)
        {
            // Diagonal terms
            for (unsigned int comp = 0; comp<TDim; comp++)
            {
                auxLeftHandSideMatrix(i*BlockSize+comp, i*BlockSize+comp) = P_gamma(i,i);
            }
            // Off-diagonal terms
            for (unsigned int j = i+1; j<TNumNodes; j++)
            {
                for (unsigned int comp = 0; comp<TDim; comp++)
                {
                    auxLeftHandSideMatrix(i*BlockSize+comp, j*BlockSize+comp) = P_gamma(i,j);
                    auxLeftHandSideMatrix(j*BlockSize+comp, i*BlockSize+comp) = P_gamma(i,j);
                }
            }
        }

        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS penalty contribution assembly
        if (this->Has(EMBEDDED_VELOCITY))
        {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; i++)
            {
                for (unsigned int comp=0; comp<TDim; comp++)
                {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            rRightHandSideVector += prod(auxLeftHandSideMatrix, aux_embedded_vel);
        }

        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol); // Residual contribution assembly

    }


    /**
    * This functions adds the slip no penetration penalty extra term level set contribution.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddSlipBoundaryConditionPenaltyContribution(MatrixType& rLeftHandSideMatrix,
                                                     VectorType& rRightHandSideVector,
                                                     const ElementDataType& rData,
                                                     const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
        }

        // Compute the consistency coefficient
        const double h = rData.h; // Characteristic element size
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);
        const double cons_coef = eff_mu/h;

        // Declare the auxiliar matrices used in the assembly
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);
        MatrixType auxLeftHandSideMatrixNoPenetration = ZeroMatrix(MatrixSize, MatrixSize);
        MatrixType auxLeftHandSideMatrixNoShearStress = ZeroMatrix(MatrixSize, MatrixSize);

        // Compute the normal projection matrix
        bounded_matrix<double, TDim, TDim> normal_projection_matrix;
        SetNormalProjectionMatrix(rSplittingData, normal_projection_matrix);

        // Compute the Voigt notation tangential projection matrix
        bounded_matrix<double, TDim, (TDim-1)*3> voigt_tang_proj_matrix;        // Set the Voigt notation tangential projection matrix
        SetVoigtTangentialProjectionMatrix2(rSplittingData, voigt_tang_proj_matrix);

        // Compute the constant contribution in the null tangential stress slip requirement
        // const bounded_matrix<double, TDim, (TDim-1)*3> tensor_C_proj = prod(voigt_tang_proj_matrix, rData.C);
        // bounded_matrix<double, (TDim-1)*3, TNumNodes*TDim> B_matrix;
        // SetStrainMatrix(rData, B_matrix);
        // const bounded_matrix<double, TDim, TNumNodes*TDim> aux_stress_matrix = prod(tensor_C_proj, B_matrix);

        // Compute the penalty LHS contribution
        for (unsigned int icut = 0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

            MatrixType shfunc_matrix_vel = ZeroMatrix(TNumNodes*TDim,TDim);
            // Fill the velocity shape functions matrix
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int j=0; j<TDim; ++j)
                {
                    shfunc_matrix_vel(i*TDim+j,j) = aux_cut(i);
                }
            }

            // No penetration condition penalty contribution computation
            const bounded_matrix<double, TNumNodes*TDim, TDim> aux_matrix_1 = prod(shfunc_matrix_vel, normal_projection_matrix);
            const bounded_matrix<double, TNumNodes*TDim, TNumNodes*TDim> aux_matrix_2 = prod(aux_matrix_1, trans(shfunc_matrix_vel));

            // No tangential stres condition penalty contribution computation
            // const bounded_matrix<double, TNumNodes*TDim, TNumNodes*TDim> aux_matrix_3 = prod(shfunc_matrix_vel, aux_stress_matrix);

            // Assemble the previous contribution
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int j=0; j<TNumNodes; ++j)
                {
                    for (unsigned int ii=0; ii<TDim; ++ii)
                    {
                        for (unsigned int jj=0; jj<TDim; ++jj)
                        {
                            auxLeftHandSideMatrixNoPenetration(i*BlockSize+ii,j*BlockSize+jj) += cons_coef*weight*aux_matrix_2(i*TDim+ii,j*TDim+jj);
                            // auxLeftHandSideMatrixNoShearStress(i*BlockSize+ii,j*BlockSize+jj) += weight*aux_matrix_3(i*TDim+ii,j*TDim+jj);
                        }
                    }
                }
            }
        }

        // Add both LHS contributions
        auxLeftHandSideMatrix = auxLeftHandSideMatrixNoPenetration + auxLeftHandSideMatrixNoShearStress;

        // Multiply the LHS contribution by the penalty coefficient
        const double pen_coef = ComputePenaltyCoefficient(rLeftHandSideMatrix, rData, rSplittingData);
        auxLeftHandSideMatrix *= pen_coef;

        // Add the LHS penalty contribution
        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS penalty contribution assembly
        if (this->Has(EMBEDDED_VELOCITY))
        {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; i++)
            {
                for (unsigned int comp=0; comp<TDim; comp++)
                {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            rRightHandSideVector += prod(auxLeftHandSideMatrixNoPenetration, aux_embedded_vel);
        }

        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol); // Residual contribution assembly

    }


    /**
    * This functions adds the level set strong boundary condition imposition contribution.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddBoundaryConditionModifiedNitcheContribution(MatrixType& rLeftHandSideMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        const ElementDataType& rData,
                                                        const ElementSplittingDataStruct& rSplittingData)
    {

        constexpr unsigned int BlockSize = TDim+1;                 // Block size
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;   // Matrix size

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
        }

        // Compute the BCs imposition matrices
        MatrixType M_gamma = ZeroMatrix(rSplittingData.n_neg, rSplittingData.n_neg);      // Outside nodes matrix (Nitche contribution)
        MatrixType N_gamma = ZeroMatrix(rSplittingData.n_neg, rSplittingData.n_pos);      // Interior nodes matrix (Nitche contribution)
        MatrixType f_gamma = ZeroMatrix(rSplittingData.n_neg, TNumNodes); // Matrix to compute the RHS (Nitche contribution)

        VectorType aux_out(rSplittingData.n_neg);
        VectorType aux_int(rSplittingData.n_pos);
        VectorType aux_cut(TNumNodes);

        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            double weight = rSplittingData.cut_edge_areas(icut);

            aux_out = row(rSplittingData.N_container_out, icut);
            aux_int = row(rSplittingData.N_container_int, icut);
            aux_cut = row(rSplittingData.N_container_cut, icut);

            M_gamma += weight*outer_prod(aux_out,aux_out);
            N_gamma += weight*outer_prod(aux_out,aux_int);
            f_gamma += weight*outer_prod(aux_out,aux_cut);
        }

        // Declare auxLeftHandSideMatrix
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // LHS outside nodes contribution assembly
        // Outer nodes contribution assembly
        for (unsigned int i = 0; i<rSplittingData.n_neg; i++)
        {
            unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

            for (unsigned int j = 0; j<rSplittingData.n_neg; j++)
            {
                unsigned int out_node_col_id = rSplittingData.out_vec_identifiers[j];

                for (unsigned int comp = 0; comp<TDim; comp++)
                {
                    auxLeftHandSideMatrix(out_node_row_id*BlockSize+comp, out_node_col_id*BlockSize+comp) = M_gamma(i, j);
                }
            }
        }

        // Interior nodes contribution assembly
        for (unsigned int i = 0; i<rSplittingData.n_neg; i++)
        {
            unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

            for (unsigned int j = 0; j<rSplittingData.n_pos; j++)
            {
                unsigned int int_node_col_id = rSplittingData.int_vec_identifiers[j];

                for (unsigned int comp = 0; comp<TDim; comp++)
                {
                    auxLeftHandSideMatrix(out_node_row_id*BlockSize+comp, int_node_col_id*BlockSize+comp) = N_gamma(i, j);
                }
            }
        }

        // LHS outside Nitche contribution assembly
        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS outside Nitche contribution assembly
        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol);

        // Compute f_gamma if level set velocity is not 0
        if (this->Has(EMBEDDED_VELOCITY))
        {
            auxLeftHandSideMatrix.clear();

            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; i++)
            {
                aux_embedded_vel(i*BlockSize) = embedded_vel(0);
                aux_embedded_vel(i*BlockSize+1) = embedded_vel(1);
                aux_embedded_vel(i*BlockSize+2) = embedded_vel(2);
            }

            // Asemble the RHS f_gamma contribution
            for (unsigned int i=0; i<rSplittingData.n_neg; i++)
            {
                unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<TNumNodes; j++)
                {
                    for (unsigned int comp = 0; comp<TDim; comp++)
                    {
                        auxLeftHandSideMatrix(out_node_row_id*BlockSize+comp, j*BlockSize+comp) = f_gamma(i,j);
                    }
                }
            }

            rRightHandSideVector += prod(auxLeftHandSideMatrix, aux_embedded_vel);

        }
    }

    /**
    * This functions adds the level set strong slip boundary no penetration condition imposition.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddSlipNoPenetrationNitcheContribution(MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                const ElementDataType& rData,
                                                const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;                 // Block size
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;   // Matrix size

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
            prev_sol(i*BlockSize+TDim) = rData.p(i);
        }

        // Compute the consistency coefficient
        const double h = rData.h; // Characteristic element size
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);
        const double cons_coef = eff_mu/h;

        // Compute the normal projection matrix
        bounded_matrix<double, TDim, TDim> normal_projection_matrix;
        SetNormalProjectionMatrix(rSplittingData, normal_projection_matrix);

        // Declare auxLeftHandSideNormalMatrix
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Compute the no penetration condition contributions
        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);

            const VectorType aux_out = row(rSplittingData.N_container_out, icut);
            const VectorType aux_int = row(rSplittingData.N_container_int, icut);

            MatrixType aux_out_mat = ZeroMatrix(TDim, rSplittingData.n_neg*TDim);
            MatrixType aux_int_mat = ZeroMatrix(TDim, rSplittingData.n_pos*TDim);

            // Fill the previous auxiliar matrices
            for (unsigned int i=0; i<TDim; ++i)
            {
                for (unsigned int j=0; j<rSplittingData.n_neg; ++j)
                {
                    aux_out_mat(i,j*TDim+i) = aux_out(j);
                }
            }

            for (unsigned int i=0; i<TDim; ++i)
            {
                for (unsigned int j=0; j<rSplittingData.n_pos; ++j)
                {
                    aux_int_mat(i,j*TDim+i) = aux_int(j);
                }
            }

            const MatrixType aux_matrix = prod(trans(aux_out_mat), normal_projection_matrix);
            const MatrixType M_gamma_normal = prod(aux_matrix, aux_out_mat);
            const MatrixType N_gamma_normal = prod(aux_matrix, aux_int_mat);

            // Outer nodes contribution assembly
            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<rSplittingData.n_neg; ++j)
                {
                    const unsigned int out_node_col_id = rSplittingData.out_vec_identifiers[j];

                    for (unsigned int ii=0; ii<TDim; ++ii)
                    {
                        for (unsigned int jj=0; jj<TDim; ++jj)
                        {
                            auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, out_node_col_id*BlockSize+jj) += cons_coef*weight*M_gamma_normal(i*TDim+ii, j*TDim+jj);
                        }
                    }
                }
            }

            // Interior nodes contribution assembly
            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<rSplittingData.n_pos; ++j)
                {
                    const unsigned int int_node_col_id = rSplittingData.int_vec_identifiers[j];

                    for (unsigned int ii=0; ii<TDim; ++ii)
                    {
                        for (unsigned int jj=0; jj<TDim; ++jj)
                        {
                            auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, int_node_col_id*BlockSize+jj) += cons_coef*weight*N_gamma_normal(i*TDim+ii, j*TDim+jj);
                        }
                    }
                }
            }

        }

        // LHS outside Nitche contribution assembly
        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS outside Nitche contribution assembly
        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol);

        // Compute f_gamma if level set velocity is not 0
        if (this->Has(EMBEDDED_VELOCITY))
        {
            auxLeftHandSideMatrix.clear();

            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; i++)
            {
                aux_embedded_vel(i*BlockSize) = embedded_vel(0);
                aux_embedded_vel(i*BlockSize+1) = embedded_vel(1);
                aux_embedded_vel(i*BlockSize+2) = embedded_vel(2);
            }

            // Compute and assemble the level set velocity contribution
            for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
            {
                const double weight = rSplittingData.cut_edge_areas(icut);

                const VectorType aux_out = row(rSplittingData.N_container_out, icut);
                const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

                MatrixType aux_out_mat = ZeroMatrix(TDim, rSplittingData.n_neg*TDim);
                MatrixType aux_cut_mat = ZeroMatrix(TDim, TNumNodes*TDim);

                // Fill the previous auxiliar matrices
                for (unsigned int i=0; i<TDim; ++i)
                {
                    for (unsigned int j=0; j<rSplittingData.n_neg; ++j)
                    {
                        aux_out_mat(i,j*TDim+i) = aux_out(j);
                    }
                }

                for (unsigned int i=0; i<TDim; ++i)
                {
                    for (unsigned int j=0; j<TNumNodes; ++j)
                    {
                        aux_cut_mat(i,j*TDim+i) = aux_cut(j);
                    }
                }

                const MatrixType aux_matrix = prod(trans(aux_out_mat), normal_projection_matrix);
                const MatrixType f_gamma_normal = prod(aux_matrix, aux_cut_mat);

                // Asemble the RHS f_gamma contribution
                for (unsigned int i=0; i<rSplittingData.n_neg; i++)
                {
                    unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                    for (unsigned int j=0; j<TNumNodes; j++)
                    {
                        for (unsigned int ii = 0; ii<TDim; ++ii)
                        {
                            for (unsigned int jj = 0; jj<TDim; ++jj)
                            {
                                auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, j*BlockSize+jj) -= cons_coef*weight*f_gamma_normal(i*TDim+ii,j*TDim+jj);
                            }
                        }
                    }
                }
            }

            rRightHandSideVector += prod(auxLeftHandSideMatrix, aux_embedded_vel);
        }
    }


    /**
    * This functions adds the level set strong slip boundary null tangential stress condition imposition.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddSlipNoTangentialStressNitcheContribution(MatrixType& rLeftHandSideMatrix,
                                                     VectorType& rRightHandSideVector,
                                                     const ElementDataType& rData,
                                                     const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        array_1d<double, TNumNodes*TDim> prev_vel = ZeroVector(TNumNodes*TDim);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
                prev_vel(i*TDim+comp) = rData.v(i,comp);
            }
            prev_sol(i*BlockSize+TDim) = rData.p(i);
        }

        // Compute the effective viscosity with the previous iteration data
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);

        // Declare auxLeftHandSideNormalMatrix
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        bounded_matrix<double, (TDim-1)*3, (TDim-1)*3> voigt_tang_proj_matrix;        // Set the Voigt notation tangential projection matrix
        SetVoigtTangentialProjectionMatrix(rSplittingData, voigt_tang_proj_matrix);

        MatrixType OuterNodesB_matrix;    // Set outer nodes strain matrix (constant since linear elements are used)
        SetOuterNodesStrainMatrix(rData, rSplittingData, OuterNodesB_matrix);

        MatrixType InnerNodesB_matrix;    // Set inner nodes strain matrix (constant since linear elements are used)
        SetInnerNodesStrainMatrix(rData, rSplittingData, InnerNodesB_matrix);

        // Compute the auxiliar matrix that is used in all contributions
        const bounded_matrix<double, (TDim-1)*3, (TDim-1)*3> aux_matrix_1 = prod(trans(rData.C), voigt_tang_proj_matrix);
        const MatrixType aux_matrix_2 = prod(trans(OuterNodesB_matrix), aux_matrix_1);
        const MatrixType aux_matrix_3 = prod(aux_matrix_2, rData.C);

        // Compute the no tangential stress condition contribution matrices
        const MatrixType Aux_M_gamma_tang = prod(aux_matrix_3, OuterNodesB_matrix);
        const MatrixType Aux_N_gamma_tang = prod(aux_matrix_3, InnerNodesB_matrix);

        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);

            // LHS outside nodes null tangential stress contribution assembly
            // Outer nodes contribution assembly
            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<rSplittingData.n_neg; ++j)
                {
                    const unsigned int out_node_col_id = rSplittingData.out_vec_identifiers[j];

                    for (unsigned int ii = 0; ii<TDim; ++ii)
                    {
                        for (unsigned int jj = 0; jj<TDim; ++jj)
                        {
                            auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, out_node_col_id*BlockSize+jj) += (1/eff_mu)*weight*Aux_M_gamma_tang(i*TDim+ii, j*TDim+jj);
                        }
                    }
                }
            }

            // Interior nodes contribution assembly
            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<rSplittingData.n_pos; ++j)
                {
                    const unsigned int int_node_col_id = rSplittingData.int_vec_identifiers[j];

                    for (unsigned int ii = 0; ii<TDim; ++ii)
                    {
                        for (unsigned int jj = 0; jj<TDim; ++jj)
                        {
                            auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, int_node_col_id*BlockSize+jj) += (1/eff_mu)*weight*Aux_N_gamma_tang(i*TDim+ii, j*TDim+jj);
                        }
                    }
                }
            }

            // Pressure contribution assembly
            // Set the current Gauss pt. pressure to Voigt notation matrix
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);
            MatrixType press_to_voigt_mat = ZeroMatrix((TDim-1)*3,BlockSize*TNumNodes);
            for (unsigned int i=0; i<TDim; ++i)
            {
                for (unsigned int j=0; j<TNumNodes; ++j)
                {
                    press_to_voigt_mat(i,j*BlockSize) = aux_cut(j);
                }
            }

            const MatrixType Aux_Pres_gamma_tang = prod(aux_matrix_2, press_to_voigt_mat);

            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<TNumNodes; ++j)
                {
                    for (unsigned int ii = 0; ii<TDim; ++ii)
                    {
                        auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, j*BlockSize) -= weight*Aux_Pres_gamma_tang(i*TDim+ii, j*BlockSize);
                    }
                }
            }
        }

        // LHS outside Nitche contribution assembly
        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS outside Nitche contribution assembly
        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol);

    }


    /**
    * This functions adds the level set strong slip boundary null tangential stress condition imposition.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddSlipNoTangentialStressNitcheContribution2(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      const ElementDataType& rData,
                                                      const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        array_1d<double, TNumNodes*TDim> prev_vel = ZeroVector(TNumNodes*TDim);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
                prev_vel(i*TDim+comp) = rData.v(i,comp);
            }
            prev_sol(i*BlockSize+TDim) = rData.p(i);
        }

        // Declare auxLeftHandSideNormalMatrix
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        bounded_matrix<double, TDim, (TDim-1)*3> voigt_tang_proj_matrix;        // Set the Voigt notation tangential projection matrix
        SetVoigtTangentialProjectionMatrix2(rSplittingData, voigt_tang_proj_matrix);

        MatrixType OuterNodesB_matrix;    // Set outer nodes strain matrix (constant since linear elements are used)
        SetOuterNodesStrainMatrix(rData, rSplittingData, OuterNodesB_matrix);

        MatrixType InnerNodesB_matrix;    // Set inner nodes strain matrix (constant since linear elements are used)
        SetInnerNodesStrainMatrix(rData, rSplittingData, InnerNodesB_matrix);

        // Compute the auxiliar matrix that is used in all contributions
        // const bounded_matrix<double, (TDim-1)*3, (TDim-1)*3> aux_matrix_1 = prod(trans(rData.C), voigt_tang_proj_matrix);
        const bounded_matrix<double, TDim, (TDim-1)*3> aux_matrix = prod(voigt_tang_proj_matrix, rData.C);
        const MatrixType aux_matrix_out = prod(aux_matrix, OuterNodesB_matrix);
        const MatrixType aux_matrix_int = prod(aux_matrix, InnerNodesB_matrix);

        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);

            const VectorType aux_out = row(rSplittingData.N_container_out, icut);
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

            MatrixType shfunc_matrix_out = ZeroMatrix(rSplittingData.n_neg*BlockSize, TDim);
            MatrixType pres_to_voigt_matrix_op = ZeroMatrix((TDim-1)*3, TNumNodes*BlockSize);

            // Fill the test function matrix
            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    shfunc_matrix_out(i*BlockSize+comp, comp) = aux_out(i);
                }
            }

            // Fill the pressure to Voigt notation operator matrix
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    pres_to_voigt_matrix_op(comp,i*BlockSize+TDim) = aux_cut(i);
                }
            }

            const bounded_matrix<double, TDim, TNumNodes*BlockSize> voigt_proj_pressure = prod(aux_matrix, pres_to_voigt_matrix_op);

            MatrixType out_nodes_contribution = prod(shfunc_matrix_out, aux_matrix_out);
            MatrixType int_nodes_contribution = prod(shfunc_matrix_out, aux_matrix_int);
            MatrixType pressure_contribution  = prod(shfunc_matrix_out, voigt_proj_pressure);

            // LHS outside nodes null tangential stress contribution assembly
            // Outer nodes contribution assembly
            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<rSplittingData.n_neg; ++j)
                {
                    const unsigned int out_node_col_id = rSplittingData.out_vec_identifiers[j];

                    for (unsigned int ii = 0; ii<TDim; ++ii)
                    {
                        for (unsigned int jj = 0; jj<TDim; ++jj)
                        {
                            auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, out_node_col_id*BlockSize+jj) += weight*out_nodes_contribution(i*BlockSize+ii, j*TDim+jj);
                        }
                    }
                }
            }

            // Interior nodes contribution assembly
            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<rSplittingData.n_pos; ++j)
                {
                    const unsigned int int_node_col_id = rSplittingData.int_vec_identifiers[j];

                    for (unsigned int ii = 0; ii<TDim; ++ii)
                    {
                        for (unsigned int jj = 0; jj<TDim; ++jj)
                        {
                            auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, int_node_col_id*BlockSize+jj) += weight*int_nodes_contribution(i*BlockSize+ii, j*TDim+jj);
                        }
                    }
                }
            }

            // Pressure contribution assembly
            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<TNumNodes; ++j)
                {
                    for (unsigned int ii = 0; ii<TDim; ++ii)
                    {
                        auxLeftHandSideMatrix(out_node_row_id*BlockSize+ii, j*BlockSize) -= weight*pressure_contribution(i*BlockSize+ii, j*BlockSize);
                    }
                }
            }
        }

        // LHS outside Nitche contribution assembly
        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS outside Nitche contribution assembly
        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol);

    }


    /**
    * This functions adds the pressure equations slip level set BC contributions.
    * The idea is to use the pressure as a Lagrange multiplier for the no penetration condition.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddSlipNoPenetrationNitchePressureContribution(MatrixType& rLeftHandSideMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        const ElementDataType& rData,
                                                        const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
            prev_sol(i*BlockSize+TDim) = rData.p(i);
        }

        // Compute the consistency coefficient
        const double h = rData.h; // Characteristic element size
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);
        const double cons_coef = eff_mu/h;

        // Declare auxLeftHandSideNormalMatrix
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);

            const VectorType aux_out = row(rSplittingData.N_container_out, icut);
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

            for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
            {
                const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

                for (unsigned int j=0; j<TNumNodes; ++j)
                {
                    for (unsigned int comp=0; comp<TDim; ++comp)
                    {
                        auxLeftHandSideMatrix(out_node_row_id*BlockSize+TDim, j*BlockSize+comp) += cons_coef*weight*aux_out(i)*aux_cut(j)*rSplittingData.intersection_normal(comp);
                    }
                }
            }
        }

        // LHS outside Nitche contribution assembly
        // rLeftHandSideMatrix += auxLeftHandSideMatrix;
        rLeftHandSideMatrix -= auxLeftHandSideMatrix;

        // RHS outside Nitche contribution assembly
        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        rRightHandSideVector += prod(auxLeftHandSideMatrix, prev_sol);

        // If level set velocity is not 0, add its contribution to the RHS
        if (this->Has(EMBEDDED_VELOCITY))
        {
            auxLeftHandSideMatrix.clear();

            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; i++)
            {
                aux_embedded_vel(i*BlockSize)   = embedded_vel(0);
                aux_embedded_vel(i*BlockSize+1) = embedded_vel(1);
                aux_embedded_vel(i*BlockSize+2) = embedded_vel(2);
            }

            //TODO: Check this! A clear has been done before
            rRightHandSideVector += prod(auxLeftHandSideMatrix, aux_embedded_vel);
        }
    }


    /**
    * This drops the outer nodes velocity constributions in both LHS and RHS matrices.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void DropOuterNodesVelocityContribution(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ElementDataType& rData,
                                            const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Set the LHS and RHS u_out rows to zero (outside nodes used to impose the BC)
        for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
        {
            const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

            for (unsigned int j=0; j<TDim; ++j)
            {
                // LHS matrix u_out zero set (note that just the velocity rows are set to 0)
                for (unsigned int col = 0; col<MatrixSize; col++)
                {
                    rLeftHandSideMatrix(out_node_row_id*BlockSize+j, col) = 0.0;
                }

                // RHS vector u_out zero set (note that just the velocity rows are set to 0)
                rRightHandSideVector(out_node_row_id*BlockSize+j) = 0.0;
            }
        }

    }


    /**
    * This drops the outer nodes velocity constributions in both LHS and RHS matrices.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void DropOuterNodesPressureContribution(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ElementDataType& rData,
                                            const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Set the LHS and RHS u_out rows to zero (outside nodes used to impose the BC)
        for (unsigned int i=0; i<rSplittingData.n_neg; ++i)
        {
            const unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

            // LHS matrix u_out zero set (note that just the velocity rows are set to 0)
            for (unsigned int col = 0; col<MatrixSize; col++)
            {
                rLeftHandSideMatrix(out_node_row_id*BlockSize+TDim, col) = 0.0;
            }

            // RHS vector u_out zero set (note that just the velocity rows are set to 0)
            rRightHandSideVector(out_node_row_id*BlockSize+TDim) = 0.0;
        }

    }


    /**
    * This function adds the penalty Nitsche term (Freund and Stenberg formulation).
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddSlipNoPenetrationFreundNitschePenaltyContribution(MatrixType& rLeftHandSideMatrix,
                                                              VectorType& rRightHandSideVector,
                                                              const ElementDataType& rData,
                                                              const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
            prev_sol(i*BlockSize+TDim) = rData.p(i);
        }

        // Nitsche coefficient
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);
        const double cons_coef = 10000.0*eff_mu/rData.h;

        // Declare auxiliar arrays
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);
        bounded_matrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Set the normal projection matrix nxn
        bounded_matrix<double, TDim, TDim> normal_projection_matrix;
        SetNormalProjectionMatrix(rSplittingData, normal_projection_matrix);

        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

            // Set the shape functions auxiliar matrices
            bounded_matrix<double, TDim, MatrixSize> N_aux = ZeroMatrix(TDim, MatrixSize);
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    N_aux(comp,i*BlockSize+comp) = aux_cut(i);
                }
            }
            const bounded_matrix<double, MatrixSize, TDim> aux_trans = trans(N_aux);

            // Compute the current cut point auxLHS contribution
            const bounded_matrix<double, MatrixSize, TDim> aux_1 = prod(aux_trans, normal_projection_matrix);
            const bounded_matrix<double, MatrixSize, MatrixSize> aux_2 = prod(aux_1, N_aux);
            auxLeftHandSideMatrix += cons_coef*weight*aux_2;
        }

        // If level set velocity is not 0, add its contribution to the RHS
        if (this->Has(EMBEDDED_VELOCITY))
        {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> embedded_vel_exp = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    embedded_vel_exp(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            auxRightHandSideVector += prod(auxLeftHandSideMatrix, embedded_vel_exp);
        }

        // LHS outside Nitche contribution assembly
        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS outside Nitche contribution assembly
        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        rRightHandSideVector += auxRightHandSideVector;
        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol);

    }


    /**
    * This function adds the Nitsche symmetric counterpart of the fluxes (Freund and Stenberg formulation).
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddSlipNoPenetrationFreundNitscheSymmetricCounterpartContribution(MatrixType& rLeftHandSideMatrix,
                                                                           VectorType& rRightHandSideVector,
                                                                           const ElementDataType& rData,
                                                                           const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
            prev_sol(i*BlockSize+TDim) = rData.p(i);
        }

        // Nitsche coefficient
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);

        // Declare auxiliar arrays
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);
        bounded_matrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Get the normal projection matrix
        bounded_matrix<double, TDim, TDim> normal_projection_matrix = ZeroMatrix(TDim, TDim);
        SetNormalProjectionMatrix(rSplittingData, normal_projection_matrix);

        // Get the normal projection matrix in Voigt notation
        bounded_matrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix = ZeroMatrix(TDim, (TDim-1)*3);
        SetVoigtNormalProjectionMatrix(rSplittingData, voigt_normal_projection_matrix);

        // Get the strain matrix (constant since linear elements are used)
        bounded_matrix<double, (TDim-1)*3, TNumNodes*TDim> B_matrix = ZeroMatrix((TDim-1)*3, TNumNodes*TDim);
        SetStrainMatrix(rData, B_matrix);

        // Expand the B matrix to set 0 in the pressure rows. Besides, transpose it.
        bounded_matrix<double, MatrixSize, (TDim-1)*3> trans_B_matrix_exp = ZeroMatrix(MatrixSize, (TDim-1)*3);
        for (unsigned int i=0; i<TNumNodes; ++i)
        {
            for (unsigned int j=0; j<TDim; ++j)
            {
                for (unsigned int k=0; k<(TDim-1)*3; ++k)
                {
                    trans_B_matrix_exp(i*BlockSize+j, k) = B_matrix(i*TDim+j,k);
                }
            }
        }

        // Compute some element constant matrices
        const bounded_matrix<double, MatrixSize, (TDim-1)*3> aux_matrix_BC = prod(trans_B_matrix_exp, trans(rData.C));

        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

            // Fill the pressure to Voigt notation operator matrix
            bounded_matrix<double, MatrixSize, (TDim-1)*3> trans_pres_to_voigt_matrix_op = ZeroMatrix(MatrixSize, (TDim-1)*3);
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    trans_pres_to_voigt_matrix_op(i*BlockSize+TDim, comp) = aux_cut(i);
                }
            }

            // Set the shape functions auxiliar matrix
            bounded_matrix<double, TDim, MatrixSize> N_aux = ZeroMatrix(TDim, MatrixSize);
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    N_aux(comp,i*BlockSize+comp) = aux_cut(i);
                }
            }

            const bounded_matrix<double, TDim, MatrixSize> N_aux_proj = prod(normal_projection_matrix, N_aux);
            const bounded_matrix<double, (TDim-1)*3, MatrixSize> N_aux_proj_voigt = prod(trans(voigt_normal_projection_matrix), N_aux_proj);

            // Contribution coming fron the shear stress operator
            //TODO: CHECK THE CONSISTENCY PARAMETERS
            auxLeftHandSideMatrix += eff_mu*weight*prod(aux_matrix_BC, N_aux_proj_voigt);

            // Contribution coming from the pressure terms
            //TODO: CHECK THE CONSISTENCY PARAMETERS
            auxLeftHandSideMatrix -= eff_mu*weight*prod(trans_pres_to_voigt_matrix_op, N_aux_proj_voigt);

        }

        // LHS outside Nitche contribution assembly
        rLeftHandSideMatrix -= auxLeftHandSideMatrix; // The minus sign comes from the Nitsche formulation

        // RHS outside Nitche contribution assembly
        // If level set velocity is not 0, add its contribution to the RHS
        if (this->Has(EMBEDDED_VELOCITY))
        {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> embedded_vel_exp = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    embedded_vel_exp(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            auxRightHandSideVector += prod(auxLeftHandSideMatrix, embedded_vel_exp);
        }

        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        rRightHandSideVector -= auxRightHandSideVector;
        rRightHandSideVector += prod(auxLeftHandSideMatrix, prev_sol); // The plus sign comes from the Nitsche formulation

    }


    /**
    * This function adds the null tangential stress contribution (Freund and Stenberg formulation).
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddSlipNoPenetrationFreundNitscheNullTangentialStressContribution(MatrixType& rLeftHandSideMatrix,
                                                                           VectorType& rRightHandSideVector,
                                                                           const ElementDataType& rData,
                                                                           const ElementSplittingDataStruct& rSplittingData)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);

        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int comp=0; comp<TDim; comp++)
            {
                prev_sol(i*BlockSize+comp) = rData.v(i,comp);
            }
            prev_sol(i*BlockSize+TDim) = rData.p(i);
        }

        // Declare auxiliar arrays
        bounded_matrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Get the normal projection matrix
        bounded_matrix<double, TDim, TDim> normal_projection_matrix = ZeroMatrix(TDim, TDim);
        SetNormalProjectionMatrix(rSplittingData, normal_projection_matrix);

        // Get the normal projection matrix in Voigt notation
        bounded_matrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix = ZeroMatrix(TDim, (TDim-1)*3);
        SetVoigtNormalProjectionMatrix(rSplittingData, voigt_normal_projection_matrix);

        // Get the strain matrix (constant since linear elements are used)
        bounded_matrix<double, (TDim-1)*3, TNumNodes*TDim> B_matrix = ZeroMatrix((TDim-1)*3, TNumNodes*TDim);
        SetStrainMatrix(rData, B_matrix);

        // Expand the B matrix to set 0 in the pressure rows.
        bounded_matrix<double, (TDim-1)*3, MatrixSize> B_matrix_exp = ZeroMatrix((TDim-1)*3, MatrixSize);
        for (unsigned int i=0; i<TNumNodes; ++i)
        {
            for (unsigned int j=0; j<TDim; ++j)
            {
                for (unsigned int k=0; k<(TDim-1)*3; ++k)
                {
                    B_matrix_exp(k, i*BlockSize+j) = B_matrix(k, i*TDim+j);
                }
            }
        }

        // Compute some element constant matrices
        const bounded_matrix<double, TDim, (TDim-1)*3> projection_matrix = prod(normal_projection_matrix, voigt_normal_projection_matrix);
        const bounded_matrix<double, (TDim-1)*3, MatrixSize> aux_matrix_BC = prod(rData.C, B_matrix_exp);

        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            const double weight = rSplittingData.cut_edge_areas(icut);
            const VectorType aux_cut = row(rSplittingData.N_container_cut, icut);

            // Fill the pressure to Voigt notation operator matrix
            bounded_matrix<double, (TDim-1)*3, MatrixSize> pres_to_voigt_matrix_op = ZeroMatrix((TDim-1)*3, MatrixSize);
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    pres_to_voigt_matrix_op(comp, i*BlockSize+TDim) = aux_cut(i);
                }
            }

            // Set the shape functions auxiliar transpose matrix
            bounded_matrix<double, MatrixSize, TDim> N_aux_trans = ZeroMatrix(MatrixSize, TDim);
            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    N_aux_trans(i*BlockSize+comp, comp) = aux_cut(i);
                }
            }

            const bounded_matrix<double, MatrixSize, (TDim-1)*3> N_aux_proj = prod(N_aux_trans, projection_matrix);

            // Contribution coming fron the shear stress operator
            auxLeftHandSideMatrix += weight*prod(N_aux_proj, aux_matrix_BC);

            // Contribution coming from the pressure terms
            auxLeftHandSideMatrix -= weight*prod(N_aux_proj, pres_to_voigt_matrix_op);

        }

        // LHS outside Nitche contribution assembly
        rLeftHandSideMatrix -= auxLeftHandSideMatrix; // The minus sign comes from the Nitsche formulation

        // RHS outside Nitche contribution assembly
        rRightHandSideVector += prod(auxLeftHandSideMatrix, prev_sol); // The plus sign comes from the Nitsche formulation

    }


    /**
    * This functions collects and adds all the level set boundary condition contributions
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddBoundaryConditionElementContribution(MatrixType& rLeftHandSideMatrix,
                                                 VectorType& rRightHandSideVector,
                                                 const array_1d<double,TNumNodes>& rDistances,
                                                 const ElementDataType& rData,
                                                 const ElementSplittingDataStruct& rSplittingData)
    {

        // TODO: CREATE A METHOD TO COMPUTE THE PREVIOUS SOLUTION ONCE

        // Add all the boundary intersection terms in those elements that are split
        if (rSplittingData.ndivisions > 1)
        {
            // Compute and assemble the boundary terms comping from the integration by parts
            AddIntersectionBoundaryTermsContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);

            // First, compute and assemble the penalty level set BC imposition contribution
            // Secondly, compute and assemble the modified Nitche method level set BC imposition contribution
            // Note that the Nitche contribution has to be computed the last since it drops the outer nodes rows previous constributions
            if (this->Is(SLIP)) // TODO: Now, all the elements are marked with the SLIP flag. Even though the contribution is only added
            {
                // Previous implementation and tests
                // AddSlipBoundaryConditionPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData); //TODO: Check and decoment this when finished
                // DropOuterNodesVelocityContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
                // // DropOuterNodesPressureContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
                // AddSlipNoPenetrationNitcheContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
                // AddSlipNoTangentialStressNitcheContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData); // Joan and Ramon
                // // AddSlipNoTangentialStressNitcheContribution2(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);   // Just testing against W
                // // AddSlipNoPenetrationNitchePressureContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);

                // Pure Nitche implementation (Freund and Stenberg)
                AddSlipNoPenetrationFreundNitschePenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
                AddSlipNoPenetrationFreundNitscheSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
                AddSlipNoPenetrationFreundNitscheNullTangentialStressContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
            }
            else
            {
                AddBoundaryConditionPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
                DropOuterNodesVelocityContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
                AddBoundaryConditionModifiedNitcheContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
            }
        }
        // Add a penalty contribution to enforce the BC imposition at the level set when the distance value is close to 0
        // Note that this is an excepcional case in where there are both positive and negative distance value but the intersection
        // utility determines that there are no intersections (this might happen if the level set is too close to a node).
        // TODO: Add slip version
        // else if (rSplittingData.ndivisions == 1)
        // {
        //     constexpr unsigned int BlockSize = TDim+1;                 // Block size
        //     constexpr unsigned int MatrixSize = TNumNodes*BlockSize;   // Matrix size
        //
        //     double diag_max = 0.0;
        //
        //     for (unsigned int i=0; i<MatrixSize; i++)
        //     {
        //         if ((fabs(rLeftHandSideMatrix(i,i)) > diag_max) && ((i+1)%BlockSize != 0))
        //         {
        //             diag_max = fabs(rLeftHandSideMatrix(i,i));
        //         }
        //     }
        //
        //     double tol_d;
        //
        //     if (TDim == 2)
        //     {
        //         tol_d = 1e-2*sqrt(rData.h);
        //     }
        //     else
        //     {
        //         tol_d = 1e-2*pow(rData.h, 1.0/3.0);
        //     }
        //
        //     double pen_coef = std::max((diag_max/(0.1*rData.h*rData.h)),(1000*rData.h*rData.h));
        //
        //     for (unsigned int i=0; i<TNumNodes; i++)
        //     {
        //         if (fabs(rDistances[i])<tol_d)
        //         {
        //             // LHS penalty contribution assembly
        //             for (unsigned int comp=0; comp<TDim; comp++)
        //             {
        //                 rLeftHandSideMatrix(i*BlockSize+comp,i*BlockSize+comp) += pen_coef;
        //             }
        //
        //             // RHS penalty contribution assembly
        //             if (this->Has(EMBEDDED_VELOCITY))
        //             {
        //                 const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        //                 for (unsigned int comp = 0; comp<TDim; comp++)
        //                 {
        //                     rRightHandSideVector(i*BlockSize+comp) += pen_coef*embedded_vel(comp);
        //                 }
        //             }
        //
        //             // RHS residual contribution assembly
        //             for (unsigned int comp = 0; comp<TDim; comp++)
        //             {
        //                 rRightHandSideVector(i*BlockSize+comp) -= pen_coef*rData.v(i,comp);
        //             }
        //         }
        //     }
        // }
    }


    /**
    * This functions sets the intersection normal matrix
    * @param rIntersection_normal: reference to intersection normal vector
    * @param rNormal_matrix: reference to the computed intersection normal matrix
    */
    void SetNormalMatrix(const VectorType& rIntersection_normal,
                         bounded_matrix<double, TDim, (TDim-1)*3>& rNormal_matrix)
    {
        rNormal_matrix.clear();

        if (TDim == 3)
        {
            rNormal_matrix(0,0) = rIntersection_normal(0);
            rNormal_matrix(0,3) = rIntersection_normal(1);
            rNormal_matrix(0,5) = rIntersection_normal(2);
            rNormal_matrix(1,1) = rIntersection_normal(1);
            rNormal_matrix(1,3) = rIntersection_normal(0);
            rNormal_matrix(1,4) = rIntersection_normal(2);
            rNormal_matrix(2,2) = rIntersection_normal(2);
            rNormal_matrix(2,4) = rIntersection_normal(1);
            rNormal_matrix(2,5) = rIntersection_normal(0);
        }
        else
        {
            rNormal_matrix(0,0) = rIntersection_normal(0);
            rNormal_matrix(0,2) = rIntersection_normal(1);
            rNormal_matrix(1,1) = rIntersection_normal(1);
            rNormal_matrix(1,2) = rIntersection_normal(0);
        }
    }


    /**
    * This functions sets the B strain matrix
    * @param rData: reference to element data structure (it contains the shape functions derivatives)
    * @param rB_matrix: reference to the computed B strain matrix
    */
    void SetStrainMatrix(const ElementDataType& rData,
                         bounded_matrix<double, (TDim-1)*3, TNumNodes*TDim>& rB_matrix)
    {
        rB_matrix.clear();

        if (TDim == 3)
        {
            for (unsigned int i=0; i<TNumNodes; i++)
            {
                rB_matrix(0,i*TDim)   = rData.DN_DX(i,0);
                rB_matrix(1,i*TDim+1) = rData.DN_DX(i,1);
                rB_matrix(2,i*TDim+2) = rData.DN_DX(i,2);
                rB_matrix(3,i*TDim)   = rData.DN_DX(i,1);
                rB_matrix(3,i*TDim+1) = rData.DN_DX(i,0);
                rB_matrix(4,i*TDim+1) = rData.DN_DX(i,2);
                rB_matrix(4,i*TDim+2) = rData.DN_DX(i,1);
                rB_matrix(5,i*TDim)   = rData.DN_DX(i,2);
                rB_matrix(5,i*TDim+2) = rData.DN_DX(i,0);
            }
        }
        else
        {
            for (unsigned int i=0; i<TNumNodes; i++)
            {
                rB_matrix(0,i*TDim)   = rData.DN_DX(i,0);
                rB_matrix(1,i*TDim+1) = rData.DN_DX(i,1);
                rB_matrix(2,i*TDim)   = rData.DN_DX(i,1);
                rB_matrix(2,i*TDim+1) = rData.DN_DX(i,0);
            }
        }
    }


    /**
    * This functions sets the B strain matrix for the outer nodes
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    * @param rOuterNodesB_matrix: reference to the computed B strain matrix
    */
    void SetOuterNodesStrainMatrix(const ElementDataType& rData,
                                   const ElementSplittingDataStruct& rSplittingData,
                                   MatrixType& rOuterNodesB_matrix)
    {
        rOuterNodesB_matrix = ZeroMatrix((TDim-1)*3, rSplittingData.n_neg*TDim);

        if (TDim == 3)
        {
            for (unsigned int i=0; i<rSplittingData.n_neg; i++)
            {
                const unsigned int out_node_id = rSplittingData.out_vec_identifiers[i];

                rOuterNodesB_matrix(0,i*TDim)   = rData.DN_DX(out_node_id,0);
                rOuterNodesB_matrix(1,i*TDim+1) = rData.DN_DX(out_node_id,1);
                rOuterNodesB_matrix(2,i*TDim+2) = rData.DN_DX(out_node_id,2);
                rOuterNodesB_matrix(3,i*TDim)   = rData.DN_DX(out_node_id,1);
                rOuterNodesB_matrix(3,i*TDim+1) = rData.DN_DX(out_node_id,0);
                rOuterNodesB_matrix(4,i*TDim+1) = rData.DN_DX(out_node_id,2);
                rOuterNodesB_matrix(4,i*TDim+2) = rData.DN_DX(out_node_id,1);
                rOuterNodesB_matrix(5,i*TDim)   = rData.DN_DX(out_node_id,2);
                rOuterNodesB_matrix(5,i*TDim+2) = rData.DN_DX(out_node_id,0);
            }
        }
        else
        {
            for (unsigned int i=0; i<rSplittingData.n_neg; i++)
            {
                const unsigned int out_node_id = rSplittingData.out_vec_identifiers[i];

                rOuterNodesB_matrix(0,i*TDim)   = rData.DN_DX(out_node_id,0);
                rOuterNodesB_matrix(1,i*TDim+1) = rData.DN_DX(out_node_id,1);
                rOuterNodesB_matrix(2,i*TDim)   = rData.DN_DX(out_node_id,1);
                rOuterNodesB_matrix(2,i*TDim+1) = rData.DN_DX(out_node_id,0);
            }
        }
    }


    /**
    * This functions sets the B strain matrix for the inner nodes
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    * @param rInnerNodesB_matrix: reference to the computed B strain matrix
    */
    void SetInnerNodesStrainMatrix(const ElementDataType& rData,
                                   const ElementSplittingDataStruct& rSplittingData,
                                   MatrixType& rInnerNodesB_matrix)
    {
        rInnerNodesB_matrix = ZeroMatrix((TDim-1)*3, rSplittingData.n_pos*TDim);

        if (TDim == 3)
        {
            for (unsigned int i=0; i<rSplittingData.n_pos; i++)
            {
                const unsigned int int_node_id = rSplittingData.int_vec_identifiers[i];

                rInnerNodesB_matrix(0,i*TDim)   = rData.DN_DX(int_node_id,0);
                rInnerNodesB_matrix(1,i*TDim+1) = rData.DN_DX(int_node_id,1);
                rInnerNodesB_matrix(2,i*TDim+2) = rData.DN_DX(int_node_id,2);
                rInnerNodesB_matrix(3,i*TDim)   = rData.DN_DX(int_node_id,1);
                rInnerNodesB_matrix(3,i*TDim+1) = rData.DN_DX(int_node_id,0);
                rInnerNodesB_matrix(4,i*TDim+1) = rData.DN_DX(int_node_id,2);
                rInnerNodesB_matrix(4,i*TDim+2) = rData.DN_DX(int_node_id,1);
                rInnerNodesB_matrix(5,i*TDim)   = rData.DN_DX(int_node_id,2);
                rInnerNodesB_matrix(5,i*TDim+2) = rData.DN_DX(int_node_id,0);
            }
        }
        else
        {
            for (unsigned int i=0; i<rSplittingData.n_pos; i++)
            {
                const unsigned int int_node_id = rSplittingData.int_vec_identifiers[i];

                rInnerNodesB_matrix(0,i*TDim)   = rData.DN_DX(int_node_id,0);
                rInnerNodesB_matrix(1,i*TDim+1) = rData.DN_DX(int_node_id,1);
                rInnerNodesB_matrix(2,i*TDim)   = rData.DN_DX(int_node_id,1);
                rInnerNodesB_matrix(2,i*TDim+1) = rData.DN_DX(int_node_id,0);
            }
        }
    }

    /**
    * This functions sets the auxiliar matrix to compute the tangential projection in Voigt notation
    * @param rSplittingData: reference to the intersection data structure
    * @param rVoigtTangProj_matrix: reference to the computed tangential projection auxiliar matrix
    */
    void SetNormalProjectionMatrix(const ElementSplittingDataStruct& rSplittingData,
                                   bounded_matrix<double, TDim, TDim>& rNormProj_matrix)
    {
        rNormProj_matrix.clear();

        if (TDim == 3)
        {
            // Fill the normal projection matrix (nxn)
            rNormProj_matrix = outer_prod(rSplittingData.intersection_normal, rSplittingData.intersection_normal);
        }
        else
        {
            // Fill the normal projection matrix (nxn)
            rNormProj_matrix(0,0) = rSplittingData.intersection_normal(0)*rSplittingData.intersection_normal(0);
            rNormProj_matrix(0,1) = rSplittingData.intersection_normal(0)*rSplittingData.intersection_normal(1);
            rNormProj_matrix(1,0) = rSplittingData.intersection_normal(1)*rSplittingData.intersection_normal(0);
            rNormProj_matrix(1,1) = rSplittingData.intersection_normal(1)*rSplittingData.intersection_normal(1);
        }
    }

    /**
    * This functions sets the auxiliar matrix to compute the tangential projection in Voigt notation
    * @param rSplittingData: reference to the intersection data structure
    * @param rVoigtTangProj_matrix: reference to the computed tangential projection auxiliar matrix
    */
    void SetVoigtNormalProjectionMatrix(const ElementSplittingDataStruct& rSplittingData,
                                        bounded_matrix<double, TDim, (TDim-1)*3>& rVoigtNormProj_matrix)
    {
        rVoigtNormProj_matrix.clear();

        if (TDim == 3)
        {
            // Fill the normal projection matrix for Voigt notation
            rVoigtNormProj_matrix(0,0) = rSplittingData.intersection_normal(0);
            rVoigtNormProj_matrix(0,3) = rSplittingData.intersection_normal(1);
            rVoigtNormProj_matrix(0,5) = rSplittingData.intersection_normal(2);
            rVoigtNormProj_matrix(1,1) = rSplittingData.intersection_normal(1);
            rVoigtNormProj_matrix(1,3) = rSplittingData.intersection_normal(0);
            rVoigtNormProj_matrix(1,4) = rSplittingData.intersection_normal(2);
            rVoigtNormProj_matrix(2,2) = rSplittingData.intersection_normal(2);
            rVoigtNormProj_matrix(2,4) = rSplittingData.intersection_normal(1);
            rVoigtNormProj_matrix(2,5) = rSplittingData.intersection_normal(0);
        }
        else
        {
            // Fill the noromal projection matrix for Voigt notation
            rVoigtNormProj_matrix(0,0) = rSplittingData.intersection_normal(0);
            rVoigtNormProj_matrix(0,2) = rSplittingData.intersection_normal(1);
            rVoigtNormProj_matrix(1,1) = rSplittingData.intersection_normal(1);
            rVoigtNormProj_matrix(1,2) = rSplittingData.intersection_normal(0);

        }
    }

    /**
    * This functions sets the auxiliar matrix to compute the tangential projection in Voigt notation
    * @param rSplittingData: reference to the intersection data structure
    * @param rVoigtTangProj_matrix: reference to the computed tangential projection auxiliar matrix
    */
    void SetVoigtTangentialProjectionMatrix(const ElementSplittingDataStruct& rSplittingData,
                                            bounded_matrix<double, (TDim-1)*3, (TDim-1)*3>& rVoigtTangProj_matrix)
    {
        rVoigtTangProj_matrix.clear();

        bounded_matrix<double, TDim, TDim> tang_vector_matrix;
        bounded_matrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix;

        if (TDim == 3)
        {
            // Fill the tangential projection matrix (I-nxn)
            const bounded_matrix<double, TDim, TDim> norm_vector_matrix = outer_prod(rSplittingData.intersection_normal,
                                                                                     rSplittingData.intersection_normal);
            tang_vector_matrix.clear();
            for (unsigned int i=0; i<TDim; ++i) {tang_vector_matrix(i,i) = 1.0;}
            tang_vector_matrix -= norm_vector_matrix;

            // Fill the Voigt notation normal projection matrix
            voigt_normal_projection_matrix.clear();
            voigt_normal_projection_matrix(0,0) = rSplittingData.intersection_normal(0);
            voigt_normal_projection_matrix(0,3) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(0,5) = rSplittingData.intersection_normal(2);
            voigt_normal_projection_matrix(1,1) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(1,3) = rSplittingData.intersection_normal(0);
            voigt_normal_projection_matrix(1,4) = rSplittingData.intersection_normal(2);
            voigt_normal_projection_matrix(2,2) = rSplittingData.intersection_normal(2);
            voigt_normal_projection_matrix(2,4) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(2,5) = rSplittingData.intersection_normal(0);
        }
        else
        {
            // Fill the tangential projection matrix (I-nxn)
            tang_vector_matrix(0,0) = 1.0 - rSplittingData.intersection_normal(0)*rSplittingData.intersection_normal(0);
            tang_vector_matrix(0,1) =     - rSplittingData.intersection_normal(0)*rSplittingData.intersection_normal(1);
            tang_vector_matrix(1,0) =     - rSplittingData.intersection_normal(1)*rSplittingData.intersection_normal(0);
            tang_vector_matrix(1,1) = 1.0 - rSplittingData.intersection_normal(1)*rSplittingData.intersection_normal(1);

            // Fill the Voigt notation normal projection matrix
            voigt_normal_projection_matrix.clear();
            voigt_normal_projection_matrix(0,0) = rSplittingData.intersection_normal(0);
            voigt_normal_projection_matrix(0,2) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(1,1) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(1,2) = rSplittingData.intersection_normal(0);
        }

        const bounded_matrix<double, TDim, (TDim-1)*3> aux_proj_voigt = prod(tang_vector_matrix, voigt_normal_projection_matrix);
        rVoigtTangProj_matrix = prod(trans(voigt_normal_projection_matrix), aux_proj_voigt);
    }

    /**
    * This functions sets the auxiliar matrix to compute the tangential projection in Voigt notation
    * @param rSplittingData: reference to the intersection data structure
    * @param rVoigtTangProj_matrix: reference to the computed tangential projection auxiliar matrix
    */
    void SetVoigtTangentialProjectionMatrix2(const ElementSplittingDataStruct& rSplittingData,
                                             bounded_matrix<double, TDim, (TDim-1)*3>& rVoigtTangProj_matrix)
    {
        rVoigtTangProj_matrix.clear();

        bounded_matrix<double, TDim, TDim> tang_vector_matrix;
        bounded_matrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix;

        if (TDim == 3)
        {
            // Fill the tangential projection matrix (I-nxn)
            const bounded_matrix<double, TDim, TDim> norm_vector_matrix = outer_prod(rSplittingData.intersection_normal,
                                                                                     rSplittingData.intersection_normal);
            tang_vector_matrix.clear();
            for (unsigned int i=0; i<TDim; ++i) {tang_vector_matrix(i,i) = 1.0;}
            tang_vector_matrix -= norm_vector_matrix;

            // Fill the Voigt notation normal projection matrix
            voigt_normal_projection_matrix.clear();
            voigt_normal_projection_matrix(0,0) = rSplittingData.intersection_normal(0);
            voigt_normal_projection_matrix(0,3) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(0,5) = rSplittingData.intersection_normal(2);
            voigt_normal_projection_matrix(1,1) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(1,3) = rSplittingData.intersection_normal(0);
            voigt_normal_projection_matrix(1,4) = rSplittingData.intersection_normal(2);
            voigt_normal_projection_matrix(2,2) = rSplittingData.intersection_normal(2);
            voigt_normal_projection_matrix(2,4) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(2,5) = rSplittingData.intersection_normal(0);
        }
        else
        {
            // Fill the tangential projection matrix (I-nxn)
            tang_vector_matrix(0,0) = 1.0 - rSplittingData.intersection_normal(0)*rSplittingData.intersection_normal(0);
            tang_vector_matrix(0,1) =     - rSplittingData.intersection_normal(0)*rSplittingData.intersection_normal(1);
            tang_vector_matrix(1,0) =     - rSplittingData.intersection_normal(1)*rSplittingData.intersection_normal(0);
            tang_vector_matrix(1,1) = 1.0 - rSplittingData.intersection_normal(1)*rSplittingData.intersection_normal(1);

            // Fill the Voigt notation normal projection matrix
            voigt_normal_projection_matrix.clear();
            voigt_normal_projection_matrix(0,0) = rSplittingData.intersection_normal(0);
            voigt_normal_projection_matrix(0,2) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(1,1) = rSplittingData.intersection_normal(1);
            voigt_normal_projection_matrix(1,2) = rSplittingData.intersection_normal(0);
        }

        rVoigtTangProj_matrix = prod(tang_vector_matrix, voigt_normal_projection_matrix);
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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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

///@}

} // namespace Kratos.

#endif // KRATOS_STOKES_ELEMENT_SYMBOLIC_INCLUDED  defined
