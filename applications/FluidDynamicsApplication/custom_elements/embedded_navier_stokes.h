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
//~ #include "utilities/enrichment_utilities_duplicate_dofs.h"  // Tetrahedra splitting
#include "utilities/split_tetrahedra_utilities.h"               // Tetrahedra splitting
//~ #include "utilities/enrich_2d_2dofs.h"                      // Triangle splitting
//~ #include "utilities/discont_utils.h"                        // Tetrahedra splitting

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

        // Allocate memory needed
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_local;
        array_1d<double, MatrixSize> rhs_local;

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
            ComputeElementAsFluid<MatrixSize>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo);
        }
        else if(nneg == TNumNodes) // All nodes belong to structure domain
        {
            rLeftHandSideMatrix.clear();
            rRightHandSideVector.clear();
        }
        else // Element intersects both fluid and structure domains
        {
            ComputeElementAsMixed<MatrixSize>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo, distances);
        }

        KRATOS_CATCH("Error in embedded Navier-Stokes symbolic element")
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
            KRATOS_THROW_ERROR(std::invalid_argument,"DISTANCE Key is 0. Check if the application was correctly registered.","");

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
    void ComputeElementAsFluid(bounded_matrix<double,MatrixSize, MatrixSize>& lhs_local,
                               array_1d<double,MatrixSize>& rhs_local,
                               MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               ElementDataType& rData,
                               ProcessInfo& rCurrentProcessInfo)
    {
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
    void ComputeElementAsMixed(bounded_matrix<double,MatrixSize, MatrixSize>& lhs_local,
                               array_1d<double,MatrixSize>& rhs_local,
                               MatrixType& rLeftHandSideMatrix,
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

        //~ std::cout << "Mixed element: " << this->Id() << " ndivisions = " << ndivisions << " distances: " << rDistances[0] << " " << rDistances[1] << " " << rDistances[2] << " " << distances[3] << " " << std::endl;

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
                ComputeElementAsFluid<MatrixSize>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, rData, rCurrentProcessInfo);
            }
        }
        else
        {
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
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;    // Matrix size

        array_1d<double, TNumNodes> aux;
        array_1d<double, MatrixSize> aux_normal;

        VectorType auxRightHandSideVector = ZeroVector(MatrixSize);
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, TNumNodes);

        // Compute the contribution of the boundary pressure term at the intersection
        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {

            aux = row(rSplittingData.N_container_cut, icut);
            aux_normal = ZeroVector(MatrixSize);

            double weight = (rSplittingData.cut_edge_areas)(icut);

            for (unsigned int i=0; i<TNumNodes; i++)
            {
                for(unsigned int j=0; j<TDim; j++)
                {
                    aux_normal(i*BlockSize+j) = aux(i)*rSplittingData.intersection_normal(j);
                }
            }

            auxLeftHandSideMatrix += weight*outer_prod(aux_normal, aux);

        }

        // Assemble the LHS contribution to the pressure columns
        for (unsigned int i=0; i<MatrixSize; i++)
        {
            for (unsigned int j=0; j<TNumNodes; j++)
            {
                rLeftHandSideMatrix(i,j*BlockSize+TDim) += auxLeftHandSideMatrix(i,j);
            }
        }

        // Assemble and compute the RHS pressure residual contribution at the intersection
        noalias(rRightHandSideVector) -= prod(auxLeftHandSideMatrix, rData.p);

        // Compute the contribution of the tangential boundary stress term at the intersection
        bounded_matrix<double, TDim, (TDim-1)*3> normal_matrix;         // Set normal matrix
        SetNormalMatrix(rSplittingData.intersection_normal, normal_matrix);

        bounded_matrix<double, (TDim-1)*3, TNumNodes*TDim> B_matrix;    // Set strain matrix
        SetStrainMatrix(rData, B_matrix);

        bounded_matrix<double, TNumNodes*TDim, TDim> w_matrix;
        bounded_matrix<double, TNumNodes*TDim, (TDim-1)*3> aux_matrix;

        auxLeftHandSideMatrix.resize(TNumNodes*TDim, TNumNodes*TDim, false);
        auxRightHandSideVector.resize(TNumNodes*TDim, false);
        auxLeftHandSideMatrix.clear();
        auxRightHandSideVector.clear();

        // Compute the tangential stress LHS contribution
        for (unsigned int icut=0; icut<rSplittingData.ncutpoints; icut++)
        {
            aux = row(rSplittingData.N_container_cut, icut);
            double weight = (rSplittingData.cut_edge_areas)(icut);

            w_matrix.clear();
            aux_matrix.clear();

            for (unsigned int i=0; i<TNumNodes; i++)
            {
                for (unsigned int j=0; j<TDim; j++)
                {
                    w_matrix(i*TDim+j,j) = aux(i);
                }
            }

            noalias(aux_matrix) = prod(w_matrix, normal_matrix);
            aux_matrix = prod(aux_matrix, rData.C);

            noalias(auxLeftHandSideMatrix) += weight*prod(aux_matrix, B_matrix);

        }

        // Assemble the LHS tangential stress contribution to the velocity rows
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int j=0; j<TNumNodes; j++)
            {
                for (unsigned int row=0; row<TDim; row++)
                {
                    for (unsigned int col=0; col<TDim; col++)
                    {
                        //~ rLeftHandSideMatrix(i*BlockSize+row,j*BlockSize+col) += auxLeftHandSideMatrix(i*TDim+row, j*TDim+col); // TODO: Check these signs
                        rLeftHandSideMatrix(i*BlockSize+row,j*BlockSize+col) -= auxLeftHandSideMatrix(i*TDim+row, j*TDim+col);
                    }
                }
            }
        }

        // Obtain the previous iteration velocity solution
        array_1d<double, TNumNodes*TDim> prev_sol = ZeroVector(TNumNodes*TDim);

        for (unsigned int i = 0; i<TNumNodes; i++)
        {
            prev_sol(i*TDim) = rData.v(i,0);
            prev_sol(i*TDim+1) = rData.v(i,1);
            prev_sol(i*TDim+2) = rData.v(i,2);
        }

        noalias(auxRightHandSideVector) = prod(auxLeftHandSideMatrix, prev_sol);

        // Assemble the RHS tangential stress contribution to the velocity rows
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int row=0; row<TDim; row++)
            {
                //~ rRightHandSideVector(i*BlockSize+row) -= auxRightHandSideVector(i*TDim+row); // TODO: Check these signs
                rRightHandSideVector(i*BlockSize+row) += auxRightHandSideVector(i*TDim+row);
            }
        }

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

        MatrixType P_gamma(TNumNodes, TNumNodes);    // Penalty matrix
        noalias(P_gamma) = ZeroMatrix(TNumNodes, TNumNodes);

        VectorType aux_cut(TNumNodes);

        double intersection_area = 0.0;

        for (unsigned int icut = 0; icut<rSplittingData.ncutpoints; icut++)
        {
            double weight = rSplittingData.cut_edge_areas(icut);
            intersection_area += weight;

            aux_cut = row(rSplittingData.N_container_cut, icut);

            P_gamma += weight*outer_prod(aux_cut,aux_cut);
        }

        // Compute the penalty coefficient
        double diag_max = 0.0;
        for (unsigned int i=0; i<MatrixSize; i++)
        {
            if ((fabs(rLeftHandSideMatrix(i,i)) > diag_max) && (i%BlockSize != 0.0))
            {
                diag_max = fabs(rLeftHandSideMatrix(i,i)); // Maximum diagonal value (associated to velocity)
            }
        }

        // TODO: Think about this value. Now is K*max(LHS(i,i))*IntArea (we integrate P_gamma over the intersection area)
        double h = rData.h;
        double denominator = std::max(0.0001*h*h, intersection_area);

        double pen_coef = 100.0*diag_max/denominator;

        // Multiply the penalty matrix by the penalty coefficient
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
    * This functions adds the level set strong boundary condition imposition contribution.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rSplittingData: reference to the intersection data structure
    */
    void AddBoundaryConditionNitcheContribution(MatrixType& rLeftHandSideMatrix,
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

        // Set the u_out rows to zero (outside nodes used to impose the BC)
        for (unsigned int i = 0; i<rSplittingData.n_neg; i++)
        {
            unsigned int out_node_row_id = rSplittingData.out_vec_identifiers[i];

            for (unsigned int j = 0; j<TDim; j++)
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

        // Add all the boundary intersection terms
        if (rSplittingData.ndivisions > 1)
        {
            // Compute and assemble the boundary terms level set contribution
            AddIntersectionBoundaryTermsContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);

            // Compute and assemble the penalty boundary condition imposition contribution
            // TODO: AddSlipBoundaryConditionPenaltyContribution
            AddBoundaryConditionPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);

            // Compute and assemble the modified Nitche method boundary condition imposition contribution
            // Note that the Nitche contribution has to be computed the last since it modifies the outer nodes contribution
            // TODO: AddSlipBoundaryConditionNitcheContribution
            AddBoundaryConditionNitcheContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rSplittingData);
        }
        // Add a penalty contribution to enforce the BC imposition at the level set when the distance value is close to 0
        // Note that this is an excepcional case in where there are both positive and negative distance value but the intersection
        // utility determines that there are no intersections (this might happen if the level set is too close to a node).
        // TODO: Add slip version
        else if (rSplittingData.ndivisions == 1)
        {
            constexpr unsigned int BlockSize = TDim+1;                 // Block size
            constexpr unsigned int MatrixSize = TNumNodes*BlockSize;   // Matrix size

            double diag_max = 0.0;

            for (unsigned int i=0; i<MatrixSize; i++)
            {
                if ((fabs(rLeftHandSideMatrix(i,i)) > diag_max) && ((i+1)%BlockSize != 0))
                {
                    diag_max = fabs(rLeftHandSideMatrix(i,i));
                }
            }

            double tol_d;

            if (TDim == 2)
            {
                tol_d = 1e-2*sqrt(rData.h);
            }
            else
            {
                tol_d = 1e-2*pow(rData.h, 1.0/3.0);
            }

            double pen_coef = std::max((diag_max/(0.1*rData.h*rData.h)),(1000*rData.h*rData.h));

            for (unsigned int i=0; i<TNumNodes; i++)
            {
                if (fabs(rDistances[i])<tol_d)
                {
                    // LHS penalty contribution assembly
                    for (unsigned int comp=0; comp<TDim; comp++)
                    {
                        rLeftHandSideMatrix(i*BlockSize+comp,i*BlockSize+comp) += pen_coef;
                    }

                    // RHS penalty contribution assembly
                    if (this->Has(EMBEDDED_VELOCITY))
                    {
                        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
                        for (unsigned int comp = 0; comp<TDim; comp++)
                        {
                            rRightHandSideVector(i*BlockSize+comp) += pen_coef*embedded_vel(comp);
                        }
                    }

                    // RHS residual contribution assembly
                    for (unsigned int comp = 0; comp<TDim; comp++)
                    {
                        rRightHandSideVector(i*BlockSize+comp) -= pen_coef*rData.v(i,comp);
                    }
                }
            }
        }
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
