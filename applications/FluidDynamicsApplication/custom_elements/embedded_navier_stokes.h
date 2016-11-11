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
//~ #include "utilities/enrichment_utilities_duplicate_dofs.h"      // Tetrahedra splitting
#include "utilities/split_tetrahedra_utilities.h"      // Tetrahedra splitting
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

// TODO: UPDATE THIS INFORMATION
/**this element is a 3D stokes element, stabilized by employing an ASGS stabilization
* formulation is described in the file: 
*    https://drive.google.com/file/d/0B_gRLnSH5vCwZ2Zxd09YUmlPZ28/view?usp=sharing
* symbolic implementation is defined in the file: 
*    https://drive.google.com/file/d/0B_gRLnSH5vCwaXRKRUpDbmx4VXM/view?usp=sharing
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
    
    typedef typename BaseType::element_data                     ElementDataType;
    
    typedef typename BaseType::VectorType                            VectorType;
    
    typedef typename BaseType::MatrixType                            MatrixType;
    
    typedef typename BaseType::IndexType                              IndexType;
    
    typedef typename BaseType::GeometryType::Pointer        GeometryPointerType;
    
    typedef typename BaseType::NodesArrayType                    NodesArrayType;
    
    typedef typename BaseType::PropertiesType::Pointer    PropertiesPointerType;
    
    
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

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, Element::PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return boost::make_shared< EmbeddedNavierStokes < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }
    
    
    Element::Pointer Create(IndexType NewId, Element::GeometryType::Pointer pGeom, Element::PropertiesType::Pointer pProperties) const
    {
        return boost::make_shared< EmbeddedNavierStokes < TDim, TNumNodes > >(NewId, pGeom, pProperties);
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);    // Matrix size
        
        if (rLeftHandSideMatrix.size1() != MatrixSize)
        {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); // false says not to preserve existing storage!!
        }
        else if (rLeftHandSideMatrix.size2() != MatrixSize)
        {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); // false says not to preserve existing storage!!
        }

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); // false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataType data;
        
        // Getting data for the given geometry
        double Volume;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, Volume);
        
        // Compute element size
        data.h = BaseType::ComputeH(data.DN_DX, Volume); // TODO: Check if h has to be recomputed in split elements
                
        // Database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        data.bdf0 = BDFVector[0];
        data.bdf1 = BDFVector[1];
        data.bdf2 = BDFVector[2];
        
        data.delta_t = rCurrentProcessInfo[DELTA_TIME];         // Only needed if the temporal dependent term is considered in the subscales
        data.dyn_tau_coeff = rCurrentProcessInfo[DYNAMIC_TAU];  // Only needed if the temporal dependent term is considered in the subscales

        array_1d<double, TNumNodes> distances;                  // Array to store the nodal value of the distance function

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const array_1d<double,3>& body_force = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& vel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double,3>& vel_n = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            const array_1d<double,3>& vel_nn = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);
            const array_1d<double,3>& vel_mesh = this->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
            const array_1d<double,3>& vel_conv = vel - vel_mesh;

            for(unsigned int k=0; k<TDim; k++)
            {
                data.v(i,k)   = vel[k];
                data.vn(i,k)  = vel_n[k];
                data.vnn(i,k) = vel_nn[k];
                data.vmesh(i,k) = vel_mesh[k];
                data.f(i,k)   = body_force[k];
            }
                        
            data.p[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
            data.rho[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
            data.nu[i] = this->GetGeometry()[i].FastGetSolutionStepValue(VISCOSITY);
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Allocate memory needed
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_local; 
        array_1d<double, MatrixSize> rhs_local;

        // Number of positive and negative distance function values 
        unsigned int npos=0, nneg=0;
        for (unsigned int i = 0; i<TNumNodes; i++)
        {
            if(distances[i] > 0)
                npos++;
            else
                nneg++;
        }
        
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize, MatrixSize);   // LHS initialization
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);              // RHS initialization

        // Element LHS and RHS contributions computation        
        if(npos == TNumNodes) // All nodes belong to fluid domain 
        {
            ComputeElementAsFluid<MatrixSize>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo);
        }
        else if(nneg == TNumNodes) // All nodes belong to structure domain
        {
            //~ ComputeElementAsFluid<MatrixSize>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo);
        }
        else // Element intersects both fluid and structure domains
        {            
            ComputeElementAsMixed<MatrixSize>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo, distances);
        }

        KRATOS_CATCH("Error in embedded Navier-Stokes symbolic element")
    }    

    
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);    // Matrix size
        
        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataType data;
        
        // Getting data for the given geometry
        double Volume;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, Volume);

        // Compute element size
        data.h = NavierStokes<TDim, TNumNodes>::ComputeH(data.DN_DX, Volume);
        
        // Database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        data.bdf0 = BDFVector[0];
        data.bdf1 = BDFVector[1];
        data.bdf2 = BDFVector[2];
        
        data.dyn_tau_coeff = rCurrentProcessInfo[DYNAMIC_TAU] * data.bdf0; // Only, needed if the temporal dependent term is considered in the subscales
        data.delta_t = rCurrentProcessInfo[DELTA_TIME];                    // Only, needed if the temporal dependent term is considered in the subscales
        
        array_1d<double, TNumNodes> distances;   // Array to store the nodal value of the distance function

        // Data collection
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            
            const array_1d<double,3>& body_force = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& vel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double,3>& vel_n = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            const array_1d<double,3>& vel_nn = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);
            const array_1d<double,3>& vel_mesh = this->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
            const array_1d<double,3>& vel_conv = vel - vel_mesh;

            for(unsigned int k=0; k<TDim; k++)
            {
                data.v(i,k)   = vel[k];
                data.vn(i,k)  = vel_n[k];
                data.vnn(i,k) = vel_nn[k];
                data.vmesh(i,k) = vel_mesh[k];
                data.f(i,k)   = body_force[k];
            }
            
            data.p[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
            data.rho[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
            data.nu[i] = this->GetGeometry()[i].FastGetSolutionStepValue(VISCOSITY);
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }
        
        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;
        
        // Number of positive and negative distance function values 
        unsigned int npos=0, nneg=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            if(distances[i] > 0)
                npos++;
            else
                nneg++;
        }
        
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize); // RHS initialization
        
        // Element LHS and RHS contributions computation
        if(npos == TNumNodes) // All nodes belong to fluid domain 
        {            
            ComputeRHSAsFluid<MatrixSize>(rhs_local, rRightHandSideVector, Volume, data, rCurrentProcessInfo);
        }
        if(nneg == TNumNodes) // All nodes belong to structure domain
        {
            //~ ComputeRHSAsFluid<MatrixSize>(rhs_local, rRightHandSideVector, Volume, data, rCurrentProcessInfo);
        }
        else // Element intersects both fluid and structure domains
        {
            ComputeRHSAsMixed<MatrixSize>(rhs_local, rRightHandSideVector, Volume, data, rCurrentProcessInfo, distances);
        }

        KRATOS_CATCH("")        
      
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
    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
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


    template<unsigned int MatrixSize>
    void ComputeElementAsFluid(bounded_matrix<double,MatrixSize, MatrixSize>& lhs_local,
                               array_1d<double,MatrixSize>& rhs_local,
                               MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const double& Volume,
                               ElementDataType& data,
                               ProcessInfo& rCurrentProcessInfo)
    {
        // Shape functions Gauss points values
        bounded_matrix<double, TNumNodes, TNumNodes> Ncontainer; // Container with the evaluation of the 4 shape functions in the 4 Gauss pts.
        BaseType::GetShapeFunctionsOnGauss(Ncontainer);
        
        // Loop on gauss point        
        for(unsigned int igauss = 0; igauss<Ncontainer.size2(); igauss++)
        {
            noalias(data.N) = row(Ncontainer, igauss);
            
            BaseType::ComputeConstitutiveResponse(data, rCurrentProcessInfo);
            
            BaseType::ComputeGaussPointLHSContribution(lhs_local, data);
            BaseType::ComputeGaussPointRHSContribution(rhs_local, data);
            
            // All the Gauss pts. have the same weight so the accumulated contributions can be multiplied by volume/n_nodes at the end
            noalias(rLeftHandSideMatrix) += lhs_local;
            noalias(rRightHandSideVector) += rhs_local;
        }
        
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);;
    }


    template<unsigned int MatrixSize>
    void ComputeElementAsMixed(bounded_matrix<double,MatrixSize, MatrixSize>& lhs_local,
                               array_1d<double,MatrixSize>& rhs_local,
                               MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const double& Volume,
                               ElementDataType& data,
                               ProcessInfo& rCurrentProcessInfo,
                               array_1d<double, TNumNodes>& distances)
    {
        constexpr unsigned int nEdges = (TDim-1)*3;                // Edges per element
        
        MatrixType Ncontainer;
        VectorType gauss_volumes;
        VectorType signs(nEdges); //ATTENTION: this shall be initialized of size 6 in 3D and 3 in 2D
        VectorType edge_areas(nEdges);
                
        // Splitting to determine the new Gauss pts.
        //~ unsigned int ndivisions = ComputeSplitting(data, Ncontainer, gauss_volumes, DNvalues, Nenriched, signs, distances, edge_areas);
        unsigned int ndivisions = ComputeSplitting(data, Ncontainer, gauss_volumes, signs, distances, edge_areas);
                        
        if(ndivisions == 1)
        {
            // Cases exist when the element is like not subdivided due to the characteristics of the provided distance.
            // Use the distance value at the central Gauss pt. to determine if the element is wether fluid or structure.
            array_1d<double, TNumNodes> Ncenter;
            for(unsigned int i=0; i<TNumNodes; i++) Ncenter[i]=0.25;
            const double dgauss = inner_prod(distances, Ncenter);

            // Gauss pt. is FLUID
            if(dgauss > 0)
            {
                ComputeElementAsFluid<MatrixSize>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo);
            }
            // Gaus pt. is STRUCTURE
            //~ else
            //~ {
                //~ std::cout << "Null Gauss pt. contribution"  << std::endl;
                //~ ComputeElementAsFluid<MatrixSize>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, Ncontainer, rCurrentProcessInfo);
            //~ }
        }
        else
        {            
            // Loop over subdivisions
            for(unsigned int division = 0; division<ndivisions; division++)            
            {
                // Loop on gauss points
                for(unsigned int igauss = 0; igauss<4; igauss++)            
                {
                    noalias(data.N) = row(Ncontainer, division*4+igauss); // Take the new Gauss pts. shape functions values

                    const double dgauss = inner_prod(distances, data.N);  // Compute the distance on the gauss points
                    
                    // Gauss pt. is FLUID
                    if(dgauss > 0) 
                    {
                        BaseType::ComputeConstitutiveResponse(data, rCurrentProcessInfo);
                        
                        BaseType::ComputeGaussPointLHSContribution(lhs_local, data);
                        BaseType::ComputeGaussPointRHSContribution(rhs_local, data);
                    
                        noalias(rLeftHandSideMatrix) += gauss_volumes[division*4+igauss]*lhs_local;
                        noalias(rRightHandSideVector) += gauss_volumes[division*4+igauss]*rhs_local;
                    }
                    // Gaus pt. is STRUCTURE
                    //~ else 
                    //~ {
                        //~ BaseType::ComputeConstitutiveResponse(data, rCurrentProcessInfo);
                        
                        //~ BaseType::ComputeGaussPointLHSContribution(lhs_local, data);
                        //~ BaseType::ComputeGaussPointRHSContribution(rhs_local, data);
                    
                        //~ noalias(rLeftHandSideMatrix) += gauss_volumes[division*4+igauss]*lhs_local;
                        //~ noalias(rRightHandSideVector) += gauss_volumes[division*4+igauss]*rhs_local;
                        
                        // std::cout << "Null Gauss pt. contribution"  << std::endl;
                    //~ }
                }
            }
            
            AddBoundaryConditionElementContribution(rLeftHandSideMatrix, rRightHandSideVector, data, distances, edge_areas);
        }
    }

    template<unsigned int MatrixSize>
    void ComputeRHSAsFluid(array_1d<double, MatrixSize>& rhs_local,
                           VectorType& rRightHandSideVector,
                           const double& Volume,
                           ElementDataType& data,
                           ProcessInfo& rCurrentProcessInfo)
    {                
        // Shape functions Gauss points values
        bounded_matrix<double, TNumNodes, TNumNodes> Ncontainer; // Container with the evaluation of the 4 shape functions in the 4 Gauss pts.
        NavierStokes<TDim, TNumNodes>::GetShapeFunctionsOnGauss(Ncontainer);
        
        // Loop on gauss point        
        for(unsigned int igauss = 0; igauss<Ncontainer.size2(); igauss++)
        {
            noalias(data.N) = row(Ncontainer, igauss);
            
            BaseType::ComputeConstitutiveResponse(data, rCurrentProcessInfo);
            
            BaseType::ComputeGaussPointRHSContribution(rhs_local, data);
            
            // All the Gauss pts. have the same weight so the accumulated contributions can be multiplied by volume/n_nodes at the end
            noalias(rRightHandSideVector) += rhs_local;
        }
        
        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);;
    }

    
    template<unsigned int MatrixSize>
    void ComputeRHSAsMixed(array_1d<double, MatrixSize>& rhs_local,
                           VectorType& rRightHandSideVector,
                           const double& Volume,
                           ElementDataType& data,
                           ProcessInfo& rCurrentProcessInfo,
                           array_1d<double, TNumNodes>& distances)
    {
        constexpr unsigned int nEdges = (TDim-1)*3;                // Edges per element
        
        MatrixType Ncontainer;
        VectorType gauss_volumes;
        VectorType signs(nEdges); //ATTENTION: this shall be initialized of size 6
        VectorType edge_areas(nEdges);
        
        // Splitting to determine the new Gauss pts.
        //~ unsigned int ndivisions = ComputeSplitting(data, Ncontainer, gauss_volumes, DNvalues, Nenriched, signs, distances, edge_areas);
        unsigned int ndivisions = ComputeSplitting(data, Ncontainer, gauss_volumes, signs, distances, edge_areas);
        
        if(ndivisions == 1)
        {
            // Cases exist when the element is like not subdivided due to the characteristics of the provided distance.
            // Use the distance value at the central Gauss pt. to determine if the element is wether fluid or structure.
            array_1d<double, TNumNodes> Ncenter;
            for(unsigned int i=0; i<TNumNodes; i++) Ncenter[i]=0.25;
            const double dgauss = inner_prod(distances, Ncenter);
            
            // Gauss pt. is FLUID
            if(dgauss > 0)
            {
                ComputeRHSAsFluid<MatrixSize>(rhs_local, rRightHandSideVector, Volume, data, rCurrentProcessInfo);                
            }
            // Gaus pt. is STRUCTURE
            //~ else
            //~ {
                //~ std::cout << "ComputeRHSAsStructure" << std::endl;
            //~ }
        }
        else
        {
            // Loop over subdivisions
            for(unsigned int division = 0; division<ndivisions; division++)            
            {
                // Loop on gauss points
                for(unsigned int igauss = 0; igauss<4; igauss++)            
                {
                    noalias(data.N) = row(Ncontainer, division*4+igauss); // Take the new Gauss pts. shape functions values

                    const double dgauss = inner_prod(distances, data.N); // Compute the distance on the gauss points
                    
                    // Gauss pt. is FLUID
                    if(dgauss > 0) 
                    {
                        BaseType::ComputeConstitutiveResponse(data, rCurrentProcessInfo);
                        
                        BaseType::ComputeGaussPointRHSContribution(rhs_local, data);
                    
                        noalias(rRightHandSideVector) += gauss_volumes[division*4+igauss]*rhs_local; 
                    }
                    // Gaus pt. is STRUCTURE
                    //~ else 
                    //~ {
                        //~ std::cout << "Null Gauss pt. contribution"  << std::endl;
                    //~ }
                }
            }
            
            // TODO: ComputeBoundaryConditionElementContribution (SIMILAR TO AddBoundaryConditionElementContribution) 
            //~ ComputeOutsideNodesRHSContribution(rLeftHandSideMatrix, rRightHandSideVector, data, distances, edge_areas);
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

    virtual std::string Info() const
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
    
    //~ unsigned int ComputeSplitting(ElementDataType& data,
                                  //~ MatrixType& rShapeFunctionValues, 
                                  //~ VectorType& rGaussVolumes,
                                  //~ std::vector< MatrixType >& rEnrGradientsValues,
                                  //~ MatrixType& rEnrShapeFunctionValues, 
                                  //~ VectorType& rPartitionSigns,
                                  //~ array_1d<double,TNumNodes>& distances,
                                  //~ VectorType& rEdgeAreas)
    
    unsigned int ComputeSplitting(const ElementDataType& data,
                                  MatrixType& rShapeFunctionValues, 
                                  VectorType& rGaussVolumes,
                                  VectorType& rPartitionSigns,
                                  const array_1d<double,TNumNodes>& distances,
                                  VectorType& rEdgeAreas)
    {
        
        unsigned int ndivisions = 1;
        VectorType NodalDistances = distances;
        MatrixType NoSplitGradients = data.DN_DX;
        
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
        
        //~ ndivisions = DiscontinuousShapeFunctionsUtilities::CalculateDiscontinuousShapeFunctions(NodalCoords, NoSplitGradients, 
                                                                                                //~ NodalDistances, rGaussVolumes, 
                                                                                                //~ rShapeFunctionValues, rPartitionSigns, 
                                                                                                //~ rEnrGradientsValues, rEnrShapeFunctionValues, 
                                                                                                //~ rEdgeAreas);
                                                                                                        
        //~ ndivisions = EnrichmentUtilitiesDuplicateDofs_new::CalculateTetrahedraEnrichedShapeFuncions_new(NodalCoords, NoSplitGradients, 
                                                                                                        //~ NodalDistances, rGaussVolumes, 
                                                                                                        //~ rShapeFunctionValues, rPartitionSigns, 
                                                                                                        //~ rEnrGradientsValues, rEnrShapeFunctionValues, 
                                                                                                        //~ rEdgeAreas);
                                                                                                        
        ndivisions = SplitTetrahedraUtilities::CalculateSplitTetrahedraShapeFuncions(NodalCoords, NoSplitGradients, 
                                                                                     NodalDistances, rGaussVolumes, 
                                                                                     rShapeFunctionValues, rPartitionSigns, 
                                                                                     rEdgeAreas);
    
        return ndivisions;
    }
    
    
    unsigned int IntersectionGeometryUtility(const array_1d<double,TNumNodes>& distances,
                                             const VectorType& edge_areas,
                                             const ElementDataType& data,
                                             unsigned int& npos,
                                             unsigned int& nneg,
                                             std::vector<unsigned int>& int_vec,
                                             std::vector<unsigned int>& out_vec,
                                             VectorType& cut_areas,
                                             array_1d<double, TDim>& intersection_normal,
                                             MatrixType& Ncontainer_cut,
                                             MatrixType& Ncontainer_int,
                                             MatrixType& Ncontainer_out)
    {
        constexpr unsigned int nDOFs = TDim+1;                     // DOFs per node
        constexpr unsigned int nEdges = (TDim-1)*3;                // Edges per element
        
        npos = 0;
        nneg = 0;
        
        int_vec.clear();
        out_vec.clear();
        
        for (unsigned int inode = 0; inode<TNumNodes; inode++)
        {
            if(distances[inode] > 0)
            {
                int_vec.push_back(inode);
                npos++; // Number of positive distance value nodes (fluid)
            }
            else
            {
                out_vec.push_back(inode);
                nneg++; // Number of negative distance value nodes (structure)
            }
        }
                        
        // Identify the cut edges
        unsigned int cut_points = 0;                    // Number of cut edges (or cut points)
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
                
                if((distances[i]*distances[j])<0)
                {
                    cut_points++;
                    cut_edges(aux_count) = 1; // Flag 1 says that this edge is cut
                }
                
                aux_count++;
            }
        }
        
        // Store the non-zero edge areas in an auxiliar vector cut_areas
        cut_areas.resize(cut_points, false);        
        
        // Compute the shape function values at the intersection points
        Ncontainer_cut = ZeroMatrix(cut_points, TNumNodes);     // Shape function values at the interface
        Ncontainer_int = ZeroMatrix(cut_points, npos);          // Interior nodes shape function values at the interface
        Ncontainer_out = ZeroMatrix(cut_points, nneg);          // Exterior nodes shape function values at the interface
        
        // Ncontainer_cut matrix fill
        unsigned int icut = 0;
        for (unsigned int i = 0; i<nEdges; i++)
        {            
            if (cut_edges(i) == 1) 
            {                
                unsigned int i_node = i_edges(i);
                unsigned int j_node = j_edges(i);
            
                // We take advantage of the fact that the sh. funct. that do not belong to the edge vanish along it (linear elements)
                double aux = fabs(distances[i_node])/(fabs(distances[i_node])+fabs(distances[j_node])+1e-30);
            
                // The columns associated to the nodes that conform the cut edge are filled. The other two columns remain as zero
                Ncontainer_cut(icut,i_node) = 1.0-aux;
                Ncontainer_cut(icut,j_node) = aux;
                
                cut_areas(icut) = edge_areas(i);
                
                icut++;
            }
        }
        
        // Ncontainer_int and Ncontainer_out matrices fill
        for (unsigned int i = 0; i<cut_points; i++)
        {
            for (unsigned int j = 0; j<npos; j++)
            {
                Ncontainer_int(i,j) = Ncontainer_cut(i,int_vec[j]);
            }
            
            for (unsigned int j = 0; j<nneg; j++)
            {
                Ncontainer_out(i,j) = Ncontainer_cut(i,out_vec[j]);
            }
        }
        
        // Compute intersection normal (gradient of the distance function)  
        noalias(intersection_normal) = -prod(trans(data.DN_DX), distances);      
        intersection_normal /= (norm_2(intersection_normal)+1e-15);
        
        return cut_points;
        
    }
    
    
    void AddIntersectionBoundaryContribution(MatrixType& rLeftHandSideMatrix, 
                                             VectorType& rRightHandSideVector,
                                             const ElementDataType& data, 
                                             const array_1d<double, TNumNodes>& distances,
                                             const VectorType& edge_areas,
                                             const VectorType& intersection_normal,
                                             const MatrixType& Ncontainer_cut)
    {
        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);    // Matrix size
        constexpr unsigned int BlockSize = TDim+1;
        
        // Get the previous iteration stored values
        //~ VectorType Stress = data.stress;
        
        array_1d<double, TNumNodes> aux;
        array_1d<double, MatrixSize> aux_normal;
        
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);
        bounded_matrix<double, MatrixSize, TNumNodes> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, TNumNodes);
        
        unsigned int cut_points = edge_areas.size();
        
        // Compute the contribution of the boundary pressure term at the intersection
        for (unsigned int icut=0; icut<cut_points; icut++)
        {            

            aux = row(Ncontainer_cut, icut);
            aux_normal = ZeroVector(MatrixSize);
            
            double weight = edge_areas(icut);
                        
            for (unsigned int i=0; i<TNumNodes; i++)
            {
                for(unsigned int j=0; j<TDim; j++)
                {
                    aux_normal(i*BlockSize+j) = aux(i)*intersection_normal(j);
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
        noalias(rRightHandSideVector) -= prod(auxLeftHandSideMatrix, data.p);
        
        // Compute the contribution of the tangential boundary stress term at the intersection
        bounded_matrix<double, TDim, (TDim-1)*3> normal_matrix;         // Set normal matrix
        SetNormalMatrix(intersection_normal, normal_matrix);
        
        bounded_matrix<double, (TDim-1)*3, TNumNodes*TDim> B_matrix;    // Set strain matrix
        SetStrainMatrix(data, B_matrix);
                
        bounded_matrix<double, TNumNodes*TDim, TDim> w_matrix;
        bounded_matrix<double, TNumNodes*TDim, (TDim-1)*3> aux_matrix;
        
        auxLeftHandSideMatrix.resize(TNumNodes*TDim, TNumNodes*TDim, false);
        auxRightHandSideVector.resize(TNumNodes*TDim);
        auxLeftHandSideMatrix.clear();
        auxRightHandSideVector.clear();
        
        // Compute the tangential stress LHS contribution
        for (unsigned int icut=0; icut<cut_points; icut++)
        {
            aux = row(Ncontainer_cut, icut);
            double weight = edge_areas(icut);
            
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
            aux_matrix = prod(aux_matrix, data.C);
            
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
                        rLeftHandSideMatrix(i*BlockSize+row,j*BlockSize+col) += auxLeftHandSideMatrix(i*TDim+row, j*TDim+col);
                    }
                }
            }
        }
        
        // Obtain the previous iteration velocity solution
        array_1d<double, TNumNodes*TDim> prev_sol = ZeroVector(TNumNodes*TDim);
        
        for (unsigned int i = 0; i<TNumNodes; i++)
        {
            prev_sol(i*TDim) = data.v(i,0);
            prev_sol(i*TDim+1) = data.v(i,1);
            prev_sol(i*TDim+2) = data.v(i,2);
        }
        
        noalias(auxRightHandSideVector) = prod(auxLeftHandSideMatrix, prev_sol);
        
        // Assemble the RHS tangential stress contribution to the velocity rows
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int row=0; row<TDim; row++)
            {
                rRightHandSideVector(i*BlockSize+row) -= auxRightHandSideVector(i*TDim+row);
            }
        }
        
    }
    
    
    void AddBoundaryConditionElementContribution(MatrixType& rLeftHandSideMatrix, 
                                                 VectorType& rRightHandSideVector,
                                                 const ElementDataType& data, 
                                                 const array_1d<double, TNumNodes>& distances,
                                                 const VectorType& edge_areas)
    {
        
        constexpr unsigned int BlockSize = TDim+1;                 // Block size
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;   // Matrix size
               
        // Intersection geometry data computation
        unsigned int npos = 0, nneg=0;
        std::vector<unsigned int> int_vec;
        std::vector<unsigned int> out_vec;
        VectorType cut_areas;
        array_1d<double, TDim> intersection_normal;
        MatrixType Ncontainer_cut, Ncontainer_int, Ncontainer_out;
                
        unsigned int cut_points = IntersectionGeometryUtility(distances, edge_areas, data, npos, nneg, 
                                                              int_vec, out_vec, cut_areas, intersection_normal,
                                                              Ncontainer_cut, Ncontainer_int, Ncontainer_out);       
                                                              
        AddIntersectionBoundaryContribution(rLeftHandSideMatrix, rRightHandSideVector, data, 
                                            distances, cut_areas, intersection_normal, Ncontainer_cut);
        
        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        
        for (unsigned int i = 0; i<TNumNodes; i++)
        {
            prev_sol(i*BlockSize) = data.v(i,0);
            prev_sol(i*BlockSize+1) = data.v(i,1);
            prev_sol(i*BlockSize+2) = data.v(i,2);
        }
        
        // Compute the BCs imposition matrices                
        MatrixType M_gamma(nneg, nneg);              // Outside nodes matrix (Nitche contribution)
        MatrixType N_gamma(nneg, npos);              // Interior nodes matrix (Nitche contribution)
        MatrixType f_gamma(nneg, TNumNodes);         // Matrix to compute the RHS (Nitche contribution)
        MatrixType P_gamma(TNumNodes, TNumNodes);    // Penalty matrix
        
        noalias(M_gamma) = ZeroMatrix(nneg, nneg);
        noalias(N_gamma) = ZeroMatrix(nneg, npos);
        noalias(f_gamma) = ZeroMatrix(nneg, npos);
        noalias(P_gamma) = ZeroMatrix(TNumNodes, TNumNodes);
        
        VectorType aux_out(nneg);
        VectorType aux_int(npos);
        VectorType aux_cut(TNumNodes);
        
        double intersection_area = 0.0;
        
        for (unsigned int icut = 0; icut<cut_points; icut++)
        {
            double weight = cut_areas(icut);
            intersection_area += weight;
            
            aux_out = row(Ncontainer_out, icut);
            aux_int = row(Ncontainer_int, icut);
            aux_cut = row(Ncontainer_cut, icut);
            
            M_gamma += weight*outer_prod(aux_out,aux_out);
            N_gamma += weight*outer_prod(aux_out,aux_int);
            f_gamma += weight*outer_prod(aux_out,aux_cut);
            P_gamma += weight*outer_prod(aux_cut,aux_cut);
        }
        
        // ADD PENALTY CONTRIBUTION
        // Compute the penalty coefficient
        double diag_max = 0.0;
        for (unsigned int i=0; i<MatrixSize; i++)
        {
            if ((rLeftHandSideMatrix(i,i) > diag_max) && (i%BlockSize != 0.0))
            {
                diag_max = rLeftHandSideMatrix(i,i); // Maximum diagonal value (associated to velocity)
            }
        }
        
        // TODO: Think about this value. Now is K*max(LHS(i,i))*IntArea (we integrate P_gamma over the intersection area)
        double h = data.h;
        double denominator = std::max(0.0001*h*h, intersection_area);
        double pen_coef = 100.0*diag_max/denominator;
        
        // Multiply the penalty matrix by the penalty coefficient
        P_gamma *= pen_coef;
        
        // Declare auxLeftHandSideMatrix (note that firstly it contains the penalty contribution. Then is reused in the Nitche's contributions)
        MatrixType auxLeftHandSideMatrix(MatrixSize, MatrixSize);
        noalias(auxLeftHandSideMatrix) = ZeroMatrix(MatrixSize, MatrixSize);
        
        VectorType auxRightHandSideVector(MatrixSize);  
        noalias(auxRightHandSideVector) = ZeroVector(MatrixSize);
        
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
                aux_embedded_vel(i*BlockSize) = embedded_vel(0);
                aux_embedded_vel(i*BlockSize+1) = embedded_vel(1);
                aux_embedded_vel(i*BlockSize+2) = embedded_vel(2);
            }
            
            rRightHandSideVector += prod(auxLeftHandSideMatrix, aux_embedded_vel);
        }
        
        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol); // Residual contribution assembly
                
        // ADD MODIFIED NITCHE METHOD CONTRIBUTION
        noalias(auxLeftHandSideMatrix) = ZeroMatrix(MatrixSize, MatrixSize); // Recall to reinitialize the auxLeftHandSideMatrix
        
        // Set the u_out rows to zero (outside nodes used to impose the BC)        
        for (unsigned int i = 0; i<nneg; i++)
        {
            unsigned int out_node_row_id = out_vec[i];
            
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
        for (unsigned int i = 0; i<nneg; i++)
        {
            unsigned int out_node_row_id = out_vec[i];
            
            for (unsigned int j = 0; j<nneg; j++)
            {
                unsigned int out_node_col_id = out_vec[j];
                
                for (unsigned int comp = 0; comp<TDim; comp++)
                {
                    auxLeftHandSideMatrix(out_node_row_id*BlockSize+comp, out_node_col_id*BlockSize+comp) = M_gamma(i, j);
                }
            }
        }
                        
        // Interior nodes contribution assembly
        for (unsigned int i = 0; i<nneg; i++)
        {
            unsigned int out_node_row_id = out_vec[i];
            
            for (unsigned int j = 0; j<npos; j++)
            {
                unsigned int int_node_col_id = int_vec[j];
                
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
            for (unsigned int i=0; i<nneg; i++)
            {
                unsigned int out_node_row_id = out_vec[i];
                
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
    
    
    void SetNormalMatrix(const VectorType& intersection_normal, 
                         bounded_matrix<double, TDim, (TDim-1)*3>& normal_matrix)
    {
        normal_matrix.clear();
        
        if (TDim == 3)
        {
            normal_matrix(0,0) = intersection_normal(0);
            normal_matrix(0,3) = intersection_normal(1);
            normal_matrix(0,5) = intersection_normal(2);
            normal_matrix(1,1) = intersection_normal(1);
            normal_matrix(1,3) = intersection_normal(0);
            normal_matrix(1,4) = intersection_normal(2);
            normal_matrix(2,2) = intersection_normal(2);
            normal_matrix(2,4) = intersection_normal(1);
            normal_matrix(2,5) = intersection_normal(0);            
        }
        else
        {
            normal_matrix(0,0) = intersection_normal(0);
            normal_matrix(0,2) = intersection_normal(1);
            normal_matrix(1,1) = intersection_normal(1);
            normal_matrix(1,2) = intersection_normal(0);
        }
    }
    
    
    void SetStrainMatrix(const ElementDataType& data, 
                         bounded_matrix<double, (TDim-1)*3, TNumNodes*TDim>& B_matrix)
    {
        B_matrix.clear();

        if (TDim == 3)
        {
            for (unsigned int i=0; i<TNumNodes; i++) 
            {
                B_matrix(0,i*TNumNodes) = data.DN_DX(i,0);
                B_matrix(1,i*TNumNodes+1) = data.DN_DX(i,1);
                B_matrix(2,i*TNumNodes+2) = data.DN_DX(i,2);
                B_matrix(3,i*TNumNodes) = data.DN_DX(i,1);
                B_matrix(3,i*TNumNodes+1) = data.DN_DX(i,0);
                B_matrix(4,i*TNumNodes+1) = data.DN_DX(i,2);
                B_matrix(4,i*TNumNodes+2) = data.DN_DX(i,1);
                B_matrix(5,i*TNumNodes) = data.DN_DX(i,2);
                B_matrix(5,i*TNumNodes+2) = data.DN_DX(i,0);
            }
        }
        else
        {
            for (unsigned int i=0; i<TNumNodes; i++)
            {
                B_matrix(0,i*TNumNodes) = data.DN_DX(i,0);
                B_matrix(1,i*TNumNodes+1) = data.DN_DX(i,1);
                B_matrix(2,i*TNumNodes) = data.DN_DX(i,1);
                B_matrix(2,i*TNumNodes+1) = data.DN_DX(i,0);
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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer)
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


