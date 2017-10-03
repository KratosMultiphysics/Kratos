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

#if !defined(KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES)
#define  KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "utilities/discont_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

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
class EmbeddedAusasNavierStokes : public Element
{
public:
    ///@name Type Definitions
    ///@{
    typedef GeometryType::IntegrationPointsArrayType                InteGrationPointsType;
    typedef GeometryType::ShapeFunctionsGradientsType   ShapeFunctionDerivativesArrayType;

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedAusasNavierStokes);

    struct ElementDataStruct
    {
        bounded_matrix<double, TNumNodes, TDim> v, vn, vnn, vmesh, f;
        array_1d<double,TNumNodes> p, pn, pnn, rho, mu;

        bounded_matrix<double, TNumNodes, TDim> DN_DX;
        array_1d<double, TNumNodes> N;

        Matrix C;             // Matrix to store the constitutive matrix in Voigt notation
        Vector stress;        // Vector to store the stress values in Voigt notation
        Vector strain;        // Vector to store the stain values in Voigt notation

        double bdf0;          // BDF2 scheme coefficient 0
        double bdf1;          // BDF2 scheme coefficient 1
        double bdf2;          // BDF2 scheme coefficient 2
        double c;             // Wave velocity (needed if artificial compressibility is considered)
        double h;             // Element size
        double dt;            // Time increment
        double dyn_tau;       // Dynamic tau considered in ASGS stabilization coefficients
    };

    struct ElementGeometryDataStruct
    {
        array_1d<double, TDim> intersection_normal;   // Intersection unit normal vector

        Matrix N_container;                           // Container with the shape functions values in each partition Gauss points
        Matrix N_positive_cut;                        // Container with the shape functions values in the edges intersection pts (positive distance side).
        Matrix N_negative_cut;                        // Container with the shape functions values in the edges intersection pts (negative distance side).
        Vector gauss_weights;                         // Vector containing the Gauss pts. weights
        Vector cut_edge_areas;                        // Vector containing the cut edges intersection areas

        ShapeFunctionDerivativesArrayType DN_DX;      // Array containing the shape functions derivatives

        unsigned int ndivisions;                      // Number of element subdivisions
        unsigned int ncutpoints;                      // Number of element intersected edges

        double total_volume;                          // In 2D: element area. In 3D: element volume
        double element_size;                          // Element size
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    EmbeddedAusasNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    EmbeddedAusasNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~EmbeddedAusasNavierStokes() override {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< EmbeddedAusasNavierStokes < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< EmbeddedAusasNavierStokes < TDim, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
        {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!
        }
        else if (rLeftHandSideMatrix.size2() != MatrixSize)
        {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!
        }

        if (rRightHandSideVector.size() != MatrixSize)
        {
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!
        }

        // Initialize LHS and RHS
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);

        // Set the elemental distance vector
        Vector& elemental_distances = this->GetValue(ELEMENTAL_DISTANCES);
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            elemental_distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Number of positive and negative distance function values
        unsigned int npos = 0, nneg = 0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            if (elemental_distances[i] > 0.0)
            {
                npos++;
            }
            else
            {
                nneg++;
            }
        }

        if (npos != 0 && nneg != 0)
        {
            this->Set(TO_SPLIT, true);
        }

        // Fill the data structure to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);
        ElementGeometryDataStruct geometry_data;
        this->FillElementGeometryData(geometry_data);

        // Element LHS and RHS contributions computation
        CalculateLocalSystemContribution(rLeftHandSideMatrix, rRightHandSideVector, data, geometry_data, rCurrentProcessInfo);

        KRATOS_CATCH("Error in embedded Ausas Navier-Stokes element CalculateLocalSystem!")

    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rRightHandSideVector.size() != MatrixSize)
        {
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!
        }

        // Initialize RHS
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Set the elemental distance vector
        Vector &elemental_distances = this->GetValue(ELEMENTAL_DISTANCES);
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            elemental_distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Number of positive and negative distance function values
        unsigned int npos = 0, nneg = 0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            if (elemental_distances[i] > 0.0)
            {
                npos++;
            }
            else
            {
                nneg++;
            }
        }

        if (npos != 0 && nneg != 0)
        {
            this->Set(TO_SPLIT, true);
        }

        // Fill the data structure to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);
        ElementGeometryDataStruct geometry_data;
        this->FillElementGeometryData(geometry_data);

        // Element LHS and RHS contributions computation
        CalculateRightHandSideContribution(rRightHandSideVector, data, geometry_data, rCurrentProcessInfo);

        KRATOS_CATCH("");
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
        KRATOS_TRY;

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        // Check that all required variables have been registered
        if(VELOCITY.Key() == 0)
            KRATOS_ERROR << "VELOCITY Key is 0. Check if the application was correctly registered.";
        if(PRESSURE.Key() == 0)
            KRATOS_ERROR << "PRESSURE Key is 0. Check if the application was correctly registered.";
        if(DENSITY.Key() == 0)
            KRATOS_ERROR << "DENSITY Key is 0. Check if the application was correctly registered.";
        if(DYNAMIC_TAU.Key() == 0)
            KRATOS_ERROR << "DYNAMIC_TAU Key is 0. Check if the application was correctly registered.";
        if(DELTA_TIME.Key() == 0)
            KRATOS_ERROR << "DELTA_TIME Key is 0. Check if the application was correctly registered.";
        if(SOUND_VELOCITY.Key() == 0)
            KRATOS_ERROR << "SOUND_VELOCITY Key is 0. Check if the application was correctly registered.";

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR << "Missing VELOCITY variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_ERROR << "Missing PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_ERROR << "Missing VELOCITY component degree of freedom on node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_ERROR << "Missing PRESSURE component degree of freedom on node " << this->GetGeometry()[i].Id();
        }

        // Check constitutive law
        if(mpConstitutiveLaw == nullptr)
            KRATOS_ERROR << "The constitutive law was not set. Cannot proceed. Call the navier_stokes.h Initialize() method needs to be called.";

        mpConstitutiveLaw->Check(GetProperties(), this->GetGeometry(), rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("");
    }


    void Calculate(const Variable<double>& rVariable,
                   double& rOutput,
                   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_CATCH("")
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
        return "EmbeddedAusasNavierStokes";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const override

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static member variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    // Constitutive law pointer
    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    // Symbolic function implementing the element
    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override;
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs, const ElementDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ElementDataStruct& data);

    ///@}
    ///@name Protected Operators
    ///@{

    EmbeddedAusasNavierStokes() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

    // Element initialization (constitutive law)
    void Initialize() override
    {
        KRATOS_TRY;

        // Initalize the constitutive law pointer
        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( GetProperties(), this->GetGeometry(), row( this->GetGeometry().ShapeFunctionsValues(), 0 ) );

        // Initialize the ELEMENTAL_DISTANCES variable (make it threadsafe)
        Vector zero_vector(TNumNodes, 0.0);
        this->SetValue(ELEMENTAL_DISTANCES, zero_vector);

        KRATOS_CATCH("");
    }

    // Auxiliar function to fill the element data structure
    void FillElementData(ElementDataStruct& rData, const ProcessInfo& rCurrentProcessInfo)
    {
        // Database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        rData.bdf0 = BDFVector[0];
        rData.bdf1 = BDFVector[1];
        rData.bdf2 = BDFVector[2];

        rData.dyn_tau = rCurrentProcessInfo[DYNAMIC_TAU];   // Only, needed if the temporal dependent term is considered in the subscales
        rData.dt = rCurrentProcessInfo[DELTA_TIME];         // Only, needed if the temporal dependent term is considered in the subscales

        rData.c = rCurrentProcessInfo[SOUND_VELOCITY];      // Wave velocity

        for (unsigned int i = 0; i < TNumNodes; i++)
        {

            const array_1d<double,3>& body_force = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& vel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double,3>& vel_n = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            const array_1d<double,3>& vel_nn = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);
            const array_1d<double,3>& vel_mesh = this->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);

            for(unsigned int k=0; k<TDim; k++)
            {
                rData.v(i,k)   = vel[k];
                rData.vn(i,k)  = vel_n[k];
                rData.vnn(i,k) = vel_nn[k];
                rData.vmesh(i,k) = vel_mesh[k];
                rData.f(i,k)   = body_force[k];
            }

            rData.p[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
            rData.pn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,1);
            rData.pnn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,2);
            rData.rho[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
            rData.mu[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
        }

    }

    // Auxiliar function to fill the element splitting data
    void FillElementGeometryData(ElementGeometryDataStruct &rGeometryData)
    {
        constexpr unsigned int MaxPartitions = 3*(TDim-1);

        GeometryType& rGeom = this->GetGeometry();

        // Getting data for the given geometry
        array_1d<double, TNumNodes> N_continuous;
        bounded_matrix<double, TNumNodes, TDim> DN_DX_continuous;
        GeometryUtils::CalculateGeometryData(rGeom, DN_DX_continuous, N_continuous, rGeometryData.total_volume);

        // Compute element size
        rGeometryData.element_size = ComputeH(DN_DX_continuous);

        if (this->Is(TO_SPLIT))
        {
            // Compute the intersection normal
            this->ComputeIntersectionNormal(rGeometryData, DN_DX_continuous);

            // Arrays initialization
            bounded_matrix<double, TNumNodes, TDim> nodal_coords;
            bounded_matrix<double, MaxPartitions, TNumNodes> shape_function_values;
            bounded_matrix<double, MaxPartitions, TNumNodes> enriched_shape_function_values;

            array_1d<double, TNumNodes> nodal_distances = this->GetValue(ELEMENTAL_DISTANCES);

            array_1d<double, MaxPartitions> partition_volumes;
            array_1d<double, MaxPartitions> partition_signs;
            array_1d<double, MaxPartitions> edge_areas;

            std::vector<Matrix> gauss_gradients(MaxPartitions);

            // Fill the nodal coordinates matrix
            this->FillNodalCoordinatesMatrix(nodal_coords);

            rGeometryData.ndivisions = DiscontinuousShapeFunctionsUtilities::CalculateDiscontinuousShapeFunctions(nodal_coords,
                                                                                                                  DN_DX_continuous,
                                                                                                                  nodal_distances,
                                                                                                                  partition_volumes,
                                                                                                                  shape_function_values,
                                                                                                                  partition_signs,
                                                                                                                  gauss_gradients,
                                                                                                                  enriched_shape_function_values,
                                                                                                                  edge_areas);

            // Compute the shape function values at the intersection points in the edges
            this->ComputeIntersectionPointsShapeFunctionValues(rGeometryData, edge_areas);

            // Resize the splitting data according to the obtained splitting pattern
            if ((rGeometryData.gauss_weights).size() != rGeometryData.ndivisions)
            {
                (rGeometryData.gauss_weights).resize(rGeometryData.ndivisions, false);
            }

            if ((rGeometryData.N_container).size1() != rGeometryData.ndivisions || (rGeometryData.N_container).size2() != TDim+1)
            {
                (rGeometryData.N_container).resize(rGeometryData.ndivisions, TDim+1, false);
            }

            if ((rGeometryData.DN_DX).size() != rGeometryData.ndivisions)
            {
                (rGeometryData.DN_DX).resize(rGeometryData.ndivisions);
            }

            // Gather the computed splitting data to the splitting data structure
            for (unsigned int division=0; division<rGeometryData.ndivisions; ++division)
            {
                for (unsigned int dof=0; dof<TDim+1; ++dof)
                {
                    rGeometryData.N_container(division, dof) = enriched_shape_function_values(division, dof);
                }

                rGeometryData.gauss_weights[division] = partition_volumes[division];
                rGeometryData.DN_DX[division] = gauss_gradients[division];
            }
        }
        else
        {
            // Fill the shape functions container
            rGeometryData.N_container = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

            // Fill the shape functions gradient container
            Vector DetJ;
            rGeom.ShapeFunctionsIntegrationPointsGradients(rGeometryData.DN_DX, DetJ, GeometryData::GI_GAUSS_2);

            // Fill the Gauss pts. weights container
            const unsigned int ngauss = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
            const InteGrationPointsType& rIntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
            rGeometryData.gauss_weights.resize(ngauss, false);
            for (unsigned int igauss = 0; igauss<ngauss; ++igauss)
            {
                rGeometryData.gauss_weights[igauss] = DetJ[igauss]*rIntegrationPoints[igauss].Weight();
            }
        }
    }

    void ComputeIntersectionPointsShapeFunctionValues(ElementGeometryDataStruct& rGeometryData,
                                                      array_1d<double, 3*(TDim-1)>& rEdgeAreas)
    {
        constexpr unsigned int nedges = (TDim-1)*3; // Edges per element

        // Identify the cut edges
        rGeometryData.ncutpoints = 0;             // Number of cut edges (or cut points)
        array_1d<unsigned int, nedges> i_edges;   // i-node of the edges vector
        array_1d<unsigned int, nedges> j_edges;   // j-node of the edges vector
        array_1d<unsigned int, nedges> cut_edges; // 0: non cut edge 1: cut edge

        const array_1d<double, TNumNodes>& nodal_distances = this->GetValue(ELEMENTAL_DISTANCES);

        unsigned int aux_count = 0;
        for (unsigned int i=0; i<TNumNodes; ++i)
        {
            for (unsigned int j=i+1; j<TNumNodes; ++j)
            {
                i_edges(aux_count) = i; // Construct the edges i-nodes vector
                j_edges(aux_count) = j; // Construct the edges j-nodes vector
                cut_edges(aux_count) = 0;

                if ((nodal_distances[i] * nodal_distances[j]) < 0)
                {
                    rGeometryData.ncutpoints++;
                    cut_edges(aux_count) = 1; // Flag 1 says that this edge is cut
                }

                aux_count++;
            }
        }

        if ((rGeometryData.cut_edge_areas).size() != rGeometryData.ncutpoints)
        {
            (rGeometryData.cut_edge_areas).resize(rGeometryData.ncutpoints, false);
        }

        // Compute the enriched shape function values at the edges intersection pts.
        rGeometryData.N_positive_cut = ZeroMatrix(rGeometryData.ncutpoints, TNumNodes);
        rGeometryData.N_negative_cut = ZeroMatrix(rGeometryData.ncutpoints, TNumNodes);

        unsigned int icut = 0;
        for (unsigned int iedge=0; iedge<nedges; ++iedge)
        {
            // Check if the edge is cut
            if (cut_edges(iedge) == 1)
            {
                // Get the nodes that compose the cut edge
                unsigned int inode = i_edges(iedge);
                unsigned int jnode = j_edges(iedge);

                // Get the intersected point edge area
                rGeometryData.cut_edge_areas(icut) = rEdgeAreas(iedge);

                // First edge node shape functions intersection values
                if (nodal_distances[inode] > 0.0)
                {
                    rGeometryData.N_positive_cut(icut, inode) = 1.0; // Fill the positive distance side shape functions array
                }
                else
                {
                    rGeometryData.N_negative_cut(icut, inode) = 1.0; // Fill the negative distance side shape functions array
                }

                // Second edge node shape functions intersection values
                if (nodal_distances[jnode] > 0.0)
                {
                    rGeometryData.N_positive_cut(icut, jnode) = 1.0; // Fill the positive distance side shape functions array
                }
                else
                {
                    rGeometryData.N_negative_cut(icut, jnode) = 1.0; // Fill the negative distance side shape functions array
                }

                icut++;
            }
        }

    }

    void CalculateLocalSystemContribution(MatrixType &rLeftHandSideMatrix,
                                          VectorType &rRightHandSideVector,
                                          ElementDataStruct &rData,
                                          ElementGeometryDataStruct &rGeometryData,
                                          ProcessInfo &rCurrentProcessInfo)
    {
        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;
        bounded_matrix<double, MatrixSize, MatrixSize> lhs_local;

        rData.h = rGeometryData.element_size;
        const unsigned int ngauss = (rGeometryData.N_container).size2();

        // Loop on gauss points
        for (unsigned int igauss = 0; igauss < ngauss; ++igauss)
        {
            // Gather in the Gauss pt. data from the geometry data structure
            const double wgauss = rGeometryData.gauss_weights[igauss];
            noalias(rData.N) = row(rGeometryData.N_container, igauss);
            noalias(rData.DN_DX) = rGeometryData.DN_DX[igauss];

            ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, rData);
            ComputeGaussPointLHSContribution(lhs_local, rData);

            noalias(rLeftHandSideMatrix) += wgauss * lhs_local;
            noalias(rRightHandSideVector) += wgauss * rhs_local;
        }

        if (this->Is(TO_SPLIT))
        {
            // Add the intersection boundary fluxes contribution comping from the integration by parts
            this->AddSystemPressureBoundaryTermsContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rGeometryData);

            // Add the normal component penalty contribution
            // this->AddSystemNormalVelocityPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rGeometryData);

            // Use the pressure as a Lagrange multiplier to enforce the no penetration condition
            this->AddSystemNormalVelocityLagrangeMultiplierContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rGeometryData);
        }
    }

    void CalculateRightHandSideContribution(VectorType &rRightHandSideVector,
                                            ElementDataStruct &rData,
                                            ElementGeometryDataStruct &rGeometryData,
                                            ProcessInfo &rCurrentProcessInfo)
    {
        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;

        rData.h = rGeometryData.element_size;
        const unsigned int ngauss = (rGeometryData.N_container).size2();

        // Loop on gauss points
        for (unsigned int igauss = 0; igauss < ngauss; ++igauss)
        {
            // Gather in the Gauss pt. data from the geometry data structure
            const double wgauss = rGeometryData.gauss_weights[igauss];
            noalias(rData.N) = row(rGeometryData.N_container, igauss);
            noalias(rData.DN_DX) = rGeometryData.DN_DX[igauss];

            ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, rData);

            noalias(rRightHandSideVector) += wgauss * rhs_local;
        }

        if (this->Is(TO_SPLIT))
        {
            // Add the intersection boundary fluxes contribution comping from the integration by parts
            this->AddRHSPressureBoundaryTermsContribution(rRightHandSideVector, rData, rGeometryData);

            // Add the normal component penalty contribution
            this->AddRHSNormalVelocityPenaltyContribution(rRightHandSideVector, rData, rGeometryData);
        }
    }

    /**
    * This function adds the local system contribution of the interface pressure boundary terms,
    * coming from the integration by parts.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rGeometryData: reference to the intersection data structure
    */
    void AddSystemPressureBoundaryTermsContribution(MatrixType &rLeftHandSideMatrix,
                                                    VectorType &rRightHandSideVector,
                                                    const ElementDataStruct &rData,
                                                    const ElementGeometryDataStruct &rGeometryData)
    {
        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Declare auxiliar arrays
        bounded_matrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Compute the LHS and RHS boundary terms contributions
        array_1d<double, TDim> side_normal; 
        array_1d<double, TNumNodes> aux_cut;

        // Contribution coming fron the shear stress operator
        // At the momoent the contribution coming from the stress operator is neglected since it is relatively low and implies the computation of the sh. function gradients in the intersection pts.

        // Contribution coming from the positive side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_positive_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = -1.0 * rGeometryData.intersection_normal;

            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        auxLeftHandSideMatrix(i * BlockSize + d, j * BlockSize + TDim) += weight * aux_cut(i) * aux_cut(j) * side_normal(d);
                    }
                }
            }
        }

        // Contribution coming from the negative side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_negative_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = rGeometryData.intersection_normal;

            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        auxLeftHandSideMatrix(i * BlockSize + d, j * BlockSize + TDim) += weight * aux_cut(i) * aux_cut(j) * side_normal(d);
                    }
                }
            }
        }

        // LHS assembly
        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS assembly
        rRightHandSideVector -= prod(auxLeftHandSideMatrix, prev_sol);
    }

    /**
    * This function adds the RHS contribution of the interface boundary terms,
    * coming from the integration by parts.
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rGeometryData: reference to the intersection data structure
    */
    void AddRHSPressureBoundaryTermsContribution(VectorType &rRightHandSideVector,
                                                 const ElementDataStruct &rData,
                                                 const ElementGeometryDataStruct &rGeometryData)
    {
        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Declare auxiliar arrays
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);

        // Compute the LHS and RHS boundary terms contributions
        array_1d<double, TDim> side_normal;
        array_1d<double, TNumNodes> aux_cut;

        // Contribution coming fron the shear stress operator
        // At the momoent the contribution coming from the stress operator is neglected since it is relatively low and implies the computation of the sh. function gradients in the intersection pts.

        // Contribution coming from the positive side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_positive_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = -1.0 * rGeometryData.intersection_normal;

            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        auxRightHandSideVector(i * BlockSize + d) += weight * aux_cut(i) * aux_cut(j) * side_normal(d) * prev_sol(j * BlockSize + TDim);
                    }
                }
            }
        }

        // Contribution coming from the negative side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_negative_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = rGeometryData.intersection_normal;

            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        auxRightHandSideVector(i * BlockSize + d) += weight * aux_cut(i) * aux_cut(j) * side_normal(d) * prev_sol(j * BlockSize + TDim);
                    }
                }
            }
        }


        // RHS assembly
        rRightHandSideVector -= auxRightHandSideVector;
    }

    /**
    * This function adds the local system contribution of the penalty no penetration imposition.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rGeometryData: reference to the intersection data structure
    */
    void AddSystemNormalVelocityPenaltyContribution(MatrixType &rLeftHandSideMatrix,
                                                    VectorType &rRightHandSideVector,
                                                    const ElementDataStruct &rData,
                                                    const ElementGeometryDataStruct &rGeometryData)
    {
        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        array_1d<double, MatrixSize> solution_jump = ZeroVector(MatrixSize);

        // Obtain the previous iteration velocity solution
        GetPreviousSolutionVector(rData, prev_sol);

        // Compute the velocity diference to penalize
        if (this->Has(EMBEDDED_VELOCITY))
        {
            const array_1d<double, 3> &embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            solution_jump = aux_embedded_vel - prev_sol;
        }
        else
        {
            solution_jump = - prev_sol;
        }

        // Compute the penalty coefficient (TODO: Implement a K depending on the element size)
        const double pen_coef = ComputePenaltyCoefficient(rLeftHandSideMatrix, rData, rGeometryData);

        // Compute the LHS and RHS penalty contributions
        array_1d<double, TDim> side_normal;
        array_1d<double, TNumNodes> aux_cut;
        bounded_matrix<double, MatrixSize, MatrixSize> P_gamma = ZeroMatrix(MatrixSize, MatrixSize);

        // Contribution coming from the positive side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_positive_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = -1.0 * rGeometryData.intersection_normal;

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int m = 0; m < TDim; ++m)
                    {
                        const unsigned int row = i * BlockSize + m;
                        for (unsigned int n = 0; n < TDim; ++n)
                        {
                            const unsigned int col = j * BlockSize + n;
                            P_gamma(row, col) += pen_coef * weight * aux_cut(i) * side_normal(m) * side_normal(n) * aux_cut(j);
                        }
                    }
                }
            }
        }

        // Contribution coming from the negative side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_negative_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = rGeometryData.intersection_normal;

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int m = 0; m < TDim; ++m)
                    {
                        const unsigned int row = i * BlockSize + m;
                        for (unsigned int n = 0; n < TDim; ++n)
                        {
                            const unsigned int col = j * BlockSize + n;
                            P_gamma(row, col) += pen_coef * weight * aux_cut(i) * side_normal(m) * side_normal(n) * aux_cut(j);
                        }
                    }
                }
            }
        }

        // LHS assembly
        rLeftHandSideMatrix += P_gamma;

        // RHS assembly
        rRightHandSideVector += prod(P_gamma, solution_jump);

    }

    /**
    * This function adds the local system contribution of the no penetration imposition,
    * by means of the pressure acting as a Lagrange multiplier.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rGeometryData: reference to the intersection data structure
    */
    void AddSystemNormalVelocityLagrangeMultiplierContribution(MatrixType &rLeftHandSideMatrix,
                                                               VectorType &rRightHandSideVector,
                                                               const ElementDataStruct &rData,
                                                               const ElementGeometryDataStruct &rGeometryData)
    {
        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        array_1d<double, MatrixSize> solution_jump = ZeroVector(MatrixSize);

        // Obtain the previous iteration velocity solution
        GetPreviousSolutionVector(rData, prev_sol);

        // Compute the velocity diference to penalize
        if (this->Has(EMBEDDED_VELOCITY))
        {
            const array_1d<double, 3> &embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            solution_jump = prev_sol - aux_embedded_vel;
        }
        else
        {
            solution_jump = prev_sol;
        }

        // Compute the LHS and RHS penalty contributions
        array_1d<double, TDim> side_normal;
        array_1d<double, TNumNodes> aux_cut;
        bounded_matrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Drop the pressure rows and columns
        for (unsigned int i=0; i<TNumNodes; ++i)
        {   
            const unsigned int aux = i*BlockSize + TDim;
            for (unsigned int j=0; j<MatrixSize; ++j)
            {
                rLeftHandSideMatrix(aux, j) = 0.0;
                rLeftHandSideMatrix(j, aux) = 0.0;
            }
        }

        // Contribution coming from the positive side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_positive_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = -1.0 * rGeometryData.intersection_normal;

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                const unsigned int row = i * BlockSize + TDim;

                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int m = 0; m < TDim; ++m)
                    {
                        const unsigned int col = j * BlockSize + m;
                        const double aux_value = weight * aux_cut(i) * side_normal(m) * aux_cut(j);
                        auxLeftHandSideMatrix(row, col) += aux_value;
                        auxLeftHandSideMatrix(col, row) += aux_value;
                    }
                }
            }
        }

        // Contribution coming from the negative side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_negative_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = rGeometryData.intersection_normal;

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                const unsigned int row = i * BlockSize + TDim;

                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int m = 0; m < TDim; ++m)
                    {
                        const unsigned int col = j * BlockSize + m;
                        const double aux_value = weight * aux_cut(i) * side_normal(m) * aux_cut(j);
                        auxLeftHandSideMatrix(row, col) += aux_value;
                        auxLeftHandSideMatrix(col, row) += aux_value;
                    }
                }
            }
        }

        // LHS assembly
        rLeftHandSideMatrix += auxLeftHandSideMatrix;

        // RHS assembly
        rRightHandSideVector -= prod(auxLeftHandSideMatrix, solution_jump);
    }

    /**
    * This function adds the RHS contribution of the penalty no penetration imposition.
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    * @param rGeometryData: reference to the intersection data structure
    */
    void AddRHSNormalVelocityPenaltyContribution(VectorType &rRightHandSideVector,
                                                 const ElementDataStruct &rData,
                                                 const ElementGeometryDataStruct &rGeometryData)
    {
        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        array_1d<double, MatrixSize> solution_jump = ZeroVector(MatrixSize);

        // Obtain the previous iteration velocity solution
        GetPreviousSolutionVector(rData, prev_sol);

        // Compute the velocity diference to penalize
        if (this->Has(EMBEDDED_VELOCITY))
        {
            const array_1d<double, 3> &embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i)
            {
                for (unsigned int comp=0; comp<TDim; ++comp)
                {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            solution_jump = aux_embedded_vel - prev_sol;
        }
        else
        {
            solution_jump = - prev_sol;
        }

        // Compute the penalty coefficient (TODO: Implement a K independent of the LHS)
        const double pen_coef = 100.0;

        // Compute the RHS penalty contributions
        array_1d<double, TDim> side_normal;
        array_1d<double, TNumNodes> aux_cut;
        array_1d<double, MatrixSize> P_gamma_RHS = ZeroVector(MatrixSize); 

        // Contribution coming from the positive side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_positive_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = -1.0 * rGeometryData.intersection_normal;

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int m = 0; m < TDim; ++m)
                    {
                        const unsigned int row = i * BlockSize + m;
                        for (unsigned int n = 0; n < TDim; ++n)
                        {
                            P_gamma_RHS(row) += pen_coef * weight * aux_cut(i) * side_normal(m) * side_normal(n) * aux_cut(j) * solution_jump(row);
                        }
                    }
                }
            }
        }

        // Contribution coming from the negative side pressure term
        for (unsigned int icut = 0; icut < rGeometryData.ncutpoints; icut++) // Consider the Gauss points as the edge intersection points
        {
            const double weight = rGeometryData.cut_edge_areas(icut);

            // Get the shape functions according to the positive or negative distance sides
            aux_cut = row(rGeometryData.N_negative_cut, icut);

            // Get the normal according to the positive or negative distance sides
            side_normal = rGeometryData.intersection_normal;

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    for (unsigned int m = 0; m < TDim; ++m)
                    {
                        const unsigned int row = i * BlockSize + m;
                        for (unsigned int n = 0; n < TDim; ++n)
                        {
                            P_gamma_RHS(row) += pen_coef * weight * aux_cut(i) * side_normal(m) * side_normal(n) * aux_cut(j) * solution_jump(row);
                        }
                    }
                }
            }
        }

        // RHS assembly
        rRightHandSideVector += P_gamma_RHS;

    }

    // Auxiliar function to fill the tetrahedra nodal coordinates
    void FillNodalCoordinatesMatrix(bounded_matrix<double,4,3>& rNodalCoordinates)
    {
        const GeometryType& rGeom = this->GetGeometry();

        for (unsigned int inode=0; inode<TNumNodes; ++inode)
        {
            rNodalCoordinates(inode,0) = rGeom[inode].X();
            rNodalCoordinates(inode,1) = rGeom[inode].Y();
            rNodalCoordinates(inode,2) = rGeom[inode].Z();
        }
    }

    // Auxiliar function to fill the triangle nodal coordinates
    void FillNodalCoordinatesMatrix(bounded_matrix<double,3,2>& rNodalCoordinates)
    {
        const GeometryType& rGeom = this->GetGeometry();

        for (unsigned int inode=0; inode<TNumNodes; ++inode)
        {
            rNodalCoordinates(inode,0) = rGeom[inode].X();
            rNodalCoordinates(inode,1) = rGeom[inode].Y();
        }
    }

    // Auxiliar function to compute the element size
    double ComputeH(bounded_matrix<double,TNumNodes, TDim>& DN_DX)
    {
        double h=0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }
        h = sqrt(h)/static_cast<double>(TNumNodes);
        return h;
    }

    /**
    * This function computes the penalty coefficient for the level set BC imposition
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rData: reference to element data structure
    * @param rGeometryData: reference to the intersection data structure
    */
    double ComputePenaltyCoefficient(MatrixType &rLeftHandSideMatrix,
                                     const ElementDataStruct &rData,
                                     const ElementGeometryDataStruct &rGeometryData)
    {
        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        // Compute the penalty coefficient as K*max(LHS(i,i))*IntArea (we integrate P_gamma over the intersection area)
        double diag_max = 0.0;
        for (unsigned int i = 0; i < MatrixSize; i++)
        {
            if ((std::abs(rLeftHandSideMatrix(i, i)) > diag_max) && ((i+1) % BlockSize != 0.0))
            {
                diag_max = std::abs(rLeftHandSideMatrix(i, i)); // Maximum diagonal value (associated to velocity)
            }
        }

        // Compute the intersection area
        double intersection_area = 0.0;
        for (unsigned int i = 0; i < rGeometryData.ncutpoints; ++i)
        {
            intersection_area += rGeometryData.cut_edge_areas(i);
        }

        // Return the penalty coefficient
        const double K = 1000.0;
        // const double denominator = std::max(0.001 * rData.h * rData.h, intersection_area);
        // const double pen_coef = K * diag_max / denominator;
        const double pen_coef = K * diag_max / intersection_area;

        return pen_coef;
    }

    /**
    * This functions sets the auxiliar matrix to compute the tangential projection in Voigt notation
    * @param rData: reference to the element data structure
    * @param rPrevSolVector: reference to the previous solution vector
    */
    void GetPreviousSolutionVector(const ElementDataStruct &rData,
                                   array_1d<double, TNumNodes *(TDim + 1)> &rPrevSolVector)
    {
        rPrevSolVector.clear();

        for (unsigned int i=0; i<TNumNodes; ++i)
        {
            for (unsigned int comp=0; comp<TDim; ++comp)
            {
                rPrevSolVector(i*(TDim+1)+comp) = rData.v(i,comp);
            }
            rPrevSolVector(i*(TDim+1)+TDim) = rData.p(i);
        }
    }

    void ComputeIntersectionNormal(ElementGeometryDataStruct& rGeometryData,
                                   const bounded_matrix<double, TNumNodes, TDim>& rDN_DX)
    {
        array_1d<double, TNumNodes> distances = this->GetValue(ELEMENTAL_DISTANCES);
        rGeometryData.intersection_normal = prod(trans(rDN_DX), distances);
        rGeometryData.intersection_normal /= norm_2(rGeometryData.intersection_normal);
    }

    // Computes the strain rate in Voigt notation
    void ComputeStrain(ElementDataStruct& rData, const unsigned int& strain_size)
    {
        const bounded_matrix<double, TNumNodes, TDim>& v = rData.v;
        const bounded_matrix<double, TNumNodes, TDim>& DN = rData.DN_DX;

        // Compute strain (B*v)
        // 3D strain computation
        if (strain_size == 6)
        {
            rData.strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
            rData.strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
            rData.strain[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
            rData.strain[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
            rData.strain[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
            rData.strain[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
        }
        // 2D strain computation
        else if (strain_size == 3)
        {
            rData.strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
            rData.strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
            rData.strain[2] = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
        }
    }

    // Call the constitutive law to get the stress value
    virtual void ComputeConstitutiveResponse(ElementDataStruct& rData, const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int strain_size = (TDim*3)-3;

        if(rData.C.size1() != strain_size)
        {
            rData.C.resize(strain_size, strain_size, false);
        }
        else if(rData.C.size2() != strain_size)
        {
            rData.C.resize(strain_size, strain_size, false);
        }

        if(rData.stress.size() != strain_size)
        {
            rData.stress.resize(strain_size,false);
        }

        if(rData.strain.size() != strain_size)
        {
            rData.strain.resize(strain_size,false);
        }

        this->ComputeStrain(rData, strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(this->GetGeometry(), GetProperties(), rCurrentProcessInfo);

        const Vector Nvec(rData.N);
        Values.SetShapeFunctionsValues(Nvec);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        Values.SetStrainVector(rData.strain);       //this is the input parameter
        Values.SetStressVector(rData.stress);       //this is an ouput parameter
        Values.SetConstitutiveMatrix(rData.C);      //this is an ouput parameter

        //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
        //this is ok under the hypothesis that no history dependent behaviour is employed
        mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    }

    virtual double ComputeEffectiveViscosity(const ElementDataStruct& rData)
    {
        // Computes the effective viscosity as the average of the lower diagonal constitutive tensor
        double mu_eff = 0.0;
        const unsigned int strain_size = (TDim-1)*3;
        for (unsigned int i=TDim; i<strain_size; ++i){mu_eff += rData.C(i,i);}
        mu_eff /= (strain_size - TDim);

        return mu_eff;
    }

    /**
    * This functions sets the auxiliar matrix to compute the tangential projection in Voigt notation
    * @param rSplittingData: reference to the intersection data structure
    * @param rVoigtNormProjMatrix: reference to the computed tangential projection auxiliar matrix
    */
    void SetVoigtNormalProjectionMatrix(const ElementGeometryDataStruct &rGeometryData,
                                        bounded_matrix<double, TDim, (TDim - 1) * 3> &rVoigtNormProjMatrix)
    {
        rVoigtNormProjMatrix.clear();

        if (TDim == 3)
        {
            // Fill the normal projection matrix for Voigt notation
            rVoigtNormProjMatrix(0, 0) = rGeometryData.intersection_normal(0);
            rVoigtNormProjMatrix(0, 3) = rGeometryData.intersection_normal(1);
            rVoigtNormProjMatrix(0, 5) = rGeometryData.intersection_normal(2);
            rVoigtNormProjMatrix(1, 1) = rGeometryData.intersection_normal(1);
            rVoigtNormProjMatrix(1, 3) = rGeometryData.intersection_normal(0);
            rVoigtNormProjMatrix(1, 4) = rGeometryData.intersection_normal(2);
            rVoigtNormProjMatrix(2, 2) = rGeometryData.intersection_normal(2);
            rVoigtNormProjMatrix(2, 4) = rGeometryData.intersection_normal(1);
            rVoigtNormProjMatrix(2, 5) = rGeometryData.intersection_normal(0);
        }
        else
        {
            // Fill the noromal projection matrix for Voigt notation
            rVoigtNormProjMatrix(0, 0) = rGeometryData.intersection_normal(0);
            rVoigtNormProjMatrix(0, 2) = rGeometryData.intersection_normal(1);
            rVoigtNormProjMatrix(1, 1) = rGeometryData.intersection_normal(1);
            rVoigtNormProjMatrix(1, 2) = rGeometryData.intersection_normal(0);
        }
    }

    /**
    * This functions sets the B strain matrix with zero value in the pressure rows
    * @param rData: reference to element data structure (it contains the shape functions derivatives)
    * @param rB_matrix: reference to the computed B strain matrix
    */
    void SetExpandedStrainMatrix(const bounded_matrix<double, TNumNodes, TDim> &rDN_DX,
                                 bounded_matrix<double, (TDim - 1) * 3, TNumNodes * TDim> &rB_matrix)
    {
        constexpr unsigned int BlockSize = TDim + 1;

        rB_matrix.clear();

        // Set the shape function derivatives values
        if (TDim == 3)
        {
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                rB_matrix(0, i * BlockSize)     = rDN_DX(i, 0);
                rB_matrix(1, i * BlockSize + 1) = rDN_DX(i, 1);
                rB_matrix(2, i * BlockSize + 2) = rDN_DX(i, 2);
                rB_matrix(3, i * BlockSize)     = rDN_DX(i, 1);
                rB_matrix(3, i * BlockSize + 1) = rDN_DX(i, 0);
                rB_matrix(4, i * BlockSize + 1) = rDN_DX(i, 2);
                rB_matrix(4, i * BlockSize + 2) = rDN_DX(i, 1);
                rB_matrix(5, i * BlockSize)     = rDN_DX(i, 2);
                rB_matrix(5, i * BlockSize + 2) = rDN_DX(i, 0);    
            }
        }
        else
        {
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                rB_matrix(0, i * BlockSize)     = rDN_DX(i, 0);
                rB_matrix(1, i * BlockSize + 1) = rDN_DX(i, 1);
                rB_matrix(2, i * BlockSize)     = rDN_DX(i, 1);
                rB_matrix(2, i * BlockSize + 1) = rDN_DX(i, 0);
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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
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

#endif // KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES_ELEMENT_INCLUDED  defined
