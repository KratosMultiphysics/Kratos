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
        Matrix  N_container;                     // Container with the shape functions values in each partition Gauss points
        Vector GaussWeights;                     // Vector containing the Gauss pts. weights
        ShapeFunctionDerivativesArrayType DN_DX; // Array containing the shape functions derivatives
        unsigned int ndivisions;                 // Number of element subdivisions
        double total_volume;                     // In 2D: element area. In 3D: element volume
        double element_size;                     // Element size
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
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

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

        // Initialize LHS and RHS 
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);

        // Element LHS and RHS contributions computation
        if (npos == TNumNodes) // All nodes belong to fluid domain
        {
            CalculateLocalSystemAsFluid<MatrixSize>(rLeftHandSideMatrix, rRightHandSideVector, data, geometry_data, rCurrentProcessInfo);
        }
        else // Element intersects both fluid and structure domains
        {
            CalculateLocalSystemAsMixed<MatrixSize>(rLeftHandSideMatrix, rRightHandSideVector, data, geometry_data, rCurrentProcessInfo);
        }

        KRATOS_CATCH("Error in embedded Ausas Navier-Stokes element!")

    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

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

        // Initialize RHS
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Element LHS and RHS contributions computation
        if (npos == TNumNodes) // All nodes belong to fluid domain
        {
            CalculateRHSAsFluid<MatrixSize>(rRightHandSideVector, data, geometry_data, rCurrentProcessInfo);
        }
        else // Element intersects both fluid and structure domains
        {
            CalculateRHSAsMixed<MatrixSize>(rRightHandSideVector, data, geometry_data, rCurrentProcessInfo);
        }

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
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
        if(PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","");
        if(DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
        if(DYNAMIC_TAU.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_TAU Key is 0. Check if the application was correctly registered.","");
        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if the application was correctly registered.","");
        if(SOUND_VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"SOUND_VELOCITY Key is 0. Check if the application was correctly registered.","");

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());

            if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
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

            // Resize the splitting data according to the obtained splitting pattern
            if ((rGeometryData.GaussWeights).size() != rGeometryData.ndivisions)
            {
                (rGeometryData.GaussWeights).resize(rGeometryData.ndivisions, false);
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

                rGeometryData.GaussWeights[division] = partition_volumes[division];
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
            rGeometryData.GaussWeights.resize(ngauss, false);
            for (unsigned int igauss = 0; igauss<ngauss; ++igauss)
            {
                rGeometryData.GaussWeights[igauss] = DetJ[igauss]*rIntegrationPoints[igauss].Weight();
            }
        }
    }

    template <unsigned int MatrixSize>
    void CalculateLocalSystemAsFluid(MatrixType &rLeftHandSideMatrix,
                                     VectorType &rRightHandSideVector,
                                     ElementDataStruct &rData,
                                     ElementGeometryDataStruct &rGeometryData,
                                     ProcessInfo &rCurrentProcessInfo)
    {
        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;
        bounded_matrix<double, MatrixSize, MatrixSize> lhs_local;

        rData.h = rGeometryData.element_size;
        const unsigned int ngauss = (rGeometryData.N_container).size2();

        // Loop on gauss points
        for (unsigned int igauss = 0; igauss < ngauss; ++igauss)
        {
            // Gather in the Gauss pt. data from the geometry data structure
            const double wgauss = rGeometryData.GaussWeights[igauss];
            noalias(rData.N) = row(rGeometryData.N_container, igauss);
            noalias(rData.DN_DX) = rGeometryData.DN_DX[igauss];

            ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, rData);
            ComputeGaussPointLHSContribution(lhs_local, rData);

            noalias(rLeftHandSideMatrix) += wgauss * lhs_local;
            noalias(rRightHandSideVector) += wgauss * rhs_local;
        }
    }

    template <unsigned int MatrixSize>
    void CalculateLocalSystemAsMixed(MatrixType &rLeftHandSideMatrix,
                                     VectorType &rRightHandSideVector,
                                     ElementDataStruct &rData,
                                     ElementGeometryDataStruct &rGeometryData,
                                     ProcessInfo &rCurrentProcessInfo)
    {
        // Allocate memory needed
        bounded_matrix<double, MatrixSize, MatrixSize> lhs_local;
        array_1d<double, MatrixSize> rhs_local;

        rData.h = rGeometryData.element_size;
        const unsigned int ngauss = (rGeometryData.N_container).size2();

        // Loop on gauss points
        for (unsigned int igauss = 0; igauss < ngauss; ++igauss)
        {
            // Gather in the Gauss pt. data from the geometry data structure
            const double wgauss = rGeometryData.GaussWeights[igauss];
            noalias(rData.N) = row(rGeometryData.N_container, igauss);
            noalias(rData.DN_DX) = rGeometryData.DN_DX[igauss];

            ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, rData);
            ComputeGaussPointLHSContribution(lhs_local, rData);

            noalias(rLeftHandSideMatrix) += wgauss * lhs_local;
            noalias(rRightHandSideVector) += wgauss * rhs_local;
        }

        // TODO: ADD THE PENALTY TERMS IN HERE
    }

    template <unsigned int MatrixSize>
    void CalculateRHSAsFluid(VectorType &rRightHandSideVector,
                             ElementDataStruct &rData,
                             ElementGeometryDataStruct &rGeometryData,
                             ProcessInfo &rCurrentProcessInfo)
    {
        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;
        bounded_matrix<double, MatrixSize, MatrixSize> lhs_local;

        rData.h = rGeometryData.element_size;
        const unsigned int ngauss = (rGeometryData.N_container).size2();

        // Loop on gauss points
        for (unsigned int igauss = 0; igauss < ngauss; ++igauss)
        {
            // Gather in the Gauss pt. data from the geometry data structure
            const double wgauss = rGeometryData.GaussWeights[igauss];
            noalias(rData.N) = row(rGeometryData.N_container, igauss);
            noalias(rData.DN_DX) = rGeometryData.DN_DX[igauss];

            ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, rData);

            noalias(rRightHandSideVector) += wgauss * rhs_local;
        }
    }

    template <unsigned int MatrixSize>
    void CalculateRHSAsMixed(VectorType &rRightHandSideVector,
                             ElementDataStruct &rData,
                             ElementGeometryDataStruct &rGeometryData,
                             ProcessInfo &rCurrentProcessInfo)
    {
        // Allocate memory needed
        bounded_matrix<double, MatrixSize, MatrixSize> lhs_local;
        array_1d<double, MatrixSize> rhs_local;

        rData.h = rGeometryData.element_size;
        const unsigned int ngauss = (rGeometryData.N_container).size2();

        // Loop on gauss points
        for (unsigned int igauss = 0; igauss < ngauss; ++igauss)
        {
            // Gather in the Gauss pt. data from the geometry data structure
            const double wgauss = rGeometryData.GaussWeights[igauss];
            noalias(rData.N) = row(rGeometryData.N_container, igauss);
            noalias(rData.DN_DX) = rGeometryData.DN_DX[igauss];

            ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, rData);

            noalias(rRightHandSideVector) += wgauss * rhs_local;
        }

        // TODO: ADD THE PENALTY TERMS IN HERE
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

    // 3D tetrahedra shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(bounded_matrix<double,4,4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
        Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
        Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
        Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    }

    // 2D triangle shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(bounded_matrix<double,3,3>& Ncontainer)
    {
        const double one_sixt = 1.0/6.0;
        const double two_third = 2.0/3.0;
        Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
        Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
        Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
    }

    // // 3D tetrahedra shape functions values at centered Gauss point
    // void GetShapeFunctionsOnUniqueGauss(boost::numeric::ublas::bounded_matrix<double,1,4>& Ncontainer)
    // {
    //     Ncontainer(0,0) = 0.25; Ncontainer(0,1) = 0.25; Ncontainer(0,2) = 0.25; Ncontainer(0,3) = 0.25;
    // }

    // // 2D triangle shape functions values at centered Gauss point
    // void GetShapeFunctionsOnUniqueGauss(boost::numeric::ublas::bounded_matrix<double,1,3>& Ncontainer)
    // {
    //     Ncontainer(0,0) = 1.0/3.0; Ncontainer(0,1) = 1.0/3.0; Ncontainer(0,2) = 1.0/3.0;
    // }

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
            rData.C.resize(strain_size,strain_size,false);
        if(rData.stress.size() != strain_size)
            rData.stress.resize(strain_size,false);
        if(rData.strain.size() != strain_size)
            rData.strain.resize(strain_size,false);

        ComputeStrain(rData, strain_size);

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


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                                    Fluid2DASGS& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                                    const Fluid2DASGS& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES_ELEMENT_INCLUDED  defined
