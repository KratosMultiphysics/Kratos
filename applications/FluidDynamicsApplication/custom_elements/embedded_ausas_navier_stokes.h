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
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

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

    typedef GeometryType::Pointer                                     GeometryPointerType;

    typedef GeometryType::IntegrationPointsArrayType                InteGrationPointsType;

    typedef GeometryType::ShapeFunctionsGradientsType         ShapeFunctionsGradientsType;

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedAusasNavierStokes);

    struct EmbeddedAusasElementDataStruct {

        bounded_matrix<double, TNumNodes, TDim> v, vn, vnn, vmesh, f;
        array_1d<double,TNumNodes> p, pn, pnn, rho, mu;

        array_1d<double, TNumNodes>                 N;        // Shape functions values on Gauss pt. container
        bounded_matrix<double, TNumNodes, TDim>     DN_DX;    // Shape functions gradients values on Gauss pt. container

        Matrix C;             // Matrix to store the constitutive matrix in Voigt notation
        Vector stress;        // Vector to store the stress values in Voigt notation
        Vector strain;        // Vector to store the stain values in Voigt notation

        double bdf0;          // BDF2 scheme coefficient 0
        double bdf1;          // BDF2 scheme coefficient 1
        double bdf2;          // BDF2 scheme coefficient 2
        double c;             // Wave velocity (needed if artificial compressibility is considered)
        double h;             // Element size
        double volume;        // In 2D: element area. In 3D: element volume
        double dt;            // Time increment
        double dyn_tau;       // Dynamic tau considered in ASGS stabilization coefficients

        // No splitted elements geometry data containers
        VectorType                      w_gauss;       // No splitted element Gauss pts. weights values
        MatrixType                      N_gauss;       // No splitted element shape function values on Gauss pts.
        ShapeFunctionsGradientsType     DN_DX_gauss;   // No splitted element shape function gradients values on Gauss pts.

        // Splitted element geometry data containers
        // Positive side geometry data
        MatrixType                  N_pos_side;             // Positive distance element side shape functions values
        std::vector<MatrixType>     DN_DX_pos_side;         // Positive distance element side shape functions gradients values
        VectorType                  w_gauss_pos_side;       // Positive distance element side Gauss pts. weights

        // Negative side geometry data
        MatrixType                  N_neg_side;             // Negative distance element side shape functions values
        std::vector<MatrixType>     DN_DX_neg_side;         // Negative distance element side shape functions gradients values
        VectorType                  w_gauss_neg_side;       // Negative distance element side Gauss pts. weights

        // Positive interface geometry data
        MatrixType                  N_pos_int;              // Positive interface Gauss pts. shape functions values
        std::vector<MatrixType>     DN_DX_pos_int;          // Positive interface Gauss pts. shape functions gradients values
        VectorType                  w_gauss_pos_int;        // Positive interface Gauss pts. weights
        std::vector<VectorType>     pos_int_unit_normals;   // Positive interface unit normal vector in each Gauss pt.

        // Negative interface geometry data
        MatrixType                  N_neg_int;              // Positive interface Gauss pts. shape functions values
        std::vector<MatrixType>     DN_DX_neg_int;          // Positive interface Gauss pts. shape functions gradients values
        VectorType                  w_gauss_neg_int;        // Positive interface Gauss pts. weights
        std::vector<VectorType>     neg_int_unit_normals;   // Positive interface unit normal vector in each Gauss pt.

        std::vector<unsigned int>   int_vec_identifiers;    // Interior (fluid) nodes identifiers
        std::vector<unsigned int>   out_vec_identifiers;    // Outside (stucture) nodes identifiers

        unsigned int    n_pos = 0;                          // Number of postivie distance nodes
        unsigned int    n_neg = 0;                          // Number of negative distance nodes

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

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override {
        KRATOS_TRY
        return boost::make_shared< EmbeddedAusasNavierStokes < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override {
        KRATOS_TRY
        return boost::make_shared< EmbeddedAusasNavierStokes < TDim, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override {

        KRATOS_TRY;

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize) {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false);
        } else if (rLeftHandSideMatrix.size2() != MatrixSize) {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false);
        }

        if (rRightHandSideVector.size() != MatrixSize) {
            rRightHandSideVector.resize(MatrixSize, false);
        }

        // Initialize LHS and RHS
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);

        // Set the elemental distances vector variable
        Vector& elemental_distances = this->GetValue(ELEMENTAL_DISTANCES);
        for (unsigned int i = 0; i < TNumNodes; i++) {
            elemental_distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Decides if the element is split and fills data structure accordingly
        EmbeddedAusasElementDataStruct data;
        this->FillEmbeddedAusasElementData(data, rCurrentProcessInfo);

        // Element LHS and RHS contributions computation
        CalculateLocalSystemContribution(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo);

        KRATOS_CATCH("Error in embedded Ausas Navier-Stokes element CalculateLocalSystem!")

    }

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override {

        KRATOS_TRY;

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rRightHandSideVector.size() != MatrixSize) {
            rRightHandSideVector.resize(MatrixSize, false);
        }

        // Initialize RHS
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Set the elemental distances vector variable
        Vector& elemental_distances = this->GetValue(ELEMENTAL_DISTANCES);
        for (unsigned int i = 0; i < TNumNodes; i++) {
            elemental_distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Decides if the element is split and fills data structure accordingly
        EmbeddedAusasElementDataStruct data;
        this->FillEmbeddedAusasElementData(data, rCurrentProcessInfo);

        // Element LHS and RHS contributions computation
        CalculateRightHandSideContribution(rRightHandSideVector, data, rCurrentProcessInfo);

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
    int Check(const ProcessInfo& rCurrentProcessInfo) override {
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
        for(unsigned int i_node = 0; i_node < this->GetGeometry().size(); ++i_node) {
            if(this->GetGeometry()[i_node].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR << "Missing VELOCITY variable on solution step data for node " << this->GetGeometry()[i_node].Id();
            if(this->GetGeometry()[i_node].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_ERROR << "Missing PRESSURE variable on solution step data for node " << this->GetGeometry()[i_node].Id();
            if(this->GetGeometry()[i_node].HasDofFor(VELOCITY_X) == false ||
               this->GetGeometry()[i_node].HasDofFor(VELOCITY_Y) == false ||
               this->GetGeometry()[i_node].HasDofFor(VELOCITY_Z) == false)
                KRATOS_ERROR << "Missing VELOCITY component degree of freedom on node " << this->GetGeometry()[i_node].Id();
            if(this->GetGeometry()[i_node].HasDofFor(PRESSURE) == false)
                KRATOS_ERROR << "Missing PRESSURE component degree of freedom on node " << this->GetGeometry()[i_node].Id();
        }

        // Check constitutive law
        if(mpConstitutiveLaw == nullptr)
            KRATOS_ERROR << "The constitutive law was not set. Cannot proceed. Call the navier_stokes.h Initialize() method needs to be called.";

        mpConstitutiveLaw->Check(GetProperties(), this->GetGeometry(), rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("");
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
    std::string Info() const override {
        return "EmbeddedAusasNavierStokes";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
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

    void ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs, const EmbeddedAusasElementDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const EmbeddedAusasElementDataStruct& data);

    ///@}
    ///@name Protected Operators
    ///@{

    EmbeddedAusasNavierStokes() : Element()
    {}

    ///@}
    ///@name Protected Operations
    ///@{

    // Element initialization (constitutive law)
    void Initialize() override {
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
    void FillEmbeddedAusasElementData(
        EmbeddedAusasElementDataStruct& rData,
        const ProcessInfo& rCurrentProcessInfo) {

        // Get the element geometry
        GeometryType& r_geom = this->GetGeometry();

        // Compute characteristic parent element size
        rData.h = ComputeH();

        // Database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        rData.bdf0 = BDFVector[0];
        rData.bdf1 = BDFVector[1];
        rData.bdf2 = BDFVector[2];

        rData.dyn_tau = rCurrentProcessInfo[DYNAMIC_TAU];   // Only, needed if the temporal dependent term is considered in the subscales
        rData.dt = rCurrentProcessInfo[DELTA_TIME];         // Only, needed if the temporal dependent term is considered in the subscales

        rData.c = rCurrentProcessInfo[SOUND_VELOCITY];      // Wave velocity

        for (unsigned int i = 0; i < TNumNodes; i++) {

            const array_1d<double,3>& body_force = r_geom[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& vel = r_geom[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double,3>& vel_n = r_geom[i].FastGetSolutionStepValue(VELOCITY,1);
            const array_1d<double,3>& vel_nn = r_geom[i].FastGetSolutionStepValue(VELOCITY,2);
            const array_1d<double,3>& vel_mesh = r_geom[i].FastGetSolutionStepValue(MESH_VELOCITY);

            for(unsigned int k=0; k<TDim; k++) {
                rData.v(i,k)   = vel[k];
                rData.vn(i,k)  = vel_n[k];
                rData.vnn(i,k) = vel_nn[k];
                rData.vmesh(i,k) = vel_mesh[k];
                rData.f(i,k)   = body_force[k];
            }

            rData.p[i] = r_geom[i].FastGetSolutionStepValue(PRESSURE);
            rData.pn[i] = r_geom[i].FastGetSolutionStepValue(PRESSURE,1);
            rData.pnn[i] = r_geom[i].FastGetSolutionStepValue(PRESSURE,2);
            rData.rho[i] = r_geom[i].FastGetSolutionStepValue(DENSITY);
            rData.mu[i] = r_geom[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
        }

        // Getting the nodal distances vector
        const array_1d<double, TNumNodes>& distances = this->GetValue(ELEMENTAL_DISTANCES);

        (rData.int_vec_identifiers).clear();
        (rData.out_vec_identifiers).clear();

        // Number of positive and negative distance function values
        for (unsigned int inode = 0; inode<TNumNodes; inode++) {
            if(distances[inode] > 0.0) {
                rData.n_pos++;
                (rData.int_vec_identifiers).push_back(inode);
            } else {
                rData.n_neg++;
                (rData.out_vec_identifiers).push_back(inode);
            }
        }

        if (rData.n_pos != 0 && rData.n_neg != 0) {
            this->Set(TO_SPLIT, true);
        }

        // If the element is split, get the modified shape functions
        if (this->Is(TO_SPLIT)) {

            GeometryPointerType p_geom = this->pGetGeometry();

            // Construct the modified shape fucntions utility
            ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = nullptr;
            if (TNumNodes == 4) {
                p_ausas_modified_sh_func = boost::make_shared<Tetrahedra3D4AusasModifiedShapeFunctions>(p_geom, distances);
            } else {
                p_ausas_modified_sh_func = boost::make_shared<Triangle2D3AusasModifiedShapeFunctions>(p_geom, distances);
            }

            // Call the positive side modified shape functions calculator
            p_ausas_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
                rData.N_pos_side,
                rData.DN_DX_pos_side,
                rData.w_gauss_pos_side,
                GeometryData::GI_GAUSS_2);

            // Call the negative side modified shape functions calculator
            p_ausas_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                rData.N_neg_side,
                rData.DN_DX_neg_side,
                rData.w_gauss_neg_side,
                GeometryData::GI_GAUSS_2);

            // Call the positive side interface modified shape functions calculator
            p_ausas_modified_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                rData.N_pos_int,
                rData.DN_DX_pos_int,
                rData.w_gauss_pos_int,
                GeometryData::GI_GAUSS_2);

            // Call the negative side interface modified shape functions calculator
            p_ausas_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                rData.N_neg_int,
                rData.DN_DX_neg_int,
                rData.w_gauss_neg_int,
                GeometryData::GI_GAUSS_2);

            // Call the positive side Gauss pts. unit normal calculator
            p_ausas_modified_sh_func->ComputePositiveSideInterfaceAreaNormals(
                rData.pos_int_unit_normals,
                GeometryData::GI_GAUSS_2);

            // Call the negative side Gauss pts. unit normal calculator
            p_ausas_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
                rData.neg_int_unit_normals,
                GeometryData::GI_GAUSS_2);

            // Normalize the obtained positive and negative sides area normals
            const double tol = std::pow(1e-3*rData.h, TDim-1); // Tolerance to avoid the unit normal to blow up
            const unsigned int n_gauss_pos = (rData.pos_int_unit_normals).size();
            const unsigned int n_gauss_neg = (rData.neg_int_unit_normals).size();

            for (unsigned int i_gauss = 0;  i_gauss < n_gauss_pos; ++i_gauss) {
                Vector& normal = rData.pos_int_unit_normals[i_gauss];
                const double n_norm = norm_2(normal);
                normal /= std::max(n_norm, tol);
            }

            for (unsigned int i_gauss = 0;  i_gauss < n_gauss_neg; ++i_gauss) {
                Vector& normal = rData.neg_int_unit_normals[i_gauss];
                const double n_norm = norm_2(normal);
                normal /= std::max(n_norm, tol);
            }

        } else {
            // Fill the shape functions container
            rData.N_gauss = r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

            // Fill the shape functions gradient container
            Vector det_jacobian;
            r_geom.ShapeFunctionsIntegrationPointsGradients(rData.DN_DX_gauss, det_jacobian, GeometryData::GI_GAUSS_2);

            // Fill the Gauss pts. weights container
            const unsigned int n_gauss = r_geom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
            const InteGrationPointsType& rIntegrationPoints = r_geom.IntegrationPoints(GeometryData::GI_GAUSS_2);
            rData.w_gauss.resize(n_gauss, false);
            for (unsigned int i_gauss = 0; i_gauss<n_gauss; ++i_gauss) {
                rData.w_gauss[i_gauss] = det_jacobian[i_gauss] * rIntegrationPoints[i_gauss].Weight();
            }
        }
    }

    void CalculateLocalSystemContribution(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        EmbeddedAusasElementDataStruct &rData,
        ProcessInfo &rCurrentProcessInfo) {

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;
        bounded_matrix<double, MatrixSize, MatrixSize> lhs_local;

        // Decide if the element is wether split or not and add the contribution accordingly
        if (this->Is(TO_SPLIT)) {

            // Add the positive side volume contribution
            const unsigned int n_pos_gauss = (rData.w_gauss_pos_side).size();
            for (unsigned int i_pos_gauss = 0; i_pos_gauss < n_pos_gauss; ++i_pos_gauss) {
                // Gather Ausas Gauss pt. positive side values from data structure
                noalias(rData.N) = row(rData.N_pos_side, i_pos_gauss);
                noalias(rData.DN_DX) = rData.DN_DX_pos_side[i_pos_gauss];
                const double w_gauss = rData.w_gauss_pos_side[i_pos_gauss];

                ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

                ComputeGaussPointRHSContribution(rhs_local, rData);
                ComputeGaussPointLHSContribution(lhs_local, rData);

                noalias(rLeftHandSideMatrix) += w_gauss * lhs_local;
                noalias(rRightHandSideVector) += w_gauss * rhs_local;
            }

            // Add the negative side volume contribution
            const unsigned int n_neg_gauss = (rData.w_gauss_neg_side).size();
            for (unsigned int i_neg_gauss = 0; i_neg_gauss < n_neg_gauss; ++i_neg_gauss) {
                // Gather Ausas Gauss pt. negative side values from data structure
                noalias(rData.N) = row(rData.N_neg_side, i_neg_gauss);
                noalias(rData.DN_DX) = rData.DN_DX_neg_side[i_neg_gauss];
                const double w_gauss = rData.w_gauss_neg_side[i_neg_gauss];

                ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

                ComputeGaussPointRHSContribution(rhs_local, rData);
                ComputeGaussPointLHSContribution(lhs_local, rData);

                noalias(rLeftHandSideMatrix) += w_gauss * lhs_local;
                noalias(rRightHandSideVector) += w_gauss * rhs_local;
            }

            // Add the intersection boundary fluxes contribution comping from the integration by parts
            this->AddSystemBoundaryTermsContribution(rLeftHandSideMatrix, rRightHandSideVector, rData);

            // Add the normal component penalty contribution
            this->AddSystemNormalVelocityPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData);

            // Use the pressure as a Lagrange multiplier to enforce the no penetration condition
            // this->AddSystemNormalVelocityLagrangeMultiplierContribution(rLeftHandSideMatrix, rRightHandSideVector, rData);

        } else {

            // If the element is not splitted, add the standard Navier-Stokes contribution
            const unsigned int n_gauss = (rData.w_gauss).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                // Gather the standard shape functions Gauss pt. values from data structure
                const double w_gauss = rData.w_gauss[i_gauss];
                noalias(rData.N) = row(rData.N_gauss, i_gauss);
                noalias(rData.DN_DX) = rData.DN_DX_gauss[i_gauss];

                ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

                ComputeGaussPointRHSContribution(rhs_local, rData);
                ComputeGaussPointLHSContribution(lhs_local, rData);

                noalias(rLeftHandSideMatrix) += w_gauss * lhs_local;
                noalias(rRightHandSideVector) += w_gauss * rhs_local;
            }
        }
    }

    void CalculateRightHandSideContribution(
        VectorType &rRightHandSideVector,
        EmbeddedAusasElementDataStruct &rData,
        ProcessInfo &rCurrentProcessInfo) {

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;

        // Decide if the element is wether split or not and add the contribution accordingly
        if (this->Is(TO_SPLIT)) {

            // Add the positive side volume contribution
            const unsigned int n_pos_gauss = (rData.w_gauss_pos_side).size();
            for (unsigned int i_pos_gauss = 0; i_pos_gauss < n_pos_gauss; ++i_pos_gauss) {
                // Gather Ausas Gauss pt. positive side values from data structure
                noalias(rData.N) = row(rData.N_pos_side, i_pos_gauss);
                noalias(rData.DN_DX) = rData.DN_DX_pos_side[i_pos_gauss];
                const double w_gauss = rData.w_gauss_pos_side[i_pos_gauss];

                ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

                ComputeGaussPointRHSContribution(rhs_local, rData);

                noalias(rRightHandSideVector) += w_gauss * rhs_local;
            }

            // Add the negative side volume contribution
            const unsigned int n_neg_gauss = (rData.w_gauss_neg_side).size();
            for (unsigned int i_neg_gauss = 0; i_neg_gauss < n_neg_gauss; ++i_neg_gauss) {
                // Gather Ausas Gauss pt. negative side values from data structure
                noalias(rData.N) = row(rData.N_neg_side, i_neg_gauss);
                noalias(rData.DN_DX) = rData.DN_DX_neg_side[i_neg_gauss];
                const double w_gauss = rData.w_gauss_neg_side[i_neg_gauss];

                ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

                ComputeGaussPointRHSContribution(rhs_local, rData);

                noalias(rRightHandSideVector) += w_gauss * rhs_local;
            }

            // Add the intersection boundary fluxes contribution comping from the integration by parts
            this->AddRHSBoundaryTermsContribution(rRightHandSideVector, rData);

            // Add the normal component penalty contribution
            this->AddRHSNormalVelocityPenaltyContribution(rRightHandSideVector, rData);

        } else {

            // If the element is not splitted, add the standard Navier-Stokes contribution
            const unsigned int n_gauss = (rData.w_gauss).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                // Gather the standard shape functions Gauss pt. values from data structure
                const double w_gauss = rData.w_gauss[i_gauss];
                noalias(rData.N) = row(rData.N_gauss, i_gauss);
                noalias(rData.DN_DX) = rData.DN_DX_gauss[i_gauss];

                ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

                ComputeGaussPointRHSContribution(rhs_local, rData);

                noalias(rRightHandSideVector) += w_gauss * rhs_local;
            }
        }
    }

    /**
    * This function adds the local system contribution of the interface pressure boundary terms,
    * coming from the integration by parts.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    */
    void AddSystemBoundaryTermsContribution(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const EmbeddedAusasElementDataStruct &rData) {

        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Declare auxiliar arrays
        bounded_matrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Contribution coming from the positive side boundary term
        const unsigned int n_int_pos_gauss = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss = 0; i_gauss < n_int_pos_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_pos_int(i_gauss);

            // Get the positive side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss);

            // Get the positive side shape functions gradients Gauss pt. values
            const bounded_matrix<double, TNumNodes, TDim> aux_DN_DX = rData.DN_DX_pos_int[i_gauss];

            // Get the positive side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.pos_int_unit_normals[i_gauss];

            // Set the auxiliar test function matrix
            bounded_matrix<double, MatrixSize, TDim> test_mat = ZeroMatrix(MatrixSize, TDim);
            this->SetTestMatrix(aux_N, test_mat);

            // Set the normal projection matrix in Voigt notation
            bounded_matrix<double, TDim, (TDim-1)*3> normal_proj_mat = ZeroMatrix(TDim, (TDim-1)*3);
            this->SetVoigtNormalProjectionMatrix(side_normal, normal_proj_mat);

            // Set the strain matrix (expanded with zeros in the pressure rows)
            bounded_matrix<double, (TDim-1)*3, MatrixSize> exp_strain_mat = ZeroMatrix((TDim-1)*3, MatrixSize);
            this->SetExpandedStrainMatrix(aux_DN_DX, exp_strain_mat);

            // Pressure contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int d = 0; d < TDim; ++d) {
                        auxLeftHandSideMatrix(i * BlockSize + d, j * BlockSize + TDim) += w_gauss * aux_N(i) * aux_N(j) * side_normal(d);
                    }
                }
            }

            // Shear stress contribution
            for (unsigned int i = 0; i < MatrixSize; ++i) {
                for (unsigned int j = 0; j < MatrixSize; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        for (unsigned int k = 0; k < (TDim-1)*3; ++k) {
                            for (unsigned int n = 0; n < (TDim-1)*3; ++n) {
                                auxLeftHandSideMatrix(i,j) -= w_gauss * test_mat(i,m) * normal_proj_mat(m,n) * rData.C(n,k) * exp_strain_mat(k,j);
                            }
                        }
                    }
                }
            }
        }

        // Contribution coming from the negative side boundary term
        const unsigned int n_int_neg_gauss = (rData.w_gauss_neg_int).size();

        for (unsigned int i_gauss = 0; i_gauss < n_int_neg_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_neg_int(i_gauss);

            // Get the negative side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_neg_int, i_gauss);

            // Get the negative side shape functions gradients Gauss pt. values
            const bounded_matrix<double, TNumNodes, TDim> aux_DN_DX = rData.DN_DX_neg_int[i_gauss];

            // Get the negative side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.neg_int_unit_normals[i_gauss];

            // Set the auxiliar test function matrix
            bounded_matrix<double, MatrixSize, TDim> test_mat = ZeroMatrix(MatrixSize, TDim);
            this->SetTestMatrix(aux_N, test_mat);

            // Set the normal projection matrix in Voigt notation
            bounded_matrix<double, TDim, (TDim-1)*3> normal_proj_mat = ZeroMatrix(TDim, (TDim-1)*3);
            this->SetVoigtNormalProjectionMatrix(side_normal, normal_proj_mat);

            // Set the strain matrix (expanded with zeros in the pressure rows)
            bounded_matrix<double, (TDim-1)*3, MatrixSize> exp_strain_mat = ZeroMatrix((TDim-1)*3, MatrixSize);
            this->SetExpandedStrainMatrix(aux_DN_DX, exp_strain_mat);

            // Pressure contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int d = 0; d < TDim; ++d) {
                        auxLeftHandSideMatrix(i * BlockSize + d, j * BlockSize + TDim) += w_gauss * aux_N(i) * aux_N(j) * side_normal(d);
                    }
                }
            }

            // Shear stress contribution
            for (unsigned int i = 0; i < MatrixSize; ++i) {
                for (unsigned int j = 0; j < MatrixSize; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        for (unsigned int k = 0; k < (TDim-1)*3; ++k) {
                            for (unsigned int n = 0; n < (TDim-1)*3; ++n) {
                                auxLeftHandSideMatrix(i,j) -= w_gauss * test_mat(i,m) * normal_proj_mat(m,n) * rData.C(n,k) * exp_strain_mat(k,j);
                            }
                        }
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
    */
    void AddRHSBoundaryTermsContribution(
        VectorType &rRightHandSideVector,
        const EmbeddedAusasElementDataStruct &rData) {

        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Declare auxiliar arrays
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);

        // Contribution coming from the positive side boundary term
        const unsigned int n_int_pos_gauss = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss = 0; i_gauss < n_int_pos_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_pos_int(i_gauss);

            // Get the positive side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss);

            // Get the positive side shape functions gradients Gauss pt. values
            const bounded_matrix<double, TNumNodes, TDim> aux_DN_DX = rData.DN_DX_pos_int[i_gauss];

            // Get the positive side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.pos_int_unit_normals[i_gauss];

            // Set the auxiliar test function matrix
            bounded_matrix<double, MatrixSize, TDim> test_mat = ZeroMatrix(MatrixSize, TDim);
            this->SetTestMatrix(aux_N, test_mat);

            // Set the normal projection matrix in Voigt notation
            bounded_matrix<double, TDim, (TDim-1)*3> normal_proj_mat = ZeroMatrix(TDim, (TDim-1)*3);
            this->SetVoigtNormalProjectionMatrix(side_normal, normal_proj_mat);

            // Set the strain matrix (expanded with zeros in the pressure rows)
            bounded_matrix<double, (TDim-1)*3, MatrixSize> exp_strain_mat = ZeroMatrix((TDim-1)*3, MatrixSize);
            this->SetExpandedStrainMatrix(aux_DN_DX, exp_strain_mat);

            // Pressure contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int d = 0; d < TDim; ++d) {
                        auxRightHandSideVector(i * BlockSize + d) += w_gauss * aux_N(i) * aux_N(j) * side_normal(d) * prev_sol(j * BlockSize + TDim);
                    }
                }
            }

            // Shear stress contribution
            for (unsigned int i = 0; i < MatrixSize; ++i) {
                for (unsigned int j = 0; j < MatrixSize; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        for (unsigned int k = 0; k < (TDim-1)*3; ++k) {
                            for (unsigned int n = 0; n < (TDim-1)*3; ++n) {
                                auxRightHandSideVector(i) -= w_gauss * test_mat(i,m) * normal_proj_mat(m,n) * rData.C(n,k) * exp_strain_mat(k,j) * prev_sol(j);
                            }
                        }
                    }
                }
            }
        }

        // Contribution coming from the negative side boundary term
        const unsigned int n_int_neg_gauss = (rData.w_gauss_neg_int).size();

        for (unsigned int i_gauss = 0; i_gauss < n_int_neg_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_neg_int(i_gauss);

            // Get the negative side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_neg_int, i_gauss);

            // Get the negative side shape functions gradients Gauss pt. values
            const bounded_matrix<double, TNumNodes, TDim> aux_DN_DX = rData.DN_DX_neg_int[i_gauss];

            // Get the negative side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.neg_int_unit_normals[i_gauss];

            // Set the auxiliar test function matrix
            bounded_matrix<double, MatrixSize, TDim> test_mat = ZeroMatrix(MatrixSize, TDim);
            this->SetTestMatrix(aux_N, test_mat);

            // Set the normal projection matrix in Voigt notation
            bounded_matrix<double, TDim, (TDim-1)*3> normal_proj_mat = ZeroMatrix(TDim, (TDim-1)*3);
            this->SetVoigtNormalProjectionMatrix(side_normal, normal_proj_mat);

            // Set the strain matrix (expanded with zeros in the pressure rows)
            bounded_matrix<double, (TDim-1)*3, MatrixSize> exp_strain_mat = ZeroMatrix((TDim-1)*3, MatrixSize);
            this->SetExpandedStrainMatrix(aux_DN_DX, exp_strain_mat);

            // Pressure contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int d = 0; d < TDim; ++d) {
                        auxRightHandSideVector(i * BlockSize + d) += w_gauss * aux_N(i) * aux_N(j) * side_normal(d) * prev_sol(j * BlockSize + TDim);
                    }
                }
            }

            // Shear stress contribution
            for (unsigned int i = 0; i < MatrixSize; ++i) {
                for (unsigned int j = 0; j < MatrixSize; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        for (unsigned int k = 0; k < (TDim-1)*3; ++k) {
                            for (unsigned int n = 0; n < (TDim-1)*3; ++n) {
                                auxRightHandSideVector(i) -= w_gauss * test_mat(i,m) * normal_proj_mat(m,n) * rData.C(n,k) * exp_strain_mat(k,j) * prev_sol(j);
                            }
                        }
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
    */
    void AddSystemNormalVelocityPenaltyContribution(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const EmbeddedAusasElementDataStruct &rData) {

        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        array_1d<double, MatrixSize> solution_jump = ZeroVector(MatrixSize);

        // Obtain the previous iteration velocity solution
        GetPreviousSolutionVector(rData, prev_sol);

        // Compute the velocity diference to penalize
        if (this->Has(EMBEDDED_VELOCITY)) {
            const array_1d<double, 3> &embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            solution_jump = aux_embedded_vel - prev_sol;
        } else {
            solution_jump = - prev_sol;
        }

        // Compute the penalty coefficient
        const double pen_coef = ComputePenaltyCoefficient(rData);

        // Compute the LHS and RHS penalty contributions
        bounded_matrix<double, MatrixSize, MatrixSize> P_gamma = ZeroMatrix(MatrixSize, MatrixSize);

        // Contribution coming from the positive side penalty term
        const unsigned int n_int_pos_gauss = (rData.w_gauss_pos_int).size();
        for (unsigned int i_gauss = 0; i_gauss < n_int_pos_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_pos_int(i_gauss);

            // Get the positive side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss);

            // Get the positive side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.pos_int_unit_normals[i_gauss];

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        const unsigned int row = i * BlockSize + m;

                        for (unsigned int n = 0; n < TDim; ++n) {
                            const unsigned int col = j * BlockSize + n;
                            P_gamma(row, col) += pen_coef * w_gauss * aux_N(i) * side_normal(m) * side_normal(n) * aux_N(j);
                        }
                    }
                }
            }
        }

        // Contribution coming from the negative side penalty term
        const unsigned int n_int_neg_gauss = (rData.w_gauss_neg_int).size();
        for (unsigned int i_gauss = 0; i_gauss < n_int_neg_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_neg_int(i_gauss);

            // Get the negative side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_neg_int, i_gauss);

            // Get the negative side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.neg_int_unit_normals[i_gauss];

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        const unsigned int row = i * BlockSize + m;

                        for (unsigned int n = 0; n < TDim; ++n) {
                            const unsigned int col = j * BlockSize + n;
                            P_gamma(row, col) += pen_coef * w_gauss * aux_N(i) * side_normal(m) * side_normal(n) * aux_N(j);
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
    * This function adds the RHS contribution of the penalty no penetration imposition.
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    */
    void AddRHSNormalVelocityPenaltyContribution(
        VectorType &rRightHandSideVector,
        const EmbeddedAusasElementDataStruct &rData) {

        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        array_1d<double, MatrixSize> solution_jump = ZeroVector(MatrixSize);

        // Obtain the previous iteration velocity solution
        GetPreviousSolutionVector(rData, prev_sol);

        // Compute the velocity diference to penalize
        if (this->Has(EMBEDDED_VELOCITY)) {
            const array_1d<double, 3> &embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            solution_jump = aux_embedded_vel - prev_sol;
        } else {
            solution_jump = - prev_sol;
        }

        // Compute the penalty coefficient
        const double pen_coef = ComputePenaltyCoefficient(rData);

        // Compute the RHS penalty contributions
        array_1d<double, MatrixSize> P_gamma_RHS = ZeroVector(MatrixSize);

        // Contribution coming from the positive side penalty term
        const unsigned int n_int_pos_gauss = (rData.w_gauss_pos_int).size();
        for (unsigned int i_gauss = 0; i_gauss < n_int_pos_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_pos_int(i_gauss);

            // Get the positive side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss);

            // Get the positive side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.pos_int_unit_normals[i_gauss];

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        const unsigned int row = i * BlockSize + m;
                        for (unsigned int n = 0; n < TDim; ++n) {
                            P_gamma_RHS(row) += pen_coef * w_gauss * aux_N(i) * side_normal(m) * side_normal(n) * aux_N(j) * solution_jump(row);
                        }
                    }
                }
            }
        }

        // Contribution coming from the negative side penalty term
        const unsigned int n_int_neg_gauss = (rData.w_gauss_neg_int).size();
        for (unsigned int i_gauss = 0; i_gauss < n_int_neg_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_neg_int(i_gauss);

            // Get the negative side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_neg_int, i_gauss);

            // Get the negative side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.neg_int_unit_normals[i_gauss];

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        const unsigned int row = i * BlockSize + m;
                        for (unsigned int n = 0; n < TDim; ++n) {
                            P_gamma_RHS(row) += pen_coef * w_gauss * aux_N(i) * side_normal(m) * side_normal(n) * aux_N(j) * solution_jump(row);
                        }
                    }
                }
            }
        }

        // RHS assembly
        rRightHandSideVector += P_gamma_RHS;
    }

    /**
    * This function adds the local system contribution of the no penetration imposition,
    * by means of the pressure acting as a Lagrange multiplier.
    * @param rLeftHandSideMatrix: reference to the LHS matrix
    * @param rRightHandSideVector: reference to the RHS vector
    * @param rData: reference to element data structure
    */
    void AddSystemNormalVelocityLagrangeMultiplierContribution(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const EmbeddedAusasElementDataStruct &rData) {

        constexpr unsigned int BlockSize = TDim + 1;
        constexpr unsigned int MatrixSize = TNumNodes * BlockSize;

        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        array_1d<double, MatrixSize> solution_jump = ZeroVector(MatrixSize);

        // Obtain the previous iteration velocity solution
        GetPreviousSolutionVector(rData, prev_sol);

        // Compute the velocity diference to penalize
        if (this->Has(EMBEDDED_VELOCITY)) {
            const array_1d<double, 3> &embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            solution_jump = prev_sol - aux_embedded_vel;
        } else {
            solution_jump = prev_sol;
        }

        // Compute the LHS and RHS penalty contributions
        bounded_matrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // Drop the pressure rows and columns
        for (unsigned int i=0; i<TNumNodes; ++i) {
            const unsigned int aux = i*BlockSize + TDim;
            for (unsigned int j=0; j<MatrixSize; ++j) {
                rLeftHandSideMatrix(aux, j) = 0.0;
                rLeftHandSideMatrix(j, aux) = 0.0;
            }
        }

        // Contribution coming from the positive side Lagrange multiplier term
        const unsigned int n_int_pos_gauss = (rData.w_gauss_pos_int).size();
        for (unsigned int i_gauss = 0; i_gauss < n_int_pos_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_pos_int(i_gauss);

            // Get the positive side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss);

            // Get the positive side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.pos_int_unit_normals[i_gauss];

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                const unsigned int row = i * BlockSize + TDim;

                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        const unsigned int col = j * BlockSize + m;
                        const double aux_value = w_gauss * aux_N(i) * side_normal(m) * aux_N(j);
                        auxLeftHandSideMatrix(row, col) += aux_value;
                        auxLeftHandSideMatrix(col, row) += aux_value;
                    }
                }
            }
        }

        // Contribution coming from the negative side Lagrange multiplier term
        const unsigned int n_int_neg_gauss = (rData.w_gauss_neg_int).size();
        for (unsigned int i_gauss = 0; i_gauss < n_int_neg_gauss; ++i_gauss) {
            const double w_gauss = rData.w_gauss_neg_int(i_gauss);

            // Get the negative side shape functions Gauss pt. values
            const array_1d<double, TNumNodes> aux_N = row(rData.N_neg_int, i_gauss);

            // Get the negative side Gauss pt. normal
            const array_1d<double, 3> side_normal = rData.neg_int_unit_normals[i_gauss];

            // Compute and assemble the LHS contribution
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                const unsigned int row = i * BlockSize + TDim;

                for (unsigned int j = 0; j < TNumNodes; ++j) {
                    for (unsigned int m = 0; m < TDim; ++m) {
                        const unsigned int col = j * BlockSize + m;
                        const double aux_value = w_gauss * aux_N(i) * side_normal(m) * aux_N(j);
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
    * Auxiliar function to compute the element size
    * @return h: characteristic element size computed using the gradients
    */
    double ComputeH() {

        double aux_volume;
        array_1d<double, TNumNodes> aux_N;
        bounded_matrix<double, TNumNodes, TDim> aux_DN_DX;

        GeometryUtils::CalculateGeometryData(this->GetGeometry(), aux_DN_DX, aux_N, aux_volume);

        double h=0.0;
        for(unsigned int i=0; i<TNumNodes; i++) {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++) {
                h_inv += aux_DN_DX(i,k)*aux_DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }

        h = sqrt(h)/static_cast<double>(TNumNodes);

        return h;
    }

    /**
    * This function computes the penalty coefficient for the level set normal velocity imposition
    * @param rData: reference to element data structure
    */
    double ComputePenaltyCoefficient(const EmbeddedAusasElementDataStruct &rData) {

        // Compute the intersection area using the Gauss pts. weights
        double intersection_area = 0.0;
        for (unsigned int i_gauss = 0; i_gauss < (rData.w_gauss_pos_int).size(); ++i_gauss) {
            intersection_area += rData.w_gauss_pos_int(i_gauss);
        }

        // Compute the element average values
        double avg_rho = 0.0;
        double avg_visc = 0.0;
        array_1d<double, 3> avg_vel = ZeroVector(3);

        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
            avg_rho += rData.rho(i_node);
            avg_visc += rData.mu(i_node);
            avg_vel += row(rData.v, i_node);
        }

        avg_rho /= TNumNodes;
        avg_visc /= TNumNodes;
        avg_vel /= TNumNodes;

        const double v_norm = norm_2(avg_vel);

        // Compute the penalty constant
        const double pen_cons = avg_rho*std::pow(rData.h, TDim)/rData.dt +
                                avg_rho*avg_visc*std::pow(rData.h,TDim-2) +
                                avg_rho*v_norm*std::pow(rData.h, TDim-1);

        // Return the penalty coefficient
        const double K = 10.0;
        const double pen_coef = K * pen_cons / intersection_area;

        return pen_coef;
    }

    /**
    * This functions sets a vector containing the element previous solution
    * @param rData: reference to the element data structure
    * @param rPrevSolVector: reference to the previous solution vector
    */
    void GetPreviousSolutionVector(
        const EmbeddedAusasElementDataStruct& rData,
        array_1d<double, TNumNodes*(TDim+1)>& rPrevSolVector) {

        rPrevSolVector.clear();

        for (unsigned int i=0; i<TNumNodes; i++) {
            for (unsigned int comp=0; comp<TDim; comp++) {
                rPrevSolVector(i*(TDim+1)+comp) = rData.v(i,comp);
            }
            rPrevSolVector(i*(TDim+1)+TDim) = rData.p(i);
        }
    }

    /**
    * This functions computes the strain rate in Voigt notation
    * @param rData: reference to the element data structure
    * @param StrainSize: strain size (3 in 2D and 6 in 3D)
    */
    void ComputeStrain(
        EmbeddedAusasElementDataStruct& rData,
        const unsigned int StrainSize) {

        const bounded_matrix<double, TNumNodes, TDim>& v = rData.v;
        const bounded_matrix<double, TNumNodes, TDim>& DN = rData.DN_DX;

        // Compute strain (B*v)
        if (StrainSize == 6) {
            // 3D strain computation
            rData.strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
            rData.strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
            rData.strain[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
            rData.strain[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
            rData.strain[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
            rData.strain[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
        } else if (StrainSize == 3) {
            // 2D strain computation
            rData.strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
            rData.strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
            rData.strain[2] = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
        }
    }

    /**
    * This function calls the constitutive law to get the stress value
    * @param rData: reference to the element data structure
    * @param rCurrentProcessInfo: reference to the ProcessInfo
    */
    void ComputeConstitutiveResponse(
        EmbeddedAusasElementDataStruct &rData,
        const ProcessInfo& rCurrentProcessInfo) {

        const unsigned int strain_size = (TDim*3)-3;

        if(rData.C.size1() != strain_size) {
            rData.C.resize(strain_size, strain_size, false);
        } else if(rData.C.size2() != strain_size) {
            rData.C.resize(strain_size, strain_size, false);
        }

        if(rData.stress.size() != strain_size) {
            rData.stress.resize(strain_size,false);
        }

        if(rData.strain.size() != strain_size) {
            rData.strain.resize(strain_size,false);
        }

        this->ComputeStrain(rData, strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(this->GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Values.SetShapeFunctionsValues(rData.N);

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

    /**
    * This function computes the effective viscosity as the average of the lower diagonal constitutive tensor
    * @param rData: reference to the element data structure
    * @return mu_eff: computed effective viscosity value
    */
    double ComputeEffectiveViscosity(const EmbeddedAusasElementDataStruct& rData) {
        double mu_eff = 0.0;
        const unsigned int strain_size = (TDim-1)*3;
        for (unsigned int i=TDim; i<strain_size; ++i){mu_eff += rData.C(i,i);}
        mu_eff /= (strain_size - TDim);

        return mu_eff;
    }

    /**
    * This functions sets the auxiliar matrix to compute the normal projection in Voigt notation
    * @param rUnitNormal: reference to Gauss pt. unit normal vector
    * @return rVoigtNormProjMatrix: reference to the computed normal projection auxiliar matrix
    */
    void SetVoigtNormalProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        bounded_matrix<double, TDim, (TDim-1)*3>& rVoigtNormProjMatrix) {

        rVoigtNormProjMatrix.clear();

        if (TDim == 3) {
            rVoigtNormProjMatrix(0,0) = rUnitNormal(0);
            rVoigtNormProjMatrix(0,3) = rUnitNormal(1);
            rVoigtNormProjMatrix(0,5) = rUnitNormal(2);
            rVoigtNormProjMatrix(1,1) = rUnitNormal(1);
            rVoigtNormProjMatrix(1,3) = rUnitNormal(0);
            rVoigtNormProjMatrix(1,4) = rUnitNormal(2);
            rVoigtNormProjMatrix(2,2) = rUnitNormal(2);
            rVoigtNormProjMatrix(2,4) = rUnitNormal(1);
            rVoigtNormProjMatrix(2,5) = rUnitNormal(0);
        } else {
            rVoigtNormProjMatrix(0,0) = rUnitNormal(0);
            rVoigtNormProjMatrix(0,2) = rUnitNormal(1);
            rVoigtNormProjMatrix(1,1) = rUnitNormal(1);
            rVoigtNormProjMatrix(1,2) = rUnitNormal(0);
        }
    }

    /**
    * This functions sets the B strain matrix (pressure columns are set to zero)
    * @param rData: reference to element data structure (it contains the shape functions derivatives)
    * @return rBmatrix: reference to the computed B strain matrix
    */
    void SetExpandedStrainMatrix(
        const bounded_matrix<double, TNumNodes, TDim> &rDN_DX,
        bounded_matrix<double, (TDim-1)*3, TNumNodes*(TDim+1)> &rBmatrix) {

        constexpr unsigned int block_size = TDim + 1;
        rBmatrix.clear();

        // Set the shape function derivatives values
        if (TDim == 3) {
            for (unsigned int i = 0; i < TNumNodes; i++) {
                rBmatrix(0, i * block_size)     = rDN_DX(i, 0);
                rBmatrix(1, i * block_size + 1) = rDN_DX(i, 1);
                rBmatrix(2, i * block_size + 2) = rDN_DX(i, 2);
                rBmatrix(3, i * block_size)     = rDN_DX(i, 1);
                rBmatrix(3, i * block_size + 1) = rDN_DX(i, 0);
                rBmatrix(4, i * block_size + 1) = rDN_DX(i, 2);
                rBmatrix(4, i * block_size + 2) = rDN_DX(i, 1);
                rBmatrix(5, i * block_size)     = rDN_DX(i, 2);
                rBmatrix(5, i * block_size + 2) = rDN_DX(i, 0);
            }
        } else {
            for (unsigned int i = 0; i < TNumNodes; i++) {
                rBmatrix(0, i * block_size)     = rDN_DX(i, 0);
                rBmatrix(1, i * block_size + 1) = rDN_DX(i, 1);
                rBmatrix(2, i * block_size)     = rDN_DX(i, 1);
                rBmatrix(2, i * block_size + 1) = rDN_DX(i, 0);
            }
        }
    }

    /**
    * This functions sets the test function matrix given the Gauss pt. shape function values
    * @param rN: shape function values on a Gauss pt.
    * @return rTestMatrix: computed test function matrix
    */
    void SetTestMatrix(
        const array_1d<double, TNumNodes> &rN,
        bounded_matrix<double, (TDim+1)*TNumNodes, TDim> &rTestMatrix) {

        constexpr unsigned int block_size = TDim + 1;
        rTestMatrix.clear();

        // Set the test function matrix using the shape functions values
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            for (unsigned int d = 0; d < TDim; ++d) {
                rTestMatrix(i * block_size + d, d) = rN(i);
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
        // TODO: Serialize constitutive law
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        // TODO: Serialize constitutive law
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
