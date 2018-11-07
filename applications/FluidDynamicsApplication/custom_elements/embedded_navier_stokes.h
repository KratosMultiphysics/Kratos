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

// Project includes
#include "custom_elements/navier_stokes.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

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

    typedef typename BaseType::GeometryType::ShapeFunctionsGradientsType
                                                    ShapeFunctionsGradientsType;

    struct EmbeddedElementDataStruct : public ElementDataType {
        // Element geometry data
        MatrixType                  N_pos_side;             // Positive distance element side shape functions values
        ShapeFunctionsGradientsType DN_DX_pos_side;         // Positive distance element side shape functions gradients values
        VectorType                  w_gauss_pos_side;       // Positive distance element side Gauss pts. weights

        // Intersection geometry data
        MatrixType                  N_pos_int;              // Positive interface Gauss pts. shape functions values
        ShapeFunctionsGradientsType DN_DX_pos_int;          // Positive interface Gauss pts. shape functions gradients values
        VectorType                  w_gauss_pos_int;        // Positive interface Gauss pts. weights
        std::vector<VectorType>     pos_int_unit_normals;   // Positive interface unit normal vector in each Gauss pt.

        std::vector<unsigned int>   int_vec_identifiers;    // Interior (fluid) nodes identifiers
        std::vector<unsigned int>   out_vec_identifiers;    // Outside (stucture) nodes identifiers

        unsigned int    n_pos = 0;                          // Number of postivie distance nodes
        unsigned int    n_neg = 0;                          // Number of negative distance nodes
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
    ~EmbeddedNavierStokes() override {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, Element::PropertiesType::Pointer pProperties) const override {
        KRATOS_TRY
        return Kratos::make_shared< EmbeddedNavierStokes < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }


    Element::Pointer Create(IndexType NewId, Element::GeometryType::Pointer pGeom, Element::PropertiesType::Pointer pProperties) const override {
        return Kratos::make_shared< EmbeddedNavierStokes < TDim, TNumNodes > >(NewId, pGeom, pProperties);
    }


    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param rThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override {
        Element::Pointer pNewElement = Create(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

        pNewElement->SetData(this->GetData());
        pNewElement->SetFlags(this->GetFlags());

        return pNewElement;
    }

    /**
     * Fill the element data structure. If the element is split,
     * calls the modified shape functions calculator.
     * @param rData reference to the element data structure
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void FillEmbeddedElementData(
        EmbeddedElementDataStruct& rData,
        const ProcessInfo& rCurrentProcessInfo) {

        // Fill the basic element data (base class call)
        BaseType::FillElementData(rData, rCurrentProcessInfo);

        // Getting the nodal distance values
        array_1d<double, TNumNodes> distances;
        for(unsigned int i=0; i<TNumNodes; i++) {
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

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

        // If the element is split, get the modified shape functions
        if (rData.n_pos != 0 && rData.n_neg != 0){

            GeometryPointerType p_geom = this->pGetGeometry();

            // Construct the modified shape fucntions utility
            ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
            if (TNumNodes == 4) {
                p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, distances);
            } else {
                p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, distances);
            }

            // Call the fluid side modified shape functions calculator
            p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
                rData.N_pos_side,
                rData.DN_DX_pos_side,
                rData.w_gauss_pos_side,
                GeometryData::GI_GAUSS_2);

            // Call the fluid side interface modified shape functions calculator
            p_modified_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                rData.N_pos_int,
                rData.DN_DX_pos_int,
                rData.w_gauss_pos_int,
                GeometryData::GI_GAUSS_2);

            // Call the fluid side Gauss pts. unit normal calculator
            p_modified_sh_func->ComputePositiveSideInterfaceAreaNormals(
                rData.pos_int_unit_normals,
                GeometryData::GI_GAUSS_2);

            // Normalize the obtained area normals
            const double tol = std::pow(1e-3*rData.h, TDim-1); // Tolerance to avoid the unit normal to blow up
            const unsigned int n_gauss = (rData.pos_int_unit_normals).size();

            for (unsigned int i_gauss = 0;  i_gauss < n_gauss; ++i_gauss) {
                Vector& normal = rData.pos_int_unit_normals[i_gauss];
                const double n_norm = norm_2(normal);
                normal /= std::max(n_norm, tol);
            }
        }
    };

    /**
     * Calculates both LHS and RHS contributions
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override {

        KRATOS_TRY;

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        // Initialize LHS and RHS
        if (rLeftHandSideMatrix.size1() != MatrixSize) {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); // false says to not preserve existing storage!!
        } else if (rLeftHandSideMatrix.size2() != MatrixSize) {
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); // false says to not preserve existing storage!!
        }

        if (rRightHandSideVector.size() != MatrixSize) {
            rRightHandSideVector.resize(MatrixSize, false);            // false says to not preserve existing storage!!
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize, MatrixSize);   // LHS initialization
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);              // RHS initialization

        // Embedded data struct filling
        EmbeddedElementDataStruct data;
        this->FillEmbeddedElementData(data, rCurrentProcessInfo);

        // Element LHS and RHS contributions computation
        if(data.n_pos == TNumNodes) {
            ComputeElementAsFluid<MatrixSize>(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo);
        } else if ((data.n_pos != 0) && (data.n_neg != 0)) {
            ComputeElementAsMixed<MatrixSize>(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo);
        }

        KRATOS_CATCH("Error in embedded Navier-Stokes element CalculateLocalSystem method.");
    }

    /**
     * Calculates both LHS and RHS elemental contributions for those cases in where
     * all the nodes belong to the fluid domain.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rData reference to element data structure
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    template<unsigned int MatrixSize>
    void ComputeElementAsFluid(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        EmbeddedElementDataStruct& rData,
        const ProcessInfo& rCurrentProcessInfo) {

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;
        BoundedMatrix<double,MatrixSize, MatrixSize> lhs_local;

        // Shape functions Gauss points values
        // TODO: CHANGE THIS, USE THE GEOMETRY WITH A QUADRATURE
        BoundedMatrix<double, TNumNodes, TNumNodes> Ncontainer; // Container with the evaluation of the n shape functions in the n Gauss pts.
        BaseType::GetShapeFunctionsOnGauss(Ncontainer);

        // Loop on gauss point
        for(unsigned int igauss = 0; igauss<Ncontainer.size1(); igauss++) {
            noalias(rData.N) = row(Ncontainer, igauss);

            BaseType::ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            BaseType::ComputeGaussPointLHSContribution(lhs_local, rData);
            BaseType::ComputeGaussPointRHSContribution(rhs_local, rData);

            // All the Gauss pts. have the same weight so the accumulated contributions can be multiplied by volume/n_nodes at the end
            noalias(rLeftHandSideMatrix) += lhs_local;
            noalias(rRightHandSideVector) += rhs_local;
        }

        rLeftHandSideMatrix *= rData.volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= rData.volume/static_cast<double>(TNumNodes);
    }

    /**
    * Calculates both LHS and RHS elemental contributions for those cases in where
    * the element has both fluid and structure nodes.
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    * @param rCurrentProcessInfo reference to the ProcessInfo
    */
    template<unsigned int MatrixSize>
    void ComputeElementAsMixed(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        EmbeddedElementDataStruct& rData,
        const ProcessInfo& rCurrentProcessInfo) {

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_local;
        BoundedMatrix<double,MatrixSize, MatrixSize> lhs_local;

        // Gauss points loop
        const unsigned int n_gauss_pts = (rData.w_gauss_pos_side).size();

        for(unsigned int i_gauss = 0; i_gauss < n_gauss_pts; i_gauss++) {
            noalias(rData.N) = row(rData.N_pos_side, i_gauss);          // Take the new Gauss pts. shape functions values
            noalias(rData.DN_DX) = rData.DN_DX_pos_side[i_gauss];       // Take the new Gauss pts. shape functions gradients values
            const double weight = rData.w_gauss_pos_side(i_gauss);      // Subvolume Gauss pt. weights

            BaseType::ComputeConstitutiveResponse(rData, rCurrentProcessInfo);

            BaseType::ComputeGaussPointLHSContribution(lhs_local, rData);
            BaseType::ComputeGaussPointRHSContribution(rhs_local, rData);

            noalias(rLeftHandSideMatrix)  += weight*lhs_local;
            noalias(rRightHandSideVector) += weight*rhs_local;
        }

        // Add level set boundary terms, penalty and modified Nitche contributions
        AddBoundaryConditionElementContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rCurrentProcessInfo);
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    int Check(const ProcessInfo &rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        // Base element check
        int error_code = NavierStokes<TDim, TNumNodes>::Check(rCurrentProcessInfo);
        if (error_code != 0){
            return error_code;
        }

        // Specific embedded element check
        if (DISTANCE.Key() == 0){
            KRATOS_ERROR << "DISTANCE Key is 0. Check if the application was correctly registered.";
        }

        for (unsigned int i = 0; i < (this->GetGeometry()).size(); ++i){
            if (this->GetGeometry()[i].SolutionStepsDataHas(DISTANCE) == false){
                KRATOS_ERROR << "missing VELOCITY variable on solution step data for node " << this->GetGeometry()[i].Id();
            }
        }

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * Given a vector variable, this function computes its value inside de element.
     * If the function has not implemented this variable computation, throws an error.
     * @param rVariable Variable to be computed.
     * @param rOutput Reference to the output array.
     * @param rCurrentProcessInfo Reference to the process info.
     */
    void Calculate(
        const Variable<array_1d<double, 3>> &rVariable,
        array_1d<double, 3> &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override {

        rOutput = ZeroVector(3);

        // If the element is split, integrate sigma·n over the interface
        // Note that in the ausas formulation, both interface sides need to be integrated
        if (rVariable == DRAG_FORCE) {

            EmbeddedElementDataStruct data;
            this->FillEmbeddedElementData(data, rCurrentProcessInfo);

            // Check if the element is split
            if (data.n_pos != 0 && data.n_neg != 0){

                // Integrate positive interface side drag
                const unsigned int n_int_pos_gauss = (data.w_gauss_pos_int).size();
                for (unsigned int i_gauss = 0; i_gauss < n_int_pos_gauss; ++i_gauss) {
                    // Get Gauss pt. data
                    const double w_gauss = data.w_gauss_pos_int(i_gauss);
                    const array_1d<double, TNumNodes> aux_N = row(data.N_pos_int, i_gauss);
                    const array_1d<double, 3> side_normal = data.pos_int_unit_normals[i_gauss];

                    // Obtain Gauss pt. pressure
                    const double p_gauss = inner_prod(aux_N, data.p);

                    // Call the constitutive law to compute the shear contribution
                    // Recall to set data.N and data.DN_DX (required by the constitutive law)
                    noalias(data.N) = aux_N;
                    noalias(data.DN_DX) = data.DN_DX_pos_int[i_gauss];
                    this->ComputeConstitutiveResponse(data, rCurrentProcessInfo);

                    // Get the Voigt notation normal projection matrix
                    BoundedMatrix<double, TDim, (TDim - 1) * 3> normal_proj_mat = ZeroMatrix(TDim, (TDim - 1) * 3);
                    this->SetVoigtNormalProjectionMatrix(side_normal, normal_proj_mat);

                    // Add the shear and pressure drag contributions
                    const array_1d<double, TDim> shear_proj = w_gauss * prod(normal_proj_mat, data.stress);
                    for (unsigned int i = 0; i < TDim ; ++i){
                        rOutput(i) -= shear_proj(i);
                    }
                    rOutput += w_gauss * p_gauss * side_normal;
                }
            }
        }
        else
        {
            KRATOS_ERROR << "Calculate method not implemented for the requested variable.";
        }
    }

    void Calculate(const Variable<double>& rVariable,
                   double& Output,
                   const ProcessInfo& rCurrentProcessInfo) override
    {}

    void Calculate(const Variable<Vector >& rVariable,
                   Vector& Output,
                   const ProcessInfo& rCurrentProcessInfo) override
    {}

    void Calculate(const Variable<Matrix >& rVariable,
                   Matrix& Output,
                   const ProcessInfo& rCurrentProcessInfo) override
    {}

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

    EmbeddedNavierStokes() : NavierStokes<TDim, TNumNodes>() {}

    ///@}
    ///@name Protected Operations
    ///@{

    /**
    * This functions adds the contribution of the boundary terms in the level set cut
    * These terms, which do not vanish at the level set since the test function is not zero
    * at the intersection points, come from the integration by parts of the stress term.
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddIntersectionBoundaryTermsContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData) {

        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Declare auxiliar arrays
        BoundedMatrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        const unsigned int n_gauss_total = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss_int = 0; i_gauss_int < n_gauss_total; ++i_gauss_int) {

            // Get the current Gauss pt. data
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss_int);                    // Shape function values
            const BoundedMatrix<double, TNumNodes, TDim> aux_DN_DX = rData.DN_DX_pos_int[i_gauss_int];     // Shape function gradient values
            const double weight = rData.w_gauss_pos_int(i_gauss_int);                                       // Intersection Gauss pt. weight
            const array_1d<double, 3> aux_unit_normal = rData.pos_int_unit_normals[i_gauss_int];            // Gauss pt. unit normal

            // Set the current Gauss pt. Voigt notation normal projection matrix
            BoundedMatrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix = ZeroMatrix(TDim, (TDim-1)*3);
            SetVoigtNormalProjectionMatrix(aux_unit_normal, voigt_normal_projection_matrix);

            // Set the current Gauss pt. strain matrix
            BoundedMatrix<double, (TDim-1)*3, MatrixSize> B_matrix = ZeroMatrix((TDim-1)*3, MatrixSize);
            SetInterfaceStrainMatrix(aux_DN_DX, B_matrix);

            // Compute some Gauss pt. auxiliar matrices
            const BoundedMatrix<double, TDim, (TDim-1)*3> aux_matrix_AC = prod(voigt_normal_projection_matrix, rData.C);
            const BoundedMatrix<double, (TDim-1)*3, MatrixSize> aux_matrix_ACB = prod(aux_matrix_AC, B_matrix);

            // Fill the pressure to Voigt notation operator matrix
            BoundedMatrix<double, (TDim-1)*3, MatrixSize> pres_to_voigt_matrix_op = ZeroMatrix((TDim-1)*3, MatrixSize);
            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    pres_to_voigt_matrix_op(comp, i*BlockSize+TDim) = aux_N(i);
                }
            }

            // Set the shape functions auxiliar transpose matrix
            BoundedMatrix<double, MatrixSize, TDim> N_aux_trans = ZeroMatrix(MatrixSize, TDim);
            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    N_aux_trans(i*BlockSize+comp, comp) = aux_N(i);
                }
            }

            // Contribution coming fron the shear stress operator
            noalias(auxLeftHandSideMatrix) += weight*prod(N_aux_trans, aux_matrix_ACB);

            // Contribution coming from the pressure terms
            const BoundedMatrix<double, MatrixSize, (TDim-1)*3> N_voigt_proj_matrix = prod(N_aux_trans, voigt_normal_projection_matrix);
            noalias(auxLeftHandSideMatrix) -= weight*prod(N_voigt_proj_matrix, pres_to_voigt_matrix_op);
        }

        // LHS assembly
        noalias(rLeftHandSideMatrix) -= auxLeftHandSideMatrix;

        // RHS assembly
        noalias(rRightHandSideVector) += prod(auxLeftHandSideMatrix, prev_sol);
    }

    /**
     * This function computes the penalty coefficient for the level set BC imposition
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rData reference to element data structure
     */
    double ComputePenaltyCoefficient(
        const EmbeddedElementDataStruct& rData,
        const ProcessInfo &rCurrentProcessInfo)
    {
        // Compute the intersection area using the Gauss pts. weights
        double intersection_area = 0.0;
        for (unsigned int i_gauss = 0; i_gauss < (rData.w_gauss_pos_int).size(); ++i_gauss) {
            intersection_area += rData.w_gauss_pos_int(i_gauss);
        }

        // Compute the element average values
        double avg_rho = 0.0;
        double avg_visc = 0.0;
        array_1d<double, TDim> avg_vel = ZeroVector(TDim);

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
        const double K = rCurrentProcessInfo[PENALTY_COEFFICIENT];
        const double pen_coef = K * pen_cons / intersection_area;

        return pen_coef;

    }

    /**
    * This functions adds the penalty extra term level set contribution.
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddBoundaryConditionPenaltyContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData,
        const ProcessInfo &rCurrentProcessInfo)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Set the penalty matrix
        MatrixType P_gamma(TNumNodes, TNumNodes);
        noalias(P_gamma) = ZeroMatrix(TNumNodes, TNumNodes);

        const unsigned int n_gauss_total = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss_int = 0; i_gauss_int < n_gauss_total; ++i_gauss_int) {
            const double weight = rData.w_gauss_pos_int(i_gauss_int);
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss_int);
            P_gamma += weight*outer_prod(aux_N,aux_N);
        }

        // Multiply the penalty matrix by the penalty coefficient
        double pen_coef = ComputePenaltyCoefficient(rData, rCurrentProcessInfo);
        P_gamma *= pen_coef;

        VectorType auxRightHandSideVector = ZeroVector(MatrixSize);
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // LHS penalty contribution assembly (symmetric mass matrix)
        for (unsigned int i = 0; i<TNumNodes; i++) {
            // Diagonal terms
            for (unsigned int comp = 0; comp<TDim; comp++) {
                auxLeftHandSideMatrix(i*BlockSize+comp, i*BlockSize+comp) = P_gamma(i,i);
            }
            // Off-diagonal terms
            for (unsigned int j = i+1; j<TNumNodes; j++) {
                for (unsigned int comp = 0; comp<TDim; comp++) {
                    auxLeftHandSideMatrix(i*BlockSize+comp, j*BlockSize+comp) = P_gamma(i,j);
                    auxLeftHandSideMatrix(j*BlockSize+comp, i*BlockSize+comp) = P_gamma(i,j);
                }
            }
        }

        noalias(rLeftHandSideMatrix) += auxLeftHandSideMatrix;

        // RHS penalty contribution assembly
        if (this->Has(EMBEDDED_VELOCITY)) {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; i++) {
                for (unsigned int comp=0; comp<TDim; comp++) {
                    aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            noalias(rRightHandSideVector) += prod(auxLeftHandSideMatrix, aux_embedded_vel);
        }

        noalias(rRightHandSideVector) -= prod(auxLeftHandSideMatrix, prev_sol); // Residual contribution assembly
    }

    /**
    * This functions adds the level set strong boundary condition imposition contribution.
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddBoundaryConditionModifiedNitcheContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData) {

        constexpr unsigned int BlockSize = TDim+1;                 // Block size
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;   // Matrix size

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Compute the BCs imposition matrices
        MatrixType M_gamma = ZeroMatrix(rData.n_neg, rData.n_neg);  // Outside nodes matrix (Nitche contribution)
        MatrixType N_gamma = ZeroMatrix(rData.n_neg, rData.n_pos);  // Interior nodes matrix (Nitche contribution)
        MatrixType f_gamma = ZeroMatrix(rData.n_neg, TNumNodes);    // Matrix to compute the RHS (Nitche contribution)

        VectorType aux_out(rData.n_neg);
        VectorType aux_int(rData.n_pos);

        const unsigned int n_gauss_total = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss_int = 0; i_gauss_int < n_gauss_total; ++i_gauss_int) {
            const double weight = rData.w_gauss_pos_int(i_gauss_int);

            const VectorType aux_cut = row(rData.N_pos_int, i_gauss_int);

            for (unsigned int i_out = 0; i_out < rData.n_neg; ++i_out) {
                const unsigned int i_out_nodeid = rData.out_vec_identifiers[i_out];
                aux_out(i_out) = aux_cut(i_out_nodeid);
            }

            for (unsigned int i_int = 0; i_int < rData.n_pos; ++i_int) {
                const unsigned int i_int_nodeid = rData.int_vec_identifiers[i_int];
                aux_int(i_int) = aux_cut(i_int_nodeid);
            }

            M_gamma += weight*outer_prod(aux_out,aux_out);
            N_gamma += weight*outer_prod(aux_out,aux_int);
            f_gamma += weight*outer_prod(aux_out,aux_cut);
        }

        // Declare auxLeftHandSideMatrix
        MatrixType auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        // LHS outside nodes contribution assembly
        // Outer nodes contribution assembly
        for (unsigned int i = 0; i<rData.n_neg; i++) {
            unsigned int out_node_row_id = rData.out_vec_identifiers[i];

            for (unsigned int j = 0; j<rData.n_neg; j++) {
                unsigned int out_node_col_id = rData.out_vec_identifiers[j];

                for (unsigned int comp = 0; comp<TDim; comp++) {
                    auxLeftHandSideMatrix(out_node_row_id*BlockSize+comp, out_node_col_id*BlockSize+comp) = M_gamma(i, j);
                }
            }
        }

        // Interior nodes contribution assembly
        for (unsigned int i = 0; i<rData.n_neg; i++) {
            unsigned int out_node_row_id = rData.out_vec_identifiers[i];

            for (unsigned int j = 0; j<rData.n_pos; j++) {
                unsigned int int_node_col_id = rData.int_vec_identifiers[j];

                for (unsigned int comp = 0; comp<TDim; comp++) {
                    auxLeftHandSideMatrix(out_node_row_id*BlockSize+comp, int_node_col_id*BlockSize+comp) = N_gamma(i, j);
                }
            }
        }

        // LHS outside Nitche contribution assembly
        noalias(rLeftHandSideMatrix) += auxLeftHandSideMatrix;

        // RHS outside Nitche contribution assembly
        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        noalias(rRightHandSideVector) -= prod(auxLeftHandSideMatrix, prev_sol);

        // Compute f_gamma if level set velocity is not 0
        if (this->Has(EMBEDDED_VELOCITY)) {
            auxLeftHandSideMatrix.clear();

            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> aux_embedded_vel = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; i++) {
                aux_embedded_vel(i*BlockSize) = embedded_vel(0);
                aux_embedded_vel(i*BlockSize+1) = embedded_vel(1);
                aux_embedded_vel(i*BlockSize+2) = embedded_vel(2);
            }

            // Asemble the RHS f_gamma contribution
            for (unsigned int i=0; i<rData.n_neg; i++) {
                unsigned int out_node_row_id = rData.out_vec_identifiers[i];

                for (unsigned int j=0; j<TNumNodes; j++) {
                    for (unsigned int comp = 0; comp<TDim; comp++) {
                        auxLeftHandSideMatrix(out_node_row_id*BlockSize+comp, j*BlockSize+comp) = f_gamma(i,j);
                    }
                }
            }

            noalias(rRightHandSideVector) += prod(auxLeftHandSideMatrix, aux_embedded_vel);
        }
    }

    /**
    * This drops the outer nodes velocity constributions in both LHS and RHS matrices.
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void DropOuterNodesVelocityContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData) {

        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Set the LHS and RHS u_out rows to zero (outside nodes used to impose the BC)
        for (unsigned int i=0; i<rData.n_neg; ++i) {
            const unsigned int out_node_row_id = rData.out_vec_identifiers[i];

            for (unsigned int j=0; j<TDim; ++j) {
                // LHS matrix u_out zero set (note that just the velocity rows are set to 0)
                for (unsigned int col = 0; col<MatrixSize; col++) {
                    rLeftHandSideMatrix(out_node_row_id*BlockSize+j, col) = 0.0;
                }

                // RHS vector u_out zero set (note that just the velocity rows are set to 0)
                rRightHandSideVector(out_node_row_id*BlockSize+j) = 0.0;
            }
        }
    }

    /**
    * This function adds the Nitsche normal component of the penalty contribution (Winter formulation).
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddSlipWinterNormalPenaltyContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData,
        const ProcessInfo &rCurrentProcessInfo)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Nitsche coefficient computation
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);

        // Compute the element average velocity norm
        double v_norm = 0.0;
        for (unsigned int comp=0; comp<TDim; ++comp) {
            double aux_vel = 0.0;
            for (unsigned int j=0; j<TNumNodes; ++j) {
                aux_vel += rData.v(j,comp);
            }
            aux_vel /= TNumNodes;
            v_norm += aux_vel*aux_vel;
        }
        v_norm = std::sqrt(v_norm);

        // Compute the element average density
        double avg_rho = 0.0;
        for (unsigned int j=0; j<TNumNodes; ++j) {
            avg_rho += rData.rho(j);
        }
        avg_rho /= TNumNodes;

        // Compute the Nitsche coefficient (considering the Winter stabilization term)
        const double penalty = 1.0/rCurrentProcessInfo[PENALTY_COEFFICIENT];
        const double cons_coef = (eff_mu + eff_mu + avg_rho*v_norm*rData.h + avg_rho*rData.h*rData.h/rData.dt)/(rData.h*penalty);

        // Declare auxiliar arrays
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);
        BoundedMatrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        const unsigned int n_gauss_total = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss_int = 0; i_gauss_int < n_gauss_total; ++i_gauss_int) {
            // Get the Gauss pt. data
            const double weight = rData.w_gauss_pos_int(i_gauss_int);
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss_int);
            const array_1d<double, 3> aux_unit_normal = rData.pos_int_unit_normals[i_gauss_int];

            // Set the shape functions auxiliar matrices
            BoundedMatrix<double, TDim, MatrixSize> N_aux = ZeroMatrix(TDim, MatrixSize);
            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    N_aux(comp,i*BlockSize+comp) = aux_N(i);
                }
            }
            const BoundedMatrix<double, MatrixSize, TDim> N_aux_trans = trans(N_aux);

            // Set the normal projection matrix nxn
            BoundedMatrix<double, TDim, TDim> normal_projection_matrix;
            SetNormalProjectionMatrix(aux_unit_normal, normal_projection_matrix);

            // Compute the current cut point auxLHS contribution
            const BoundedMatrix<double, MatrixSize, TDim> aux_1 = prod(N_aux_trans, normal_projection_matrix);
            const BoundedMatrix<double, MatrixSize, MatrixSize> aux_2 = prod(aux_1, N_aux);
            noalias(auxLeftHandSideMatrix) += cons_coef*weight*aux_2;
        }

        // If level set velocity is not 0, add its contribution to the RHS
        if (this->Has(EMBEDDED_VELOCITY)) {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> embedded_vel_exp = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    embedded_vel_exp(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            noalias(auxRightHandSideVector) += prod(auxLeftHandSideMatrix, embedded_vel_exp);
        }

        // LHS outside Nitche contribution assembly
        noalias(rLeftHandSideMatrix) += auxLeftHandSideMatrix;

        // RHS outside Nitche contribution assembly
        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        noalias(rRightHandSideVector) += auxRightHandSideVector;
        noalias(rRightHandSideVector) -= prod(auxLeftHandSideMatrix, prev_sol);
    }

    /**
    * This function adds the Nitsche normal component of the symmetric counterpart of the fluxes (Winter formulation).
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddSlipWinterNormalSymmetricCounterpartContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData) {

        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Set if the shear stress term is adjoint consistent (1.0) or not (-1.0)
        const double adjoint_consistency_term = -1.0;

        // Declare auxiliar arrays
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);
        BoundedMatrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);

        const unsigned int n_gauss_total = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss_int = 0; i_gauss_int < n_gauss_total; ++i_gauss_int) {
            // Get the Gauss pt. data
            const double weight = rData.w_gauss_pos_int(i_gauss_int);
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss_int);
            const BoundedMatrix<double, TNumNodes, TDim> aux_DN_DX = rData.DN_DX_pos_int[i_gauss_int];
            const array_1d<double, 3> aux_unit_normal = rData.pos_int_unit_normals[i_gauss_int];

            // Fill the pressure to Voigt notation operator normal projected matrix
            BoundedMatrix<double, MatrixSize, TDim> trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(MatrixSize, TDim);
            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    trans_pres_to_voigt_matrix_normal_op(i*BlockSize+TDim, comp) = aux_N(i)*aux_unit_normal(comp);
                }
            }

            // Set the shape functions auxiliar matrix
            BoundedMatrix<double, TDim, MatrixSize> N_aux = ZeroMatrix(TDim, MatrixSize);
            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    N_aux(comp,i*BlockSize+comp) = aux_N(i);
                }
            }

            // Get the strain matrix
            BoundedMatrix<double, (TDim-1)*3, MatrixSize> B_matrix = ZeroMatrix((TDim-1)*3, MatrixSize);
            SetInterfaceStrainMatrix(aux_DN_DX, B_matrix);

            // Get the normal projection matrix
            BoundedMatrix<double, TDim, TDim> normal_projection_matrix = ZeroMatrix(TDim, TDim);
            SetNormalProjectionMatrix(aux_unit_normal, normal_projection_matrix);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix = ZeroMatrix(TDim, (TDim-1)*3);
            SetVoigtNormalProjectionMatrix(aux_unit_normal, voigt_normal_projection_matrix);

            // Compute some Gauss pt. auxiliar matrices
            const BoundedMatrix<double, MatrixSize, (TDim-1)*3> aux_matrix_BC = prod(trans(B_matrix), trans(rData.C));
            const BoundedMatrix<double, (TDim-1)*3, TDim> aux_matrix_APnorm = prod(trans(voigt_normal_projection_matrix), normal_projection_matrix);
            const BoundedMatrix<double, MatrixSize, TDim> aux_matrix_BCAPnorm = prod(aux_matrix_BC, aux_matrix_APnorm);

            // Contribution coming fron the shear stress operator
            noalias(auxLeftHandSideMatrix) += adjoint_consistency_term*weight*prod(aux_matrix_BCAPnorm, N_aux);

            // Contribution coming from the pressure terms
            const BoundedMatrix<double, MatrixSize, TDim> aux_matrix_VPnorm = prod(trans_pres_to_voigt_matrix_normal_op, normal_projection_matrix);
            noalias(auxLeftHandSideMatrix) += weight*prod(aux_matrix_VPnorm, N_aux);

        }

        // LHS outside Nitche contribution assembly
        noalias(rLeftHandSideMatrix) -= auxLeftHandSideMatrix; // The minus sign comes from the Nitsche formulation

        // RHS outside Nitche contribution assembly
        // If level set velocity is not 0, add its contribution to the RHS
        if (this->Has(EMBEDDED_VELOCITY)) {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> embedded_vel_exp = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    embedded_vel_exp(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            noalias(auxRightHandSideVector) += prod(auxLeftHandSideMatrix, embedded_vel_exp);
        }

        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        noalias(rRightHandSideVector) -= auxRightHandSideVector;
        noalias(rRightHandSideVector) += prod(auxLeftHandSideMatrix, prev_sol);
    }

    /**
    * This function adds the Nitsche tangential component of the penalty contribution (Winter formulation).
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddSlipWinterTangentialPenaltyContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData,
        const ProcessInfo &rCurrentProcessInfo)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Compute the penalty coefficients
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);
        const double penalty = 1.0/rCurrentProcessInfo[PENALTY_COEFFICIENT];
        const double slip_length = rCurrentProcessInfo[SLIP_LENGTH];
        const double coeff_1 = slip_length / (slip_length + penalty*rData.h);
        const double coeff_2 = eff_mu / (slip_length + penalty*rData.h);

        // Declare auxiliar arrays
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);
        BoundedMatrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix_1 = ZeroMatrix(MatrixSize, MatrixSize); // Adds the contribution coming from the tangential component of the Cauchy stress vector
        BoundedMatrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix_2 = ZeroMatrix(MatrixSize, MatrixSize); // Adds the contribution generated by the viscous shear force generated by the velocity

        const unsigned int n_gauss_total = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss_int = 0; i_gauss_int < n_gauss_total; ++i_gauss_int) {
            // Get the Gauss pt. data
            const double weight = rData.w_gauss_pos_int(i_gauss_int);
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss_int);
            const BoundedMatrix<double, TNumNodes, TDim> aux_DN_DX = rData.DN_DX_pos_int[i_gauss_int];
            const array_1d<double, 3> aux_unit_normal = rData.pos_int_unit_normals[i_gauss_int];

            // Set the shape functions auxiliar matrix
            BoundedMatrix<double, MatrixSize, TDim> N_aux_trans = ZeroMatrix(MatrixSize, TDim);
            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    N_aux_trans(i*BlockSize+comp, comp) = aux_N(i);
                }
            }

            // Get the strain matrix
            BoundedMatrix<double, (TDim-1)*3, MatrixSize> B_matrix = ZeroMatrix((TDim-1)*3, MatrixSize);
            SetInterfaceStrainMatrix(aux_DN_DX, B_matrix);

            // Get the normal projection matrix
            BoundedMatrix<double, TDim, TDim> tangential_projection_matrix = ZeroMatrix(TDim, TDim);
            SetTangentialProjectionMatrix(aux_unit_normal, tangential_projection_matrix);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix = ZeroMatrix(TDim, (TDim-1)*3);
            SetVoigtNormalProjectionMatrix(aux_unit_normal, voigt_normal_projection_matrix);

            // Compute some Gauss pt. auxiliar matrices
            const BoundedMatrix<double, (TDim-1)*3, MatrixSize> aux_matrix_CB = prod(rData.C, B_matrix);
            const BoundedMatrix<double, (TDim-1)*3, TDim> aux_matrix_PtangA = prod(tangential_projection_matrix, voigt_normal_projection_matrix);
            const BoundedMatrix<double, MatrixSize, TDim> aux_matrix_PtangACB = prod(aux_matrix_PtangA, aux_matrix_CB);

            // Contribution coming from the traction vector tangencial component
            noalias(auxLeftHandSideMatrix_1) += coeff_1*weight*prod(N_aux_trans, aux_matrix_PtangACB);

            // Contribution coming from the shear force generated by the velocity jump
            const BoundedMatrix<double, MatrixSize, TDim> aux_matrix_N_trans_tang = prod(N_aux_trans, tangential_projection_matrix);
            noalias(auxLeftHandSideMatrix_2) += coeff_2*weight*prod(aux_matrix_N_trans_tang, trans(N_aux_trans));
        }

        // LHS outside Nitche contribution assembly
        noalias(rLeftHandSideMatrix) += auxLeftHandSideMatrix_1;
        noalias(rLeftHandSideMatrix) += auxLeftHandSideMatrix_2;

        // RHS outside Nitche contribution assembly
        // If level set velocity is not 0, add its contribution to the RHS
        if (this->Has(EMBEDDED_VELOCITY)) {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> embedded_vel_exp = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    embedded_vel_exp(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            noalias(auxRightHandSideVector) += prod(auxLeftHandSideMatrix_2, embedded_vel_exp);
        }

        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        noalias(rRightHandSideVector) += auxRightHandSideVector;
        noalias(rRightHandSideVector) -= prod(auxLeftHandSideMatrix_1, prev_sol);
        noalias(rRightHandSideVector) -= prod(auxLeftHandSideMatrix_2, prev_sol);
    }

    /**
    * This function adds the Nitsche tangential component of the symmetric counterpart of the fluxes (Winter formulation).
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddSlipWinterTangentialSymmetricCounterpartContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData,
        const ProcessInfo &rCurrentProcessInfo)
    {
        constexpr unsigned int BlockSize = TDim+1;
        constexpr unsigned int MatrixSize = TNumNodes*BlockSize;

        // Obtain the previous iteration velocity solution
        array_1d<double, MatrixSize> prev_sol = ZeroVector(MatrixSize);
        GetPreviousSolutionVector(rData, prev_sol);

        // Set if the shear stress term is adjoint consistent (1.0) or not (-1.0)
        const double adjoint_consistency_term = -1.0;

        // Compute the coefficients
        const double eff_mu = BaseType::ComputeEffectiveViscosity(rData);
        const double penalty = 1.0/rCurrentProcessInfo[PENALTY_COEFFICIENT];
        const double slip_length = rCurrentProcessInfo[SLIP_LENGTH];
        const double coeff_1 = slip_length*penalty*rData.h / (slip_length + penalty*rData.h);
        const double coeff_2 = eff_mu*penalty*rData.h / (slip_length + penalty*rData.h);

        // Declare auxiliar arrays
        array_1d<double, MatrixSize> auxRightHandSideVector = ZeroVector(MatrixSize);
        BoundedMatrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix_1 = ZeroMatrix(MatrixSize, MatrixSize); // Adds the contribution coming from the tangential component of the Cauchy stress vector
        BoundedMatrix<double, MatrixSize, MatrixSize> auxLeftHandSideMatrix_2 = ZeroMatrix(MatrixSize, MatrixSize); // Adds the contribution generated by the viscous shear force generated by the velocity

        const unsigned int n_gauss_total = (rData.w_gauss_pos_int).size();

        for (unsigned int i_gauss_int = 0; i_gauss_int < n_gauss_total; ++i_gauss_int) {
            // Get the Gauss pt. data
            const double weight = rData.w_gauss_pos_int(i_gauss_int);
            const array_1d<double, TNumNodes> aux_N = row(rData.N_pos_int, i_gauss_int);
            const BoundedMatrix<double, TNumNodes, TDim> aux_DN_DX = rData.DN_DX_pos_int[i_gauss_int];
            const array_1d<double, 3> aux_unit_normal = rData.pos_int_unit_normals[i_gauss_int];

            // Set the shape functions auxiliar matrix
            BoundedMatrix<double, TDim, MatrixSize> N_aux = ZeroMatrix(TDim, MatrixSize);
            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    N_aux(comp, i*BlockSize+comp) = aux_N(i);
                }
            }

            // Get the strain matrix
            BoundedMatrix<double, (TDim-1)*3, MatrixSize> B_matrix = ZeroMatrix((TDim-1)*3, MatrixSize);
            SetInterfaceStrainMatrix(aux_DN_DX, B_matrix);

            // Get the normal projection matrix
            BoundedMatrix<double, TDim, TDim> tangential_projection_matrix = ZeroMatrix(TDim, TDim);
            SetTangentialProjectionMatrix(aux_unit_normal, tangential_projection_matrix);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, TDim, (TDim-1)*3> voigt_normal_projection_matrix = ZeroMatrix(TDim, (TDim-1)*3);
            SetVoigtNormalProjectionMatrix(aux_unit_normal, voigt_normal_projection_matrix);

            // Compute some Gauss pt. auxiliar matrices
            const BoundedMatrix<double, MatrixSize, TDim> aux_matrix_BtransAtrans = prod(trans(B_matrix), trans(voigt_normal_projection_matrix));
            const BoundedMatrix<double, MatrixSize, TDim> aux_matrix_BtransAtransPtan = prod(aux_matrix_BtransAtrans, tangential_projection_matrix);
            const BoundedMatrix<double, (TDim-1)*3, MatrixSize> aux_matrix_CB = prod(rData.C, B_matrix);
            const BoundedMatrix<double, TDim, MatrixSize> aux_matrix_ACB = prod(voigt_normal_projection_matrix, aux_matrix_CB);
            const BoundedMatrix<double, MatrixSize, MatrixSize> aux_matrix_BtransAtransPtanACB = prod(aux_matrix_BtransAtransPtan, aux_matrix_ACB);

            // Contribution coming from the traction vector tangencial component
            noalias(auxLeftHandSideMatrix_1) += adjoint_consistency_term*coeff_1*weight*aux_matrix_BtransAtransPtanACB;

            // Contribution coming from the shear force generated by the velocity jump
            noalias(auxLeftHandSideMatrix_2) += adjoint_consistency_term*coeff_2*weight*prod(aux_matrix_BtransAtransPtan, N_aux);

        }

        // LHS outside Nitche contribution assembly
        noalias(rLeftHandSideMatrix) -= auxLeftHandSideMatrix_1;
        noalias(rLeftHandSideMatrix) -= auxLeftHandSideMatrix_2;

        // RHS outside Nitche contribution assembly
        // If level set velocity is not 0, add its contribution to the RHS
        if (this->Has(EMBEDDED_VELOCITY)) {
            const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            array_1d<double, MatrixSize> embedded_vel_exp = ZeroVector(MatrixSize);

            for (unsigned int i=0; i<TNumNodes; ++i) {
                for (unsigned int comp=0; comp<TDim; ++comp) {
                    embedded_vel_exp(i*BlockSize+comp) = embedded_vel(comp);
                }
            }

            noalias(auxRightHandSideVector) += prod(auxLeftHandSideMatrix_2, embedded_vel_exp);
        }

        // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
        noalias(rRightHandSideVector) -= auxRightHandSideVector;
        noalias(rRightHandSideVector) += prod(auxLeftHandSideMatrix_1, prev_sol);
        noalias(rRightHandSideVector) += prod(auxLeftHandSideMatrix_2, prev_sol);
    }

    /**
    * This functions collects and adds all the level set boundary condition contributions
    * @param rLeftHandSideMatrix reference to the LHS matrix
    * @param rRightHandSideVector reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddBoundaryConditionElementContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const EmbeddedElementDataStruct& rData,
        const ProcessInfo &rCurrentProcessInfo)
    {
        // Compute and assemble the boundary terms comping from the integration by parts
        AddIntersectionBoundaryTermsContribution(rLeftHandSideMatrix, rRightHandSideVector, rData);

        if (this->Is(SLIP)) {
            // Winter Navier-slip condition
            AddSlipWinterNormalPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rCurrentProcessInfo);
            AddSlipWinterNormalSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, rData);
            AddSlipWinterTangentialPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rCurrentProcessInfo);
            AddSlipWinterTangentialSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rCurrentProcessInfo);

        } else {
            // First, compute and assemble the penalty level set BC imposition contribution
            // Secondly, compute and assemble the modified Nitche method level set BC imposition contribution (Codina and Baiges, 2009)
            // Note that the Nitche contribution has to be computed the last since it drops the outer nodes rows previous constributions
            AddBoundaryConditionPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, rData, rCurrentProcessInfo);
            DropOuterNodesVelocityContribution(rLeftHandSideMatrix, rRightHandSideVector, rData);
            AddBoundaryConditionModifiedNitcheContribution(rLeftHandSideMatrix, rRightHandSideVector, rData);
        }
    }

    /**
    * This functions sets the B strain matrix (pressure columns are set to zero)
    * @param rDN_DX reference to the current Gauss pt. shape function gradients
    * @param rB_matrix reference to the computed B strain matrix
    */
    void SetInterfaceStrainMatrix(
        const BoundedMatrix<double, TNumNodes, TDim>& rDN_DX,
        BoundedMatrix<double, (TDim-1)*3, TNumNodes*(TDim+1)>& rB_matrix) {

        constexpr unsigned int BlockSize = TDim+1;
        rB_matrix.clear();

        if (TDim == 3) {
            for (unsigned int i=0; i<TNumNodes; i++) {
                rB_matrix(0,i*BlockSize)   = rDN_DX(i,0);
                rB_matrix(1,i*BlockSize+1) = rDN_DX(i,1);
                rB_matrix(2,i*BlockSize+2) = rDN_DX(i,2);
                rB_matrix(3,i*BlockSize)   = rDN_DX(i,1);
                rB_matrix(3,i*BlockSize+1) = rDN_DX(i,0);
                rB_matrix(4,i*BlockSize+1) = rDN_DX(i,2);
                rB_matrix(4,i*BlockSize+2) = rDN_DX(i,1);
                rB_matrix(5,i*BlockSize)   = rDN_DX(i,2);
                rB_matrix(5,i*BlockSize+2) = rDN_DX(i,0);
            }
        } else {
            for (unsigned int i=0; i<TNumNodes; i++) {
                rB_matrix(0,i*BlockSize)   = rDN_DX(i,0);
                rB_matrix(1,i*BlockSize+1) = rDN_DX(i,1);
                rB_matrix(2,i*BlockSize)   = rDN_DX(i,1);
                rB_matrix(2,i*BlockSize+1) = rDN_DX(i,0);
            }
        }
    }

    /**
    * This functions sets the normal projection matrix nxn
    * @param rUnitNormal reference to Gauss pt. unit normal vector
    * @param rNormProjMatrix reference to the computed normal projection matrix
    */
    void SetNormalProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        BoundedMatrix<double, TDim, TDim>& rNormProjMatrix) {

        rNormProjMatrix.clear();

        // Fill the normal projection matrix (nxn)
        if (TDim == 3) {
            noalias(rNormProjMatrix) = outer_prod(rUnitNormal, rUnitNormal);
        } else {
            rNormProjMatrix(0,0) = rUnitNormal(0)*rUnitNormal(0);
            rNormProjMatrix(0,1) = rUnitNormal(0)*rUnitNormal(1);
            rNormProjMatrix(1,0) = rUnitNormal(1)*rUnitNormal(0);
            rNormProjMatrix(1,1) = rUnitNormal(1)*rUnitNormal(1);
        }
    }

    /**
    * This functions sets the tangential projection matrix I - nxn
    * @param rUnitNormal reference to Gauss pt. unit normal vector
    * @param rTangProjMatrix reference to the computed tangential projection matrix
    */
    void SetTangentialProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        BoundedMatrix<double, TDim, TDim>& rTangProjMatrix) {

        rTangProjMatrix.clear();

        // Fill the tangential projection matrix (I - nxn)
        if (TDim == 3) {
            #ifdef KRATOS_USE_AMATRIX
            BoundedMatrix<double,3,3> id_matrix = IdentityMatrix(TDim);
            #else
            BoundedMatrix<double,3,3> id_matrix = IdentityMatrix(TDim,TDim);
            #endif
            noalias(rTangProjMatrix) = id_matrix - outer_prod(rUnitNormal, rUnitNormal);
        } else {
            rTangProjMatrix(0,0) = 1.0 - rUnitNormal(0)*rUnitNormal(0);
            rTangProjMatrix(0,1) = - rUnitNormal(0)*rUnitNormal(1);
            rTangProjMatrix(1,0) = - rUnitNormal(1)*rUnitNormal(0);
            rTangProjMatrix(1,1) = 1.0 - rUnitNormal(1)*rUnitNormal(1);
        }
    }

    /**
    * This functions sets the auxiliar matrix to compute the normal projection in Voigt notation
    * @param rUnitNormal reference to Gauss pt. unit normal vector
    * @param rVoigtNormProjMatrix reference to the computed normal projection auxiliar matrix
    */
    void SetVoigtNormalProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        BoundedMatrix<double, TDim, (TDim-1)*3>& rVoigtNormProjMatrix) {

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
    * This functions sets a vector containing the element previous solution
    * @param rData reference to the element data structure
    * @param rPrevSolVector reference to the previous solution vector
    */
    void GetPreviousSolutionVector(
        const ElementDataType& rData,
        array_1d<double, TNumNodes*(TDim+1)>& rPrevSolVector) {

        rPrevSolVector.clear();

        for (unsigned int i=0; i<TNumNodes; i++) {
            for (unsigned int comp=0; comp<TDim; comp++) {
                rPrevSolVector(i*(TDim+1)+comp) = rData.v(i,comp);
            }
            rPrevSolVector(i*(TDim+1)+TDim) = rData.p(i);
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

    void save(Serializer& rSerializer) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, NavierStokes<TDim>);
    }

    void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, NavierStokes<TDim>);
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

#endif // KRATOS_EMBEDDED_NAVIER_STOKES  defined
