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

#ifndef KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES_WALL_CONDITION_H
#define KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES_WALL_CONDITION_H

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
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

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
template <unsigned int TDim, unsigned int TNumNodes = TDim>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) EmbeddedAusasNavierStokesWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef GeometryType::Pointer GeometryPointerType;

    typedef GeometryType::PointsArrayType NodesArrayType;

    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    typedef Element::WeakPointer ElementWeakPointerType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    /// Pointer definition of EmbeddedAusasNavierStokesWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedAusasNavierStokesWallCondition);

    struct ConditionDataStruct
    {
        double charVel;                // Problem characteristic velocity (used in the outlet inflow prevention)
        double delta;                  // Non-dimensional positive sufficiently small constant (used in the outlet inflow prevention)

        // Data required in the RHS and LHS calculation
        double wGauss;                              // Gauss point weight
        array_1d<double, 3> Normal;                 // Condition normal
        array_1d<double, TNumNodes> N;              // Gauss point shape functions values
        BoundedMatrix<double, TNumNodes, TDim> v;  // Current step velocity

        // Data containers for no-split faces
        MatrixType N_container;
        VectorType w_gauss_container;
        std::vector<VectorType> area_normals_container;

        // Data containers for split faces
        // Positive face geometry data
        MatrixType N_pos_face;                          // Positive interface Gauss pts. shape functions values
        ShapeFunctionsGradientsType DN_DX_pos_face;     // Positive interface Gauss pts. shape functions gradients values
        VectorType w_gauss_pos_face;                    // Positive interface Gauss pts. weights
        std::vector<VectorType> pos_face_area_normals;  // Positive interface unit normal vector in each Gauss pt.

        // Negative face geometry data
        MatrixType N_neg_face;                          // Positive interface Gauss pts. shape functions values
        ShapeFunctionsGradientsType DN_DX_neg_face;     // Positive interface Gauss pts. shape functions gradients values
        VectorType w_gauss_neg_face;                    // Positive interface Gauss pts. weights
        std::vector<VectorType> neg_face_area_normals;  // Positive interface unit normal vector in each Gauss pt.

        unsigned int n_pos = 0;     // Number of positive distance nodes
        unsigned int n_neg = 0;     // Number of negative distance nodes
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
     @param NewId Index for the new
     */
    EmbeddedAusasNavierStokesWallCondition(IndexType NewId = 0) : Condition(NewId) {}

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    EmbeddedAusasNavierStokesWallCondition(IndexType NewId, const NodesArrayType& ThisNodes) :
        Condition(NewId,ThisNodes) {}

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    EmbeddedAusasNavierStokesWallCondition(IndexType NewId, GeometryType::Pointer pGeometry) :
        Condition(NewId,pGeometry) {}

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    EmbeddedAusasNavierStokesWallCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Condition(NewId,pGeometry,pProperties) {}

    /// Copy constructor.
    EmbeddedAusasNavierStokesWallCondition(EmbeddedAusasNavierStokesWallCondition const& rOther) :
        Condition(rOther) {}

    /// Destructor.
    ~EmbeddedAusasNavierStokesWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    EmbeddedAusasNavierStokesWallCondition & operator=(EmbeddedAusasNavierStokesWallCondition const& rOther) {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new EmbeddedAusasNavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override {
        return Kratos::make_shared<EmbeddedAusasNavierStokesWallCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new EmbeddedAusasNavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override {
        return Kratos::make_shared< EmbeddedAusasNavierStokesWallCondition >(NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override {
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    /**
     * If the condition is split, finds the condition parent element.
     * Note that this needs to be done at each time step for that cases
     * in where the distance function varies.
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override {
        KRATOS_TRY;

        // Set a reference to the current condition geometry
        GeometryType &r_geometry = this->GetGeometry();

        // Check if the condition is split
        unsigned int n_pos = 0, n_neg = 0;
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
            const double aux_dist = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
            if (aux_dist < 0) {
                n_neg++;
            } else {
                n_pos++;
            }
        }

        // If the condition is split, save a pointer to its parent element
        if (n_pos != 0 && n_neg != 0){

            // Get all the possible element candidates
            WeakPointerVector<Element> element_candidates;
            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
                WeakPointerVector<Element> &r_node_element_candidates = r_geometry[i_node].GetValue(NEIGHBOUR_ELEMENTS);
                for (unsigned int j = 0; j < r_node_element_candidates.size(); j++) {
                    element_candidates.push_back(r_node_element_candidates(j));
                }
            }

            // Check that the condition has candidate parent elements
            KRATOS_ERROR_IF(element_candidates.size() == 0) <<
                "Condition " << this->Id() << " has no candidate parent elements.\n" <<
                "Check that the FindNodalNeighboursProcess has been executed.";

            // Get a sort array with the current condition nodal ids
            std::vector<unsigned int> node_ids(TNumNodes), element_nodes_ids;

            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
                node_ids[i_node] = r_geometry[i_node].Id();
            }
            std::sort(node_ids.begin(), node_ids.end());

            // Iterate the candidate elements
            for (unsigned int i_candidate = 0; i_candidate < element_candidates.size(); ++i_candidate) {
                GeometryType &r_elem_geom = element_candidates[i_candidate].GetGeometry();
                const unsigned int n_elem_nodes = r_elem_geom.PointsNumber();

                // Get a sort array with the iterated candidate element nodal ids
                element_nodes_ids.resize(n_elem_nodes);
                for (unsigned int j = 0; j < n_elem_nodes; ++j) {
                    element_nodes_ids[j] = r_elem_geom[j].Id();
                }
                std::sort(element_nodes_ids.begin(), element_nodes_ids.end());

                // Check if the current condition ids are included in the iterated candidate element nodal ids
                if (std::includes(element_nodes_ids.begin(), element_nodes_ids.end(), node_ids.begin(), node_ids.end())) {
                    // Save a pointer to the parent element
                    mpParentElement = element_candidates(i_candidate);

                    // Save the parent element local ids. corresponding to the condition nodes
                    mParentElementIds.resize(TNumNodes);

                    std::vector<unsigned int > aux_elem_ids(n_elem_nodes);
                    for (unsigned int j = 0; j < n_elem_nodes; ++j) {
                        aux_elem_ids[j] = r_elem_geom[j].Id();
                    }

                    for (unsigned int i = 0; i < TNumNodes; ++i) {
                        const unsigned int aux_id = r_geometry[i].Id();
                        const std::vector<unsigned int >::iterator aux_it = std::find(aux_elem_ids.begin(), aux_elem_ids.end(), aux_id);
                        mParentElementIds[i] = std::distance(aux_elem_ids.begin(), aux_it);
                    }

                    // Leave the parent element search
                    return;
                }
            }

            KRATOS_ERROR << "Condition " << this->Id() << " cannot find parent element.";
        }

        KRATOS_CATCH("Error in EmbeddedAusasNavierStokesWallCondition InitializeSolutionStep() method.");
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
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false);

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false);

        // Struct to pass around the data
        ConditionDataStruct data;
        this->FillConditionData(data);

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_gauss;
        BoundedMatrix<double,MatrixSize, MatrixSize> lhs_gauss;

        // LHS and RHS contributions initialization
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Store the outlet inflow prevention constants in the data structure
        data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
        const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // const to avoid race conditions on data_value_container access/initialization
        data.charVel = rProcessInfo[CHARACTERISTIC_VELOCITY];

        // Loop on gauss points
        if (data.n_pos != 0 && data.n_neg != 0){

            // Positive side Gauss pts. loop
            const unsigned int n_gauss_pos = (data.w_gauss_pos_face).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss_pos; ++i_gauss) {

                Vector aux_N = row(data.N_pos_face, i_gauss);
                for (unsigned int i = 0; i < TNumNodes; ++i) {
                    data.N(i) = aux_N(mParentElementIds[i]);
                }
                data.wGauss = data.w_gauss_pos_face(i_gauss);
                data.Normal = data.pos_face_area_normals[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointRHSContribution(rhs_gauss, data);
                ComputeGaussPointLHSContribution(lhs_gauss, data);

                noalias(rLeftHandSideMatrix) += lhs_gauss;
                noalias(rRightHandSideVector) += rhs_gauss;
            }

            // Negative side Gauss pts. loop
            const unsigned int n_gauss_neg = (data.w_gauss_neg_face).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss_neg; ++i_gauss) {

                Vector aux_N = row(data.N_neg_face, i_gauss);
                for (unsigned int i = 0; i < TNumNodes; ++i) {
                    data.N(i) = aux_N(mParentElementIds[i]);
                }
                data.wGauss = data.w_gauss_neg_face(i_gauss);
                data.Normal = data.neg_face_area_normals[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointRHSContribution(rhs_gauss, data);
                ComputeGaussPointLHSContribution(lhs_gauss, data);

                noalias(rLeftHandSideMatrix) += lhs_gauss;
                noalias(rRightHandSideVector) += rhs_gauss;
            }
        } else {
            const unsigned int n_gauss = (data.w_gauss_container).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

                data.N = row(data.N_container, i_gauss);
                data.wGauss = data.w_gauss_container(i_gauss);
                data.Normal = data.area_normals_container[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointRHSContribution(rhs_gauss, data);
                ComputeGaussPointLHSContribution(lhs_gauss, data);

                noalias(rLeftHandSideMatrix) += lhs_gauss;
                noalias(rRightHandSideVector) += rhs_gauss;
            }
        }

        KRATOS_CATCH("")
    }

    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;
        this->FillConditionData(data);

        // Allocate memory needed
        BoundedMatrix<double, MatrixSize, MatrixSize> lhs_gauss;

        // LHS contributions initialization
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize, MatrixSize);

        // Loop on gauss points
        if (data.n_pos != 0 && data.n_neg != 0){

            // Positive side Gauss pts. loop
            const unsigned int n_gauss_pos = (data.w_gauss_pos_face).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss_pos; ++i_gauss) {

                Vector aux_N = row(data.N_pos_face, i_gauss);
                for (unsigned int i = 0; i < TNumNodes; ++i) {
                    data.N(i) = aux_N(mParentElementIds[i]);
                }
                data.wGauss = data.w_gauss_pos_face(i_gauss);
                data.Normal = data.pos_face_area_normals[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointLHSContribution(lhs_gauss, data);

                noalias(rLeftHandSideMatrix) += lhs_gauss;
            }

            // Negative side Gauss pts. loop
            const unsigned int n_gauss_neg = (data.w_gauss_neg_face).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss_neg; ++i_gauss) {

                Vector aux_N = row(data.N_neg_face, i_gauss);
                for (unsigned int i = 0; i < TNumNodes; ++i) {
                    data.N(i) = aux_N(mParentElementIds[i]);
                }
                data.wGauss = data.w_gauss_neg_face(i_gauss);
                data.Normal = data.neg_face_area_normals[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointLHSContribution(lhs_gauss, data);

                noalias(rLeftHandSideMatrix) += lhs_gauss;
            }
        } else {
            const unsigned int n_gauss = (data.w_gauss_container).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

                data.N = row(data.N_container, i_gauss);
                data.wGauss = data.w_gauss_container(i_gauss);
                data.Normal = data.area_normals_container[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointLHSContribution(lhs_gauss, data);

                noalias(rLeftHandSideMatrix) += lhs_gauss;
            }
        }

        KRATOS_CATCH("")
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

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;
        this->FillConditionData(data);

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_gauss;

        // RHS contributions initialization
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Store the outlet inflow prevention constants in the data structure
        data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
        const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // const to avoid race conditions on data_value_container access/initialization
        data.charVel = rProcessInfo[CHARACTERISTIC_VELOCITY];

        // Loop on gauss points
        if (data.n_pos != 0 && data.n_neg != 0){

            // Positive side Gauss pts. loop
            const unsigned int n_gauss_pos = (data.w_gauss_pos_face).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss_pos; ++i_gauss) {

                Vector aux_N = row(data.N_pos_face, i_gauss);
                for (unsigned int i = 0; i < TNumNodes; ++i) {
                    data.N(i) = aux_N(mParentElementIds[i]);
                }
                data.wGauss = data.w_gauss_pos_face(i_gauss);
                data.Normal = data.pos_face_area_normals[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointRHSContribution(rhs_gauss, data);

                noalias(rRightHandSideVector) += rhs_gauss;
            }

            // Negative side Gauss pts. loop
            const unsigned int n_gauss_neg = (data.w_gauss_neg_face).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss_neg; ++i_gauss) {

                Vector aux_N = row(data.N_neg_face, i_gauss);
                for (unsigned int i = 0; i < TNumNodes; ++i) {
                    data.N(i) = aux_N(mParentElementIds[i]);
                }
                data.wGauss = data.w_gauss_neg_face(i_gauss);
                data.Normal = data.neg_face_area_normals[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointRHSContribution(rhs_gauss, data);

                noalias(rRightHandSideVector) += rhs_gauss;
            }
        } else {
            const unsigned int n_gauss = (data.w_gauss_container).size();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

                data.N = row(data.N_container, i_gauss);
                data.wGauss = data.w_gauss_container(i_gauss);
                data.Normal = data.area_normals_container[i_gauss];
                data.Normal /= norm_2(data.Normal); // Normalize the area normal

                ComputeGaussPointRHSContribution(rhs_gauss, data);

                noalias(rRightHandSideVector) += rhs_gauss;
            }
        }

        KRATOS_CATCH("")
    }


    /// Condition check
    /**
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const GeometryType& rGeom = this->GetGeometry();

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0) {
            return Check;
        } else {
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
            if(DISTANCE.Key() == 0)
                KRATOS_ERROR << "DISTANCE Key is 0. Check if the application was correctly registered.";

            // Checks on nodes
            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i = 0; i < rGeom.size(); ++i) {
                if(rGeom[i].SolutionStepsDataHas(VELOCITY) == false)
                    KRATOS_ERROR << "Missing VELOCITY variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(PRESSURE) == false)
                    KRATOS_ERROR << "Missing PRESSURE variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                    KRATOS_ERROR << "Missing MESH_VELOCITY variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(ACCELERATION) == false)
                    KRATOS_ERROR << "Missing ACCELERATION variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(EXTERNAL_PRESSURE) == false)
                    KRATOS_ERROR << "Missing EXTERNAL_PRESSURE variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(DISTANCE) == false)
                    KRATOS_ERROR << "Missing DISTANCE variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].HasDofFor(VELOCITY_X) == false || rGeom[i].HasDofFor(VELOCITY_Y) == false || rGeom[i].HasDofFor(VELOCITY_Z) == false)
                    KRATOS_ERROR << "Missing VELOCITY component degree of freedom on node " << rGeom[i].Id();
                if(rGeom[i].HasDofFor(PRESSURE) == false)
                    KRATOS_ERROR << "Missing PRESSURE component degree of freedom on node " << rGeom[i].Id();
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
    std::string Info() const override {
        std::stringstream buffer;
        buffer << "EmbeddedAusasNavierStokesWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info() << "\nCondition id: " << Id();
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

    Element::Pointer pGetElement() {
		return mpParentElement.lock();
	}

    std::vector<unsigned int > GetParentElementIds() {
        return mParentElementIds;
    }

    void ComputeGaussPointLHSContribution(BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs, const ConditionDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);

    void ComputeRHSNeumannContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);
    void ComputeRHSOutletInflowContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& rData);

    // Auxiliar function to fill the element data structure
    void FillConditionData(ConditionDataStruct &rData)
    {
        const GeometryType& r_geometry = this->GetGeometry();

        // Check if the condition is split
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
            const double aux_dist = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
            if (aux_dist < 0) {
                rData.n_neg++;
            } else {
                rData.n_pos++;
            }
        }

        // If the element is split, take the values from the parent element modified shape functions utility
        // Otherwise, take the values from the current condition geometry
        if (rData.n_pos != 0 && rData.n_neg != 0){
            // Get the parent element nodal distances
            Element::Pointer p_parent_element = this->pGetElement();
            GeometryPointerType p_parent_geometry = p_parent_element->pGetGeometry();
            const Vector &distances = p_parent_element->GetValue(ELEMENTAL_DISTANCES);
            const unsigned int n_parent_nodes = p_parent_geometry->PointsNumber();

            // Construct the modified shape functions utility with the parent element pointer
            ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = nullptr;
            if (n_parent_nodes == 4) {
                p_ausas_modified_sh_func = Kratos::make_shared<Tetrahedra3D4AusasModifiedShapeFunctions>(p_parent_geometry, distances);
            }
            else if (n_parent_nodes == 3) {
                p_ausas_modified_sh_func = Kratos::make_shared<Triangle2D3AusasModifiedShapeFunctions>(p_parent_geometry, distances);
            } else {
                KRATOS_ERROR << "Asking for a non-implemented geometry modified shape functions utility.";
            }

            // Get the current condition global ids
            std::vector<unsigned int> cond_ids(TNumNodes);
            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
                cond_ids[i_node] = r_geometry[i_node].Id();
            }
            std::sort(cond_ids.begin(), cond_ids.end());

            DenseMatrix<unsigned int> elem_face_loc_ids;
            DenseVector<unsigned int> elem_nodes_in_face;
            p_parent_geometry->NodesInFaces(elem_face_loc_ids);
            p_parent_geometry->NumberNodesInFaces(elem_nodes_in_face);
            const unsigned int n_elem_faces = elem_face_loc_ids.size2();

            // Iterate the element faces to find the condition correspondent one
            unsigned int face_id = n_elem_faces + 1;
            for (unsigned int i_face = 0; i_face < n_elem_faces; ++i_face) {
                const unsigned int n_face_nodes = elem_nodes_in_face(i_face);

                // Get the element face local nodal ids
                // Note that the first index represent the node out of the face
                // (check this in case quads or tetras are used).
                std::vector<unsigned int> face_loc_ids(n_face_nodes);
                for (unsigned int i_node = 0; i_node < n_face_nodes; ++i_node) {
                    face_loc_ids[i_node] = elem_face_loc_ids(i_node + 1, i_face);
                }

                // Get the element face global nodal ids
                std::vector<unsigned int> face_glob_ids(n_face_nodes);
                for (unsigned int i_node = 0; i_node < n_face_nodes; ++i_node) {
                    const int aux_loc_id = face_loc_ids[i_node];
                    face_glob_ids[i_node] = (*p_parent_geometry)[aux_loc_id].Id();
                }
                std::sort(face_glob_ids.begin(), face_glob_ids.end());

                // Check if the element face global ids correspond with the current condition ones
                if (std::includes(face_glob_ids.begin(), face_glob_ids.end(), cond_ids.begin(), cond_ids.end())) {
                    face_id = i_face;
                    break;
                }
            }

            KRATOS_ERROR_IF(face_id == n_elem_faces + 1) <<
                "No parent element face found for condition " << this->Id() << " and parent element " << p_parent_element->Id();

            // Call the positive and negative sides modified shape functions face utilities
            p_ausas_modified_sh_func->ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                rData.N_pos_face,
                rData.DN_DX_pos_face,
                rData.w_gauss_pos_face,
                face_id,
                GeometryData::GI_GAUSS_2);

            p_ausas_modified_sh_func->ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                rData.N_neg_face,
                rData.DN_DX_neg_face,
                rData.w_gauss_neg_face,
                face_id,
                GeometryData::GI_GAUSS_2);

            p_ausas_modified_sh_func->ComputePositiveExteriorFaceAreaNormals(
                rData.pos_face_area_normals,
                face_id,
                GeometryData::GI_GAUSS_2);

            p_ausas_modified_sh_func->ComputeNegativeExteriorFaceAreaNormals(
                rData.neg_face_area_normals,
                face_id,
                GeometryData::GI_GAUSS_2);

        } else {
            // If the condition is not split, take the geometry shape function values
            GeometryType::IntegrationPointsArrayType integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
            const unsigned int n_gauss = integration_points.size();

            // Get the condition geometry shape functions values
            rData.N_container = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

            // Compute each Gauss pt. weight
            Vector gauss_pts_J_det(n_gauss);
            r_geometry.DeterminantOfJacobian(gauss_pts_J_det, GeometryData::GI_GAUSS_2);
            (rData.w_gauss_container).resize(n_gauss);
            for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                rData.w_gauss_container(i_gauss) = integration_points[i_gauss].Weight() * gauss_pts_J_det(i_gauss);
            }

            // Compute each Gauss pt. area normal
            (rData.area_normals_container).clear();
            for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                const CoordinatesArrayType& gauss_pt_loc_coords = integration_points[i_gauss].Coordinates();
                (rData.area_normals_container).push_back(r_geometry.Normal(gauss_pt_loc_coords));
            }
        }

        // Fill the nodal velocity array
        for (unsigned int i = 0; i < TNumNodes; i++) {
            const array_1d<double, 3> &vel = r_geometry[i].FastGetSolutionStepValue(VELOCITY);

            for (unsigned int k = 0; k < TDim; k++) {
                rData.v(i, k) = vel[k];
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

    ElementWeakPointerType mpParentElement;
    std::vector<unsigned int> mParentElementIds;

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

}; // Class EmbeddedAusasNavierStokesWallCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, EmbeddedAusasNavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const EmbeddedAusasNavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES_WALL_CONDITION_H
