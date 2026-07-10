//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined(KRATOS_CALCULATE_EMBEDDED_VARIABLE_FROM_SKIN_PROCESS_INCLUDED )
#define  KRATOS_CALCULATE_EMBEDDED_VARIABLE_FROM_SKIN_PROCESS_INCLUDED

// System includes


// External includes


// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "includes/kratos_flags.h"
#include "factories/linear_solver_factory.h"
#include "elements/embedded_nodal_variable_calculation_element_simplex.h"
#include "processes/process.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/intersection_utilities.h"
#include "utilities/variable_utils.h"

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

template< class TVarType >
class EmbeddedNodalVariableFromSkinTypeHelperClass
{
public:

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of EmbeddedNodalVariableFromSkinTypeHelperClass
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedNodalVariableFromSkinTypeHelperClass);

    ///@}
    ///@name Life Cycle
    ///@{


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Get the Unknown Variable object
     * This method returns a reference to the unknown variable. For double embedded nodal
     * variables this is a reference to NODAL_MAUX. In case of array type embedded nodal
     * variables, it returns a reference to NODAL_VAUX. These are the variables used when
     * solving the embedded nodal values least squares minimization problem.
     * @return const Variable<TVarType>& Reference to the unknown variable
     */
    static inline const Variable<TVarType> &GetUnknownVariable();

    /**
     * @brief Add the unknown variable to a model part
     * This method adds the unknown variable to the model part of interest.
     * @param rModelPart Reference to the model part to which the variable is added
     */
    static inline void AddUnknownVariable(ModelPart &rModelPart);

    /**
     * @brief Add the unknown variable DOFs to a model part
     * This method adds the unknown variable DOFs to the model part of interest
     * @param rModelPart Reference to the model part to which the variable DOFs are added
     */
    static inline void AddUnknownVariableDofs(ModelPart &rModelPart);

    ///@}
};

template <>
inline const Variable<double> &EmbeddedNodalVariableFromSkinTypeHelperClass<double>::GetUnknownVariable()
{
    return KratosComponents<Variable<double>>::Get("NODAL_MAUX");
}

template <>
inline const Variable<array_1d<double,3>> &EmbeddedNodalVariableFromSkinTypeHelperClass<array_1d<double,3>>::GetUnknownVariable()
{
    return KratosComponents<Variable<array_1d<double, 3>>>::Get("NODAL_VAUX");
}

template <>
inline void EmbeddedNodalVariableFromSkinTypeHelperClass<double>::AddUnknownVariable(ModelPart &rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(NODAL_MAUX);
}

template <>
inline void EmbeddedNodalVariableFromSkinTypeHelperClass<array_1d<double,3>>::AddUnknownVariable(ModelPart &rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(NODAL_VAUX);
}

template <>
inline void EmbeddedNodalVariableFromSkinTypeHelperClass<double>::AddUnknownVariableDofs(ModelPart &rModelPart)
{
    VariableUtils().AddDof(NODAL_MAUX, rModelPart);
}

template <>
inline void EmbeddedNodalVariableFromSkinTypeHelperClass<array_1d<double, 3>>::AddUnknownVariableDofs(ModelPart &rModelPart)
{
    VariableUtils().AddDof(NODAL_VAUX_X, rModelPart);
    VariableUtils().AddDof(NODAL_VAUX_Y, rModelPart);
    VariableUtils().AddDof(NODAL_VAUX_Z, rModelPart);
}

template <class TVarType, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class CalculateEmbeddedNodalVariableFromSkinProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    typedef typename TLinearSolver::Pointer LinearSolverPointerType;
    typedef typename Scheme<TSparseSpace,TDenseSpace>::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverPointerType;
    typedef typename ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::UniquePointer SolvingStrategyPointerType;
    typedef typename FindIntersectedGeometricalObjectsProcess::UniquePointer FindIntersectedGeometricalObjectsProcessPointerType;

    typedef std::unordered_set<std::pair<std::size_t, std::size_t>, PairHasher<std::size_t, std::size_t>, PairComparor<std::size_t, std::size_t>> EdgesSetType;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of CalculateEmbeddedNodalVariableFromSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateEmbeddedNodalVariableFromSkinProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Calculate Embedded Nodal Variable From Skin Process object
     * Constructor with model and json settings
     * @param rModel Model container
     * @param rSettings Settings json string
     */
    CalculateEmbeddedNodalVariableFromSkinProcess(
        Model &rModel,
        Parameters rSettings)
        : CalculateEmbeddedNodalVariableFromSkinProcess(
            rModel.GetModelPart(rSettings["base_model_part_name"].GetString()),
            rModel.GetModelPart(rSettings["skin_model_part_name"].GetString()),
            [] (Parameters x) -> Parameters {x.ValidateAndAssignDefaults(StaticGetDefaultParameters()); return x;} (rSettings))
    {
    }

    /**
     * @brief Construct a new Calculate Embedded Nodal Variable From Skin Process object
     *
     * @param rBaseModelPart Background mesh model part reference
     * @param rSkinModelPart Embedded skin model part reference
     * @param LinearSolverSettings Linear solver json settings
     * @param rSkinVariable Skin variable to take the values from
     * @param rEmbeddedNodalVariable Background mesh destination variable
     * @param LevelSetType Level set type (continuous or discontinuous)
     * @param BufferPosition Position in the buffer to take and save the values
     * @param AuxPartName Auxiliary intersections model part name
     */
    CalculateEmbeddedNodalVariableFromSkinProcess(
        ModelPart &rBaseModelPart,
        ModelPart &rSkinModelPart,
        Parameters LinearSolverSettings,
        const Variable<TVarType> &rSkinVariable,
        const Variable<TVarType> &rEmbeddedNodalVariable,
        const double GradientPenaltyCoefficient = 0.0,
        const unsigned int BufferPosition = 0,
        const std::string& AuxPartName = "IntersectedElementsModelPart",
        const std::size_t EchoLevel = 0)
        : Process()
        , mEchoLevel(EchoLevel)
        , mBufferPosition(BufferPosition),
          mAuxModelPartName(AuxPartName),
          mGradientPenaltyCoefficient(GradientPenaltyCoefficient),
          mrBaseModelPart(rBaseModelPart),
          mrSkinModelPart(rSkinModelPart),
          mrSkinVariable(rSkinVariable),
          mrEmbeddedNodalVariable(rEmbeddedNodalVariable)
    {
        KRATOS_TRY

        // Check the process settings
        KRATOS_ERROR_IF(!(mBufferPosition < rBaseModelPart.GetBufferSize())) <<
            "Asked for buffer position " << mBufferPosition << " buf base model part buffer size is " << rBaseModelPart.GetBufferSize() << std::endl;
        KRATOS_ERROR_IF(!(mBufferPosition < rSkinModelPart.GetBufferSize())) <<
            "Asked for buffer position " << mBufferPosition << " buf skin model part buffer size is " << rSkinModelPart.GetBufferSize() << std::endl;

        // Check that there is at least one element and node in the model
        int n_loc_mesh_nodes = mrBaseModelPart.GetCommunicator().pLocalMesh()->NumberOfNodes();
        int n_loc_mesh_elements = mrBaseModelPart.GetCommunicator().pLocalMesh()->NumberOfElements();
        KRATOS_ERROR_IF(mrBaseModelPart.GetCommunicator().GetDataCommunicator().SumAll(n_loc_mesh_nodes) == 0) << "The base model part has no nodes." << std::endl;
        KRATOS_ERROR_IF(mrBaseModelPart.GetCommunicator().GetDataCommunicator().SumAll(n_loc_mesh_elements) == 0) << "The base model Part has no elements." << std::endl;

        // Check that the base model part is conformed by simplex elements
        const auto &r_aux_geom = (mrBaseModelPart.ElementsBegin())->GetGeometry();
        const unsigned int dim = r_aux_geom.WorkingSpaceDimension();
        if(dim == 2){
            KRATOS_ERROR_IF(r_aux_geom.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle) <<
                "In 2D the element type is expected to be a triangle." << std::endl;
        } else if(dim == 3) {
            KRATOS_ERROR_IF(r_aux_geom.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra) <<
                "In 3D the element type is expected to be a tetrahedron" << std::endl;
        } else {
            KRATOS_ERROR << "Wrong geometry WorkingSpaceDimension(). Expected 2 or 3 and obtained: " << dim;
        }

        // Construct the linear solver pointer
        LinearSolverFactory<TSparseSpace, TDenseSpace> linear_solver_factory;
        mpLinearSolver = linear_solver_factory.Create(LinearSolverSettings);

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~CalculateEmbeddedNodalVariableFromSkinProcess() override
    {
        Model& current_model = mrBaseModelPart.GetModel();
        if(current_model.HasModelPart(mAuxModelPartName)) {
            current_model.DeleteModelPart(mAuxModelPartName);
        }
    };

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    ModelPart &GetIntersectedEdgesModelPart() const
    {
        Model &current_model = mrBaseModelPart.GetModel();
        return current_model.GetModelPart(mAuxModelPartName);
    }

    void Execute() override
    {
        KRATOS_TRY;
        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        this->GenerateIntersectedEdgesElementsModelPart();
        
        // Set the linear strategy to solve the regression problem
        this->SetLinearStrategy();

        // Solve the regression problem
        mpSolvingStrategy->Solve();

        // Copy the obtained values from the unknown variable to the user-defined variable
        this->SetObtainedEmbeddedNodalValues();

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        Model& current_model = mrBaseModelPart.GetModel();
        ModelPart& r_intersected_edges_model_part = current_model.GetModelPart( mAuxModelPartName );
        r_intersected_edges_model_part.Nodes().clear();
        r_intersected_edges_model_part.Elements().clear();
        r_intersected_edges_model_part.Conditions().clear();

        mpSolvingStrategy->Clear();
    }

    /**
     * @brief Get the Default Settings object
     * This method returns the default parameters for this process.
     * Note that it is required to be static since it is called during
     * the construction of the object so no instantiation exists yet.
     * @return Parameters Default parameters json string
     */
    static Parameters StaticGetDefaultParameters()
    {
        Parameters default_settings(R"(
        {
            "echo_level" : 0,
            "base_model_part_name": "",
            "skin_model_part_name": "",
            "skin_variable_name": "",
            "embedded_nodal_variable_name": "",
            "buffer_position": 0,
            "gradient_penalty_coefficient": 0.0,
            "aux_model_part_name": "IntersectedElementsModelPart",
            "linear_solver_settings": {
                "preconditioner_type": "amg",
                "solver_type": "amgcl",
                "smoother_type": "ilu0",
                "krylov_type": "cg",
                "max_iteration": 1000,
                "verbosity": 0,
                "tolerance": 1e-8,
                "scaling": false,
                "block_size": 1,
                "use_block_matrices_if_possible": true
            }
        }
        )");

        return default_settings;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        return StaticGetDefaultParameters();
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
        return "CalculateEmbeddedNodalVariableFromSkinProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateEmbeddedNodalVariableFromSkinProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
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

    const std::size_t mEchoLevel;
    const unsigned int mBufferPosition;
    const std::string mAuxModelPartName;
    const double mGradientPenaltyCoefficient;

    ModelPart& mrBaseModelPart;
    ModelPart& mrSkinModelPart;

    const Variable<TVarType> &mrSkinVariable;
    const Variable<TVarType> &mrEmbeddedNodalVariable;

    LinearSolverPointerType mpLinearSolver = nullptr;
    SolvingStrategyPointerType mpSolvingStrategy = nullptr;

    FindIntersectedGeometricalObjectsProcessPointerType mpFindIntersectedGeometricalObjectsProcess;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    virtual void GenerateIntersectedEdgesElementsModelPart()
    {
        KRATOS_TRY

        // Compute element intersections
        this->CalculateIntersections();

        Model& current_model = mrBaseModelPart.GetModel();
        if(current_model.HasModelPart(mAuxModelPartName)) {
            current_model.DeleteModelPart(mAuxModelPartName);
        }

        // Generate the auxiliary model part
        ModelPart& r_int_elems_model_part = current_model.CreateModelPart(mAuxModelPartName);

        r_int_elems_model_part.Nodes().clear();
        r_int_elems_model_part.Elements().clear();
        r_int_elems_model_part.Conditions().clear();

        r_int_elems_model_part.SetBufferSize(1);
        r_int_elems_model_part.CreateNewProperties(0, 0);

        // Set the gradient penalty coefficient in the auxiliary model part process info
        r_int_elems_model_part.GetProcessInfo()[GRADIENT_PENALTY_COEFFICIENT] = mGradientPenaltyCoefficient;

        // Add the minimization problem auxiliary variables
        this->AddIntersectedElementsVariables(r_int_elems_model_part);

        // Add intersected elements
        this->AddIntersectedElementsModelPartElements(r_int_elems_model_part);

        // Add DOFs to intersected elements model part
        this->AddIntersectedElementsModelPartDOFs(r_int_elems_model_part);

        KRATOS_CATCH("")
    }

    void SetObtainedEmbeddedNodalValues() const
    {
        const auto &rUnknownVariable = EmbeddedNodalVariableFromSkinTypeHelperClass<TVarType>::GetUnknownVariable();
        const auto &r_int_elems_model_part = (mrBaseModelPart.GetModel()).GetModelPart(mAuxModelPartName);

        block_for_each(r_int_elems_model_part.Nodes(), [&](Node& rNode){
            auto &r_emb_nod_val = (mrBaseModelPart.GetNode(rNode.Id())).FastGetSolutionStepValue(mrEmbeddedNodalVariable, mBufferPosition);
            r_emb_nod_val = rNode.FastGetSolutionStepValue(rUnknownVariable);
        });
    }

    inline void AddIntersectedElementsVariables(ModelPart &rModelPart) const
    {
        EmbeddedNodalVariableFromSkinTypeHelperClass<TVarType>::AddUnknownVariable(rModelPart);
    }

    void AddIntersectedElementsModelPartDOFs(ModelPart &rModelPart) const
    {
        EmbeddedNodalVariableFromSkinTypeHelperClass<TVarType>::AddUnknownVariableDofs(rModelPart);
    }

    void AddIntersectedElementsModelPartElements(ModelPart &rModelPart) const
    {
        // Initialize the VISITED flag in the origin model part
        // It will be used to mark the nodes already added to the intersected elements model part
        VariableUtils().SetFlag(VISITED, false, mrBaseModelPart.Nodes());

        // Initialize the INTERFACE flag in the origin model part
        // It will be used to mark the elements that have any intersection with the skin model part
        VariableUtils().SetFlag(INTERFACE, false, mrBaseModelPart.Elements());

        // Create element edges map
        EdgesSetType edges_set;

        // Get the base model part intersections
        auto &r_int_obj_vect = mpFindIntersectedGeometricalObjectsProcess->GetIntersections();

        // Get the unknown variable from Kratos components
        const auto &rUnknownVariable = EmbeddedNodalVariableFromSkinTypeHelperClass<TVarType>::GetUnknownVariable();

        // Temporary container of nodes
        // This is intentionally done to add the nodes at once and avoid the sort at each CreateNewNode call
        std::unordered_map<unsigned int, Node::Pointer> map_of_nodes;

        // Loop the base model part elements
        std::size_t new_elem_id = 1;
        for (unsigned int i_elem = 0; i_elem < mrBaseModelPart.NumberOfElements(); ++i_elem) {
            auto it_elem = mrBaseModelPart.ElementsBegin() + i_elem;
            // Check if the current element has intersections
            if (r_int_obj_vect[i_elem].size() != 0) {
                // Initialize the element values
                auto &r_geom = it_elem->GetGeometry();
                const auto edges = r_geom.GenerateEdges();

                // Loop the edges
                for (unsigned int i_edge = 0; i_edge < r_geom.EdgesNumber(); ++i_edge) {
                    // Check if the current edge is already stored
                    auto &r_i_edge_geom = edges[i_edge];
                    auto i_edge_pair = this->SetEdgePair(r_i_edge_geom);

                    if (edges_set.find(i_edge_pair) == edges_set.end()) {
                        // Initialize edge values
                        double i_edge_d = 0.0; // Average normalized distance from lower id. node
                        unsigned int n_int_obj = 0; // Number edge of intersecting entities
                        TVarType i_edge_val = mrEmbeddedNodalVariable.Zero(); // Average edge variable value

                        // Check the edge intersection against all the candidates
                        for (auto &r_int_obj : r_int_obj_vect[i_elem]) {
                            Point intersection_point;
                            const bool is_intersected = this->ComputeEdgeIntersection(
                                r_int_obj.GetGeometry(),
                                r_i_edge_geom[0],
                                r_i_edge_geom[1],
                                intersection_point);

                            // Compute the variable value in the intersection point
                            if (is_intersected) {
                                n_int_obj++;
                                Vector int_obj_N;
                                array_1d<double,3> local_coords;
                                r_int_obj.GetGeometry().PointLocalCoordinates(local_coords, intersection_point);
                                r_int_obj.GetGeometry().ShapeFunctionsValues(int_obj_N, local_coords);
                                for (unsigned int i_node = 0; i_node < r_int_obj.GetGeometry().PointsNumber(); ++i_node) {
                                    i_edge_val += r_int_obj.GetGeometry()[i_node].FastGetSolutionStepValue(mrSkinVariable, mBufferPosition) * int_obj_N[i_node];
                                }
                                i_edge_d += intersection_point.Distance(r_i_edge_geom[0]) / r_i_edge_geom.Length();
                            }
                        }

                        // Check if the edge is intersected
                        if (n_int_obj != 0) {
                            // Flag the edge parent element if the edge is intersected by any entity
                            it_elem->Set(INTERFACE, true);

                            // Add the average edge value (there might exist cases in where
                            // more than one geometry intersects the edge of interest).
                            i_edge_d /= n_int_obj;
                            i_edge_val /= n_int_obj;

                            // If not added yet, add the edge nodes
                            this->AddEdgeNodes(r_i_edge_geom, rModelPart, map_of_nodes);

                            // Create a new element with the intersected edge geometry and fake properties
                            auto p_element = Kratos::make_intrusive<EmbeddedNodalVariableCalculationElementSimplex<TVarType>>(
                                new_elem_id,
                                this->pSetEdgeElementGeometry(map_of_nodes, r_i_edge_geom, i_edge_pair),
                                rModelPart.pGetProperties(0));

                            // Save the edge values in the new element
                            p_element->SetValue(DISTANCE, i_edge_d);
                            p_element->SetValue(rUnknownVariable, i_edge_val);

                            // Update the id. counter
                            new_elem_id++;

                            // Add the new edge element to the hash map
                            edges_set.insert(i_edge_pair);

                            // Add the new edge element to the intersected elements model part
                            rModelPart.Elements().push_back(p_element);
                        }
                    }
                }
            }
        }


        // Populate the modelpart with all the nodes in NodesMap
        // Note that a temporary vector is created from the set to add all nodes at once
        PointerVectorSet<Node> tmp;
        tmp.reserve(rModelPart.NumberOfElements()*2);
        for(auto& item: map_of_nodes){
            tmp.push_back(item.second);
        }
        rModelPart.AddNodes(tmp.begin(), tmp.end());    
    }

    void SetLinearStrategy()
    {
        // Create the linear strategy
        SchemePointerType p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>>();

        bool calculate_norm_dx = false;
        bool calculate_reactions = false;
        bool reform_dof_at_each_iteration = false;
        BuilderSolverPointerType p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>>(mpLinearSolver);

        Model &current_model = mrBaseModelPart.GetModel();
        ModelPart &r_aux_model_part = current_model.GetModelPart(mAuxModelPartName);

        mpSolvingStrategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>>(
            r_aux_model_part,
            p_scheme,
            p_builder_and_solver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_dx);

        mpSolvingStrategy->Check();
        mpSolvingStrategy->SetEchoLevel(mEchoLevel);
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

    /**
     * @brief Construct a new Calculate Embedded Nodal Variable From Skin Process object
     * Constructor with background and skin model parts as well as json settings. This
     * constructor is intentionally protected to avoid exposing it to the user since it
     * is intended to serve as an auxiliar constructor to bridge from the model and
     * parameters one, which checks the provided settings with the defaults, to the "old
     * fashioned" one. This allows keeping the member variables as const as well as to
     * have a unique implementation of the constructor required checks and operations.
     * @param rBaseModelPart Background mesh model part reference
     * @param rSkinModelPart Embedded skin model part reference
     * @param rSettings Settings json string
     */
    CalculateEmbeddedNodalVariableFromSkinProcess(
        ModelPart &rBaseModelPart,
        ModelPart &rSkinModelPart,
        Parameters rSettings)
        : CalculateEmbeddedNodalVariableFromSkinProcess(
            rBaseModelPart,
            rSkinModelPart,
            rSettings["linear_solver_settings"],
            KratosComponents<Variable<TVarType>>::Get(rSettings["skin_variable_name"].GetString()),
            KratosComponents<Variable<TVarType>>::Get(rSettings["embedded_nodal_variable_name"].GetString()),
            rSettings["gradient_penalty_coefficient"].GetDouble(),
            rSettings["buffer_position"].GetInt(),
            rSettings["aux_model_part_name"].GetString(),
            rSettings["echo_level"].GetInt())
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateIntersections()
    {
        mpFindIntersectedGeometricalObjectsProcess = Kratos::make_unique<FindIntersectedGeometricalObjectsProcess>(mrBaseModelPart, mrSkinModelPart);
        mpFindIntersectedGeometricalObjectsProcess->ExecuteInitialize();
        mpFindIntersectedGeometricalObjectsProcess->FindIntersections();
    }

    void ClearIntersections()
    {
        mpFindIntersectedGeometricalObjectsProcess->Clear();
    }

	bool ComputeEdgeIntersection(
		const Element::GeometryType& rIntObjGeometry,
		const Element::NodeType& rEdgePoint1,
		const Element::NodeType& rEdgePoint2,
		Point& rIntersectionPoint) const
	{
		bool intersection_flag = false;
		const unsigned int work_dim = rIntObjGeometry.WorkingSpaceDimension();
		if (work_dim == 2){
			const unsigned int intersection_status = IntersectionUtilities::ComputeLineLineIntersection<Element::GeometryType>(
				rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
            if (intersection_status == 1 || intersection_status == 3) {
                intersection_flag = true;
            }
		} else if (work_dim == 3){
			const unsigned int intersection_status = IntersectionUtilities::ComputeTriangleLineIntersection<Element::GeometryType>(
				rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
            if (intersection_status == 1) {
                intersection_flag = true;
            }
		} else {
			KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
		}

        return intersection_flag;
    }

    void AddEdgeNodes(
        const Geometry<Node> &rEdgeGeometry,
        ModelPart &rModelPart,
        std::unordered_map<unsigned int, Node::Pointer>& rNodesMap
        ) const
    {
        const auto& rp_var_list = rModelPart.pGetNodalSolutionStepVariablesList();
        unsigned int buffer_size = rModelPart.GetBufferSize();
        
        // Loop the edge nodes
        for (std::size_t i = 0; i < 2; ++i) {
            auto p_i_node = rEdgeGeometry(i);
            // Check if the node has been already added
            if (!p_i_node->Is(VISITED)) {
                p_i_node->Set(VISITED, true);
                auto p_node_copy = Kratos::make_intrusive< Node >(
                    p_i_node->Id(), 
                    p_i_node->Coordinates());
                p_node_copy->SetSolutionStepVariablesList(rp_var_list);
                p_node_copy->SetBufferSize(buffer_size);

                rNodesMap[p_i_node->Id()] = p_node_copy;
            }
        }
    }

    Element::GeometryType::Pointer pSetEdgeElementGeometry(
        std::unordered_map<unsigned int, Node::Pointer>& rNodesMap,
        const Element::GeometryType &rCurrentEdgeGeometry,
        const std::pair<std::size_t, std::size_t> NewEdgeIds) const
    {
        Element::GeometryType::PointsArrayType points_array;
        points_array.push_back(rNodesMap[std::get<0>(NewEdgeIds)]);
        points_array.push_back(rNodesMap[std::get<1>(NewEdgeIds)]);
        return rCurrentEdgeGeometry.Create(points_array);
    }

    inline std::pair<std::size_t, std::size_t> SetEdgePair(const Geometry<Node> &rEdgeGeom) const
    {
        std::pair<std::size_t, std::size_t> edge_pair(
            (rEdgeGeom[0].Id() < rEdgeGeom[1].Id()) ? rEdgeGeom[0].Id() : rEdgeGeom[1].Id(),
            (rEdgeGeom[0].Id() > rEdgeGeom[1].Id()) ? rEdgeGeom[0].Id() : rEdgeGeom[1].Id());
        return edge_pair;
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
    CalculateEmbeddedNodalVariableFromSkinProcess& operator=(CalculateEmbeddedNodalVariableFromSkinProcess const& rOther) = delete;

    /// Copy constructor.
    CalculateEmbeddedNodalVariableFromSkinProcess(CalculateEmbeddedNodalVariableFromSkinProcess const& rOther) = delete;

    ///@}
}; // Class CalculateEmbeddedNodalVariableFromSkinProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template< class TVarType, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator >> (
    std::istream& rIStream,
    CalculateEmbeddedNodalVariableFromSkinProcess<TVarType, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template< class TVarType, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CalculateEmbeddedNodalVariableFromSkinProcess<TVarType, TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_EMBEDDED_VARIABLE_FROM_SKIN_PROCESS_INCLUDED  defined
