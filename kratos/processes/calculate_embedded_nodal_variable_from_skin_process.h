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
#include "includes/kratos_flags.h"
#include "elements/embedded_nodal_variable_calculation_element_simplex.h"
#include "linear_solvers/cg_solver.h"
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
class KRATOS_API(KRATOS_CORE) CalculateEmbeddedNodalVariableFromSkinProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    typedef typename Element::NodeType::Pointer NodePointerType;
    typedef typename Scheme<TSparseSpace,TDenseSpace>::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverPointerType;
    typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::UniquePointer SolvingStrategyPointerType;
    typedef typename FindIntersectedGeometricalObjectsProcess::Pointer FindIntersectedGeometricalObjectsProcessPointerType;

    struct Hash
    {
        std::size_t operator()(const std::pair<unsigned int,bool>& k) const
        {
            std::size_t h1 = std::hash<unsigned int>()(std::get<0>(k));
            std::size_t h2 = std::hash<bool>()(std::get<1>(k));
            return h1 ^ (h2 << 1);
        }
    };

    struct KeyEqual
    {
        bool operator()(const std::pair<std::size_t, std::size_t>& lhs, const std::pair<std::size_t, std::size_t>& rhs) const
        {
            return ((std::get<0>(lhs) == std::get<0>(rhs)) && (std::get<1>(lhs) == std::get<1>(rhs)));
        }
    };

    typedef std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t, Hash, KeyEqual> EdgesMapType;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of CalculateEmbeddedNodalVariableFromSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateEmbeddedNodalVariableFromSkinProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    CalculateEmbeddedNodalVariableFromSkinProcess(
        ModelPart &rBaseModelPart,
        ModelPart &rSkinModelPart,
        const Variable<TVarType> &rSkinVariable,
        const Variable<TVarType> &rEmbeddedNodalVariable,
        const std::string LevelSetType = "continuous",
        const std::string AuxPartName = "IntersectedElementsModelPart")
        : mLevelSetType(LevelSetType),
          mAuxModelPartName(AuxPartName),
          mrBaseModelPart(rBaseModelPart),
          mrSkinModelPart(rSkinModelPart),
          mrSkinVariable(rSkinVariable),
          mrEmbeddedNodalVariable(rEmbeddedNodalVariable)
    {
        KRATOS_TRY

        // Check that there is at least one element and node in the model
        int n_loc_mesh_nodes = mrBaseModelPart.GetCommunicator().pLocalMesh()->NumberOfNodes();
        int n_loc_mesh_elements = mrBaseModelPart.GetCommunicator().pLocalMesh()->NumberOfElements();
        KRATOS_ERROR_IF(mrBaseModelPart.GetCommunicator().SumAll(n_loc_mesh_nodes) == 0) << "The base model part has no nodes." << std::endl;
        KRATOS_ERROR_IF(mrBaseModelPart.GetCommunicator().SumAll(n_loc_mesh_elements) == 0) << "The base model Part has no elements." << std::endl;

        // Check that the base model part is conformed by simplex elements
        const auto &r_aux_geom = (mrBaseModelPart.ElementsBegin())->GetGeometry();
        const unsigned int dim = r_aux_geom.Dimension();
        if(dim == 2){
            KRATOS_ERROR_IF(r_aux_geom.GetGeometryFamily() != GeometryData::Kratos_Triangle) <<
                "In 2D the element type is expected to be a triangle." << std::endl;
        } else if(dim == 3) {
            KRATOS_ERROR_IF(r_aux_geom.GetGeometryFamily() != GeometryData::Kratos_Tetrahedra) <<
                "In 3D the element type is expected to be a tetrahedron" << std::endl;
        } else {
            KRATOS_ERROR << "Wrong geometry Dimension(). Expected 2 or 3 and obtained: " << dim;
        }

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

    void Execute() override
    {
        KRATOS_TRY;

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        this->GenerateIntersectedEdgesElementsModelPart();

        // Set the linear strategy to solve the regression problem
        this->SetLinearStrategy();

        // Solve the regression problem
        mpSolvingStrategy->Solve();

        // Move the obtained values to the user-defined variable
        this->SetObtainedEmbeddedNodalValues();

        // Free the memory
        this->Clear();

        KRATOS_CATCH("")
    }

    virtual void Clear()
    {
        Model& current_model = mrBaseModelPart.GetModel();
        ModelPart& r_distance_model_part = current_model.GetModelPart( mAuxModelPartName );
        r_distance_model_part.Nodes().clear();
        r_distance_model_part.Elements().clear();
        r_distance_model_part.Conditions().clear();

        mpSolvingStrategy->Clear();
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

    const std::string mLevelSetType;
    const std::string mAuxModelPartName;

    ModelPart& mrBaseModelPart;
    ModelPart& mrSkinModelPart;

    const Variable<TVarType> &mrSkinVariable;
    const Variable<TVarType> &mrEmbeddedNodalVariable;

    SolvingStrategyPointerType mpSolvingStrategy;

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
        r_int_elems_model_part.SetProcessInfo(mrBaseModelPart.pGetProcessInfo());

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
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(r_int_elems_model_part.NumberOfNodes()); ++i_node) {
            const auto it_node = r_int_elems_model_part.NodesBegin() + i_node;
            auto &r_emb_nod_val = (mrBaseModelPart.GetNode(it_node->Id())).FastGetSolutionStepValue(mrEmbeddedNodalVariable);
            r_emb_nod_val = it_node->FastGetSolutionStepValue(rUnknownVariable);
        }
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

        // Create element edges map
        EdgesMapType edges_map;

        // Get the base model part intersections
        auto &r_int_obj_vect = mpFindIntersectedGeometricalObjectsProcess->GetIntersections();

        // Get the unknown variable from Kratos components
        const auto &rUnknownVariable = EmbeddedNodalVariableFromSkinTypeHelperClass<TVarType>::GetUnknownVariable();

        // Loop the base model part elements
        std::size_t new_elem_id = 1;
        for (unsigned int i_elem = 0; i_elem < mrBaseModelPart.NumberOfElements(); ++i_elem) {
            auto it_elem = mrBaseModelPart.ElementsBegin() + i_elem;
            const auto elem_dist = this->SetDistancesVector(it_elem);
            // Check if the current element is split
            if (IsSplit(elem_dist)) {
                if (r_int_obj_vect[i_elem].size() != 0) {
                    // Initialize the element values
                    auto &r_geom = it_elem->GetGeometry();
                    const auto edges = r_geom.Edges();

                    // Loop the edges
                    for (unsigned int i_edge = 0; i_edge < r_geom.EdgesNumber(); ++i_edge) {
                        // Check if the current edge is already stored
                        auto &r_i_edge_geom = edges[i_edge];
                        auto i_edge_pair = this->SetEdgePair(r_i_edge_geom);

                        if (edges_map.find(i_edge_pair) == edges_map.end()) {
                            // Initialize edge values
                            double i_edge_d = 0.0; // Average normalized distance from lower id. node
                            unsigned int n_int_obj = 0; // Number edge of intersecting entities
                            TVarType i_edge_val = mrEmbeddedNodalVariable.Zero(); // Average edge variable value

                            // Check the edge intersection against all the candidates
                            for (auto &r_int_obj : r_int_obj_vect[i_elem]) {
                                Point intersection_point;
                                const int is_intersected = this->ComputeEdgeIntersection(
                                    r_int_obj.GetGeometry(),
                                    r_i_edge_geom[0],
                                    r_i_edge_geom[1],
                                    intersection_point);

                                // Compute the variable value in the intersection point
                                if (is_intersected == 1) {
                                    n_int_obj++;
                                    Vector int_obj_N;
                                    array_1d<double,3> local_coords;
                                    r_int_obj.GetGeometry().PointLocalCoordinates(local_coords, intersection_point);
                                    r_int_obj.GetGeometry().ShapeFunctionsValues(int_obj_N, local_coords);
                                    for (unsigned int i_node = 0; i_node < r_int_obj.GetGeometry().PointsNumber(); ++i_node) {
                                        i_edge_val += r_int_obj.GetGeometry()[i_node].FastGetSolutionStepValue(mrSkinVariable) * int_obj_N[i_node];
                                    }
                                    i_edge_d = norm_2(intersection_point - r_i_edge_geom[0]) / r_i_edge_geom.Length();
                                }
                            }

                            // Check if the edge is intersected
                            if (n_int_obj != 0) {
                                // Add the average edge value (there might exist cases in where
                                // more than one geometry intersects the edge of interest).
                                i_edge_d /= n_int_obj;
                                i_edge_val /= n_int_obj;

                                // If not added yet, add the edge nodes
                                this->AddEdgeNodes(r_i_edge_geom, rModelPart);

                                // Create a new element with the intersected edge geometry and fake properties
                                auto p_element = Kratos::make_shared<EmbeddedNodalVariableCalculationElementSimplex<TVarType>>(
                                    new_elem_id,
                                    this->pSetEdgeElementGeometry(rModelPart, r_i_edge_geom, i_edge_pair),
                                    rModelPart.pGetProperties(0));

                                // Save the edge values in the new element
                                p_element->SetValue(DISTANCE, i_edge_d);
                                p_element->SetValue(rUnknownVariable, i_edge_val);

                                // Update the id. counter
                                new_elem_id++;

                                // Add the new edge element to the hash map
                                edges_map.insert(std::make_pair(i_edge_pair, new_elem_id));

                                // Add the new edge element to the intersected elements model part
                                rModelPart.Elements().push_back(p_element);
                            }
                        }
                    }
                }
            }
        }
    }

    void SetLinearStrategy()
    {
        // Create the CG linear solver (take advantage of the symmetry of the problem)
        Parameters linear_solver_parameters(R"({
            "solver_type": "cg_solver",
            "tolerance": 1.0e-8,
            "max_iteration": 200,
            "preconditioner_type": "none",
            "scaling": false
        })");
        auto p_linear_solver = Kratos::make_shared<CGSolver<TSparseSpace, TDenseSpace>>(linear_solver_parameters);

        // Create the linear strategy
        SchemePointerType p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>>();

        bool calculate_norm_dx = false;
        bool calculate_reactions = false;
        bool reform_dof_at_each_iteration = false;
        BuilderSolverPointerType p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>>(p_linear_solver);

        Model &current_model = mrBaseModelPart.GetModel();
        ModelPart &r_aux_model_part = current_model.GetModelPart(mAuxModelPartName);

        mpSolvingStrategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>>(
            r_aux_model_part,
            p_scheme,
            p_linear_solver,
            p_builder_and_solver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_dx);

        mpSolvingStrategy->Check();
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
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateIntersections()
    {
        mpFindIntersectedGeometricalObjectsProcess = Kratos::make_shared<FindIntersectedGeometricalObjectsProcess>(mrBaseModelPart, mrSkinModelPart);
        mpFindIntersectedGeometricalObjectsProcess->Initialize();
        mpFindIntersectedGeometricalObjectsProcess->FindIntersections();
    }

    void ClearIntersections()
    {
        mpFindIntersectedGeometricalObjectsProcess->Clear();
    }

    const Vector SetDistancesVector(ModelPart::ElementIterator ItElem) const
    {
        const auto &r_geom = ItElem->GetGeometry();
        Vector nodal_distances(r_geom.PointsNumber());

        if (mLevelSetType == "continuous"){
            // Continuous nodal distance function case
            for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                nodal_distances[i_node] = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
            }
        } else if (mLevelSetType == "discontinuous") {
            // Discontinuous elemental distance function case
            nodal_distances = ItElem->GetValue(ELEMENTAL_DISTANCES);
        } else {
            KRATOS_ERROR << "Level set type must be either 'continuous' or 'discontinuous'. Got " << mLevelSetType;
        }

        return nodal_distances;
    }

    inline bool IsSplit(const Vector &rDistances) const
    {
        unsigned int n_pos = 0, n_neg = 0;

        for (double dist : rDistances) {
            if(dist >= 0) {
                ++n_pos;
            } else {
                ++n_neg;
            }
        }

        if (n_pos > 0 && n_neg > 0) {
            return true;
        }

        return false;
    }

	int ComputeEdgeIntersection(
		const Element::GeometryType& rIntObjGeometry,
		const Element::NodeType& rEdgePoint1,
		const Element::NodeType& rEdgePoint2,
		Point& rIntersectionPoint) const
	{
		int intersection_flag = 0;
		const unsigned int work_dim = rIntObjGeometry.WorkingSpaceDimension();
		if (work_dim == 2){
			intersection_flag = IntersectionUtilities::ComputeLineLineIntersection<Element::GeometryType>(
				rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
		} else if (work_dim == 3){
			intersection_flag = IntersectionUtilities::ComputeTriangleLineIntersection<Element::GeometryType>(
				rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
		} else {
			KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
		}

        return intersection_flag;
    }

    void AddEdgeNodes(
        const Geometry<Node<3>> &rEdgeGeometry,
        ModelPart &rModelPart) const
    {
        // Loop the edge nodes
        for (std::size_t i = 0; i < 2; ++i) {
            auto p_i_node = rEdgeGeometry(i);
            // Check if the node has been already added
            if (!p_i_node->Is(VISITED)) {
                p_i_node->Set(VISITED, true);
                rModelPart.CreateNewNode(p_i_node->Id(), *p_i_node);
            }
        }
    }

    Element::GeometryType::Pointer pSetEdgeElementGeometry(
        ModelPart &rModelPart,
        const Element::GeometryType &rCurrentEdgeGeometry,
        const std::pair<std::size_t, std::size_t> NewEdgeIds) const
    {
        Element::GeometryType::PointsArrayType points_array;
        points_array.push_back(rModelPart.pGetNode(std::get<0>(NewEdgeIds)));
        points_array.push_back(rModelPart.pGetNode(std::get<1>(NewEdgeIds)));
        return rCurrentEdgeGeometry.Create(points_array);
    }

    inline std::pair<std::size_t, std::size_t> SetEdgePair(const Geometry<Node<3>> &rEdgeGeom) const
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
