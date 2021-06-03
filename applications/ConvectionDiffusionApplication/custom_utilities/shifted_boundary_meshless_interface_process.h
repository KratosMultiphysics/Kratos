// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_SHIFTED_BOUNDARY_INTERFACE_UTILITIES )
#define  KRATOS_SHIFTED_BOUNDARY_INTERFACE_UTILITIES

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "containers/pointer_vector.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "processes/process.h"
#include "utilities/element_size_calculator.h"
#include "utilities/mls_shape_functions_utility.h"

// Application includes

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

/// Short class definition.
class ShiftedBoundaryMeshlessInterfaceProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    using IndexType = ModelPart::IndexType;

    using NodeType = ModelPart::NodeType;

    using GeometryType = ModelPart::GeometryType;

    using ModifiedShapeFunctionsFactoryType = std::function<ModifiedShapeFunctions::UniquePointer(const GeometryType::Pointer, const Vector&)>;

    using MLSShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&, Matrix&)>;

    using ElementSizeFunctionType = std::function<double(const GeometryType&)>;

    using NodesCloudSetType = std::unordered_set<NodeType::Pointer, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of ShiftedBoundaryMeshlessInterfaceProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShiftedBoundaryMeshlessInterfaceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    ShiftedBoundaryMeshlessInterfaceProcess(
        Model& rModel,
        Parameters ThisParameters)
        : Process()
    {
        const std::string model_part_name = ThisParameters["model_part_name"].GetString();
        const std::string boundary_sub_model_part_name = ThisParameters["boundary_sub_model_part_name"].GetString();
        mpModelPart = &rModel.GetModelPart(model_part_name);
        if (mpModelPart->HasSubModelPart(boundary_sub_model_part_name)) {
            mpBoundarySubModelPart = &(mpModelPart->GetSubModelPart(boundary_sub_model_part_name));
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfNodes() != 0) << "Provided SBM model part has nodes." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfElements() != 0) << "Provided SBM model part has elements." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfConditions() != 0) << "Provided SBM model part has conditions." << std::endl;
        } else {
            mpBoundarySubModelPart = &(mpModelPart->CreateSubModelPart(boundary_sub_model_part_name));
        }
    };

    /// Destructor.
    virtual ~ShiftedBoundaryMeshlessInterfaceProcess() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute()
    {
        // Set the required interface flags
        SetInterfaceFlags();

        // Set the modified shape functions factory
        // Note that unique geometry in the mesh is assumed
        const auto& r_begin_geom = mpModelPart->ElementsBegin()->GetGeometry();
        auto p_mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_begin_geom);

        // Get the MLS shape functions function
        auto p_mls_sh_func = GetMLSShapeFunctionsFunction();

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        auto p_element_size_func = GetElementSizeFunction(r_begin_geom);

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mpModelPart->Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        // Loop the elements to find the intersected ones
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto p_geom = rElement.pGetGeometry();

            // If the element is split, get the interface Gauss pt. data
            // If the element is negative, it is deactivated to not be assembled
            if (IsSplit(*p_geom)) {
                // Set the meshless cloud of point support for the MLS
                Matrix cloud_nodes_coordinates;
                PointerVector<NodeType> cloud_nodes;
                SetSplitElementSupportCloud(rElement, cloud_nodes, cloud_nodes_coordinates);

                // Set up the distances vector
                const auto& r_geom = *p_geom;
                const std::size_t n_nodes = r_geom.PointsNumber();
                Vector nodal_distances(n_nodes);
                SetNodalDistancesVector(r_geom, nodal_distances);

                // Set the modified shape functions pointer and calculate the positive interface data
                auto p_mod_sh_func = p_mod_sh_func_factory(p_geom, nodal_distances);
                Vector pos_int_w;
                Matrix pos_int_N;
                std::vector<Vector> pos_int_n;
                typename ModifiedShapeFunctions::ShapeFunctionsGradientsType pos_int_DN_DX;
                //TODO: Add a method without the interface gradients
                p_mod_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(pos_int_N, pos_int_DN_DX, pos_int_w, GeometryData::GI_GAUSS_2);
                p_mod_sh_func->ComputePositiveSideInterfaceAreaNormals(pos_int_n, GeometryData::GI_GAUSS_2);

                // Calculate parent element size for the SBM BC imposition
                const double h = p_element_size_func(*p_geom);

                // Iterate the interface Gauss pts.
                DenseVector<double> i_g_N;
                array_1d<double,3> i_g_coords;
                const std::size_t n_int_pts = pos_int_w.size();
                for (std::size_t i_g = 0; i_g < n_int_pts; ++i_g) {
                    // Calculate Gauss pt. coordinates
                    i_g_N = row(pos_int_N, i_g);
                    noalias(i_g_coords) = ZeroVector(3);
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        noalias(i_g_coords) += i_g_N[i_node] * r_geom[i_node].Coordinates();
                    }

                    // Create a new condition with a geometry made up with the basis nodes
                    auto p_prop = rElement.pGetProperties();
                    auto p_cond = Kratos::make_intrusive<LaplacianShiftedBoundaryCondition>(++max_cond_id, cloud_nodes);
                    p_cond->SetProperties(p_prop); //TODO: Think if we want properties in these conditions or not
                    p_cond->Set(ACTIVE, true);
                    mpBoundarySubModelPart->AddCondition(p_cond);

                    // Store the SBM BC data in the condition database
                    p_cond->SetValue(ELEMENT_H, h);
                    p_cond->SetValue(NORMAL, pos_int_n[i_g]);
                    p_cond->SetValue(INTEGRATION_WEIGHT, pos_int_w[i_g]);
                    p_cond->SetValue(INTEGRATION_COORDINATES, i_g_coords);

                    // Calculate the MLS shape functions and gradients and save in the database
                    Vector N_container;
                    Matrix DN_DX_container;
                    const double mls_kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, i_g_coords);
                    p_mls_sh_func(cloud_nodes_coordinates, i_g_coords, mls_kernel_rad, N_container, DN_DX_container);
                    //FIXME: Find variables for these
                    p_cond->SetValue(BDF_COEFFICIENTS, N_container);
                    p_cond->SetValue(LOCAL_AXES_MATRIX, DN_DX_container);
                }
            }
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
    virtual std::string Info() const {
        return "ShiftedBoundaryMeshlessInterfaceProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {
        rOStream << "ShiftedBoundaryMeshlessInterfaceProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
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

    ///@}
    ///@name Protected Operations
    ///@{

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

    ModelPart* mpModelPart = nullptr;
    ModelPart* mpBoundarySubModelPart = nullptr;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void SetInterfaceFlags()
    {
        // Initialize flags to false
        block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
            rNode.Set(ACTIVE, false);
            rNode.Set(BOUNDARY, false);
            rNode.Set(INTERFACE, false);
        });
        block_for_each(mpModelPart->Elements(), [](Element& rElement){
            rElement.Set(ACTIVE, false);
            rElement.Set(BOUNDARY, false);
            rElement.Set(INTERFACE, false);
        });

        // Find active regions
        for (auto& rElement : mpModelPart->Elements()) {
            auto& r_geom = rElement.GetGeometry();
            if (IsSplit(r_geom)) {
                // Mark the intersected elements as BOUNDARY
                // The intersected elements are also flagged as non ACTIVE to avoid assembling them
                // Note that the split elements BC is applied by means of the extension operators
                rElement.Set(ACTIVE, false);
                rElement.Set(BOUNDARY, true);
                // Mark the positive intersected nodes as BOUNDARY
                for (auto& rNode : r_geom) {
                    if (!(rNode.FastGetSolutionStepValue(DISTANCE) < 0.0)) {
                        rNode.Set(BOUNDARY, true);
                    }
                }
            } else if (IsNegative(r_geom)) {
                // Mark the negative element as non ACTIVE to avoid assembling it
                rElement.Set(ACTIVE, false);
            } else {
                // Mark the positive element as ACTIVE to assemble it
                rElement.Set(ACTIVE, true);
                // Mark the active element nodes as ACTIVE
                for (auto& rNode : r_geom) {
                    rNode.Set(ACTIVE, true);
                }
            }
        }

        // Find the surrogate boundary
        // Note that we rely on the fact that the neighbours are sorted according to the faces
        for (auto& rElement : mpModelPart->Elements()) {
            if (rElement.Is(BOUNDARY)) {
                const auto& r_geom = rElement.GetGeometry();
                const std::size_t n_faces = r_geom.FacesNumber();
                auto& r_neigh_elems = rElement.GetValue(NEIGHBOUR_ELEMENTS);
                for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
                    // The neighbour corresponding to the current face is ACTIVE means that the current face is surrogate boundary
                    // Flag the current neighbour owning the surrogate face as INTERFACE as well as its nodes, which will form the MLS cloud
                    auto& r_neigh_elem = r_neigh_elems[i_face];
                    if (r_neigh_elem.Is(ACTIVE)) {
                        r_neigh_elem.Set(INTERFACE, true);
                        auto& r_neigh_geom = r_neigh_elem.GetGeometry();
                        for (auto& rNode : r_neigh_geom) {
                            rNode.Set(INTERFACE, true);
                        }
                    }
                }
            }
        }
    }

    bool IsSplit(const GeometryType& rGeometry)
    {
        std::size_t n_neg = 0;
        std::size_t n_pos = 0;
        for (const auto& r_node : rGeometry) {
            if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
                n_neg++;
            } else {
                n_pos++;
            }
        }
        return (n_pos != 0 && n_neg != 0);
    }

    bool IsNegative(const GeometryType& rGeometry)
    {
        std::size_t n_neg = 0;
        for (const auto& r_node : rGeometry) {
            if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
                n_neg++;
            }
        }
        return (n_neg == rGeometry.PointsNumber());
    }

    void SetNodalDistancesVector(
        const GeometryType& rGeometry,
        Vector& rNodalDistances)
    {
        const std::size_t n_nodes = rGeometry.PointsNumber();
        if (rNodalDistances.size() != n_nodes) {
            rNodalDistances.resize(n_nodes);
        }
        std::size_t i = 0;
        for (const auto& r_node : rGeometry) {
            rNodalDistances[i++] = r_node.FastGetSolutionStepValue(DISTANCE);
        }
    }

    ModifiedShapeFunctionsFactoryType GetStandardModifiedShapeFunctionsFactory(const GeometryType& rGeometry)
    {
        switch (rGeometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return [](const GeometryType::Pointer pGeometry, const Vector& rNodalDistances)->ModifiedShapeFunctions::UniquePointer{
                    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rNodalDistances);};
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return [](const GeometryType::Pointer pGeometry, const Vector& rNodalDistances)->ModifiedShapeFunctions::UniquePointer{
                    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rNodalDistances);};
            default:
                KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
        }
    }

    MLSShapeFunctionsFunctionType GetMLSShapeFunctionsFunction()
    {
        switch (mpModelPart->GetProcessInfo()[DOMAIN_SIZE]) {
            case 2:
                return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, Matrix& rDN_DX){
                    MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2>(rPoints, rX, h, rN, rDN_DX);};
            case 3:
                return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, Matrix& rDN_DX){
                    MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3>(rPoints, rX, h, rN, rDN_DX);};
            default:
                KRATOS_ERROR << "Wrong domain size. MLS shape functions utility cannot be set.";
        }
    }

    ElementSizeFunctionType GetElementSizeFunction(const GeometryType& rGeometry)
    {
        switch (rGeometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return [](const GeometryType& rGeometry)->double{return ElementSizeCalculator<2,3>::AverageElementSize(rGeometry);};
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return [](const GeometryType& rGeometry)->double{return ElementSizeCalculator<3,4>::AverageElementSize(rGeometry);};
            default:
                KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
        }
    }

    void SetSplitElementSupportCloud(
        const Element& rSplitElement,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates)
    {
        // Find the positive side support cloud of nodes
        // Note that we use an unordered_set to ensure that these are unique
        const auto& r_split_geom = rSplitElement.GetGeometry();
        const std::size_t n_nodes = r_split_geom.PointsNumber();
        NodesCloudSetType aux_set;
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            auto p_node = r_split_geom(i_node);
            if (p_node->Is(INTERFACE)) {
                // Add current positive node to map
                aux_set.insert(p_node);

                // Add positive neighbours to map
                auto& r_pos_node_neigh = p_node->GetValue(NEIGHBOUR_NODES);
                const std::size_t n_neigh = r_pos_node_neigh.size();
                for (std::size_t i_neigh = 0; i_neigh < n_neigh; ++i_neigh) {
                    auto& r_neigh = r_pos_node_neigh[i_neigh];
                    if (r_neigh.Is(INTERFACE)) {
                        NodeType::Pointer p_neigh = &r_neigh;
                        aux_set.insert(p_neigh);
                    }
                }
            }
        }

        // Check that the current nodes are enough to perform the MLS calculation (3 in 2D and 4 in 3D)
        // If not sufficient, add the nodal neighbours of the current nodes to the cloud of points
        const std::size_t n_cloud_nodes_temp = aux_set.size();
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        KRATOS_ERROR_IF(n_cloud_nodes_temp == 0) << "Degenerated case with no neighbours. Check if the element " << rSplitElement.Id() << " is intersected and isolated." << std::endl;
        if (n_cloud_nodes_temp < n_dim + 1) {
            NodesCloudSetType aux_extra_set(aux_set);
            for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
                auto& r_it_set_neighs = (*it_set)->GetValue(NEIGHBOUR_NODES);
                const std::size_t n_neigh = r_it_set_neighs.size();
                for (std::size_t i_neigh = 0; i_neigh < n_neigh; ++i_neigh) {
                    auto& r_neigh = r_it_set_neighs[i_neigh];
                    if (r_neigh.Is(ACTIVE)) {
                        NodeType::Pointer p_neigh = &r_neigh;
                        aux_extra_set.insert(p_neigh);
                    }
                }
            }
            aux_set = aux_extra_set;
        }

        // Sort the obtained nodes by id
        const std::size_t n_cloud_nodes = aux_set.size();
        rCloudNodes.resize(n_cloud_nodes);
        std::size_t aux_i = 0;
        for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
            rCloudNodes(aux_i++) = *it_set;
        }
        std::sort(rCloudNodes.ptr_begin(), rCloudNodes.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});

        // Fill the coordinates matrix
        array_1d<double,3> aux_coord;
        rCloudCoordinates.resize(n_cloud_nodes, 3);
        for (std::size_t i_node = 0; i_node < n_cloud_nodes; ++i_node) {
            noalias(aux_coord) = rCloudNodes[i_node].Coordinates();
            rCloudCoordinates(i_node, 0) = aux_coord[0];
            rCloudCoordinates(i_node, 1) = aux_coord[1];
            rCloudCoordinates(i_node, 2) = aux_coord[2];
        }
    }

    double CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin)
    {
        const std::size_t n_nodes = rCloudCoordinates.size1();
        const double squared_rad = IndexPartition<std::size_t>(n_nodes).for_each<MaxReduction<double>>([&](std::size_t I){
            return std::pow(rCloudCoordinates(I,0) - rOrigin(0),2) + std::pow(rCloudCoordinates(I,1) - rOrigin(1),2) + std::pow(rCloudCoordinates(I,2) - rOrigin(2),2);
        });
        return std::sqrt(squared_rad);
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
    // ShiftedBoundaryMeshlessInterfaceProcess& operator=(ShiftedBoundaryMeshlessInterfaceProcess const& rOther);

    /// Copy constructor.
    //ShiftedBoundaryMeshlessInterfaceProcess(ShiftedBoundaryMeshlessInterfaceProcess const& rOther);

    ///@}

}; // Class ShiftedBoundaryMeshlessInterfaceProcess

} // namespace Kratos.

#endif // KRATOS_SHIFTED_BOUNDARY_INTERFACE_UTILITIES  defined
