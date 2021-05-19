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
        , mrModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))
    {
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
        // Set the modified shape functions factory
        // Note that unique geometry in the mesh is assumed
        const auto& r_begin_geom = mrModelPart.ElementsBegin()->GetGeometry();
        auto p_mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_begin_geom);

        // Get the MLS shape functions function
        const double mls_kernel_rad = 0.2;
        auto p_mls_sh_func = GetMLSShapeFunctionsFunction();

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mrModelPart.Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        //FIXME: This is a temporary solution until we finally decide what to do with the already existent conditions
        // Deactivate the already existent conditions
        block_for_each(mrModelPart.Conditions(),[](Condition& rCondition){rCondition.Set(ACTIVE, false);});

        // Loop the elements to find the intersected ones
        for (auto& rElement : mrModelPart.Elements()) {
            // Check if the element is split
            const auto p_geom = rElement.pGetGeometry();

            // If the element is split, get the interface Gauss pt. data
            // If the element is negative, it is deactivated to not be assembled
            if (IsSplit(*p_geom)) {
                // Set the meshless cloud of point support for the MLS
                Matrix cloud_nodes_coordinates;
                PointerVector<NodeType> cloud_nodes;
                // std::vector<NodeType::Pointer> cloud_nodes;
                SetSplitElementSupportCloud(*p_geom, cloud_nodes, cloud_nodes_coordinates);

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

                    // Store the Gauss pt. weight and normal in the database
                    p_cond->SetValue(NORMAL, pos_int_n[i_g]);
                    p_cond->SetValue(INTEGRATION_WEIGHT, pos_int_w[i_g]);

                    // Calculate the MLS shape functions and gradients and save in the database
                    Vector N_container;
                    Matrix DN_DX_container;
                    p_mls_sh_func(cloud_nodes_coordinates, i_g_coords, mls_kernel_rad, N_container, DN_DX_container);
                    //FIXME: Find variables for these
                    p_cond->SetValue(BDF_COEFFICIENTS, N_container);
                    p_cond->SetValue(LOCAL_AXES_MATRIX, DN_DX_container);
                }

                // Deactivate the split element to avoid assembling it
                // Note that the split elements BC is applied by means of the extension operators
                rElement.Set(ACTIVE, false);
            } else if (IsNegative(*p_geom)) {
                // Deactivate the negative element to avoid assembling it
                rElement.Set(ACTIVE, false);
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

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
        switch (mrModelPart.GetProcessInfo()[DOMAIN_SIZE]) {
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

    void SetSplitElementSupportCloud(
        const GeometryType& rSplitGeometry,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates)
    {
        // Find the positive side support cloud of nodes
        // Note that we use an unordered_set to ensure that these are unique
        const std::size_t n_nodes = rSplitGeometry.PointsNumber();
        NodesCloudSetType aux_set;
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            auto p_node = rSplitGeometry(i_node);
            if (!(p_node->FastGetSolutionStepValue(DISTANCE) < 0.0)) {
                // Add current positive node to map
                aux_set.insert(p_node);

                // Add positive neighbours to map
                auto& r_pos_node_neigh = p_node->GetValue(NEIGHBOUR_NODES);
                const std::size_t n_neigh = r_pos_node_neigh.size();
                for (std::size_t i_neigh = 0; i_neigh < n_neigh; ++i_neigh) {
                    auto& r_neigh = r_pos_node_neigh[i_neigh];
                    if (!(r_neigh.FastGetSolutionStepValue(DISTANCE) < 0.0)) {
                        NodeType::Pointer p_neigh = &r_neigh;
                        aux_set.insert(p_neigh);
                    }
                }
            }
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
