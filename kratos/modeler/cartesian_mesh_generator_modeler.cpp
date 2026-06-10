//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//

// System includes

// External includes

// Project includes
#include "utilities/timer.h"
#include "modeler/cartesian_mesh_generator_modeler.h"

namespace Kratos
{

CartesianMeshGeneratorModeler::CartesianMeshGeneratorModeler()
    : Modeler()
    , mpModel(nullptr)
    , mpSourceModelPart(nullptr)
    , mElementSize(0.0)
{
}

/***********************************************************************************/
/***********************************************************************************/

CartesianMeshGeneratorModeler::CartesianMeshGeneratorModeler(
    Model& rModel,
    Parameters ModelerParameters
    )
    : Modeler(rModel, ModelerParameters)
    , mpModel(&rModel)
    , mpSourceModelPart(nullptr)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    mElementSize = mParameters["element_size"].GetDouble();
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters CartesianMeshGeneratorModeler::GetDefaultParameters() const
{
    return Parameters(R"({
        "input_model_part_name"  : "",
        "output_model_part_name" : "",
        "element_name"           : "Element3D4N",
        "element_size"           : 1.0
    })");
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::SetupModelPart()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mpModel == nullptr)
        << "CartesianMeshGeneratorModeler::SetupModelPart called without a Model." << std::endl;

    const std::string input_name   = mParameters["input_model_part_name"].GetString();
    const std::string output_name  = mParameters["output_model_part_name"].GetString();
    const std::string element_name = mParameters["element_name"].GetString();

    KRATOS_ERROR_IF(input_name.empty())
        << "CartesianMeshGeneratorModeler: \"input_model_part_name\" must not be empty." << std::endl;
    KRATOS_ERROR_IF(output_name.empty())
        << "CartesianMeshGeneratorModeler: \"output_model_part_name\" must not be empty." << std::endl;

    mpSourceModelPart = &mpModel->GetModelPart(input_name);

    ModelPart& r_output = mpModel->HasModelPart(output_name)
                              ? mpModel->GetModelPart(output_name)
                              : mpModel->CreateModelPart(output_name);

    if (!r_output.RecursivelyHasProperties(0))
        r_output.CreateNewProperties(0);

    const Element& r_ref_element = KratosComponents<Element>::Get(element_name);
    GenerateCartesianMesh(r_output, r_ref_element);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::GenerateCartesianMesh(
    ModelPart& rModelPart,
    Element const& rReferenceElement
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mpSourceModelPart == nullptr)
        << "CartesianMeshGeneratorModeler: source model part is not set." << std::endl;

    const unsigned int dimension = rReferenceElement.GetGeometry().WorkingSpaceDimension();

    KRATOS_INFO("CartesianMeshGeneratorModeler") << "Generating mesh for dimension: " << dimension << std::endl;

    Timer::Start("Generating Mesh");

    CalculateBoundingBox(*mpSourceModelPart, mMinPoint, mMaxPoint);
    CalculateDivisionNumbers();

    unsigned int start_node_id = mpSourceModelPart->NumberOfNodes() + 1;
    unsigned int start_element_id = 1;

    if (dimension == 3) {
        // 3D Cartesian tetrahedral mesh using Freudenthal decomposition (6 tets per hex cell)
        Timer::Start("Generating Nodes");

        const double x0 = mMinPoint.X();
        const double y0 = mMinPoint.Y();
        const double z0 = mMinPoint.Z();
        const unsigned int nx = mDivisionsNumber[0];
        const unsigned int ny = mDivisionsNumber[1];
        const unsigned int nz = mDivisionsNumber[2];

        // 3D index with correct strides (nx, nx*ny)
        auto NodeIdx3D = [&](unsigned int i, unsigned int j, unsigned int k) -> unsigned int {
            return i + nx * j + nx * ny * k;
        };

        const unsigned int total_nodes = nx * ny * nz;
        std::vector<Node::Pointer> node_ptrs(total_nodes);
        ModelPart::NodesContainerType::ContainerType& r_nodes_array = rModelPart.NodesArray();
        r_nodes_array.resize(total_nodes);

        auto p_variables_list = rModelPart.pGetNodalSolutionStepVariablesList();
        const SizeType buffer_size = rModelPart.GetBufferSize();
        for (unsigned int k = 0; k < nz; k++) {
            for (unsigned int j = 0; j < ny; j++) {
                for (unsigned int i = 0; i < nx; i++) {
                    auto p_node = Node::Pointer(new Node(
                        start_node_id++,
                        x0 + i * mElementSize,
                        y0 + j * mElementSize,
                        z0 + k * mElementSize));
                    p_node->SetSolutionStepVariablesList(p_variables_list);
                    p_node->SetBufferSize(buffer_size);
                    node_ptrs[NodeIdx3D(i,j,k)] = p_node;
                    r_nodes_array[NodeIdx3D(i,j,k)] = p_node;
                }
            }
        }

        Timer::Stop("Generating Nodes");
        Timer::Start("Generating Elements");

        const unsigned int sx = mSegmentsNumber[0];
        const unsigned int sy = mSegmentsNumber[1];
        const unsigned int sz = mSegmentsNumber[2];
        ModelPart::ElementsContainerType::ContainerType& r_elements_array = rModelPart.ElementsArray();
        r_elements_array.resize(6 * sx * sy * sz);

        Element::NodesArrayType tet_nodes(4);
        unsigned int elem_counter = 0;

        for (unsigned int k = 0; k < sz; k++) {
            for (unsigned int j = 0; j < sy; j++) {
                for (unsigned int i = 0; i < sx; i++) {
                    auto c000 = node_ptrs[NodeIdx3D(i,   j,   k  )];
                    auto c100 = node_ptrs[NodeIdx3D(i+1, j,   k  )];
                    auto c010 = node_ptrs[NodeIdx3D(i,   j+1, k  )];
                    auto c110 = node_ptrs[NodeIdx3D(i+1, j+1, k  )];
                    auto c001 = node_ptrs[NodeIdx3D(i,   j,   k+1)];
                    auto c101 = node_ptrs[NodeIdx3D(i+1, j,   k+1)];
                    auto c011 = node_ptrs[NodeIdx3D(i,   j+1, k+1)];
                    auto c111 = node_ptrs[NodeIdx3D(i+1, j+1, k+1)];

                    // Freudenthal decomposition: 6 tets all sharing the c000-c111 body diagonal
                    tet_nodes(0)=c000; tet_nodes(1)=c100; tet_nodes(2)=c110; tet_nodes(3)=c111;
                    r_elements_array[elem_counter++] = rReferenceElement.Create(start_element_id++, tet_nodes, rReferenceElement.pGetProperties());

                    tet_nodes(0)=c000; tet_nodes(1)=c100; tet_nodes(2)=c101; tet_nodes(3)=c111;
                    r_elements_array[elem_counter++] = rReferenceElement.Create(start_element_id++, tet_nodes, rReferenceElement.pGetProperties());

                    tet_nodes(0)=c000; tet_nodes(1)=c010; tet_nodes(2)=c110; tet_nodes(3)=c111;
                    r_elements_array[elem_counter++] = rReferenceElement.Create(start_element_id++, tet_nodes, rReferenceElement.pGetProperties());

                    tet_nodes(0)=c000; tet_nodes(1)=c010; tet_nodes(2)=c011; tet_nodes(3)=c111;
                    r_elements_array[elem_counter++] = rReferenceElement.Create(start_element_id++, tet_nodes, rReferenceElement.pGetProperties());

                    tet_nodes(0)=c000; tet_nodes(1)=c001; tet_nodes(2)=c101; tet_nodes(3)=c111;
                    r_elements_array[elem_counter++] = rReferenceElement.Create(start_element_id++, tet_nodes, rReferenceElement.pGetProperties());

                    tet_nodes(0)=c000; tet_nodes(1)=c001; tet_nodes(2)=c011; tet_nodes(3)=c111;
                    r_elements_array[elem_counter++] = rReferenceElement.Create(start_element_id++, tet_nodes, rReferenceElement.pGetProperties());
                }
            }
        }

        Timer::Stop("Generating Elements");
        Timer::Stop("Generating Mesh");
        return;
    }

    KRATOS_ERROR_IF(dimension != 2) << "Only 2D and 3D meshes are supported!" << std::endl;

    // 2D path: quadrilateral elements with inside-test via ray casting
    const unsigned int number_of_nodes = mDivisionsNumber[0] * mDivisionsNumber[1] * mDivisionsNumber[2];
    const unsigned int number_of_elements = mSegmentsNumber[0] * mSegmentsNumber[1] * mSegmentsNumber[2];

    KRATOS_INFO("CartesianMeshGeneratorModeler") << "Number of nodes: " << number_of_nodes << std::endl;
    KRATOS_INFO("CartesianMeshGeneratorModeler") << "Number of elements: " << number_of_elements << std::endl;

    CalculateBoundaryIntersections(*mpSourceModelPart);
    CalculateIsInside(*mpSourceModelPart);
    CalculateNormals();

    Timer::Start("Generating Nodes");

    const double x0 = mMinPoint.X();
    const double y0 = mMinPoint.Y();
    const double z0 = mMinPoint.Z();

    ModelPart::NodesContainerType::ContainerType& r_nodes_array = rModelPart.NodesArray();
    ModelPart::NodesContainerType::ContainerType temp_nodes_array(number_of_nodes);

    ModelPart::ElementsContainerType::ContainerType& r_elements_array = rModelPart.ElementsArray();

    auto p_variables_list = rModelPart.pGetNodalSolutionStepVariablesList();
    const SizeType buffer_size = rModelPart.GetBufferSize();
    for (unsigned int k = 0; k < mDivisionsNumber[2]; k++)
        for (unsigned int j = 0; j < mDivisionsNumber[1]; j++)
            for (unsigned int i = 0; i < mDivisionsNumber[0]; i++) {
                auto p_node = Node::Pointer(new Node(
                    start_node_id++,
                    x0 + i * mElementSize,
                    y0 + j * mElementSize,
                    z0 + k * mElementSize));
                p_node->SetSolutionStepVariablesList(p_variables_list);
                p_node->SetBufferSize(buffer_size);
                temp_nodes_array[NodeIndex(i,j,k)] = p_node;
            }

    unsigned int number_of_active_nodes = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
        if (mIsInside[i])
            number_of_active_nodes++;

    r_nodes_array.resize(number_of_active_nodes);

    unsigned int index = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
        if (mIsInside[i])
            r_nodes_array[index++] = temp_nodes_array[i];

    Timer::Stop("Generating Nodes");
    Timer::Start("Generating Elements");

    Element::NodesArrayType element_nodes(4);

    unsigned int number_of_active_elements = 0;
    for (unsigned int j = 0; j < mSegmentsNumber[1]; j++)
        for (unsigned int i = 0; i < mSegmentsNumber[0]; i++)
            if (mIsInside[NodeIndex(i,j,0)] & mIsInside[NodeIndex(i+1,j,0)] &
                mIsInside[NodeIndex(i+1,j+1,0)] & mIsInside[NodeIndex(i,j+1,0)])
                number_of_active_elements++;

    r_elements_array.resize(number_of_active_elements);

    unsigned int counter = 0;
    for (unsigned int j = 0; j < mSegmentsNumber[1]; j++) {
        for (unsigned int i = 0; i < mSegmentsNumber[0]; i++) {
            if (mIsInside[NodeIndex(i,j,0)] & mIsInside[NodeIndex(i+1,j,0)] &
                mIsInside[NodeIndex(i+1,j+1,0)] & mIsInside[NodeIndex(i,j+1,0)]) {
                element_nodes(0) = temp_nodes_array[NodeIndex(i,   j,   0)];
                element_nodes(1) = temp_nodes_array[NodeIndex(i+1, j,   0)];
                element_nodes(2) = temp_nodes_array[NodeIndex(i+1, j+1, 0)];
                element_nodes(3) = temp_nodes_array[NodeIndex(i,   j+1, 0)];
                r_elements_array[counter++] = rReferenceElement.Create(start_element_id++, element_nodes, rReferenceElement.pGetProperties());
            }
        }
    }

    for (auto& r_elem : mpSourceModelPart->Elements()) {
        auto& r_geometry = r_elem.GetGeometry();
        array_1d<double,3> coords1, coords2;
        noalias(coords1) = r_geometry[0].Coordinates() + mNormals[r_geometry[0].Id()-1];
        noalias(coords2) = r_geometry[1].Coordinates() + mNormals[r_geometry[1].Id()-1];
        Point point1(coords1[0], coords1[1], coords1[2]);
        Point point2(coords2[0], coords2[1], coords2[2]);
        unsigned int index1 = FindNearestNodeIndex(point1, mNormals[r_geometry[0].Id()-1]);
        unsigned int index2 = FindNearestNodeIndex(point2, mNormals[r_geometry[1].Id()-1]);
        element_nodes(0) = r_geometry(0);
        element_nodes(1) = r_geometry(1);
        element_nodes(2) = temp_nodes_array[index2];
        element_nodes(3) = temp_nodes_array[index1];
        r_elements_array.push_back(rReferenceElement.Create(start_element_id++, element_nodes, rReferenceElement.pGetProperties()));
    }

    Timer::Stop("Generating Elements");
    Timer::Stop("Generating Mesh");

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::CalculateNormals()
{
    const array_1d<double,3> zero = ZeroVector(3);

    if (mNormals.size() != mpSourceModelPart->NumberOfNodes())
        mNormals.resize(mpSourceModelPart->NumberOfNodes(), zero);
    else
        std::fill(mNormals.begin(), mNormals.end(), zero);

    const double coefficient = mElementSize / 2.0;

    for (auto& r_elem : mpSourceModelPart->Elements()) {
        auto& r_geometry = r_elem.GetGeometry();
        array_1d<double,3> normal;
        normal[0] =   r_geometry[1].Y() - r_geometry[0].Y();
        normal[1] = -(r_geometry[1].X() - r_geometry[0].X());
        normal[2] = 0.0;
        normal *= coefficient / r_geometry.Length();
        mNormals[r_geometry[0].Id()-1] += normal;
        mNormals[r_geometry[1].Id()-1] += normal;
    }

    for (auto& r_node : mpSourceModelPart->Nodes())
        noalias(r_node.FastGetSolutionStepValue(NORMAL)) = mNormals[r_node.Id()-1];
}

/***********************************************************************************/
/***********************************************************************************/

unsigned int CartesianMeshGeneratorModeler::FindNearestNodeIndex(
    Point& rThisPoint,
    array_1d<double,3>& rNormal
    )
{
    const double x = (rThisPoint.X() - mMinPoint.X()) / mElementSize;
    const double y = (rThisPoint.Y() - mMinPoint.Y()) / mElementSize;

    unsigned int i = static_cast<unsigned int>(x);
    unsigned int j = static_cast<unsigned int>(y);

    if (mIsInside[NodeIndex(i,j,0)])
        return NodeIndex(i,j,0);

    if (rNormal[0] >= 0.0) {
        if (rNormal[1] >= 0.0) {
            if (rNormal[0] > rNormal[1]) i++;
            else                          j++;
        } else {
            if (rNormal[0] > -rNormal[1]) i++;
            else                           j--;
        }
    } else {
        if (rNormal[1] >= 0.0) {
            if (-rNormal[0] > rNormal[1]) i--;
            else                           j++;
        } else {
            if (-rNormal[0] > -rNormal[1]) i--;
            else                            j--;
        }
    }

    return NodeIndex(i,j,0);
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::CalculateIsInside(ModelPart& rModelPart)
{
    const unsigned int number_of_nodes = mDivisionsNumber[0] * mDivisionsNumber[1] * mDivisionsNumber[2];

    if (mIsInside.size() != number_of_nodes)
        mIsInside.resize(number_of_nodes, 0);
    else
        std::fill(mIsInside.begin(), mIsInside.end(), 0);

    const int size = static_cast<int>(mSegmentsNumber[1]) + 1;

    for (int j = 0; j < size; j++) {
        std::vector<double>& r_j_intersections = mIntersections[j];

        for (auto j_x = r_j_intersections.begin(); j_x != r_j_intersections.end(); j_x++) {
            const auto start = j_x++;
            const unsigned int i_start = static_cast<unsigned int>(*start / mElementSize);
            const unsigned int i_end   = static_cast<unsigned int>(*j_x   / mElementSize);
            for (unsigned int i = i_start + 1; i < i_end; i++)
                mIsInside[NodeIndex(i, j, 0)] = true;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::CalculateBoundaryIntersections(ModelPart& rModelPart)
{
    if (mIntersections.size() != mSegmentsNumber[1] + 1)
        mIntersections.resize(mSegmentsNumber[1] + 1);

    for (auto& r_elem : rModelPart.Elements()) {
        const auto& r_geometry = r_elem.GetGeometry();

        double x1, x2, y1, y2;
        if (r_geometry[0].Y() < r_geometry[1].Y()) {
            x1 = r_geometry[0].X() - mMinPoint.X();
            x2 = r_geometry[1].X() - mMinPoint.X();
            y1 = r_geometry[0].Y() - mMinPoint.Y();
            y2 = r_geometry[1].Y() - mMinPoint.Y();
        } else {
            x1 = r_geometry[1].X() - mMinPoint.X();
            x2 = r_geometry[0].X() - mMinPoint.X();
            y1 = r_geometry[1].Y() - mMinPoint.Y();
            y2 = r_geometry[0].Y() - mMinPoint.Y();
        }

        unsigned int i_start = static_cast<unsigned int>(y1 / mElementSize);
        const unsigned int i_end = static_cast<unsigned int>(y2 / mElementSize);

        if (i_start * mElementSize < y1)
            i_start++;

        const double m       = (y1 != y2) ? (x2 - x1) / (y2 - y1) : 0.0;
        const double delta_x = m * mElementSize;
        double x = x1 + (i_start * mElementSize - y1) * m;

        for (unsigned int i = i_start; i <= i_end; i++) {
            mIntersections[i].push_back(x);
            x += delta_x;
        }
    }

    for (auto& r_row : mIntersections)
        std::sort(r_row.begin(), r_row.end());
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::CalculateBoundingBox(
    ModelPart& rModelPart,
    Point& rMinPoint,
    Point& rMaxPoint
    )
{
    if (rModelPart.NumberOfElements() == 0 && rModelPart.NumberOfNodes() == 0) {
        rMinPoint = Point();
        rMaxPoint = Point();
        return;
    }

    // Fall back to node-based BB when no elements are present (e.g. condition-only model parts from STL)
    if (rModelPart.NumberOfElements() == 0) {
        rMinPoint = *rModelPart.NodesBegin();
        rMaxPoint = rMinPoint;
        for (const auto& r_node : rModelPart.Nodes()) {
            for (unsigned int i = 0; i < 3; i++) {
                if (rMinPoint[i] > r_node[i]) rMinPoint[i] = r_node[i];
                if (rMaxPoint[i] < r_node[i]) rMaxPoint[i] = r_node[i];
            }
        }
        return;
    }

    if (rModelPart.ElementsBegin()->GetGeometry().empty()) {
        rMinPoint = Point();
        rMaxPoint = Point();
        return;
    }

    rMinPoint = rModelPart.ElementsBegin()->GetGeometry()[0];
    rMaxPoint = rMinPoint;

    for (const auto& r_elem : rModelPart.Elements()) {
        for (const auto& r_point : r_elem.GetGeometry()) {
            for (unsigned int i = 0; i < Point::Dimension(); i++) {
                if (rMinPoint[i] > r_point[i]) rMinPoint[i] = r_point[i];
                if (rMaxPoint[i] < r_point[i]) rMaxPoint[i] = r_point[i];
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::CalculateDivisionNumbers()
{
    if (mElementSize == 0.0)
        return;

    for (unsigned int i = 0; i < Point::Dimension(); i++) {
        const double delta = mMaxPoint[i] - mMinPoint[i];
        int segments_number = static_cast<int>(delta / mElementSize);

        if ((segments_number * mElementSize) < delta)
            segments_number++;

        mSegmentsNumber[i] = segments_number;
        mDivisionsNumber[i] = segments_number + 1;

        if (mSegmentsNumber[i] == 0)
            mSegmentsNumber[i]++;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::string CartesianMeshGeneratorModeler::Info() const
{
    return "CartesianMeshGeneratorModeler";
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

/***********************************************************************************/
/***********************************************************************************/

void CartesianMeshGeneratorModeler::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos