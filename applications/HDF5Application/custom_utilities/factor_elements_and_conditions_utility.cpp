// System includes
#include <sstream>
#include <utility>

// External includes
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif

// Project includes
#include "utilities/compare_elements_and_conditions_utility.h"

// Application includes
#include "custom_utilities/factor_elements_and_conditions_utility.h"

namespace Kratos
{

namespace
{

// Return the set of element or condition names found in the container.
template<class TContainerType>
// ElementContainerType || ConditionContainerType
std::vector<std::string> GetLocalComponentNames(TContainerType const& rElements)
{
    KRATOS_TRY;

    if (rElements.size() == 0)
        return {};

    std::vector<typename TContainerType::data_type const*> components(1);
    // Add first component.
    components[0] = &rElements.front();

    // Find all remaining components.
    components.reserve(5);
    unsigned pos = 0;
    for (auto it = rElements.begin() + 1; it != rElements.end(); ++it)
    {
        if (GeometricalObject::IsSame(*it, *components[pos]) == false)
        {
            bool found = false;
            for (unsigned k = 0; k < components.size(); ++k)
            {
                if (GeometricalObject::IsSame(*it, *components[k]))
                {
                    found = true;
                    pos = k;
                    break;
                }
            }

            if (!found)
            {
                components.push_back(&(*it));
                pos = components.size() - 1;
            }
        }
    }

    // Return component names.
    std::vector<std::string> component_names(components.size());
    for (unsigned i = 0; i < components.size(); ++i)
        CompareElementsAndConditionsUtility::GetRegisteredName(*components[i], component_names[i]);
    return component_names;

    KRATOS_CATCH("");
}

// Return the set of element or condition names found across all partitions.
template<class TContainerType>
// ElementContainerType || ConditionContainerType
std::vector<std::string> GetComponentNames(const TContainerType& rElements)
{
    KRATOS_TRY;

    std::vector<std::string> local_component_names =
        GetLocalComponentNames<TContainerType>(rElements);

    // Construct global vector of component names across all processes if MPI is used.
    std::vector<std::string> component_names;
#ifdef KRATOS_USING_MPI
    int mpi_is_initialized, ierr;
    ierr = MPI_Initialized(&mpi_is_initialized);
    KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Initialized failed." << std::endl;

    if (mpi_is_initialized)
    {
        constexpr std::size_t max_size = 1000;
        const int root = 0;
        int comm_rank, comm_size, ierr;
        MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

        // Fill the send buffer with local component names.
        std::stringstream ss;
        for (const auto& r_name : local_component_names)
            ss << r_name << '\n';
        const std::string send_buf = ss.str();

        // Gather components from each process and construct the global set of
        // component names.
        for (int pid = 0; pid < comm_size; ++pid)
        {
            if (pid == root)
                continue;

            if (comm_rank == pid)
            {
                const int send_size = send_buf.size() + 1; // Include '\0'.
                ierr = MPI_Send(send_buf.c_str(), send_size, MPI_CHAR, root, 0, MPI_COMM_WORLD);
                KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Send failed." << std::endl;
            }
            else if (comm_rank == root)
            {
                std::stringstream recv_components;
                char recv_buf[max_size];
                MPI_Status status;
                ierr = MPI_Recv(recv_buf, max_size, MPI_CHAR, pid, 0,
                                MPI_COMM_WORLD, &status);
                KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Recv failed." << std::endl;
                recv_components.str(recv_buf);

                std::string name;
                while (std::getline(recv_components, name))
                {
                    if (ss.str().find(name) == std::string::npos)
                        ss << name << '\n';
                }
            }
        }

        // Scatter the global set of component names.
        if (comm_rank == root)
        {
            std::string send_components = ss.str();
            char send_buf[max_size];
            const std::size_t length = send_components.copy(send_buf, max_size);
            KRATOS_ERROR_IF(length == max_size) << "Max string size exceeded." << std::endl;
            send_buf[length] = '\0';
            int ssize = length + 1;
            ierr = MPI_Bcast(&ssize, 1, MPI_INT, root, MPI_COMM_WORLD);
            KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Bcast failed." << std::endl;
            ierr = MPI_Bcast(send_buf, ssize, MPI_CHAR, root, MPI_COMM_WORLD);
            KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Bcast failed." << std::endl;

            std::string name;
            while (std::getline(ss, name))
                component_names.push_back(name);
        }
        else
        {
            std::stringstream recv_components;
            char recv_buf[max_size];
            int ssize = 0;
            ierr = MPI_Bcast(&ssize, 1, MPI_INT, root, MPI_COMM_WORLD);
            KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Bcast failed." << std::endl;
            ierr = MPI_Bcast(recv_buf, ssize, MPI_CHAR, root, MPI_COMM_WORLD);
            KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Bcast failed." << std::endl;
            recv_components.str(recv_buf);

            std::string name;
            while (std::getline(recv_components, name))
                component_names.push_back(name);
        }
    }
    else
    {
        component_names = std::move(local_component_names);
    }

#else /* KRATOS_USING_MPI */
    component_names = std::move(local_component_names);
#endif

    return component_names;

    KRATOS_CATCH("");
}

template <class TContainerType>
void FactorEntities(TContainerType const& rEntities,
                    std::vector<std::string>& rNames,
                    std::vector<TContainerType>& rFactoredEntities)
{
    KRATOS_TRY;

    rNames = GetComponentNames<TContainerType>(rEntities);

    rFactoredEntities.resize(rNames.size());
    if (rEntities.empty())
        return;

    std::vector<typename TContainerType::data_type const*> components(rNames.size());

    for (unsigned i = 0; i < rNames.size(); ++i)
        components[i] = &KratosComponents<typename TContainerType::data_type>::Get(rNames[i]);

    int pos = 0; // If !rEntities.empty(), then components.size() > 0.
    for (auto it = rEntities.ptr_begin(); it != rEntities.ptr_end(); ++it)
    {
        if (GeometricalObject::IsSame(**it, *components[pos]) == false)
        {
            // Find the new position.
            bool found = false;
            for (unsigned k = 0; k < components.size(); ++k)
            {
                if (GeometricalObject::IsSame(**it, *components[k]))
                {
                    pos = k;
                    found = true;
                    break;
                }
            }

            KRATOS_ERROR_IF_NOT(found) << "Did not find entity with id #" << (**it).Id()
                                       << " in components." << std::endl;
        }

        rFactoredEntities[pos].push_back(*it);
    }

    KRATOS_CATCH("");
}

} // unnamed namespace

std::vector<ElementsContainerType> FactorElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    std::vector<std::string> names;
    std::vector<ElementsContainerType> factored_elements;
    FactorElements(rElements, names, factored_elements);
    return factored_elements;

    KRATOS_CATCH("");
}

void FactorElements(ElementsContainerType const& rElements,
                    std::vector<std::string>& rNames,
                    std::vector<ElementsContainerType>& rFactoredElements)
{
    KRATOS_TRY;

    FactorEntities(rElements, rNames, rFactoredElements);

    KRATOS_CATCH("");
}

std::vector<ConditionsContainerType> FactorConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    std::vector<std::string> names;
    std::vector<ConditionsContainerType> factored_conditions;
    FactorConditions(rConditions, names, factored_conditions);
    return factored_conditions;

    KRATOS_CATCH("");
}

void FactorConditions(ConditionsContainerType const& rConditions,
                      std::vector<std::string>& rNames,
                      std::vector<ConditionsContainerType>& rFactoredConditions)
{
    KRATOS_TRY;

    FactorEntities(rConditions, rNames, rFactoredConditions);

    KRATOS_CATCH("");
}

} // namespace Kratos.
