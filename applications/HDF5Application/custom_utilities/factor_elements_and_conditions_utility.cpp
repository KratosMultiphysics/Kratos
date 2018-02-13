// System includes
#include <sstream>

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

template<class TContainerType>
// ElementContainerType || ConditionContainerType
std::vector<std::string> GetLocalComponentNames(const TContainerType& rElements)
{
    KRATOS_TRY;

    if (rElements.size() == 0)
        return std::vector<std::string>();

    std::vector<std::pair<std::string, const typename TContainerType::data_type *>> components(1);
    std::string name;
    // Add first component.
    CompareElementsAndConditionsUtility::GetRegisteredName(rElements.front(), name);
    components[0] = std::pair<std::string, const typename TContainerType::data_type *>(
        name, &rElements.front());

    // Find all remaining components.
    components.reserve(5);
    unsigned pos = 0;
    for (typename TContainerType::const_iterator it = rElements.begin() + 1;
         it != rElements.end(); ++it)
    {
        if (GeometricalObject::IsSame(*it, *components[pos].second) == false)
        {
            bool found = false;
            for (unsigned k = 0; k < components.size(); ++k)
            {
                if (GeometricalObject::IsSame(*it, *components[k].second))
                {
                    found = true;
                    pos = k;
                    break;
                }
            }

            if (!found)
            {
                CompareElementsAndConditionsUtility::GetRegisteredName(*it, name);
                components.push_back(std::pair<std::string, const typename TContainerType::data_type*>(
                    name, &(*it)));
                pos = components.size() - 1;
            }
        }
    }

    // Return component names.
    std::vector<std::string> component_names(components.size());
    for (unsigned i = 0; i < components.size(); ++i)
        component_names[i] = components[i].first;
    return component_names;

    KRATOS_CATCH("");
}

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
        for (int pid = 0; pid != root && pid < comm_size; ++pid)
        {
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
                ierr = MPI_Recv(recv_buf, max_size, MPI_CHAR, root, 0,
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
        component_names = local_component_names;
    }

#else /* KRATOS_USING_MPI */
    component_names = local_component_names;
#endif

    return component_names;

    KRATOS_CATCH("");
}

FactorElementsUtility::FactorElementsUtility(const ElementsContainerType& rElements)
{
    KRATOS_TRY;

    if (rElements.size() == 0)
        return;

    std::vector<std::string> component_names =
        GetComponentNames<ElementsContainerType>(rElements);

    std::vector<const Element*> components(
        component_names.size());
    mFactoredElements.resize(component_names.size());

    for (unsigned i = 0; i < component_names.size(); ++i)
    {
        mFactoredElements[i].first = component_names[i];
        components[i] = &KratosComponents<Element>::Get(component_names[i]);
    }

    int pos = 0;
    for (auto it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
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

            KRATOS_ERROR_IF_NOT(found) << "Did not find element #" << (**it).Id()
                                       << " in components." << std::endl;
        }

        mFactoredElements[pos].second.push_back(*it);
    }

    KRATOS_CATCH("");
}

FactorConditionsUtility::FactorConditionsUtility(const ConditionsContainerType& rConditions)
{
    KRATOS_TRY;

    if (rConditions.size() == 0)
        return;

    std::vector<std::string> component_names =
        GetComponentNames<ConditionsContainerType>(rConditions);

    std::vector<const Condition*> components(
        component_names.size());
    mFactoredConditions.resize(component_names.size());

    for (unsigned i = 0; i < component_names.size(); ++i)
    {
        mFactoredConditions[i].first = component_names[i];
        components[i] = &KratosComponents<Condition>::Get(component_names[i]);
    }

    int pos = 0;
    for (auto it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
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

            KRATOS_ERROR_IF_NOT(found) << "Did not find condition #" << (**it).Id()
                                       << " in components." << std::endl;
        }

        mFactoredConditions[pos].second.push_back(*it);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
