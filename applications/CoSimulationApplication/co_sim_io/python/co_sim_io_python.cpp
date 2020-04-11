#include <pybind11/pybind11.h>

#include "../co_sim_io.h"

PYBIND11_MODULE(CoSimIO, m)
{
    void (*ConnectWithSettingsFileName)(const std::string&, const std::string&) = &CoSimIO::Connect;
    void (*ConnectWithSettings)(const std::string&, CoSimIO::SettingsType) = &CoSimIO::Connect;

    m.def("Connect", ConnectWithSettingsFileName);
    m.def("Connect", ConnectWithSettings);

    m.def("Disconnect", CoSimIO::Disconnect);

    // m.def("__str__", "ssssss");

}