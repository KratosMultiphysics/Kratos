#include <iostream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <set>
#include <algorithm>
#include <vexcl/devlist.hpp>

using namespace std;

template <class Extensions>
void show_extensions(const Extensions &extensions) {
    cout << "    " << left << setw(32) << "CL_DEVICE_EXTENSIONS" << " = ";
    size_t w = 40;
    for(auto s = extensions.begin(); s != extensions.end(); ++s) {
        w += s->length() + 1;
        if (w > 80) {
            cout << endl << setw(w = 8) << "";
            w += s->length() + 1;
        }
        cout << *s << " ";
    }
}

int main() {

    try {
        auto dev = vex::backend::device_list(vex::Filter::Any);

#if defined(VEXCL_BACKEND_OPENCL)
        cout << "OpenCL devices:" << endl << endl;
        for (auto d = dev.begin(); d != dev.end(); d++) {
            cout << "  " << d->getInfo<CL_DEVICE_NAME>() << endl
                 << "    " << left << setw(32) << "CL_PLATFORM_NAME" << " = "
                 << cl::Platform(d->getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_NAME>()
                 << endl;

#define SHOW_DEVPROP(name) \
            cout << "    " << left << setw(32) << #name << " = " \
                      << d->getInfo< name >() << endl

            SHOW_DEVPROP(CL_DEVICE_TYPE);
            SHOW_DEVPROP(CL_DEVICE_VENDOR);
            SHOW_DEVPROP(CL_DEVICE_VERSION);
            SHOW_DEVPROP(CL_DEVICE_MAX_COMPUTE_UNITS);
            SHOW_DEVPROP(CL_DEVICE_HOST_UNIFIED_MEMORY);
            SHOW_DEVPROP(CL_DEVICE_GLOBAL_MEM_SIZE);
            SHOW_DEVPROP(CL_DEVICE_LOCAL_MEM_SIZE);
            SHOW_DEVPROP(CL_DEVICE_MAX_MEM_ALLOC_SIZE);
            SHOW_DEVPROP(CL_DEVICE_ADDRESS_BITS);
            SHOW_DEVPROP(CL_DEVICE_MAX_CLOCK_FREQUENCY);

#undef SHOW_DEVPROP

            istringstream extensions(d->getInfo<CL_DEVICE_EXTENSIONS>());
            show_extensions( set<string>(
                        istream_iterator<string>(extensions),
                        istream_iterator<string>()
                        ) );

            cout << endl << endl;
        }
#elif defined(VEXCL_BACKEND_COMPUTE)
        cout << "Compute devices:" << endl << endl;
        for(auto d = dev.begin(); d != dev.end(); d++) {
            cout << "  " << d->name() << endl
                 << "    " << left << setw(32) << "CL_PLATFORM_NAME" << " = "
                 << d->platform().name() << endl;

#define SHOW_DEVPROP(name, method) \
            cout << "    " << left << setw(32) << #name << " = " << d->method() << endl

            SHOW_DEVPROP(CL_DEVICE_VENDOR, vendor);
            SHOW_DEVPROP(CL_DEVICE_VERSION, version);
            SHOW_DEVPROP(CL_DEVICE_MAX_COMPUTE_UNITS, compute_units);
            SHOW_DEVPROP(CL_DEVICE_GLOBAL_MEM_SIZE, global_memory_size);
            SHOW_DEVPROP(CL_DEVICE_LOCAL_MEM_SIZE, local_memory_size);
            SHOW_DEVPROP(CL_DEVICE_MAX_MEM_ALLOC_SIZE, max_memory_alloc_size);
            SHOW_DEVPROP(CL_DEVICE_ADDRESS_BITS, address_bits);
            SHOW_DEVPROP(CL_DEVICE_MAX_CLOCK_FREQUENCY, clock_frequency);

#undef SHOW_DEVPROP

            show_extensions( d->extensions() );
            cout << endl << endl;
        }
#elif defined(VEXCL_BACKEND_CUDA)
        cout << "CUDA devices:" << endl << endl;
        unsigned pos = 0;
        for(auto d = dev.begin(); d != dev.end(); d++)
            cout << ++pos << ". " << *d << endl;
#elif defined(VEXCL_BACKEND_JIT)
        cout << "Devices:" << endl << endl;
        unsigned pos = 0;
        for(auto d = dev.begin(); d != dev.end(); d++)
            cout << ++pos << ". " << d->name() << endl;
#else
#error Unsupported backend
#endif
    } catch (const vex::error &e) {
        std::cerr << e << std::endl;
    }
}

// vim: et
