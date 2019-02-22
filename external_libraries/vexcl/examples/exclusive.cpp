#include <vexcl/devlist.hpp>
using namespace vex;

int main() {
    // Get exclusive access to compute devices.
    Context ctx( Filter::Exclusive(Filter::Env) );

    if (ctx.size()) {
        std::cout
            << "Locked devices:" << std::endl
            << ctx << std::endl
            << "Press any key to exit: " << std::endl;
        std::cin.get();
    } else {
        std::cout << "No available devices found." << std::endl;
    }
}

// vim: et
