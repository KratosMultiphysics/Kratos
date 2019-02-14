#include <vexcl/vexcl.hpp>

int main() {
    vex::Context ctx(vex::Filter::All);
    std::cout << ctx << std::endl;
}
