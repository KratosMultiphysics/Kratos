#include "queso/includes/parameters.h"

namespace queso {

std::ostream& operator<< (std::ostream& rOStream, const Parameters& rThis){
    rThis.PrintInfo(rOStream);
    return rOStream;
}

} // End queso namespace