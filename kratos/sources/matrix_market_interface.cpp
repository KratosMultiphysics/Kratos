#include "includes/matrix_market_interface.h"

namespace Kratos {

void TestMatrixRead() {
    const char * harcoded_path = "D:\\a\\Kratos\\Kratos\\kratos\\tests\\test_files\\matrix.mm";
    FILE *f = fopen(harcoded_path, "r");
    fclose(f);
}

}