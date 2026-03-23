#include "xc.h"

int main() {
    xc_func_type func;
    double rho[1] = {1.0};
    double ex[1], vx[1];
    int version_major;

    version_major = (int)XC_VERSION_MAJOR;

    xc_func_init(&func, XC_FUNC_TYPE_LDA, XC_POLARIZATION_NONE, 1, "LDA_X");
    xc_lda(&func, 1, rho, ex, vx);
        
    xc_func_end(&func);
}
