#include "ambit-stochastics.h"

int main(int argc, char *argv[]) {
    double tbl[] = {0, 1, 4, 9};
    double a = 0, b = 3;
    int n = 4;
    double y0 = 0;
    for (double x = -1; x < b + 1; x += 0.25)
        printf(" % 4.3f % 4.3f\n", x, linear_interpolation(x, a, b, n, tbl, y0));
    
}
