#include "ambit-stochastics.h"

int
main (int argc, char *argv[]) {
    int err = 0;

    struct ambit_dense_array *x = ambit_dense_array_malloc_1d(11);
    struct ambit_dense_array *k = ambit_dense_array_malloc_1d(4);


    for (int i0 = 0; i0 < x->dim[0]; ++i0)
        x->data[ambit_dense_array_linear_index_1d(x, i0)] = sqrt((double)i0);
    for (int i0 = 0; i0 < k->dim[0]; ++i0)
        k->data[ambit_dense_array_linear_index_1d(k, i0)] = (i0 % 2 ? 1.0 * i0 : -1.0 * i0);
    
    err = ambit_circular_convolution_inplace(k, x);
    printf("err = %d\n", err);
    ambit_dense_array_fprintf(stdout, x, true);
    
    return 0;
}
