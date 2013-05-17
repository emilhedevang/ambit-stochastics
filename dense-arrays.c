#include "ambit-stochastics.h"

//
// Dense arrays
//
//
// Naming scheme: ambit_dense_array_SOMETHING_MORE
// Exceptions are helper functions not to be used elsewhere.

//
// Subscripts
//

inline size_t
wrap_subscript(size_t sub, size_t dim) {
    return sub >= dim ? sub - dim : sub;
}

inline void
increment_subscript(size_t *sub, int rank, size_t *dim) {
    for (int i = rank - 1; i >= 0; --i)
        if (++sub[i] < dim[i])
            return;
        else 
            sub[i] = 0;
}

//
// Indexing functions
//

inline size_t
ambit_dense_array_linear_index(struct ambit_dense_array *array, size_t *sub) {
    size_t linear_index = 0;
    for (int i = 0; i < array->rank; ++i)
        linear_index += array->strideembed[i] * wrap_subscript(sub[i], array->dim[i]);
    return linear_index;
}

inline size_t
ambit_dense_array_linear_index_1d(struct ambit_dense_array *array, size_t sub0) {
    return sub0;
}

inline size_t
ambit_dense_array_linear_index_2d(struct ambit_dense_array *array, size_t sub0, size_t sub1) {
    return array->strideembed[0] * sub0 + sub1;
}

inline size_t
ambit_dense_array_linear_index_3d(struct ambit_dense_array *array, 
                                  size_t sub0, size_t sub1, size_t sub2) {
    return array->strideembed[0] * sub0 + array->strideembed[1] * sub1 + sub2;
}

inline R
ambit_dense_array_get(struct ambit_dense_array *array, size_t *sub) {
    return array->data[ambit_dense_array_linear_index(array, sub)];
}

inline R
ambit_dense_array_get_1d(struct ambit_dense_array *array, size_t sub0) {
    return array->data[ambit_dense_array_linear_index_1d(array, sub0)];
}

inline R
ambit_dense_array_get_2d(struct ambit_dense_array *array, size_t sub0, size_t sub1) {
    return array->data[ambit_dense_array_linear_index_2d(array, sub0, sub1)];
}

inline R
ambit_dense_array_get_3d(struct ambit_dense_array *array, size_t sub0, size_t sub1, size_t sub2) {
    return array->data[ambit_dense_array_linear_index_3d(array, sub0, sub1, sub2)];
}

inline void
ambit_dense_array_set(struct ambit_dense_array *array, size_t *sub, R value) {
    array->data[ambit_dense_array_linear_index(array, sub)] = value;
}

//
// Utilities
//

void
ambit_dense_array_parallel_fill(struct ambit_dense_array *a,
                                R (*func)(const void *, const void *),
                                const void  *p,  size_t p_size,
                                const void **pp, size_t pp_size) {
    size_t *sub = NULL;
    void *p0    = NULL;
    void *pp0   = NULL;
    size_t inner_n;
#pragma omp parallel      \
    private(sub, inner_n, p0, pp0)    \
    shared(a) 
    {
        inner_n = a->n / a->dim[0];
        sub     = malloc(a->rank * sizeof(size_t));
        assert(sub);
        
        if (p_size) {
            p0 = malloc(p_size);
            assert(p0);
            memcpy(p0, p, p_size);
        }
        else 
            p0 = NULL;
        
        if (pp_size) {
            pp0 = malloc(pp_size);
            assert(pp0);
            memcpy(pp0, pp[omp_get_thread_num()], pp_size);
        }
        else 
            pp0 = NULL;
        
#pragma omp for 
        for (size_t i0 = 0; i0 < a->dim[0]; ++i0) {
            sub[0] = i0;
            for (int j = 1; j < a->rank; ++j)
                sub[j] = 0;
            for (size_t i = 0; i < inner_n; ++i) {
                ambit_dense_array_set(a, sub, func(p0, pp0));
                increment_subscript(sub, a->rank, a->dim);
            }
        }
        free(sub);
        free(p0);
        free(pp0);
    }    
}


void 
ambit_dense_array_set_all(struct ambit_dense_array *a, R c) {
    size_t *sub0 = NULL;
    size_t lin_idx0;
    size_t inner_n = a->n / a->dim[0];
    R c0;
#pragma omp parallel            \
    private(sub0, lin_idx0, c0) \
    shared(a) 
    {
        sub0 = malloc(a->rank * sizeof(size_t));
        c0   = c;
#pragma omp for 
        for (size_t i0 = 0; i0 < a->dim[0]; ++i0) {
            sub0[0] = i0;
            for (int j = 1; j < a->rank; ++j)
                sub0[j] = 0;
            for (size_t i = 0; i < inner_n; ++i) {
                lin_idx0 = ambit_dense_array_linear_index(a, sub0);
                a->data[lin_idx0] = c0;
                increment_subscript(sub0, a->rank, a->dim);
            }
        }
        free(sub0);
    }    
}

//
// Checks
//

inline bool 
ambit_dense_array_valid (struct ambit_dense_array *a) {
    if (!a || !a->dim || !a->dimembed || !a->strideembed)
        return false;
    for (int i = 0; i < a->rank; ++i)
        if (a->dim[i] == 0 || a->dimembed[i] < a->dim[i])
            return false;
    size_t n = 1;
    size_t nembed = 1;
    for (int i = 0; i < a->rank; ++i) {
        n *= a->dim[i];
        nembed *= a->dimembed[i];
    }
    if (n != a->n || nembed != a->nembed) 
        return false;
    for (int i = a->rank - 1; i >= 0; --i)  
        if (a->strideembed[i] != (i == a->rank - 1 ? 1 : a->dimembed[i + 1] * a->strideembed[i + 1]))
            return false;
    return true;
}

inline bool
ambit_dense_array_ranks_eq (struct ambit_dense_array *a, struct ambit_dense_array *b) {
    return (ambit_dense_array_valid(a) && 
            ambit_dense_array_valid(b) && 
            a->rank == b->rank);
}

inline bool
ambit_dense_array_dims_eq (struct ambit_dense_array *a, struct ambit_dense_array *b) {
    if (!ambit_dense_array_ranks_eq(a, b))
        return false;
    bool ok = true;
    for (int i = 0; i < a->rank; ++i)
        ok = ok && a->dim[i] == b->dim[i];
    return ok;
}

inline bool
ambit_dense_array_dims_leq (struct ambit_dense_array *a, struct ambit_dense_array *b) {
    if (!ambit_dense_array_ranks_eq(a, b))
        return false;
    bool ok = true;
    for (int i = 0; i < a->rank; ++i)
        ok = ok && a->dim[i] <= b->dim[i];
    return ok;
}


inline bool 
ambit_dense_array_embedded_tightly(struct ambit_dense_array *a) {
    if (!ambit_dense_array_valid(a))
        return false;
    for (int i = 0; i < a->rank; ++i) 
        if (a->dim[i] != a->dimembed[i]) 
            return false;
    return true;
}

inline bool 
ambit_dense_array_fftw_r2c_embedded_tightly(struct ambit_dense_array *a) {
    if (!ambit_dense_array_valid(a))
        return false;
    for (int i = 0; i < a->rank; ++i)
        if (a->dimembed[i] != (i < a->rank - 1 ? a->dim[i] : 2 * (a->dim[i] / 2 + 1)))
            return false;
    return true;
}

inline bool 
ambit_dense_array_fftw_r2c_embedded(struct ambit_dense_array *a) {
    if (!ambit_dense_array_valid(a))
        return false;
    for (int i = 0; i < a->rank; ++i)
        if (a->dimembed[i] < (i < a->rank - 1 ? a->dim[i] : 2 * (a->dim[i] / 2 + 1)))
            return false;
    return true;
}

//
// Circular embedding of a into b
//
int
ambit_dense_array_circular_embed (struct ambit_dense_array *a,
                                  struct ambit_dense_array *b,
                                  size_t *offset) {
    if (!ambit_dense_array_dims_leq(a, b)) {
        fprintf(stderr, "Error: First array does not fit into second array.\n");
        return -1;
    }
    size_t  rank = a->rank;
    size_t  n    = a->n;
    size_t *adim = a->dim;
    size_t *asub = malloc(rank * sizeof(size_t));
    size_t *bsub = malloc(rank * sizeof(size_t));
    if (!asub || !bsub) {
        fprintf(stderr, "Error: Could not allocate memory.\n");
        return -1;
    }
    for (int i = 0; i < rank; ++i) 
        asub[i] = 0;
    ptrdiff_t aindex = 0;
    ptrdiff_t bindex = 0;
    for (ptrdiff_t j = 0; j < n; ++j) {
        for (int i = 0; i < rank; ++i)
            bsub[i] = asub[i] + offset[i];
        aindex = ambit_dense_array_linear_index(a, asub);
        bindex = ambit_dense_array_linear_index(b, bsub);
//        printf(" --> %zd %zd\n", aindex, bindex);
        b->data[bindex] = a->data[aindex];
        increment_subscript(asub, rank, adim);
    }
    return 0;
};
    
//
// Memory management
//

struct ambit_dense_array *
ambit_dense_array_malloc (int rank, size_t *dim, size_t *dimembed) {
    
    struct ambit_dense_array *array = NULL;
    R      *data2        = NULL;
    size_t *dim2         = NULL;
    size_t *dimembed2    = NULL;
    size_t *strideembed2 = NULL;
    
    size_t size   = sizeof(R);
    size_t n      = 1;
    size_t nembed = 1;

    if (rank < 0 || !dim || !dimembed) {
        fprintf(stderr, "%s: Invalid input.\n", __func__);
        goto cleanup;
    }
    
    for (int i = 0; i < rank; ++i)
        size *= dimembed[i];
    size = size % AMBIT_ALIGNMENT ? size + AMBIT_ALIGNMENT - (size % AMBIT_ALIGNMENT) : size;

    array = malloc(sizeof(struct ambit_dense_array));
    // R *data = aligned_alloc(AMBIT_ALIGNMENT, size); // Requires GLIBC 2.16
    posix_memalign((void **)&data2, AMBIT_ALIGNMENT, size);       
    dim2         = malloc(rank * sizeof(size_t));
    dimembed2    = malloc(rank * sizeof(size_t));
    strideembed2 = malloc(rank * sizeof(size_t));
    
    if (!array || !data2 || !dim2 || !dimembed2 || !strideembed2) {
        fprintf(stderr, "%s: Could not allocate memory.\n", __func__);
        goto cleanup;
    }

    for (int i = rank - 1; i >= 0; --i)  
        strideembed2[i] = (i == rank - 1 ? 1 : dimembed[i + 1] * strideembed2[i + 1]);

    for (int i = 0; i < rank; ++i) {
        n      *= dim[i];
        nembed *= dimembed[i];
        dim2[i]      = dim[i];
        dimembed2[i] = dimembed[i];
    }

    for (size_t i = 0; i < n; ++i)
        data2[i] = 0.0;
        
    array->rank     = rank;
    array->dim      = dim2;
    array->dimembed = dimembed2;
    array->data     = data2;
    array->strideembed = strideembed2;
    array->n      = n;
    array->nembed = nembed;
    return array;
    
  cleanup:
    free(array);
    free(data2);
    free(dim2);
    free(dimembed2);
    free(strideembed2);
    return NULL;
}

struct ambit_dense_array *
ambit_dense_array_malloc_1d(size_t dim0) {
    size_t *dim      = malloc(1 * sizeof(size_t));
    size_t *dimembed = malloc(1 * sizeof(size_t));
    dim[0]      = dim0;
    dimembed[0] = 2 * (dim0 / 2 + 1);
    return ambit_dense_array_malloc(1, dim, dimembed);
}

struct ambit_dense_array *
ambit_dense_array_malloc_2d(size_t dim0, size_t dim1) {
    size_t *dim      = malloc(2 * sizeof(size_t));
    size_t *dimembed = malloc(2 * sizeof(size_t));
    dim[0]      = dim0;
    dim[1]      = dim1;
    dimembed[0] = dim0;
    dimembed[1] = 2 * (dim1 / 2 + 1);
    return ambit_dense_array_malloc(2, dim, dimembed);
}

struct ambit_dense_array *
ambit_dense_array_malloc_3d(size_t dim0, size_t dim1, size_t dim2) {
    size_t *dim      = malloc(3 * sizeof(size_t));
    size_t *dimembed = malloc(3 * sizeof(size_t));
    dim[0]      = dim0;
    dim[1]      = dim1;
    dim[2]      = dim2;
    dimembed[0] = dim0;
    dimembed[1] = dim1;
    dimembed[2] = 2 * (dim2 / 2 + 1);
    return ambit_dense_array_malloc(3, dim, dimembed);
}

void 
ambit_dense_array_free (struct ambit_dense_array *array) {
    if (array) {
        free(array->dim);
        free(array->dimembed);
        free(array->strideembed);
        free(array->data);
        free(array);
    }
}

//
// IO
//

void ambit_dense_array_fprintf (FILE *stream, struct ambit_dense_array *array, bool also_data) {
    if (ambit_dense_array_valid(array)) {
        fprintf(stream, "rank     = %d\n", array->rank);
        fprintf(stream, "dim      =");
        for (int i = 0; i < array->rank; ++i)
            fprintf(stream, " %zd", array->dim[i]);
        fprintf(stream, "\n");
        fprintf(stream, "dimembed =");
        for (int i = 0; i < array->rank; ++i)
            fprintf(stream, " %zd", array->dimembed[i]);
        fprintf(stream, "\n");
        fprintf(stream, "n        = %zd\n", array->n);
        fprintf(stream, "nembed   = %zd\n", array->nembed);
        if (also_data) {
            fprintf(stream, "data     =");
            size_t n = 1;
            for (int i = 0; i < array->rank; ++i)
                n *= array->dimembed[i];
            for (size_t i = 0; i < n; ++i)
                fprintf(stream, " %g", array->data[i]);
            fprintf(stream, "\n");
        }
    }
    else {
        fprintf(stderr, "%s: Invalid array.\n", __func__);
    }    
}


int 
ambit_dense_array_hdf5_write(struct ambit_dense_array *a, 
                             char *filename, char *datasetname, bool append) {
    herr_t   err       = 0;
    
    hid_t    file      = 0;
    hid_t    dataset   = 0;
    hid_t    filespace = 0;
    hid_t    memspace  = 0;
    
    hsize_t *offset               = NULL;
    unsigned long *dimembed_ulong = NULL;
    
    if (!ambit_dense_array_valid(a)) {
        fprintf(stderr, "%s: Invalid dense array.\n", __func__);
        err = -1;
        goto cleanup;
    }

    if (append) 
        file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    else
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "%s: Could not open or create \"%s\".\n", __func__, filename);
        err = -1;
        goto cleanup;
    }

    if (H5LTfind_dataset(file, datasetname)) {
        fprintf(stderr, "%s: Dataset \"%s\" already exists.\n", __func__, datasetname);
        err = -1;
        goto cleanup;
    }

    memspace = H5Screate_simple(a->rank, (const hsize_t *)a->dimembed, NULL);
    assert(memspace >= 0);
    
    offset = malloc(a->rank * sizeof(hsize_t));
    assert(offset);
    for (int i = 0; i < a->rank; ++i)
        offset[i] = 0;

    err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, (const hsize_t *)a->dim, NULL);
    assert(err >= 0);
    
    filespace = H5Screate_simple(a->rank, (const hsize_t *)a->dim, NULL);
    assert(filespace >= 0);
    
    dataset = H5Dcreate(file, datasetname, H5T_NATIVE_DOUBLE, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset >= 0);
    
    err = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, a->data);
    if (err < 0) {
        fprintf(stderr, "%s: Could not write dataset.\n", __func__);
        err = -1;
        goto cleanup;
    }

    dimembed_ulong = malloc(a->rank * sizeof(unsigned long));
    assert(dimembed_ulong);
    for (int i = 0; i < a->rank; ++i)
        dimembed_ulong[i] = (unsigned long)a->dimembed[i];
    err = H5LTset_attribute_ulong(file, datasetname, "dimembed", dimembed_ulong, a->rank);
    assert(err >= 0);

  cleanup:
    if (memspace  >= 0) H5Sclose(memspace);
    if (filespace >= 0) H5Sclose(filespace);
    if (dataset   >= 0) H5Dclose(dataset);
    if (file      >= 0) H5Fclose(file);
    free(offset);
    free(dimembed_ulong);
    return (int)err;
}



struct ambit_dense_array *
ambit_dense_array_hdf5_read(char *filename, char *datasetname, 
                            size_t *dimembed, bool use_attr_dimembed) {
    herr_t err = 0;

    struct ambit_dense_array *a = NULL;

    int     rank2     = 0;
    size_t *dim2      = NULL;
    size_t *dimembed2 = NULL;
    unsigned long *dimembed_ulong = NULL;
    hsize_t *offset = NULL;
    
    hid_t file      = 0;
    hid_t dataset   = 0;
    hid_t filespace = 0;
    hid_t memspace  = 0;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "%s: Could not open \"%s\".\n", __func__, filename);
        err = -1;
        goto cleanup;
    }
    if (!H5LTfind_dataset(file, datasetname)) {
        fprintf(stderr, "%s: Dataset \"%s\" does not exist.\n", __func__, datasetname);
        err = -1;
        goto cleanup;
    }
    
    dataset = H5Dopen(file, datasetname, H5P_DEFAULT);
    assert(dataset >= 0);
    
    filespace = H5Dget_space(dataset);
    assert(filespace >= 0);
    
    rank2 = H5Sget_simple_extent_ndims(filespace);
    assert(rank2 >= 0);
    
    dim2 = malloc(rank2 * sizeof(size_t));
    dimembed2 = malloc(rank2 * sizeof(size_t));
    assert(dim2 && dimembed2);
    
    err = H5Sget_simple_extent_dims(filespace, (hsize_t *)dim2, NULL);
    assert(err >= 0);
    
    if (dimembed) {
        for (int i = 0; i < rank2; ++i)
            dimembed2[i] = dimembed[i];
    }
    else {
        if (use_attr_dimembed && H5LTfind_attribute(dataset, "dimembed")) {
            dimembed_ulong = malloc(rank2 * sizeof(unsigned long));
            assert(dimembed_ulong);
            err = H5LTget_attribute_ulong(file, datasetname, "dimembed", dimembed_ulong);
            assert(err >= 0);
            for (int i = 0; i < rank2; ++i)
                dimembed2[i] = (size_t)dimembed_ulong[i];
        }
        else {
            for (int i = 0; i < rank2; ++i)
                dimembed2[i] = dim2[i];
        }
    }
    
    a = ambit_dense_array_malloc(rank2, dim2, dimembed2);
    assert(a);

    memspace = H5Screate_simple(a->rank, (const hsize_t *)a->dimembed, NULL);
    assert(memspace >= 0);
    
    offset = malloc(a->rank * sizeof(hsize_t));
    assert(offset);
    for (int i = 0; i < a->rank; ++i)
        offset[i] = 0;
    
    err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, (const hsize_t *)a->dim, NULL);
    assert(err >= 0);
    
    err = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, a->data);
    assert(err >= 0);

    assert(ambit_dense_array_valid(a));
    
  cleanup:
    if (memspace  >= 0) H5Sclose(memspace);
    if (filespace >= 0) H5Sclose(filespace);
    if (dataset   >= 0) H5Dclose(dataset);
    if (file      >= 0) H5Fclose(file);
    free(dim2);
    free(dimembed2);
    free(dimembed_ulong);
    return a;
}
