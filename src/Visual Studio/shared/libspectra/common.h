#pragma once

#define WHICH_LM 0
#define WHICH_LR 1
#define WHICH_LI 2
#define WHICH_LA 3
#define WHICH_SM 4
#define WHICH_SR 5
#define WHICH_SI 6
#define WHICH_SA 7
#define WHICH_BE 8

#define SELECTION_RULE_FROM_INT(i)                      \
switch(i)                                               \
{                                                       \
    case WHICH_LM :                                     \
        selection = SortRule::LargestMagn;              \
        break;                                          \
    case WHICH_LR :                                     \
        selection = SortRule::LargestReal;              \
        break;                                          \
    case WHICH_LI :                                     \
        selection = SortRule::LargestImag;              \
        break;                                          \
    case WHICH_LA :                                     \
        selection = SortRule::LargestAlge;              \
        break;                                          \
    case WHICH_SM :                                     \
        selection = SortRule::SmallestMagn;             \
        break;                                          \
    case WHICH_SR :                                     \
        selection = SortRule::SmallestReal;             \
        break;                                          \
    case WHICH_SI :                                     \
        selection = SortRule::SmallestImag;             \
        break;                                          \
    case WHICH_SA :                                     \
        selection = SortRule::SmallestAlge;             \
        break;                                          \
    case WHICH_BE :                                     \
        selection = SortRule::BothEnds;                 \
        break;                                          \
    default:                                            \
        return -1001;                                   \
}

/* Sparse matrix in column compressed format. */
typedef struct spectra_spmat_t {
	/* Number of rows/columns. */
    int  n;
	/* Array of nonzero values. */
    void *x;
	/* Array of column indices. */
    int  *i;
	/* Array of row pointers. */
    int  *p;
	/* Number of nonzeros in the matrix. */
    int  nnz;
} spectra_spmat;

typedef struct spectra_result_t {
	/* Eigenvalues. */
    void *eigval;
	/* Eigenvectors. */
    void *eigvec;
	/* Number of iterations taken. */
    int  iterations;
	/* Error info. */
    int  info;
} spectra_result;
