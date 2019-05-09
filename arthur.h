#ifndef ARTHUR_HEADER
#define ARTHUR_HEADER

#include "_hypre_parcsr_ls.h"


HYPRE_Int HYPRE_IJMatrixRead_binary(const char     *filename,
                             MPI_Comm        comm,
                             HYPRE_Int       type,
                             HYPRE_IJMatrix *matrix);


// HYPRE_Int HYPRE_ParCSRPCGSetupArthur2 ( HYPRE_Solver solver , HYPRE_ParCSRMatrix A , HYPRE_ParCSRMatrix K , HYPRE_ParVector b , HYPRE_ParVector x );
// HYPRE_Int HYPRE_ParCSRPCGSolveArthur2 ( HYPRE_Solver solver , HYPRE_ParCSRMatrix A , HYPRE_ParCSRMatrix K , HYPRE_ParVector b , HYPRE_ParVector x );

// HYPRE_Int HYPRE_ParCSRPCGSetupArthurTRUE ( HYPRE_Solver solver , HYPRE_ParCSRMatrix A , HYPRE_ParCSRMatrix G , HYPRE_ParVector f , HYPRE_ParVector h , HYPRE_ParVector vel, HYPRE_ParVector pres );
// HYPRE_Int HYPRE_ParCSRPCGSolveArthurTRUE ( HYPRE_Solver solver , HYPRE_ParCSRMatrix A , HYPRE_ParCSRMatrix G , HYPRE_ParVector f , HYPRE_ParVector h , HYPRE_ParVector vel, HYPRE_ParVector pres );

int setArrayZeros(double *xVector, int Size);
double dotprod (const double *x,const double *y, int Size);
int copyVector (const double *xVector, double *yVector, int Sized);
void writingSplitVector (double *x, int Sized, int nProc);
int sizeFile(char *str);
int sizeVectorFile (char *str);
int loadVector(char *str, double *xVector, const int nSized);
int load2SplitVector(char *str, char *sstr, double *xVector, const int nSized);
int scalarProduct( double *xVector, int scalar, int Sized);
int incrementOp(double *yVector, const double *xVector,const double scalar,const int Sized);
int D2incrementOp(double *xVector, const double *yVector, const double scalar, const int Sized);



#ifdef __cplusplus
}
#endif

#endif

