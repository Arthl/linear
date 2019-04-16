#ifndef ARTHUR_HEADER
#define ARTHUR_HEADER

#include "_hypre_parcsr_ls.h"

typedef struct
{
   char *       (*CAlloc)        ( size_t count, size_t elt_size );
   HYPRE_Int    (*Free)          ( char *ptr );
   HYPRE_Int    (*CommInfo)      ( void  *A, HYPRE_Int   *my_id,
                                   HYPRE_Int   *num_procs );
   void *       (*CreateVector)  ( void *vector );
   HYPRE_Int    (*DestroyVector) ( void *vector );
   void *       (*MatvecCreate)  ( void *A, void *x );
   HYPRE_Int    (*Matvec)        ( void *matvec_data, HYPRE_Complex alpha, void *A,
                                   void *x, HYPRE_Complex beta, void *y );
   HYPRE_Int    (*MatvecT)       ( void *matvec_data, HYPRE_Complex alpha, void *A,
                                   void *x, HYPRE_Complex beta, void *y );
   HYPRE_Int    (*MatvecDestroy) ( void *matvec_data );
   HYPRE_Real   (*InnerProd)     ( void *x, void *y );
   HYPRE_Int    (*CopyVector)    ( void *x, void *y );
   HYPRE_Int    (*ClearVector)   ( void *x );
   HYPRE_Int    (*ScaleVector)   ( HYPRE_Complex alpha, void *x );
   HYPRE_Int    (*Axpy)          ( HYPRE_Complex alpha, void *x, void *y );

   HYPRE_Int    (*precond)();
   HYPRE_Int    (*precond_setup)();

} hypre_PCGFunctions2;


typedef struct
{
   HYPRE_Real   tol;
   HYPRE_Real   atolf;
   HYPRE_Real   cf_tol;
   HYPRE_Real   a_tol;
   HYPRE_Real   rtol;
   HYPRE_Int    max_iter;
   HYPRE_Int    two_norm;
   HYPRE_Int    rel_change;
   HYPRE_Int    recompute_residual;
   HYPRE_Int    recompute_residual_p;
   HYPRE_Int    stop_crit;
   HYPRE_Int    converged;

   void    *A;
   void    *p;
   void    *s;
   void    *r; /* ...contains the residual.  This is currently kept permanently.
                  If that is ever changed, it still must be kept if logging>1 */

   HYPRE_Int  owns_matvec_data;  /* normally 1; if 0, don't delete it */
   void      *matvec_data;
   void      *precond_data;

   hypre_PCGFunctions2 * functions;

   /* log info (always logged) */
   HYPRE_Int    num_iterations;
   HYPRE_Real   rel_residual_norm;

   HYPRE_Int    print_level; /* printing when print_level>0 */
   HYPRE_Int    logging;  /* extra computations for logging when logging>0 */
   HYPRE_Real  *norms;
   HYPRE_Real  *rel_norms;

} hypre_PCGData2;

#define hypre_PCGDataOwnsMatvecData(pcgdata)  ((pcgdata) -> owns_matvec_data)

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @name generic PCG Solver
 *
 * Description...
 **/
/*@{*/

/**
 * Description...
 *
 * @param param [IN] ...
 **/

hypre_PCGFunctions2 *
hypre_PCGFunctionsCreate2(
   char *       (*CAlloc)        ( size_t count, size_t elt_size ),
   HYPRE_Int    (*Free)          ( char *ptr ),
   HYPRE_Int    (*CommInfo)      ( void  *A, HYPRE_Int   *my_id,
                                   HYPRE_Int   *num_procs ),
   void *       (*CreateVector)  ( void *vector ),
   HYPRE_Int    (*DestroyVector) ( void *vector ),
   void *       (*MatvecCreate)  ( void *A, void *x ),
   HYPRE_Int    (*Matvec)        ( void *matvec_data, HYPRE_Complex alpha, void *A,
                                   void *x, HYPRE_Complex beta, void *y ),
   HYPRE_Int    (*MatvecT)       ( void *matvec_data, HYPRE_Complex alpha, void *A,
                                   void *x, HYPRE_Complex beta, void *y ),
   HYPRE_Int    (*MatvecDestroy) ( void *matvec_data ),
   HYPRE_Real   (*InnerProd)     ( void *x, void *y ),
   HYPRE_Int    (*CopyVector)    ( void *x, void *y ),
   HYPRE_Int    (*ClearVector)   ( void *x ),
   HYPRE_Int    (*ScaleVector)   ( HYPRE_Complex alpha, void *x ),
   HYPRE_Int    (*Axpy)          ( HYPRE_Complex alpha, void *x, void *y ),
   HYPRE_Int    (*PrecondSetup)  ( void *vdata, void *A, void *b, void *x ),
   HYPRE_Int    (*Precond)       ( void *vdata, void *A, void *b, void *x )
   );

/**
 * Description...
 *
 * @param param [IN] ...
 **/

void *
hypre_PCGCreate2( hypre_PCGFunctions2 *pcg_functions2 );

#ifdef __cplusplus
}
#endif


void *hypre_PCGCreate2 ( hypre_PCGFunctions2 *pcg_functions2 );

HYPRE_Int HYPRE_ParCSRPCGSetupArthur2 ( HYPRE_Solver solver , HYPRE_ParCSRMatrix A , HYPRE_ParCSRMatrix K , HYPRE_ParVector b , HYPRE_ParVector x );
HYPRE_Int HYPRE_ParCSRPCGSolveArthur2 ( HYPRE_Solver solver , HYPRE_ParCSRMatrix A , HYPRE_ParCSRMatrix K , HYPRE_ParVector b , HYPRE_ParVector x );

HYPRE_Int HYPRE_ParCSRPCGSetupArthurTRUE ( HYPRE_Solver solver , HYPRE_ParCSRMatrix A , HYPRE_ParCSRMatrix G , HYPRE_ParVector f , HYPRE_ParVector h , HYPRE_ParVector vel, HYPRE_ParVector pres );
HYPRE_Int HYPRE_ParCSRPCGSolveArthurTRUE ( HYPRE_Solver solver , HYPRE_ParCSRMatrix A , HYPRE_ParCSRMatrix G , HYPRE_ParVector f , HYPRE_ParVector h , HYPRE_ParVector vel, HYPRE_ParVector pres );


#ifdef __cplusplus
}
#endif

#endif

