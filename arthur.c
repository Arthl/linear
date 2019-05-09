#include "arthur.h"

#include <stdio.h>

void arthur() {
}


HYPRE_Int 
HYPRE_IJMatrixRead_binary( const char     *filename,
                           MPI_Comm        comm,
                           HYPRE_Int       type,
                           HYPRE_IJMatrix *matrix_ptr )
{
   HYPRE_IJMatrix  matrix;
   HYPRE_Int       ilower, iupper, jlower, jupper;
   HYPRE_Int       nnz, iVal;
   HYPRE_Int       ncols, I, J;
   HYPRE_Real      value, realI, realJ; /* no complex number allowed */
   HYPRE_Int       myid, ret;
   char            new_filename[255];
   FILE           *file;

   hypre_MPI_Comm_rank(comm, &myid);
   
   hypre_sprintf(new_filename,"%s.%05d", filename, myid);

   if ((file = fopen(new_filename, "rb")) == NULL)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }

   fread(&ilower, sizeof(ilower), 1, file);
   fread(&iupper, sizeof(iupper), 1, file);
   fread(&jlower, sizeof(jlower), 1, file);
   fread(&jupper, sizeof(jupper), 1, file);
   fread(&nnz,    sizeof(nnz),    1, file); /* added by sz */

   HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &matrix);

   HYPRE_IJMatrixSetObjectType(matrix, type);
   HYPRE_IJMatrixInitialize(matrix);

   ncols = 1;
   while( !feof(file)) 
   {
      ret = fread(&realI,     sizeof(realI),     1, file);
      ret = fread(&realJ,     sizeof(realJ),     1, file);
      ret = fread(&value, sizeof(value), 1, file);

      I = (HYPRE_Int) realI;
      J = (HYPRE_Int) realJ;

      if (I < ilower || I > iupper)
         HYPRE_IJMatrixAddToValues(matrix, 1, &ncols, &I, &J, &value);
      else
         HYPRE_IJMatrixSetValues(matrix, 1, &ncols, &I, &J, &value);
   }

   HYPRE_IJMatrixAssemble(matrix);

   fclose(file);

   *matrix_ptr = matrix;

   return hypre_error_flag;
}


/* Functions Implemented for Matrix-Vector Operations */
int setArrayZeros(double *xVector, int Size)
{
	int i;
	i = 0;
    for (i = 0; i < Size; i++)
    {
        xVector[i]= 0;
    }
    return ( 1 );
}
/* dot product (x,y) */
double dotprod (const double *x,const double *y, int Size)
{
	int i ;
    double prod = 0;
    for (i = 0 ; i < Size ; i++)
    {
    	prod += x[i]*y[i];
	}
    return ( prod ) ;
}

int copyVector (const double *xVector, double *yVector, int Sized)
{
	int i ;
	for (i = 0 ; i < Sized ; i++)
	{
		yVector[i] = xVector[i];
	}
	return ( 1 ) ;
}

void writingSplitVector (double *xVector, int Sized, int nProc)
{
  int j;
  int partSize;
  partSize = Sized/nProc;
  for(j=0; j < nProc; j++)
  {
    char buf[0x100];
    int len;
    len = sizeof(buf);
    snprintf(buf, len, "myfile.0.%1d", j);
    FILE * fp;
    fp = fopen (buf,"w");
    fprintf (fp, "%d\n", partSize);
    /* write lines of text into the file stream*/
    int z;
    for(z = 0; z < partSize ;z++)
    {
      fprintf (fp, "%.15e\n", xVector[z + partSize*j]);
    }
    fclose (fp);
  }
}

// What about instead ?
//    if (myid==0) 
  //  {
  //     ACTION
  //  }
  
int sizeFile(char *str)
{
	FILE *in_file  = fopen(str, "r"); // read only 
  	if (in_file == NULL) 
  	{   
      printf("Error! Could not open file\n"); 
      exit(-1); // must include stdlib.h 
    } 
	int Sized;
	fscanf(in_file, "%d", &Sized);
	fclose(in_file);
	return (Sized) ;
  }  
  
int sizeVectorFile (char *str)
{
	FILE *in_file  = fopen(str, "r"); // read only 
  	if (in_file == NULL) 
  	{   
      printf("Error! Could not open file\n"); 
      exit(-1); // must include stdlib.h 
    } 
    int Sized;
    fscanf(in_file, "%d", &Sized);
    return ( Sized );
  }  
  
int loadVector(char *str, double *xVector, const int nSized)
{
	FILE *in_file  = fopen(str, "r"); // read only 
  	if (in_file == NULL) 
  	{   
      printf("Error! Could not open file\n"); 
      exit(-1); // must include stdlib.h 
    } 
    int i;
    int Sized;
    fscanf(in_file, "%d", &Sized);
    if (Sized != nSized)
    {
    	printf("Error! Wrong size of the Vector\n"); 
	}
	for (i = 0; i < Sized; i++)
    {
		fscanf(in_file, "%lf", &xVector[i]);
    }
	fclose(in_file);
	return ( 1 ) ;
}


int load2SplitVector(char *str, char *sstr, double *xVector, const int nSized)
{
	FILE *in_file  = fopen(str, "r"); // read only 
  	if (in_file == NULL) 
  	{   
      printf("Error! Could not open file\n"); 
      exit(-1); // must include stdlib.h 
    } 
    int i;
    int Sized;
    fscanf(in_file, "%d", &Sized);
    if (Sized != nSized/2)
    {
    	printf("Error! Wrong size of the Vector\n"); 
	}
	for (i = 0; i < Sized; i++)
    {
		fscanf(in_file, "%lf", &xVector[i]);
    }
	fclose(in_file);

	FILE *in2_file  = fopen(sstr, "r"); // read only 
  	if (in2_file == NULL) 
  	{   
      printf("Error! Could not open file\n"); 
      exit(-1); // must include stdlib.h 
    } 
    fscanf(in2_file, "%d", &Sized);
    if (Sized != nSized/2)
    {
    	printf("Error! Wrong size of the Vector\n"); 
	}
	for (i = 0; i < Sized; i++)
    {
		fscanf(in2_file, "%lf", &xVector[i+Sized]);
    }
	fclose(in2_file);
	return ( 1 ) ;
}


int scalarProduct( double *xVector, int scalar, int Sized)
{
	int i;
	for(i = 0; i < Sized; i++)
  {
    xVector[i] *= scalar;
  }
	return ( 1 );
}

int incrementOp(double *yVector, const double *xVector,const double scalar,const int Sized)
{
	int i;
	for(i = 0; i < Sized; i++)
  {
    yVector[i] += scalar*xVector[i];
  }
	return ( 1 );
}

int D2incrementOp(double *xVector, const double *yVector, const double scalar, const int Sized)
{
	int i;
	for(i = 0; i < Sized; i++)
  {
    xVector[i] = yVector[i] + scalar*xVector[i];
  }
	return ( 1 );
}







//LINES should be into krylov.h, lines 954-1070, computes matvec with transpose
//   HYPRE_Int    (*MatvecT)       ( void *matvec_data, HYPRE_Complex alpha, void *A,
//                                   void *x, HYPRE_Complex beta, void *y );
//TO COMPUTE (*(pcg_functions->MatvecT))(matvec_data, -1.0, A, x, 1.0, r);


// HYPRE_Int
// hypre_PCGSetupArthurTRUE( void *pcg_vdata,
//                 void *K,
//                 void *G,
//                 void *f,
//                 void *h,
//                 void *vel,
//                 void *pres         )
// {
//    hypre_PCGData *pcg_data =  (hypre_PCGData *)pcg_vdata;
//    hypre_PCGFunctions *pcg_functions = pcg_data->functions;
//    HYPRE_Int            max_iter         = (pcg_data -> max_iter);
//    HYPRE_Int          (*precond_setup)(void*,void*,void*,void*) = (pcg_functions -> precond_setup);
//    void          *precond_data     = (pcg_data -> precond_data);


//    (pcg_data -> A) = K; //TROUBLE HERE //NOT NECESSARY ????

//    /*--------------------------------------------------
//     * The arguments for CreateVector are important to
//     * maintain consistency between the setup and
//     * compute phases of matvec and the preconditioner.
//     *--------------------------------------------------*/

//    if ( pcg_data -> p != NULL )
//       (*(pcg_functions->DestroyVector))(pcg_data -> p);
//    (pcg_data -> p) = (*(pcg_functions->CreateVector))(vel); //not correct vel

//    if ( pcg_data -> s != NULL )
//       (*(pcg_functions->DestroyVector))(pcg_data -> s);
//    (pcg_data -> s) = (*(pcg_functions->CreateVector))(vel); //not correct vel

//    if ( pcg_data -> r != NULL )
//       (*(pcg_functions->DestroyVector))(pcg_data -> r);
//    (pcg_data -> r) = (*(pcg_functions->CreateVector))(f); //not correct vel

//    if ( pcg_data -> matvec_data != NULL && pcg_data->owns_matvec_data )
//       (*(pcg_functions->MatvecDestroy))(pcg_data -> matvec_data);
//    (pcg_data -> matvec_data) = (*(pcg_functions->MatvecCreate))(K, vel);

//    precond_setup(precond_data, K, f, vel); //DONT KNOW WHAT TO DO HERE

//    /*-----------------------------------------------------
//     * Allocate space for log info
//     *-----------------------------------------------------*/

//    if ( (pcg_data->logging)>0  || (pcg_data->print_level)>0 ) 
//    {
//       if ( (pcg_data -> norms) != NULL )
//          hypre_TFreeF( pcg_data -> norms, pcg_functions );
//       (pcg_data -> norms)     = hypre_CTAllocF( HYPRE_Real, max_iter + 1,
//                                                 pcg_functions);

//       if ( (pcg_data -> rel_norms) != NULL )
//          hypre_TFreeF( pcg_data -> rel_norms, pcg_functions );
//       (pcg_data -> rel_norms) = hypre_CTAllocF( HYPRE_Real, max_iter + 1,
//                                                 pcg_functions );
//    }

//    return hypre_error_flag;
// }


// HYPRE_Int 
// HYPRE_PCGSetupArthurTRUE( HYPRE_Solver solver,
//                 HYPRE_Matrix K,
//                 HYPRE_Matrix G,
//                 HYPRE_Vector f,
//                 HYPRE_Vector h,
//                 HYPRE_Vector vel,
//                 HYPRE_Vector pres )
// {
//    return( hypre_PCGSetupArthurTRUE( solver,
//                            K,
//                            G,
//                            f,
//                            h,
//                            vel,
//                            pres ) );
// }

// HYPRE_Int 
// HYPRE_ParCSRPCGSetupArthurTRUE( HYPRE_Solver solver,
//                       HYPRE_ParCSRMatrix K,
//                       HYPRE_ParCSRMatrix G,
//                       HYPRE_ParVector f,
//                       HYPRE_ParVector h,
//                       HYPRE_ParVector vel,
//                       HYPRE_ParVector pres      )
// {
//    return( HYPRE_PCGSetupArthurTRUE( solver,
//                            (HYPRE_Matrix) K,
//                            (HYPRE_Matrix) G,
//                            (HYPRE_Vector) f,
//                            (HYPRE_Vector) h,
//                            (HYPRE_Vector) vel,
//                            (HYPRE_Vector) pres ) );
// }


// /*--------------------------------------------------------------------------
//  * hypre_PCGSolveArthurTRUE
//  *--------------------------------------------------------------------------
//  *
//  * We use the following convergence test as the default (see Ashby, Holst,
//  * Manteuffel, and Saylor):
//  *
//  *       ||e||_A                           ||r||_C
//  *       -------  <=  [kappa_A(C*A)]^(1/2) -------  < tol
//  *       ||x||_A                           ||b||_C
//  *
//  * where we let (for the time being) kappa_A(CA) = 1.
//  * We implement the test as:
//  *
//  *       gamma = <C*r,r>/<C*b,b>  <  (tol^2) = eps
//  *
//  *--------------------------------------------------------------------------*/

// HYPRE_Int
// hypre_PCGSolveArthurTRUE( void *pcg_vdata,
//                void *A,
//                void *G,
//                void *x,
//                void *b,
//                void *vel,
//                void *pres         )
// {
//    hypre_PCGData  *pcg_data     =  (hypre_PCGData *)pcg_vdata;
//    hypre_PCGFunctions *pcg_functions = pcg_data->functions;
//    // hypre_PCGFunctions2 *pcg_functions2 = pcg_data->functions;
//    // And line 491

//    HYPRE_Real      r_tol        = (pcg_data -> tol);
//    HYPRE_Real      a_tol        = (pcg_data -> a_tol);
//    HYPRE_Real      atolf        = (pcg_data -> atolf);
//    HYPRE_Real      cf_tol       = (pcg_data -> cf_tol);
//    HYPRE_Real      rtol         = (pcg_data -> rtol);
//    HYPRE_Int             max_iter     = (pcg_data -> max_iter);
//    HYPRE_Int             two_norm     = (pcg_data -> two_norm);
//    HYPRE_Int             rel_change   = (pcg_data -> rel_change);
//    HYPRE_Int             recompute_residual = (pcg_data -> recompute_residual);
//    HYPRE_Int             recompute_residual_p = (pcg_data -> recompute_residual_p);
//    HYPRE_Int             stop_crit    = (pcg_data -> stop_crit);
// /*
//    HYPRE_Int             converged    = (pcg_data -> converged);
// */
//    void           *p            = (pcg_data -> p);
//    void           *s            = (pcg_data -> s);
//    void           *r            = (pcg_data -> r);
//    void           *matvec_data  = (pcg_data -> matvec_data);
//    HYPRE_Int           (*precond)(void*,void*,void*,void*)   = (pcg_functions -> precond);
//    void           *precond_data = (pcg_data -> precond_data);
//    HYPRE_Int             print_level  = (pcg_data -> print_level);
//    HYPRE_Int             logging      = (pcg_data -> logging);
//    HYPRE_Real     *norms        = (pcg_data -> norms);
//    HYPRE_Real     *rel_norms    = (pcg_data -> rel_norms);
                
//    HYPRE_Real      alpha, beta;
//    HYPRE_Real      gamma, gamma_old;
//    HYPRE_Real      bi_prod, eps;
//    HYPRE_Real      pi_prod, xi_prod;
//    HYPRE_Real      ieee_check = 0.;
                
//    HYPRE_Real      i_prod = 0.0;
//    HYPRE_Real      i_prod_0 = 0.0;
//    HYPRE_Real      cf_ave_0 = 0.0;
//    HYPRE_Real      cf_ave_1 = 0.0;
//    HYPRE_Real      weight;
//    HYPRE_Real      ratio;

//    HYPRE_Real      guard_zero_residual, sdotp;
//    HYPRE_Int             tentatively_converged = 0;
//    HYPRE_Int             recompute_true_residual = 0;

//    HYPRE_Int             i = 0;
//    HYPRE_Int             my_id, num_procs;

//    (pcg_data -> converged) = 0;

//    (*(pcg_functions->CommInfo))(A,&my_id,&num_procs);

//    /*-----------------------------------------------------------------------
//     * With relative change convergence test on, it is possible to attempt
//     * another iteration with a zero residual. This causes the parameter
//     * alpha to go NaN. The guard_zero_residual parameter is to circumvent
//     * this. Perhaps it should be set to something non-zero (but small).
//     *-----------------------------------------------------------------------*/

//    guard_zero_residual = 0.0;

//    /*-----------------------------------------------------------------------
//     * Start pcg solve
//     *-----------------------------------------------------------------------*/

//    /* compute eps */
//    if (two_norm)
//    {
//       /* bi_prod = <b,b> */
//       bi_prod = (*(pcg_functions->InnerProd))(b, b);
//       if (print_level > 1 && my_id == 0) 
//           hypre_printf("<b,b>: %e\n",bi_prod);
//    }
//    else
//    {
//       /* bi_prod = <C*b,b> */
//       (*(pcg_functions->ClearVector))(p);
//       precond(precond_data, A, b, p);
//       bi_prod = (*(pcg_functions->InnerProd))(p, b);
//       if (print_level > 1 && my_id == 0)
//           hypre_printf("<C*b,b>: %e\n",bi_prod);
//    };

//    /* Since it is does not diminish performance, attempt to return an error flag
//       and notify users when they supply bad input. */
//    if (bi_prod != 0.) ieee_check = bi_prod/bi_prod; /* INF -> NaN conversion */
//    if (ieee_check != ieee_check)
//    {
//       /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
//          for ieee_check self-equality works on all IEEE-compliant compilers/
//          machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
//          by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
//          found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
//       if (print_level > 0 || logging > 0)
//       {
//         hypre_printf("\n\nERROR detected by Hypre ...  BEGIN\n");
//         hypre_printf("ERROR -- hypre_PCGSolve: INFs and/or NaNs detected in input.\n");
//         hypre_printf("User probably placed non-numerics in supplied b.\n");
//         hypre_printf("Returning error flag += 101.  Program not terminated.\n");
//         hypre_printf("ERROR detected by Hypre ...  END\n\n\n");
//       }
//       hypre_error(HYPRE_ERROR_GENERIC);
//       return hypre_error_flag;
//    }

//    eps = r_tol*r_tol; /* note: this may be re-assigned below */
//    if ( bi_prod > 0.0 ) {
//       if ( stop_crit && !rel_change && atolf<=0 ) {  /* pure absolute tolerance */
//          eps = eps / bi_prod;
//          /* Note: this section is obsolete.  Aside from backwards comatability
//             concerns, we could delete the stop_crit parameter and related code,
//             using tol & atolf instead. */
//       }
//       else if ( atolf>0 )  /* mixed relative and absolute tolerance */
//          bi_prod += atolf;
//       else /* DEFAULT (stop_crit and atolf exist for backwards compatibilty
//               and are not in the reference manual) */
//       {
//         /* convergence criteria:  <C*r,r>  <= max( a_tol^2, r_tol^2 * <C*b,b> )
//             note: default for a_tol is 0.0, so relative residual criteria is used unless
//             user specifies a_tol, or sets r_tol = 0.0, which means absolute
//             tol only is checked  */
//          eps = hypre_max(r_tol*r_tol, a_tol*a_tol/bi_prod);
         
//       }
//    }
//    else    /* bi_prod==0.0: the rhs vector b is zero */
//    {
//       /* Set x equal to zero and return */
//       (*(pcg_functions->CopyVector))(b, x);
//       if (logging>0 || print_level>0)
//       {
//          norms[0]     = 0.0;
//          rel_norms[i] = 0.0;
//       }

//       return hypre_error_flag;
//       /* In this case, for the original parcsr pcg, the code would take special
//          action to force iterations even though the exact value was known. */
//    };

//    /* r = b - Ax */
//    (*(pcg_functions->CopyVector))(b, r);
//    (*(pcg_functions->Matvec))(matvec_data, -1.0, A, x, 1.0, r);
 
//    /* p = C*r */
//    (*(pcg_functions->ClearVector))(p);
//    precond(precond_data, A, r, p);

//    /* gamma = <r,p> */
//    gamma = (*(pcg_functions->InnerProd))(r,p);

//    /* Since it is does not diminish performance, attempt to return an error flag
//       and notify users when they supply bad input. */
//    if (gamma != 0.) ieee_check = gamma/gamma; /* INF -> NaN conversion */
//    if (ieee_check != ieee_check)
//    {
//       /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
//          for ieee_check self-equality works on all IEEE-compliant compilers/
//          machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
//          by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
//          found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
//       if (print_level > 0 || logging > 0)
//       {
//         hypre_printf("\n\nERROR detected by Hypre ...  BEGIN\n");
//         hypre_printf("ERROR -- hypre_PCGSolve: INFs and/or NaNs detected in input.\n");
//         hypre_printf("User probably placed non-numerics in supplied A or x_0.\n");
//         hypre_printf("Returning error flag += 101.  Program not terminated.\n");
//         hypre_printf("ERROR detected by Hypre ...  END\n\n\n");
//       }
//       hypre_error(HYPRE_ERROR_GENERIC);
//       return hypre_error_flag;
//    }

//    /* Set initial residual norm */
//    if ( logging>0 || print_level > 0 || cf_tol > 0.0 )
//    {
//       if (two_norm)
//          i_prod_0 = (*(pcg_functions->InnerProd))(r,r);
//       else
//          i_prod_0 = gamma;

//       if ( logging>0 || print_level>0 ) norms[0] = sqrt(i_prod_0);
//    }
//    if ( print_level > 1 && my_id==0 )
//    {
//       hypre_printf("\n\n");
//       if (two_norm)
//       {
//          if ( stop_crit && !rel_change && atolf==0 ) {  /* pure absolute tolerance */
//             hypre_printf("Iters       ||r||_2     conv.rate\n");
//             hypre_printf("-----    ------------   ---------\n");
//          }
//          else {
//             hypre_printf("Iters       ||r||_2     conv.rate  ||r||_2/||b||_2\n");
//             hypre_printf("-----    ------------   ---------  ------------ \n");
//          }
//       }
//       else  /* !two_norm */
//       {
//          hypre_printf("Iters       ||r||_C     conv.rate  ||r||_C/||b||_C\n");
//          hypre_printf("-----    ------------    ---------  ------------ \n");
//       }
//    }

//    while ((i+1) <= max_iter)
//    {
//       /*--------------------------------------------------------------------
//        * the core CG calculations...
//        *--------------------------------------------------------------------*/
//       i++;

//       /* At user request, periodically recompute the residual from the formula
//          r = b - A x (instead of using the recursive definition). Note that this
//          is potentially expensive and can lead to degraded convergence (since it
//          essentially a "restarted CG"). */
//       recompute_true_residual = recompute_residual_p && !(i%recompute_residual_p);

//       /* s = A*p */
//       (*(pcg_functions->Matvec))(matvec_data, 1.0, A, p, 0.0, s);

//       /* alpha = gamma / <s,p> */
//       sdotp = (*(pcg_functions->InnerProd))(s, p);
//       if ( sdotp==0.0 )
//       {
//          /* ++ierr;*/
//          if (i==1) i_prod=i_prod_0;
//          break;
//       }
//       alpha = gamma / sdotp;

//       gamma_old = gamma;

//       /* x = x + alpha*p */
//       (*(pcg_functions->Axpy))(alpha, p, x);

//       /* r = r - alpha*s */
//       if ( !recompute_true_residual )
//       {
//          (*(pcg_functions->Axpy))(-alpha, s, r);
//       }
//       else
//       {
//          if (print_level > 1 && my_id == 0)
//          {
//             hypre_printf("Recomputing the residual...\n");
//          }
//          (*(pcg_functions->CopyVector))(b, r);
//          (*(pcg_functions->Matvec))(matvec_data, -1.0, A, x, 1.0, r);
//       }

//       /* residual-based stopping criteria: ||r_new-r_old|| < rtol ||b|| */
//       if (rtol && two_norm)
//       {
//          /* use that r_new-r_old = alpha * s */
//          HYPRE_Real drob2 = alpha*alpha*(*(pcg_functions->InnerProd))(s,s)/bi_prod;
//          if ( drob2 < rtol*rtol )
//          {
//             if (print_level > 1 && my_id == 0)
//             {
//                hypre_printf("\n\n||r_old-r_new||/||b||: %e\n", sqrt(drob2));
//             }
//             break;
//          }
//       }

//       /* s = C*r */
//       (*(pcg_functions->ClearVector))(s);
//       precond(precond_data, A, r, s);

//       /* gamma = <r,s> */
//       gamma = (*(pcg_functions->InnerProd))(r, s);

//       /* residual-based stopping criteria: ||r_new-r_old||_C < rtol ||b||_C */
//       if (rtol && !two_norm)
//       {
//          /* use that ||r_new-r_old||_C^2 = (r_new ,C r_new) + (r_old, C r_old) */
//          HYPRE_Real r2ob2 = (gamma + gamma_old)/bi_prod;
//          if ( r2ob2 < rtol*rtol)
//          {
//             if (print_level > 1 && my_id == 0)
//             {
//                hypre_printf("\n\n||r_old-r_new||_C/||b||_C: %e\n", sqrt(r2ob2));
//             }
//             break;
//          }
//       }

//       /* set i_prod for convergence test */
//       if (two_norm)
//          i_prod = (*(pcg_functions->InnerProd))(r,r);
//       else
//          i_prod = gamma;

//       /*--------------------------------------------------------------------
//        * optional output
//        *--------------------------------------------------------------------*/
// #if 0
//       if (two_norm)
//          hypre_printf("Iter (%d): ||r||_2 = %e, ||r||_2/||b||_2 = %e\n",
//                 i, sqrt(i_prod), (bi_prod ? sqrt(i_prod/bi_prod) : 0));
//       else
//          hypre_printf("Iter (%d): ||r||_C = %e, ||r||_C/||b||_C = %e\n",
//                 i, sqrt(i_prod), (bi_prod ? sqrt(i_prod/bi_prod) : 0));
// #endif
 
//       /* print norm info */
//       if ( logging>0 || print_level>0 )
//       {
//          norms[i]     = sqrt(i_prod);
//          rel_norms[i] = bi_prod ? sqrt(i_prod/bi_prod) : 0;
//       }
//       if ( print_level > 1 && my_id==0 )
//       {
//          if (two_norm)
//          {
//             if ( stop_crit && !rel_change && atolf==0 ) {  /* pure absolute tolerance */
//                hypre_printf("% 5d    %e    %f\n", i, norms[i],
//                       norms[i]/norms[i-1] );
//             }
//             else 
//             {
//                hypre_printf("% 5d    %e    %f    %e\n", i, norms[i],
//                       norms[i]/norms[i-1], rel_norms[i] );
//             }
//          }
//          else 
//          {
//                hypre_printf("% 5d    %e    %f    %e\n", i, norms[i],
//                       norms[i]/norms[i-1], rel_norms[i] );
//          }
//       }


//       /*--------------------------------------------------------------------
//        * check for convergence
//        *--------------------------------------------------------------------*/
//       if (i_prod / bi_prod < eps)  /* the basic convergence test */
//             tentatively_converged = 1;
//       if ( tentatively_converged && recompute_residual )
//          /* At user request, don't trust the convergence test until we've recomputed
//             the residual from scratch.  This is expensive in the usual case where an
//             the norm is the energy norm.
//             This calculation is coded on the assumption that r's accuracy is only a
//             concern for problems where CG takes many iterations. */
//       {
//          /* r = b - Ax */
//          (*(pcg_functions->CopyVector))(b, r);
//          (*(pcg_functions->Matvec))(matvec_data, -1.0, A, x, 1.0, r);
//          // (*(pcg_functions->MatvecT))(matvec_data, -1.0, A, x, 1.0, r);

//          /* set i_prod for convergence test */
//          if (two_norm)
//             i_prod = (*(pcg_functions->InnerProd))(r,r);
//          else
//          {
//             /* s = C*r */
//             (*(pcg_functions->ClearVector))(s);
//             precond(precond_data, A, r, s);
//             /* iprod = gamma = <r,s> */
//             i_prod = (*(pcg_functions->InnerProd))(r, s);
//          }
//          if (i_prod / bi_prod >= eps) tentatively_converged = 0;
//       }
//       if ( tentatively_converged && rel_change && (i_prod > guard_zero_residual ))
//          /* At user request, don't treat this as converged unless x didn't change
//             much in the last iteration. */
//       {
//             pi_prod = (*(pcg_functions->InnerProd))(p,p); 
//             xi_prod = (*(pcg_functions->InnerProd))(x,x);
//             ratio = alpha*alpha*pi_prod/xi_prod;
//             if (ratio >= eps) tentatively_converged = 0;
//       }
//       if ( tentatively_converged )
//          /* we've passed all the convergence tests, it's for real */
//       {
//          (pcg_data -> converged) = 1;
//          break;
//       }

//       if ( (gamma<1.0e-292) && ((-gamma)<1.0e-292) ) {
//          /* ierr = 1;*/
//          hypre_error(HYPRE_ERROR_CONV);
         
//          break;
//       }
//       /* ... gamma should be >=0.  IEEE subnormal numbers are < 2**(-1022)=2.2e-308
//          (and >= 2**(-1074)=4.9e-324).  So a gamma this small means we're getting
//          dangerously close to subnormal or zero numbers (usually if gamma is small,
//          so will be other variables).  Thus further calculations risk a crash.
//          Such small gamma generally means no hope of progress anyway. */

//       /*--------------------------------------------------------------------
//        * Optional test to see if adequate progress is being made.
//        * The average convergence factor is recorded and compared
//        * against the tolerance 'cf_tol'. The weighting factor is  
//        * intended to pay more attention to the test when an accurate
//        * estimate for average convergence factor is available.  
//        *--------------------------------------------------------------------*/

//       if (cf_tol > 0.0)
//       {
//          cf_ave_0 = cf_ave_1;
//          if ( i_prod_0<1.0e-292 ) {
//             /* i_prod_0 is zero, or (almost) subnormal, yet i_prod wasn't small
//                enough to pass the convergence test.  Therefore initial guess was good,
//                and we're just calculating garbage - time to bail out before the
//                next step, which will be a divide by zero (or close to it). */
//             /* ierr = 1; */
//             hypre_error(HYPRE_ERROR_CONV);
            
//             break;
//          }
//          cf_ave_1 = pow( i_prod / i_prod_0, 1.0/(2.0*i)); 

//          weight   = fabs(cf_ave_1 - cf_ave_0);
//          weight   = weight / hypre_max(cf_ave_1, cf_ave_0);
//          weight   = 1.0 - weight;
// #if 0
//          hypre_printf("I = %d: cf_new = %e, cf_old = %e, weight = %e\n",
//                 i, cf_ave_1, cf_ave_0, weight );
// #endif
//          if (weight * cf_ave_1 > cf_tol) break;
//       }

//       /*--------------------------------------------------------------------
//        * back to the core CG calculations
//        *--------------------------------------------------------------------*/

//       /* beta = gamma / gamma_old */
//       beta = gamma / gamma_old;

//       /* p = s + beta p */
//       if ( !recompute_true_residual )
//       {
//          (*(pcg_functions->ScaleVector))(beta, p);
//          (*(pcg_functions->Axpy))(1.0, s, p);
//       }
//       else
//          (*(pcg_functions->CopyVector))(s, p);
//    }

//    /*--------------------------------------------------------------------
//     * Finish up with some outputs.
//     *--------------------------------------------------------------------*/

//    if ( print_level > 1 && my_id==0 )
//       hypre_printf("\n\n");

//    (pcg_data -> num_iterations) = i;
//    if (bi_prod > 0.0)
//       (pcg_data -> rel_residual_norm) = sqrt(i_prod/bi_prod);
//    else /* actually, we'll never get here... */
//       (pcg_data -> rel_residual_norm) = 0.0;

//    return hypre_error_flag;
// }

// HYPRE_Int 
// HYPRE_PCGSolveArthurTRUE( HYPRE_Solver solver,
//                 HYPRE_Matrix K,
//                 HYPRE_Matrix G,
//                 HYPRE_Vector f,
//                 HYPRE_Vector h,
//                 HYPRE_Vector vel,
//                 HYPRE_Vector pres )
// {
//    return( hypre_PCGSolveArthurTRUE( (void *) solver,
//                            (void *) K,
//                            (void *) G,
//                            (void *) f,
//                            (void *) h,
//                            (void *) vel,
//                            (void *) pres ) );
// }

// HYPRE_Int 
// HYPRE_ParCSRPCGSolveArthurTRUE( HYPRE_Solver solver,
//                       HYPRE_ParCSRMatrix K,
//                       HYPRE_ParCSRMatrix G,
//                       HYPRE_ParVector f,
//                       HYPRE_ParVector h,
//                       HYPRE_ParVector vel,
//                       HYPRE_ParVector pres      )
// {
//    return( HYPRE_PCGSolveArthurTRUE( solver,
//                            (HYPRE_Matrix) K,
//                            (HYPRE_Matrix) G,
//                            (HYPRE_Vector) f,
//                            (HYPRE_Vector) h,
//                            (HYPRE_Vector) vel,
//                            (HYPRE_Vector) pres )  );
// }









// // IN CASE 
// HYPRE_Int
// hypre_PCGSetupArthur( void *pcg_vdata,
//                 void *A,
//                 void *K,
//                 void *b,
//                 void *x         )
// {
//    hypre_PCGData *pcg_data =  (hypre_PCGData *)pcg_vdata;
//    hypre_PCGFunctions *pcg_functions = pcg_data->functions;
//    HYPRE_Int            max_iter         = (pcg_data -> max_iter);
//    HYPRE_Int          (*precond_setup)(void*,void*,void*,void*) = (pcg_functions -> precond_setup);
//    void          *precond_data     = (pcg_data -> precond_data);


//    (pcg_data -> A) = A;

//    /*--------------------------------------------------
//     * The arguments for CreateVector are important to
//     * maintain consistency between the setup and
//     * compute phases of matvec and the preconditioner.
//     *--------------------------------------------------*/

//    if ( pcg_data -> p != NULL )
//       (*(pcg_functions->DestroyVector))(pcg_data -> p);
//    (pcg_data -> p) = (*(pcg_functions->CreateVector))(x);

//    if ( pcg_data -> s != NULL )
//       (*(pcg_functions->DestroyVector))(pcg_data -> s);
//    (pcg_data -> s) = (*(pcg_functions->CreateVector))(x);

//    if ( pcg_data -> r != NULL )
//       (*(pcg_functions->DestroyVector))(pcg_data -> r);
//    (pcg_data -> r) = (*(pcg_functions->CreateVector))(b);

//    if ( pcg_data -> matvec_data != NULL && pcg_data->owns_matvec_data )
//       (*(pcg_functions->MatvecDestroy))(pcg_data -> matvec_data);
//    (pcg_data -> matvec_data) = (*(pcg_functions->MatvecCreate))(A, x);

//    precond_setup(precond_data, A, b, x);

//    /*-----------------------------------------------------
//     * Allocate space for log info
//     *-----------------------------------------------------*/

//    if ( (pcg_data->logging)>0  || (pcg_data->print_level)>0 ) 
//    {
//       if ( (pcg_data -> norms) != NULL )
//          hypre_TFreeF( pcg_data -> norms, pcg_functions );
//       (pcg_data -> norms)     = hypre_CTAllocF( HYPRE_Real, max_iter + 1,
//                                                 pcg_functions);

//       if ( (pcg_data -> rel_norms) != NULL )
//          hypre_TFreeF( pcg_data -> rel_norms, pcg_functions );
//       (pcg_data -> rel_norms) = hypre_CTAllocF( HYPRE_Real, max_iter + 1,
//                                                 pcg_functions );
//    }

//    return hypre_error_flag;
// }

// HYPRE_Int 
// HYPRE_PCGSetupArthur( HYPRE_Solver solver,
//                 HYPRE_Matrix A,
//                 HYPRE_Matrix K,
//                 HYPRE_Vector b,
//                 HYPRE_Vector x )
// {
//    return( hypre_PCGSetupArthur( solver,
//                            A,
//                            K,
//                            b,
//                            x ) );
// }

// HYPRE_Int 
// HYPRE_ParCSRPCGSetupArthur2( HYPRE_Solver solver,
//                       HYPRE_ParCSRMatrix A,
//                       HYPRE_ParCSRMatrix K,
//                       HYPRE_ParVector b,
//                       HYPRE_ParVector x      )
// {
//    return( HYPRE_PCGSetupArthur( solver,
//                            (HYPRE_Matrix) A,
//                            (HYPRE_Matrix) K,
//                            (HYPRE_Vector) b,
//                            (HYPRE_Vector) x ) );
// }


// /*--------------------------------------------------------------------------
//  * hypre_PCGSolve
//  *--------------------------------------------------------------------------
//  *
//  * We use the following convergence test as the default (see Ashby, Holst,
//  * Manteuffel, and Saylor):
//  *
//  *       ||e||_A                           ||r||_C
//  *       -------  <=  [kappa_A(C*A)]^(1/2) -------  < tol
//  *       ||x||_A                           ||b||_C
//  *
//  * where we let (for the time being) kappa_A(CA) = 1.
//  * We implement the test as:
//  *
//  *       gamma = <C*r,r>/<C*b,b>  <  (tol^2) = eps
//  *
//  *--------------------------------------------------------------------------*/

// HYPRE_Int
// hypre_PCGSolveArthur( void *pcg_vdata,
//                 void *A,
//                 void *K,
//                 void *b,
//                 void *x         )
// //void *K,
// //void *G,
// //void *f,
// //void *h,
// //void *vel,
// //void *pres
// {
//    hypre_PCGData  *pcg_data     =  (hypre_PCGData *)pcg_vdata;
//    hypre_PCGFunctions *pcg_functions = pcg_data->functions;

//    HYPRE_Real      r_tol        = (pcg_data -> tol);
//    HYPRE_Real      a_tol        = (pcg_data -> a_tol);
//    HYPRE_Real      atolf        = (pcg_data -> atolf);
//    HYPRE_Real      cf_tol       = (pcg_data -> cf_tol);
//    HYPRE_Real      rtol         = (pcg_data -> rtol);
//    HYPRE_Int             max_iter     = (pcg_data -> max_iter);
//    HYPRE_Int             two_norm     = (pcg_data -> two_norm);
//    HYPRE_Int             rel_change   = (pcg_data -> rel_change);
//    HYPRE_Int             recompute_residual = (pcg_data -> recompute_residual);
//    HYPRE_Int             recompute_residual_p = (pcg_data -> recompute_residual_p);
//    HYPRE_Int             stop_crit    = (pcg_data -> stop_crit);
// /*
//    HYPRE_Int             converged    = (pcg_data -> converged);
// */
//    void           *p            = (pcg_data -> p);
//    void           *s            = (pcg_data -> s);
//    void           *r            = (pcg_data -> r);
//    void           *matvec_data  = (pcg_data -> matvec_data);
//    HYPRE_Int           (*precond)(void*,void*,void*,void*)   = (pcg_functions -> precond);
//    void           *precond_data = (pcg_data -> precond_data);
//    HYPRE_Int             print_level  = (pcg_data -> print_level);
//    HYPRE_Int             logging      = (pcg_data -> logging);
//    HYPRE_Real     *norms        = (pcg_data -> norms);
//    HYPRE_Real     *rel_norms    = (pcg_data -> rel_norms);
                
//    HYPRE_Real      alpha, beta;
//    HYPRE_Real      gamma, gamma_old;
//    HYPRE_Real      bi_prod, eps;
//    HYPRE_Real      pi_prod, xi_prod;
//    HYPRE_Real      ieee_check = 0.;
                
//    HYPRE_Real      i_prod = 0.0;
//    HYPRE_Real      i_prod_0 = 0.0;
//    HYPRE_Real      cf_ave_0 = 0.0;
//    HYPRE_Real      cf_ave_1 = 0.0;
//    HYPRE_Real      weight;
//    HYPRE_Real      ratio;

//    HYPRE_Real      guard_zero_residual, sdotp;
//    HYPRE_Int             tentatively_converged = 0;
//    HYPRE_Int             recompute_true_residual = 0;

//    HYPRE_Int             i = 0;
//    HYPRE_Int             my_id, num_procs;

//    (pcg_data -> converged) = 0;

//    (*(pcg_functions->CommInfo))(A,&my_id,&num_procs);

//    /*-----------------------------------------------------------------------
//     * With relative change convergence test on, it is possible to attempt
//     * another iteration with a zero residual. This causes the parameter
//     * alpha to go NaN. The guard_zero_residual parameter is to circumvent
//     * this. Perhaps it should be set to something non-zero (but small).
//     *-----------------------------------------------------------------------*/

//    guard_zero_residual = 0.0;

//    /*-----------------------------------------------------------------------
//     * Start pcg solve
//     *-----------------------------------------------------------------------*/

//    /* compute eps */
//    if (two_norm)
//    {
//       /* bi_prod = <b,b> */
//       bi_prod = (*(pcg_functions->InnerProd))(b, b);
//       if (print_level > 1 && my_id == 0) 
//           hypre_printf("<b,b>: %e\n",bi_prod);
//    }
//    else
//    {
//       /* bi_prod = <C*b,b> */
//       (*(pcg_functions->ClearVector))(p);
//       precond(precond_data, A, b, p);
//       bi_prod = (*(pcg_functions->InnerProd))(p, b);
//       if (print_level > 1 && my_id == 0)
//           hypre_printf("<C*b,b>: %e\n",bi_prod);
//    };

//    /* Since it is does not diminish performance, attempt to return an error flag
//       and notify users when they supply bad input. */
//    if (bi_prod != 0.) ieee_check = bi_prod/bi_prod; /* INF -> NaN conversion */
//    if (ieee_check != ieee_check)
//    {
//       /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
//          for ieee_check self-equality works on all IEEE-compliant compilers/
//          machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
//          by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
//          found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
//       if (print_level > 0 || logging > 0)
//       {
//         hypre_printf("\n\nERROR detected by Hypre ...  BEGIN\n");
//         hypre_printf("ERROR -- hypre_PCGSolve: INFs and/or NaNs detected in input.\n");
//         hypre_printf("User probably placed non-numerics in supplied b.\n");
//         hypre_printf("Returning error flag += 101.  Program not terminated.\n");
//         hypre_printf("ERROR detected by Hypre ...  END\n\n\n");
//       }
//       hypre_error(HYPRE_ERROR_GENERIC);
//       return hypre_error_flag;
//    }

//    eps = r_tol*r_tol; /* note: this may be re-assigned below */
//    if ( bi_prod > 0.0 ) {
//       if ( stop_crit && !rel_change && atolf<=0 ) {  /* pure absolute tolerance */
//          eps = eps / bi_prod;
//          /* Note: this section is obsolete.  Aside from backwards comatability
//             concerns, we could delete the stop_crit parameter and related code,
//             using tol & atolf instead. */
//       }
//       else if ( atolf>0 )  /* mixed relative and absolute tolerance */
//          bi_prod += atolf;
//       else /* DEFAULT (stop_crit and atolf exist for backwards compatibilty
//               and are not in the reference manual) */
//       {
//         /* convergence criteria:  <C*r,r>  <= max( a_tol^2, r_tol^2 * <C*b,b> )
//             note: default for a_tol is 0.0, so relative residual criteria is used unless
//             user specifies a_tol, or sets r_tol = 0.0, which means absolute
//             tol only is checked  */
//          eps = hypre_max(r_tol*r_tol, a_tol*a_tol/bi_prod);
         
//       }
//    }
//    else    /* bi_prod==0.0: the rhs vector b is zero */
//    {
//       /* Set x equal to zero and return */
//       (*(pcg_functions->CopyVector))(b, x);
//       if (logging>0 || print_level>0)
//       {
//          norms[0]     = 0.0;
//          rel_norms[i] = 0.0;
//       }

//       return hypre_error_flag;
//       /* In this case, for the original parcsr pcg, the code would take special
//          action to force iterations even though the exact value was known. */
//    };

//    /* r = b - Ax */
//    (*(pcg_functions->CopyVector))(b, r);
//    (*(pcg_functions->Matvec))(matvec_data, -1.0, A, x, 1.0, r);
 
//    /* p = C*r */
//    (*(pcg_functions->ClearVector))(p);
//    precond(precond_data, A, r, p);

//    /* gamma = <r,p> */
//    gamma = (*(pcg_functions->InnerProd))(r,p);

//    /* Since it is does not diminish performance, attempt to return an error flag
//       and notify users when they supply bad input. */
//    if (gamma != 0.) ieee_check = gamma/gamma; /* INF -> NaN conversion */
//    if (ieee_check != ieee_check)
//    {
//       /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
//          for ieee_check self-equality works on all IEEE-compliant compilers/
//          machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
//          by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
//          found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
//       if (print_level > 0 || logging > 0)
//       {
//         hypre_printf("\n\nERROR detected by Hypre ...  BEGIN\n");
//         hypre_printf("ERROR -- hypre_PCGSolve: INFs and/or NaNs detected in input.\n");
//         hypre_printf("User probably placed non-numerics in supplied A or x_0.\n");
//         hypre_printf("Returning error flag += 101.  Program not terminated.\n");
//         hypre_printf("ERROR detected by Hypre ...  END\n\n\n");
//       }
//       hypre_error(HYPRE_ERROR_GENERIC);
//       return hypre_error_flag;
//    }

//    /* Set initial residual norm */
//    if ( logging>0 || print_level > 0 || cf_tol > 0.0 )
//    {
//       if (two_norm)
//          i_prod_0 = (*(pcg_functions->InnerProd))(r,r);
//       else
//          i_prod_0 = gamma;

//       if ( logging>0 || print_level>0 ) norms[0] = sqrt(i_prod_0);
//    }
//    if ( print_level > 1 && my_id==0 )
//    {
//       hypre_printf("\n\n");
//       if (two_norm)
//       {
//          if ( stop_crit && !rel_change && atolf==0 ) {  /* pure absolute tolerance */
//             hypre_printf("Iters       ||r||_2     conv.rate\n");
//             hypre_printf("-----    ------------   ---------\n");
//          }
//          else {
//             hypre_printf("Iters       ||r||_2     conv.rate  ||r||_2/||b||_2\n");
//             hypre_printf("-----    ------------   ---------  ------------ \n");
//          }
//       }
//       else  /* !two_norm */
//       {
//          hypre_printf("Iters       ||r||_C     conv.rate  ||r||_C/||b||_C\n");
//          hypre_printf("-----    ------------    ---------  ------------ \n");
//       }
//    }

//    while ((i+1) <= max_iter)
//    {
//       /*--------------------------------------------------------------------
//        * the core CG calculations...
//        *--------------------------------------------------------------------*/
//       i++;

//       /* At user request, periodically recompute the residual from the formula
//          r = b - A x (instead of using the recursive definition). Note that this
//          is potentially expensive and can lead to degraded convergence (since it
//          essentially a "restarted CG"). */
//       recompute_true_residual = recompute_residual_p && !(i%recompute_residual_p);

//       /* s = A*p */
//       (*(pcg_functions->Matvec))(matvec_data, 1.0, A, p, 0.0, s);

//       /* alpha = gamma / <s,p> */
//       sdotp = (*(pcg_functions->InnerProd))(s, p);
//       if ( sdotp==0.0 )
//       {
//          /* ++ierr;*/
//          if (i==1) i_prod=i_prod_0;
//          break;
//       }
//       alpha = gamma / sdotp;

//       gamma_old = gamma;

//       /* x = x + alpha*p */
//       (*(pcg_functions->Axpy))(alpha, p, x);

//       /* r = r - alpha*s */
//       if ( !recompute_true_residual )
//       {
//          (*(pcg_functions->Axpy))(-alpha, s, r);
//       }
//       else
//       {
//          if (print_level > 1 && my_id == 0)
//          {
//             hypre_printf("Recomputing the residual...\n");
//          }
//          (*(pcg_functions->CopyVector))(b, r);
//          (*(pcg_functions->Matvec))(matvec_data, -1.0, A, x, 1.0, r);
//       }

//       /* residual-based stopping criteria: ||r_new-r_old|| < rtol ||b|| */
//       if (rtol && two_norm)
//       {
//          /* use that r_new-r_old = alpha * s */
//          HYPRE_Real drob2 = alpha*alpha*(*(pcg_functions->InnerProd))(s,s)/bi_prod;
//          if ( drob2 < rtol*rtol )
//          {
//             if (print_level > 1 && my_id == 0)
//             {
//                hypre_printf("\n\n||r_old-r_new||/||b||: %e\n", sqrt(drob2));
//             }
//             break;
//          }
//       }

//       /* s = C*r */
//       (*(pcg_functions->ClearVector))(s);
//       precond(precond_data, A, r, s);

//       /* gamma = <r,s> */
//       gamma = (*(pcg_functions->InnerProd))(r, s);

//       /* residual-based stopping criteria: ||r_new-r_old||_C < rtol ||b||_C */
//       if (rtol && !two_norm)
//       {
//          /* use that ||r_new-r_old||_C^2 = (r_new ,C r_new) + (r_old, C r_old) */
//          HYPRE_Real r2ob2 = (gamma + gamma_old)/bi_prod;
//          if ( r2ob2 < rtol*rtol)
//          {
//             if (print_level > 1 && my_id == 0)
//             {
//                hypre_printf("\n\n||r_old-r_new||_C/||b||_C: %e\n", sqrt(r2ob2));
//             }
//             break;
//          }
//       }

//       /* set i_prod for convergence test */
//       if (two_norm)
//          i_prod = (*(pcg_functions->InnerProd))(r,r);
//       else
//          i_prod = gamma;

//       /*--------------------------------------------------------------------
//        * optional output
//        *--------------------------------------------------------------------*/
// #if 0
//       if (two_norm)
//          hypre_printf("Iter (%d): ||r||_2 = %e, ||r||_2/||b||_2 = %e\n",
//                 i, sqrt(i_prod), (bi_prod ? sqrt(i_prod/bi_prod) : 0));
//       else
//          hypre_printf("Iter (%d): ||r||_C = %e, ||r||_C/||b||_C = %e\n",
//                 i, sqrt(i_prod), (bi_prod ? sqrt(i_prod/bi_prod) : 0));
// #endif
 
//       /* print norm info */
//       if ( logging>0 || print_level>0 )
//       {
//          norms[i]     = sqrt(i_prod);
//          rel_norms[i] = bi_prod ? sqrt(i_prod/bi_prod) : 0;
//       }
//       if ( print_level > 1 && my_id==0 )
//       {
//          if (two_norm)
//          {
//             if ( stop_crit && !rel_change && atolf==0 ) {  /* pure absolute tolerance */
//                hypre_printf("% 5d    %e    %f\n", i, norms[i],
//                       norms[i]/norms[i-1] );
//             }
//             else 
//             {
//                hypre_printf("% 5d    %e    %f    %e\n", i, norms[i],
//                       norms[i]/norms[i-1], rel_norms[i] );
//             }
//          }
//          else 
//          {
//                hypre_printf("% 5d    %e    %f    %e\n", i, norms[i],
//                       norms[i]/norms[i-1], rel_norms[i] );
//          }
//       }


//       /*--------------------------------------------------------------------
//        * check for convergence
//        *--------------------------------------------------------------------*/
//       if (i_prod / bi_prod < eps)  /* the basic convergence test */
//             tentatively_converged = 1;
//       if ( tentatively_converged && recompute_residual )
//          /* At user request, don't trust the convergence test until we've recomputed
//             the residual from scratch.  This is expensive in the usual case where an
//             the norm is the energy norm.
//             This calculation is coded on the assumption that r's accuracy is only a
//             concern for problems where CG takes many iterations. */
//       {
//          /* r = b - Ax */
//          (*(pcg_functions->CopyVector))(b, r);
//          (*(pcg_functions->Matvec))(matvec_data, -1.0, A, x, 1.0, r);

//          /* set i_prod for convergence test */
//          if (two_norm)
//             i_prod = (*(pcg_functions->InnerProd))(r,r);
//          else
//          {
//             /* s = C*r */
//             (*(pcg_functions->ClearVector))(s);
//             precond(precond_data, A, r, s);
//             /* iprod = gamma = <r,s> */
//             i_prod = (*(pcg_functions->InnerProd))(r, s);
//          }
//          if (i_prod / bi_prod >= eps) tentatively_converged = 0;
//       }
//       if ( tentatively_converged && rel_change && (i_prod > guard_zero_residual ))
//          /* At user request, don't treat this as converged unless x didn't change
//             much in the last iteration. */
//       {
//             pi_prod = (*(pcg_functions->InnerProd))(p,p); 
//             xi_prod = (*(pcg_functions->InnerProd))(x,x);
//             ratio = alpha*alpha*pi_prod/xi_prod;
//             if (ratio >= eps) tentatively_converged = 0;
//       }
//       if ( tentatively_converged )
//          /* we've passed all the convergence tests, it's for real */
//       {
//          (pcg_data -> converged) = 1;
//          break;
//       }

//       if ( (gamma<1.0e-292) && ((-gamma)<1.0e-292) ) {
//          /* ierr = 1;*/
//          hypre_error(HYPRE_ERROR_CONV);
         
//          break;
//       }
//       /* ... gamma should be >=0.  IEEE subnormal numbers are < 2**(-1022)=2.2e-308
//          (and >= 2**(-1074)=4.9e-324).  So a gamma this small means we're getting
//          dangerously close to subnormal or zero numbers (usually if gamma is small,
//          so will be other variables).  Thus further calculations risk a crash.
//          Such small gamma generally means no hope of progress anyway. */

//       /*--------------------------------------------------------------------
//        * Optional test to see if adequate progress is being made.
//        * The average convergence factor is recorded and compared
//        * against the tolerance 'cf_tol'. The weighting factor is  
//        * intended to pay more attention to the test when an accurate
//        * estimate for average convergence factor is available.  
//        *--------------------------------------------------------------------*/

//       if (cf_tol > 0.0)
//       {
//          cf_ave_0 = cf_ave_1;
//          if ( i_prod_0<1.0e-292 ) {
//             /* i_prod_0 is zero, or (almost) subnormal, yet i_prod wasn't small
//                enough to pass the convergence test.  Therefore initial guess was good,
//                and we're just calculating garbage - time to bail out before the
//                next step, which will be a divide by zero (or close to it). */
//             /* ierr = 1; */
//             hypre_error(HYPRE_ERROR_CONV);
            
//             break;
//          }
//          cf_ave_1 = pow( i_prod / i_prod_0, 1.0/(2.0*i)); 

//          weight   = fabs(cf_ave_1 - cf_ave_0);
//          weight   = weight / hypre_max(cf_ave_1, cf_ave_0);
//          weight   = 1.0 - weight;
// #if 0
//          hypre_printf("I = %d: cf_new = %e, cf_old = %e, weight = %e\n",
//                 i, cf_ave_1, cf_ave_0, weight );
// #endif
//          if (weight * cf_ave_1 > cf_tol) break;
//       }

//       /*--------------------------------------------------------------------
//        * back to the core CG calculations
//        *--------------------------------------------------------------------*/

//       /* beta = gamma / gamma_old */
//       beta = gamma / gamma_old;

//       /* p = s + beta p */
//       if ( !recompute_true_residual )
//       {
//          (*(pcg_functions->ScaleVector))(beta, p);
//          (*(pcg_functions->Axpy))(1.0, s, p);
//       }
//       else
//          (*(pcg_functions->CopyVector))(s, p);
//    }

//    /*--------------------------------------------------------------------
//     * Finish up with some outputs.
//     *--------------------------------------------------------------------*/

//    if ( print_level > 1 && my_id==0 )
//       hypre_printf("\n\n");

//    (pcg_data -> num_iterations) = i;
//    if (bi_prod > 0.0)
//       (pcg_data -> rel_residual_norm) = sqrt(i_prod/bi_prod);
//    else /* actually, we'll never get here... */
//       (pcg_data -> rel_residual_norm) = 0.0;

//    return hypre_error_flag;
// }


// HYPRE_Int 
// HYPRE_PCGSolveArthur( HYPRE_Solver solver,
//                 HYPRE_Matrix A,
//                 HYPRE_Matrix K,
//                 HYPRE_Vector b,
//                 HYPRE_Vector x )
// {
//    return( hypre_PCGSolveArthur( (void *) solver,
//                            (void *) A,
//                            (void *) K,
//                            (void *) b,
//                            (void *) x ) );
// }

// HYPRE_Int 
// HYPRE_ParCSRPCGSolveArthur2( HYPRE_Solver solver,
//                       HYPRE_ParCSRMatrix A,
//                       HYPRE_ParCSRMatrix K,
//                       HYPRE_ParVector b,
//                       HYPRE_ParVector x      )
// {
//    return( HYPRE_PCGSolveArthur( solver,
//                            (HYPRE_Matrix) A,
//                            (HYPRE_Matrix) K,
//                            (HYPRE_Vector) b,
//                            (HYPRE_Vector) x ) );
// }
