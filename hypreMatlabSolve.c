/*
 Read a linear system assembled in Matlab
 Solve the linear system using hypre
 Interface: Linear-Algebraic (IJ)
 Available solvers: 0  - AMG (default)
                    1  - AMG-PCG
                    8  - ParaSails-PCG
                    50 - PCG
                    61 - AMG-FlexGMRES

 Copyright: M. Giacomini (2017)
 (based on ex5.c in the hypre examples folder)
*/

#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

#include "vis.c"

#include "arthur.h"

int main (int argc, char *argv[])
{
   HYPRE_Int i;
   HYPRE_Int iErr = 0;
    
   HYPRE_Int myid, num_procs, dummy;
   HYPRE_Int solver_id;
    
   HYPRE_Int first_local_row, last_local_row;//, local_num_rows;
   HYPRE_Int first_local_col, last_local_col, local_num_cols;
    
   HYPRE_Real *values;

   HYPRE_IJMatrix ij_A;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJMatrix ij_G;
   HYPRE_ParCSRMatrix parcsr_G;
   //HYPRE_IJMatrix ij_Gprime;
   //HYPRE_ParCSRMatrix parcsr_Gprime;
   HYPRE_IJVector ij_f;
   HYPRE_ParVector par_f;
   HYPRE_IJVector ij_vel;
   HYPRE_ParVector par_vel;
   //ADD g, and pres with dimension of G when G can be readable.

   char saveMatsAs;
   
   FILE *fset;
    
   void *object;
    
   HYPRE_Solver solver;//, precond;
    

    /* Read dimension of the system and chosen solver from file */
    fset = fopen("solverPar.dat","r"); // read mode
    if( fset == NULL )
    {
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
    hypre_fscanf(fset, "%d\n", &solver_id);
    hypre_fscanf(fset, "%d\n", &num_procs);
    hypre_fscanf(fset, "%d\n", &dummy);
    hypre_fscanf(fset, "%c\n", &saveMatsAs);
    fclose(fset);

    
    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    
   /* Read the matrix previously assembled
      <filename>  = IJ.A.out to read in what has been printed out (processor numbers are omitted). */
    iErr = HYPRE_IJMatrixRead( "matrixA", MPI_COMM_WORLD, HYPRE_PARCSR, &ij_A );


    if (iErr) 
    {
        hypre_printf("ERROR: Problem reading in the system matrix!\n");
        exit(1);
    }
    /* Get dimension info */
    iErr = HYPRE_IJMatrixGetLocalRange( ij_A,
                                       &first_local_row, &last_local_row ,
                                       &first_local_col, &last_local_col );
    
    //local_num_rows = last_local_row - first_local_row + 1;
    local_num_cols = last_local_col - first_local_col + 1;
   /* Get the parcsr matrix object to use */
   iErr += HYPRE_IJMatrixGetObject( ij_A, &object);
   parcsr_A = (HYPRE_ParCSRMatrix) object;

   /* Read the matrix previously assembled
      <filename>  = IJ.A.out to read in what has been printed out (processor numbers are omitted). */
    iErr = HYPRE_IJMatrixRead( "matrixA", MPI_COMM_WORLD, HYPRE_PARCSR, &ij_G );


    if (iErr) 
    {
        hypre_printf("ERROR: Problem reading in the system matrix!\n");
        exit(1);
    }
    /* Get dimension info */
    iErr = HYPRE_IJMatrixGetLocalRange( ij_G,
                                       &first_local_row, &last_local_row ,
                                       &first_local_col, &last_local_col );
    
    //local_num_rows = last_local_row - first_local_row + 1;
    local_num_cols = last_local_col - first_local_col + 1;
   /* Get the parcsr matrix object to use */
   iErr += HYPRE_IJMatrixGetObject( ij_G, &object);
   parcsr_G = (HYPRE_ParCSRMatrix) object;
    

   /* Read the matrix previously assembled
      <filename>  = IJ.A.out to read in what has been printed out (processor numbers are omitted). */
   //  iErr = HYPRE_IJMatrixRead( "matrixA", MPI_COMM_WORLD, HYPRE_PARCSR, &ij_Gprime );


   //  if (iErr) 
   //  {
   //      hypre_printf("ERROR: Problem reading in the system matrix!\n");
   //      exit(1);
   //  }
   //  /* Get dimension info */
   //  iErr = HYPRE_IJMatrixGetLocalRange( ij_Gprime,
   //                                     &first_local_row, &last_local_row ,
   //                                     &first_local_col, &last_local_col );
    
   //  //local_num_rows = last_local_row - first_local_row + 1;
   //  local_num_cols = last_local_col - first_local_col + 1;
   // /* Get the parcsr matrix object to use */
   // iErr += HYPRE_IJMatrixGetObject( ij_Gprime, &object);
   // parcsr_Gprime = (HYPRE_ParCSRMatrix) object;
    
    
   /*  Read the RHS previously assembled */
   iErr = HYPRE_ParVectorRead(MPI_COMM_WORLD, "vectorB.0", &par_f);
   if (iErr)
   {
       hypre_printf("ERROR: Problem reading in the right-hand-side!\n");
       exit(1);
   }
   ij_f = NULL;
    
    
   /* Create the initial solution and set it to zero */
   HYPRE_IJVectorCreate(MPI_COMM_WORLD, first_local_col, last_local_col, &ij_vel);
   HYPRE_IJVectorSetObjectType(ij_vel, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(ij_vel);
   /* Initialize the guess vector */
   values = hypre_CTAlloc(HYPRE_Real, local_num_cols);
   for (i = 0; i < local_num_cols; i++)
      values[i] = 0.;
   HYPRE_IJVectorSetValues(ij_vel, local_num_cols, NULL, values);
   hypre_TFree(values);
   /* Get the parcsr vector object to use */
   iErr = HYPRE_IJVectorGetObject( ij_vel, &object );
   par_vel = (HYPRE_ParVector) object;


   /* Choose a solver and solve the system */


   /* PCG */
   if (solver_id == 50)
   {
      int num_iterations;
      double final_res_norm;

      /* Create solver */
      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

      /* Set some parameters (See Reference Manual for more parameters) */
      HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
      HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
      HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
      HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
      HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

      /* Now setup and solve! */
      HYPRE_ParCSRPCGSetupArthur2(solver, parcsr_A, parcsr_G, par_f, par_vel);
      HYPRE_ParCSRPCGSolveArthur2(solver, parcsr_A, parcsr_G, par_f, par_vel);

      //HYPRE_ParCSRPCGSetupArthurTRUE(solver, parcsr_A, parcsr_G, par_vel, par_f, par_vel, par_vel);
      //HYPRE_ParCSRPCGSolveArthurTRUE(solver, parcsr_A, parcsr_G, par_vel, par_f, par_vel, par_vel);

      /* Run info - needed logging turned on */
      HYPRE_PCGGetNumIterations(solver, &num_iterations);
      HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
      if (myid == 0)
      {
         printf("\n");
         printf("Iterations = %d\n", num_iterations);
         printf("Final Relative Residual Norm = %e\n", final_res_norm);
         printf("\n");
      }

      /* Destroy solver */
      HYPRE_ParCSRPCGDestroy(solver);
   }
   
    
   /* Save the solution to file */
   HYPRE_ParVectorPrint(par_vel, "solution.0");

    
   /* Clean up */
   HYPRE_IJMatrixDestroy(ij_A);
   HYPRE_IJMatrixDestroy(ij_G);
   //HYPRE_IJMatrixDestroy(ij_Gprime);
   HYPRE_IJVectorDestroy(ij_f);
   HYPRE_IJVectorDestroy(ij_vel);


   /* Finalize MPI*/
   if (myid == 0)
   {
       printf("\n Linear system correctly solved.");
   }
   MPI_Finalize();

    
   return(0);
}
