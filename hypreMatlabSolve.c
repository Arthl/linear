/*
 Read a linear system assembled in Matlab
 Solve the linear system using hypre
 Interface: Linear-Algebraic (IJ)
 Available solvers: 0  - AMG (default)
                    1  - AMG-PCG
                    8  - ParaSails-PCG
                    50 - PCG
                    61 - AMG-FlexGMRES

 Copyright: M. Giacomini (2017) / S. Zlotnik / A. Lustman (2019)
 (based on ex5.c in the hypre examples folder)
*/

#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

#include "csparse.h"

#include "vis.c"

#include "arthur.h"


cs *loadMatrix(char *str)
  {
  	 cs *T;
  	 cs *matrix;
  	FILE *in_file  = fopen(str, "r"); // read only 
  	if (in_file == NULL) 
  	{   
      printf("Error! Could not open file\n"); 
      exit(-1); // must include stdlib.h 
    } 
    T = cs_load ( in_file );	
  	matrix = cs_triplet ( T );
  	cs_spfree ( T );
  	return ( matrix ) ;
  }

int main (int argc, char *argv[])
{
   HYPRE_Int i;
   HYPRE_Int iErr = 0;
    
   HYPRE_Int myid, num_procs, dummy;
   HYPRE_Int solver_id;
    
   HYPRE_Int first_local_row, last_local_row;//, local_num_rows;
   HYPRE_Int first_local_col, last_local_col, local_num_cols;

  //  HYPRE_Int first_G_row, last_G_row;//, local_num_rows;
  //  HYPRE_Int first_G_col, last_G_col, local_G_cols; //local_G_cols used for initiating the pressure vector
    
   HYPRE_Real *values;

   HYPRE_IJMatrix ij_A;
   HYPRE_ParCSRMatrix parcsr_A;
  //  HYPRE_IJMatrix ij_G;
  //  HYPRE_ParCSRMatrix parcsr_G;
   HYPRE_IJVector ij_b;
   HYPRE_ParVector par_b;
  //  HYPRE_IJVector ij_h;
  //  HYPRE_ParVector par_h;
   HYPRE_IJVector ij_vec;
   HYPRE_ParVector par_vec;
   //HYPRE_IJVector ij_pres;
   //HYPRE_ParVector par_pres;
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
    if (saveMatsAs=='a')
    {
      iErr = HYPRE_IJMatrixRead( "matrixA", MPI_COMM_WORLD, HYPRE_PARCSR, &ij_A );
    }
    else
    {
      iErr = HYPRE_IJMatrixRead_binary( "matrixA", MPI_COMM_WORLD, HYPRE_PARCSR, &ij_A );
    }


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
  //     <filename>  = IJ.A.out to read in what has been printed out (processor numbers are omitted). */
  //   iErr = HYPRE_IJMatrixRead( "matrixA", MPI_COMM_WORLD, HYPRE_PARCSR, &ij_G );
  //   //iErr = HYPRE_IJMatrixRead( "matrixA", MPI_COMM_WORLD, HYPRE_PARCSR, &ij_G );



  //   if (iErr) 
  //   {
  //       hypre_printf("ERROR: Problem reading in the system matrix!\n");
  //       exit(1);
  //   }
  //   /* Get dimension info */
  //   iErr = HYPRE_IJMatrixGetLocalRange( ij_G,
  //                                      &first_G_row, &last_G_row ,
  //                                      &first_G_col, &last_G_col );
    
  //   //local_num_rows = last_local_row - first_local_row + 1;
  //   local_G_cols = last_G_col - first_G_col + 1;
  //  /* Get the parcsr matrix object to use */
  //  iErr += HYPRE_IJMatrixGetObject( ij_G, &object);
  //  parcsr_G = (HYPRE_ParCSRMatrix) object;
    
    
   /*  Read the RHS previously assembled */
   iErr = HYPRE_ParVectorRead(MPI_COMM_WORLD, "vectorF.0", &par_b);
   if (iErr)
   {
       hypre_printf("ERROR: Problem reading in the right-hand-side!\n");
       exit(1);
   }
   ij_b = NULL;


   /*  Read the RHS previously assembled */
  //  iErr = HYPRE_ParVectorRead(MPI_COMM_WORLD, "vectorH.0", &par_h);
  //  if (iErr)
  //  {
  //      hypre_printf("ERROR: Problem reading in the right-hand-side!\n");
  //      exit(1);
  //  }
  //  ij_h = NULL;
    
  /* Reading of the Matrix  */
  cs *Gmatrix;
  Gmatrix = loadMatrix("matrixG.00000");

  /* Get the transpose matrix G*/
  cs *GTmatrix;
  GTmatrix = cs_transpose(Gmatrix,1);

  /* Load for both Vectors*/
  int sizeVectorF;
  sizeVectorF = sizeVectorFile("vectorF");
  double fVector[sizeVectorF];
  loadVector("vectorF", fVector, sizeVectorF);

  int sizeVectorH;
  sizeVectorH = sizeVectorFile("vectorH");
  double hVector[sizeVectorH];
  loadVector("vectorH", hVector, sizeVectorH);


   /* Create the initial solution and set it to zero */
   HYPRE_IJVectorCreate(MPI_COMM_WORLD, first_local_col, last_local_col, &ij_vec);
   HYPRE_IJVectorSetObjectType(ij_vec, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(ij_vec);
   /* Initialize the guess vector */
   values = hypre_CTAlloc(HYPRE_Real, local_num_cols);
   for (i = 0; i < local_num_cols; i++)
      values[i] = 0.;
   HYPRE_IJVectorSetValues(ij_vec, local_num_cols, NULL, values);
   hypre_TFree(values);
   /* Get the parcsr vector object to use */
   iErr = HYPRE_IJVectorGetObject( ij_vec, &object );
   par_vec = (HYPRE_ParVector) object;


  // Create the assigned variables for the resolution
  double bi_prod, alpha, beta, gamma, gamma_old;
  double d_1vector[sizeVectorF];//, uSol[sizeVectorF]; 
  double d_2vector[sizeVectorH], pSol[sizeVectorH], residual[sizeVectorH];
  setArrayZeros(d_1vector, sizeVectorF);
  // setArrayZeros(uSol, sizeVectorF);
  setArrayZeros(d_2vector, sizeVectorH); //not necessary since it is copied late
  setArrayZeros(pSol, sizeVectorH);
  setArrayZeros(residual, sizeVectorH);

  // Usage of temporary arrays and values, never declared in the Uzawa algorithm
  double temporary[sizeVectorF]; // array used for the Hypre communication
  double transFirst[sizeVectorF], transSecond[sizeVectorH]; // array used for Gmatrix*d_2vector
  double abba; // used for temporary values to assign
  setArrayZeros(temporary, sizeVectorF); // technically not needed because done in the loop
  setArrayZeros(transFirst, sizeVectorF);
  setArrayZeros(transSecond, sizeVectorH);

  /* Variables used for the loop */
  int it, max_it;
  it = 0;
  max_it = 50;
  double tol_r;
  tol_r = 1e-8;

  // CheckUp of the dimensions
  // // // // // // HOW TO PRINT THE DIMENSIONS OF Gmatrix ??????
  printf("Size of F array is %ld\n", sizeof(fVector)/sizeof(double));
  printf("Size of H array is %ld\n", sizeof(hVector)/sizeof(double));

  /* Before the loop*/
  // // // // // // SOLVE PROBLEM Kmatrix*uSol = fVector;
  int num_iterations;
  double final_res_norm;

  HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

  /* Set some parameters (See Reference Manual for more parameters) */
  HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
  HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
  HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
  HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
  HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

  HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_vec);
  HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_vec);

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

  /* Prints the solution*/
  HYPRE_ParVectorPrint(par_vec, "initSolution.0");
  
  // Initialize the Vector uSol from the solution computed
  int sizeVectorFverif;
  sizeVectorFverif = sizeVectorFile("initSolution.0.0");
  double uSol[sizeVectorF];
  load2SplitVector("initSolution.0.0", "initSolution.0.1" , uSol, sizeVectorFverif*2);

 // // EVERYTHING IS FINE TILL HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 // command r = GT*u + r */ int cs_gaxpy (const cs *A, const double *x, double *y)
 // never put the same array in there ! always pass by another array
 // CAUTION : if a dimension error is made it is not shown !!!
 bi_prod = dotprod(fVector,fVector,sizeVectorF) + dotprod(hVector,hVector,sizeVectorH);
 printf("Bi prod of RHS is %f\n", bi_prod); // NEVER NEEDED
 cs_gaxpy(GTmatrix, uSol, residual);
 copyVector(residual, d_2vector, sizeVectorH);
 gamma = dotprod(residual, residual, sizeVectorH);
 
 printf("Initial residual %f\n", pow(gamma,0.5));

 // Loop Starts Here
 for (it = 0; it < 1; it++) 
  {
    if (myid==0)
    {

      setArrayZeros(temporary, sizeVectorF);
      cs_gaxpy(Gmatrix,d_2vector, temporary);
      writingSplitVector(temporary, sizeVectorF, num_procs);
    }
    // MPI_Barrier(MPI_COMM_WORLD);

    // To read vector for solver!!!!!!!!
    iErr = HYPRE_ParVectorRead(MPI_COMM_WORLD, "myfile.0", &par_b);
    if (iErr)
    {
        hypre_printf("ERROR: Problem reading in the splitted Vector!\n");
        exit(1);
    }
    ij_b = NULL;


    /* PCG */
    if (solver_id == 50)
      {

        /* Create solver */
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
        HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
        HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
        HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
        HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

        /* Now setup and solve! */
        HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_vec);
        HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_vec);

        //HYPRE_PCGGetNumIterations(solver, &num_iterations);
        //HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (myid == 0)
        {
            // printf("\n");
            //printf("Iterations = %d\n", num_iterations);
            // printf("Final Relative Residual Norm = %e\n", final_res_norm);
            // printf("\n");
        }

        /* Destroy solver */
        HYPRE_ParCSRPCGDestroy(solver);

      }
        
        /* Save the solution to file */
        HYPRE_ParVectorPrint(par_vec, "solution.0");

    if (myid==0)
      {
        sizeVectorFverif = sizeVectorFile("solution.0.0");
        load2SplitVector("solution.0.0", "solution.0.1" , d_1vector, sizeVectorFverif*2);
	printf("Hey, vector dot product of d_1 is %.5e\n", dotprod(d_1vector,d_1vector,sizeVectorF));

        // transFirst = Gmatrix * d_2vector;
        setArrayZeros(transFirst, sizeVectorF);
        cs_gaxpy(Gmatrix, d_2vector, transFirst);

        abba = dotprod(transFirst, d_1vector, sizeVectorF);
        alpha = gamma / abba;
        gamma_old = gamma;

        // pSol = pSol + alpha * d_2vector;
        incrementOp(pSol, d_2vector, alpha, sizeVectorH);
        //uSol = uSol - alpha * d_1vector;
        incrementOp(uSol, d_1vector, -alpha, sizeVectorF);

        // transSecond = Gmatrix.transpose() * d_1vector;
        setArrayZeros(transSecond, sizeVectorH);
        cs_gaxpy(GTmatrix, d_1vector, transSecond);

        // residual = residual - alpha * transSecond;
        incrementOp(residual, transSecond, -alpha, sizeVectorH);
        
        abba = dotprod(residual, residual, sizeVectorH);
        abba = pow(abba,0.5);
        if (abba<tol_r)
          {
            printf("You got it fam\n");
            printf("Congrats on the convergence, %e\n", abba);
            printf("Congrats on the convergence, %e\n", tol_r);
	    writingSplitVector(uSol, sizeVectorF, 1);
	    //writingSplitVector(pSol, sizeVectorH, 1);
            printf("Dot prod of velocity, %e\n", dotprod(uSol,uSol, sizeVectorF));
            printf("Dot prod of the pressure, %e\n", dotprod(pSol, pSol, sizeVectorH));
            return( -1 );
          }

        // it++;
      
        abba = dotprod(residual,residual,sizeVectorH);
        printf("Iteration num %d, ", it);
        printf("norm of residual is %f\n", pow(abba,0.5));

        gamma = dotprod(residual, residual, sizeVectorH);
        beta = gamma / gamma_old;

        // d_2vector = transSecond + beta * d_2vector;
        D2incrementOp(d_2vector, transSecond, beta, sizeVectorH);

	printf("Norm of the residual %e\n", pow(dotprod(residual,residual,sizeVectorH),0.5));
	printf("Norm of the d_2 direction %e\n", pow(dotprod(d_2vector,d_2vector,sizeVectorH),0.5));
	printf("Dot prod of velocity, %e\n", dotprod(uSol,uSol, sizeVectorF));
	printf("Dot prod of the pressure, %e\n", dotprod(pSol, pSol, sizeVectorH));
      }

 }

	printf("Norm of the residual %e\n", pow(dotprod(residual,residual,sizeVectorH),0.5));
	printf("Norm of the d_2 direction %e\n", dotprod(d_2vector,d_2vector,sizeVectorH));
	printf("Dot prod of velocity, %e\n", dotprod(uSol,uSol, sizeVectorF));
	printf("Dot prod of the pressure, %e\n", dotprod(pSol, pSol, sizeVectorH));

   /* Clean up */
   HYPRE_IJMatrixDestroy(ij_A);
  //  HYPRE_IJMatrixDestroy(ij_G);
  //  HYPRE_IJVectorDestroy(ij_f);
  //  HYPRE_IJVectorDestroy(ij_h);
  //  HYPRE_IJVectorDestroy(ij_vel);
   // HYPRE_IJVectorDestroy(ij_pres);
   HYPRE_IJVectorDestroy(ij_b);
   HYPRE_IJVectorDestroy(ij_vec);
   

   /* Finalize MPI*/
   if (myid == 0)
   {
       printf("\n Linear system correctly solved.");
   }
   MPI_Finalize();

    
   return(0);
}
