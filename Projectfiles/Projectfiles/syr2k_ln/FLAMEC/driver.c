#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0
#define BLOCK_SZ 2

char *strings_unb[6] = {
    "data_unb_var1( %d, 1:3 ) = [ %d %le %le];\n",
    "data_unb_var2( %d, 1:3 ) = [ %d %le %le];\n",
    "data_unb_var3( %d, 1:3 ) = [ %d %le %le];\n",
    "data_unb_var4( %d, 1:3 ) = [ %d %le %le];\n",
    "data_unb_var5( %d, 1:3 ) = [ %d %le %le];\n",
    "data_unb_var9( %d, 1:3 ) = [ %d %le %le];\n"};

char* strings_blk[6] = {
    "data_blk_var1( %d, 1:3 ) = [ %d %le %le];\n",
    "data_blk_var2( %d, 1:3 ) = [ %d %le %le];\n",
    "data_blk_var3( %d, 1:3 ) = [ %d %le %le];\n",
    "data_blk_var4( %d, 1:3 ) = [ %d %le %le];\n",
    "data_blk_var5( %d, 1:3 ) = [ %d %le %le];\n",
    "data_blk_var9( %d, 1:3 ) = [ %d %le %le];\n"};

int syr2k_ln_unb_var1(FLA_Obj A, FLA_Obj B, FLA_Obj C);
int syr2k_ln_unb_var2(FLA_Obj A, FLA_Obj B, FLA_Obj C);
int syr2k_ln_unb_var3(FLA_Obj A, FLA_Obj B, FLA_Obj C);
int syr2k_ln_unb_var4(FLA_Obj A, FLA_Obj B, FLA_Obj C);
int syr2k_ln_unb_var5(FLA_Obj A, FLA_Obj B, FLA_Obj C);
int syr2k_ln_unb_var9(FLA_Obj A, FLA_Obj B, FLA_Obj C);
int syr2k_ln_blk_var1(FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg);
int syr2k_ln_blk_var2(FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg);
int syr2k_ln_blk_var3(FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg);
int syr2k_ln_blk_var4(FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg);
int syr2k_ln_blk_var5(FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg);
int syr2k_ln_blk_var9(FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg);

int (*functions_blk[6])(FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg) = {
    syr2k_ln_blk_var1,
    syr2k_ln_blk_var2,
    syr2k_ln_blk_var3,
    syr2k_ln_blk_var4,
    syr2k_ln_blk_var5,
    syr2k_ln_blk_var9};

int (*functions_unb[6])(FLA_Obj A, FLA_Obj B, FLA_Obj C) = {
    syr2k_ln_unb_var1,
    syr2k_ln_unb_var2,
    syr2k_ln_unb_var3,
    syr2k_ln_unb_var4,
    syr2k_ln_unb_var5,
    syr2k_ln_unb_var9};


int main(int argc, char *argv[])
{
  int n, nfirst, nlast, ninc, i, irep, nrepeats;

  double
      dtime,
      dtime_best,
      diff;

  FLA_Obj Aobj, Bobj, Cobj, Aold, Bold, Cold, Cref;

  /* Initialize FLAME. */
  FLA_Init();

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf("%% number of repeats:");
  scanf("%d", &nrepeats);
  printf("%% %d\n", nrepeats);

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed. */
  printf("%% enter nfirst, nlast, ninc:");
  scanf("%d%d%d", &nfirst, &nlast, &ninc);
  printf("%% %d %d %d \n", nfirst, nlast, ninc);
  fflush(stdout);

  i = 1;
  for (n = nfirst; n <= nlast; n += ninc)
  {

    /* Allocate space for the matrices and vectors */
    FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Aobj);
    FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Bobj);
    FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Cobj);
    FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Aold);
    FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Bold);
    FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Cold);
    FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Cref);

    /* Generate random matrix A, and vectors x, and y */
    FLA_Random_matrix(Aold);
    FLA_Random_matrix(Bold);
    FLA_Random_symm_matrix(FLA_LOWER_TRIANGULAR, Cold);

    for (irep = 0; irep < nrepeats; irep++)
    {
      /* Time reference implementation (from libflame) */
      FLA_Copy(Cold, Cref);
      FLA_Copy(Aold, Aobj);
      FLA_Copy(Bold, Bobj);
      /* start clock */
      dtime = FLA_Clock();

      /* Compute Cref = AB^T + BA^T + Cref where Cref is lower triangular stored in the
	 lower triangular part of array Cref, by calling FLA_Syr2k.  The
	 result ends up in Cref, which we will consider to be the
	 correct result. */
      FLA_Syr2k(FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, Aobj, Bobj, FLA_ONE, Cref);

      /* stop clock */
      dtime = FLA_Clock() - dtime;

      if (irep == 0)
        dtime_best = dtime;
      else
        dtime_best = (dtime < dtime_best ? dtime : dtime_best);
    }

    printf("data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best);
    fflush(stdout);

    /* Time your unblocked Variant 1 */
    for (int variant = 0; variant < 6; variant++) 
    {
      // blocked
      for (irep = 0; irep < nrepeats; irep++)
      {
        /* Copy vector yold to y */
        FLA_Copy(Cold, Cobj);
        FLA_Copy(Aold, Aobj);
        FLA_Copy(Bold, Bobj);
        
        /* start clock */
        dtime = FLA_Clock();

        /* Comment out the below call and call your routine instead */
        functions_blk[variant](Aobj, Bobj, Cobj, BLOCK_SZ);
        //syr2k_ln_blk_var1(Aobj, Bobj, Cobj, BLOCK_SZ);

        /* stop clock */
        dtime = FLA_Clock() - dtime;

        if (irep == 0)
          dtime_best = dtime;
        else
          dtime_best = (dtime < dtime_best ? dtime : dtime_best);
      }

      diff = FLA_Max_elemwise_diff(Cobj, Cref);

      printf(strings_blk[variant], i, n,
             dtime_best, diff);

      fflush(stdout);

      // unblocked
      for (irep = 0; irep < nrepeats; irep++)
      {
        /* Copy vector yold to y */
        FLA_Copy(Cold, Cobj);
        FLA_Copy(Aold, Aobj);
        FLA_Copy(Bold, Bobj);

        /* start clock */
        dtime = FLA_Clock();

        /* Comment out the below call and call your routine instead */
        functions_unb[variant](Aobj, Bobj, Cobj);
        //syr2k_ln_unb_var1(Aobj, Bobj, Cobj);

        /* stop clock */
        dtime = FLA_Clock() - dtime;

        if (irep == 0)
          dtime_best = dtime;
        else
          dtime_best = (dtime < dtime_best ? dtime : dtime_best);
      }

      diff = FLA_Max_elemwise_diff(Cobj, Cref);

      printf(strings_unb[variant], i, n,
             dtime_best, diff);

      fflush(stdout);
    }
    FLA_Obj_free(&Aobj);
    FLA_Obj_free(&Bobj);
    FLA_Obj_free(&Cobj);
    FLA_Obj_free(&Cref);
    FLA_Obj_free(&Aold);
    FLA_Obj_free(&Bold);
    FLA_Obj_free(&Cold);

    i++;
  }
  FLA_Finalize();

  exit(0);
}
