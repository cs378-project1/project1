#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

void trmm_llnn_unb_var1( FLA_Obj, FLA_Obj );
int Syr2k_unb_var2( FLA_Obj A, FLA_Obj B, FLA_Obj C );

int main(int argc, char *argv[])
{
  int n, nfirst, nlast, ninc, i, irep, nrepeats;

  double
    dtime, dtime_best, 
    diff;

  dtime_best = 0.0;
  
  FLA_Obj Aobj, Bobj, Cobj, Cold, Cref;

  /* Initialize FLAME. */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed. */
  printf( "%% enter nfirst, nlast, ninc:" );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  printf( "%% %d %d %d \n", nfirst, nlast, ninc );
  fflush( stdout );

  i = 1;
  for ( n=nfirst; n<= nlast; n+=ninc ){

    /* Allocate space for the matrices and vectors */
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Aobj );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Bobj );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Cobj );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Cold );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Cref );

    /* Generate random matrix A, and vectors x, and y */
    FLA_Random_matrix( Aobj );
    FLA_Random_matrix( Bobj );
    FLA_Random_symm_matrix( FLA_LOWER_TRIANGULAR, Cold );



    for ( irep=0; irep<nrepeats; irep++ ) {
    /* Time reference implementation (from libflame) */
      FLA_Copy( Cold, Cref );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute Bref = L Bref where L is lower triangular stored in the
	 lower triangular part of array L, by calling FLA_Trmm.  The
	 result ends up in Bref, which we will consider to be the
	 correct result. */
      FLA_Syr2k(FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, Aobj, Bobj, FLA_ONE, Cref);

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best );
    fflush( stdout );

    /* Time your unblocked Variant 1 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( Cold, Cobj );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      //      FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, Lobj, Bobj );
      Syr2k_unb_var2( Aobj, Bobj, Cobj );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( Cobj, Cref );
  
    printf( "data_unb_var1( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    fflush( stdout );

    FLA_Obj_free( &Aobj );
    FLA_Obj_free( &Bobj );
    FLA_Obj_free( &Cobj );
    FLA_Obj_free( &Cref );
    FLA_Obj_free( &Cold );

    i++;
  }
  FLA_Finalize( );

  exit( 0 );
}
