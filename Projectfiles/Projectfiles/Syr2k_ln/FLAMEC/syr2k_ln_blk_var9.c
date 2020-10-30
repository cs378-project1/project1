/* Copyright 2020 The University of Texas at Austin
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html

   Programmed by: Daimu Iwata
                  di2937@utexas.edu
                                                                     */

#include "FLAME.h"

int syr2k_ln_blk_var9( FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  int b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    b = min( FLA_Obj_width( AR ), nb_alg );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    /* C := A1 * B1' + B1 * A1' + C */
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A1, B1, FLA_ONE, C );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, B1, A1, FLA_ONE, C );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
