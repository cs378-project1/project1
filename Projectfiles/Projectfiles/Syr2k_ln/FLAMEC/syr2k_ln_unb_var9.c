/* Copyright 2020 The University of Texas at Austin
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html

   Programmed by: Daimu Iwata
                  di2937@utexas.edu
                                                                     */

#include "FLAME.h"

int syr2k_ln_unb_var9( FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &b1, &B2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    /* C := a1 * b1' + b1 * a1' + C */
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, a1, b1, FLA_ONE, C );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, b1, a1, FLA_ONE, C );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, b1, /**/ B2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
