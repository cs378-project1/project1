/* Copyright 2020 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: Name of author
                  Email of author
                                                                     */

#include "FLAME.h"

int Syr2k_blk_var5( FLA_Obj C, FLA_Obj A, FLA_Obj B, int nb_alg )
{
  FLA_Obj CTL,   CTR,      C00, C01, C02, 
          CBL,   CBR,      C10, C11, C12,
                           C20, C21, C22;

  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  int b;

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_BR );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( CBR ) < FLA_Obj_length( C ) ){

    b = min( FLA_Obj_length( CTL ), nb_alg );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00, &C01, /**/ &C02,
                                                &C10, &C11, /**/ &C12,
                        /* ************* */   /* ******************** */
                           CBL, /**/ CBR,       &C20, &C21, /**/ &C22,
                           b, b, FLA_TL );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &A1, 
                        /* ** */            /* ** */
                           AB,                &A2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &B1, 
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /*                       update line 1                        */
    /*                             :                              */
    /*                       update line n                        */

    FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A1, B1, FLA_ONE, C11);
    FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, B1, A1, FLA_ONE, C11);

    FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A1, B0, FLA_ONE, C10);
    FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, B2, A1, FLA_ONE, C21);
        

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00, /**/ C01, C02,
                            /* ************** */  /* ****************** */
                                                     C10, /**/ C11, C12,
                              &CBL, /**/ &CBR,       C20, /**/ C21, C22,
                              FLA_BR );

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* ** */
                                                  A1, 
                              &AB,                A2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* ** */
                                                  B1, 
                              &BB,                B2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}