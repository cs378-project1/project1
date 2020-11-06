/* Copyright 2020 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: Name of author
                  Email of author
                                                                     */

#include "FLAME.h"

int Syr2k_unb_var5(FLA_Obj A, FLA_Obj B, FLA_Obj C)
{
    FLA_Obj AT, A0,
        AB, a1t,
        A2;

    FLA_Obj BT, B0,
        BB, b1t,
        B2;

    FLA_Obj CTL, CTR, C00, c01, C02,
        CBL, CBR, c10t, gamma11, c12t,
        C20, c21, C22;

    FLA_Part_2x1(A, &AT,
                 &AB, 0, FLA_BOTTOM);

    FLA_Part_2x1(B, &BT,
                 &BB, 0, FLA_BOTTOM);

    FLA_Part_2x2(C, &CTL, &CTR,
                 &CBL, &CBR, 0, 0, FLA_BR);

    while (FLA_Obj_length(CBR) < FLA_Obj_length(C))
    {

        // 5a

        FLA_Repart_2x1_to_3x1(AT, &A0,
                              &a1t,
                              /* ** */ /* *** */
                              AB, &A2, 1, FLA_TOP);

        FLA_Repart_2x1_to_3x1(BT, &B0,
                              &b1t,
                              /* ** */ /* *** */
                              BB, &B2, 1, FLA_TOP);

        FLA_Repart_2x2_to_3x3(CTL, /**/ CTR, &C00, /**/ &c01, &C02,
                              /* ************* */ /* ************************** */
                              &c10t, /**/ &gamma11, &c12t,
                              CBL, /**/ CBR, &C20, /**/ &c21, &C22,
                              1, 1, FLA_TL);

        /*------------------------------------------------------------*/

        /*                       update line 1                        */
        /*                             :                              */
        /*                       update line n                        */
        
        FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, a1t, b1t, FLA_ONE, gamma11);
        FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, b1t, a1t, FLA_ONE, gamma11);

        FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, a1t, B0, FLA_ONE, c10t);

        FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, B2, a1t, FLA_ONE, c21);
        
        /*------------------------------------------------------------*/

        // 5b

        FLA_Cont_with_3x1_to_2x1(&AT, A0,
                                 /* ** */ /* *** */
                                 a1t,
                                 &AB, A2, FLA_BOTTOM);

        FLA_Cont_with_3x1_to_2x1(&BT, B0,
                                 /* ** */ /* *** */
                                 b1t,
                                 &BB, B2, FLA_BOTTOM);

        FLA_Cont_with_3x3_to_2x2(&CTL, /**/ &CTR, C00, c01, /**/ C02,
                                 c10t, gamma11, /**/ c12t,
                                 /* ************** */ /* ************************ */
                                 &CBL, /**/ &CBR, C20, c21, /**/ C22,
                                 FLA_BR);
    }

    return FLA_SUCCESS;
}
