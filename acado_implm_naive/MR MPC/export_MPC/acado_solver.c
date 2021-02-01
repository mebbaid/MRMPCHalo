/*
 *    This file was auto-generated using the ACADO Toolkit.
 *    
 *    While ACADO Toolkit is free software released under the terms of
 *    the GNU Lesser General Public License (LGPL), the generated code
 *    as such remains the property of the user who used ACADO Toolkit
 *    to generate this code. In particular, user dependent data of the code
 *    do not inherit the GNU LGPL license. On the other hand, parts of the
 *    generated code that are a direct copy of source code from the
 *    ACADO Toolkit or the software tools it is based on, remain, as derived
 *    work, automatically covered by the LGPL license.
 *    
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *    
 */


#include "acado_common.h"




/******************************************************************************/
/*                                                                            */
/* ACADO code generation                                                      */
/*                                                                            */
/******************************************************************************/


int acado_modelSimulation(  )
{
int ret;

int lRun1;
ret = 0;
for (lRun1 = 0; lRun1 < 10; ++lRun1)
{
acadoWorkspace.state[0] = acadoVariables.x[lRun1 * 6];
acadoWorkspace.state[1] = acadoVariables.x[lRun1 * 6 + 1];
acadoWorkspace.state[2] = acadoVariables.x[lRun1 * 6 + 2];
acadoWorkspace.state[3] = acadoVariables.x[lRun1 * 6 + 3];
acadoWorkspace.state[4] = acadoVariables.x[lRun1 * 6 + 4];
acadoWorkspace.state[5] = acadoVariables.x[lRun1 * 6 + 5];

acadoWorkspace.state[60] = acadoVariables.u[lRun1 * 3];
acadoWorkspace.state[61] = acadoVariables.u[lRun1 * 3 + 1];
acadoWorkspace.state[62] = acadoVariables.u[lRun1 * 3 + 2];
acadoWorkspace.state[63] = acadoVariables.od[lRun1 * 4];
acadoWorkspace.state[64] = acadoVariables.od[lRun1 * 4 + 1];
acadoWorkspace.state[65] = acadoVariables.od[lRun1 * 4 + 2];
acadoWorkspace.state[66] = acadoVariables.od[lRun1 * 4 + 3];

ret = acado_integrate(acadoWorkspace.state, 1);

acadoWorkspace.d[lRun1 * 6] = acadoWorkspace.state[0] - acadoVariables.x[lRun1 * 6 + 6];
acadoWorkspace.d[lRun1 * 6 + 1] = acadoWorkspace.state[1] - acadoVariables.x[lRun1 * 6 + 7];
acadoWorkspace.d[lRun1 * 6 + 2] = acadoWorkspace.state[2] - acadoVariables.x[lRun1 * 6 + 8];
acadoWorkspace.d[lRun1 * 6 + 3] = acadoWorkspace.state[3] - acadoVariables.x[lRun1 * 6 + 9];
acadoWorkspace.d[lRun1 * 6 + 4] = acadoWorkspace.state[4] - acadoVariables.x[lRun1 * 6 + 10];
acadoWorkspace.d[lRun1 * 6 + 5] = acadoWorkspace.state[5] - acadoVariables.x[lRun1 * 6 + 11];

acadoWorkspace.evGx[lRun1 * 36] = acadoWorkspace.state[6];
acadoWorkspace.evGx[lRun1 * 36 + 1] = acadoWorkspace.state[7];
acadoWorkspace.evGx[lRun1 * 36 + 2] = acadoWorkspace.state[8];
acadoWorkspace.evGx[lRun1 * 36 + 3] = acadoWorkspace.state[9];
acadoWorkspace.evGx[lRun1 * 36 + 4] = acadoWorkspace.state[10];
acadoWorkspace.evGx[lRun1 * 36 + 5] = acadoWorkspace.state[11];
acadoWorkspace.evGx[lRun1 * 36 + 6] = acadoWorkspace.state[12];
acadoWorkspace.evGx[lRun1 * 36 + 7] = acadoWorkspace.state[13];
acadoWorkspace.evGx[lRun1 * 36 + 8] = acadoWorkspace.state[14];
acadoWorkspace.evGx[lRun1 * 36 + 9] = acadoWorkspace.state[15];
acadoWorkspace.evGx[lRun1 * 36 + 10] = acadoWorkspace.state[16];
acadoWorkspace.evGx[lRun1 * 36 + 11] = acadoWorkspace.state[17];
acadoWorkspace.evGx[lRun1 * 36 + 12] = acadoWorkspace.state[18];
acadoWorkspace.evGx[lRun1 * 36 + 13] = acadoWorkspace.state[19];
acadoWorkspace.evGx[lRun1 * 36 + 14] = acadoWorkspace.state[20];
acadoWorkspace.evGx[lRun1 * 36 + 15] = acadoWorkspace.state[21];
acadoWorkspace.evGx[lRun1 * 36 + 16] = acadoWorkspace.state[22];
acadoWorkspace.evGx[lRun1 * 36 + 17] = acadoWorkspace.state[23];
acadoWorkspace.evGx[lRun1 * 36 + 18] = acadoWorkspace.state[24];
acadoWorkspace.evGx[lRun1 * 36 + 19] = acadoWorkspace.state[25];
acadoWorkspace.evGx[lRun1 * 36 + 20] = acadoWorkspace.state[26];
acadoWorkspace.evGx[lRun1 * 36 + 21] = acadoWorkspace.state[27];
acadoWorkspace.evGx[lRun1 * 36 + 22] = acadoWorkspace.state[28];
acadoWorkspace.evGx[lRun1 * 36 + 23] = acadoWorkspace.state[29];
acadoWorkspace.evGx[lRun1 * 36 + 24] = acadoWorkspace.state[30];
acadoWorkspace.evGx[lRun1 * 36 + 25] = acadoWorkspace.state[31];
acadoWorkspace.evGx[lRun1 * 36 + 26] = acadoWorkspace.state[32];
acadoWorkspace.evGx[lRun1 * 36 + 27] = acadoWorkspace.state[33];
acadoWorkspace.evGx[lRun1 * 36 + 28] = acadoWorkspace.state[34];
acadoWorkspace.evGx[lRun1 * 36 + 29] = acadoWorkspace.state[35];
acadoWorkspace.evGx[lRun1 * 36 + 30] = acadoWorkspace.state[36];
acadoWorkspace.evGx[lRun1 * 36 + 31] = acadoWorkspace.state[37];
acadoWorkspace.evGx[lRun1 * 36 + 32] = acadoWorkspace.state[38];
acadoWorkspace.evGx[lRun1 * 36 + 33] = acadoWorkspace.state[39];
acadoWorkspace.evGx[lRun1 * 36 + 34] = acadoWorkspace.state[40];
acadoWorkspace.evGx[lRun1 * 36 + 35] = acadoWorkspace.state[41];

acadoWorkspace.evGu[lRun1 * 18] = acadoWorkspace.state[42];
acadoWorkspace.evGu[lRun1 * 18 + 1] = acadoWorkspace.state[43];
acadoWorkspace.evGu[lRun1 * 18 + 2] = acadoWorkspace.state[44];
acadoWorkspace.evGu[lRun1 * 18 + 3] = acadoWorkspace.state[45];
acadoWorkspace.evGu[lRun1 * 18 + 4] = acadoWorkspace.state[46];
acadoWorkspace.evGu[lRun1 * 18 + 5] = acadoWorkspace.state[47];
acadoWorkspace.evGu[lRun1 * 18 + 6] = acadoWorkspace.state[48];
acadoWorkspace.evGu[lRun1 * 18 + 7] = acadoWorkspace.state[49];
acadoWorkspace.evGu[lRun1 * 18 + 8] = acadoWorkspace.state[50];
acadoWorkspace.evGu[lRun1 * 18 + 9] = acadoWorkspace.state[51];
acadoWorkspace.evGu[lRun1 * 18 + 10] = acadoWorkspace.state[52];
acadoWorkspace.evGu[lRun1 * 18 + 11] = acadoWorkspace.state[53];
acadoWorkspace.evGu[lRun1 * 18 + 12] = acadoWorkspace.state[54];
acadoWorkspace.evGu[lRun1 * 18 + 13] = acadoWorkspace.state[55];
acadoWorkspace.evGu[lRun1 * 18 + 14] = acadoWorkspace.state[56];
acadoWorkspace.evGu[lRun1 * 18 + 15] = acadoWorkspace.state[57];
acadoWorkspace.evGu[lRun1 * 18 + 16] = acadoWorkspace.state[58];
acadoWorkspace.evGu[lRun1 * 18 + 17] = acadoWorkspace.state[59];
}
return ret;
}

void acado_evaluateLSQ(const real_t* in, real_t* out)
{
const real_t* xd = in;
const real_t* u = in + 6;

/* Compute outputs: */
out[0] = xd[0];
out[1] = xd[1];
out[2] = xd[2];
out[3] = xd[3];
out[4] = xd[4];
out[5] = xd[5];
out[6] = u[0];
out[7] = u[1];
out[8] = u[2];
}

void acado_evaluateLSQEndTerm(const real_t* in, real_t* out)
{
const real_t* xd = in;

/* Compute outputs: */
out[0] = xd[0];
out[1] = xd[1];
out[2] = xd[2];
out[3] = xd[3];
out[4] = xd[4];
out[5] = xd[5];
}

void acado_setObjQ1Q2( real_t* const tmpObjS, real_t* const tmpQ1, real_t* const tmpQ2 )
{
tmpQ2[0] = +tmpObjS[0];
tmpQ2[1] = +tmpObjS[1];
tmpQ2[2] = +tmpObjS[2];
tmpQ2[3] = +tmpObjS[3];
tmpQ2[4] = +tmpObjS[4];
tmpQ2[5] = +tmpObjS[5];
tmpQ2[6] = +tmpObjS[6];
tmpQ2[7] = +tmpObjS[7];
tmpQ2[8] = +tmpObjS[8];
tmpQ2[9] = +tmpObjS[9];
tmpQ2[10] = +tmpObjS[10];
tmpQ2[11] = +tmpObjS[11];
tmpQ2[12] = +tmpObjS[12];
tmpQ2[13] = +tmpObjS[13];
tmpQ2[14] = +tmpObjS[14];
tmpQ2[15] = +tmpObjS[15];
tmpQ2[16] = +tmpObjS[16];
tmpQ2[17] = +tmpObjS[17];
tmpQ2[18] = +tmpObjS[18];
tmpQ2[19] = +tmpObjS[19];
tmpQ2[20] = +tmpObjS[20];
tmpQ2[21] = +tmpObjS[21];
tmpQ2[22] = +tmpObjS[22];
tmpQ2[23] = +tmpObjS[23];
tmpQ2[24] = +tmpObjS[24];
tmpQ2[25] = +tmpObjS[25];
tmpQ2[26] = +tmpObjS[26];
tmpQ2[27] = +tmpObjS[27];
tmpQ2[28] = +tmpObjS[28];
tmpQ2[29] = +tmpObjS[29];
tmpQ2[30] = +tmpObjS[30];
tmpQ2[31] = +tmpObjS[31];
tmpQ2[32] = +tmpObjS[32];
tmpQ2[33] = +tmpObjS[33];
tmpQ2[34] = +tmpObjS[34];
tmpQ2[35] = +tmpObjS[35];
tmpQ2[36] = +tmpObjS[36];
tmpQ2[37] = +tmpObjS[37];
tmpQ2[38] = +tmpObjS[38];
tmpQ2[39] = +tmpObjS[39];
tmpQ2[40] = +tmpObjS[40];
tmpQ2[41] = +tmpObjS[41];
tmpQ2[42] = +tmpObjS[42];
tmpQ2[43] = +tmpObjS[43];
tmpQ2[44] = +tmpObjS[44];
tmpQ2[45] = +tmpObjS[45];
tmpQ2[46] = +tmpObjS[46];
tmpQ2[47] = +tmpObjS[47];
tmpQ2[48] = +tmpObjS[48];
tmpQ2[49] = +tmpObjS[49];
tmpQ2[50] = +tmpObjS[50];
tmpQ2[51] = +tmpObjS[51];
tmpQ2[52] = +tmpObjS[52];
tmpQ2[53] = +tmpObjS[53];
tmpQ1[0] = + tmpQ2[0];
tmpQ1[1] = + tmpQ2[1];
tmpQ1[2] = + tmpQ2[2];
tmpQ1[3] = + tmpQ2[3];
tmpQ1[4] = + tmpQ2[4];
tmpQ1[5] = + tmpQ2[5];
tmpQ1[6] = + tmpQ2[9];
tmpQ1[7] = + tmpQ2[10];
tmpQ1[8] = + tmpQ2[11];
tmpQ1[9] = + tmpQ2[12];
tmpQ1[10] = + tmpQ2[13];
tmpQ1[11] = + tmpQ2[14];
tmpQ1[12] = + tmpQ2[18];
tmpQ1[13] = + tmpQ2[19];
tmpQ1[14] = + tmpQ2[20];
tmpQ1[15] = + tmpQ2[21];
tmpQ1[16] = + tmpQ2[22];
tmpQ1[17] = + tmpQ2[23];
tmpQ1[18] = + tmpQ2[27];
tmpQ1[19] = + tmpQ2[28];
tmpQ1[20] = + tmpQ2[29];
tmpQ1[21] = + tmpQ2[30];
tmpQ1[22] = + tmpQ2[31];
tmpQ1[23] = + tmpQ2[32];
tmpQ1[24] = + tmpQ2[36];
tmpQ1[25] = + tmpQ2[37];
tmpQ1[26] = + tmpQ2[38];
tmpQ1[27] = + tmpQ2[39];
tmpQ1[28] = + tmpQ2[40];
tmpQ1[29] = + tmpQ2[41];
tmpQ1[30] = + tmpQ2[45];
tmpQ1[31] = + tmpQ2[46];
tmpQ1[32] = + tmpQ2[47];
tmpQ1[33] = + tmpQ2[48];
tmpQ1[34] = + tmpQ2[49];
tmpQ1[35] = + tmpQ2[50];
}

void acado_setObjR1R2( real_t* const tmpObjS, real_t* const tmpR1, real_t* const tmpR2 )
{
tmpR2[0] = +tmpObjS[54];
tmpR2[1] = +tmpObjS[55];
tmpR2[2] = +tmpObjS[56];
tmpR2[3] = +tmpObjS[57];
tmpR2[4] = +tmpObjS[58];
tmpR2[5] = +tmpObjS[59];
tmpR2[6] = +tmpObjS[60];
tmpR2[7] = +tmpObjS[61];
tmpR2[8] = +tmpObjS[62];
tmpR2[9] = +tmpObjS[63];
tmpR2[10] = +tmpObjS[64];
tmpR2[11] = +tmpObjS[65];
tmpR2[12] = +tmpObjS[66];
tmpR2[13] = +tmpObjS[67];
tmpR2[14] = +tmpObjS[68];
tmpR2[15] = +tmpObjS[69];
tmpR2[16] = +tmpObjS[70];
tmpR2[17] = +tmpObjS[71];
tmpR2[18] = +tmpObjS[72];
tmpR2[19] = +tmpObjS[73];
tmpR2[20] = +tmpObjS[74];
tmpR2[21] = +tmpObjS[75];
tmpR2[22] = +tmpObjS[76];
tmpR2[23] = +tmpObjS[77];
tmpR2[24] = +tmpObjS[78];
tmpR2[25] = +tmpObjS[79];
tmpR2[26] = +tmpObjS[80];
tmpR1[0] = + tmpR2[6];
tmpR1[1] = + tmpR2[7];
tmpR1[2] = + tmpR2[8];
tmpR1[3] = + tmpR2[15];
tmpR1[4] = + tmpR2[16];
tmpR1[5] = + tmpR2[17];
tmpR1[6] = + tmpR2[24];
tmpR1[7] = + tmpR2[25];
tmpR1[8] = + tmpR2[26];
}

void acado_setObjQN1QN2( real_t* const tmpObjSEndTerm, real_t* const tmpQN1, real_t* const tmpQN2 )
{
tmpQN2[0] = +tmpObjSEndTerm[0];
tmpQN2[1] = +tmpObjSEndTerm[1];
tmpQN2[2] = +tmpObjSEndTerm[2];
tmpQN2[3] = +tmpObjSEndTerm[3];
tmpQN2[4] = +tmpObjSEndTerm[4];
tmpQN2[5] = +tmpObjSEndTerm[5];
tmpQN2[6] = +tmpObjSEndTerm[6];
tmpQN2[7] = +tmpObjSEndTerm[7];
tmpQN2[8] = +tmpObjSEndTerm[8];
tmpQN2[9] = +tmpObjSEndTerm[9];
tmpQN2[10] = +tmpObjSEndTerm[10];
tmpQN2[11] = +tmpObjSEndTerm[11];
tmpQN2[12] = +tmpObjSEndTerm[12];
tmpQN2[13] = +tmpObjSEndTerm[13];
tmpQN2[14] = +tmpObjSEndTerm[14];
tmpQN2[15] = +tmpObjSEndTerm[15];
tmpQN2[16] = +tmpObjSEndTerm[16];
tmpQN2[17] = +tmpObjSEndTerm[17];
tmpQN2[18] = +tmpObjSEndTerm[18];
tmpQN2[19] = +tmpObjSEndTerm[19];
tmpQN2[20] = +tmpObjSEndTerm[20];
tmpQN2[21] = +tmpObjSEndTerm[21];
tmpQN2[22] = +tmpObjSEndTerm[22];
tmpQN2[23] = +tmpObjSEndTerm[23];
tmpQN2[24] = +tmpObjSEndTerm[24];
tmpQN2[25] = +tmpObjSEndTerm[25];
tmpQN2[26] = +tmpObjSEndTerm[26];
tmpQN2[27] = +tmpObjSEndTerm[27];
tmpQN2[28] = +tmpObjSEndTerm[28];
tmpQN2[29] = +tmpObjSEndTerm[29];
tmpQN2[30] = +tmpObjSEndTerm[30];
tmpQN2[31] = +tmpObjSEndTerm[31];
tmpQN2[32] = +tmpObjSEndTerm[32];
tmpQN2[33] = +tmpObjSEndTerm[33];
tmpQN2[34] = +tmpObjSEndTerm[34];
tmpQN2[35] = +tmpObjSEndTerm[35];
tmpQN1[0] = + tmpQN2[0];
tmpQN1[1] = + tmpQN2[1];
tmpQN1[2] = + tmpQN2[2];
tmpQN1[3] = + tmpQN2[3];
tmpQN1[4] = + tmpQN2[4];
tmpQN1[5] = + tmpQN2[5];
tmpQN1[6] = + tmpQN2[6];
tmpQN1[7] = + tmpQN2[7];
tmpQN1[8] = + tmpQN2[8];
tmpQN1[9] = + tmpQN2[9];
tmpQN1[10] = + tmpQN2[10];
tmpQN1[11] = + tmpQN2[11];
tmpQN1[12] = + tmpQN2[12];
tmpQN1[13] = + tmpQN2[13];
tmpQN1[14] = + tmpQN2[14];
tmpQN1[15] = + tmpQN2[15];
tmpQN1[16] = + tmpQN2[16];
tmpQN1[17] = + tmpQN2[17];
tmpQN1[18] = + tmpQN2[18];
tmpQN1[19] = + tmpQN2[19];
tmpQN1[20] = + tmpQN2[20];
tmpQN1[21] = + tmpQN2[21];
tmpQN1[22] = + tmpQN2[22];
tmpQN1[23] = + tmpQN2[23];
tmpQN1[24] = + tmpQN2[24];
tmpQN1[25] = + tmpQN2[25];
tmpQN1[26] = + tmpQN2[26];
tmpQN1[27] = + tmpQN2[27];
tmpQN1[28] = + tmpQN2[28];
tmpQN1[29] = + tmpQN2[29];
tmpQN1[30] = + tmpQN2[30];
tmpQN1[31] = + tmpQN2[31];
tmpQN1[32] = + tmpQN2[32];
tmpQN1[33] = + tmpQN2[33];
tmpQN1[34] = + tmpQN2[34];
tmpQN1[35] = + tmpQN2[35];
}

void acado_evaluateObjective(  )
{
int runObj;
for (runObj = 0; runObj < 10; ++runObj)
{
acadoWorkspace.objValueIn[0] = acadoVariables.x[runObj * 6];
acadoWorkspace.objValueIn[1] = acadoVariables.x[runObj * 6 + 1];
acadoWorkspace.objValueIn[2] = acadoVariables.x[runObj * 6 + 2];
acadoWorkspace.objValueIn[3] = acadoVariables.x[runObj * 6 + 3];
acadoWorkspace.objValueIn[4] = acadoVariables.x[runObj * 6 + 4];
acadoWorkspace.objValueIn[5] = acadoVariables.x[runObj * 6 + 5];
acadoWorkspace.objValueIn[6] = acadoVariables.u[runObj * 3];
acadoWorkspace.objValueIn[7] = acadoVariables.u[runObj * 3 + 1];
acadoWorkspace.objValueIn[8] = acadoVariables.u[runObj * 3 + 2];
acadoWorkspace.objValueIn[9] = acadoVariables.od[runObj * 4];
acadoWorkspace.objValueIn[10] = acadoVariables.od[runObj * 4 + 1];
acadoWorkspace.objValueIn[11] = acadoVariables.od[runObj * 4 + 2];
acadoWorkspace.objValueIn[12] = acadoVariables.od[runObj * 4 + 3];

acado_evaluateLSQ( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.Dy[runObj * 9] = acadoWorkspace.objValueOut[0];
acadoWorkspace.Dy[runObj * 9 + 1] = acadoWorkspace.objValueOut[1];
acadoWorkspace.Dy[runObj * 9 + 2] = acadoWorkspace.objValueOut[2];
acadoWorkspace.Dy[runObj * 9 + 3] = acadoWorkspace.objValueOut[3];
acadoWorkspace.Dy[runObj * 9 + 4] = acadoWorkspace.objValueOut[4];
acadoWorkspace.Dy[runObj * 9 + 5] = acadoWorkspace.objValueOut[5];
acadoWorkspace.Dy[runObj * 9 + 6] = acadoWorkspace.objValueOut[6];
acadoWorkspace.Dy[runObj * 9 + 7] = acadoWorkspace.objValueOut[7];
acadoWorkspace.Dy[runObj * 9 + 8] = acadoWorkspace.objValueOut[8];

acado_setObjQ1Q2( acadoVariables.W, &(acadoWorkspace.Q1[ runObj * 36 ]), &(acadoWorkspace.Q2[ runObj * 54 ]) );

acado_setObjR1R2( acadoVariables.W, &(acadoWorkspace.R1[ runObj * 9 ]), &(acadoWorkspace.R2[ runObj * 27 ]) );

}
acadoWorkspace.objValueIn[0] = acadoVariables.x[60];
acadoWorkspace.objValueIn[1] = acadoVariables.x[61];
acadoWorkspace.objValueIn[2] = acadoVariables.x[62];
acadoWorkspace.objValueIn[3] = acadoVariables.x[63];
acadoWorkspace.objValueIn[4] = acadoVariables.x[64];
acadoWorkspace.objValueIn[5] = acadoVariables.x[65];
acadoWorkspace.objValueIn[6] = acadoVariables.od[40];
acadoWorkspace.objValueIn[7] = acadoVariables.od[41];
acadoWorkspace.objValueIn[8] = acadoVariables.od[42];
acadoWorkspace.objValueIn[9] = acadoVariables.od[43];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );

acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0];
acadoWorkspace.DyN[1] = acadoWorkspace.objValueOut[1];
acadoWorkspace.DyN[2] = acadoWorkspace.objValueOut[2];
acadoWorkspace.DyN[3] = acadoWorkspace.objValueOut[3];
acadoWorkspace.DyN[4] = acadoWorkspace.objValueOut[4];
acadoWorkspace.DyN[5] = acadoWorkspace.objValueOut[5];

acado_setObjQN1QN2( acadoVariables.WN, acadoWorkspace.QN1, acadoWorkspace.QN2 );

}

void acado_multGxGu( real_t* const Gx1, real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = + Gx1[0]*Gu1[0] + Gx1[1]*Gu1[3] + Gx1[2]*Gu1[6] + Gx1[3]*Gu1[9] + Gx1[4]*Gu1[12] + Gx1[5]*Gu1[15];
Gu2[1] = + Gx1[0]*Gu1[1] + Gx1[1]*Gu1[4] + Gx1[2]*Gu1[7] + Gx1[3]*Gu1[10] + Gx1[4]*Gu1[13] + Gx1[5]*Gu1[16];
Gu2[2] = + Gx1[0]*Gu1[2] + Gx1[1]*Gu1[5] + Gx1[2]*Gu1[8] + Gx1[3]*Gu1[11] + Gx1[4]*Gu1[14] + Gx1[5]*Gu1[17];
Gu2[3] = + Gx1[6]*Gu1[0] + Gx1[7]*Gu1[3] + Gx1[8]*Gu1[6] + Gx1[9]*Gu1[9] + Gx1[10]*Gu1[12] + Gx1[11]*Gu1[15];
Gu2[4] = + Gx1[6]*Gu1[1] + Gx1[7]*Gu1[4] + Gx1[8]*Gu1[7] + Gx1[9]*Gu1[10] + Gx1[10]*Gu1[13] + Gx1[11]*Gu1[16];
Gu2[5] = + Gx1[6]*Gu1[2] + Gx1[7]*Gu1[5] + Gx1[8]*Gu1[8] + Gx1[9]*Gu1[11] + Gx1[10]*Gu1[14] + Gx1[11]*Gu1[17];
Gu2[6] = + Gx1[12]*Gu1[0] + Gx1[13]*Gu1[3] + Gx1[14]*Gu1[6] + Gx1[15]*Gu1[9] + Gx1[16]*Gu1[12] + Gx1[17]*Gu1[15];
Gu2[7] = + Gx1[12]*Gu1[1] + Gx1[13]*Gu1[4] + Gx1[14]*Gu1[7] + Gx1[15]*Gu1[10] + Gx1[16]*Gu1[13] + Gx1[17]*Gu1[16];
Gu2[8] = + Gx1[12]*Gu1[2] + Gx1[13]*Gu1[5] + Gx1[14]*Gu1[8] + Gx1[15]*Gu1[11] + Gx1[16]*Gu1[14] + Gx1[17]*Gu1[17];
Gu2[9] = + Gx1[18]*Gu1[0] + Gx1[19]*Gu1[3] + Gx1[20]*Gu1[6] + Gx1[21]*Gu1[9] + Gx1[22]*Gu1[12] + Gx1[23]*Gu1[15];
Gu2[10] = + Gx1[18]*Gu1[1] + Gx1[19]*Gu1[4] + Gx1[20]*Gu1[7] + Gx1[21]*Gu1[10] + Gx1[22]*Gu1[13] + Gx1[23]*Gu1[16];
Gu2[11] = + Gx1[18]*Gu1[2] + Gx1[19]*Gu1[5] + Gx1[20]*Gu1[8] + Gx1[21]*Gu1[11] + Gx1[22]*Gu1[14] + Gx1[23]*Gu1[17];
Gu2[12] = + Gx1[24]*Gu1[0] + Gx1[25]*Gu1[3] + Gx1[26]*Gu1[6] + Gx1[27]*Gu1[9] + Gx1[28]*Gu1[12] + Gx1[29]*Gu1[15];
Gu2[13] = + Gx1[24]*Gu1[1] + Gx1[25]*Gu1[4] + Gx1[26]*Gu1[7] + Gx1[27]*Gu1[10] + Gx1[28]*Gu1[13] + Gx1[29]*Gu1[16];
Gu2[14] = + Gx1[24]*Gu1[2] + Gx1[25]*Gu1[5] + Gx1[26]*Gu1[8] + Gx1[27]*Gu1[11] + Gx1[28]*Gu1[14] + Gx1[29]*Gu1[17];
Gu2[15] = + Gx1[30]*Gu1[0] + Gx1[31]*Gu1[3] + Gx1[32]*Gu1[6] + Gx1[33]*Gu1[9] + Gx1[34]*Gu1[12] + Gx1[35]*Gu1[15];
Gu2[16] = + Gx1[30]*Gu1[1] + Gx1[31]*Gu1[4] + Gx1[32]*Gu1[7] + Gx1[33]*Gu1[10] + Gx1[34]*Gu1[13] + Gx1[35]*Gu1[16];
Gu2[17] = + Gx1[30]*Gu1[2] + Gx1[31]*Gu1[5] + Gx1[32]*Gu1[8] + Gx1[33]*Gu1[11] + Gx1[34]*Gu1[14] + Gx1[35]*Gu1[17];
}

void acado_moveGuE( real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = Gu1[0];
Gu2[1] = Gu1[1];
Gu2[2] = Gu1[2];
Gu2[3] = Gu1[3];
Gu2[4] = Gu1[4];
Gu2[5] = Gu1[5];
Gu2[6] = Gu1[6];
Gu2[7] = Gu1[7];
Gu2[8] = Gu1[8];
Gu2[9] = Gu1[9];
Gu2[10] = Gu1[10];
Gu2[11] = Gu1[11];
Gu2[12] = Gu1[12];
Gu2[13] = Gu1[13];
Gu2[14] = Gu1[14];
Gu2[15] = Gu1[15];
Gu2[16] = Gu1[16];
Gu2[17] = Gu1[17];
}

void acado_multBTW1( real_t* const Gu1, real_t* const Gu2, int iRow, int iCol )
{
acadoWorkspace.H[(iRow * 90) + (iCol * 3)] = + Gu1[0]*Gu2[0] + Gu1[3]*Gu2[3] + Gu1[6]*Gu2[6] + Gu1[9]*Gu2[9] + Gu1[12]*Gu2[12] + Gu1[15]*Gu2[15];
acadoWorkspace.H[(iRow * 90) + (iCol * 3 + 1)] = + Gu1[0]*Gu2[1] + Gu1[3]*Gu2[4] + Gu1[6]*Gu2[7] + Gu1[9]*Gu2[10] + Gu1[12]*Gu2[13] + Gu1[15]*Gu2[16];
acadoWorkspace.H[(iRow * 90) + (iCol * 3 + 2)] = + Gu1[0]*Gu2[2] + Gu1[3]*Gu2[5] + Gu1[6]*Gu2[8] + Gu1[9]*Gu2[11] + Gu1[12]*Gu2[14] + Gu1[15]*Gu2[17];
acadoWorkspace.H[(iRow * 90 + 30) + (iCol * 3)] = + Gu1[1]*Gu2[0] + Gu1[4]*Gu2[3] + Gu1[7]*Gu2[6] + Gu1[10]*Gu2[9] + Gu1[13]*Gu2[12] + Gu1[16]*Gu2[15];
acadoWorkspace.H[(iRow * 90 + 30) + (iCol * 3 + 1)] = + Gu1[1]*Gu2[1] + Gu1[4]*Gu2[4] + Gu1[7]*Gu2[7] + Gu1[10]*Gu2[10] + Gu1[13]*Gu2[13] + Gu1[16]*Gu2[16];
acadoWorkspace.H[(iRow * 90 + 30) + (iCol * 3 + 2)] = + Gu1[1]*Gu2[2] + Gu1[4]*Gu2[5] + Gu1[7]*Gu2[8] + Gu1[10]*Gu2[11] + Gu1[13]*Gu2[14] + Gu1[16]*Gu2[17];
acadoWorkspace.H[(iRow * 90 + 60) + (iCol * 3)] = + Gu1[2]*Gu2[0] + Gu1[5]*Gu2[3] + Gu1[8]*Gu2[6] + Gu1[11]*Gu2[9] + Gu1[14]*Gu2[12] + Gu1[17]*Gu2[15];
acadoWorkspace.H[(iRow * 90 + 60) + (iCol * 3 + 1)] = + Gu1[2]*Gu2[1] + Gu1[5]*Gu2[4] + Gu1[8]*Gu2[7] + Gu1[11]*Gu2[10] + Gu1[14]*Gu2[13] + Gu1[17]*Gu2[16];
acadoWorkspace.H[(iRow * 90 + 60) + (iCol * 3 + 2)] = + Gu1[2]*Gu2[2] + Gu1[5]*Gu2[5] + Gu1[8]*Gu2[8] + Gu1[11]*Gu2[11] + Gu1[14]*Gu2[14] + Gu1[17]*Gu2[17];
}

void acado_multBTW1_R1( real_t* const R11, real_t* const Gu1, real_t* const Gu2, int iRow )
{
acadoWorkspace.H[iRow * 93] = + Gu1[0]*Gu2[0] + Gu1[3]*Gu2[3] + Gu1[6]*Gu2[6] + Gu1[9]*Gu2[9] + Gu1[12]*Gu2[12] + Gu1[15]*Gu2[15] + R11[0];
acadoWorkspace.H[iRow * 93 + 1] = + Gu1[0]*Gu2[1] + Gu1[3]*Gu2[4] + Gu1[6]*Gu2[7] + Gu1[9]*Gu2[10] + Gu1[12]*Gu2[13] + Gu1[15]*Gu2[16] + R11[1];
acadoWorkspace.H[iRow * 93 + 2] = + Gu1[0]*Gu2[2] + Gu1[3]*Gu2[5] + Gu1[6]*Gu2[8] + Gu1[9]*Gu2[11] + Gu1[12]*Gu2[14] + Gu1[15]*Gu2[17] + R11[2];
acadoWorkspace.H[iRow * 93 + 30] = + Gu1[1]*Gu2[0] + Gu1[4]*Gu2[3] + Gu1[7]*Gu2[6] + Gu1[10]*Gu2[9] + Gu1[13]*Gu2[12] + Gu1[16]*Gu2[15] + R11[3];
acadoWorkspace.H[iRow * 93 + 31] = + Gu1[1]*Gu2[1] + Gu1[4]*Gu2[4] + Gu1[7]*Gu2[7] + Gu1[10]*Gu2[10] + Gu1[13]*Gu2[13] + Gu1[16]*Gu2[16] + R11[4];
acadoWorkspace.H[iRow * 93 + 32] = + Gu1[1]*Gu2[2] + Gu1[4]*Gu2[5] + Gu1[7]*Gu2[8] + Gu1[10]*Gu2[11] + Gu1[13]*Gu2[14] + Gu1[16]*Gu2[17] + R11[5];
acadoWorkspace.H[iRow * 93 + 60] = + Gu1[2]*Gu2[0] + Gu1[5]*Gu2[3] + Gu1[8]*Gu2[6] + Gu1[11]*Gu2[9] + Gu1[14]*Gu2[12] + Gu1[17]*Gu2[15] + R11[6];
acadoWorkspace.H[iRow * 93 + 61] = + Gu1[2]*Gu2[1] + Gu1[5]*Gu2[4] + Gu1[8]*Gu2[7] + Gu1[11]*Gu2[10] + Gu1[14]*Gu2[13] + Gu1[17]*Gu2[16] + R11[7];
acadoWorkspace.H[iRow * 93 + 62] = + Gu1[2]*Gu2[2] + Gu1[5]*Gu2[5] + Gu1[8]*Gu2[8] + Gu1[11]*Gu2[11] + Gu1[14]*Gu2[14] + Gu1[17]*Gu2[17] + R11[8];
acadoWorkspace.H[iRow * 93] += 1.0000000000000001e-05;
acadoWorkspace.H[iRow * 93 + 31] += 1.0000000000000001e-05;
acadoWorkspace.H[iRow * 93 + 62] += 1.0000000000000001e-05;
}

void acado_multGxTGu( real_t* const Gx1, real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = + Gx1[0]*Gu1[0] + Gx1[6]*Gu1[3] + Gx1[12]*Gu1[6] + Gx1[18]*Gu1[9] + Gx1[24]*Gu1[12] + Gx1[30]*Gu1[15];
Gu2[1] = + Gx1[0]*Gu1[1] + Gx1[6]*Gu1[4] + Gx1[12]*Gu1[7] + Gx1[18]*Gu1[10] + Gx1[24]*Gu1[13] + Gx1[30]*Gu1[16];
Gu2[2] = + Gx1[0]*Gu1[2] + Gx1[6]*Gu1[5] + Gx1[12]*Gu1[8] + Gx1[18]*Gu1[11] + Gx1[24]*Gu1[14] + Gx1[30]*Gu1[17];
Gu2[3] = + Gx1[1]*Gu1[0] + Gx1[7]*Gu1[3] + Gx1[13]*Gu1[6] + Gx1[19]*Gu1[9] + Gx1[25]*Gu1[12] + Gx1[31]*Gu1[15];
Gu2[4] = + Gx1[1]*Gu1[1] + Gx1[7]*Gu1[4] + Gx1[13]*Gu1[7] + Gx1[19]*Gu1[10] + Gx1[25]*Gu1[13] + Gx1[31]*Gu1[16];
Gu2[5] = + Gx1[1]*Gu1[2] + Gx1[7]*Gu1[5] + Gx1[13]*Gu1[8] + Gx1[19]*Gu1[11] + Gx1[25]*Gu1[14] + Gx1[31]*Gu1[17];
Gu2[6] = + Gx1[2]*Gu1[0] + Gx1[8]*Gu1[3] + Gx1[14]*Gu1[6] + Gx1[20]*Gu1[9] + Gx1[26]*Gu1[12] + Gx1[32]*Gu1[15];
Gu2[7] = + Gx1[2]*Gu1[1] + Gx1[8]*Gu1[4] + Gx1[14]*Gu1[7] + Gx1[20]*Gu1[10] + Gx1[26]*Gu1[13] + Gx1[32]*Gu1[16];
Gu2[8] = + Gx1[2]*Gu1[2] + Gx1[8]*Gu1[5] + Gx1[14]*Gu1[8] + Gx1[20]*Gu1[11] + Gx1[26]*Gu1[14] + Gx1[32]*Gu1[17];
Gu2[9] = + Gx1[3]*Gu1[0] + Gx1[9]*Gu1[3] + Gx1[15]*Gu1[6] + Gx1[21]*Gu1[9] + Gx1[27]*Gu1[12] + Gx1[33]*Gu1[15];
Gu2[10] = + Gx1[3]*Gu1[1] + Gx1[9]*Gu1[4] + Gx1[15]*Gu1[7] + Gx1[21]*Gu1[10] + Gx1[27]*Gu1[13] + Gx1[33]*Gu1[16];
Gu2[11] = + Gx1[3]*Gu1[2] + Gx1[9]*Gu1[5] + Gx1[15]*Gu1[8] + Gx1[21]*Gu1[11] + Gx1[27]*Gu1[14] + Gx1[33]*Gu1[17];
Gu2[12] = + Gx1[4]*Gu1[0] + Gx1[10]*Gu1[3] + Gx1[16]*Gu1[6] + Gx1[22]*Gu1[9] + Gx1[28]*Gu1[12] + Gx1[34]*Gu1[15];
Gu2[13] = + Gx1[4]*Gu1[1] + Gx1[10]*Gu1[4] + Gx1[16]*Gu1[7] + Gx1[22]*Gu1[10] + Gx1[28]*Gu1[13] + Gx1[34]*Gu1[16];
Gu2[14] = + Gx1[4]*Gu1[2] + Gx1[10]*Gu1[5] + Gx1[16]*Gu1[8] + Gx1[22]*Gu1[11] + Gx1[28]*Gu1[14] + Gx1[34]*Gu1[17];
Gu2[15] = + Gx1[5]*Gu1[0] + Gx1[11]*Gu1[3] + Gx1[17]*Gu1[6] + Gx1[23]*Gu1[9] + Gx1[29]*Gu1[12] + Gx1[35]*Gu1[15];
Gu2[16] = + Gx1[5]*Gu1[1] + Gx1[11]*Gu1[4] + Gx1[17]*Gu1[7] + Gx1[23]*Gu1[10] + Gx1[29]*Gu1[13] + Gx1[35]*Gu1[16];
Gu2[17] = + Gx1[5]*Gu1[2] + Gx1[11]*Gu1[5] + Gx1[17]*Gu1[8] + Gx1[23]*Gu1[11] + Gx1[29]*Gu1[14] + Gx1[35]*Gu1[17];
}

void acado_multQEW2( real_t* const Q11, real_t* const Gu1, real_t* const Gu2, real_t* const Gu3 )
{
Gu3[0] = + Q11[0]*Gu1[0] + Q11[1]*Gu1[3] + Q11[2]*Gu1[6] + Q11[3]*Gu1[9] + Q11[4]*Gu1[12] + Q11[5]*Gu1[15] + Gu2[0];
Gu3[1] = + Q11[0]*Gu1[1] + Q11[1]*Gu1[4] + Q11[2]*Gu1[7] + Q11[3]*Gu1[10] + Q11[4]*Gu1[13] + Q11[5]*Gu1[16] + Gu2[1];
Gu3[2] = + Q11[0]*Gu1[2] + Q11[1]*Gu1[5] + Q11[2]*Gu1[8] + Q11[3]*Gu1[11] + Q11[4]*Gu1[14] + Q11[5]*Gu1[17] + Gu2[2];
Gu3[3] = + Q11[6]*Gu1[0] + Q11[7]*Gu1[3] + Q11[8]*Gu1[6] + Q11[9]*Gu1[9] + Q11[10]*Gu1[12] + Q11[11]*Gu1[15] + Gu2[3];
Gu3[4] = + Q11[6]*Gu1[1] + Q11[7]*Gu1[4] + Q11[8]*Gu1[7] + Q11[9]*Gu1[10] + Q11[10]*Gu1[13] + Q11[11]*Gu1[16] + Gu2[4];
Gu3[5] = + Q11[6]*Gu1[2] + Q11[7]*Gu1[5] + Q11[8]*Gu1[8] + Q11[9]*Gu1[11] + Q11[10]*Gu1[14] + Q11[11]*Gu1[17] + Gu2[5];
Gu3[6] = + Q11[12]*Gu1[0] + Q11[13]*Gu1[3] + Q11[14]*Gu1[6] + Q11[15]*Gu1[9] + Q11[16]*Gu1[12] + Q11[17]*Gu1[15] + Gu2[6];
Gu3[7] = + Q11[12]*Gu1[1] + Q11[13]*Gu1[4] + Q11[14]*Gu1[7] + Q11[15]*Gu1[10] + Q11[16]*Gu1[13] + Q11[17]*Gu1[16] + Gu2[7];
Gu3[8] = + Q11[12]*Gu1[2] + Q11[13]*Gu1[5] + Q11[14]*Gu1[8] + Q11[15]*Gu1[11] + Q11[16]*Gu1[14] + Q11[17]*Gu1[17] + Gu2[8];
Gu3[9] = + Q11[18]*Gu1[0] + Q11[19]*Gu1[3] + Q11[20]*Gu1[6] + Q11[21]*Gu1[9] + Q11[22]*Gu1[12] + Q11[23]*Gu1[15] + Gu2[9];
Gu3[10] = + Q11[18]*Gu1[1] + Q11[19]*Gu1[4] + Q11[20]*Gu1[7] + Q11[21]*Gu1[10] + Q11[22]*Gu1[13] + Q11[23]*Gu1[16] + Gu2[10];
Gu3[11] = + Q11[18]*Gu1[2] + Q11[19]*Gu1[5] + Q11[20]*Gu1[8] + Q11[21]*Gu1[11] + Q11[22]*Gu1[14] + Q11[23]*Gu1[17] + Gu2[11];
Gu3[12] = + Q11[24]*Gu1[0] + Q11[25]*Gu1[3] + Q11[26]*Gu1[6] + Q11[27]*Gu1[9] + Q11[28]*Gu1[12] + Q11[29]*Gu1[15] + Gu2[12];
Gu3[13] = + Q11[24]*Gu1[1] + Q11[25]*Gu1[4] + Q11[26]*Gu1[7] + Q11[27]*Gu1[10] + Q11[28]*Gu1[13] + Q11[29]*Gu1[16] + Gu2[13];
Gu3[14] = + Q11[24]*Gu1[2] + Q11[25]*Gu1[5] + Q11[26]*Gu1[8] + Q11[27]*Gu1[11] + Q11[28]*Gu1[14] + Q11[29]*Gu1[17] + Gu2[14];
Gu3[15] = + Q11[30]*Gu1[0] + Q11[31]*Gu1[3] + Q11[32]*Gu1[6] + Q11[33]*Gu1[9] + Q11[34]*Gu1[12] + Q11[35]*Gu1[15] + Gu2[15];
Gu3[16] = + Q11[30]*Gu1[1] + Q11[31]*Gu1[4] + Q11[32]*Gu1[7] + Q11[33]*Gu1[10] + Q11[34]*Gu1[13] + Q11[35]*Gu1[16] + Gu2[16];
Gu3[17] = + Q11[30]*Gu1[2] + Q11[31]*Gu1[5] + Q11[32]*Gu1[8] + Q11[33]*Gu1[11] + Q11[34]*Gu1[14] + Q11[35]*Gu1[17] + Gu2[17];
}

void acado_macATw1QDy( real_t* const Gx1, real_t* const w11, real_t* const w12, real_t* const w13 )
{
w13[0] = + Gx1[0]*w11[0] + Gx1[6]*w11[1] + Gx1[12]*w11[2] + Gx1[18]*w11[3] + Gx1[24]*w11[4] + Gx1[30]*w11[5] + w12[0];
w13[1] = + Gx1[1]*w11[0] + Gx1[7]*w11[1] + Gx1[13]*w11[2] + Gx1[19]*w11[3] + Gx1[25]*w11[4] + Gx1[31]*w11[5] + w12[1];
w13[2] = + Gx1[2]*w11[0] + Gx1[8]*w11[1] + Gx1[14]*w11[2] + Gx1[20]*w11[3] + Gx1[26]*w11[4] + Gx1[32]*w11[5] + w12[2];
w13[3] = + Gx1[3]*w11[0] + Gx1[9]*w11[1] + Gx1[15]*w11[2] + Gx1[21]*w11[3] + Gx1[27]*w11[4] + Gx1[33]*w11[5] + w12[3];
w13[4] = + Gx1[4]*w11[0] + Gx1[10]*w11[1] + Gx1[16]*w11[2] + Gx1[22]*w11[3] + Gx1[28]*w11[4] + Gx1[34]*w11[5] + w12[4];
w13[5] = + Gx1[5]*w11[0] + Gx1[11]*w11[1] + Gx1[17]*w11[2] + Gx1[23]*w11[3] + Gx1[29]*w11[4] + Gx1[35]*w11[5] + w12[5];
}

void acado_macBTw1( real_t* const Gu1, real_t* const w11, real_t* const U1 )
{
U1[0] += + Gu1[0]*w11[0] + Gu1[3]*w11[1] + Gu1[6]*w11[2] + Gu1[9]*w11[3] + Gu1[12]*w11[4] + Gu1[15]*w11[5];
U1[1] += + Gu1[1]*w11[0] + Gu1[4]*w11[1] + Gu1[7]*w11[2] + Gu1[10]*w11[3] + Gu1[13]*w11[4] + Gu1[16]*w11[5];
U1[2] += + Gu1[2]*w11[0] + Gu1[5]*w11[1] + Gu1[8]*w11[2] + Gu1[11]*w11[3] + Gu1[14]*w11[4] + Gu1[17]*w11[5];
}

void acado_macQSbarW2( real_t* const Q11, real_t* const w11, real_t* const w12, real_t* const w13 )
{
w13[0] = + Q11[0]*w11[0] + Q11[1]*w11[1] + Q11[2]*w11[2] + Q11[3]*w11[3] + Q11[4]*w11[4] + Q11[5]*w11[5] + w12[0];
w13[1] = + Q11[6]*w11[0] + Q11[7]*w11[1] + Q11[8]*w11[2] + Q11[9]*w11[3] + Q11[10]*w11[4] + Q11[11]*w11[5] + w12[1];
w13[2] = + Q11[12]*w11[0] + Q11[13]*w11[1] + Q11[14]*w11[2] + Q11[15]*w11[3] + Q11[16]*w11[4] + Q11[17]*w11[5] + w12[2];
w13[3] = + Q11[18]*w11[0] + Q11[19]*w11[1] + Q11[20]*w11[2] + Q11[21]*w11[3] + Q11[22]*w11[4] + Q11[23]*w11[5] + w12[3];
w13[4] = + Q11[24]*w11[0] + Q11[25]*w11[1] + Q11[26]*w11[2] + Q11[27]*w11[3] + Q11[28]*w11[4] + Q11[29]*w11[5] + w12[4];
w13[5] = + Q11[30]*w11[0] + Q11[31]*w11[1] + Q11[32]*w11[2] + Q11[33]*w11[3] + Q11[34]*w11[4] + Q11[35]*w11[5] + w12[5];
}

void acado_macASbar( real_t* const Gx1, real_t* const w11, real_t* const w12 )
{
w12[0] += + Gx1[0]*w11[0] + Gx1[1]*w11[1] + Gx1[2]*w11[2] + Gx1[3]*w11[3] + Gx1[4]*w11[4] + Gx1[5]*w11[5];
w12[1] += + Gx1[6]*w11[0] + Gx1[7]*w11[1] + Gx1[8]*w11[2] + Gx1[9]*w11[3] + Gx1[10]*w11[4] + Gx1[11]*w11[5];
w12[2] += + Gx1[12]*w11[0] + Gx1[13]*w11[1] + Gx1[14]*w11[2] + Gx1[15]*w11[3] + Gx1[16]*w11[4] + Gx1[17]*w11[5];
w12[3] += + Gx1[18]*w11[0] + Gx1[19]*w11[1] + Gx1[20]*w11[2] + Gx1[21]*w11[3] + Gx1[22]*w11[4] + Gx1[23]*w11[5];
w12[4] += + Gx1[24]*w11[0] + Gx1[25]*w11[1] + Gx1[26]*w11[2] + Gx1[27]*w11[3] + Gx1[28]*w11[4] + Gx1[29]*w11[5];
w12[5] += + Gx1[30]*w11[0] + Gx1[31]*w11[1] + Gx1[32]*w11[2] + Gx1[33]*w11[3] + Gx1[34]*w11[4] + Gx1[35]*w11[5];
}

void acado_expansionStep( real_t* const Gx1, real_t* const Gu1, real_t* const U1, real_t* const w11, real_t* const w12 )
{
w12[0] += + Gx1[0]*w11[0] + Gx1[1]*w11[1] + Gx1[2]*w11[2] + Gx1[3]*w11[3] + Gx1[4]*w11[4] + Gx1[5]*w11[5];
w12[1] += + Gx1[6]*w11[0] + Gx1[7]*w11[1] + Gx1[8]*w11[2] + Gx1[9]*w11[3] + Gx1[10]*w11[4] + Gx1[11]*w11[5];
w12[2] += + Gx1[12]*w11[0] + Gx1[13]*w11[1] + Gx1[14]*w11[2] + Gx1[15]*w11[3] + Gx1[16]*w11[4] + Gx1[17]*w11[5];
w12[3] += + Gx1[18]*w11[0] + Gx1[19]*w11[1] + Gx1[20]*w11[2] + Gx1[21]*w11[3] + Gx1[22]*w11[4] + Gx1[23]*w11[5];
w12[4] += + Gx1[24]*w11[0] + Gx1[25]*w11[1] + Gx1[26]*w11[2] + Gx1[27]*w11[3] + Gx1[28]*w11[4] + Gx1[29]*w11[5];
w12[5] += + Gx1[30]*w11[0] + Gx1[31]*w11[1] + Gx1[32]*w11[2] + Gx1[33]*w11[3] + Gx1[34]*w11[4] + Gx1[35]*w11[5];
w12[0] += + Gu1[0]*U1[0] + Gu1[1]*U1[1] + Gu1[2]*U1[2];
w12[1] += + Gu1[3]*U1[0] + Gu1[4]*U1[1] + Gu1[5]*U1[2];
w12[2] += + Gu1[6]*U1[0] + Gu1[7]*U1[1] + Gu1[8]*U1[2];
w12[3] += + Gu1[9]*U1[0] + Gu1[10]*U1[1] + Gu1[11]*U1[2];
w12[4] += + Gu1[12]*U1[0] + Gu1[13]*U1[1] + Gu1[14]*U1[2];
w12[5] += + Gu1[15]*U1[0] + Gu1[16]*U1[1] + Gu1[17]*U1[2];
}

void acado_copyHTH( int iRow, int iCol )
{
acadoWorkspace.H[(iRow * 90) + (iCol * 3)] = acadoWorkspace.H[(iCol * 90) + (iRow * 3)];
acadoWorkspace.H[(iRow * 90) + (iCol * 3 + 1)] = acadoWorkspace.H[(iCol * 90 + 30) + (iRow * 3)];
acadoWorkspace.H[(iRow * 90) + (iCol * 3 + 2)] = acadoWorkspace.H[(iCol * 90 + 60) + (iRow * 3)];
acadoWorkspace.H[(iRow * 90 + 30) + (iCol * 3)] = acadoWorkspace.H[(iCol * 90) + (iRow * 3 + 1)];
acadoWorkspace.H[(iRow * 90 + 30) + (iCol * 3 + 1)] = acadoWorkspace.H[(iCol * 90 + 30) + (iRow * 3 + 1)];
acadoWorkspace.H[(iRow * 90 + 30) + (iCol * 3 + 2)] = acadoWorkspace.H[(iCol * 90 + 60) + (iRow * 3 + 1)];
acadoWorkspace.H[(iRow * 90 + 60) + (iCol * 3)] = acadoWorkspace.H[(iCol * 90) + (iRow * 3 + 2)];
acadoWorkspace.H[(iRow * 90 + 60) + (iCol * 3 + 1)] = acadoWorkspace.H[(iCol * 90 + 30) + (iRow * 3 + 2)];
acadoWorkspace.H[(iRow * 90 + 60) + (iCol * 3 + 2)] = acadoWorkspace.H[(iCol * 90 + 60) + (iRow * 3 + 2)];
}

void acado_multRDy( real_t* const R2, real_t* const Dy1, real_t* const RDy1 )
{
RDy1[0] = + R2[0]*Dy1[0] + R2[1]*Dy1[1] + R2[2]*Dy1[2] + R2[3]*Dy1[3] + R2[4]*Dy1[4] + R2[5]*Dy1[5] + R2[6]*Dy1[6] + R2[7]*Dy1[7] + R2[8]*Dy1[8];
RDy1[1] = + R2[9]*Dy1[0] + R2[10]*Dy1[1] + R2[11]*Dy1[2] + R2[12]*Dy1[3] + R2[13]*Dy1[4] + R2[14]*Dy1[5] + R2[15]*Dy1[6] + R2[16]*Dy1[7] + R2[17]*Dy1[8];
RDy1[2] = + R2[18]*Dy1[0] + R2[19]*Dy1[1] + R2[20]*Dy1[2] + R2[21]*Dy1[3] + R2[22]*Dy1[4] + R2[23]*Dy1[5] + R2[24]*Dy1[6] + R2[25]*Dy1[7] + R2[26]*Dy1[8];
}

void acado_multQDy( real_t* const Q2, real_t* const Dy1, real_t* const QDy1 )
{
QDy1[0] = + Q2[0]*Dy1[0] + Q2[1]*Dy1[1] + Q2[2]*Dy1[2] + Q2[3]*Dy1[3] + Q2[4]*Dy1[4] + Q2[5]*Dy1[5] + Q2[6]*Dy1[6] + Q2[7]*Dy1[7] + Q2[8]*Dy1[8];
QDy1[1] = + Q2[9]*Dy1[0] + Q2[10]*Dy1[1] + Q2[11]*Dy1[2] + Q2[12]*Dy1[3] + Q2[13]*Dy1[4] + Q2[14]*Dy1[5] + Q2[15]*Dy1[6] + Q2[16]*Dy1[7] + Q2[17]*Dy1[8];
QDy1[2] = + Q2[18]*Dy1[0] + Q2[19]*Dy1[1] + Q2[20]*Dy1[2] + Q2[21]*Dy1[3] + Q2[22]*Dy1[4] + Q2[23]*Dy1[5] + Q2[24]*Dy1[6] + Q2[25]*Dy1[7] + Q2[26]*Dy1[8];
QDy1[3] = + Q2[27]*Dy1[0] + Q2[28]*Dy1[1] + Q2[29]*Dy1[2] + Q2[30]*Dy1[3] + Q2[31]*Dy1[4] + Q2[32]*Dy1[5] + Q2[33]*Dy1[6] + Q2[34]*Dy1[7] + Q2[35]*Dy1[8];
QDy1[4] = + Q2[36]*Dy1[0] + Q2[37]*Dy1[1] + Q2[38]*Dy1[2] + Q2[39]*Dy1[3] + Q2[40]*Dy1[4] + Q2[41]*Dy1[5] + Q2[42]*Dy1[6] + Q2[43]*Dy1[7] + Q2[44]*Dy1[8];
QDy1[5] = + Q2[45]*Dy1[0] + Q2[46]*Dy1[1] + Q2[47]*Dy1[2] + Q2[48]*Dy1[3] + Q2[49]*Dy1[4] + Q2[50]*Dy1[5] + Q2[51]*Dy1[6] + Q2[52]*Dy1[7] + Q2[53]*Dy1[8];
}

void acado_condensePrep(  )
{
/* Column: 0 */
acado_moveGuE( acadoWorkspace.evGu, acadoWorkspace.E );
acado_multGxGu( &(acadoWorkspace.evGx[ 36 ]), acadoWorkspace.E, &(acadoWorkspace.E[ 18 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 72 ]), &(acadoWorkspace.E[ 18 ]), &(acadoWorkspace.E[ 36 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.E[ 36 ]), &(acadoWorkspace.E[ 54 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.E[ 54 ]), &(acadoWorkspace.E[ 72 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.E[ 72 ]), &(acadoWorkspace.E[ 90 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.E[ 90 ]), &(acadoWorkspace.E[ 108 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.E[ 108 ]), &(acadoWorkspace.E[ 126 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.E[ 126 ]), &(acadoWorkspace.E[ 144 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 144 ]), &(acadoWorkspace.E[ 162 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 162 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 144 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 126 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.W1, 7, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 108 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.W1, 6, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 90 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 90 ]), acadoWorkspace.W1, 5, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 180 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 72 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 72 ]), acadoWorkspace.W1, 4, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 144 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.E[ 54 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 54 ]), acadoWorkspace.W1, 3, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 108 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 108 ]), &(acadoWorkspace.E[ 36 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 36 ]), acadoWorkspace.W1, 2, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 72 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 72 ]), &(acadoWorkspace.E[ 18 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 18 ]), acadoWorkspace.W1, 1, 0 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 36 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 36 ]), acadoWorkspace.E, acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( acadoWorkspace.R1, acadoWorkspace.evGu, acadoWorkspace.W1, 0 );

/* Column: 1 */
acado_moveGuE( &(acadoWorkspace.evGu[ 18 ]), &(acadoWorkspace.E[ 180 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 72 ]), &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.E[ 198 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.E[ 216 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.E[ 216 ]), &(acadoWorkspace.E[ 234 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.E[ 234 ]), &(acadoWorkspace.E[ 252 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.E[ 252 ]), &(acadoWorkspace.E[ 270 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.E[ 288 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.E[ 306 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.E[ 324 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 324 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 1 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 306 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8, 1 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 288 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.W1, 7, 1 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 270 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.W1, 6, 1 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 252 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 90 ]), acadoWorkspace.W1, 5, 1 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 180 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 234 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 72 ]), acadoWorkspace.W1, 4, 1 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 144 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.E[ 216 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 54 ]), acadoWorkspace.W1, 3, 1 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 108 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 108 ]), &(acadoWorkspace.E[ 198 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 36 ]), acadoWorkspace.W1, 2, 1 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 72 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 72 ]), &(acadoWorkspace.E[ 180 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 9 ]), &(acadoWorkspace.evGu[ 18 ]), acadoWorkspace.W1, 1 );

/* Column: 2 */
acado_moveGuE( &(acadoWorkspace.evGu[ 36 ]), &(acadoWorkspace.E[ 342 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.E[ 342 ]), &(acadoWorkspace.E[ 360 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.E[ 360 ]), &(acadoWorkspace.E[ 378 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.E[ 396 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.E[ 414 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.E[ 432 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.E[ 450 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 450 ]), &(acadoWorkspace.E[ 468 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 468 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 450 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8, 2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 432 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.W1, 7, 2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 414 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.W1, 6, 2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 396 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 90 ]), acadoWorkspace.W1, 5, 2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 180 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 378 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 72 ]), acadoWorkspace.W1, 4, 2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 144 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.E[ 360 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 54 ]), acadoWorkspace.W1, 3, 2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 108 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 108 ]), &(acadoWorkspace.E[ 342 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 18 ]), &(acadoWorkspace.evGu[ 36 ]), acadoWorkspace.W1, 2 );

/* Column: 3 */
acado_moveGuE( &(acadoWorkspace.evGu[ 54 ]), &(acadoWorkspace.E[ 486 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.E[ 486 ]), &(acadoWorkspace.E[ 504 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.E[ 522 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.E[ 540 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.E[ 558 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.E[ 576 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.E[ 594 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 594 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 3 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 576 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8, 3 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 558 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.W1, 7, 3 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 540 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.W1, 6, 3 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 522 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 90 ]), acadoWorkspace.W1, 5, 3 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 180 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 504 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 72 ]), acadoWorkspace.W1, 4, 3 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 144 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.E[ 486 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 27 ]), &(acadoWorkspace.evGu[ 54 ]), acadoWorkspace.W1, 3 );

/* Column: 4 */
acado_moveGuE( &(acadoWorkspace.evGu[ 72 ]), &(acadoWorkspace.E[ 612 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.E[ 612 ]), &(acadoWorkspace.E[ 630 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.E[ 630 ]), &(acadoWorkspace.E[ 648 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.E[ 666 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.E[ 684 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.E[ 702 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 702 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 4 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 684 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8, 4 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 666 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.W1, 7, 4 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 648 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.W1, 6, 4 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 630 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 90 ]), acadoWorkspace.W1, 5, 4 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 180 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 612 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 36 ]), &(acadoWorkspace.evGu[ 72 ]), acadoWorkspace.W1, 4 );

/* Column: 5 */
acado_moveGuE( &(acadoWorkspace.evGu[ 90 ]), &(acadoWorkspace.E[ 720 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.E[ 738 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.E[ 756 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.E[ 756 ]), &(acadoWorkspace.E[ 774 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 774 ]), &(acadoWorkspace.E[ 792 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 792 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 5 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 774 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8, 5 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 756 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.W1, 7, 5 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 738 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.W1, 6, 5 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 720 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 45 ]), &(acadoWorkspace.evGu[ 90 ]), acadoWorkspace.W1, 5 );

/* Column: 6 */
acado_moveGuE( &(acadoWorkspace.evGu[ 108 ]), &(acadoWorkspace.E[ 810 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.E[ 828 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.E[ 846 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.E[ 864 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 864 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 6 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 846 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8, 6 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 828 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.W1, 7, 6 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 810 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 54 ]), &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.W1, 6 );

/* Column: 7 */
acado_moveGuE( &(acadoWorkspace.evGu[ 126 ]), &(acadoWorkspace.E[ 882 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.E[ 900 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.E[ 918 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 918 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 7 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 900 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8, 7 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 882 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 63 ]), &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.W1, 7 );

/* Column: 8 */
acado_moveGuE( &(acadoWorkspace.evGu[ 144 ]), &(acadoWorkspace.E[ 936 ]) );
acado_multGxGu( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.E[ 936 ]), &(acadoWorkspace.E[ 954 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 954 ]), acadoWorkspace.W1 );
acado_multBTW1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9, 8 );
acado_multGxTGu( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 936 ]), acadoWorkspace.W2, acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 72 ]), &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.W1, 8 );

/* Column: 9 */
acado_moveGuE( &(acadoWorkspace.evGu[ 162 ]), &(acadoWorkspace.E[ 972 ]) );

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 972 ]), acadoWorkspace.W1 );
acado_multBTW1_R1( &(acadoWorkspace.R1[ 81 ]), &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.W1, 9 );

acado_copyHTH( 0, 1 );
acado_copyHTH( 0, 2 );
acado_copyHTH( 1, 2 );
acado_copyHTH( 0, 3 );
acado_copyHTH( 1, 3 );
acado_copyHTH( 2, 3 );
acado_copyHTH( 0, 4 );
acado_copyHTH( 1, 4 );
acado_copyHTH( 2, 4 );
acado_copyHTH( 3, 4 );
acado_copyHTH( 0, 5 );
acado_copyHTH( 1, 5 );
acado_copyHTH( 2, 5 );
acado_copyHTH( 3, 5 );
acado_copyHTH( 4, 5 );
acado_copyHTH( 0, 6 );
acado_copyHTH( 1, 6 );
acado_copyHTH( 2, 6 );
acado_copyHTH( 3, 6 );
acado_copyHTH( 4, 6 );
acado_copyHTH( 5, 6 );
acado_copyHTH( 0, 7 );
acado_copyHTH( 1, 7 );
acado_copyHTH( 2, 7 );
acado_copyHTH( 3, 7 );
acado_copyHTH( 4, 7 );
acado_copyHTH( 5, 7 );
acado_copyHTH( 6, 7 );
acado_copyHTH( 0, 8 );
acado_copyHTH( 1, 8 );
acado_copyHTH( 2, 8 );
acado_copyHTH( 3, 8 );
acado_copyHTH( 4, 8 );
acado_copyHTH( 5, 8 );
acado_copyHTH( 6, 8 );
acado_copyHTH( 7, 8 );
acado_copyHTH( 0, 9 );
acado_copyHTH( 1, 9 );
acado_copyHTH( 2, 9 );
acado_copyHTH( 3, 9 );
acado_copyHTH( 4, 9 );
acado_copyHTH( 5, 9 );
acado_copyHTH( 6, 9 );
acado_copyHTH( 7, 9 );
acado_copyHTH( 8, 9 );

acadoWorkspace.sbar[6] = acadoWorkspace.d[0];
acadoWorkspace.sbar[7] = acadoWorkspace.d[1];
acadoWorkspace.sbar[8] = acadoWorkspace.d[2];
acadoWorkspace.sbar[9] = acadoWorkspace.d[3];
acadoWorkspace.sbar[10] = acadoWorkspace.d[4];
acadoWorkspace.sbar[11] = acadoWorkspace.d[5];
acadoWorkspace.sbar[12] = acadoWorkspace.d[6];
acadoWorkspace.sbar[13] = acadoWorkspace.d[7];
acadoWorkspace.sbar[14] = acadoWorkspace.d[8];
acadoWorkspace.sbar[15] = acadoWorkspace.d[9];
acadoWorkspace.sbar[16] = acadoWorkspace.d[10];
acadoWorkspace.sbar[17] = acadoWorkspace.d[11];
acadoWorkspace.sbar[18] = acadoWorkspace.d[12];
acadoWorkspace.sbar[19] = acadoWorkspace.d[13];
acadoWorkspace.sbar[20] = acadoWorkspace.d[14];
acadoWorkspace.sbar[21] = acadoWorkspace.d[15];
acadoWorkspace.sbar[22] = acadoWorkspace.d[16];
acadoWorkspace.sbar[23] = acadoWorkspace.d[17];
acadoWorkspace.sbar[24] = acadoWorkspace.d[18];
acadoWorkspace.sbar[25] = acadoWorkspace.d[19];
acadoWorkspace.sbar[26] = acadoWorkspace.d[20];
acadoWorkspace.sbar[27] = acadoWorkspace.d[21];
acadoWorkspace.sbar[28] = acadoWorkspace.d[22];
acadoWorkspace.sbar[29] = acadoWorkspace.d[23];
acadoWorkspace.sbar[30] = acadoWorkspace.d[24];
acadoWorkspace.sbar[31] = acadoWorkspace.d[25];
acadoWorkspace.sbar[32] = acadoWorkspace.d[26];
acadoWorkspace.sbar[33] = acadoWorkspace.d[27];
acadoWorkspace.sbar[34] = acadoWorkspace.d[28];
acadoWorkspace.sbar[35] = acadoWorkspace.d[29];
acadoWorkspace.sbar[36] = acadoWorkspace.d[30];
acadoWorkspace.sbar[37] = acadoWorkspace.d[31];
acadoWorkspace.sbar[38] = acadoWorkspace.d[32];
acadoWorkspace.sbar[39] = acadoWorkspace.d[33];
acadoWorkspace.sbar[40] = acadoWorkspace.d[34];
acadoWorkspace.sbar[41] = acadoWorkspace.d[35];
acadoWorkspace.sbar[42] = acadoWorkspace.d[36];
acadoWorkspace.sbar[43] = acadoWorkspace.d[37];
acadoWorkspace.sbar[44] = acadoWorkspace.d[38];
acadoWorkspace.sbar[45] = acadoWorkspace.d[39];
acadoWorkspace.sbar[46] = acadoWorkspace.d[40];
acadoWorkspace.sbar[47] = acadoWorkspace.d[41];
acadoWorkspace.sbar[48] = acadoWorkspace.d[42];
acadoWorkspace.sbar[49] = acadoWorkspace.d[43];
acadoWorkspace.sbar[50] = acadoWorkspace.d[44];
acadoWorkspace.sbar[51] = acadoWorkspace.d[45];
acadoWorkspace.sbar[52] = acadoWorkspace.d[46];
acadoWorkspace.sbar[53] = acadoWorkspace.d[47];
acadoWorkspace.sbar[54] = acadoWorkspace.d[48];
acadoWorkspace.sbar[55] = acadoWorkspace.d[49];
acadoWorkspace.sbar[56] = acadoWorkspace.d[50];
acadoWorkspace.sbar[57] = acadoWorkspace.d[51];
acadoWorkspace.sbar[58] = acadoWorkspace.d[52];
acadoWorkspace.sbar[59] = acadoWorkspace.d[53];
acadoWorkspace.sbar[60] = acadoWorkspace.d[54];
acadoWorkspace.sbar[61] = acadoWorkspace.d[55];
acadoWorkspace.sbar[62] = acadoWorkspace.d[56];
acadoWorkspace.sbar[63] = acadoWorkspace.d[57];
acadoWorkspace.sbar[64] = acadoWorkspace.d[58];
acadoWorkspace.sbar[65] = acadoWorkspace.d[59];
acadoWorkspace.lb[0] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[0];
acadoWorkspace.lb[1] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[1];
acadoWorkspace.lb[2] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[2];
acadoWorkspace.lb[3] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[3];
acadoWorkspace.lb[4] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[4];
acadoWorkspace.lb[5] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[5];
acadoWorkspace.lb[6] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[6];
acadoWorkspace.lb[7] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[7];
acadoWorkspace.lb[8] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[8];
acadoWorkspace.lb[9] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[9];
acadoWorkspace.lb[10] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[10];
acadoWorkspace.lb[11] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[11];
acadoWorkspace.lb[12] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[12];
acadoWorkspace.lb[13] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[13];
acadoWorkspace.lb[14] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[14];
acadoWorkspace.lb[15] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[15];
acadoWorkspace.lb[16] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[16];
acadoWorkspace.lb[17] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[17];
acadoWorkspace.lb[18] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[18];
acadoWorkspace.lb[19] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[19];
acadoWorkspace.lb[20] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[20];
acadoWorkspace.lb[21] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[21];
acadoWorkspace.lb[22] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[22];
acadoWorkspace.lb[23] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[23];
acadoWorkspace.lb[24] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[24];
acadoWorkspace.lb[25] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[25];
acadoWorkspace.lb[26] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[26];
acadoWorkspace.lb[27] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[27];
acadoWorkspace.lb[28] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[28];
acadoWorkspace.lb[29] = (real_t)-1.0000000000000000e+12 - acadoVariables.u[29];
acadoWorkspace.ub[0] = (real_t)1.0000000000000000e+12 - acadoVariables.u[0];
acadoWorkspace.ub[1] = (real_t)1.0000000000000000e+12 - acadoVariables.u[1];
acadoWorkspace.ub[2] = (real_t)1.0000000000000000e+12 - acadoVariables.u[2];
acadoWorkspace.ub[3] = (real_t)1.0000000000000000e+12 - acadoVariables.u[3];
acadoWorkspace.ub[4] = (real_t)1.0000000000000000e+12 - acadoVariables.u[4];
acadoWorkspace.ub[5] = (real_t)1.0000000000000000e+12 - acadoVariables.u[5];
acadoWorkspace.ub[6] = (real_t)1.0000000000000000e+12 - acadoVariables.u[6];
acadoWorkspace.ub[7] = (real_t)1.0000000000000000e+12 - acadoVariables.u[7];
acadoWorkspace.ub[8] = (real_t)1.0000000000000000e+12 - acadoVariables.u[8];
acadoWorkspace.ub[9] = (real_t)1.0000000000000000e+12 - acadoVariables.u[9];
acadoWorkspace.ub[10] = (real_t)1.0000000000000000e+12 - acadoVariables.u[10];
acadoWorkspace.ub[11] = (real_t)1.0000000000000000e+12 - acadoVariables.u[11];
acadoWorkspace.ub[12] = (real_t)1.0000000000000000e+12 - acadoVariables.u[12];
acadoWorkspace.ub[13] = (real_t)1.0000000000000000e+12 - acadoVariables.u[13];
acadoWorkspace.ub[14] = (real_t)1.0000000000000000e+12 - acadoVariables.u[14];
acadoWorkspace.ub[15] = (real_t)1.0000000000000000e+12 - acadoVariables.u[15];
acadoWorkspace.ub[16] = (real_t)1.0000000000000000e+12 - acadoVariables.u[16];
acadoWorkspace.ub[17] = (real_t)1.0000000000000000e+12 - acadoVariables.u[17];
acadoWorkspace.ub[18] = (real_t)1.0000000000000000e+12 - acadoVariables.u[18];
acadoWorkspace.ub[19] = (real_t)1.0000000000000000e+12 - acadoVariables.u[19];
acadoWorkspace.ub[20] = (real_t)1.0000000000000000e+12 - acadoVariables.u[20];
acadoWorkspace.ub[21] = (real_t)1.0000000000000000e+12 - acadoVariables.u[21];
acadoWorkspace.ub[22] = (real_t)1.0000000000000000e+12 - acadoVariables.u[22];
acadoWorkspace.ub[23] = (real_t)1.0000000000000000e+12 - acadoVariables.u[23];
acadoWorkspace.ub[24] = (real_t)1.0000000000000000e+12 - acadoVariables.u[24];
acadoWorkspace.ub[25] = (real_t)1.0000000000000000e+12 - acadoVariables.u[25];
acadoWorkspace.ub[26] = (real_t)1.0000000000000000e+12 - acadoVariables.u[26];
acadoWorkspace.ub[27] = (real_t)1.0000000000000000e+12 - acadoVariables.u[27];
acadoWorkspace.ub[28] = (real_t)1.0000000000000000e+12 - acadoVariables.u[28];
acadoWorkspace.ub[29] = (real_t)1.0000000000000000e+12 - acadoVariables.u[29];

}

void acado_condenseFdb(  )
{
acadoWorkspace.Dx0[0] = acadoVariables.x0[0] - acadoVariables.x[0];
acadoWorkspace.Dx0[1] = acadoVariables.x0[1] - acadoVariables.x[1];
acadoWorkspace.Dx0[2] = acadoVariables.x0[2] - acadoVariables.x[2];
acadoWorkspace.Dx0[3] = acadoVariables.x0[3] - acadoVariables.x[3];
acadoWorkspace.Dx0[4] = acadoVariables.x0[4] - acadoVariables.x[4];
acadoWorkspace.Dx0[5] = acadoVariables.x0[5] - acadoVariables.x[5];
acadoWorkspace.Dy[0] -= acadoVariables.y[0];
acadoWorkspace.Dy[1] -= acadoVariables.y[1];
acadoWorkspace.Dy[2] -= acadoVariables.y[2];
acadoWorkspace.Dy[3] -= acadoVariables.y[3];
acadoWorkspace.Dy[4] -= acadoVariables.y[4];
acadoWorkspace.Dy[5] -= acadoVariables.y[5];
acadoWorkspace.Dy[6] -= acadoVariables.y[6];
acadoWorkspace.Dy[7] -= acadoVariables.y[7];
acadoWorkspace.Dy[8] -= acadoVariables.y[8];
acadoWorkspace.Dy[9] -= acadoVariables.y[9];
acadoWorkspace.Dy[10] -= acadoVariables.y[10];
acadoWorkspace.Dy[11] -= acadoVariables.y[11];
acadoWorkspace.Dy[12] -= acadoVariables.y[12];
acadoWorkspace.Dy[13] -= acadoVariables.y[13];
acadoWorkspace.Dy[14] -= acadoVariables.y[14];
acadoWorkspace.Dy[15] -= acadoVariables.y[15];
acadoWorkspace.Dy[16] -= acadoVariables.y[16];
acadoWorkspace.Dy[17] -= acadoVariables.y[17];
acadoWorkspace.Dy[18] -= acadoVariables.y[18];
acadoWorkspace.Dy[19] -= acadoVariables.y[19];
acadoWorkspace.Dy[20] -= acadoVariables.y[20];
acadoWorkspace.Dy[21] -= acadoVariables.y[21];
acadoWorkspace.Dy[22] -= acadoVariables.y[22];
acadoWorkspace.Dy[23] -= acadoVariables.y[23];
acadoWorkspace.Dy[24] -= acadoVariables.y[24];
acadoWorkspace.Dy[25] -= acadoVariables.y[25];
acadoWorkspace.Dy[26] -= acadoVariables.y[26];
acadoWorkspace.Dy[27] -= acadoVariables.y[27];
acadoWorkspace.Dy[28] -= acadoVariables.y[28];
acadoWorkspace.Dy[29] -= acadoVariables.y[29];
acadoWorkspace.Dy[30] -= acadoVariables.y[30];
acadoWorkspace.Dy[31] -= acadoVariables.y[31];
acadoWorkspace.Dy[32] -= acadoVariables.y[32];
acadoWorkspace.Dy[33] -= acadoVariables.y[33];
acadoWorkspace.Dy[34] -= acadoVariables.y[34];
acadoWorkspace.Dy[35] -= acadoVariables.y[35];
acadoWorkspace.Dy[36] -= acadoVariables.y[36];
acadoWorkspace.Dy[37] -= acadoVariables.y[37];
acadoWorkspace.Dy[38] -= acadoVariables.y[38];
acadoWorkspace.Dy[39] -= acadoVariables.y[39];
acadoWorkspace.Dy[40] -= acadoVariables.y[40];
acadoWorkspace.Dy[41] -= acadoVariables.y[41];
acadoWorkspace.Dy[42] -= acadoVariables.y[42];
acadoWorkspace.Dy[43] -= acadoVariables.y[43];
acadoWorkspace.Dy[44] -= acadoVariables.y[44];
acadoWorkspace.Dy[45] -= acadoVariables.y[45];
acadoWorkspace.Dy[46] -= acadoVariables.y[46];
acadoWorkspace.Dy[47] -= acadoVariables.y[47];
acadoWorkspace.Dy[48] -= acadoVariables.y[48];
acadoWorkspace.Dy[49] -= acadoVariables.y[49];
acadoWorkspace.Dy[50] -= acadoVariables.y[50];
acadoWorkspace.Dy[51] -= acadoVariables.y[51];
acadoWorkspace.Dy[52] -= acadoVariables.y[52];
acadoWorkspace.Dy[53] -= acadoVariables.y[53];
acadoWorkspace.Dy[54] -= acadoVariables.y[54];
acadoWorkspace.Dy[55] -= acadoVariables.y[55];
acadoWorkspace.Dy[56] -= acadoVariables.y[56];
acadoWorkspace.Dy[57] -= acadoVariables.y[57];
acadoWorkspace.Dy[58] -= acadoVariables.y[58];
acadoWorkspace.Dy[59] -= acadoVariables.y[59];
acadoWorkspace.Dy[60] -= acadoVariables.y[60];
acadoWorkspace.Dy[61] -= acadoVariables.y[61];
acadoWorkspace.Dy[62] -= acadoVariables.y[62];
acadoWorkspace.Dy[63] -= acadoVariables.y[63];
acadoWorkspace.Dy[64] -= acadoVariables.y[64];
acadoWorkspace.Dy[65] -= acadoVariables.y[65];
acadoWorkspace.Dy[66] -= acadoVariables.y[66];
acadoWorkspace.Dy[67] -= acadoVariables.y[67];
acadoWorkspace.Dy[68] -= acadoVariables.y[68];
acadoWorkspace.Dy[69] -= acadoVariables.y[69];
acadoWorkspace.Dy[70] -= acadoVariables.y[70];
acadoWorkspace.Dy[71] -= acadoVariables.y[71];
acadoWorkspace.Dy[72] -= acadoVariables.y[72];
acadoWorkspace.Dy[73] -= acadoVariables.y[73];
acadoWorkspace.Dy[74] -= acadoVariables.y[74];
acadoWorkspace.Dy[75] -= acadoVariables.y[75];
acadoWorkspace.Dy[76] -= acadoVariables.y[76];
acadoWorkspace.Dy[77] -= acadoVariables.y[77];
acadoWorkspace.Dy[78] -= acadoVariables.y[78];
acadoWorkspace.Dy[79] -= acadoVariables.y[79];
acadoWorkspace.Dy[80] -= acadoVariables.y[80];
acadoWorkspace.Dy[81] -= acadoVariables.y[81];
acadoWorkspace.Dy[82] -= acadoVariables.y[82];
acadoWorkspace.Dy[83] -= acadoVariables.y[83];
acadoWorkspace.Dy[84] -= acadoVariables.y[84];
acadoWorkspace.Dy[85] -= acadoVariables.y[85];
acadoWorkspace.Dy[86] -= acadoVariables.y[86];
acadoWorkspace.Dy[87] -= acadoVariables.y[87];
acadoWorkspace.Dy[88] -= acadoVariables.y[88];
acadoWorkspace.Dy[89] -= acadoVariables.y[89];
acadoWorkspace.DyN[0] -= acadoVariables.yN[0];
acadoWorkspace.DyN[1] -= acadoVariables.yN[1];
acadoWorkspace.DyN[2] -= acadoVariables.yN[2];
acadoWorkspace.DyN[3] -= acadoVariables.yN[3];
acadoWorkspace.DyN[4] -= acadoVariables.yN[4];
acadoWorkspace.DyN[5] -= acadoVariables.yN[5];

acado_multRDy( acadoWorkspace.R2, acadoWorkspace.Dy, acadoWorkspace.g );
acado_multRDy( &(acadoWorkspace.R2[ 27 ]), &(acadoWorkspace.Dy[ 9 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 54 ]), &(acadoWorkspace.Dy[ 18 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 81 ]), &(acadoWorkspace.Dy[ 27 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 108 ]), &(acadoWorkspace.Dy[ 36 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 135 ]), &(acadoWorkspace.Dy[ 45 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 162 ]), &(acadoWorkspace.Dy[ 54 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 189 ]), &(acadoWorkspace.Dy[ 63 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 216 ]), &(acadoWorkspace.Dy[ 72 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 243 ]), &(acadoWorkspace.Dy[ 81 ]), &(acadoWorkspace.g[ 27 ]) );

acado_multQDy( acadoWorkspace.Q2, acadoWorkspace.Dy, acadoWorkspace.QDy );
acado_multQDy( &(acadoWorkspace.Q2[ 54 ]), &(acadoWorkspace.Dy[ 9 ]), &(acadoWorkspace.QDy[ 6 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 108 ]), &(acadoWorkspace.Dy[ 18 ]), &(acadoWorkspace.QDy[ 12 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 162 ]), &(acadoWorkspace.Dy[ 27 ]), &(acadoWorkspace.QDy[ 18 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 216 ]), &(acadoWorkspace.Dy[ 36 ]), &(acadoWorkspace.QDy[ 24 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 270 ]), &(acadoWorkspace.Dy[ 45 ]), &(acadoWorkspace.QDy[ 30 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 324 ]), &(acadoWorkspace.Dy[ 54 ]), &(acadoWorkspace.QDy[ 36 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 378 ]), &(acadoWorkspace.Dy[ 63 ]), &(acadoWorkspace.QDy[ 42 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 432 ]), &(acadoWorkspace.Dy[ 72 ]), &(acadoWorkspace.QDy[ 48 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 486 ]), &(acadoWorkspace.Dy[ 81 ]), &(acadoWorkspace.QDy[ 54 ]) );

acadoWorkspace.QDy[60] = + acadoWorkspace.QN2[0]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[1]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[2]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[3]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[4]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[5]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[61] = + acadoWorkspace.QN2[6]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[7]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[8]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[9]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[10]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[11]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[62] = + acadoWorkspace.QN2[12]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[13]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[14]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[15]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[16]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[17]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[63] = + acadoWorkspace.QN2[18]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[19]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[20]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[21]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[22]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[23]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[64] = + acadoWorkspace.QN2[24]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[25]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[26]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[27]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[28]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[29]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[65] = + acadoWorkspace.QN2[30]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[31]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[32]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[33]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[34]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[35]*acadoWorkspace.DyN[5];

acadoWorkspace.sbar[0] = acadoWorkspace.Dx0[0];
acadoWorkspace.sbar[1] = acadoWorkspace.Dx0[1];
acadoWorkspace.sbar[2] = acadoWorkspace.Dx0[2];
acadoWorkspace.sbar[3] = acadoWorkspace.Dx0[3];
acadoWorkspace.sbar[4] = acadoWorkspace.Dx0[4];
acadoWorkspace.sbar[5] = acadoWorkspace.Dx0[5];
acado_macASbar( acadoWorkspace.evGx, acadoWorkspace.sbar, &(acadoWorkspace.sbar[ 6 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 36 ]), &(acadoWorkspace.sbar[ 6 ]), &(acadoWorkspace.sbar[ 12 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 72 ]), &(acadoWorkspace.sbar[ 12 ]), &(acadoWorkspace.sbar[ 18 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.sbar[ 18 ]), &(acadoWorkspace.sbar[ 24 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.sbar[ 24 ]), &(acadoWorkspace.sbar[ 30 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.sbar[ 30 ]), &(acadoWorkspace.sbar[ 36 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.sbar[ 36 ]), &(acadoWorkspace.sbar[ 42 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.sbar[ 42 ]), &(acadoWorkspace.sbar[ 48 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.sbar[ 48 ]), &(acadoWorkspace.sbar[ 54 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.sbar[ 54 ]), &(acadoWorkspace.sbar[ 60 ]) );

acadoWorkspace.w1[0] = + acadoWorkspace.QN1[0]*acadoWorkspace.sbar[60] + acadoWorkspace.QN1[1]*acadoWorkspace.sbar[61] + acadoWorkspace.QN1[2]*acadoWorkspace.sbar[62] + acadoWorkspace.QN1[3]*acadoWorkspace.sbar[63] + acadoWorkspace.QN1[4]*acadoWorkspace.sbar[64] + acadoWorkspace.QN1[5]*acadoWorkspace.sbar[65] + acadoWorkspace.QDy[60];
acadoWorkspace.w1[1] = + acadoWorkspace.QN1[6]*acadoWorkspace.sbar[60] + acadoWorkspace.QN1[7]*acadoWorkspace.sbar[61] + acadoWorkspace.QN1[8]*acadoWorkspace.sbar[62] + acadoWorkspace.QN1[9]*acadoWorkspace.sbar[63] + acadoWorkspace.QN1[10]*acadoWorkspace.sbar[64] + acadoWorkspace.QN1[11]*acadoWorkspace.sbar[65] + acadoWorkspace.QDy[61];
acadoWorkspace.w1[2] = + acadoWorkspace.QN1[12]*acadoWorkspace.sbar[60] + acadoWorkspace.QN1[13]*acadoWorkspace.sbar[61] + acadoWorkspace.QN1[14]*acadoWorkspace.sbar[62] + acadoWorkspace.QN1[15]*acadoWorkspace.sbar[63] + acadoWorkspace.QN1[16]*acadoWorkspace.sbar[64] + acadoWorkspace.QN1[17]*acadoWorkspace.sbar[65] + acadoWorkspace.QDy[62];
acadoWorkspace.w1[3] = + acadoWorkspace.QN1[18]*acadoWorkspace.sbar[60] + acadoWorkspace.QN1[19]*acadoWorkspace.sbar[61] + acadoWorkspace.QN1[20]*acadoWorkspace.sbar[62] + acadoWorkspace.QN1[21]*acadoWorkspace.sbar[63] + acadoWorkspace.QN1[22]*acadoWorkspace.sbar[64] + acadoWorkspace.QN1[23]*acadoWorkspace.sbar[65] + acadoWorkspace.QDy[63];
acadoWorkspace.w1[4] = + acadoWorkspace.QN1[24]*acadoWorkspace.sbar[60] + acadoWorkspace.QN1[25]*acadoWorkspace.sbar[61] + acadoWorkspace.QN1[26]*acadoWorkspace.sbar[62] + acadoWorkspace.QN1[27]*acadoWorkspace.sbar[63] + acadoWorkspace.QN1[28]*acadoWorkspace.sbar[64] + acadoWorkspace.QN1[29]*acadoWorkspace.sbar[65] + acadoWorkspace.QDy[64];
acadoWorkspace.w1[5] = + acadoWorkspace.QN1[30]*acadoWorkspace.sbar[60] + acadoWorkspace.QN1[31]*acadoWorkspace.sbar[61] + acadoWorkspace.QN1[32]*acadoWorkspace.sbar[62] + acadoWorkspace.QN1[33]*acadoWorkspace.sbar[63] + acadoWorkspace.QN1[34]*acadoWorkspace.sbar[64] + acadoWorkspace.QN1[35]*acadoWorkspace.sbar[65] + acadoWorkspace.QDy[65];
acado_macBTw1( &(acadoWorkspace.evGu[ 162 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 27 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 54 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.sbar[ 54 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 24 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 48 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.sbar[ 48 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 126 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 21 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 42 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.sbar[ 42 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 18 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 36 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.sbar[ 36 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 90 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 15 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 180 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 30 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.sbar[ 30 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 72 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 12 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 144 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 24 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.sbar[ 24 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 54 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 9 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 108 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 18 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 108 ]), &(acadoWorkspace.sbar[ 18 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 36 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 6 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 72 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 12 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 72 ]), &(acadoWorkspace.sbar[ 12 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 18 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 3 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 36 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 6 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 36 ]), &(acadoWorkspace.sbar[ 6 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( acadoWorkspace.evGu, acadoWorkspace.w1, acadoWorkspace.g );


}

void acado_expand(  )
{
acadoVariables.u[0] += acadoWorkspace.x[0];
acadoVariables.u[1] += acadoWorkspace.x[1];
acadoVariables.u[2] += acadoWorkspace.x[2];
acadoVariables.u[3] += acadoWorkspace.x[3];
acadoVariables.u[4] += acadoWorkspace.x[4];
acadoVariables.u[5] += acadoWorkspace.x[5];
acadoVariables.u[6] += acadoWorkspace.x[6];
acadoVariables.u[7] += acadoWorkspace.x[7];
acadoVariables.u[8] += acadoWorkspace.x[8];
acadoVariables.u[9] += acadoWorkspace.x[9];
acadoVariables.u[10] += acadoWorkspace.x[10];
acadoVariables.u[11] += acadoWorkspace.x[11];
acadoVariables.u[12] += acadoWorkspace.x[12];
acadoVariables.u[13] += acadoWorkspace.x[13];
acadoVariables.u[14] += acadoWorkspace.x[14];
acadoVariables.u[15] += acadoWorkspace.x[15];
acadoVariables.u[16] += acadoWorkspace.x[16];
acadoVariables.u[17] += acadoWorkspace.x[17];
acadoVariables.u[18] += acadoWorkspace.x[18];
acadoVariables.u[19] += acadoWorkspace.x[19];
acadoVariables.u[20] += acadoWorkspace.x[20];
acadoVariables.u[21] += acadoWorkspace.x[21];
acadoVariables.u[22] += acadoWorkspace.x[22];
acadoVariables.u[23] += acadoWorkspace.x[23];
acadoVariables.u[24] += acadoWorkspace.x[24];
acadoVariables.u[25] += acadoWorkspace.x[25];
acadoVariables.u[26] += acadoWorkspace.x[26];
acadoVariables.u[27] += acadoWorkspace.x[27];
acadoVariables.u[28] += acadoWorkspace.x[28];
acadoVariables.u[29] += acadoWorkspace.x[29];
acadoWorkspace.sbar[0] = acadoWorkspace.Dx0[0];
acadoWorkspace.sbar[1] = acadoWorkspace.Dx0[1];
acadoWorkspace.sbar[2] = acadoWorkspace.Dx0[2];
acadoWorkspace.sbar[3] = acadoWorkspace.Dx0[3];
acadoWorkspace.sbar[4] = acadoWorkspace.Dx0[4];
acadoWorkspace.sbar[5] = acadoWorkspace.Dx0[5];
acadoWorkspace.sbar[6] = acadoWorkspace.d[0];
acadoWorkspace.sbar[7] = acadoWorkspace.d[1];
acadoWorkspace.sbar[8] = acadoWorkspace.d[2];
acadoWorkspace.sbar[9] = acadoWorkspace.d[3];
acadoWorkspace.sbar[10] = acadoWorkspace.d[4];
acadoWorkspace.sbar[11] = acadoWorkspace.d[5];
acadoWorkspace.sbar[12] = acadoWorkspace.d[6];
acadoWorkspace.sbar[13] = acadoWorkspace.d[7];
acadoWorkspace.sbar[14] = acadoWorkspace.d[8];
acadoWorkspace.sbar[15] = acadoWorkspace.d[9];
acadoWorkspace.sbar[16] = acadoWorkspace.d[10];
acadoWorkspace.sbar[17] = acadoWorkspace.d[11];
acadoWorkspace.sbar[18] = acadoWorkspace.d[12];
acadoWorkspace.sbar[19] = acadoWorkspace.d[13];
acadoWorkspace.sbar[20] = acadoWorkspace.d[14];
acadoWorkspace.sbar[21] = acadoWorkspace.d[15];
acadoWorkspace.sbar[22] = acadoWorkspace.d[16];
acadoWorkspace.sbar[23] = acadoWorkspace.d[17];
acadoWorkspace.sbar[24] = acadoWorkspace.d[18];
acadoWorkspace.sbar[25] = acadoWorkspace.d[19];
acadoWorkspace.sbar[26] = acadoWorkspace.d[20];
acadoWorkspace.sbar[27] = acadoWorkspace.d[21];
acadoWorkspace.sbar[28] = acadoWorkspace.d[22];
acadoWorkspace.sbar[29] = acadoWorkspace.d[23];
acadoWorkspace.sbar[30] = acadoWorkspace.d[24];
acadoWorkspace.sbar[31] = acadoWorkspace.d[25];
acadoWorkspace.sbar[32] = acadoWorkspace.d[26];
acadoWorkspace.sbar[33] = acadoWorkspace.d[27];
acadoWorkspace.sbar[34] = acadoWorkspace.d[28];
acadoWorkspace.sbar[35] = acadoWorkspace.d[29];
acadoWorkspace.sbar[36] = acadoWorkspace.d[30];
acadoWorkspace.sbar[37] = acadoWorkspace.d[31];
acadoWorkspace.sbar[38] = acadoWorkspace.d[32];
acadoWorkspace.sbar[39] = acadoWorkspace.d[33];
acadoWorkspace.sbar[40] = acadoWorkspace.d[34];
acadoWorkspace.sbar[41] = acadoWorkspace.d[35];
acadoWorkspace.sbar[42] = acadoWorkspace.d[36];
acadoWorkspace.sbar[43] = acadoWorkspace.d[37];
acadoWorkspace.sbar[44] = acadoWorkspace.d[38];
acadoWorkspace.sbar[45] = acadoWorkspace.d[39];
acadoWorkspace.sbar[46] = acadoWorkspace.d[40];
acadoWorkspace.sbar[47] = acadoWorkspace.d[41];
acadoWorkspace.sbar[48] = acadoWorkspace.d[42];
acadoWorkspace.sbar[49] = acadoWorkspace.d[43];
acadoWorkspace.sbar[50] = acadoWorkspace.d[44];
acadoWorkspace.sbar[51] = acadoWorkspace.d[45];
acadoWorkspace.sbar[52] = acadoWorkspace.d[46];
acadoWorkspace.sbar[53] = acadoWorkspace.d[47];
acadoWorkspace.sbar[54] = acadoWorkspace.d[48];
acadoWorkspace.sbar[55] = acadoWorkspace.d[49];
acadoWorkspace.sbar[56] = acadoWorkspace.d[50];
acadoWorkspace.sbar[57] = acadoWorkspace.d[51];
acadoWorkspace.sbar[58] = acadoWorkspace.d[52];
acadoWorkspace.sbar[59] = acadoWorkspace.d[53];
acadoWorkspace.sbar[60] = acadoWorkspace.d[54];
acadoWorkspace.sbar[61] = acadoWorkspace.d[55];
acadoWorkspace.sbar[62] = acadoWorkspace.d[56];
acadoWorkspace.sbar[63] = acadoWorkspace.d[57];
acadoWorkspace.sbar[64] = acadoWorkspace.d[58];
acadoWorkspace.sbar[65] = acadoWorkspace.d[59];
acado_expansionStep( acadoWorkspace.evGx, acadoWorkspace.evGu, acadoWorkspace.x, acadoWorkspace.sbar, &(acadoWorkspace.sbar[ 6 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 36 ]), &(acadoWorkspace.evGu[ 18 ]), &(acadoWorkspace.x[ 3 ]), &(acadoWorkspace.sbar[ 6 ]), &(acadoWorkspace.sbar[ 12 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 72 ]), &(acadoWorkspace.evGu[ 36 ]), &(acadoWorkspace.x[ 6 ]), &(acadoWorkspace.sbar[ 12 ]), &(acadoWorkspace.sbar[ 18 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.evGu[ 54 ]), &(acadoWorkspace.x[ 9 ]), &(acadoWorkspace.sbar[ 18 ]), &(acadoWorkspace.sbar[ 24 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.evGu[ 72 ]), &(acadoWorkspace.x[ 12 ]), &(acadoWorkspace.sbar[ 24 ]), &(acadoWorkspace.sbar[ 30 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.evGu[ 90 ]), &(acadoWorkspace.x[ 15 ]), &(acadoWorkspace.sbar[ 30 ]), &(acadoWorkspace.sbar[ 36 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.evGu[ 108 ]), &(acadoWorkspace.x[ 18 ]), &(acadoWorkspace.sbar[ 36 ]), &(acadoWorkspace.sbar[ 42 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.evGu[ 126 ]), &(acadoWorkspace.x[ 21 ]), &(acadoWorkspace.sbar[ 42 ]), &(acadoWorkspace.sbar[ 48 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.evGu[ 144 ]), &(acadoWorkspace.x[ 24 ]), &(acadoWorkspace.sbar[ 48 ]), &(acadoWorkspace.sbar[ 54 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.evGu[ 162 ]), &(acadoWorkspace.x[ 27 ]), &(acadoWorkspace.sbar[ 54 ]), &(acadoWorkspace.sbar[ 60 ]) );
acadoVariables.x[0] += acadoWorkspace.sbar[0];
acadoVariables.x[1] += acadoWorkspace.sbar[1];
acadoVariables.x[2] += acadoWorkspace.sbar[2];
acadoVariables.x[3] += acadoWorkspace.sbar[3];
acadoVariables.x[4] += acadoWorkspace.sbar[4];
acadoVariables.x[5] += acadoWorkspace.sbar[5];
acadoVariables.x[6] += acadoWorkspace.sbar[6];
acadoVariables.x[7] += acadoWorkspace.sbar[7];
acadoVariables.x[8] += acadoWorkspace.sbar[8];
acadoVariables.x[9] += acadoWorkspace.sbar[9];
acadoVariables.x[10] += acadoWorkspace.sbar[10];
acadoVariables.x[11] += acadoWorkspace.sbar[11];
acadoVariables.x[12] += acadoWorkspace.sbar[12];
acadoVariables.x[13] += acadoWorkspace.sbar[13];
acadoVariables.x[14] += acadoWorkspace.sbar[14];
acadoVariables.x[15] += acadoWorkspace.sbar[15];
acadoVariables.x[16] += acadoWorkspace.sbar[16];
acadoVariables.x[17] += acadoWorkspace.sbar[17];
acadoVariables.x[18] += acadoWorkspace.sbar[18];
acadoVariables.x[19] += acadoWorkspace.sbar[19];
acadoVariables.x[20] += acadoWorkspace.sbar[20];
acadoVariables.x[21] += acadoWorkspace.sbar[21];
acadoVariables.x[22] += acadoWorkspace.sbar[22];
acadoVariables.x[23] += acadoWorkspace.sbar[23];
acadoVariables.x[24] += acadoWorkspace.sbar[24];
acadoVariables.x[25] += acadoWorkspace.sbar[25];
acadoVariables.x[26] += acadoWorkspace.sbar[26];
acadoVariables.x[27] += acadoWorkspace.sbar[27];
acadoVariables.x[28] += acadoWorkspace.sbar[28];
acadoVariables.x[29] += acadoWorkspace.sbar[29];
acadoVariables.x[30] += acadoWorkspace.sbar[30];
acadoVariables.x[31] += acadoWorkspace.sbar[31];
acadoVariables.x[32] += acadoWorkspace.sbar[32];
acadoVariables.x[33] += acadoWorkspace.sbar[33];
acadoVariables.x[34] += acadoWorkspace.sbar[34];
acadoVariables.x[35] += acadoWorkspace.sbar[35];
acadoVariables.x[36] += acadoWorkspace.sbar[36];
acadoVariables.x[37] += acadoWorkspace.sbar[37];
acadoVariables.x[38] += acadoWorkspace.sbar[38];
acadoVariables.x[39] += acadoWorkspace.sbar[39];
acadoVariables.x[40] += acadoWorkspace.sbar[40];
acadoVariables.x[41] += acadoWorkspace.sbar[41];
acadoVariables.x[42] += acadoWorkspace.sbar[42];
acadoVariables.x[43] += acadoWorkspace.sbar[43];
acadoVariables.x[44] += acadoWorkspace.sbar[44];
acadoVariables.x[45] += acadoWorkspace.sbar[45];
acadoVariables.x[46] += acadoWorkspace.sbar[46];
acadoVariables.x[47] += acadoWorkspace.sbar[47];
acadoVariables.x[48] += acadoWorkspace.sbar[48];
acadoVariables.x[49] += acadoWorkspace.sbar[49];
acadoVariables.x[50] += acadoWorkspace.sbar[50];
acadoVariables.x[51] += acadoWorkspace.sbar[51];
acadoVariables.x[52] += acadoWorkspace.sbar[52];
acadoVariables.x[53] += acadoWorkspace.sbar[53];
acadoVariables.x[54] += acadoWorkspace.sbar[54];
acadoVariables.x[55] += acadoWorkspace.sbar[55];
acadoVariables.x[56] += acadoWorkspace.sbar[56];
acadoVariables.x[57] += acadoWorkspace.sbar[57];
acadoVariables.x[58] += acadoWorkspace.sbar[58];
acadoVariables.x[59] += acadoWorkspace.sbar[59];
acadoVariables.x[60] += acadoWorkspace.sbar[60];
acadoVariables.x[61] += acadoWorkspace.sbar[61];
acadoVariables.x[62] += acadoWorkspace.sbar[62];
acadoVariables.x[63] += acadoWorkspace.sbar[63];
acadoVariables.x[64] += acadoWorkspace.sbar[64];
acadoVariables.x[65] += acadoWorkspace.sbar[65];
}

int acado_preparationStep(  )
{
int ret;

ret = acado_modelSimulation();
acado_evaluateObjective(  );
acado_condensePrep(  );
return ret;
}

int acado_feedbackStep(  )
{
int tmp;

acado_condenseFdb(  );

tmp = acado_solve( );

acado_expand(  );
return tmp;
}

int acado_initializeSolver(  )
{
int ret;

/* This is a function which must be called once before any other function call! */


ret = 0;

memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
return ret;
}

void acado_initializeNodesByForwardSimulation(  )
{
int index;
for (index = 0; index < 10; ++index)
{
acadoWorkspace.state[0] = acadoVariables.x[index * 6];
acadoWorkspace.state[1] = acadoVariables.x[index * 6 + 1];
acadoWorkspace.state[2] = acadoVariables.x[index * 6 + 2];
acadoWorkspace.state[3] = acadoVariables.x[index * 6 + 3];
acadoWorkspace.state[4] = acadoVariables.x[index * 6 + 4];
acadoWorkspace.state[5] = acadoVariables.x[index * 6 + 5];
acadoWorkspace.state[60] = acadoVariables.u[index * 3];
acadoWorkspace.state[61] = acadoVariables.u[index * 3 + 1];
acadoWorkspace.state[62] = acadoVariables.u[index * 3 + 2];
acadoWorkspace.state[63] = acadoVariables.od[index * 4];
acadoWorkspace.state[64] = acadoVariables.od[index * 4 + 1];
acadoWorkspace.state[65] = acadoVariables.od[index * 4 + 2];
acadoWorkspace.state[66] = acadoVariables.od[index * 4 + 3];

acado_integrate(acadoWorkspace.state, index == 0);

acadoVariables.x[index * 6 + 6] = acadoWorkspace.state[0];
acadoVariables.x[index * 6 + 7] = acadoWorkspace.state[1];
acadoVariables.x[index * 6 + 8] = acadoWorkspace.state[2];
acadoVariables.x[index * 6 + 9] = acadoWorkspace.state[3];
acadoVariables.x[index * 6 + 10] = acadoWorkspace.state[4];
acadoVariables.x[index * 6 + 11] = acadoWorkspace.state[5];
}
}

void acado_shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd )
{
int index;
for (index = 0; index < 10; ++index)
{
acadoVariables.x[index * 6] = acadoVariables.x[index * 6 + 6];
acadoVariables.x[index * 6 + 1] = acadoVariables.x[index * 6 + 7];
acadoVariables.x[index * 6 + 2] = acadoVariables.x[index * 6 + 8];
acadoVariables.x[index * 6 + 3] = acadoVariables.x[index * 6 + 9];
acadoVariables.x[index * 6 + 4] = acadoVariables.x[index * 6 + 10];
acadoVariables.x[index * 6 + 5] = acadoVariables.x[index * 6 + 11];
}

if (strategy == 1 && xEnd != 0)
{
acadoVariables.x[60] = xEnd[0];
acadoVariables.x[61] = xEnd[1];
acadoVariables.x[62] = xEnd[2];
acadoVariables.x[63] = xEnd[3];
acadoVariables.x[64] = xEnd[4];
acadoVariables.x[65] = xEnd[5];
}
else if (strategy == 2) 
{
acadoWorkspace.state[0] = acadoVariables.x[60];
acadoWorkspace.state[1] = acadoVariables.x[61];
acadoWorkspace.state[2] = acadoVariables.x[62];
acadoWorkspace.state[3] = acadoVariables.x[63];
acadoWorkspace.state[4] = acadoVariables.x[64];
acadoWorkspace.state[5] = acadoVariables.x[65];
if (uEnd != 0)
{
acadoWorkspace.state[60] = uEnd[0];
acadoWorkspace.state[61] = uEnd[1];
acadoWorkspace.state[62] = uEnd[2];
}
else
{
acadoWorkspace.state[60] = acadoVariables.u[27];
acadoWorkspace.state[61] = acadoVariables.u[28];
acadoWorkspace.state[62] = acadoVariables.u[29];
}
acadoWorkspace.state[63] = acadoVariables.od[40];
acadoWorkspace.state[64] = acadoVariables.od[41];
acadoWorkspace.state[65] = acadoVariables.od[42];
acadoWorkspace.state[66] = acadoVariables.od[43];

acado_integrate(acadoWorkspace.state, 1);

acadoVariables.x[60] = acadoWorkspace.state[0];
acadoVariables.x[61] = acadoWorkspace.state[1];
acadoVariables.x[62] = acadoWorkspace.state[2];
acadoVariables.x[63] = acadoWorkspace.state[3];
acadoVariables.x[64] = acadoWorkspace.state[4];
acadoVariables.x[65] = acadoWorkspace.state[5];
}
}

void acado_shiftControls( real_t* const uEnd )
{
int index;
for (index = 0; index < 9; ++index)
{
acadoVariables.u[index * 3] = acadoVariables.u[index * 3 + 3];
acadoVariables.u[index * 3 + 1] = acadoVariables.u[index * 3 + 4];
acadoVariables.u[index * 3 + 2] = acadoVariables.u[index * 3 + 5];
}

if (uEnd != 0)
{
acadoVariables.u[27] = uEnd[0];
acadoVariables.u[28] = uEnd[1];
acadoVariables.u[29] = uEnd[2];
}
}

real_t acado_getKKT(  )
{
real_t kkt;

int index;
real_t prd;

kkt = + acadoWorkspace.g[0]*acadoWorkspace.x[0] + acadoWorkspace.g[1]*acadoWorkspace.x[1] + acadoWorkspace.g[2]*acadoWorkspace.x[2] + acadoWorkspace.g[3]*acadoWorkspace.x[3] + acadoWorkspace.g[4]*acadoWorkspace.x[4] + acadoWorkspace.g[5]*acadoWorkspace.x[5] + acadoWorkspace.g[6]*acadoWorkspace.x[6] + acadoWorkspace.g[7]*acadoWorkspace.x[7] + acadoWorkspace.g[8]*acadoWorkspace.x[8] + acadoWorkspace.g[9]*acadoWorkspace.x[9] + acadoWorkspace.g[10]*acadoWorkspace.x[10] + acadoWorkspace.g[11]*acadoWorkspace.x[11] + acadoWorkspace.g[12]*acadoWorkspace.x[12] + acadoWorkspace.g[13]*acadoWorkspace.x[13] + acadoWorkspace.g[14]*acadoWorkspace.x[14] + acadoWorkspace.g[15]*acadoWorkspace.x[15] + acadoWorkspace.g[16]*acadoWorkspace.x[16] + acadoWorkspace.g[17]*acadoWorkspace.x[17] + acadoWorkspace.g[18]*acadoWorkspace.x[18] + acadoWorkspace.g[19]*acadoWorkspace.x[19] + acadoWorkspace.g[20]*acadoWorkspace.x[20] + acadoWorkspace.g[21]*acadoWorkspace.x[21] + acadoWorkspace.g[22]*acadoWorkspace.x[22] + acadoWorkspace.g[23]*acadoWorkspace.x[23] + acadoWorkspace.g[24]*acadoWorkspace.x[24] + acadoWorkspace.g[25]*acadoWorkspace.x[25] + acadoWorkspace.g[26]*acadoWorkspace.x[26] + acadoWorkspace.g[27]*acadoWorkspace.x[27] + acadoWorkspace.g[28]*acadoWorkspace.x[28] + acadoWorkspace.g[29]*acadoWorkspace.x[29];
kkt = fabs( kkt );
for (index = 0; index < 30; ++index)
{
prd = acadoWorkspace.y[index];
if (prd > 1e-12)
kkt += fabs(acadoWorkspace.lb[index] * prd);
else if (prd < -1e-12)
kkt += fabs(acadoWorkspace.ub[index] * prd);
}
return kkt;
}

real_t acado_getObjective(  )
{
real_t objVal;

int lRun1;
/** Row vector of size: 9 */
real_t tmpDy[ 9 ];

/** Row vector of size: 6 */
real_t tmpDyN[ 6 ];

for (lRun1 = 0; lRun1 < 10; ++lRun1)
{
acadoWorkspace.objValueIn[0] = acadoVariables.x[lRun1 * 6];
acadoWorkspace.objValueIn[1] = acadoVariables.x[lRun1 * 6 + 1];
acadoWorkspace.objValueIn[2] = acadoVariables.x[lRun1 * 6 + 2];
acadoWorkspace.objValueIn[3] = acadoVariables.x[lRun1 * 6 + 3];
acadoWorkspace.objValueIn[4] = acadoVariables.x[lRun1 * 6 + 4];
acadoWorkspace.objValueIn[5] = acadoVariables.x[lRun1 * 6 + 5];
acadoWorkspace.objValueIn[6] = acadoVariables.u[lRun1 * 3];
acadoWorkspace.objValueIn[7] = acadoVariables.u[lRun1 * 3 + 1];
acadoWorkspace.objValueIn[8] = acadoVariables.u[lRun1 * 3 + 2];
acadoWorkspace.objValueIn[9] = acadoVariables.od[lRun1 * 4];
acadoWorkspace.objValueIn[10] = acadoVariables.od[lRun1 * 4 + 1];
acadoWorkspace.objValueIn[11] = acadoVariables.od[lRun1 * 4 + 2];
acadoWorkspace.objValueIn[12] = acadoVariables.od[lRun1 * 4 + 3];

acado_evaluateLSQ( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.Dy[lRun1 * 9] = acadoWorkspace.objValueOut[0] - acadoVariables.y[lRun1 * 9];
acadoWorkspace.Dy[lRun1 * 9 + 1] = acadoWorkspace.objValueOut[1] - acadoVariables.y[lRun1 * 9 + 1];
acadoWorkspace.Dy[lRun1 * 9 + 2] = acadoWorkspace.objValueOut[2] - acadoVariables.y[lRun1 * 9 + 2];
acadoWorkspace.Dy[lRun1 * 9 + 3] = acadoWorkspace.objValueOut[3] - acadoVariables.y[lRun1 * 9 + 3];
acadoWorkspace.Dy[lRun1 * 9 + 4] = acadoWorkspace.objValueOut[4] - acadoVariables.y[lRun1 * 9 + 4];
acadoWorkspace.Dy[lRun1 * 9 + 5] = acadoWorkspace.objValueOut[5] - acadoVariables.y[lRun1 * 9 + 5];
acadoWorkspace.Dy[lRun1 * 9 + 6] = acadoWorkspace.objValueOut[6] - acadoVariables.y[lRun1 * 9 + 6];
acadoWorkspace.Dy[lRun1 * 9 + 7] = acadoWorkspace.objValueOut[7] - acadoVariables.y[lRun1 * 9 + 7];
acadoWorkspace.Dy[lRun1 * 9 + 8] = acadoWorkspace.objValueOut[8] - acadoVariables.y[lRun1 * 9 + 8];
}
acadoWorkspace.objValueIn[0] = acadoVariables.x[60];
acadoWorkspace.objValueIn[1] = acadoVariables.x[61];
acadoWorkspace.objValueIn[2] = acadoVariables.x[62];
acadoWorkspace.objValueIn[3] = acadoVariables.x[63];
acadoWorkspace.objValueIn[4] = acadoVariables.x[64];
acadoWorkspace.objValueIn[5] = acadoVariables.x[65];
acadoWorkspace.objValueIn[6] = acadoVariables.od[40];
acadoWorkspace.objValueIn[7] = acadoVariables.od[41];
acadoWorkspace.objValueIn[8] = acadoVariables.od[42];
acadoWorkspace.objValueIn[9] = acadoVariables.od[43];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0] - acadoVariables.yN[0];
acadoWorkspace.DyN[1] = acadoWorkspace.objValueOut[1] - acadoVariables.yN[1];
acadoWorkspace.DyN[2] = acadoWorkspace.objValueOut[2] - acadoVariables.yN[2];
acadoWorkspace.DyN[3] = acadoWorkspace.objValueOut[3] - acadoVariables.yN[3];
acadoWorkspace.DyN[4] = acadoWorkspace.objValueOut[4] - acadoVariables.yN[4];
acadoWorkspace.DyN[5] = acadoWorkspace.objValueOut[5] - acadoVariables.yN[5];
objVal = 0.0000000000000000e+00;
for (lRun1 = 0; lRun1 < 10; ++lRun1)
{
tmpDy[0] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[0];
tmpDy[1] = + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[10];
tmpDy[2] = + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[20];
tmpDy[3] = + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[30];
tmpDy[4] = + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[40];
tmpDy[5] = + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[50];
tmpDy[6] = + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[60];
tmpDy[7] = + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[70];
tmpDy[8] = + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[80];
objVal += + acadoWorkspace.Dy[lRun1 * 9]*tmpDy[0] + acadoWorkspace.Dy[lRun1 * 9 + 1]*tmpDy[1] + acadoWorkspace.Dy[lRun1 * 9 + 2]*tmpDy[2] + acadoWorkspace.Dy[lRun1 * 9 + 3]*tmpDy[3] + acadoWorkspace.Dy[lRun1 * 9 + 4]*tmpDy[4] + acadoWorkspace.Dy[lRun1 * 9 + 5]*tmpDy[5] + acadoWorkspace.Dy[lRun1 * 9 + 6]*tmpDy[6] + acadoWorkspace.Dy[lRun1 * 9 + 7]*tmpDy[7] + acadoWorkspace.Dy[lRun1 * 9 + 8]*tmpDy[8];
}

tmpDyN[0] = + acadoWorkspace.DyN[0]*acadoVariables.WN[0];
tmpDyN[1] = + acadoWorkspace.DyN[1]*acadoVariables.WN[7];
tmpDyN[2] = + acadoWorkspace.DyN[2]*acadoVariables.WN[14];
tmpDyN[3] = + acadoWorkspace.DyN[3]*acadoVariables.WN[21];
tmpDyN[4] = + acadoWorkspace.DyN[4]*acadoVariables.WN[28];
tmpDyN[5] = + acadoWorkspace.DyN[5]*acadoVariables.WN[35];
objVal += + acadoWorkspace.DyN[0]*tmpDyN[0] + acadoWorkspace.DyN[1]*tmpDyN[1] + acadoWorkspace.DyN[2]*tmpDyN[2] + acadoWorkspace.DyN[3]*tmpDyN[3] + acadoWorkspace.DyN[4]*tmpDyN[4] + acadoWorkspace.DyN[5]*tmpDyN[5];

objVal *= 0.5;
return objVal;
}

