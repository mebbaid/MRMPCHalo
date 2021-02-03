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


void acado_rhs(const real_t* in, real_t* out)
{
const real_t* xd = in;
const real_t* u = in + 6;
/* Vector of auxiliary variables; number of elements: 30. */
real_t* a = acadoWorkspace.rhs_aux;

/* Compute intermediate quantities: */
a[0] = (pow(((real_t)(-9.8785100000000003e-01)+xd[0]),2));
a[1] = ((xd[1])*(xd[1]));
a[2] = ((xd[2])*(xd[2]));
a[3] = (sqrt(((a[0]+a[1])+a[2])));
a[4] = (pow(a[3],3));
a[5] = (pow(((real_t)(1.2298413424696638e-02)+xd[0]),2));
a[6] = ((xd[1])*(xd[1]));
a[7] = ((xd[2])*(xd[2]));
a[8] = (sqrt(((a[5]+a[6])+a[7])));
a[9] = (pow(a[8],3));
a[10] = (pow(((real_t)(-9.8785100000000003e-01)+xd[0]),2));
a[11] = ((xd[1])*(xd[1]));
a[12] = ((xd[2])*(xd[2]));
a[13] = (sqrt(((a[10]+a[11])+a[12])));
a[14] = (pow(a[13],3));
a[15] = (pow(((real_t)(1.2298413424696638e-02)+xd[0]),2));
a[16] = ((xd[1])*(xd[1]));
a[17] = ((xd[2])*(xd[2]));
a[18] = (sqrt(((a[15]+a[16])+a[17])));
a[19] = (pow(a[18],3));
a[20] = (pow(((real_t)(1.2298413424696638e-02)+xd[0]),2));
a[21] = ((xd[1])*(xd[1]));
a[22] = ((xd[2])*(xd[2]));
a[23] = (sqrt(((a[20]+a[21])+a[22])));
a[24] = (pow(a[23],3));
a[25] = (pow(((real_t)(-9.8785100000000003e-01)+xd[0]),2));
a[26] = ((xd[1])*(xd[1]));
a[27] = ((xd[2])*(xd[2]));
a[28] = (sqrt(((a[25]+a[26])+a[27])));
a[29] = (pow(a[28],3));

/* Compute outputs: */
out[0] = xd[3];
out[1] = xd[4];
out[2] = xd[5];
out[3] = ((((((((real_t)(0.0000000000000000e+00)-((real_t)(-9.8785100000000003e-01)+xd[0]))*(real_t)(1.2149000000000000e-02))/a[4])-((((real_t)(1.2149000000000000e-02)+xd[0])*(real_t)(9.8785100000000003e-01))/a[9]))+((real_t)(2.0000000000000000e+00)*xd[4]))+u[0])+xd[0]);
out[4] = ((((((real_t)(-2.0000000000000000e+00)*xd[3])-(((real_t)(1.2149000000000000e-02)/a[14])*xd[1]))-(((real_t)(9.8785100000000003e-01)/a[19])*xd[1]))+u[1])+xd[1]);
out[5] = (((((real_t)(-9.8785100000000003e-01)/a[24])*xd[2])-(((real_t)(1.2149000000000000e-02)/a[29])*xd[2]))+u[2]);
}



void acado_diffs(const real_t* in, real_t* out)
{
const real_t* xd = in;
/* Vector of auxiliary variables; number of elements: 114. */
real_t* a = acadoWorkspace.rhs_aux;

/* Compute intermediate quantities: */
a[0] = (pow(((real_t)(-9.8785100000000003e-01)+xd[0]),2));
a[1] = ((xd[1])*(xd[1]));
a[2] = ((xd[2])*(xd[2]));
a[3] = (sqrt(((a[0]+a[1])+a[2])));
a[4] = (pow(a[3],3));
a[5] = ((real_t)(1.0000000000000000e+00)/a[4]);
a[6] = ((real_t)(3.0000000000000000e+00)*((a[3])*(a[3])));
a[7] = ((real_t)(2.0000000000000000e+00)*((real_t)(-9.8785100000000003e-01)+xd[0]));
a[8] = (1.0/sqrt(((a[0]+a[1])+a[2])));
a[9] = (a[8]*(real_t)(5.0000000000000000e-01));
a[10] = (a[7]*a[9]);
a[11] = (a[6]*a[10]);
a[12] = (a[5]*a[5]);
a[13] = (pow(((real_t)(1.2298413424696638e-02)+xd[0]),2));
a[14] = ((xd[1])*(xd[1]));
a[15] = ((xd[2])*(xd[2]));
a[16] = (sqrt(((a[13]+a[14])+a[15])));
a[17] = (pow(a[16],3));
a[18] = ((real_t)(1.0000000000000000e+00)/a[17]);
a[19] = ((real_t)(3.0000000000000000e+00)*((a[16])*(a[16])));
a[20] = ((real_t)(2.0000000000000000e+00)*((real_t)(1.2298413424696638e-02)+xd[0]));
a[21] = (1.0/sqrt(((a[13]+a[14])+a[15])));
a[22] = (a[21]*(real_t)(5.0000000000000000e-01));
a[23] = (a[20]*a[22]);
a[24] = (a[19]*a[23]);
a[25] = (a[18]*a[18]);
a[26] = ((real_t)(2.0000000000000000e+00)*xd[1]);
a[27] = (a[26]*a[9]);
a[28] = (a[6]*a[27]);
a[29] = ((real_t)(2.0000000000000000e+00)*xd[1]);
a[30] = (a[29]*a[22]);
a[31] = (a[19]*a[30]);
a[32] = ((real_t)(2.0000000000000000e+00)*xd[2]);
a[33] = (a[32]*a[9]);
a[34] = (a[6]*a[33]);
a[35] = ((real_t)(2.0000000000000000e+00)*xd[2]);
a[36] = (a[35]*a[22]);
a[37] = (a[19]*a[36]);
a[38] = (pow(((real_t)(-9.8785100000000003e-01)+xd[0]),2));
a[39] = ((xd[1])*(xd[1]));
a[40] = ((xd[2])*(xd[2]));
a[41] = (sqrt(((a[38]+a[39])+a[40])));
a[42] = ((real_t)(3.0000000000000000e+00)*((a[41])*(a[41])));
a[43] = ((real_t)(2.0000000000000000e+00)*((real_t)(-9.8785100000000003e-01)+xd[0]));
a[44] = (1.0/sqrt(((a[38]+a[39])+a[40])));
a[45] = (a[44]*(real_t)(5.0000000000000000e-01));
a[46] = (a[43]*a[45]);
a[47] = (a[42]*a[46]);
a[48] = (pow(a[41],3));
a[49] = ((real_t)(1.0000000000000000e+00)/a[48]);
a[50] = (a[49]*a[49]);
a[51] = (pow(((real_t)(1.2298413424696638e-02)+xd[0]),2));
a[52] = ((xd[1])*(xd[1]));
a[53] = ((xd[2])*(xd[2]));
a[54] = (sqrt(((a[51]+a[52])+a[53])));
a[55] = ((real_t)(3.0000000000000000e+00)*((a[54])*(a[54])));
a[56] = ((real_t)(2.0000000000000000e+00)*((real_t)(1.2298413424696638e-02)+xd[0]));
a[57] = (1.0/sqrt(((a[51]+a[52])+a[53])));
a[58] = (a[57]*(real_t)(5.0000000000000000e-01));
a[59] = (a[56]*a[58]);
a[60] = (a[55]*a[59]);
a[61] = (pow(a[54],3));
a[62] = ((real_t)(1.0000000000000000e+00)/a[61]);
a[63] = (a[62]*a[62]);
a[64] = ((real_t)(2.0000000000000000e+00)*xd[1]);
a[65] = (a[64]*a[45]);
a[66] = (a[42]*a[65]);
a[67] = ((real_t)(2.0000000000000000e+00)*xd[1]);
a[68] = (a[67]*a[58]);
a[69] = (a[55]*a[68]);
a[70] = ((real_t)(2.0000000000000000e+00)*xd[2]);
a[71] = (a[70]*a[45]);
a[72] = (a[42]*a[71]);
a[73] = ((real_t)(2.0000000000000000e+00)*xd[2]);
a[74] = (a[73]*a[58]);
a[75] = (a[55]*a[74]);
a[76] = (pow(((real_t)(1.2298413424696638e-02)+xd[0]),2));
a[77] = ((xd[1])*(xd[1]));
a[78] = ((xd[2])*(xd[2]));
a[79] = (sqrt(((a[76]+a[77])+a[78])));
a[80] = ((real_t)(3.0000000000000000e+00)*((a[79])*(a[79])));
a[81] = ((real_t)(2.0000000000000000e+00)*((real_t)(1.2298413424696638e-02)+xd[0]));
a[82] = (1.0/sqrt(((a[76]+a[77])+a[78])));
a[83] = (a[82]*(real_t)(5.0000000000000000e-01));
a[84] = (a[81]*a[83]);
a[85] = (a[80]*a[84]);
a[86] = (pow(a[79],3));
a[87] = ((real_t)(1.0000000000000000e+00)/a[86]);
a[88] = (a[87]*a[87]);
a[89] = (pow(((real_t)(-9.8785100000000003e-01)+xd[0]),2));
a[90] = ((xd[1])*(xd[1]));
a[91] = ((xd[2])*(xd[2]));
a[92] = (sqrt(((a[89]+a[90])+a[91])));
a[93] = ((real_t)(3.0000000000000000e+00)*((a[92])*(a[92])));
a[94] = ((real_t)(2.0000000000000000e+00)*((real_t)(-9.8785100000000003e-01)+xd[0]));
a[95] = (1.0/sqrt(((a[89]+a[90])+a[91])));
a[96] = (a[95]*(real_t)(5.0000000000000000e-01));
a[97] = (a[94]*a[96]);
a[98] = (a[93]*a[97]);
a[99] = (pow(a[92],3));
a[100] = ((real_t)(1.0000000000000000e+00)/a[99]);
a[101] = (a[100]*a[100]);
a[102] = ((real_t)(2.0000000000000000e+00)*xd[1]);
a[103] = (a[102]*a[83]);
a[104] = (a[80]*a[103]);
a[105] = ((real_t)(2.0000000000000000e+00)*xd[1]);
a[106] = (a[105]*a[96]);
a[107] = (a[93]*a[106]);
a[108] = ((real_t)(2.0000000000000000e+00)*xd[2]);
a[109] = (a[108]*a[83]);
a[110] = (a[80]*a[109]);
a[111] = ((real_t)(2.0000000000000000e+00)*xd[2]);
a[112] = (a[111]*a[96]);
a[113] = (a[93]*a[112]);

/* Compute outputs: */
out[0] = (real_t)(0.0000000000000000e+00);
out[1] = (real_t)(0.0000000000000000e+00);
out[2] = (real_t)(0.0000000000000000e+00);
out[3] = (real_t)(1.0000000000000000e+00);
out[4] = (real_t)(0.0000000000000000e+00);
out[5] = (real_t)(0.0000000000000000e+00);
out[6] = (real_t)(0.0000000000000000e+00);
out[7] = (real_t)(0.0000000000000000e+00);
out[8] = (real_t)(0.0000000000000000e+00);
out[9] = (real_t)(0.0000000000000000e+00);
out[10] = (real_t)(0.0000000000000000e+00);
out[11] = (real_t)(0.0000000000000000e+00);
out[12] = (real_t)(0.0000000000000000e+00);
out[13] = (real_t)(1.0000000000000000e+00);
out[14] = (real_t)(0.0000000000000000e+00);
out[15] = (real_t)(0.0000000000000000e+00);
out[16] = (real_t)(0.0000000000000000e+00);
out[17] = (real_t)(0.0000000000000000e+00);
out[18] = (real_t)(0.0000000000000000e+00);
out[19] = (real_t)(0.0000000000000000e+00);
out[20] = (real_t)(0.0000000000000000e+00);
out[21] = (real_t)(0.0000000000000000e+00);
out[22] = (real_t)(0.0000000000000000e+00);
out[23] = (real_t)(1.0000000000000000e+00);
out[24] = (real_t)(0.0000000000000000e+00);
out[25] = (real_t)(0.0000000000000000e+00);
out[26] = (real_t)(0.0000000000000000e+00);
out[27] = ((((((real_t)(-1.2149000000000000e-02))*a[5])-(((((real_t)(0.0000000000000000e+00)-((real_t)(-9.8785100000000003e-01)+xd[0]))*(real_t)(1.2149000000000000e-02))*a[11])*a[12]))-(((real_t)(9.8785100000000003e-01)*a[18])-(((((real_t)(1.2149000000000000e-02)+xd[0])*(real_t)(9.8785100000000003e-01))*a[24])*a[25])))+(real_t)(1.0000000000000000e+00));
out[28] = (((real_t)(0.0000000000000000e+00)-(((((real_t)(0.0000000000000000e+00)-((real_t)(-9.8785100000000003e-01)+xd[0]))*(real_t)(1.2149000000000000e-02))*a[28])*a[12]))-((real_t)(0.0000000000000000e+00)-(((((real_t)(1.2149000000000000e-02)+xd[0])*(real_t)(9.8785100000000003e-01))*a[31])*a[25])));
out[29] = (((real_t)(0.0000000000000000e+00)-(((((real_t)(0.0000000000000000e+00)-((real_t)(-9.8785100000000003e-01)+xd[0]))*(real_t)(1.2149000000000000e-02))*a[34])*a[12]))-((real_t)(0.0000000000000000e+00)-(((((real_t)(1.2149000000000000e-02)+xd[0])*(real_t)(9.8785100000000003e-01))*a[37])*a[25])));
out[30] = (real_t)(0.0000000000000000e+00);
out[31] = (real_t)(2.0000000000000000e+00);
out[32] = (real_t)(0.0000000000000000e+00);
out[33] = (real_t)(1.0000000000000000e+00);
out[34] = (real_t)(0.0000000000000000e+00);
out[35] = (real_t)(0.0000000000000000e+00);
out[36] = (((real_t)(0.0000000000000000e+00)-(((real_t)(0.0000000000000000e+00)-(((real_t)(1.2149000000000000e-02)*a[47])*a[50]))*xd[1]))-(((real_t)(0.0000000000000000e+00)-(((real_t)(9.8785100000000003e-01)*a[60])*a[63]))*xd[1]));
out[37] = ((((real_t)(0.0000000000000000e+00)-((((real_t)(0.0000000000000000e+00)-(((real_t)(1.2149000000000000e-02)*a[66])*a[50]))*xd[1])+((real_t)(1.2149000000000000e-02)/a[48])))-((((real_t)(0.0000000000000000e+00)-(((real_t)(9.8785100000000003e-01)*a[69])*a[63]))*xd[1])+((real_t)(9.8785100000000003e-01)/a[61])))+(real_t)(1.0000000000000000e+00));
out[38] = (((real_t)(0.0000000000000000e+00)-(((real_t)(0.0000000000000000e+00)-(((real_t)(1.2149000000000000e-02)*a[72])*a[50]))*xd[1]))-(((real_t)(0.0000000000000000e+00)-(((real_t)(9.8785100000000003e-01)*a[75])*a[63]))*xd[1]));
out[39] = (real_t)(-2.0000000000000000e+00);
out[40] = (real_t)(0.0000000000000000e+00);
out[41] = (real_t)(0.0000000000000000e+00);
out[42] = (real_t)(0.0000000000000000e+00);
out[43] = (real_t)(1.0000000000000000e+00);
out[44] = (real_t)(0.0000000000000000e+00);
out[45] = ((((real_t)(0.0000000000000000e+00)-(((real_t)(-9.8785100000000003e-01)*a[85])*a[88]))*xd[2])-(((real_t)(0.0000000000000000e+00)-(((real_t)(1.2149000000000000e-02)*a[98])*a[101]))*xd[2]));
out[46] = ((((real_t)(0.0000000000000000e+00)-(((real_t)(-9.8785100000000003e-01)*a[104])*a[88]))*xd[2])-(((real_t)(0.0000000000000000e+00)-(((real_t)(1.2149000000000000e-02)*a[107])*a[101]))*xd[2]));
out[47] = (((((real_t)(0.0000000000000000e+00)-(((real_t)(-9.8785100000000003e-01)*a[110])*a[88]))*xd[2])+((real_t)(-9.8785100000000003e-01)/a[86]))-((((real_t)(0.0000000000000000e+00)-(((real_t)(1.2149000000000000e-02)*a[113])*a[101]))*xd[2])+((real_t)(1.2149000000000000e-02)/a[99])));
out[48] = (real_t)(0.0000000000000000e+00);
out[49] = (real_t)(0.0000000000000000e+00);
out[50] = (real_t)(0.0000000000000000e+00);
out[51] = (real_t)(0.0000000000000000e+00);
out[52] = (real_t)(0.0000000000000000e+00);
out[53] = (real_t)(1.0000000000000000e+00);
}



void acado_solve_dim6_triangular( real_t* const A, real_t* const b )
{

b[5] = b[5]/A[35];
b[4] -= + A[29]*b[5];
b[4] = b[4]/A[28];
b[3] -= + A[23]*b[5];
b[3] -= + A[22]*b[4];
b[3] = b[3]/A[21];
b[2] -= + A[17]*b[5];
b[2] -= + A[16]*b[4];
b[2] -= + A[15]*b[3];
b[2] = b[2]/A[14];
b[1] -= + A[11]*b[5];
b[1] -= + A[10]*b[4];
b[1] -= + A[9]*b[3];
b[1] -= + A[8]*b[2];
b[1] = b[1]/A[7];
b[0] -= + A[5]*b[5];
b[0] -= + A[4]*b[4];
b[0] -= + A[3]*b[3];
b[0] -= + A[2]*b[2];
b[0] -= + A[1]*b[1];
b[0] = b[0]/A[0];
}

real_t acado_solve_dim6_system( real_t* const A, real_t* const b, int* const rk_perm )
{
real_t det;

int i;
int j;
int k;

int indexMax;

int intSwap;

real_t valueMax;

real_t temp;

for (i = 0; i < 6; ++i)
{
rk_perm[i] = i;
}
det = 1.0000000000000000e+00;
for( i=0; i < (5); i++ ) {
	indexMax = i;
	valueMax = fabs(A[i*6+i]);
	for( j=(i+1); j < 6; j++ ) {
		temp = fabs(A[j*6+i]);
		if( temp > valueMax ) {
			indexMax = j;
			valueMax = temp;
		}
	}
	if( indexMax > i ) {
for (k = 0; k < 6; ++k)
{
	acadoWorkspace.rk_dim6_swap = A[i*6+k];
	A[i*6+k] = A[indexMax*6+k];
	A[indexMax*6+k] = acadoWorkspace.rk_dim6_swap;
}
	acadoWorkspace.rk_dim6_swap = b[i];
	b[i] = b[indexMax];
	b[indexMax] = acadoWorkspace.rk_dim6_swap;
	intSwap = rk_perm[i];
	rk_perm[i] = rk_perm[indexMax];
	rk_perm[indexMax] = intSwap;
	}
	det *= A[i*6+i];
	for( j=i+1; j < 6; j++ ) {
		A[j*6+i] = -A[j*6+i]/A[i*6+i];
		for( k=i+1; k < 6; k++ ) {
			A[j*6+k] += A[j*6+i] * A[i*6+k];
		}
		b[j] += A[j*6+i] * b[i];
	}
}
det *= A[35];
det = fabs(det);
acado_solve_dim6_triangular( A, b );
return det;
}

void acado_solve_dim6_system_reuse( real_t* const A, real_t* const b, int* const rk_perm )
{

acadoWorkspace.rk_dim6_bPerm[0] = b[rk_perm[0]];
acadoWorkspace.rk_dim6_bPerm[1] = b[rk_perm[1]];
acadoWorkspace.rk_dim6_bPerm[2] = b[rk_perm[2]];
acadoWorkspace.rk_dim6_bPerm[3] = b[rk_perm[3]];
acadoWorkspace.rk_dim6_bPerm[4] = b[rk_perm[4]];
acadoWorkspace.rk_dim6_bPerm[5] = b[rk_perm[5]];
acadoWorkspace.rk_dim6_bPerm[1] += A[6]*acadoWorkspace.rk_dim6_bPerm[0];

acadoWorkspace.rk_dim6_bPerm[2] += A[12]*acadoWorkspace.rk_dim6_bPerm[0];
acadoWorkspace.rk_dim6_bPerm[2] += A[13]*acadoWorkspace.rk_dim6_bPerm[1];

acadoWorkspace.rk_dim6_bPerm[3] += A[18]*acadoWorkspace.rk_dim6_bPerm[0];
acadoWorkspace.rk_dim6_bPerm[3] += A[19]*acadoWorkspace.rk_dim6_bPerm[1];
acadoWorkspace.rk_dim6_bPerm[3] += A[20]*acadoWorkspace.rk_dim6_bPerm[2];

acadoWorkspace.rk_dim6_bPerm[4] += A[24]*acadoWorkspace.rk_dim6_bPerm[0];
acadoWorkspace.rk_dim6_bPerm[4] += A[25]*acadoWorkspace.rk_dim6_bPerm[1];
acadoWorkspace.rk_dim6_bPerm[4] += A[26]*acadoWorkspace.rk_dim6_bPerm[2];
acadoWorkspace.rk_dim6_bPerm[4] += A[27]*acadoWorkspace.rk_dim6_bPerm[3];

acadoWorkspace.rk_dim6_bPerm[5] += A[30]*acadoWorkspace.rk_dim6_bPerm[0];
acadoWorkspace.rk_dim6_bPerm[5] += A[31]*acadoWorkspace.rk_dim6_bPerm[1];
acadoWorkspace.rk_dim6_bPerm[5] += A[32]*acadoWorkspace.rk_dim6_bPerm[2];
acadoWorkspace.rk_dim6_bPerm[5] += A[33]*acadoWorkspace.rk_dim6_bPerm[3];
acadoWorkspace.rk_dim6_bPerm[5] += A[34]*acadoWorkspace.rk_dim6_bPerm[4];


acado_solve_dim6_triangular( A, acadoWorkspace.rk_dim6_bPerm );
b[0] = acadoWorkspace.rk_dim6_bPerm[0];
b[1] = acadoWorkspace.rk_dim6_bPerm[1];
b[2] = acadoWorkspace.rk_dim6_bPerm[2];
b[3] = acadoWorkspace.rk_dim6_bPerm[3];
b[4] = acadoWorkspace.rk_dim6_bPerm[4];
b[5] = acadoWorkspace.rk_dim6_bPerm[5];
}



/** Column vector of size: 1 */
static const real_t acado_Ah_mat[ 1 ] = 
{ 1.6666666666666666e-02 };


/* Fixed step size:0.0333333 */
int acado_integrate( real_t* const rk_eta, int resetIntegrator )
{
int error;

int i;
int j;
int k;
int run;
int run1;
int tmp_index1;
int tmp_index2;

real_t det;

acadoWorkspace.rk_ttt = 0.0000000000000000e+00;
acadoWorkspace.rk_xxx[6] = rk_eta[60];
acadoWorkspace.rk_xxx[7] = rk_eta[61];
acadoWorkspace.rk_xxx[8] = rk_eta[62];
acadoWorkspace.rk_xxx[9] = rk_eta[63];
acadoWorkspace.rk_xxx[10] = rk_eta[64];
acadoWorkspace.rk_xxx[11] = rk_eta[65];
acadoWorkspace.rk_xxx[12] = rk_eta[66];

for (run = 0; run < 3; ++run)
{
if( run > 0 ) {
for (i = 0; i < 6; ++i)
{
acadoWorkspace.rk_diffsPrev2[i * 9] = rk_eta[i * 6 + 6];
acadoWorkspace.rk_diffsPrev2[i * 9 + 1] = rk_eta[i * 6 + 7];
acadoWorkspace.rk_diffsPrev2[i * 9 + 2] = rk_eta[i * 6 + 8];
acadoWorkspace.rk_diffsPrev2[i * 9 + 3] = rk_eta[i * 6 + 9];
acadoWorkspace.rk_diffsPrev2[i * 9 + 4] = rk_eta[i * 6 + 10];
acadoWorkspace.rk_diffsPrev2[i * 9 + 5] = rk_eta[i * 6 + 11];
acadoWorkspace.rk_diffsPrev2[i * 9 + 6] = rk_eta[i * 3 + 42];
acadoWorkspace.rk_diffsPrev2[i * 9 + 7] = rk_eta[i * 3 + 43];
acadoWorkspace.rk_diffsPrev2[i * 9 + 8] = rk_eta[i * 3 + 44];
}
}
if( resetIntegrator ) {
for (i = 0; i < 1; ++i)
{
for (run1 = 0; run1 < 1; ++run1)
{
for (j = 0; j < 6; ++j)
{
acadoWorkspace.rk_xxx[j] = rk_eta[j];
tmp_index1 = j;
acadoWorkspace.rk_xxx[j] += + acado_Ah_mat[run1]*acadoWorkspace.rk_kkk[tmp_index1];
}
acado_diffs( acadoWorkspace.rk_xxx, &(acadoWorkspace.rk_diffsTemp2[ run1 * 54 ]) );
for (j = 0; j < 6; ++j)
{
tmp_index1 = (run1 * 6) + (j);
acadoWorkspace.rk_A[tmp_index1 * 6] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 1] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 1)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 2] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 2)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 3] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 3)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 4] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 4)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 5] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 5)];
if( 0 == run1 ) acadoWorkspace.rk_A[(tmp_index1 * 6) + (j)] -= 1.0000000000000000e+00;
}
acado_rhs( acadoWorkspace.rk_xxx, acadoWorkspace.rk_rhsTemp );
acadoWorkspace.rk_b[run1 * 6] = acadoWorkspace.rk_kkk[run1] - acadoWorkspace.rk_rhsTemp[0];
acadoWorkspace.rk_b[run1 * 6 + 1] = acadoWorkspace.rk_kkk[run1 + 1] - acadoWorkspace.rk_rhsTemp[1];
acadoWorkspace.rk_b[run1 * 6 + 2] = acadoWorkspace.rk_kkk[run1 + 2] - acadoWorkspace.rk_rhsTemp[2];
acadoWorkspace.rk_b[run1 * 6 + 3] = acadoWorkspace.rk_kkk[run1 + 3] - acadoWorkspace.rk_rhsTemp[3];
acadoWorkspace.rk_b[run1 * 6 + 4] = acadoWorkspace.rk_kkk[run1 + 4] - acadoWorkspace.rk_rhsTemp[4];
acadoWorkspace.rk_b[run1 * 6 + 5] = acadoWorkspace.rk_kkk[run1 + 5] - acadoWorkspace.rk_rhsTemp[5];
}
det = acado_solve_dim6_system( acadoWorkspace.rk_A, acadoWorkspace.rk_b, acadoWorkspace.rk_dim6_perm );
for (j = 0; j < 1; ++j)
{
acadoWorkspace.rk_kkk[j] += acadoWorkspace.rk_b[j * 6];
acadoWorkspace.rk_kkk[j + 1] += acadoWorkspace.rk_b[j * 6 + 1];
acadoWorkspace.rk_kkk[j + 2] += acadoWorkspace.rk_b[j * 6 + 2];
acadoWorkspace.rk_kkk[j + 3] += acadoWorkspace.rk_b[j * 6 + 3];
acadoWorkspace.rk_kkk[j + 4] += acadoWorkspace.rk_b[j * 6 + 4];
acadoWorkspace.rk_kkk[j + 5] += acadoWorkspace.rk_b[j * 6 + 5];
}
}
}
for (i = 0; i < 5; ++i)
{
for (run1 = 0; run1 < 1; ++run1)
{
for (j = 0; j < 6; ++j)
{
acadoWorkspace.rk_xxx[j] = rk_eta[j];
tmp_index1 = j;
acadoWorkspace.rk_xxx[j] += + acado_Ah_mat[run1]*acadoWorkspace.rk_kkk[tmp_index1];
}
acado_rhs( acadoWorkspace.rk_xxx, acadoWorkspace.rk_rhsTemp );
acadoWorkspace.rk_b[run1 * 6] = acadoWorkspace.rk_kkk[run1] - acadoWorkspace.rk_rhsTemp[0];
acadoWorkspace.rk_b[run1 * 6 + 1] = acadoWorkspace.rk_kkk[run1 + 1] - acadoWorkspace.rk_rhsTemp[1];
acadoWorkspace.rk_b[run1 * 6 + 2] = acadoWorkspace.rk_kkk[run1 + 2] - acadoWorkspace.rk_rhsTemp[2];
acadoWorkspace.rk_b[run1 * 6 + 3] = acadoWorkspace.rk_kkk[run1 + 3] - acadoWorkspace.rk_rhsTemp[3];
acadoWorkspace.rk_b[run1 * 6 + 4] = acadoWorkspace.rk_kkk[run1 + 4] - acadoWorkspace.rk_rhsTemp[4];
acadoWorkspace.rk_b[run1 * 6 + 5] = acadoWorkspace.rk_kkk[run1 + 5] - acadoWorkspace.rk_rhsTemp[5];
}
acado_solve_dim6_system_reuse( acadoWorkspace.rk_A, acadoWorkspace.rk_b, acadoWorkspace.rk_dim6_perm );
for (j = 0; j < 1; ++j)
{
acadoWorkspace.rk_kkk[j] += acadoWorkspace.rk_b[j * 6];
acadoWorkspace.rk_kkk[j + 1] += acadoWorkspace.rk_b[j * 6 + 1];
acadoWorkspace.rk_kkk[j + 2] += acadoWorkspace.rk_b[j * 6 + 2];
acadoWorkspace.rk_kkk[j + 3] += acadoWorkspace.rk_b[j * 6 + 3];
acadoWorkspace.rk_kkk[j + 4] += acadoWorkspace.rk_b[j * 6 + 4];
acadoWorkspace.rk_kkk[j + 5] += acadoWorkspace.rk_b[j * 6 + 5];
}
}
for (run1 = 0; run1 < 1; ++run1)
{
for (j = 0; j < 6; ++j)
{
acadoWorkspace.rk_xxx[j] = rk_eta[j];
tmp_index1 = j;
acadoWorkspace.rk_xxx[j] += + acado_Ah_mat[run1]*acadoWorkspace.rk_kkk[tmp_index1];
}
acado_diffs( acadoWorkspace.rk_xxx, &(acadoWorkspace.rk_diffsTemp2[ run1 * 54 ]) );
for (j = 0; j < 6; ++j)
{
tmp_index1 = (run1 * 6) + (j);
acadoWorkspace.rk_A[tmp_index1 * 6] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 1] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 1)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 2] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 2)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 3] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 3)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 4] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 4)];
acadoWorkspace.rk_A[tmp_index1 * 6 + 5] = + acado_Ah_mat[run1]*acadoWorkspace.rk_diffsTemp2[(run1 * 54) + (j * 9 + 5)];
if( 0 == run1 ) acadoWorkspace.rk_A[(tmp_index1 * 6) + (j)] -= 1.0000000000000000e+00;
}
}
for (run1 = 0; run1 < 6; ++run1)
{
for (i = 0; i < 1; ++i)
{
acadoWorkspace.rk_b[i * 6] = - acadoWorkspace.rk_diffsTemp2[(i * 54) + (run1)];
acadoWorkspace.rk_b[i * 6 + 1] = - acadoWorkspace.rk_diffsTemp2[(i * 54) + (run1 + 9)];
acadoWorkspace.rk_b[i * 6 + 2] = - acadoWorkspace.rk_diffsTemp2[(i * 54) + (run1 + 18)];
acadoWorkspace.rk_b[i * 6 + 3] = - acadoWorkspace.rk_diffsTemp2[(i * 54) + (run1 + 27)];
acadoWorkspace.rk_b[i * 6 + 4] = - acadoWorkspace.rk_diffsTemp2[(i * 54) + (run1 + 36)];
acadoWorkspace.rk_b[i * 6 + 5] = - acadoWorkspace.rk_diffsTemp2[(i * 54) + (run1 + 45)];
}
if( 0 == run1 ) {
det = acado_solve_dim6_system( acadoWorkspace.rk_A, acadoWorkspace.rk_b, acadoWorkspace.rk_dim6_perm );
}
 else {
acado_solve_dim6_system_reuse( acadoWorkspace.rk_A, acadoWorkspace.rk_b, acadoWorkspace.rk_dim6_perm );
}
for (i = 0; i < 1; ++i)
{
acadoWorkspace.rk_diffK[i] = acadoWorkspace.rk_b[i * 6];
acadoWorkspace.rk_diffK[i + 1] = acadoWorkspace.rk_b[i * 6 + 1];
acadoWorkspace.rk_diffK[i + 2] = acadoWorkspace.rk_b[i * 6 + 2];
acadoWorkspace.rk_diffK[i + 3] = acadoWorkspace.rk_b[i * 6 + 3];
acadoWorkspace.rk_diffK[i + 4] = acadoWorkspace.rk_b[i * 6 + 4];
acadoWorkspace.rk_diffK[i + 5] = acadoWorkspace.rk_b[i * 6 + 5];
}
for (i = 0; i < 6; ++i)
{
acadoWorkspace.rk_diffsNew2[(i * 9) + (run1)] = (i == run1-0);
acadoWorkspace.rk_diffsNew2[(i * 9) + (run1)] += + acadoWorkspace.rk_diffK[i]*(real_t)3.3333333333333333e-02;
}
}
for (run1 = 0; run1 < 3; ++run1)
{
for (i = 0; i < 1; ++i)
{
for (j = 0; j < 6; ++j)
{
tmp_index1 = (i * 6) + (j);
tmp_index2 = (run1) + (j * 9);
acadoWorkspace.rk_b[tmp_index1] = - acadoWorkspace.rk_diffsTemp2[(i * 54) + (tmp_index2 + 6)];
}
}
acado_solve_dim6_system_reuse( acadoWorkspace.rk_A, acadoWorkspace.rk_b, acadoWorkspace.rk_dim6_perm );
for (i = 0; i < 1; ++i)
{
acadoWorkspace.rk_diffK[i] = acadoWorkspace.rk_b[i * 6];
acadoWorkspace.rk_diffK[i + 1] = acadoWorkspace.rk_b[i * 6 + 1];
acadoWorkspace.rk_diffK[i + 2] = acadoWorkspace.rk_b[i * 6 + 2];
acadoWorkspace.rk_diffK[i + 3] = acadoWorkspace.rk_b[i * 6 + 3];
acadoWorkspace.rk_diffK[i + 4] = acadoWorkspace.rk_b[i * 6 + 4];
acadoWorkspace.rk_diffK[i + 5] = acadoWorkspace.rk_b[i * 6 + 5];
}
for (i = 0; i < 6; ++i)
{
acadoWorkspace.rk_diffsNew2[(i * 9) + (run1 + 6)] = + acadoWorkspace.rk_diffK[i]*(real_t)3.3333333333333333e-02;
}
}
rk_eta[0] += + acadoWorkspace.rk_kkk[0]*(real_t)3.3333333333333333e-02;
rk_eta[1] += + acadoWorkspace.rk_kkk[1]*(real_t)3.3333333333333333e-02;
rk_eta[2] += + acadoWorkspace.rk_kkk[2]*(real_t)3.3333333333333333e-02;
rk_eta[3] += + acadoWorkspace.rk_kkk[3]*(real_t)3.3333333333333333e-02;
rk_eta[4] += + acadoWorkspace.rk_kkk[4]*(real_t)3.3333333333333333e-02;
rk_eta[5] += + acadoWorkspace.rk_kkk[5]*(real_t)3.3333333333333333e-02;
if( run == 0 ) {
for (i = 0; i < 6; ++i)
{
for (j = 0; j < 6; ++j)
{
tmp_index2 = (j) + (i * 6);
rk_eta[tmp_index2 + 6] = acadoWorkspace.rk_diffsNew2[(i * 9) + (j)];
}
for (j = 0; j < 3; ++j)
{
tmp_index2 = (j) + (i * 3);
rk_eta[tmp_index2 + 42] = acadoWorkspace.rk_diffsNew2[(i * 9) + (j + 6)];
}
}
}
else {
for (i = 0; i < 6; ++i)
{
for (j = 0; j < 6; ++j)
{
tmp_index2 = (j) + (i * 6);
rk_eta[tmp_index2 + 6] = + acadoWorkspace.rk_diffsNew2[i * 9]*acadoWorkspace.rk_diffsPrev2[j];
rk_eta[tmp_index2 + 6] += + acadoWorkspace.rk_diffsNew2[i * 9 + 1]*acadoWorkspace.rk_diffsPrev2[j + 9];
rk_eta[tmp_index2 + 6] += + acadoWorkspace.rk_diffsNew2[i * 9 + 2]*acadoWorkspace.rk_diffsPrev2[j + 18];
rk_eta[tmp_index2 + 6] += + acadoWorkspace.rk_diffsNew2[i * 9 + 3]*acadoWorkspace.rk_diffsPrev2[j + 27];
rk_eta[tmp_index2 + 6] += + acadoWorkspace.rk_diffsNew2[i * 9 + 4]*acadoWorkspace.rk_diffsPrev2[j + 36];
rk_eta[tmp_index2 + 6] += + acadoWorkspace.rk_diffsNew2[i * 9 + 5]*acadoWorkspace.rk_diffsPrev2[j + 45];
}
for (j = 0; j < 3; ++j)
{
tmp_index2 = (j) + (i * 3);
rk_eta[tmp_index2 + 42] = acadoWorkspace.rk_diffsNew2[(i * 9) + (j + 6)];
rk_eta[tmp_index2 + 42] += + acadoWorkspace.rk_diffsNew2[i * 9]*acadoWorkspace.rk_diffsPrev2[j + 6];
rk_eta[tmp_index2 + 42] += + acadoWorkspace.rk_diffsNew2[i * 9 + 1]*acadoWorkspace.rk_diffsPrev2[j + 15];
rk_eta[tmp_index2 + 42] += + acadoWorkspace.rk_diffsNew2[i * 9 + 2]*acadoWorkspace.rk_diffsPrev2[j + 24];
rk_eta[tmp_index2 + 42] += + acadoWorkspace.rk_diffsNew2[i * 9 + 3]*acadoWorkspace.rk_diffsPrev2[j + 33];
rk_eta[tmp_index2 + 42] += + acadoWorkspace.rk_diffsNew2[i * 9 + 4]*acadoWorkspace.rk_diffsPrev2[j + 42];
rk_eta[tmp_index2 + 42] += + acadoWorkspace.rk_diffsNew2[i * 9 + 5]*acadoWorkspace.rk_diffsPrev2[j + 51];
}
}
}
resetIntegrator = 0;
acadoWorkspace.rk_ttt += 3.3333333333333331e-01;
}
for (i = 0; i < 6; ++i)
{
}
if( det < 1e-12 ) {
error = 2;
} else if( det < 1e-6 ) {
error = 1;
} else {
error = 0;
}
return error;
}



