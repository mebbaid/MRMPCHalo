/*
*    This file is part of ACADO Toolkit.
*
*    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
*    Copyright (C) 2008-2009 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
*    Developed within the Optimization in Engineering Center (OPTEC) under
*    supervision of Moritz Diehl. All rights reserved.
*
*    ACADO Toolkit is free software; you can redistribute it and/or
*    modify it under the terms of the GNU Lesser General Public
*    License as published by the Free Software Foundation; either
*    version 3 of the License, or (at your option) any later version.
*
*    ACADO Toolkit is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public
*    License along with ACADO Toolkit; if not, write to the Free Software
*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*
*/


/**
*    Author David Ariens, Rien Quirynen
*    Date 2009-2013
*    http://www.acadotoolkit.org/matlab 
*/

#include <acado_optimal_control.hpp>
#include <acado_toolkit.hpp>
#include <acado/utils/matlab_acado_utils.hpp>

USING_NAMESPACE_ACADO

#include <mex.h>


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
 { 
 
    MatlabConsoleStreamBuf mybuf;
    RedirectStream redirect(std::cout, mybuf);
    clearAllStaticCounters( ); 
 
    mexPrintf("\nACADO Toolkit for Matlab - Developed by David Ariens and Rien Quirynen, 2009-2013 \n"); 
    mexPrintf("Support available at http://www.acadotoolkit.org/matlab \n \n"); 

    if (nrhs != 0){ 
      mexErrMsgTxt("This problem expects 0 right hand side argument(s) since you have defined 0 MexInput(s)");
    } 
 
    DifferentialState xx;
    DifferentialState xy;
    DifferentialState xz;
    DifferentialState xxd;
    DifferentialState xyd;
    DifferentialState xzd;
    Control ux;
    Control uy;
    Control uz;
    OnlineData z1; 
    OnlineData z2; 
    OnlineData z3; 
    OnlineData z4; 
    BMatrix acadodata_M1;
    acadodata_M1.read( "nmpc_data_acadodata_M1.txt" );
    BMatrix acadodata_M2;
    acadodata_M2.read( "nmpc_data_acadodata_M2.txt" );
    Function acadodata_f2;
    acadodata_f2 << xx;
    acadodata_f2 << xy;
    acadodata_f2 << xz;
    acadodata_f2 << xxd;
    acadodata_f2 << xyd;
    acadodata_f2 << xzd;
    acadodata_f2 << ux;
    acadodata_f2 << uy;
    acadodata_f2 << uz;
    Function acadodata_f3;
    acadodata_f3 << xx;
    acadodata_f3 << xy;
    acadodata_f3 << xz;
    acadodata_f3 << xxd;
    acadodata_f3 << xyd;
    acadodata_f3 << xzd;
    OCP ocp1(0, 4, 40);
    ocp1.minimizeLSQ(acadodata_M1, acadodata_f2);
    ocp1.minimizeLSQEndTerm(acadodata_M2, acadodata_f3);
    DifferentialEquation acadodata_f4;
    acadodata_f4 << dot(xx) == xxd;
    acadodata_f4 << dot(xy) == xyd;
    acadodata_f4 << dot(xz) == xzd;
    acadodata_f4 << dot(xxd) == (-(-1.00000000000000000000e+00+1.21490000000000000768e-02+xx)*1.21490000000000000768e-02/pow(sqrt((pow((-9.87851000000000034618e-01+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)-(1.21490000000000000768e-02+xx)*9.87851000000000034618e-01/pow(sqrt((pow((-(-1.22984134246966383269e-02)+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)+2.00000000000000000000e+00*xyd+ux+xx);
    acadodata_f4 << dot(xyd) == ((-2.00000000000000000000e+00)*xxd-1.21490000000000000768e-02/pow(sqrt((pow((-9.87851000000000034618e-01+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)*xy-9.87851000000000034618e-01/pow(sqrt((pow((-(-1.22984134246966383269e-02)+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)*xy+uy+xy);
    acadodata_f4 << dot(xzd) == ((-9.87851000000000034618e-01)/pow(sqrt((pow((-(-1.22984134246966383269e-02)+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)*xz-1.21490000000000000768e-02/pow(sqrt((pow((-9.87851000000000034618e-01+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)*xz+uz);

    ocp1.setModel( acadodata_f4 );


    ocp1.setNU( 3 );
    ocp1.setNP( 0 );
    ocp1.setNOD( 4 );
    OCPexport ExportModule2( ocp1 );
    ExportModule2.set( GENERATE_MATLAB_INTERFACE, 1 );
    uint options_flag;
    options_flag = ExportModule2.set( HESSIAN_APPROXIMATION, GAUSS_NEWTON );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: HESSIAN_APPROXIMATION");
    options_flag = ExportModule2.set( DISCRETIZATION_TYPE, MULTIPLE_SHOOTING );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: DISCRETIZATION_TYPE");
    options_flag = ExportModule2.set( SPARSE_QP_SOLUTION, FULL_CONDENSING_N2 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: SPARSE_QP_SOLUTION");
    options_flag = ExportModule2.set( LEVENBERG_MARQUARDT, 1.000000E-05 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: LEVENBERG_MARQUARDT");
    options_flag = ExportModule2.set( INTEGRATOR_TYPE, INT_RK4 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: INTEGRATOR_TYPE");
    options_flag = ExportModule2.set( NUM_INTEGRATOR_STEPS, 120 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: NUM_INTEGRATOR_STEPS");
    options_flag = ExportModule2.set( QP_SOLVER, QP_QPOASES );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: QP_SOLVER");
    uint export_flag;
    export_flag = ExportModule2.exportCode( "export_MPC" );
    if(export_flag != 0) mexErrMsgTxt("ACADO export failed because of the above error(s)!");


    clearAllStaticCounters( ); 
 
} 

