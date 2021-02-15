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
 
    TIME autotime;
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
    SIMexport ExportModule1( 1, 0.0001 );
    ExportModule1.set( GENERATE_MATLAB_INTERFACE, 1 );
    uint options_flag;
    options_flag = ExportModule1.set( INTEGRATOR_TYPE, INT_IRK_GL2 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: INTEGRATOR_TYPE");
    options_flag = ExportModule1.set( NUM_INTEGRATOR_STEPS, 3 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: NUM_INTEGRATOR_STEPS");
    DifferentialEquation acadodata_f1;
    acadodata_f1 << dot(xx) == xxd;
    acadodata_f1 << dot(xy) == xyd;
    acadodata_f1 << dot(xz) == xzd;
    acadodata_f1 << dot(xxd) == (-((-2.00000000000000000000e+00)*xx+(-2.00000000000000000000e+00)*xyd)*z1-(-(-1.22984134246966383269e-02)+xx)*9.87851000000000034618e-01/pow(sqrt((pow((-(-1.22984134246966383269e-02)+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)-(-1.00000000000000000000e+00)*xx*z2-(-1.00000000000000000000e+00)*xy*z3-(-2.00000000000000000000e+00)*xyd-(-9.87851000000000034618e-01+xx)*1.21490000000000000768e-02/pow(sqrt((pow((-9.87851000000000034618e-01+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)+ux+xx);
    acadodata_f1 << dot(xyd) == (-((-2.00000000000000000000e+00)*xy+2.00000000000000000000e+00*xxd)*z1-(-1.00000000000000000000e+00)*xy*z2-1.21490000000000000768e-02/pow(sqrt((pow((-9.87851000000000034618e-01+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)*xy-2.00000000000000000000e+00*xxd-9.87851000000000034618e-01/pow(sqrt((pow((-(-1.22984134246966383269e-02)+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)*xy+uy-xx*z3+xy);
    acadodata_f1 << dot(xzd) == (-1.21490000000000000768e-02/pow(sqrt((pow((-9.87851000000000034618e-01+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)*xz-9.87851000000000034618e-01/pow(sqrt((pow((-(-1.22984134246966383269e-02)+xx),2.00000000000000000000e+00)+pow(xy,2.00000000000000000000e+00)+pow(xz,2.00000000000000000000e+00))),3.00000000000000000000e+00)*xz+uz);

    ExportModule1.setModel( acadodata_f1 );

    uint export_flag = 0;
    ExportModule1.setTimingSteps( 0 );
    export_flag = ExportModule1.exportCode( "export_SIM" );
    if(export_flag != 0) mexErrMsgTxt("ACADO export failed because of the above error(s)!");


    clearAllStaticCounters( ); 
 
} 

