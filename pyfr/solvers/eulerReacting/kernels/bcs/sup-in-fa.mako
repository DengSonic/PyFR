# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% N = int(c['component']) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, ploc, t'>
    ur[0] = ${c['rho']};
  //  ur[0] = 1.813122272491455;

    ur[${nvars - 1}] = ur[0]*0.01;   
    ur[${nvars - 2}] = ur[0]*0.99;   
 
    //9,11                            
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = (${c['rho']})*(${c[v]});
% endfor
  //  ur[1] = -113.5*ur[0];
  //  ur[2] = 0.0;

    ur[${nvars - 3}] = ${c['p']}/${c['gamma'] - 1} + 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
    
  // ur[${nvars - 3}] = 156980.0/(1.4-1.0)+0.5*1.813122272491455*113.5*113.5;
</%pyfr:macro>
