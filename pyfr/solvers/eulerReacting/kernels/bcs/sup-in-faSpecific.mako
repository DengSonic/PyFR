# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% N = int(c['component']) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, ploc, t'>

    ur[0] = 2.41595;

    ur[${nvars - 1}] = ur[0]*0.01;   
    ur[${nvars - 2}] = ur[0]*0.99;   
 
    ur[1] = -230.27*ur[0];
    ur[2] = 0.0;

    
    ur[${nvars - 3}] = 249090.0/(1.4-1.0)+0.5*2.41595*230.27*230.27;

</%pyfr:macro>
