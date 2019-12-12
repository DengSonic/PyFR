# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokesReacting.kernels.bcs.common'/>
   //7, 11
<% N = int(c['component']) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, ploc, t'>
    ur[0] = ul[0];
   // ur[${nvars - 1}] = ul[${nvars - 1}];

% for i in range(nvars):
    ur[${i}] = ul[${i}];;
% endfor 

% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = -ul[${i + 1}] + 2*${c[v]}*ul[0];
% endfor
    ur[${nvars - N-1}] = ${c['cpTw']/c['gamma']}*ur[0]
                     + 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur, ploc, t'>
    ur[0] = ul[0];
   // ur[${nvars - 1}] = ul[${nvars - 1}];

% for i in range(nvars):
    ur[${i}] = ul[${i}];;
% endfor

% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = ${c[v]}*ul[0];
% endfor

    ur[${nvars - N-1}] = ${c['cpTw']/c['gamma']}*ur[0]
                     + 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
