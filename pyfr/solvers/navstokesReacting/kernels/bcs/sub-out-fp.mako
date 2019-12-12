# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokesReacting.kernels.bcs.common'/>
  //6,9,12
<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, ploc, t'>
% for i in range(nvars - 2):
    ur[${i}] = ul[${i}];
% endfor
    ur[${nvars - 2}] = ${c['p']}/${c['gamma'] - 1}
                     + 0.5*(1.0/ul[0])*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))};
    ur[${nvars - 1}] = ul[${nvars - 1}];    
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
