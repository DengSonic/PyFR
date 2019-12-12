# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokesReacting.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, ploc, t'>
    ur[0] = ${c['rho']};
    ur[${nvars - 1}] = ${c['rhoY1']};
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = (${c['rho']}) * (${c[v]});
% endfor
    ur[${nvars - 2}] = ul[${nvars - 2}]
                     - 0.5*(1.0/ul[0])*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))}
                     + 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
