# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.eulerReacting.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.navstokesReacting.kernels.bcs.common'/>
<%include file='pyfr.solvers.navstokesReacting.kernels.flux'/>
//15 16

<% N = int(c['component']) %>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur, ploc, t'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
    ur[0] = ul[0];

% for i in range(nvars):
    ur[${i}] = ul[${i}];;
% endfor  

% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor

//    ur[${nvars - 1}] = ul[${nvars - 1}];
//    ur[${nvars - 2}] = ul[${nvars - 2}];
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, artviscl, nl, magnl, ploc, t'>
    // Ghost state r
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur', 'ploc', 't')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*(ficomm[${i}]);
% endfor
</%pyfr:macro>
