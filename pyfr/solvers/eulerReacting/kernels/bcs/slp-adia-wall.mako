# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% N = int(c['component']) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, ploc, t'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
    ur[0] = ul[0];
% for i in range(nvars):
    ur[${i}] = ul[${i}];;
% endfor  

    //14
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor

  //  ur[${nvars - 1}] = ul[${nvars - 1}];      
  //  ur[${nvars - N - 1}] = ul[${nvars - N - 1}];   
</%pyfr:macro>
