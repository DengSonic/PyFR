# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% N = int(c['component']) %>

    
<%pyfr:macro name='inviscid_flux' params='s, f, p, v, a'>

    //fpdtype_t N = ${int(c['component'])};
    //N = 1;

    fpdtype_t invrho = 1.0/s[0];
    fpdtype_t E = s[${nvars - N - 1}];

    //fpdtype_t y=invrho*s[${nvars - 1}];
    //define and calculate mass fraction
    fpdtype_t y[${N}];
% for i in range(N):
    y[${i}] = invrho*s[${nvars - N + i}];
% endfor                   

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
% endfor



    // Compute the pressure
    // int N = ${int(c['component'])};

    fpdtype_t R= ${c['Rgas']};

    fpdtype_t Cp[${N}];
    fpdtype_t Wmol[${N}];
    fpdtype_t E_ref[${N}];
% for i in range(N):
    Cp[${i}] = ${c['Cp'+str(i)]};
    Wmol[${i}] = ${c['Wmol'+str(i)]};
    E_ref[${i}] = ${c['E_ref'+str(i)]};
% endfor

    fpdtype_t CpBar = 0.0;
    fpdtype_t RBar = 0.0;
    fpdtype_t E_refBar = 0.0;
% for i in range(N):
    CpBar = CpBar + max(y[${i}],1.0e-8)*Cp[${i}];
    RBar = RBar + R*(max(y[${i}],1.0e-8)/Wmol[${i}]);
    E_refBar = E_refBar + max(y[${i}],1.0e-8)*E_ref[${i}];
% endfor

    fpdtype_t gamma = CpBar/(CpBar-RBar);

    p = (gamma - 1)*((E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)})-(s[0]*E_refBar));

    //simplified method
    a=sqrt(max(gamma*(p/s[0]),1.0e-8));

   // a=sqrt(max(gamma*(p/sumRho),1.0e-8));

    // p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)}-s[0]*E_refBar);

    // Density and energy fluxes, fl[${ndims}][${nvars}], add flux for species
    //f[${i}][${nvars - 1}] = rhov[${i}]*y;
% for i in range(ndims):
    f[${i}][0] = rhov[${i}];                       
    f[${i}][${nvars - N - 1}] = (E + p)*v[${i}];
%   for j in range(N):
      f[${i}][${nvars - N + j}] = rhov[${i}]*y[${j}]; 
%   endfor
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = rhov[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor
</%pyfr:macro>
