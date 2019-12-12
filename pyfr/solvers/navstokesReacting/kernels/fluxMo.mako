# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% N = int(c['component']) %>

% if ndims == 2:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>

    // correct sum rho later

    fpdtype_t rho = uin[0], rhou = uin[1], rhov = uin[2], E = uin[3];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov;

    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];

    // Velocity derivatives (rho*grad[u,v])
    fpdtype_t u_x = grad_uin[0][1] - u*rho_x;
    fpdtype_t u_y = grad_uin[1][1] - u*rho_y;
    fpdtype_t v_x = grad_uin[0][2] - v*rho_x;
    fpdtype_t v_y = grad_uin[1][2] - v*rho_y;

    fpdtype_t E_x = grad_uin[0][3];
    fpdtype_t E_y = grad_uin[1][3];

% if visc_corr == 'sutherland':
    // Compute the temperature and viscosity
    fpdtype_t cpT = ${c['gamma']}*(rcprho*E - 0.5*(u*u + v*v));
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}*Trat*sqrt(Trat)
                   / (cpT + ${c['cpTs']});
% else:
    fpdtype_t mu_c = ${c['mu']};
% endif

    // Compute temperature derivatives (c_v*dT/d[x,y])
    fpdtype_t T_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x));
    fpdtype_t T_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y));

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*mu_c*rcprho*(u_x - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_yy = -2*mu_c*rcprho*(v_y - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_xy = -mu_c*rcprho*(v_x + u_y);

    fout[0][1] += t_xx;     fout[1][1] += t_xy;
    fout[0][2] += t_xy;     fout[1][2] += t_yy;

    fout[0][3] += u*t_xx + v*t_xy + -mu_c*${c['gamma']/c['Pr']}*T_x;
    fout[1][3] += u*t_xy + v*t_yy + -mu_c*${c['gamma']/c['Pr']}*T_y;


    fpdtype_t Y[${N}];
% for i in range(N):
    Y[${i}] = rcprho*uin[${4+i}];
% endfor  

    //mass fraction derivatives
    fpdtype_t Y_x[${N}];
    fpdtype_t Y_y[${N}];

% for i in range(N):
    Y_x[${i}] = rcprho*(grad_uin[0][${4+i}]-Y[${i}]*rho_x);
    Y_y[${i}] = rcprho*(grad_uin[1][${4+i}]-Y[${i}]*rho_y);
% endfor     

    fpdtype_t D[${N}];
% for i in range(N):
    D[${i}] = ${c['D'+str(i)]};
% endfor

    fpdtype_t sumDY_x = 0.0;
    fpdtype_t sumDY_y = 0.0;

% for i in range(N):
    sumDY_x = sumDY_x+D[${i}]*Y_x[${i}];
    sumDY_y = sumDY_y+D[${i}]*Y_y[${i}];
% endfor

    fpdtype_t rhoYV_x[${N}];
    fpdtype_t rhoYV_y[${N}];

% for i in range(N):
    rhoYV_x[${i}] = rho*D[${i}]*Y_x[${i}]+rho*Y[${i}]*sumDY_x;
    rhoYV_y[${i}] = rho*D[${i}]*Y_y[${i}]+rho*Y[${i}]*sumDY_y;
% endfor 

    //+,- 
% for i in range(N):
    fout[0][${4+i}] += -rhoYV_x[${i}];
    fout[1][${4+i}] += -rhoYV_y[${i}];
% endfor 

    //correct energy equation due to diffusion
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
    CpBar = CpBar + max(Y[${i}],1.0e-8)*Cp[${i}];
    RBar = RBar + R*(max(Y[${i}],1.0e-8)/Wmol[${i}]);
    E_refBar = E_refBar + max(Y[${i}],1.0e-8)*E_ref[${i}];
% endfor

    fpdtype_t gamma = CpBar/(CpBar-RBar);

    fpdtype_t p = (gamma - 1)*((E - 0.5*rho*(u*u + v*v))-(rho*E_refBar));

    fpdtype_t T =p/rho/RBar;


    fpdtype_t sumEc_x = 0.0;
    fpdtype_t sumEc_y = 0.0;
    //(E_ref[${i}]+Cp[${i}]*T)
% for i in range(N):
    sumEc_x = sumEc_x+rhoYV_x[${i}]*(E_ref[${i}]);
    sumEc_y = sumEc_y+rhoYV_y[${i}]*(E_ref[${i}]);
% endfor

    fout[0][3] += sumEc_x;
    fout[1][3] += sumEc_y;

</%pyfr:macro>
% elif ndims == 3:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
    fpdtype_t rho  = uin[0];
    fpdtype_t rhou = uin[1], rhov = uin[2], rhow = uin[3];
    fpdtype_t E    = uin[4];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov, w = rcprho*rhow;

    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];
    fpdtype_t rho_z = grad_uin[2][0];

    // Velocity derivatives (rho*grad[u,v,w])
    fpdtype_t u_x = grad_uin[0][1] - u*rho_x;
    fpdtype_t u_y = grad_uin[1][1] - u*rho_y;
    fpdtype_t u_z = grad_uin[2][1] - u*rho_z;
    fpdtype_t v_x = grad_uin[0][2] - v*rho_x;
    fpdtype_t v_y = grad_uin[1][2] - v*rho_y;
    fpdtype_t v_z = grad_uin[2][2] - v*rho_z;
    fpdtype_t w_x = grad_uin[0][3] - w*rho_x;
    fpdtype_t w_y = grad_uin[1][3] - w*rho_y;
    fpdtype_t w_z = grad_uin[2][3] - w*rho_z;

    fpdtype_t E_x = grad_uin[0][4];
    fpdtype_t E_y = grad_uin[1][4];
    fpdtype_t E_z = grad_uin[2][4];

% if visc_corr == 'sutherland':
    // Compute the temperature and viscosity
    fpdtype_t cpT = ${c['gamma']}*(rcprho*E - 0.5*(u*u + v*v + w*w));
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}*Trat*sqrt(Trat)
                   / (cpT + ${c['cpTs']});
% else:
    fpdtype_t mu_c = ${c['mu']};
% endif

    // Compute temperature derivatives (c_v*dT/d[x,y,z])
    fpdtype_t T_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x + w*w_x));
    fpdtype_t T_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y + w*w_y));
    fpdtype_t T_z = rcprho*(E_z - (rcprho*rho_z*E + u*u_z + v*v_z + w*w_z));

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*mu_c*rcprho*(u_x - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_yy = -2*mu_c*rcprho*(v_y - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_zz = -2*mu_c*rcprho*(w_z - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_xy = -mu_c*rcprho*(v_x + u_y);
    fpdtype_t t_xz = -mu_c*rcprho*(u_z + w_x);
    fpdtype_t t_yz = -mu_c*rcprho*(w_y + v_z);

    fout[0][1] += t_xx;     fout[1][1] += t_xy;     fout[2][1] += t_xz;
    fout[0][2] += t_xy;     fout[1][2] += t_yy;     fout[2][2] += t_yz;
    fout[0][3] += t_xz;     fout[1][3] += t_yz;     fout[2][3] += t_zz;

    fout[0][4] += u*t_xx + v*t_xy + w*t_xz + -mu_c*${c['gamma']/c['Pr']}*T_x;
    fout[1][4] += u*t_xy + v*t_yy + w*t_yz + -mu_c*${c['gamma']/c['Pr']}*T_y;
    fout[2][4] += u*t_xz + v*t_yz + w*t_zz + -mu_c*${c['gamma']/c['Pr']}*T_z;
</%pyfr:macro>
% endif
