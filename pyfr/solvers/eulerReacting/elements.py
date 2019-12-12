# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements


class BaseFluidReactingElements(object):
    formulations = ['std', 'dual']

    privarmap = {2: ['rho', 'u', 'v', 'p','Y0','Y1'],
                 3: ['rho', 'u', 'v', 'w', 'p','Y0','Y1']}


    convarmap = {2: ['rho', 'rhou', 'rhov', 'E','rhoY0','rhoY1'],
                 3: ['rho', 'rhou', 'rhov', 'rhow', 'E','rhoY0','rhoY1']}


    dualcoeffs = convarmap

    visvarmap = {
        2: [('density', ['rho']),
            ('velocity', ['u', 'v']),
            ('pressure', ['p']),
            ('massFraction0', ['Y0']),
            ('massFraction1', ['Y1'])],
        3: [('density', ['rho']),
            ('velocity', ['u', 'v', 'w']),
            ('pressure', ['p']),
            ('massFraction0',['Y0']),
            ('massFraction1',['Y1'])]
    }


    @staticmethod
    def pri_to_con(pris, cfg):   #????

        N = cfg.getint('constants', 'component')
        Npe = N+1
        #rho, Y1 = pris[0], pris[-1]
        rho = pris[0]
        
        p = pris[-Npe]

        #Y1 = pris[-1]
        Y = [1.0*y for y in pris[-N:]]
        #for i in range(N)
        #    Y[i]=pris[-N]

        # Multiply velocity components by rho
        rhovs = [rho*c for c in pris[1:-Npe]]  #from 1 to -1 not include -1

        # Compute the energy
        Cp = [None] * N
        Wmol = [None] * N
        for i in range(N):
          Cp[i]=cfg.getfloat('constants', 'Cp'+str(i))
          Wmol[i]=cfg.getfloat('constants', 'Wmol'+str(i))

        R=cfg.getfloat('constants', 'Rgas')

        CpBar = 0.0
        RBar =0.0
        for i in range(N):
          CpBar=CpBar+Y[i]*Cp[i]
          RBar=RBar+R*(Y[i]/Wmol[i])        

        #gamma = cfg.getfloat('constants', 'gamma'+str(0))
        gamma = CpBar/(CpBar-RBar)

        #reference energy
        E_ref = [None] * N
        for i in range(N):
          E_ref[i]=cfg.getfloat('constants', 'E_ref'+str(i))

        E_refBar = 0.0
        for i in range(N):
          E_refBar=E_refBar+Y[i]*E_ref[i]        
        
        E = p/(gamma - 1) + 0.5*rho*sum(c*c for c in pris[1:-Npe])+rho*E_refBar  
        #rhoY1=rho*Y1
        #rhoY = [rho*y for y in Y]
        rhoY = [rho*y for y in pris[-N:]]

        return [rho] + rhovs + [E] + rhoY

    @staticmethod
    def con_to_pri(cons, cfg):

        N = cfg.getint('constants', 'component')
        Npe = N+1
        #rho, rhoY1 = cons[0], cons[-1]

        #rho = cons[0]
        #recaculate density
        rho =0.0;
        for i in range(N):
          rho=rho+cons[-N+i]
        
        E = cons[-Npe]
        # Divide momentum components by rho
        vs = [rhov/rho for rhov in cons[1:-Npe]]

        Y = [rhoy/rho for rhoy in cons[-N:]]

        # Compute the pressure
        Cp = [None] * N
        Wmol = [None] * N
        for i in range(N):
          Cp[i]=cfg.getfloat('constants', 'Cp'+str(i))
          Wmol[i]=cfg.getfloat('constants', 'Wmol'+str(i))

        R=cfg.getfloat('constants', 'Rgas')

        CpBar = 0.0
        RBar =0.0
        for i in range(N):
          CpBar=CpBar+Y[i]*Cp[i]
          RBar=RBar+R*(Y[i]/Wmol[i])        

        #gamma = cfg.getfloat('constants', 'gamma'+str(0))
        gamma = CpBar/(CpBar-RBar)

        #reference energy
        E_ref = [None] * N
        for i in range(N):
          E_ref[i]=cfg.getfloat('constants', 'E_ref'+str(i))

        E_refBar = 0.0
        for i in range(N):
          E_refBar=E_refBar+Y[i]*E_ref[i] 

        p = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in vs)-rho*E_refBar)
  
        # Y1 = rhoY1/rho

        

        return [rho] + vs + [p] + Y



class EulerReactingElements(BaseFluidReactingElements, BaseAdvectionElements):
    def set_backend(self, backend, nscalupts, nonce):
        super().set_backend(backend, nscalupts, nonce)

        # Register our flux kernel
        backend.pointwise.register('pyfr.solvers.eulerReacting.kernels.tflux')

        # Template parameters for the flux kernel
        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.cfg.items_as('constants', float))

        if 'flux' in self.antialias:
            self.kernels['tdisf'] = lambda: backend.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts
            )
        else:
            self.kernels['tdisf'] = lambda: backend.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts
            )
