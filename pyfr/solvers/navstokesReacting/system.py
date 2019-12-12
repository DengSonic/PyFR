# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.navstokesReacting.elements import NavierStokesReactingElements
from pyfr.solvers.navstokesReacting.inters import (NavierStokesReactingBaseBCInters,
                                           NavierStokesReactingIntInters,
                                           NavierStokesReactingMPIInters)


class NavierStokesReactingSystem(BaseAdvectionDiffusionSystem):
    name = 'navier-stokesReacting'

    elementscls = NavierStokesReactingElements
    intinterscls = NavierStokesReactingIntInters
    mpiinterscls = NavierStokesReactingMPIInters
    bbcinterscls = NavierStokesReactingBaseBCInters
