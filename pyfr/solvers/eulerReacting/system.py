# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.eulerReacting.elements import EulerReactingElements
from pyfr.solvers.eulerReacting.inters import (EulerReactingIntInters,EulerReactingMPIInters,EulerReactingBaseBCInters)
#from module import object(function, variable, class)
#import ... as..   : make module another name

class EulerReactingSystem(BaseAdvectionSystem):
    name = 'eulerReacting'

    elementscls = EulerReactingElements
    intinterscls = EulerReactingIntInters
    mpiinterscls = EulerReactingMPIInters
    bbcinterscls = EulerReactingBaseBCInters
