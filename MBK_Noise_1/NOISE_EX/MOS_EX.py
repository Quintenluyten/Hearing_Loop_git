# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 10:32:19 2023

@author: moelb
"""

from SLiCAP import *

ini.factor=False

prj = initProject('MOS noise EX')

i1 = instruction()

i1.setCircuit('noise_1.cir')
i1.setSimType('numeric')



f      = sp.Symbol('f', positive=True)
g_m    = sp.Symbol('g_m', positive=True)

f_T    = sp.Symbol('f_T', positive=True)
alpha  = sp.Symbol('alpha', positive=True)
f_L    = sp.Symbol('f_ell', positive=True)
Gamma  = sp.Symbol('Gamma', positive=True)
Gamman = i1.getParValue('Gamma')
n      = sp.Symbol('n', positive=True)
nn     = i1.getParValue('n')
k      = sp.Symbol('k', positive=True)
kn     = i1.getParValue('k')
T      = sp.Symbol('T', positive=True)
Tn     = i1.getParValue('T')
f_T    = sp.Symbol('tau_i', positive=True)
c_iss  = sp.Symbol('c_iss', positive=True)

f_max = sp.Symbol('f_max', positive=True)
f_min = sp.Symbol('f_min', positive=True)
fmax  = 6e3
fmin  = 20


substDict = {f_max: fmax,
             f_min: fmin,
             alpha: 0.0033,
             f_T: 50E9,
             Gamma: Gamman,
             n: nn,
             T: Tn,
             k: kn}

i1.setSimType('symbolic')
i1.setGainType('vi')
i1.setDataType('noise')
i1.setSource('V1')
i1.setDetector('V_out')
noiseResult = i1.execute()
Showstopper=10.6e-6

symOnoise1 = assumePosParams(noiseResult.onoise)
htmlPage("Symbolic noise analysis")
noise2html(noiseResult, label='symNoise')

i1.setSimType('numeric')
noiseResult = i1.execute()
symOnoise1 = assumePosParams(noiseResult.onoise)

symOnoise2 = sp.simplify(symOnoise1.subs(f_T, g_m/(2*sp.pi*c_iss)))
symOnoise3 = sp.simplify(symOnoise1.subs(f_L, alpha*f_T))

#istn't the showstopper value in RMS and sym0noise is in V^2/hz ?
htmlPage("Numeric noise analysis")
eqn2html('symOnoise1', symOnoise1, label = 'noise', labelText = 'noise')
eqn2html('symOnoise2', symOnoise2, label = 'noise', labelText = 'noise')
eqn2html('symOnoise3', symOnoise3, label = 'noise', labelText = 'noise')







