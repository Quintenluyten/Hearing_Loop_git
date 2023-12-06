#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 20:31:27 2023

@author: anton
"""

from SLiCAP import *

ini.factor=False

prj = initProject('MOS noise Rs')

i1 = instruction()

i1.setCircuit('MOSnoiseR.cir')

i1.setSimType('numeric')

##############################################################################
f      = sp.Symbol('f', positive=True)
g_m    = sp.Symbol('g_m', positive=True)
c_iss  = sp.Symbol('c_iss', positive=True)
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

f_max = sp.Symbol('f_max', positive=True)
f_min = sp.Symbol('f_min', positive=True)
fmax  = 6E3
fmin  = 20

substDict = {f_max: fmax,
             f_min: fmin,
             alpha: 0.0033,
             f_T: 50E9,
             Gamma: Gamman,
             n: nn,
             T: Tn,
             k: kn}
##############################################################################

i1.setSimType('symbolic')
i1.setGainType('vi')
i1.setDataType('noise')
i1.setSource('V1')
i1.setDetector('V_out')
noiseResult = i1.execute()

symOnoise1 = assumePosParams(noiseResult.onoise)

# Study the result of different substitutions:
i1.setSimType('numeric')

symOnoise2 = sp.simplify(symOnoise1.subs(f_T, g_m/(2*sp.pi*c_iss)))

symOnoise3 = sp.expand(symOnoise1.subs(f_L, alpha*f_T))

# The last expression shows there must be an optimum g_m:

diff_SymOnoise3_g_m = sp.diff(symOnoise3, g_m)
g_m_opt = sp.solve(diff_SymOnoise3_g_m , g_m)[0]

# At this optimum the output noise equals:
symOnoise3 = symOnoise3
print("bip")
S_opt = sp.simplify(symOnoise3.subs(g_m, g_m_opt))
print("boep")
# Calculate the noise figure
SRCnoise = assumePosParams(noiseResult.onoiseTerms['I_noise_R1'])

try:
    NF_min = sp.integrate(S_opt, f), fmin, fmax/sp.integrate(SRCnoise, f, fmin, fmax)
except:
    # Numeric integration
    S_opt_num = S_opt.subs(substDict)
    SRCnoise_num = SRCnoise.subs(substDict)
    S_opt_num = S_opt_num/(SRCnoise_num*(fmax-fmin))
    S_opt_n_func = sp.lambdify(f, S_opt_num)
    NF_min = quad(S_opt_n_func, fmin, fmax)[0]

NFdBmin    = sp.N(10*sp.log(NF_min)/sp.log(10),2)

print("Minimum noise figure:", NFdBmin, "dB")