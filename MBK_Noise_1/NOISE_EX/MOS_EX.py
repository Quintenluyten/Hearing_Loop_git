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
R_s  = sp.Symbol('R_s', positive=True)
R_sn = i1.getParValue('R_s')
R_t  = sp.Symbol('R_t', positive=True)
R_tn = i1.getParValue('R_t')
L_s  = sp.Symbol('L_s', positive=True)
L_sn = i1.getParValue('L_s')
n      = sp.Symbol('n', positive=True)
nn     = i1.getParValue('n')
k      = sp.Symbol('k', positive=True)
kn     = i1.getParValue('k')
T      = sp.Symbol('T', positive=True)
Tn     = i1.getParValue('T')
tau_i    = sp.Symbol('tau_i', positive=True)
tau_i_n    = i1.getParValue('tau_i')
c_iss  = sp.Symbol('c_iss', positive=True)

f_max = sp.Symbol('f_max', positive=True)
f_min = sp.Symbol('f_min', positive=True)
fmax  = 6e3
fmin  = 20


substDict_inclkT = {f_max: fmax,
             f_min: fmin,
             alpha: 0.0033,
             f_T: 50E9,
             Gamma: Gamman,
             n: nn,
             R_s: R_sn,
             R_t: R_tn,
             L_s: L_sn,
             T: Tn,
             k: kn,
             tau_i: tau_i_n}

substDict_exclkT = {f_max: fmax,
             f_min: fmin,
             alpha: 0.0033,
             f_T: 50E9,
             Gamma: Gamman,
             n: nn,
             R_s: R_sn,
             R_t: R_tn,
             L_s: L_sn,
             tau_i: tau_i_n}

Showstopper=10.6e-6

i1.setSimType('symbolic')
i1.setGainType('vi')
i1.setDataType('noise')
i1.setSource('V1')
i1.setDetector('V_out')
noiseResult = i1.execute()


symOnoise1 = assumePosParams(noiseResult.onoise)
htmlPage("Symbolic noise analysis")
noise2html(noiseResult, label='symNoise')


symOnoise1 = assumePosParams(noiseResult.onoise)

symOnoise11 = sp.simplify(symOnoise1.subs(substDict_exclkT))
htmlPage("Simplified symbolic expression")
eqn2html('symOnoise11', symOnoise11, label = 'symplified noise', labelText = 'test')

i1.setSimType('numeric')




symOnoise2 = sp.simplify(symOnoise11.subs(f_T, g_m/(2*sp.pi*c_iss)))
symOnoise3 = sp.simplify(symOnoise11.subs(f_L, alpha*f_T))


diff_SymOnoise3_g_m = sp.diff(symOnoise3, g_m)
g_m_opt = sp.solve(diff_SymOnoise3_g_m , g_m)[0]

#istn't the showstopper value in RMS and sym0noise is in V^2/hz ?
htmlPage("Numeric noise analysis")
eqn2html('symOnoise1', symOnoise1, label = 'noise', labelText = 'noise')
eqn2html('symOnoise2', symOnoise2, label = 'noise', labelText = 'noise')
eqn2html('symOnoise3', symOnoise3, label = 'noise', labelText = 'noise')
eqn2html('g_m_opt', g_m_opt, label = 'noise', labelText = 'noise')


symOnoise3 = symOnoise3
S_opt = symOnoise3.subs(g_m, g_m_opt)

eqn2html('S_opt', S_opt, label = 'noise', labelText = 'noise')




try:
    RMS = sp.integrate(S_opt, f), fmin, fmax
except:
    # Numeric integration
    S_opt_num = S_opt.subs(substDict)
    S_opt_n_func = sp.lambdify(f, S_opt_num)
    RMS = quad(S_opt_n_func, fmin, fmax)[0]


eqn2html('RMS', RMS, label = 'noise', labelText = 'noise')

