# # -*- coding: utf-8 -*-
from __future__ import unicode_literals

import sys
import random
import time

sys.path.append('.')
from GF5 import *
from IntegerFunctions import *

# Force UTF-8 output: console, pipe
import codecs
sys.stdout = codecs.getwriter('utf8')(sys.stdout)

print('\u03b1\u03b2\u03b3\u03b4\u03b5\u03b6\u03b7\u03b8')
print('\u03b9\u03ba\u03bb\u03bc\u03bd\u03be\u03bf\u03c1')
print('\u03c2\u03c3\u03c4\u03c5\u03c6\u03c7\u03c8\u03c9')
print('\N{greek capital letter gamma} \N{greek capital letter psi}')
print('\u2713\u2718\N{heavy check mark}\N{heavy ballot X}')
print('')
 
start = time.time()

F7 = IntegersModP(7)
print('5/6 = %s (mod 7)' % F7.div(5, 6))
print('')

GF(IntegersModP(3), 3, 'u', 'random').info()

F2 = IntegersModP(2)
GF4 = GF(F2, 2, '\N{greek small letter phi}', [1,1,1])
GF4.info()
GF4_4_4 = GF(GF4, 2, '\N{greek small letter mu}', [(0,1),(0,1),(1,0)])
GF4_4_4.info('extension')

def field2pds(gf, q):
    pds =[]
    for i in range(1, q *(q +1) +1):
        if gf.alog[i][0] == gf.cf.zero:
            pds.append(i -1)
    return pds

def base(b, n):
    # express n in base b form, least sig digit at list index 0
    # e.g. 14, 15, 16 base 3 -- [2,1,1] [0,2,1] [1,2,1]
    # for walking through all reduced polynomial terms
    digits=[]
    while n >0:
        digits.append(n % b)
        n = n //b
    return digits

def genPDS(q):
    pf = primeFactors(q)
    assert pf[0] ==pf[-1]
    p = pf[0]
    d = len(pf)
    if d ==1:
        gf = GF(IntegersModP(p), 3, 'x', 'random')
    else:
        gf = GF(GF(IntegersModP(p), d), 3)
    pds = field2pds(gf, q)
    print('genPDS: q=%d, %s' % (q, str(pds)))

def genPDS2(q):
    pf = primeFactors(q)
    assert pf[0] ==pf[-1]
    p = pf[0]
    d = len(pf)
    if d ==1:
        cf = IntegersModP(p)
    else:
        cf = GF(IntegersModP(p), d)
    for n in range(cf.order **3):
        terms = (tuple(cf[d] for d in base(cf.order, n)) +(cf[0], cf[0], cf[0],))[:3] +(cf[1],)
        try:
            g = GF(cf, 3, 'x', terms)
        except:
            continue
        print('%s/%s' % (cf, terms))
        print('   %s' %(str(field2pds(g, q))))

print('genPDS()')
for q in [2, 3, 4, 5, 17, 25]: #7, 17, 49, 81, 125]:
    genPDS(q)
print('')

print('genPDS2() list all GF of a given order and the derived PDS')
for q in [2, 3, 4, 5, 7]:
    genPDS2(q)
    print('')

g = GF(IntegersModP(2), 2, 'x', (1,1,1))
g.info()

print('Done.  %.2fsec'% (time.time() -start))
