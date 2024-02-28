import os
import numpy as np

import ase.io
from ase.calculators.espresso import Espresso
from ase.io.espresso import write_espresso_in
from ase.io import extxyz
from ase.atoms import Atoms

from pathlib import Path
import re
import glob
from natsort import natsorted, ns

import itertools
import random
import math
import struct

#################################################################
# Convert global (fitted) coefficients to per-atom coefficients
#################################################################

f = open("HfO2FPOD_coefficients.pod", "r")
coeffs = f.read().splitlines()[1::]
f.close()

coeffs = [float(x) for x in coeffs]

Ne = 2
# 1-body
nd1 = Ne
# 2-body
nrbf2 = 12
# 3-body
nrbf3 = 8
three_body_angular_degree = 4
nabf3 = three_body_angular_degree + 1
# 4-body
nrbf4 = 5
four_body_angular_degree = 2
nb = [1,     2,     4,     7,    11,    16,    23]
nabf4 = nb[four_body_angular_degree]
nrbf23 = 0
nabf23 = 0

# 5-body
nrbf33 = 0
nabf33 = 0
# 6-body
nrbf34 = 0
nabf34 = 0
nabf43 = 0
# 7-body
nrbf44 = 0
nabf44 = 0
############
nl1 = int(nd1/Ne)
nl2 = nrbf2*Ne
nl3 = nabf3*nrbf3*Ne*(Ne+1)/2
nl4 = nabf4*nrbf4*Ne*(Ne+1)*(Ne+2)/6

nd2 = int(nrbf2*Ne*(Ne+1)/2)
nd3 = int(nabf3*nrbf3*Ne*Ne*(Ne+1)/2)
nd4 = int(nabf4*nrbf4*Ne*Ne*(Ne+1)*(Ne+2)/6)

n23 = nrbf23*Ne
n32 = nabf23*nrbf23*Ne*(Ne+1)/2
n33 = nabf33*nrbf33*Ne*(Ne+1)/2
n34 = nabf34*nrbf34*Ne*(Ne+1)/2
n43 = nabf43*nrbf34*Ne*(Ne+1)*(Ne+2)/6
n44 = nabf44*nrbf44*Ne*(Ne+1)*(Ne+2)/6

nl23 = n23*n32
nl33 = n33*(n33+1)/2
nl34 = n34*n43
nl44 = n44*(n44+1)/2

nd23 = nl23*Ne
nd33 = nl33*Ne
nd34 = nl34*Ne
nd44 = nl44*Ne

nl = nl1 + nl2 + nl3 + nl4 + nl23 + nl33 + nl34 + nl44
nd = nd1 + nd2 + nd3 + nd4 + nd23 + nd33 + nd34 + nd44


elemindex = np.zeros((Ne,Ne), dtype=np.uint)
for i in range(Ne):
    for j in range(Ne):
        elemindex[i,j] = i+j

elemindex = list(elemindex.flatten())

def twobodycoeff(newcoeffs, coeffs):
    for i in range(Ne):
        for j in range(Ne):
            for m in range(nrbf2):
                newcoeffs[nd1 + m + nrbf2*j + nrbf2*Ne*i] = coeffs[nd1 + m + nrbf2*elemindex[i + Ne*j]]
    
    return newcoeffs
    
def mknewcoeff(coeffs):
    newcoeffs = [None] * int(nl)*Ne
    for n in range(nd1):
        newcoeffs[n] = coeffs[n]
    
    if (nd2>0):
        twobodycoeff(newcoeffs, coeffs)

    for n in range(int(nd3+nd4+nd23+nd33+nd34+nd44)):
        newcoeffs[(nl1 + nl2)*Ne + n] = coeffs[nd1+nd2+n]

    return newcoeffs

newcoeffs = mknewcoeff(coeffs)

elemCoeffs = [[None]*int(nl) for _ in range(Ne)]

for k in range(len(Nl)):
    buffer = int(sum(Nl[:k]*Ne))
    for ne in range(Ne):
        for nli in range(Nl[k]):
            elemCoeffs[ne][int(buffer/Ne) + nli] = newcoeffs[buffer + Nl[k]*ne + nli]
    if k==4: break


