#!/usr/bin/python
import numpy
import math

## Fiber Parameters
D = 16.7 *pow(10,-9)/(pow(10,-9)*pow(10,3)) # Local dispersion parameter
gamma = 1.3 *pow(10,-3) # Nonlinear coefficient
a = pow(10,0.22/10)*pow(10,-3) # Attenuation parameter
c = 299792458 # Speed of light

## System Parameters
N_slots = 31 # Number of slots
L = 100 *pow(10,3) # Span length
L_eff = (1-math.exp(-a*L))/a # Effective length
N_s = 10 # Number of spans
F_c = 193.4 *pow(10,12) # Central frequency
B_slot = 12.5 *pow(10,9) # Bandwidth of slot
CR = 0 # Dispersion management scheme

## Calculations
Freqs = numpy.linspace(F_c-B_slot*(N_slots-1)/2, \
        F_c+B_slot*(N_slots-1)/2, N_slots)

lambd = [c/index for index in Freqs]

P1 = []
for index in range(1,N_slots+1):
    x1 = -(pow(lambd[index],2)*pow(B_slot,2)*D*L*(1-CR)*pow(N_slots,2))/(16*c)
    x2 = (pow(lambd[index],2)*pow(B_slot,2)*D*L*(1-CR)*pow(N_slots,2))/(2*c)
    x1 = int(math.floor(x1))
    x2 = int(math.floor(x2))

    s = 0
    for k in range(x1,x2):
        s = s + 1/(1+pow(2*k*math.pi/a*L*(1-CR),2))

    P1.append(4*pow(gamma,2)*pow(L_eff,2)*pow(N_s,2)*c) \
    /(pow(lambd[index],2)*pow(B_slot,2)*D*math.sqrt(pow(2/a,2)+ \
    (2*pow(L,2)*pow(1-CR,2)*(pow(N_s,2)-1))/pow(s,2)))

## Results
