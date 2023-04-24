# Air & JET-A1 # GASES #

#####GASES######
M_H = 1 #g/mol
M_O = 16 #g/mol
M_N = 14 #g/mol
M_C = 12 #g/mol
#################

#Air Composition
nN = 0.78 * 2
nO = 0.22 * 2
M_air = 14*nN + 16*nO
#print('El aire está compuesto por', M_air, 'M')

Ru = 8.3144 # J/mol/K
R_air = Ru/M_air*1000 # J/kg/K

T0 = 288.15 # K (ISA)
p0 = 101325 # Pa (ISA)

delta_T = 0
rho0 = 1.225 # Kg/m^3
alpha = -6.5e-3 # ºC/m

#JET-A1 Composition (25% C8H18; 25% C7H16; 50% C14H30)################################

#C8H18
M_C8 = 8 * M_C + 18* M_H
#C7H16
M_C7 = 7 * M_C + 16 * M_H
#C14H30
M_C14 = 14 * M_C + 30 * M_H

f_est_C8 = M_C8 / ((8 + 18/4)*4.76*M_air)
f_est_C7 = M_C7 / ((7 + 16/4)*4.76*M_air)
f_est_C14 = M_C14 / ((14 + 30/4)*4.76*M_air)

f_est = 0.25 * f_est_C8 + 0.25 * f_est_C7 + 0.5 * f_est_C14

#####################################################################################

g = 9.81 #m/s^2
#####Conversión######
def ft_m(ft):
    m = ft / 3.2808
    return m
def kts_ms(kts):
    mps = 0.5144 * kts
    return mps
#####################

def T_isa(h, delta_T=0):
    T = T0 + alpha * h + delta_T
    return T
def p_isa(h):
    p = p0 * (1 + alpha * h/T0) ** (-g/R_air/alpha)
    return p
def cp (T):
    if T < 1000:
        a0 = 3.56839620e+00 # Independent coefficient [kg*K^2/J]
        a1 = -6.78729429e-04 # 1st order [kg/J/K]
        a2 = 1.55371476e-06 # 2nd order [kg/J/K**2]
        a3 = -3.29937060e-12 # 3rd order [kg/J/K**3]
        a4 = -4.66395387e-13 # 4th order [kg/J/K**4]

        return R_air * (a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4)
    #
    if T >= 1000:
        a0 = 3.08792717e+00 # Independent coefficient
        a1 = 1.24597184e-03 # 1st order [kg/J/K]
        a2 = -4.23718945e-07 # 2nd order [kg/J/K**2]
        a3 = 6.74774789e-11 # 3rd order [kg/J/K**3]
        a4 = -3.97076972e-15 # 4th order [kg/J/K**4]

        return R_air * (a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4)
    
def gamma (cp):
    gamma = cp/(cp-R_air)
    return gamma

def gamma_T (T):
    if T < 1000:
        a0 = 3.56839620e+00 # Independent coefficient [kg*K^2/J]
        a1 = -6.78729429e-04 # 1st order [kg/J/K]
        a2 = 1.55371476e-06 # 2nd order [kg/J/K**2]
        a3 = -3.29937060e-12 # 3rd order [kg/J/K**3]
        a4 = -4.66395387e-13 # 4th order [kg/J/K**4]

        cp = R_air * (a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4)
        gamma = (cp/(cp-R_air))
        
        return gamma
    if T >= 1000:
        a0 = 3.08792717e+00 # Independent coefficient
        a1 = 1.24597184e-03 # 1st order [kg/J/K]
        a2 = -4.23718945e-07 # 2nd order [kg/J/K**2]
        a3 = 6.74774789e-11 # 3rd order [kg/J/K**3]
        a4 = -3.97076972e-15 # 4th order [kg/J/K**4]

        cp = R_air * (a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4)
        gamma = cp/(cp-R_air)
        
        return gamma

def r_e(T,M):
    gamma = gamma_T(T)
    Tt_T = 1 +(gamma-1)/2*M**2
    return Tt_T
def T_r(T,M):
    Tt = T * r_e(T,M)
    return Tt
def p_r(p0,T,M):
    gamma = gamma_T(T)
    pt = p0 * (r_e(T,M))**(gamma/(gamma-1))
    return pt