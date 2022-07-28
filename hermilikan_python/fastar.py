import numpy as np

b = 6  # breidd
h = 3  # haed

m_e = 80  # massi einnar einingar af ramma
aukamassi = 150  # massi frá öðru en rammanum sjálfum

Ixx = 2*(1/12*((b/3)*m_e)*b**2+((b/3)*m_e)*(h/2)**2) + \
    2*(1/12*(b/3)*m_e*h**2+(h/3)*m_e*(b/2)**2)

Iyy = 2*(1/12*(h/3)*m_e*h**2)+2*(((b/3)*m_e)*(h/2)**2)

Izz = 2*(1/12*((b/3)*m_e)*b**2)+2*((h/3)*m_e*(b/2)**2)

Ig = np.array([[Ixx, 0, 0, 0, Iyy, 0, 0, 0, Izz]]).reshape(3, 3)


m = 2*(b/3)*m_e + 2*(h/3)*m_e + aukamassi

MRB = np.block([[m*np.eye(3), np.zeros((3, 3))],
                [np.zeros((3, 3)), Ig]])

MA = 0.15*MRB  # + 15 # added mass

M = MRB + MA

inv_M = np.linalg.inv(M)

breidd = b
haed = h
A_varpa = breidd*haed  # flatarmal op vorpu m**2

# fleiri fastar
mu = 66         # Lattitude for Trondheim, Norway (deg)
g_mu = 9.82   # gravity vector (m/s2)
rho = 1026             # density of water (kg/m3)

##########################################################################

# togvir fastar
E = 105*10**9            # Youngs studull Pa (stal cable)
d_vir = 0.05             # thvermal togvir m
Av = np.pi*(d_vir/2)**2  # flatarmal togvir m**2
rho_togvir = 7800        # kg/m**3
k_damp = 1.5             # fra modeling and control of trawl systems
dempun_togvir = k_damp*Av*np.sqrt(E*rho_togvir)
#dempun_togvir = 85000

corr_f = 0  # -0.5    # stillir adeins upphafsstellingu vorpu

# rammi hnit
# haegri ofan
R_HO = np.array([corr_f, -breidd/2, -haed/2]).reshape(3, 1)
# haegri nedan
R_HN = np.array([0, -breidd/2, haed/2]).reshape(3, 1)
# vinstri ofan
R_VO = np.array([corr_f, breidd/2, -haed/2]).reshape(3, 1)
# vinstri nedan
R_VN = np.array([0, breidd/2, haed/2]).reshape(3, 1)


# vaengir
L_vaengur = 2.2  # m
b_vaengur = 0.27  # m
A_vaengur = L_vaengur*b_vaengur  # m**2


# haegri vaengur fjarlaegd fra CO (center of body coordinate system)
hv = np.array([corr_f, -breidd/2+0.4+1.1, -haed/2]).reshape(3, 1)

# vinstri vaengur fjarlaegd fra CO
vv = np.array([corr_f, breidd/2-0.4-1.1, -haed/2]).reshape(3, 1)

################## stadsetningar togvira #####################


# stadsetning TP i byrjun m.v. o_b
Bo = np.array([9, 0, -4]).reshape(3, 1)

# haegri ofan
HO = np.array([corr_f, -breidd/4, -haed/2]).reshape(3, 1)
# haegri nedan
HN = np.array([0, -breidd/2, 0]).reshape(3, 1)
# vinstri ofan
VO = np.array([corr_f, breidd/4, -haed/2]).reshape(3, 1)
# vinstri nedan
VN = np.array([0, breidd/2, 0]).reshape(3, 1)

# upphafleg lengd efri togvira
L0_O = np.linalg.norm(Bo-HO)

# upphafleg lengd nedri togvira
L0_N = np.linalg.norm(Bo-HN)

# togvira stifni
KO = 1.82*10**7  # E*Av/L0_O
KN = 1.82*10**7  # E*Av/L0_N

####### straumfraedi studlar ######

Cd = 0.45     # dragstudulll net og co.


TOW_hradi = 1.25        # m/s i +x stefnu

d_m = np.array([-2, 0, 0]).reshape(3, 1)  # dragmidja
r_bg = np.array([0, 0, 0]).reshape(3, 1)          # CG w.r.t. to the CO
r_bb = np.array([0, 0, -haed/2]).reshape(3, 1)    # CB w.r.t. to the CO

BP = 0.612              # hlutfall flotkraftur/thyngd

Vc = 0                    # straumhradi
betaVc = 0                # horn straums vorpu

########  (initial conditions) ########
x0 = np.array([TOW_hradi,  # u0
               0,  # v0
               0,  # w0
               0,  # p0
               0,  # q0
               0,  # r0
               0,  # x0
               0,  # y0
               0,  # z0
               0,  # phi0
               0,  # theta
               0,  # psi0
               0]).reshape(13, 1) # error integration

ui = np.array([0*np.pi/180,  # u1
               0*np.pi/180]).reshape(2, 1)  # u2

Us = 0

roll = 0*np.pi/180

kp = 1
ki = 0
kd = 0
N = 1  # filter coefficient

error_listi = [0]