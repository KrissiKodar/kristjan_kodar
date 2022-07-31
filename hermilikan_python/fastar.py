import numpy as np
breidd = 6.0  # breidd ramma
haed = 3.0    # haed ramma

m_eining = 80.0    # massi einnar einingar af ramma
aukamassi = 150.0  # massi frá öðru en rammanum sjálfum


### ath aukamassinn er ekki med (ma baeta vid)
Ixx = 2*(1/12*((breidd/3)*m_eining)*breidd**2+((breidd/3)*m_eining)*(haed/2)**2) + \
    2*(1/12*(breidd/3)*m_eining*haed**2+(haed/3)*m_eining*(breidd/2)**2)

Iyy = 2*(1/12*(haed/3)*m_eining*haed**2)+2*(((breidd/3)*m_eining)*(haed/2)**2)

Izz = 2*(1/12*((breidd/3)*m_eining)*breidd**2)+2*((haed/3)*m_eining*(breidd/2)**2)

hverfitregduthinur = np.array([[Ixx, 0, 0, 0, Iyy, 0, 0, 0, Izz]]).reshape(3, 3)


massi_heild = 2*(breidd/3)*m_eining + 2*(haed/3)*m_eining + aukamassi

MRB = np.block([[massi_heild*np.eye(3), np.zeros((3, 3))],
                [np.zeros((3, 3)), hverfitregduthinur]])

MA = 0.15*MRB  # + 15 # added mass

M = MRB + MA

inv_M = np.linalg.inv(M)


flatarmal_op_vorpu = breidd*haed  # flatarmal op vorpu m**2

# fleiri fastar
gravity_a = 9.82   # gravitational acceleration (m/s2)
rho_vatn = 1026.0             # density of water (kg/m3)

##########################################################################

# togvir fastar
togvir_young_modulus = 105*10**9   # Youngs studull Pa (stal cable)
tvhermal_togvir = 0.05             # thvermal togvir m
thverskurdsflatarmal_togvirs = np.pi*(tvhermal_togvir/2)**2  # flatarmal togvir m**2
rho_togvir = 7800.0      # edlismassi togvir kg/m**3
k_damp = 1.5             # fra modeling and control of trawl systems
dempun_togvir = k_damp*thverskurdsflatarmal_togvirs*np.sqrt(togvir_young_modulus*rho_togvir)
#dempun_togvir = 85000

corr_f = 0.0  # -0.5    # stillir adeins upphafsstellingu vorpu

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
lengd_vaengur = 2.2  # m
breidd_vaengur = 0.27  # m
A_vaengur = lengd_vaengur*breidd_vaengur  # m**2


# haegri vaengur fjarlaegd fra CO (center of body coordinate system)
haegri_vaengur_stadsetning_b = np.array([corr_f, -breidd/2+0.4+1.1, -haed/2]).reshape(3, 1)

# vinstri vaengur fjarlaegd fra CO
vinstri_vaengur_statsetning_b = np.array([corr_f, breidd/2-0.4-1.1, -haed/2]).reshape(3, 1)

################## stadsetningar togvira #####################


# stadsetning TP i byrjun m.v. o_b
stadsetning_hanafots_i_byrjun_b = np.array([9, 0, -4]).reshape(3, 1)

# haegri ofan
stadsetning_togvir_haegri_ofan_b = np.array([corr_f, -breidd/4, -haed/2]).reshape(3, 1)
# haegri nedan
stadsetning_togvir_haegri_nedan_b = np.array([0, -breidd/2, 0]).reshape(3, 1)
# vinstri ofan
stadsetning_togvir_vinstri_ofan_b = np.array([corr_f, breidd/4, -haed/2]).reshape(3, 1)
# vinstri nedan
stadsetning_togvir_vinstri_nedan_b = np.array([0, breidd/2, 0]).reshape(3, 1)

# upphafleg lengd efri togvira
upphafleg_L_efri_togvira = np.linalg.norm(stadsetning_hanafots_i_byrjun_b-stadsetning_togvir_haegri_ofan_b)

# upphafleg lengd nedri togvira
upphafleg_L_nedri_togvira = np.linalg.norm(stadsetning_hanafots_i_byrjun_b-stadsetning_togvir_haegri_nedan_b)

# togvira stifni
efri_togvir_stifni = 1.82*10**7  # E*Av/L0_O
nedri_togvir_stifni = 1.82*10**7  # E*Av/L0_N

####### straumfraedi studlar ######

Cd_net_og_rest = 0.45     # dragstudulll net og co.


#TOW_hradi = 1.25        # m/s i +x stefnu

dragmidja_b = np.array([-2, 0, 0]).reshape(3, 1)  # dragmidja
CGravity_wrt_CO = np.array([0, 0, 0]).reshape(3, 1)          # CG w.r.t. to the CO
CBuoyancy_wrt_CO = np.array([0, 0, -haed/2]).reshape(3, 1)    # CB w.r.t. to the CO

ratio_Buoy_weight = 0.612          # hlutfall flotkraftur/thyngd

straumhradi = 0.0                    # straumhradi
sidehorn_straums_wrt_co = 0.0        # horn straums vorpu


# Low-speed linear damping matrix parameter (kafbatur)
T1 = 20.0  # 20            # time constant in surge (s)
T2 = 20.0  # 20            # time constant in sway (s)
zeta4 = 8.0                # relative damping ratio in roll
zeta5 = 1.0                # relative damping ratio in pitch
T6 = 5.0                   # time constant in yaw (s)

########  (initial conditions) ########
""" x0 = np.array([TOW_hradi,  # u0
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
               0  # psi0
               ]).reshape(12, 1) # error integration """

TOW_hradi = 1.25