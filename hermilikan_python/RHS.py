from http.client import error
import numpy as np
from fastar import *
from foll import *


def vectorfield(x0, t, kp, ki, kd):
    """ x0 : vector of state variables
        t  : time
    """

    # State vectors and control inputs
    x = np.array(x0).reshape(13, 1)

    nu = x[0:6]
    eta = x[6:12]

    # tegrud skekkja fyrir PID
    i_err = x[12, 0]


    # u1      # haegri vaengur
    # u2      # vinstri vaengur
    #u1 = 0.1
    #u2 = 0.2
    ##########

    """ if 37 <= t <= 150:
        sp = 5
    elif 150 <= t <= 191:
        sp = 10
    elif 191 <= t <= 236:
        sp = 15
    else:
        sp = 0 """

    """ sp = 5*(np.heaviside(t-37, 0.5)-np.heaviside(t-150, 0.5)) \
    + 10*(np.heaviside(t-150, 0.5)-np.heaviside(t-191, 0.5)) \
    + 15*(np.heaviside(t-191, 0.5)-np.heaviside(t-236, 0.5))  """
    
    # setpoint
    sp = 5
    
    # her er reiknad skekkjuna milli setpoint og alvoru gildi a phi
    err = sp*np.pi/180 - eta[3, 0]
    d_err = err - error_listi[-1]
    error_listi.append(err)
    # her er PID styringin
    styrimerki = kp*err+ki*i_err
    
    u1 = styrimerki
    u2 = u1
    ##########

    ######################## saturation fyrir vaeng ###########################
    max_afallshorn = 16

    u1 = u1*180/np.pi
    # u1min = -max_afallshorn*pi/180
    # u1max = max_afallshorn*pi/180
    #
    u2 = u2*180/np.pi
    # u2min = -max_afallshorn*pi/180
    # u2max = max_afallshorn*pi/180

    u1min = -max_afallshorn
    u1max = max_afallshorn
    u2min = -max_afallshorn
    u2max = max_afallshorn

    if u1 > u1max:
        u1 = u1max
    elif u1 < u1min:
        u1 = u1min

    if u2 > u2max:
        u2 = u2max
    elif u2 < u2min:
        u2 = u2min

    ##########################################################################
    # Ocean currents expressed in BODY

    u_c = Vc * np.cos(betaVc - eta[5, 0])
    v_c = Vc * np.sin(betaVc - eta[5, 0])

    nu_c = np.array([u_c, v_c, 0, 0, 0, 0]).reshape(
        6, 1)                # ocean current velocities
    Dnu_c = np.array([nu[5, 0]*v_c, -nu[5, 0]*u_c, 0, 0, 0, 0]
                     ).reshape(6, 1)   # time derivative of nu_c

    # Relative velocities/speed, angle of attack and vehicle speed
    nu_r = nu - nu_c                                 # relative velocity

    # np.sqrt(nu_r[0]**2 + nu_r[1]**2 + nu_r[2]**2 )  # relative speed (m/s)
    U_r = np.linalg.norm(nu_r)
    # np.sqrt(nu[0]**2 + nu[1]**2 + nu[2]**2 )         # speed (m/s)
    U = np.linalg.norm(nu)

    alpha = np.arctan2(nu_r[2, 0], nu_r[0, 0])   # angle of attack (rad)

    if U_r == 0:
        beta = 0
    else:
        beta = np.arcsin(nu_r[1, 0]/U_r)         # sideslip (rad)

    # dynamic pressure

    q = 1/2 * rho * U_r**2

    # Her koma massafylki og hverfitregður

    CRB = m2c(MRB, nu)
    CA = 0.15*CRB

    CA[4, 0] = 0
    CA[4, 3] = 0
    CA[5, 0] = 0
    CA[5, 1] = 0

    #M = MRB + MA
    C = CRB + CA
    # m = MRB[0,0]
    W = m*g_mu
    B = BP*W  # ef B = W tha "neutrally buoyant"

    # Low-speed linear damping matrix parameter (kafbatur)

    T1 = 20  # 20            # time constant in surge (s)
    T2 = 20  # 20            # time constant in sway (s)
    zeta4 = 8                # relative damping ratio in roll
    zeta5 = 1                # relative damping ratio in pitch
    T6 = 5                   # time constant in yaw (s)

    # Dissipative forces and moments
    D = Dmtrx([T1, T2, T6], [zeta4, zeta5], MRB, MA, [W, r_bg, r_bb])
    D[0, 0] = D[0, 0] * np.exp(-3*U_r)   # vanish at high speed where quadratic
    D[1, 1] = D[1, 1] * np.exp(-3*U_r)   # drag and lift forces dominates
    D[5, 5] = D[5, 5] * np.exp(-3*U_r)

    ################### straumfraedilegir kraftar a net #######################

    # gefum okkur ad midad vid flaedi se bara dragkraftur
    tau_liftdrag = np.array([-q*(A_varpa)*Cd, 0, 0]).reshape(3, 1)  # flow axis

    tau_liftdrag = C_bf(alpha, beta)@tau_liftdrag  # breyta i body axis

    # verkar i dragmidju
    M_fra_neti = np.cross(d_m, tau_liftdrag, axis=0)

    tau_liftdrag = np.vstack((tau_liftdrag, M_fra_neti))

    ###########################################################################

    # Kinematics (adal snuningsfylki)
    [J, R, T_tf] = eulerang(eta[3, 0], eta[4, 0], eta[5, 0])

    # Restoring forces and moments (thyngd, flotkraftur, massamidja, flotmidja)

    g = gRvect(W, B, R, r_bg, r_bb)

    ###########################################################################

    ################# kraftar og vaegi fra vaengjum (i b-frame) ###############
    ################# gera betur her )  ######################################

    # beta = 0

    # haegri vaengur relative speed
    v_haegri = np.cross(nu_r[3:6], hv, axis=0)+nu_r[0:3]
    # v_haegri[0:2]
    U_r_haegri = np.linalg.norm(v_haegri[0:2])
    q_hv = 1/2 * rho * U_r_haegri**2

    hv_alpha = u1 + alpha*180/np.pi  # heildar angle of attack haegri vaengur
    u_v = u1

    # vinstri vaengur relative speed
    v_vinstri = np.cross(nu_r[3:6], vv, axis=0)+nu_r[0:3]
    U_r_vinstri = np.linalg.norm(v_vinstri[0:2])
    q_vv = 1/2 * rho * U_r_vinstri**2

    vv_alpha = -u2 + alpha*180/np.pi  # heildar angle of attack vinstri vaengur

    CD_hv = 0.01+0.00113333*abs(hv_alpha)-0.00024 * \
        abs(hv_alpha)**2+0.0000426667*abs(hv_alpha)**3
    C_L_alph_hv = 0.1*hv_alpha

    CD_vv = 0.01+0.00113333*abs(vv_alpha)-0.00024 * \
        abs(vv_alpha)**2+0.0000426667*abs(vv_alpha)**3
    C_L_alph_vv = 0.1*vv_alpha

    F1 = np.array([-q_hv * A_vaengur * CD_hv,  # drag force
                   0,
                   -q_hv * A_vaengur * C_L_alph_hv]).reshape(3, 1)  # lift force

    F2 = np.array([-q_vv * A_vaengur * CD_vv,       # drag force
                   0,
                   -q_vv * A_vaengur * C_L_alph_vv]).reshape(3, 1)  # lift force

    # transformation fra FLOW til BODY

    F1 = C_bf(hv_alpha*np.pi/180, beta)@F1
    F2 = C_bf(vv_alpha*np.pi/180, beta)@F2

    Fv = F1 + F2

    # vaegi fra vaengjum (a ad vera um cg, center of gravity)

    # haegri
    lever_vh = hv

    # vinstri
    lever_vv = vv

    Mh = np.cross(lever_vh, F1, axis=0)
    Mv = np.cross(lever_vv, F2, axis=0)

    ############################## togvira kraftar ############################

    # stadsetning o_b
    o_b = eta[0:3]

    HO_X = o_b + R@HO
    HN_X = o_b + R@HN
    VO_X = o_b + R@VO
    VN_X = o_b + R@VN

    xyz_hnit = np.vstack((HO_X, HN_X, VO_X, VN_X))

    rammi_xyz_hnit = np.vstack(
        (o_b + R@R_HO, o_b + R@R_HN,  o_b + R@R_VO, o_b + R@R_VN))

    HO_U = R@(np.cross(nu_r[3:6], HO, axis=0)+nu_r[0:3])
    HN_U = R@(np.cross(nu_r[3:6], HN, axis=0)+nu_r[0:3])
    VO_U = R@(np.cross(nu_r[3:6], VO, axis=0)+nu_r[0:3])
    VN_U = R@(np.cross(nu_r[3:6], VN, axis=0)+nu_r[0:3])

    pos_tow = TOW_hradi*t

    pos_X = np.array([pos_tow, 0, 0]).reshape(3, 1) + Bo

    DeltaX1 = pos_X - HO_X
    DeltaX2 = pos_X - HN_X
    DeltaX3 = pos_X - VO_X
    DeltaX4 = pos_X - VN_X

    TOW_U = np.array([TOW_hradi, 0, 0]).reshape(3, 1)

    DeltaU1 = TOW_U - HO_U
    DeltaU2 = TOW_U - HN_U
    DeltaU3 = TOW_U - VO_U
    DeltaU4 = TOW_U - VN_U

    ####################### togkraftar i cable frame #########################
    if np.linalg.norm(DeltaX1) > L0_O:
        Tow_1 = np.array([KO*(np.linalg.norm(DeltaX1)-L0_O)+np.squeeze(dempun_togvir*(
            np.dot(DeltaU1.T, (DeltaX1/np.linalg.norm(DeltaX1))))), 0, 0]).reshape(3, 1)
    else:
        Tow_1 = np.array([0, 0, 0]).reshape(3, 1)

    if np.linalg.norm(DeltaX2) > L0_N:
        Tow_2 = np.array([KN*(np.linalg.norm(DeltaX2)-L0_N)+np.squeeze(dempun_togvir*(
            np.dot(DeltaU2.T, (DeltaX2/np.linalg.norm(DeltaX2))))), 0, 0]).reshape(3, 1)
    else:
        Tow_2 = np.array([0, 0, 0]).reshape(3, 1)

    if np.linalg.norm(DeltaX3) > L0_O:
        Tow_3 = np.array([KO*(np.linalg.norm(DeltaX3)-L0_O)+np.squeeze(dempun_togvir*(
            np.dot(DeltaU3.T, (DeltaX3/np.linalg.norm(DeltaX3))))), 0, 0]).reshape(3, 1)
    else:
        Tow_3 = np.array([0, 0, 0]).reshape(3, 1)

    if np.linalg.norm(DeltaX4) > L0_N:
        Tow_4 = np.array([KN*(np.linalg.norm(DeltaX4)-L0_N)+np.squeeze(dempun_togvir*(
            np.dot(DeltaU4.T, (DeltaX4/np.linalg.norm(DeltaX4))))), 0, 0]).reshape(3, 1)
    else:
        Tow_4 = np.array([0, 0, 0]).reshape(3, 1)

    ####### Cable to inertial rotation matrix reikningar #######
    gamma1 = np.arctan(DeltaX1[2, 0]/DeltaX1[0, 0])
    sigma1 = np.arcsin(DeltaX1[1, 0]/np.linalg.norm(DeltaX1))

    gamma2 = np.arctan(DeltaX2[2, 0]/DeltaX2[0, 0])
    sigma2 = np.arcsin(DeltaX2[1, 0]/np.linalg.norm(DeltaX2))

    gamma3 = np.arctan(DeltaX3[2, 0]/DeltaX3[0, 0])
    sigma3 = np.arcsin(DeltaX3[1, 0]/np.linalg.norm(DeltaX3))

    gamma4 = np.arctan(DeltaX4[2, 0]/DeltaX4[0, 0])
    sigma4 = np.arcsin(DeltaX4[1, 0]/np.linalg.norm(DeltaX4))

    RCableToInertial1 = R_I_tog(gamma1, sigma1)
    RCableToInertial2 = R_I_tog(gamma2, sigma2)
    RCableToInertial3 = R_I_tog(gamma3, sigma3)
    RCableToInertial4 = R_I_tog(gamma4, sigma4)

    #############################################################

    ######## kraftar i togvirum yfirfaerdir i body frame ########
    Tow_1 = R.T@(RCableToInertial1@Tow_1)
    Tow_2 = R.T@(RCableToInertial2@Tow_2)
    Tow_3 = R.T@(RCableToInertial3@Tow_3)
    Tow_4 = R.T@(RCableToInertial4@Tow_4)
    ################### vaegi reikningar #########################
    M_TOW_1 = np.cross(HO, Tow_1, axis=0)
    M_TOW_2 = np.cross(HN, Tow_2, axis=0)
    M_TOW_3 = np.cross(VO, Tow_3, axis=0)
    M_TOW_4 = np.cross(VN, Tow_4, axis=0)

    TOW = Tow_1 + Tow_2 + Tow_3 + Tow_4
    M_TOW = M_TOW_1 + M_TOW_2 + M_TOW_3 + M_TOW_4

    ##### thyngd togvira #####
    mo = 0.98*L0_O/2  # deili med 2, skipt a milli tp og ramma
    mn = 0.98*L0_N/2

    thyngd_o = R.T@np.array([0, 0, mo*g_mu]).reshape(3, 1)
    thyngd_n = R.T@np.array([0, 0, mn*g_mu]).reshape(3, 1)

    TOW = TOW + 2*thyngd_o + 2*thyngd_n

    m_HO_th = np.cross(HO, thyngd_o, axis=0)
    m_HN_th = np.cross(HN, thyngd_n, axis=0)
    m_VO_th = np.cross(VO, thyngd_o, axis=0)
    m_VN_th = np.cross(VN, thyngd_n, axis=0)

    M_TOW = M_TOW + m_HO_th + m_HN_th + m_VO_th + m_VN_th
    #################### allt plusad saman her ################################

    F_heild = Fv + TOW

    M_heild = Mh + Mv + M_TOW

    tau = np.vstack((F_heild, M_heild))

    ##########################################################################

    # Generalized force vector

    allt = tau + tau_liftdrag - C@nu_r - D@nu_r - g

    # State-space model
    xdot = np.vstack((Dnu_c + inv_M@allt, J@nu, err)).flatten().tolist()

    return xdot
