import numpy as n
from fastar import *
from foll import *

#from numba import jit

#@jit
def varpa_dynamics(t, x0, u, params):
    """ x0 : vector of state variables 12x1
    nu
    #  u:       surge velocity          (m/s)
    #  v:       sway velocity           (m/s)
    #  w:       heave velocity          (m/s)
    #  p:       roll rate               (rad/s)
    #  q:       pitch rate              (rad/s)
    #  r:       yaw rate                (rad/s)
    
    eta
    #  x:       North position          (m)
    #  y:       East position           (m)
    #  z:       downwards position      (m)
    #  phi:     roll angle              (rad)       
    #  theta:   pitch angle             (rad)
    #  psi:     yaw angle               (rad)

        t  : time
    """
    # State vectors and control inputs
    x = np.array(x0).reshape(12, 1)

    nu = x[0:6]
    eta = x[6:12]

    # u1      # haegri vaengur
    # u2      # vinstri vaengur

    u1 = u*180/np.pi
    u2 = u*180/np.pi

    max_afallshorn = 16
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

    u_c = straumhradi * np.cos(sidehorn_straums_wrt_co - eta[5, 0])
    v_c = straumhradi * np.sin(sidehorn_straums_wrt_co - eta[5, 0])

    nu_c = np.array([u_c, v_c, 0, 0, 0, 0]).reshape(6, 1)  # ocean current velocities
    Dnu_c = np.array([nu[5, 0]*v_c, -nu[5, 0]*u_c, 0, 0, 0, 0]).reshape(6, 1) # time derivative of nu_c

    # Relative velocities/speed, angle of attack and vehicle speed
    nu_r = nu - nu_c                                 # relative velocity


    U_r = np.linalg.norm(nu_r) # relative speed (m/s)
    #U = np.linalg.norm(nu)

    alpha = np.arctan2(nu_r[2, 0], nu_r[0, 0])   # angle of attack (rad)

    if U_r == 0:
        beta = 0
    else:
        beta = np.arcsin(nu_r[1, 0]/U_r)         # sideslip (rad)

    # dynamic pressure
    q = 1/2 * rho_vatn * U_r**2

    # Her koma massafylki og hverfitregður
    CRB = m2c(MRB, nu)
    CA = 0.15*CRB

    CA[4, 0] = 0.0
    CA[4, 3] = 0.0
    CA[5, 0] = 0.0
    CA[5, 1] = 0.0

    C = CRB + CA
    W = massi_heild*gravity_a
    B = ratio_Buoy_weight*W  # ef B = W tha "neutrally buoyant"

    # Dissipative forces and moments
    D = Dmtrx([T1, T2, T6], [zeta4, zeta5], MRB, MA, [W, CGravity_wrt_CO, CBuoyancy_wrt_CO])
    D[0, 0] = D[0, 0] * np.exp(-3*U_r)   # vanish at high speed where quadratic
    D[1, 1] = D[1, 1] * np.exp(-3*U_r)   # drag and lift forces dominates
    D[5, 5] = D[5, 5] * np.exp(-3*U_r)

    ################### straumfraedilegir kraftar a net #######################

    # gefum okkur ad midad vid flaedi se bara dragkraftur (-x att i flow axis)
    tau_liftdrag = np.array([-q*(flatarmal_op_vorpu)*Cd_net_og_rest, 0, 0]).reshape(3, 1)  # flow axis

    # breyta i body axis
    tau_liftdrag = C_bf(alpha, beta)@tau_liftdrag

    # straumfraedilegi krafturinn latinn verka i dragmidju, thad kemur thvi vaegi fra honum
    M_fra_neti = np.cross(dragmidja_b, tau_liftdrag, axis=0)

    # straumfraedilegi kraftur og vaegi sett i einn vigur
    tau_liftdrag = np.vstack((tau_liftdrag, M_fra_neti))

    ###########################################################################

    # Kinematics (adal snuningsfylki)
    [J, R, T_tf] = eulerang(eta[3, 0], eta[4, 0], eta[5, 0])

    # Restoring forces and moments (thyngd, flotkraftur, massamidja, flotmidja)
    g = gRvect(W, B, R, CGravity_wrt_CO, CBuoyancy_wrt_CO)

    ###########################################################################

    ################# kraftar og vaegi fra vaengjum (i b-frame) ###############

    # haegri vaengur relative speed
    v_haegri = np.cross(nu_r[3:6], haegri_vaengur_stadsetning_b, axis=0)+nu_r[0:3]

    U_r_haegri = np.linalg.norm(v_haegri[0:2])
    q_hv = 1/2 * rho_vatn * U_r_haegri**2

    hv_alpha = u1 + alpha*180/np.pi  # heildar angle of attack haegri vaengur
    u_v = u1

    # vinstri vaengur relative speed
    v_vinstri = np.cross(nu_r[3:6], vinstri_vaengur_statsetning_b, axis=0)+nu_r[0:3]
    U_r_vinstri = np.linalg.norm(v_vinstri[0:2])
    q_vv = 1/2 * rho_vatn * U_r_vinstri**2

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

    Mh = np.cross(haegri_vaengur_stadsetning_b, F1, axis=0)
    Mv = np.cross(vinstri_vaengur_statsetning_b, F2, axis=0)


    ############################## togvira kraftar ############################

    # stadsetning o_b (origin of body coordinate system in inertial frame)
    o_b = eta[0:3]
    pos_tow = TOW_hradi*t
    pos_X = np.array([pos_tow, 0, 0]).reshape(3, 1) + stadsetning_hanafots_i_byrjun_b
    TOW_U = np.array([TOW_hradi, 0, 0]).reshape(3, 1)

    ALL_TOW_FORCES = np.zeros(3).reshape(3,1)
    ALL_TOW_MOMENTS = np.zeros(3).reshape(3,1)
    for i in range(len(TOGVIRAR)):
        # TOGVIRAR[i][0] stadsetning togvirs vid ramma i body-frame
        # TOGVIRAR[i][i] upphafleg lengd togvirs

        ### stadsetning og hradi tengipunkt togvirs a ramma
        H_X = o_b + R@(TOGVIRAR[i][0])
        H_U = R@(np.cross(nu_r[3:6], TOGVIRAR[i][0], axis=0)+nu_r[0:3])
        
        # munur a stadsetningu togvirs a ramma og hanastels
        DeltaX = pos_X - H_X
        # munur a hrada togvirs a ramma og hanastels
        DeltaU = TOW_U - H_U

        # reikna krafta fra togvirum
        Tow = togkraftur(DeltaX,DeltaU,togvir_stifni,dempun_togvir,TOGVIRAR[i][1])
        gamma, sigma = gamma_sigma(DeltaX)
        RCableToInertial = R_I_tog(gamma, sigma)

        ######## kraftar i togvirum yfirfaerdir i body frame ########
        Tow = R.T@(RCableToInertial@Tow)
        M_Tow = np.cross(TOGVIRAR[i][0], Tow, axis=0)

        ##### thyngd togvira #####
        m_togvir = 0.98*TOGVIRAR[i][1]/2  # deili med 2, skipt a milli tp og ramma
        thyngd_togvir = R.T@np.array([0, 0, m_togvir*gravity_a]).reshape(3, 1)
        vaegi_thyngd_togvir = np.cross(TOGVIRAR[i][0], thyngd_togvir, axis=0)

        ALL_TOW_FORCES += Tow + thyngd_togvir
        ALL_TOW_MOMENTS += M_Tow + vaegi_thyngd_togvir
    
    
    #################### allt plusad saman her ################################

    F_heild = Fv + ALL_TOW_FORCES

    M_heild = Mh + Mv + ALL_TOW_MOMENTS

    tau = np.vstack((F_heild, M_heild))

    ##########################################################################

    # Generalized force vector

    allt = tau + tau_liftdrag - C@nu_r - D@nu_r - g

    # State-space model
    xdot = np.vstack((Dnu_c + inv_M@allt, J@nu)).flatten()

    return xdot
