import numpy as np

"""
Morg af follunum her ad nedan eru fengin beint fra Thor I. Fossen, eda Marine System Simulator
fyrir Matlab
https://uk.mathworks.com/matlabcentral/fileexchange/86393-marine-systems-simulator-mss
"""


def Smtrx(a):
    """ 
    S = Smtrx(a) computes the 3x3 vector skew-symmetric matrix S(a) = -S(a)'.
    The corss product satisfies: a x b = S(a)b. The inverse can be computed
    as: vex(S(a)) = a 
    """
    a = a.flatten()
    S = np.array([   0,    -a[2],    a[1],
                  a[2],        0,   -a[0],
                 -a[1],     a[0],       0 ]).reshape(3,3)
    
    return S
 
def m2c(M,nu):
    """ 
    C = m2c(M,nu) computes the Coriolis-centripetal matrix C(nu) from the
    the system inertia matrix M > 0 for varying velocity nu. 
    If M is a 6x6 matrix and nu = [u, v, w, p, q, r]', the output is a 6x6 C matrix

    Examples: CRB = m2c(MRB,nu)     
            CA  = m2c(MA, nu)
    Output:
    C         - Coriolis-centripetal matrix C = C(nu) 

    Inputs:
    M        - 6x6 or 3x3 rigid-body MRB or added mass MA system marix 
    nu       - nu = [u, v, w, p, q, r]' or nu = [u, v, r]'

    The Coriolis and centripetal matrix depends on nu1 = [u,v,w]' and nu2 =
    [p,q,r]' as shown in Fossen (2021, Theorem 3.2). It is possible to
    compute C = C(nu2) where nu2 = [p,q,r]' using the linear velocity-
    independent representation, see 
    """
 
    M = 0.5 * (M + M.T)      #symmetrization of the inertia matrix
    M11 = M[0:3,0:3]
    M12 = M[0:3,3:6]
    M21 = M12.T
    M22 = M[3:6,3:6]
    
    nu1 = nu[0:3]
    nu2 = nu[3:6]
    dt_dnu1 = M11@nu1 + M12@nu2
    dt_dnu2 = M21@nu1 + M22@nu2

    C = np.block([[np.zeros((3,3))   , -Smtrx(dt_dnu1)],
                  [-Smtrx(dt_dnu1), -Smtrx(dt_dnu2) ]])
    
    return C

def Dmtrx(T_126,zeta_45,MRB,MA,hydrostatics):
    """ 
    D = Dmtrx([T1, T2, T6],[zeta4,zeta5],MRB,MA,hydrostatics)
     computes the 6x6 linear damping matrix for marine craft (submerged and
     floating) by specifying the time constants [T1, T2, T6] in DOFs 1,2 and 6. 
     The time constants can be found by open-loop step responses. For roll and
     pitch the relative damping ratios are specified using [zeta4, zeta5]. 
     For floating vessels it is assumed that zeta3 = 0.2 in heave, while
     submerged vehicles are assumed to be neutrally buoyant, W = B, with equal 
     time constants in heave and sway, that is T3 = T2.
    
     Inputs: T_126 = [T1, T2, T6]: time constants for DOFs 1, 2 and 6
             zeta_45 = [zeta4, zeta5]: relative damping ratios in DOFs 4 and 5
             MRB: 6x6 rigid-body system matrix (see rbody.m)
             MA:  6x6 hydrodynamic added mass system matrix
             hydrostatics = G for surface craft (see Gmtrx.m)
             hydrostatics = [W r_bg' r_bb'] for neutrally buoyant submerged 
               vehicles where W = m*g, r_bg = [xg,yg,zg]' and r_bb = [xb,yb,zb]'
    
     Output: D: 6x6 diagonal linear damping matrix
    """
    
    M = MRB + MA
    T1 = T_126[0]
    T2 = T_126[1]
    T6 = T_126[2]
    zeta4 = zeta_45[0]
    zeta5 = zeta_45[1]
    
    W = hydrostatics[0]       
    r_bg = hydrostatics[1]
    r_bb = hydrostatics[2]
    
    T3 = T2       # for AUVs, assume same time constant in sway and heave
    w4 = np.sqrt( W * (r_bg[2,0]-r_bb[2,0]) / M[3,3] )
    
    w5 = np.sqrt( W * (r_bg[2,0]-r_bb[2,0]) / M[4,4] )
    D  = np.diag(np.array([M[0,0]/T1, M[1,1]/T2, M[2,2]/T3, M[3,3]*2*zeta4*w4,  M[4,4]*2*zeta5*w5, M[5,5]/T6])) 
    
    return D

def C_bf(alpha,beta):
    # 
    """
    reiknar transformation matrix fra FLOW til BODY

    Args:
        alpha (float): angle of attack
        beta (float): sideslip angle

    Returns:
        3x3 matrix: snuningsfylki 
    """
    C_bf = np.array([np.cos(beta)*np.cos(alpha),  -np.sin(beta)*np.cos(alpha),      -np.sin(alpha), 
                     np.sin(beta),                 np.cos(beta),                    0,
                     np.cos(beta)*np.sin(alpha),  -np.sin(beta)*np.sin(alpha),      np.cos(alpha)]).reshape(3,3)
    
    return C_bf

def Rzyx(phi,theta,psi):
    # R = Rzyx(phi,theta,psi) computes the Euler angle
    # rotation matrix R in SO[2] using the zyx convention

    cphi = np.cos(phi)
    sphi = np.sin(phi)
    cth  = np.cos(theta)
    sth  = np.sin(theta)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    
    
    R = np.array([cpsi*cth,  -spsi*cphi+cpsi*sth*sphi,  spsi*sphi+cpsi*cphi*sth,
                spsi*cth,  cpsi*cphi+sphi*sth*spsi,   -cpsi*sphi+sth*spsi*cphi,
                -sth,      cth*sphi,                  cth*cphi]).reshape(3,3)

    return R


def Tzyx(phi,theta):
    # T = Tzyx(phi,theta) computes the Euler angle
    # transformation matrix T for attitude using the zyx convention

    cphi = np.cos(phi)
    sphi = np.sin(phi)
    cth  = np.cos(theta)
    sth  = np.sin(theta)
    
    if cth==0:
        print('Tzyx is singular for theta = +-90 degrees')

    T = np.array([1,  sphi*sth/cth,  cphi*sth/cth,
                0,  cphi,          -sphi,
                0,  sphi/cth,      cphi/cth]).reshape(3,3)
    return T

def eulerang(phi,theta,psi):
    """
    computes the Euler angle transformation matrix

    Args:
        phi (float): roll
        theta (float): pitch
        psi (float): yaw

    Returns:
        3x3 fylki: snuningsfylki
    """
    J1 = Rzyx(phi,theta,psi)
    J2 = Tzyx(phi,theta)

    J = np.block([[J1   , np.zeros((3,3))],
                  [np.zeros((3,3)), J2 ]])

    return [J, J1, J2]

def gRvect(W,B,R,r_bg,r_bb):
    """ 
     Computes the 6x1 vector of restoring 
     forces about an arbitrarily point CO for a submerged body using the 
     rotation matrix R as input. Use g = gvect(W,B,theta,phi,r_bg,r_bb) for 
     Euler angle inputs. For floating vessels, use Gmtrx.m. The following 
     examples show how to use gvect.m and gRvect.m.
    
     Euler angles: 
      R = Rzyx(phi,theta,psi)
      g = gRvect(W,B,R,r_bg,r_bb)            input: rotation matrix R
      g = gvect(W,B,theta,phi,r_bg,r_bb)     input: Euler angles phi, theta
     Unit quaternions:
      Rq = Rquat(q)
      g = gRvect(W,B,Rq,r_bg,r_bb)           input: rotation matrix Rq
    
     Inputs: 
      W, B: weight and buoyancy (kg)
      R: rotation matrix Rzyx (Euler angles) or Rquat (unit quaternions)
      r_bg = [x_g y_g z_g]: location of the CG with respect to the CO (m)
      r_bb = [x_b y_b z_b]: location of the CB with respect to th CO (m)
    
    """

    g=np.array([-(W-B) * R[2,0],
                -(W-B) * R[2,1],
                -(W-B) * R[2,2],
                -(r_bg[1,0]*W - r_bb[1,0]*B) * R[2,2] + (r_bg[2,0]*W - r_bb[2,0]*B) * R[2,1],
                -(r_bg[2,0]*W - r_bb[2,0]*B) * R[2,0] + (r_bg[0,0]*W - r_bb[0,0]*B) * R[2,2],
                -(r_bg[0,0]*W - r_bb[0,0]*B) * R[2,1] + (r_bg[1,0]*W - r_bb[1,0]*B) * R[2,0]]).reshape(6,1)
    return g

def gamma_sigma(DeltaX):
    """
    reiknar hornin fyrir R_I_tog

    Args:
        DeltaX (3x1 vigur): munur a stadsetningu togvirs a ramma og hanastels

    Returns:
        float, float: horn sem lysa stefnu togvirs sja skyrslu
    """
    gamma = np.arctan(DeltaX[2, 0]/DeltaX[0, 0])
    sigma = np.arcsin(DeltaX[1, 0]/np.linalg.norm(DeltaX))
    return gamma, sigma

def R_I_tog(g,s):
    """
    coordinate transformation matrix fra togvir til ramma

    Args:
        g (float): horn sem lysir stefnu togvirs sja skyrslu
        s (float): horn sem lysir stefnu togvirs sja skyrslu

    Returns:
        3x3 fylki: snuningsfylki
    """
    cg = np.cos(g)
    cs = np.cos(s)

    sg = np.sin(g)
    ss = np.sin(s)

    R_IT = np.array([
                    cg*cs,  -cg*ss,  -sg,
                    ss,       cs,     0,
                    cs*sg,  -sg*ss,  cg]).reshape(3,3)

    return R_IT


def togkraftur(DeltaX, DeltaU, stifni, dempun, upphafleg_lengd):
    """
    reiknar togkraft i togvir i hnitakerfi togvirs
    
    Args:
        DeltaX (3x1 vigur): munur a stadsetningu togvirs a ramma og hanastels
        DeltaU (3x1 vigur): munur a hrada togvirs a ramma og hanastels
        stifni (float): stifni togvirs
        dempun (float): dempun togvirs
        upphafleg_lengd (float): upphafleg lengd togvirs

    Returns:
        3x1 vigur: togkraftu i hnitakerfi togvirs (stefni i "x"-att i thvi hnitakerfi)
    """
    nuverandi_lengd = np.linalg.norm(DeltaX)
    if nuverandi_lengd > upphafleg_lengd:
        Tow = np.array([stifni*(nuverandi_lengd-upphafleg_lengd)+np.squeeze(dempun*(
        np.dot(DeltaU.T, (DeltaX/nuverandi_lengd)))), 0, 0]).reshape(3, 1)
    else:
        Tow = np.array([0, 0, 0]).reshape(3, 1)
    return Tow

