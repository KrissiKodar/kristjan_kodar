import numpy as np

def Smtrx(a):
    a = a.flatten()
    S = np.array([   0,    -a[2],    a[1],
                  a[2],        0,   -a[0],
                 -a[1],     a[0],       0 ]).reshape(3,3)
    
    return S
 
def m2c(M,nu):
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
    # transformation fra FLOW til BODY
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
    # [J,Rzyx,Tzyx] = eulerang(phi,theta,psi) computes the Euler angle
    # transformation matrix
    #
    #  J = [ Rzyx     0
    #           0  Tzyx ]
    #
    J1 = Rzyx(phi,theta,psi)
    J2 = Tzyx(phi,theta)

    J = np.block([[J1   , np.zeros((3,3))],
                  [np.zeros((3,3)), J2 ]])

    return [J, J1, J2]

def gRvect(W,B,R,r_bg,r_bb):
    # g = gRvect(W,B,R,r_bg,r_bb) computes the 6x1 vector of restoring 
    # forces about an arbitrarily point CO for a submerged body using the 
    # rotation matrix R as input. Use g = gvect(W,B,theta,phi,r_bg,r_bb) for 
    # Euler angle inputs. For floating vessels, use Gmtrx.m. The following 
    # examples show how to use gvect.m and gRvect.m.
    #
    # Euler angles: 
    #  R = Rzyx(phi,theta,psi)
    #  g = gRvect(W,B,R,r_bg,r_bb)            input: rotation matrix R
    #  g = gvect(W,B,theta,phi,r_bg,r_bb)     input: Euler angles phi, theta
    # Unit quaternions:
    #  Rq = Rquat(q)
    #  g = gRvect(W,B,Rq,r_bg,r_bb)           input: rotation matrix Rq
    #
    # Inputs: 
    #  W, B: weight and buoyancy (kg)
    #  R: rotation matrix Rzyx (Euler angles) or Rquat (unit quaternions)
    #  r_bg = [x_g y_g z_g]: location of the CG with respect to the CO (m)
    #  r_bb = [x_b y_b z_b]: location of the CB with respect to th CO (m)

    g=np.array([-(W-B) * R[2,0],
                -(W-B) * R[2,1],
                -(W-B) * R[2,2],
                -(r_bg[1,0]*W - r_bb[1,0]*B) * R[2,2] + (r_bg[2,0]*W - r_bb[2,0]*B) * R[2,1],
                -(r_bg[2,0]*W - r_bb[2,0]*B) * R[2,0] + (r_bg[0,0]*W - r_bb[0,0]*B) * R[2,2],
                -(r_bg[0,0]*W - r_bb[0,0]*B) * R[2,1] + (r_bg[1,0]*W - r_bb[1,0]*B) * R[2,0]]).reshape(6,1)
    return g

def R_I_tog(g,s):

    cg = np.cos(g)
    cs = np.cos(s)

    sg = np.sin(g)
    ss = np.sin(s)

    R_IT = np.array([
                    cg*cs,  -cg*ss,  -sg,
                    ss,       cs,     0,
                    cs*sg,  -sg*ss,  cg]).reshape(3,3)

    return R_IT

# reikna togkraft i togvir
def togkraftur(DeltaX, DeltaU, stifni, dempun, upphafleg_lengd):
    nuverandi_lengd = np.linalg.norm(DeltaX)
    if nuverandi_lengd > upphafleg_lengd:
        Tow = np.array([stifni*(nuverandi_lengd-upphafleg_lengd)+np.squeeze(dempun*(
        np.dot(DeltaU.T, (DeltaX/nuverandi_lengd)))), 0, 0]).reshape(3, 1)
    else:
        Tow = np.array([0, 0, 0]).reshape(3, 1)
    return Tow

def gamma_sigma(DeltaX):
    gamma = np.arctan(DeltaX[2, 0]/DeltaX[0, 0])
    sigma = np.arcsin(DeltaX[1, 0]/np.linalg.norm(DeltaX))
    return gamma, sigma
