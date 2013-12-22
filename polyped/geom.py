"""
Useful geometry functions

(c) Sam Burden, Berkeley and Shai Revzen, UPenn 2010
"""

import numpy as np
import numpy.linalg as la
import scipy.optimize as opt

def euler(a):
    """
    R = geom.euler  generate rotation matrix from euler angles

    INPUTS
      a - 1 x 3 - euler angles, i.e. rotations around x,y,z axes

    OUTPUTS
      R - 3 x 3 - rotation matrix as Rx * Ry * Rz
    """
    a = np.asarray(a).flatten()

    # Trig functions of rotations
    c = np.cos(a[0:3])
    s = np.sin(a[0:3])

    # Rotation matrix as Rx * Ry * Rz
    R = np.asarray([[c[2]*c[1], s[2]*c[0]+c[2]*s[1]*s[0], s[2]*s[0]-c[2]*s[1]*c[0]],
                    [-s[2]*c[1], c[2]*c[0]-s[2]*s[1]*s[0], c[2]*s[0]+s[2]*s[1]*c[0]],
                    [s[1], -c[1]*s[0], c[1]*c[0]]])

    return R

def fit(p0, ij, d):
    """
    [p, err, res] = geom.fit  fits 3D geometry to pairwise distances

    INPUTS
      p0 - M x D - initial guess for M points in D dimensional space
      ij - N x 2 - indices of points whose distance is in d
      d  - N x 1 - distance vector

    OUTPUTS
      p - M x D - solution for M points in D dimensional space
      err - scalar - norm of residual error
      res - N x 1 - residual for each distance measurement
    """

    def cost(x, dat):
        """
        err = geom.fit.cost  computes error in geometry fit

        INPUTS
          x - M x D - current guess for M points in D dimensional space
          dat - dict - problem data

        OUTPUTS
          err - N x 1 - vector of errors
        """
        dp = np.reshape(x,dat['sz'])

        dij = dp[dat['ij'][:,0],:] - dp[dat['ij'][:,1],:] + dat['dp0']

        ofs = np.sqrt((dij**2).sum(1)) - dat['d']

        mx = dp.mean(0)

        rot = dij.flatten()[dat['ind']]

        err = np.hstack((ofs, mx, rot))

        return err
         

    dp0 = p0[ij[:,0],:] - p0[ij[:,1],:]

    ind = np.triu(np.ones(dp0.shape)).flatten().nonzero()

    sz = p0.shape
    dat = {'p0'  : p0, 
           'dp0' : dp0,
           'ind' : ind,
           'ij'  : ij,
           'd'   : d.flatten(),
           'sz'  : sz}

    p, Cp, info, msg, flag = opt.leastsq(cost, p0.flatten(), args=dat, full_output=1)
    p = np.reshape(p, sz)

    return p, (Cp, info, msg), flag

def plane(p):
  """
  n = plane  fit plane to data that is not parallel to last axis

  INPUTS
    p - N x D - N samples of D-dimensional points on (i.e. near) the plane

  OUTPUTS
    n - 1 x D - normal vector to plane
  """
  q = p - p.mean(axis=0)
  n = np.hstack((-la.lstsq(q[:,:-1],q[:,-1])[0],1.))
  return n

def orient2(x,y):
    """
    R = orient  rotation matrix bringing two orthogonal vectors in 
                line with [1,0,0] and [0,1,0]

    INPUTS
      x, y - 3 x 1 - orthogonal vectors to align with [1,0,0] and [0,1,0]

    OUTPUTS
      R - 3 x 3 - rotation matrix which orients coordinate system
    """
    x = x.reshape((3,1))
    y = y.reshape((3,1))
    if (x*y).sum() > 1e-6:
        raise ValueError,'The provided vectors are not orthogonal'
    R1 = euler(np.array([0,0,np.arctan2(x[1,0],x[0,0])]))
    x = np.dot(R1,x)
    y = np.dot(R1,y)
    R2 = euler(np.array([0,-np.arctan2(x[2,0],x[0,0]),0]))
    y = np.dot(R2,y)
    R3 = euler(np.array([np.arctan2(y[2,0],y[1,0]),0,0]))

    return np.dot(R3,np.dot(R2,R1))

def orient(z):
    """
    R = orient  rotation matrix bringing vector in line with [0,0,1]

    INPUTS
      z - 3 x 1 - vector to align with [0,0,1]

    OUTPUTS
      R - 3 x 3 - rogation matrix which orients coordinate system
    """
    x0 = z.reshape((3,1))
    R1 = euler(np.array([0,np.arctan2(x0[0,0],x0[2,0]),0]))
    x1 = np.dot(R1,x0)
    R2 = euler(np.array([-np.arctan2(x1[1,0],x1[2,0]),0,0]))
    x2 = np.dot(R2,x1)

    return np.dot(R2,R1)

