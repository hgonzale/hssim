"""
Hybrid Dynamical System numerical integration package

Primary contents:
  HDS - class - hybrid dynamical system
  Euler - func -  numerical integration

Helpful routines:
  stack - func - aggregate data from trajectory snippets returned from Euler
  trans - func - aggregate data from trajectory snippets returned from Euler
  fp - func - find fixed point of a map
  jac - func - numerically approximate Jacobian of a map
  armijo - func - Armijo stepsize rule
  descent1 - func - gradient descent
  proj - func - project given point onto constraint manifold
  tgt - func - approximate tangent space to manifold
  armijoproj - func - Armijo stepsize rule, projected onto manifold
  descent1proj - func - gradient descent, projected onto manifold

by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011
"""

import numpy as np
import pylab as plt
import scipy as sp
import scipy.optimize as op

try:
  from IPython.Debugger import Tracer; dbg = Tracer()
except:
  dbg = lambda : 1/0
import sys

class Struct(object):
    """ 
    Struct  struct-like class

    Acts like a dict, but keys are members of the object.

    Also, its string representation is Python code which 
    yields the object when executed.

    >>> a = Struct(foo='bar', goo=7)
    >>> b = a.copy()
    >>> b.hoo = [1,2,3]
    >>> print a
    Struct(**{'goo': 7, 'foo': 'bar'})
    >>> print b
    Struct(**{'hoo': [1, 2, 3], 'goo': 7, 'foo': 'bar'})
    >>> c = eval(str(a))
    >>> a.ioo = []
    >>> print a
    Struct(**{'ioo': [], 'goo': 7, 'foo': 'bar'})
    >>> print c
    Struct(**{'goo': 7, 'foo': 'bar'})
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
    def copy(self):
        return Struct(**self.__dict__.copy())
    def __repr__(self):
        return 'Struct(**'+str(self.__dict__)+')'
    def __getstate__(self):
        return self.__dict__
    def __setstate__(self, kwds):
        self.__dict__.update(kwds)

class HDS(object):

  def __init__(self):
    """
    hds = HDS(debug)

    Hybrid dynamical system superclass

    Full hybrid systems must override:
      .F, .R, .E, .P, .G, .O;
    though trivial implementations are provided.

    by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011
    """

  def P(self, debug=False):
    """
    Parameters
      debug - bool - flag for printing debugging info
    """
    return Struct(j=None,debug=False)

  def G(self,t,x,p):
    """
    g = G(x,p)  guard function

    g > 0 : guard inactive
    g = 0 : guard 
    """
    return 0

  def F(self,t,z,p):
    """
    dz = F(t,z,p)  vector field
    """
    return 0.*z

  def R(self,t,z,p):
    """
    z,p = R(t,z,p)  reset map
    """
    return t,z,p

  def E(self,t,z,p,v):
    """
    z = E(t,z,p,v)  retraction
    """
    return z

  def O(self,t,x,p):
    """
    o = O(t,x,p)  observation
    """
    return x

  def dist(self, x,y):
    """
    .dist  distance in state space
    """
    return np.max(np.abs(x-y))

class BB(HDS):

  def __init__(self):
    """
    bb = BB()

    Bouncing ball hybrid dynamical model

    Outputs:
    bb - struct - bouncing ball simulation object
      .F - vector field   dx  = F(t,x,p) dt
      .R - reset map      t,x,p = R(t,x,p)
      .E - retraction     x   = E(t,x,p,v)
      .P - parameters     p   = P(g,c,s,debug)
        .j - discrete mode index
      .G - guard function g = G(t,x,p)
      .O - observations   o = O(t,x,p)

    by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011
    """
    HDS.__init__(self)

  def P(self, g=1., c=1., debug=False):
    """
    Parameters:
      g - scalar - gravitational constant
      c - scalar - coefficient of restitution
      debug - bool - flag for printing debugging info
    """
    return Struct(j=1,g=g,c=c,debug=debug)

  def G(self,t,x,p):
    """
    guard function
    g > 0 : guard inactive
    g = 0 : guard 
    """
    if p.j == 1:
      g = x[0]

    else:
      raise RuntimeError,'(G) unknown discrete mode'

    return g

  def F(self,t,z,p):
    """
    dz = F(t,z,p)
    """

    if p.j == 1:
      h,dh = z
      
      dz = np.array([dh, -p.g])

    return dz  

  def R(self,t,z,p):
    """
    z,p = R(t,z,p)
    """
    if p.j == 1:
      h,dh = z

      z = np.array([0, -p.c*dh])
      p = p

    return t,z,p

  def E(self,t,z,p,v):
    """
    z = E(t,z,p,v)
    """
    z = z + v

    return z

  def O(self,t,x,p):
    """
    o = O(x,p)
    """
    if p.j == 1:
      h,dh = x.T

      h.shape  = (h.size,1)
      dh.shape = h.shape

      z = np.zeros(h.shape)

      m = 1
      KE = 0.5*m*dh**2
      PE = m*p.g*h
      F = [];
      for tt,xx in zip(t,x):
        F += [self.F(tt,xx,p)]
      Fn = np.sqrt(np.sum(np.array(F)**2,axis=1).reshape(-1,1))

      o = np.hstack((h,dh,KE,PE,Fn))

    return o

def Euler(t, x, p, hds, h, rx, n, debug=False, Zeno=False):
  """
  trjs = Euler(t, x, p, hds, h, rx, n, debug, Zeno)

  Compute Euler approximation to relaxed hybrid execution

  Inputs:
  t   - 1 x 2 - time range
  x   - 1 x n - initial state
  p   - Struct - model paramters
  hds - hybrid model
  h   - scalar - Euler integration stepsize
  rx  - scalar - relaxation parameter
  n   - int - maximum number of hybrid transitions
  debug - flag for printing debugging info
  Zeno - flag for quitting executions with short modes

  Outputs:
  trjs - list of trajectory dicts
  trj - trajectory dict
    .t - times
    .x - states
    .p - parameters
    .hds - hybrid model

  by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011
  """
  tf = t[1] # store final time

  h0 = h # store given stepsize 

  p.debug = debug # enable model debugging

  t = np.array([t[0]])
  x = x.copy()
  if len(x.shape) == 1:
    x.shape = (1,x.size)

  trj = Struct(t=t,x=x,p=p,hds=hds)
  trjs = []

  # Euler loop
  while trj.t[-1] <= tf and not (trj.p.j == None) and len(trjs) < n:

    # do single Euler step
    t0 = trj.t[-1]
    x0 = trj.x[-1,:]
    p0 = trj.p
    dxdt = trj.hds.F(t0, x0, p0)
    dx = h*dxdt

    j = trj.p.j
    t = t0 + h
    x = trj.hds.E(t0, x0, p0, dx)
    g = trj.hds.G(t,  x,  p0)

    # halve step size until trajectory doesn't jump over strip
    k = 0
    kmax = 50
    while np.any(g < -rx) and (k <= kmax):
      #if debug:
      #  print 'RX:  jump over strip #%d' % k
      h  = h/2.
      t  = t0 + h
      dx = h*dxdt
      x  = trj.hds.E(t0, x0, p0, dx)
      g  = trj.hds.G(t,  x,p0)
      k += 1

    if (k >= kmax):
      raise RuntimeError,'(euler)  strip iterations exceeded'

    # append state to trj
    trj.t = np.append(trj.t, [t], axis=1)
    trj.x = np.append(trj.x, [x], axis=0)

    if debug:
      print '  :  j = %s, h = %0.2e, t = %0.3f, g = %s, x = %s' % (j,h,t,g,str(x))

    # if state is on strip
    if np.any(g < 0):

      # spend time on the strip 
      t = t + (rx + g.min())
      trj.t = np.append(trj.t, [t], axis=1)
      # can't move state without analytical strip
      trj.x = np.append(trj.x, [x], axis=0)

      if debug:
        print 'RX:  j = %s, t = %0.3f, g = %s, x = %s' % (j,t,g,str(x))

      # append trj to trjs
      trjs += [trj]
      trj = trj.copy()

      if Zeno and (len(trj.t) <= 4):

        print '(euler)  possible Zeno @ stepsize h = %0.6f' % h0
        print 'RX:  j = %s  t = %0.3f\nx = %s' % (j,t,str(x))
        return trjs

      # apply reset to modify trj
      t,x,p = hds.R(t,x,p0)
      if len(x.shape) == 1:
        x.shape = (1,x.size)

      t = np.array([t])
      trj.t = t
      trj.x = x
      trj.p = p

      # re-initialize step size
      h = h0

      if debug:
        j = p.j
        g = trj.hds.G(t[0],  x[0],  p)
        print 'RX:  j = %s, t = %0.3f, g = %s, x = %s' % (j,t,g,str(x[0,:]))

  trjs += [trj]
  trj = trj.copy()

  return trjs

def EulerBisect(t, x, p, hds, h, rx, n, debug=False, Zeno=False):
  """
  trjs = EulerBisect(t, x, p, hds, h, rx, n, debug, Zeno)

  Compute Euler approximation to relaxed hybrid execution

  Inputs:
  t   - 1 x 2 - time range
  x   - 1 x n - initial state
  p   - Struct - model paramters
  hds - hybrid model
  h   - scalar - Euler integration stepsize
  rx  - scalar - relaxation parameter
  n   - int - maximum number of hybrid transitions
  debug - flag for printing debugging info
  Zeno - flag for quitting executions with short modes

  Outputs:
  trjs - list of trajectory dicts
  trj - trajectory dict
    .t - times
    .x - states
    .p - parameters
    .hds - hybrid model

  by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011
  """
  tf = t[1] # store final time

  h0 = h # store given stepsize 

  p.debug = debug # enable model debugging

  t = np.array([t[0]])
  x = x.copy()
  if len(x.shape) == 1:
    x.shape = (1,x.size)

  trj = Struct(t=t,x=x,p=p,hds=hds)
  trjs = []

  # Euler loop
  while trj.t[-1] <= tf and not (trj.p.j == None) and len(trjs) < n:

    # initial data
    t0 = trj.t[-1]
    x0 = trj.x[-1,:]
    p0 = trj.p
    dx = trj.hds.F(t0, x0, p0)
    # Euler step
    j = trj.p.j
    t = t0 + h
    x = trj.hds.E(t0, x0, p0, h*dx)
    g = trj.hds.G(t,  x,  p0)

    # if transition occurs
    if ( g < -rx ):
      # find step size that lands trajectory on strip using bisection
      h = bisect( lambda s : trj.hds.G(t0+s,trj.hds.E(t0,x0,p0,s*dx),p0),
                  (0.,h), tol=rx )
      # debug
      if np.isnan(h):
        raise RuntimeError,'(euler)  cannot step to strip'
      # Euler step
      t = t0 + h
      x = trj.hds.E(t0, x0, p0, h*dx)
      g = trj.hds.G(t,  x,  p0)

    # append state to trj
    trj.t = np.append(trj.t, [t], axis=1)
    trj.x = np.append(trj.x, [x], axis=0)

    if debug:
      print '  :  j = %d  g = %0.2e  h = %0.2e  t = %0.3f  x = %s' % (j,g,h,t,str(x))

    # if state is on strip
    if (g < 0):

      # spend time on the strip 
      t = t + (rx + g)
      trj.t = np.append(trj.t, [t], axis=1)
      # can't move state without analytical strip
      trj.x = np.append(trj.x, [x], axis=0)

      if debug:
        print 'RX:  j = %d  g = %0.2e  h = %0.2e  t = %0.3f  x = %s' % (j,g,h,t,str(x))

      # append trj to trjs
      trjs += [trj]
      trj = trj.copy()

      if Zeno and (len(trj.t) <= 4):

        print '(euler)  possible Zeno @ stepsize h = %0.6f' % h0
        print 'RX:  j = %d  g = %0.2e  h = %0.2e  t = %0.3f  x = %s' % (j,g,h,t,str(x))
        return trjs

      # apply reset to modify trj
      t,x,p = hds.R(t,x,p0)
      if len(x.shape) == 1:
        x.shape = (1,x.size)

      t = np.array([t])
      trj.t = t
      trj.x = x
      trj.p = p

      # re-initialize step size
      h = h0

      if debug:
        print 'RX:  j = %d  g = %0.2e  h = %0.2e  t = %0.3f  x = %s' % (j,g,h,t,str(x[0,:]))

  trjs += [trj]
  trj = trj.copy()

  return trjs

def EulerPlanar(t, x, p, hds, h, n, debug=False, Zeno=False):
  """
  trjs = Euler(t, x, p, hds, h, n, debug, Zeno)

  *** UNFINISHED ***

  Compute Euler approximation to hybrid execution assuming 
  guards are coordinate functions
  (thus transition time can be determined exactly)

  Inputs:
  t   - 1 x 2 - time range
  x   - 1 x n - initial state
  p   - Struct - model paramters
  hds - hybrid model
  h   - scalar - Euler integration stepsize
  n   - int - maximum number of hybrid transitions
  debug - flag for printing debugging info
  Zeno - flag for quitting executions with short modes

  Outputs:
  trjs - list of trajectory dicts
  trj - trajectory dict
    .t - times
    .x - states
    .p - parameters
    .hds - hybrid model

  by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011
  """
  tf = t[1] # store final time

  h0 = h # store given stepsize 

  p.debug = debug # enable model debugging

  t = np.array([t[0]])
  x = x.copy()
  if len(x.shape) == 1:
    x.shape = (1,x.size)

  trj = Struct(t=t,x=x,p=p,hds=hds)
  trjs = []

  # Euler loop
  while trj.t[-1] <= tf and not (trj.p.j == None) and len(trjs) < n:

    # initial data
    t0 = trj.t[-1]
    x0 = trj.x[-1,:]
    p0 = trj.p
    g0 = trj.hds.G(t, x,  p0)
    # vector field
    dx = trj.hds.F(t0, x0, p0)
    # Euler step
    j = trj.p.j
    t = t0 + h
    x = trj.hds.E(t0, x0, p0, h*dx)
    g = trj.hds.G(t,  x,  p0)
    # if guard is triggered
    if g < 0:
      # move state exactly to transition
      s = g0 / trjs.hds.G(t, h*dx, p0)


    # halve step size until trajectory doesn't jump over strip
    k = 0
    kmax = 50
    while (g < -rx) and (k <= kmax):
      if debug:
        print 'RX:  jump over strip #%d' % k
      h  = h/2
      t  = t0 + h
      dx = h*dxdt + np.sqrt(h)*np.random.randn()*dxdw
      x  = trj.hds.E(t, x0, p0, dx)
      g  = trj.hds.G(t, x,p0)
      k += 1

    if (k >= kmax):
      raise RuntimeError,'(euler)  strip iterations exceeded'

    # append state to trj
    trj.t = np.append(trj.t, [t], axis=1)
    trj.x = np.append(trj.x, [x], axis=0)

    if debug:
      print '  :  j = %d  t = %0.3f  x = %s' % (j,t,str(x))

    # if state is on strip
    if (g < 0):

      # spend time on the strip 
      t = t + (rx + g)
      trj.t = np.append(trj.t, [t], axis=1)
      # can't move state without analytical strip
      trj.x = np.append(trj.x, [x], axis=0)

      if debug:
        print 'RX:  j = %d  t = %0.3f  x = %s' % (j,t,str(x))

      # append trj to trjs
      trjs += [trj]
      trj = trj.copy()

      if Zeno and (len(trj.t) <= 4):

        print '(euler)  possible Zeno @ stepsize h = %0.6f' % h0
        print 'RX:  j = %d  t = %0.3f  x = %s' % (j,t,str(x))
        return trjs

      # apply reset to modify trj
      t,x,p = hds.R(t,x,p0)
      if len(x.shape) == 1:
        x.shape = (1,x.size)

      t = np.array([t])
      trj.t = t
      trj.x = x
      trj.p = p

      # re-initialize step size
      h = h0

      if debug:
        print 'RX:  j = %d  t = %0.3f  x = %s' % (j,t,str(x[0,:]))

  trjs += [trj]
  trj = trj.copy()

  return trjs

def stack(trjs, J=None):
  """
  t,x,o,p = stack(trjs,J)

  Contatenate times & observations from trajectories

  Inputs:
    trjs - struct array - trajectories
    J - array - discrete modes to include (OPTIONAL)

  Outputs:
    t - list of 1 x Nt - observation times
    x - list of 1 x Ns - states
    o - list of 1 x No - observations
    p - list of Struct - parameters

  by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011
  """
  T = []
  X = []
  O = []
  P = []
  if J is not None:
    J = np.array(J)

  for k in range(len(trjs)):
    if J is None or np.any(trjs[k].p.j == J):
      t = trjs[k].t
      x = trjs[k].x
      p = trjs[k].p
      o = trjs[k].hds.O(t,x,p)

      T += [t]; X += [x]; O += [o]; P += [p]
 
  return T,X,O,P

def trans(trjs, J=None):
  """
  t,x,o,p = trans(trjs,J)

  Contatenate times & observations from trajectories

  Inputs:
    trjs - struct array - trajectories
    J - array - discrete modes to include (OPTIONAL)

  Outputs:
    t - list of 1 x Nt - observation times
    x - list of 1 x Ns - states
    o - list of 1 x No - observations
    p - list of Struct - parameters

  by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011
  """
  T = []
  O = []
  X = []
  P = []
  if J is not None:
    J = np.array(J)

  for k in range(len(trjs)):
    if ( J is None or np.any(trjs[k].p.j == J) ) and len(trjs[k].t) >= 3:
      t = trjs[k].t.take([0,-2,-1])
      x = trjs[k].x.take([0,-2,-1],axis=0)
      p = trjs[k].p
      o = trjs[k].hds.O(t,x,p)

      T += [t]; X += [x]; O += [o]; P += [p]
     
  return T,X,O,P

def fp(f,x0,eps=1e-6,modes=[1,2],Ni=4,N=10,
       dist=lambda x,y : np.max(np.abs(x-y)), debug=False ):
  """
  .fp  Find fixed point of f near x0 to within tolerance eps

  The algorithm has three modes:
    1. iterate f; keep result if error decreases initially
    2. run fsolve on f(x)-x; keep result if non-nan
    3. run descent1 on |f(x) - x|^2; keep result if non-nan

  Inputs:
    f : R^n --> R^n
    x0 - n vector - initial guess for fixed point  f(x0) ~= x0
    modes - list - fixed point finding algorithms to run

  Outputs:
    x - n vector - fixed point  f(x) ~= x
  """
  # compute initial error
  xx = f(x0); e = dist(xx,x0)
  suc = False
  # 1. iterate f; keep result if error decreases initially
  if 1 in modes:
    # Iterate orbit map several times, compute error
    x = reduce(lambda x, y : f(x), [x0] + range(Ni))
    xx = f(x); e = dist(xx,x); e0 = dist(x,x0)
    # If converging to fixed point
    if e < e0:
      suc = True
      # Iterate orbit map
      n = 0
      while n < N-Ni and e > eps:
        n = n+1; x = xx; xx = f(x)
        e = dist(xx,x)
      x0 = xx
  # 2. run fsolve on f(x)-x; keep result if non-nan
  if 2 in modes:
    x = x0
    # Try to find fixed point using op.fsolve
    xx = op.fsolve(lambda x : f(x)-x, x)
    # If op.fsolve succeeded
    if not np.isnan(xx).any():
      suc = True
      x0 = xx
  # 3. run descent1 on |f(x) - x|^2; keep result if non-nan
  if 3 in modes:
    x = x0
    xx,_ = descent1(lambda x : np.sum((f(x) - x)**2),x0,eps=eps,debug=debug)
    # If descent1 succeeded
    if not np.isnan(xx).any():
      suc = True
      x0 = xx
  # if all methods failed, return nan
  if not suc:
    x0 = np.nan*x0
  
  return x0

def central(f,x,fx,d):
  """
  df = central()  compute central difference 

  df = 0.5*(f(x+d) - f(x-d))/np.linalg.norm(d)
  """
  return 0.5*(f(x+d) - f(x-d))/np.linalg.norm(d)

def forward(f,x,fx,d):
  """
  df = forward()  compute forward difference 

  df = (f(x+d) - fx)/np.linalg.norm(d)
  """
  return (f(x+d) - fx)/np.linalg.norm(d)

def jac(f, x, fx=None, d=1e-6, D=None, diff=forward):
  """
  .jac  Numerically approximate Jacobian to f at x

  Inputs:
    f : R^n --> R^m
    x - n vector
    d - scalar or (1 x n) - displacement in each coordinate
  (optional)
    fx - m vector - f(x)
    D - k x n - directions to differentiate (assumes D.T D invertible)
    diff - func - numerical differencing method

  Outputs:
    Df - m x n - Jacobian of f at x
  """
  if fx is None:
    fx = f(x)
  if D is None:
    J = map(lambda dd : diff(f,x,fx,dd), list(d*np.identity(len(x))))
  else:
    J = map(lambda dd : diff(f,x,fx,dd), list(D))
    dbg()

  return np.array(J).T

def bisect(f, ab, tol=1e-12, nmax=50):
  """
  .bisect  Find root of scalar function using bisection
           (searches from the right, i.e. from b to a)

  Inputs
    f : R --> R
    ab - (a,b) - interval to search
  (optional)
    tol - terminate when  (b-a)/2 < tol
    nmax - no more than  nmax  bisections

  Outputs
    c - root, i.e.  f(c) ~= 0
  """
  a,b = ab
  for n in range(nmax):
    c = (a + b) / 2.
    if ( f(c) == 0 ) or ( (b - a) / 2. < tol ):
      return c
    elif ( (b - a) / 2. < tol ):
      return b
    if np.sign(f(c)) == np.sign(f(b)):
      b = c
    else:
      a = c
  return np.nan

def armijo(f, x, d, fx=None, Dfx=None, 
           alpha=0.1, beta=0.5, gamma=1., kmax=100):
  """
  .armijo  Compute stepsize using Armijo rule

  i.e.   a = beta**l , where
  l = min{ k <= kmax : f(x) - f(x + beta**k d) >= alpha * beta**k Dfx.T * d }

  P. Wolfe.  Convergence conditions for ascent methods.
  SIAM Review, pgs 226-235, 1969.

  Inputs
    f : R^n --> R
    x - n vector
    d - n vector - descent direction
  (optional)
    Dfx - 1 x n - Jacobian of f at x
    alpha - Armijo constant # 1
    beta  - Armijo constant # 2
    kmax  - maximum power of beta

  Outputs
    a - scalar - stepsize
  """
  # initialize stepsize
  a = np.nan
  # evaluate objective function if necessary
  if fx is None:
    fx = f(x)
  # compute gradient of f if necessary
  if Dfx is None:
    Dfx = jac(f,x,fx=fx)
  # loop over Armijo exponents
  for k in range(1,kmax+1):
    betak = gamma * beta**k
    fxd = f(x + betak * d)
    # skip step if objective function undefined
    if np.isnan(fxd):
      continue
    # terminate if Armijo condition is satisfied
    if (fx - fxd) >= -alpha * betak * np.dot(Dfx,d):
      a = betak
      break
  # return Armijo stepsize
  return a

def descent1(f, x0, D=lambda f,x,fx : jac(f,x,fx=fx), 
             eps=1e-4, ftol=-np.inf, d=None, S=None, a=armijo, 
             norm=lambda x : np.sqrt(np.sum(x * x)), 
             jmax=100, xmax=1e6, debug=False):
  """
  .descent1  Gradient descent to find stationary points

  Inputs
    f : R^n --> R
    x0 - n vector - initial guess for minimum
  (optional)
    D  : (f,x)   |--> (1 x n) - gradient
    eps - scalar - threshold on norm of gradient for stationary point
    d  : (f,x)   |--> (n x 1) - descent direction
    a  : (f,x,d) |--> scalar  - stepsize
    jmax - integer - maximum number of iterations
    xmax - scalar  - largest allowable state
    debug - bool - debugging flag

  Outputs
    x - n vector - stationary point, i.e.  || D(f,x) || < eps
    j - int - number of iterations
  """
  # initial point
  x = x0
  # evaluate function
  fx = f(x)
  # compute gradient
  Dfx = D(f,x,fx)
  # default scaling is identity
  if S is None:
    S = np.identity(len(x))
  elif ( len(S.shape) == 1 ):
    S = np.diag(S)
  if debug:
    print '(descent1)'
    print '    x = %s' % (x)
    print ' f(x) = %0.2e' % (fx)
    print 'Df(x) = %s' % (Dfx)
    sys.stdout.flush()
  # loop over descent iterations
  for j in range(jmax):
    # fail if function or gradient are undefined
    if np.isnan(fx):
      if debug:
        print '(descent1)  f(x) is nan'
        dbg()
      return x,j
    if np.any(np.isnan(Dfx)):
      if debug:
        print '(descent1)  Df(x) is nan'
        dbg()
      return x,j
    # succeed at stationary point
    if norm(Dfx) < eps or fx < ftol:
      return x,j
    # descend gradient
    if d is None:
      dd = -np.dot(S,Dfx.T)
    # descend given direction
    else:
      dd = d(f,x)
    # compute stepsize
    aa = a(f,x,dd,fx=fx,Dfx=Dfx)
    # fail if stepsize cannot be chosen
    if np.isnan(aa):
      if debug:
        print '(descent1)  stepsize is nan'
        dbg()
      return x,j
    # descend
    x = x + aa * dd
    # fail if point is unbounded
    if norm(x) >= xmax:
      if debug:
        print '(descent1)  norm(x) >= xmax'
        dbg()
      return x,j
    # evaluate function
    fx = f(x)
    # compute gradient
    Dfx = D(f,x,fx)
    if debug:
      print 'j = %3d, a = %0.2e' % (j,aa)
      print '   x = %s' % x
      print 'f(x) = %0.2e' % (f(x))
      print '   d = %s' % (dd)
      sys.stdout.flush()
      #dbg()
  # descent failed to converge
  if debug:
    print '(descent1)  j >= jmax'
    dbg()
  return x,j

def proj(v,U):
  """
  u = proj  Project vector orthogonally onto column span

  Inputs
    v - n x 1 - vector to project
    U - m x n - columns span subspace 

  Outputs
    u := U (U^T U)^{-1} U^T v
  """
  U = np.matrix(U)
  n,m = U.shape
  v = np.matrix(v).reshape(n,1)
  u = np.array((U*((U.T*U).I)*U.T)*v).flatten()

  return u

def tgt(f,x,d=1e-6,Map=map):
  """
  U = tgt  Approximate tangent space

  Inputs
    f : R^n --> R^m - submanifold coordinates
    x \in R^n
  (optional)
    d - scalar - diameter of `tangent' neighborhood

  Outputs
    U - m x n - columns span T_f(x) f(R^n)
  """
  x = np.array(x).flatten(); n = x.size;
  # apply f to points near x
  fx = f(x); m = fx.size
  X = x + np.vstack((d*np.identity(n),
                     -d*np.identity(n)))
  Z = np.array(Map(f,X))
  # find orthonormal basis for tangent space T_f(x) f(R^m)
  U,_,_ = np.linalg.svd(np.cov(Z.T))
  U = U[:,:n]

  return U

def armijoproj(f, x, d, fx=None, Dfx=None, 
               alpha=0.1, beta=0.8, gamma=1., kmax=100,
               TM=lambda x : np.identity(x.size), Pi=lambda x : x):
  """
  .armijoproj  Compute stepsize using Armijo rule

  i.e.   a = beta**l , where
  l = min{ k <= kmax : f(x) - f(x + beta**k d) >= alpha * beta**k Dfx.T * d }

  P. Wolfe.  Convergence conditions for ascent methods.
  SIAM Review, pgs 226-235, 1969.

  Inputs
    f : R^n --> R
    x - n vector
    d - n vector - descent direction
  (optional)
    fx - scalar - f(x)
    Dfx - 1 x n - Jacobian of f at x
    alpha - Armijo constant # 1
    beta  - Armijo constant # 2
    kmax  - maximum power of beta
    Pi : R^n --> M - projection onto constraint manifold
    TM : R^n --> R^(n x m) - tangent space to M

  Outputs
    a - scalar - stepsize
  """
  # initialize stepsize
  a = np.nan
  TxM = TM(x)
  # evaluate objective function if necessary
  if fx is None:
    fx = f(x)
  # compute & project gradient of f onto tangent space if necessary
  if Dfx is None:
    Dfx = proj(D(f,x),TxM)
  # loop over Armijo exponents
  for k in range(1,kmax+1):
    betak = gamma * beta**k
    fxd = f( Pi(x + betak * d) )
    # skip step if objective function undefined
    if np.isnan(fxd):
      continue
    # terminate if Armijo condition is satisfied
    if (fx - fxd) >= -alpha * betak * np.dot(Dfx,d):
      a = betak
      break
  # return Armijo stepsize
  return a

def descent1proj(f, x0, D=lambda f,x : jac(f,x), eps=1e-4, d=None, 
                    a=armijoproj, 
                    norm=lambda x : np.sqrt(np.sum(x * x)), 
                    jmax=100, xmax=1e6, debug=False, 
                    TM=lambda x : np.identity(x.size), Pi=lambda x : x,
                    alpha=0.1, beta=0.8, gamma=1.):
  """
  .descent1proj  Gradient descent to find stationary points, 
                 reprojecting into constraint manifold M

  Inputs
    f : R^n --> R
    x0 - n vector - initial guess for minimum
  (optional)
    D  : (f,x)   |--> (1 x n) - gradient
    eps - scalar - threshold on norm of gradient for stationary point
    d  : (f,x)   |--> (n x 1) - descent direction
    a  : (f,x,d) |--> scalar  - stepsize
    Pi : R^n --> M - projection onto constraint manifold
    TM : R^n --> R^(n x m) - tangent space to M
    jmax - integer - maximum number of iterations
    xmax - scalar  - largest allowable state
    debug - bool - debugging flag

  Outputs
    x - n vector - stationary point, i.e.  || D(f,x) || < eps
  """
  # initial point
  x = Pi(x0)
  TxM = TM(x)
  # evaluate function
  fx = f(x)
  # compute & project gradient onto tangent space
  Dfx = proj(D(f,x),TxM)
  # loop over descent iterations
  for j in range(jmax):
    # fail if function or gradient are undefined
    if np.isnan(fx):
      print '(descent1proj)  f(x) is nan'
      return x,j
    if np.any(np.isnan(Dfx)):
      print '(descent1proj)  Df(x) is nan'
      return x,j
    # succeed at stationary point
    if norm(Dfx) < eps:
      return x,j
    # descend gradient
    if d is None:
      dd = -Dfx.T
    # descend given direction
    else:
      dd = d(f,x)
    # compute stepsize
    aa = a(f,x,dd,fx=fx,Dfx=Dfx,TM=lambda x : TxM,Pi=Pi,
           alpha=alpha,beta=beta,gamma=gamma)
    # fail if stepsize cannot be chosen
    if np.isnan(aa):
      print '(descent1proj)  stepsize is nan'
      return x,j
    # descend
    x = Pi(x + aa * dd)
    TxM = TM(x)
    # fail if point is unbounded
    if norm(x) >= xmax:
      print '(descent1proj)  norm(x) >= xmax'
      return x,j
    # evaluate function
    fx = f(x)
    # compute & project gradient onto tangent space
    Dfx = proj(D(f,x),TxM)
    if debug:
      print 'j = %3d, a = %0.2e, x = %s, f(x) = %0.2e' % (j,aa,x,f(x))
      print '                       d = %s' % (dd)
      sys.stdout.flush()
  # descent failed to converge
  print '(descent1proj)  j >= jmax'
  return x,j


if __name__ == "__main__":

  t = (0.,10.)
  #t = (0.,0.)
  x = np.array([1., 0.])
  g = 10.
  c = 1.
  #c = 0.1
  #c = 0
  bb = BB()
  p = bb.P(g, c)
  h = 0.01
  rx = 0.01
  n = np.infty
  #debug = True
  debug = False

  trjs = Euler(t, x, p, bb, h, rx, n, debug);

  t,x,o,p = stack(trjs)
  t = np.hstack(t)
  o = np.vstack(o)

  h,dh,KE,PE,Fn = o.T

  plt.figure(1)
  plt.clf()

  plt.subplot(2,1,1)
  plt.plot(dh,h,'.-')
  plt.plot(dh[:1],h[:1],'go')
  plt.plot(dh[-1:],h[-1:],'ro')
  plt.xlabel('dh')
  plt.ylabel('h')

  plt.subplot(2,1,2)
  plt.plot(t,h,'k.-')
  plt.plot(t,dh,'b.-')
  plt.plot(t,KE,'r.-')
  plt.plot(t,PE,'g.-')
  plt.legend(('h','dh','KE','PE'))
  plt.xlabel('t')
  plt.ylabel('')

