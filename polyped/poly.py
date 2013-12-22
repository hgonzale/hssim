"""
Polyped model

by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2011 -- 2013
"""

# system
import os
import time
# profiling
import cProfile as profile
# numerics, plotting
import numpy as np
import pylab as plt
import matplotlib as mpl
import scipy as sp
font = {'family' : 'serif',
        'size'   : 14}
mpl.rc('font', **font)
np.set_printoptions(precision=2)
# simulation
import relax as rx
# util
from util import Struct
import geom
import opt
# dbg
try:
  from IPython.Debugger import Tracer; dbg = Tracer()
except:
  dbg = lambda : 1/0

def skew( v ):
  """
  Convert a 3-vector to a skew matrix such that 
    dot(skew(x),y) = cross(x,y)
  
  The function is vectorized, such that:
  INPUT:
    v -- 3 x ... x N -- input vectors
  OUTPUT:
    3 x 3 x ... x N
  """
  v = np.asarray(v)
  z = np.zeros_like(v[0,...])
  return np.array([[ z,       -v[2,...], v[1,...]],
                   [ v[2,...], z,       -v[0,...]],
                   [-v[1,...], v[0,...], z,      ]]).T

def Adjoint(p,R):
  """
  Ad = Adjoint(p,R)

  Inputs:
    p - n x 1 - position
    R - n x n - orientation

  Outputs:
    Ad - 2n x 2n - adjoint transformation
       = ( R skew(p)*R )
         ( 0     R     )
  """
  n = p.size
  assert R.shape[0] == R.shape[1], "R is square"
  assert R.shape[0] == n, "R and p have compatible shapes"
  p = p.flatten()
  Ad = np.r_[np.c_[R, np.dot(skew(p),R)], np.c_[np.zeros((n,n)), R]]
  return Ad

def homog(p,R):
  """
  g = homog(p,R)
  
  Inputs:
    p - n x 1 - position
    R - n x n - orientation

  Outputs:
    g - (n+1) x (n+1) - homogeneous representation
      = ( R p )
        ( 0 1 ) 
  """
  n = p.size
  assert R.shape[0] == R.shape[1], "R is square"
  assert R.shape[0] == n, "R and p have compatible shapes"
  p = p.flatten()[:,np.newaxis]
  zo = np.r_[np.zeros(n),1.]
  g = np.r_[ np.c_[R,p], [zo]]
  return g

def unhomog(g):
  """
  p,R = unhomog(g)
  
  Inputs:
    g - ( R p )
        ( 0 1 ) - homogeneous representation

  Outputs:
    p - n x 1 - position
    R - n x n - orientation
  """
  assert g.shape[0] == g.shape[1], "g is square"
  n = g.shape[0]-1
  g = g[:n]
  p = g[:,n]
  R = g[:n,:n]
  return p,R

class Poly(rx.HDS):

  def __init__(self):
    """
    poly = Poly()

    poly - HDS - polyped hybrid model
      .F - vector field   dx  = F(t,x,p) dt
      .R - reset map      t,x,p = R(t,x,p)
      .E - retraction     x   = E(t,x,p,v)
      .P - parameters     p   = P(q,debug)
        .j - discrete mode index
      .G - guard function g = G(t,x,p)
      .O - observations   o = O(t,x,p)
      .plot(trjs) - generate  plot from trajectory

    by Sam Burden, Humberto Gonzalez, Ram Vasudevan, Berkeley 2013
    """
    rx.HDS.__init__(self)

  def Rot(self,t,x,p):
    """
    rot = Rot(...)

    rot - 3 x 3 - body rotation matrix in lab frame
    """
    th = x[3:6] # x,y,z rotation angles
    rot = geom.euler(th)
    return rot

  def Hip(self,t,x,p,pos=None,rot=None):
    """
    hip = Hip(...)

    hip - n x 3 - hip locations in world frame for each leg
    """
    # body pos,rot in world frame
    if pos is None:
      pos = x[0:3]
    if rot is None:
      rot = self.Rot(t,x,p)
    hip = pos + np.dot(p.hip,rot.T)
    return hip

  def Foot(self,t,x,p,pos=None,rot=None,debug=False):
    """
    foot = Foot(...)

    foot - n x 3 - foot locations in world frame for each leg
    """
    # body pos,rot in world frame
    if pos is None:
      pos = x[0:3]
    if rot is None:
      rot = self.Rot(t,x,p)
    # stance legs
    foot = p.foot.copy()
    # swing legs
    for l in (1.-p.b).nonzero()[0]:
      R = geom.euler([0.,-p.psi[l],-p.zeta[l]]).T
      f = p.hip[l] - np.dot( p.lamb[l]*p.ez, R.T )
      foot[l] = pos + np.dot( f, rot.T ) 
      if debug:
        print R, f, pos, rot
    #1/0
    return foot

  def P(self, debug=False,mu=np.identity(6),g=9.81,n=1,N=0,N_=np.infty, 
              a=np.ones(1), b=np.zeros(1), 
              hip=np.zeros((1,3)), foot=np.nan*np.ones((1,3)), 
              lamb=1.*np.ones(1), kapp=1e3*np.ones(1), 
              psi=np.zeros(1), zeta=np.zeros(1)):
    """
    Parameters:
      mu - 6x6 - SE(3) inertia tensor
      g - scalar - gravitational constant
      n - int - number of legs
      N - int - step count; negative indicates descent; zero disables count
      a - 1 x n - binary vector indicating leg active
      b - 1 x n - binary vector indicating leg stance
      hip - n x 3 - hip locations in body frame
      foot - n x 3 - foot stance locations in world frame
      lamb - 1 x n - leg lengths
      kapp - 1 x n - leg stiffness
      psi,zeta - 1 x n - leg touchdown angle is 
        dot( expm(zeta*skew(e_z)), expm(psi*skew(e_y)) )
      debug - bool - flag for printing debugging info
    """
    ex,ey,ez = np.identity(3)
    # TODO: sanitize inputs (make np.array, broadcast to correct size)
    return Struct(j=None,debug=debug,mu=mu,muinv=np.linalg.inv(mu),
                  g=g,n=n,N=N,N_=N_,a=a,b=b,
                  hip=hip,foot=foot,lamb=lamb,kapp=kapp,psi=psi,zeta=zeta,
                  ex=ex,ey=ey,ez=ez)

  def F(self,t,x,p):
    """
    dx = F(...)

    dx - 1 x 12 - vector field
    """
    q,dq = x[:6],x[6:]
    # leg data in world frame
    pos = q[0:3]; rot = self.Rot(t,x,p)
    hip  = self.Hip(t,x,p,pos=pos,rot=rot)
    foot = self.Foot(t,x,p,pos=pos,rot=rot)
    len  = np.sqrt( np.sum( (hip - foot)**2, axis=1 ) ).flatten()
    # hip wrenches on COM in world frame
    W = np.zeros((p.n,6))
    # active stance legs
    for l in (p.a * p.b).nonzero()[0]:
      f = p.kapp[l] * (p.lamb[l] - len[l]) * (hip[l] - foot[l]) / len[l]
      tau = np.zeros(3)
      w = np.r_[f,tau]
      #g_ch = ( h, inv( R ) ); g_hc = inv( g_ch )
      #g_hc = homog( -np.dot(rot, p.hip[l]), rot )
      Ad = Adjoint( -np.dot(rot, p.hip[l]), rot )
      W[l] = np.dot(w,Ad) # W = Ad^T * w
    # net wrench
    ddq = np.sum(W,axis=0) - p.g*np.dot(np.r_[p.ez,np.zeros(3)],p.mu)
    # vector field
    #dx = np.r_[ dq, np.dot( ddq, np.linalg.inv(p.mu) ) ]
    dx = np.r_[ dq, np.dot( ddq, p.muinv ) ]
    return np.array(dx)

  def G(self,t,x,p):
    """
    g = G(...)

    g > 0 : guard inactive
    g = 0 : guard 
    """
    g = np.zeros(p.n+2)
    # leg data in world frame
    q,dq = x[:6],x[6:]
    pos = q[0:3]; rot = self.Rot(t,x,p)
    hip  = self.Hip(t,x,p,pos=pos,rot=rot)
    foot = self.Foot(t,x,p,pos=pos,rot=rot)
    len  = np.sqrt( np.sum( (hip - foot)**2, axis=1 ) ).flatten()
    # foot heights
    hei = foot[:,2]
    # active stance legs
    for l in (p.a * p.b).nonzero()[0]:
      g[l] = p.lamb[l] - len[l]
    # active swing legs
    for l in (p.a * (1.-p.b)).nonzero()[0]:
      g[l] = hei[l]
    # test for local extrema in height above ground
    g[p.n] = np.sign( p.N ) * dq[2]
    # don't allow COM to touch ground
    g[p.n+1] = q[2]
    return g

  def R(self,t,x,p):
    """
    t,x,p = R(...)
    """
    t = t.copy(); x = x.copy(); p = p.copy()
    g = self.G(t,x,p)
    # leg data in world frame
    q,dq = x[:6],x[6:]
    pos = q[0:3]; rot = self.Rot(t,x,p)
    hip  = self.Hip(t,x,p,pos=pos,rot=rot)
    foot = self.Foot(t,x,p,pos=pos,rot=rot)
    # active legs undergoing transition
    for l in ((g[:p.n] < 0) * p.a).nonzero()[0]:
      # stance -> swing
      if p.b[l] > .5:
        # deactivate leg
        p.a[l] = 0.
        # leg liftoff
        p.b[l] = 0.
        # VIZ HACK: change leg angle, length to liftoff values
        z = np.sum( (foot[l] - hip[l]) * [1.,0.,1.j] ) * 1.j
        p.lamb[l] = np.abs(z)
        p.psi[l] = np.angle(z)
      # swing -> stance
      elif p.b[l] < .5:
        # leg touchdown
        p.b[l] = 1.
        # foot location in world frame
        p.foot[l] = foot[l]
    # count steps
    if g[p.n] < 0:
      p.N = -p.N - (1 + np.sign(p.N))/2
    # mode index
    if g[p.n+1] < 0 or np.abs(p.N) > p.N_:
      p.j = None
    else:
      p.j = ( p.N, np.sum( (2**np.arange(p.n))*p.b*p.a ) )
    return t,x,p

  def E(self,t,x,p,v):
    """
    x_ = E(...)
    """
    x_ = x + v
    return x_

  def O(self,T,X,p):
    """
    o = O(...)
    """
    n = T.size
    q,dq = X[:,:6],X[:,6:]
    pos = np.nan*np.zeros((n,3))
    rot = np.nan*np.zeros((n,3,3))
    hip  = np.nan*np.zeros((n,p.n,3))
    foot = np.nan*np.zeros((n,p.n,3))
    len = np.nan*np.zeros((n,p.n))
    E = np.nan*np.zeros((n,4))
    for k,(t,x) in enumerate(zip(T,X)):
      # leg data in world frame
      pos[k] = x[0:3]; rot[k,...] = self.Rot(t,x,p)
      hip[k,...]  = self.Hip(t,x,p,pos=pos[k],rot=rot[k,...])
      foot[k,...] = self.Foot(t,x,p,pos=pos[k],rot=rot[k,...])
      len[k,...] = np.sqrt(np.sum((hip[k,...]-foot[k,...])**2,axis=1)).flatten()
      # kinetic energy
      E[k,0:1] = .5*np.sum( np.dot(dq[k,...], p.mu) * dq[k,...] )
      # potential energy from gravity and active stance legs
      E[k,1:2] = ( p.mu[2,2]*p.g*q[k,2] 
                   +.5*np.sum(p.a*p.b*p.kapp*(p.lamb-len[k,...])**2) )
      E[k,2:3] = p.mu[2,2]*p.g*q[k,2] 
      E[k,3:4] = .5*np.sum(p.a*p.b*p.kapp*(p.lamb-len[k,...])**2)
    return dict(q=q,dq=dq,len=len,E=E,pos=pos,rot=rot,hip=hip,foot=foot)

  def saggital(self,X):
    """
    X_ = saggital(...)
    """
    q,dq = X[:6],X[6:]
    x_,_,z_,_,thy_,thz = q
    dx,dy,dz_,_,dthy_,_ = dq
    rot = geom.euler([0.,0.,thz])
    dx_,_,_ = np.dot([dx,dy,0.],rot)
    X_ = np.r_[x_,z_,thy_,dx_,dz_,dthy_]
    return X_

  def unsaggital(self,X_,y=0.,dy=0.,thx=0.,dthx=0.,thz=0.,dthz=0.):
    """
    X_ = saggital(...)
    """
    x,z,thy,dx,dz,dthy = X_
    X = np.r_[x,y,z,thx,thy,thz,dx,dy,dz,dthx,dthy,dthz]
    return X

  def intrinsic(self,X):
    """
    X_ = intrinsic(...)
    """
    q,dq = X[:6],X[6:]
    x,y,z_,thx_,thy_,thz = q
    dx,dy,dz_,dthx_,dthy_,dthz_ = dq
    rot = geom.euler([0.,0.,thz])
    dx_,dy_,_ = np.dot([dx,dy,0.],rot)
    X_ = np.r_[z_,thx_,thy_,dx_,dy_,dz_,dthx_,dthy_,dthz_]
    return X_

  def extrinsic(self,X_,x=0.,y=0.,thz=0.):
    """
    X = extrinsic(X_)
    """
    z,thx,thy,dx_,dy_,dz,dthx,dthy,dthz = X_
    rot = geom.euler([0.,0.,thz])
    dx,dy,_ = np.dot([dx_,dy_,0.],rot.T)
    X = np.r_[x,y,z,thx,thy,thz,dx,dy,dz,dthx,dthy,dthz]
    return X

  def plot(self,trjs,fn=None,dpi=None,pxw=1000,pxh=600,offv=0,figN=1,
           plots=['pos','rot','eng'],
           alpha=1.,col=('b','g','r'),lp='',lw=4,slab='',
           clf=True,legend=True,fill=True,legloc='lower left'):
    """
    Plot trajectory

    Inputs:
      trjs - list - Poly trajectory

    Outputs:
      (T,Z,O,P),(Te,Ze,Oe,Pe) - observations
    """

    def fill(Te,Pe,ylim):
      for te,pe in zip(Te,Pe):
        if fill:
          fc = 1 - .33*np.sum( pe.b ) / pe.n
          ax.fill([te[0],te[0],te[1],te[1]],
                  [ylim[1],ylim[0],ylim[0],ylim[1]],
                  fc=np.ones(3)*fc,ec='none',zorder=-1)
                  #fc=np.ones(3)*fc,ec='k',zorder=-1)
          #ax.text(.5*(te[0]+te[1]),.5*(ylim[0]+ylim[1]),np.sum( pe.b ),
          #        horizontalalignment='center',verticalalignment='center')

    Te,Xe,Oe,Pe = rx.trans(trjs)
    T,X,O,P = rx.stack(trjs)
    vars = O[0].keys()
    o = Struct(**dict([( v, np.vstack([oo[v] for oo in O]) ) for v in vars]))

    legprop = {'size':12};
    ms = 20.; mew = 2.; ma = 0.5

    t = np.hstack(T)
    xlim = (min(t),max(t))

    fig = plt.figure(figN)
    if clf:
      fig.clf()

    sps = len(plots); sp = 1;

    if 'pos' in plots:
      ax = plt.subplot(sps,1,sp); sp += 1; plt.grid('on')
      H = (ax.plot(t,o.q[:,0],col[0]+lp,label='$x'+slab+'$',alpha=alpha),
           ax.plot(t,o.q[:,1],col[1]+lp,label='$y'+slab+'$',alpha=alpha),
           ax.plot(t,o.q[:,2],col[2]+lp,label='$z'+slab+'$',alpha=alpha))
      [[hh.set_linewidth(lw) for hh in h] for h in H]
      ax.plot(t,o.hip[:,:,0],'k:')
      ax.plot(t,o.hip[:,:,1],'k:')
      ax.plot(t,o.hip[:,:,2],'k:')
      ax.plot(t,o.foot[:,:,0],'k--')
      ax.plot(t,o.foot[:,:,1],'k--')
      ax.plot(t,o.foot[:,:,2],'k--')

      ylim = np.array(ax.get_ylim())#-0.25
      fill(Te,Pe,ylim)
      if clf:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
      ax.set_ylabel('position')

      if legend:
        leg = ax.legend(ncol=6,loc=legloc,prop=legprop)

    if 'rot' in plots:
      ax = plt.subplot(sps,1,sp); sp += 1; plt.grid('on')
      H = (ax.plot(t,o.q[:,3],col[0]+lp,label='$\\theta_x'+slab+'$',alpha=alpha),
           ax.plot(t,o.q[:,4],col[1]+lp,label='$\\theta_y'+slab+'$',alpha=alpha),
           ax.plot(t,o.q[:,5],col[2]+lp,label='$\\theta_z'+slab+'$',alpha=alpha))
      [[hh.set_linewidth(lw) for hh in h] for h in H]

      ylim = np.array(ax.get_ylim())#-0.25
      fill(Te,Pe,ylim)
      if clf:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
      ax.set_ylabel('rotation')

      if legend:
        leg = ax.legend(ncol=6,loc=legloc,prop=legprop)

    if 'eng' in plots:
      ax = plt.subplot(sps,1,sp); sp += 1; plt.grid('on')
      H = (ax.plot(t,o.E[:,0],col[0]+lp,label='$T'+slab+'$',alpha=alpha),
           ax.plot(t,o.E[:,3],col[1]+lp,label='$V_k'+slab+'$',alpha=alpha),
           ax.plot(t,o.E[:,2],col[2]+lp,label='$V_g'+slab+'$',alpha=alpha),
           ax.plot(t,o.E[:,0]+o.E[:,1],'k'+lp,label='$L'+slab+'$',alpha=alpha))
      [[hh.set_linewidth(lw) for hh in h] for h in H]

      ylim = np.array(ax.get_ylim())#-0.25
      fill(Te,Pe,ylim)
      if clf:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
      ax.set_ylabel('energy')

      if legend:
        leg = ax.legend(ncol=6,loc=legloc,prop=legprop)

    ax.set_xlabel('time (t)')

    if fn is not None:
      plt.savefig(fn,bbox_inches='tight',pad_inches=0.1,dpi=dpi)
    
    return o

  def sch(self,t,x,p,figN=2,text=False,fn='',exts=['pdf'],plots=['sch','grd'],
          clf=True,figs={},axs={},xlim=None,ylim=None,alpha=1.,lw=6,
          fill=True,light=0.,time=False,offz=0,offx=0.,offt=+1.,
          position=[0.2,0.2,0.60,0.6],dpi=80,pxh=600,offv=0,pad_inches=.1):
    """
    Plot schematic

    Inputs:
      t,x,p - time, state, parameters
    """

    def circle(z0=0.,r=1.,th=np.linspace(0.,2*np.pi,100)):
      """
      z = circle(...)

      Inputs:
        z0 - complex - circle center
        r  - scalar - circle radius
        th - 1 x n - sweep angles

      Outputs:
        z - 1 x n complex - points on circle
      """
      z = z0 + r*np.exp(1.j*th)
      return z

    def zigzag(a=.2,b=.6,c=.2,s=.2,p=4,N=100,z0=None,z1=None):
      """
      z = zigzag(...)

      Inputs:
        (optional)
        a - float - length before zigzag
        b - float - length of zigzag
        c - float - length after zigzag
        s - float - size of zigzags
        p - int - number of zigzags
        N - int - number of samples 
        z0,z1 - complex - endpoints for zigzag; neither or both must be given

      Outputs:
        z - complex N vector - x + 1.j*y pairs along zigzag
            x[0] = y[0] = 0; x[N-1] = a+b+c; y[N-1] = 0
      """
      x = np.linspace(0.,a+b+c,N); y = 0.*x
      mb = np.round(N*a/(a+b+c)); Mb = np.round(N*(a+b)/(a+b+c))
      y[mb:Mb] = s*( np.mod( np.linspace(0.,p-.01,Mb-mb), 1. ) - 0.5 )
      z = x + 1.j*y
      if z0 is not None and z1 is not None:
        z_ = z
        z = z0 + (z1-z0)*z
        #1/0
      return z

    lighten = lambda c : c + (1. - np.asarray(c)) * light

    o = Struct(**self.O(np.array([t]),np.array([x]),p))

    ms = 40.; mew = 2.; ma = 0.5
    fs = 22.
    fd = {'family':'serif','size':fs}
    mpl.rc('font', **font)
    #mpl.rc('text.latex', preamble='\usepackage[helvet]{sfmath}')

    # leg data in world frame
    q,dq = o.q[0],o.dq[0]
    pos = q[0:3]; rot = self.Rot(t,x,p)
    hip  = self.Hip(t,x,p,pos=pos,rot=rot)
    foot = self.Foot(t,x,p,pos=pos,rot=rot)
    len  = np.sqrt( np.sum( (hip - foot)**2, axis=1 ) ).flatten()

    e3 = np.identity(3)[0]
    e6 = np.identity(6)[0]
    pos += offx*e3
    hip += offx*e3; foot += offx*e3

    mw = .5*p.lamb.max();
    mh = .2*mw;

    if xlim is None:
      xlim = q[0] + np.array([-1.5*mw,2.5*mw])
    if ylim is None:
      ylim = (-0.1,1.6*p.lamb.max())
    dxl = np.diff(xlim)[0]; dyl = np.diff(ylim)[0]

    pxw = int(pxh*np.diff(xlim) / np.diff(ylim))
    pxw = 2*int(pxw/2)
    pxh = 2*int(pxh/2)

    if 'sch' in plots:
      if clf:
        fig = plt.figure(figN)
        fig.clf()
        ax = plt.axes(position)
        ax.grid(zorder=-1)
      else:
        fig = figs['sch']
        ax = axs['sch']

      gc = lighten(np.array([139.,69.,19.])/255.)
      hc = lighten(np.array([1.0,0.4,0.4]))
      mc = lighten(np.array([0.4,0.4,1.]))
      lc = lighten(np.array([0.,0.5,0.]))

      # m
      z_ = np.asarray([-mw,-mw,mw,mw,-mw])+1.j*np.asarray([-mh,mh,mh,-mh,-mh])
      z = (q[0] + 1.j*q[2]) + z_*np.exp(1.j*q[4])
      if fill:
        ax.fill(z.real,z.imag,lw=lw,ec=.5*mc,fc=mc,alpha=alpha,zorder=5+offz)
      else:
        ax.plot(z.real,z.imag,lw=lw,color=.5*mc,alpha=alpha,zorder=5+offz)
      if time:
        ax.plot(q[0],q[2],'.',mfc=mc,mec=.5*mc,mew=lw,ms=ms,zorder=6+offz,alpha=alpha)
      if time and not( light == 1. ):
        # time
        ax.text(q[0],ylim[1]-.30*mw,'$t=%0.1f$'%t,ha='center',va='bottom',zorder=10+offz,backgroundcolor='w',**fd)#bbox=dict(facecolor='w',edgecolor='k',lw=1.),**fd)

        ax.plot([q[0],q[0]],[ylim[1]-.25*mw,q[2]],'--',color=.3*np.ones(3),lw=lw/2.,zorder=1)
      if text:
        # m
        ax.text(q[0],q[2],'$m,I$',ha='center',va='center',zorder=10,**fd)
        # d
        z = hip[[1,0],0] + 1.j*hip[[1,0],2]; dz = np.diff(z)[0]; 
        z += 1.3j*.3*mw*dz
        z_ = z.mean()
        ax.arrow(z_.real,z_.imag,.4*dz.real,.4*dz.imag,
                 head_width=0.06,lw=.5*lw,fc='k')
        ax.arrow(z_.real,z_.imag,-.4*dz.real,-.4*dz.imag,
                 head_width=0.06,lw=.5*lw,fc='k')
        ax.text(z_.real,z_.imag+.2*mh,'$d$',ha='center',va='bottom',**fd)
        # rotation
        lamb = 1.2*p.lamb.max()
        z0 = q[0]+1.j*q[2]; z1 = z0 + lamb; z = np.asarray([z0,z1])
        th = np.linspace(+1.0*np.pi/2**7,q[4]-1.25*np.pi/2**5,10)
        a = circle(z0=z0,r=1.02*lamb,th=th); da = 1.j*np.exp(1.j*th)*.05
        ax.plot(z.real,z.imag,'k--',lw=.5*lw,zorder=4)
        z_ = np.asarray([z0,z0+(z1-z0)*np.exp(1.j*q[4])])
        ax.plot(z_.real,z_.imag,'k--',lw=.5*lw,zorder=4)
        ax.plot(a.real,a.imag,'k',lw=.5*lw,zorder=4)
        ax.arrow(a[-1].real,a[-1].imag,da[-1].real,da[-1].imag,
                 head_width=0.06,lw=.5*lw,fc='k')
        ax.text(a[a.size/2].real-.125*lamb,a[-1].imag,'$\\theta$',va='center',**fd)
        # g
        ax.arrow(xlim[0]+.05*dxl,ylim[1]-.1*dyl,0.,-3*mh,
                 head_width=0.06,lw=lw,fc='k')
        ax.text(xlim[0]+.03*dxl,ylim[1]-.1*dyl-mh,'$g$',ha='right',va='center',**fd)
      for l in range(p.n):
        # hip
        if True:#fill:
          ax.plot(hip[l,0],hip[l,2],'.',mfc=hc,mec=.5*hc,
                  mew=lw,ms=ms,zorder=6+offz,alpha=alpha)
        else:
          ax.plot(hip[l,0],hip[l,2],'.',mfc='none',
                  mec=.5*hc,mew=lw,ms=ms,zorder=6+offz,alpha=alpha)
        if p.b[l]:
          # foot
          ax.plot(foot[l,0],foot[l,2],'.',mfc=hc,
                  mec=.5*hc,mew=lw,ms=ms,zorder=6+offz)
        #if p.a[l]:
        #  alpha = alpha*1.
        #else:
        #  alpha = alpha*.5

        z0 = hip[l,0] + 1.j*hip[l,2]
        z1 = foot[l,0] + 1.j*foot[l,2]
        s = .1*p.lamb[l] / np.abs(z1-z0)
        z = zigzag(s=s,z0=z0,z1=z1,p=6)
        ax.plot(z.real,z.imag,lw=1.0*lw,color=lc,alpha=alpha,zorder=5+offz)
        if text and l == 1:
          # leg spring
          ax.text(z.real.mean()+.06*dxl,z.imag.mean()+.01*dyl,'$k, \\ell$',va='center',**fd)
          # leg angle
          th = q[4] + np.linspace(np.pi/2**6,p.psi[l]-np.pi/2**4,10)-np.pi/2.
          a = circle(z0=z[0],r=p.lamb[l],th=th); da = 1.j*np.exp(1.j*th)*.05
          z_ = np.asarray([z0,z0+(z1-z0)*np.exp(-1.j*p.psi[l])])
          ax.plot(z_.real,z_.imag,'k--',lw=.5*lw,zorder=4)
          ax.plot(a.real,a.imag,'k',lw=.5*lw,zorder=4)
          ax.arrow(a[-1].real,a[-1].imag,da[-1].real,da[-1].imag,
                   head_width=0.06,lw=.5*lw,fc='k')
          ax.text(a[a.size/3].real,a[a.size/3].imag+.05*p.lamb[l],'$\\psi$',va='bottom',**fd)

      if clf:
        ax.axis('equal')
        if pad_inches > 0:
          ax.set_yticks([0.,.5,1.,1.5])
          ax.set_xlabel('$x$',**fd)
          ax.set_ylabel('$z$',rotation='horizontal',**fd)
          ax.tick_params(labelsize=fs)
        else:
          ax.set_xticklabels([])
          ax.set_yticklabels([])
          #ax.set_xticks([])
          #ax.set_yticks([])
          fig.set_size_inches(pxw/dpi,pxh/dpi)
        ax.set_xlim(xlim); 
        ax.set_ylim(ylim); 

        # ground
        ax.fill([xlim[0]-.1*dxl,xlim[0]-.1*dxl,xlim[1]+.1*dxl,xlim[1]+.1*dxl],
                [0.,ylim[0]-.1*dyl,ylim[0]-.1*dyl,0.],
                lw=lw,ec=0.5*gc,fc=gc,zorder=5)

      if not fn == '':
        di,fi = os.path.split(fn)
        if not os.path.exists(di):
          os.mkdir(di)
        for ext in exts:
          fig.savefig(fn+'_sch.'+ext,bbox_inches='tight',pad_inches=pad_inches,dpi=dpi)

      figs['sch'] = fig; axs['sch'] = ax

    if 'grd' in plots:
      # grd
      lamb = p.lamb.max(); psi = p.psi.min()

      xlim = (-np.pi/8-np.pi/12,np.pi/8+np.pi/12)
      ylim = (0.2,1.25)
      dxl = np.diff(xlim)[0]; dyl = np.diff(ylim)[0]
      pxw = int(pxh*np.diff(xlim) / np.diff(ylim))

      fig = plt.figure(3); fig.clf(); 
      ax = plt.axes([0.2,0.2,0.6,0.6])
      ax.grid(zorder=-1)

      li = .5
      n = 400
      th = np.linspace(-2.0*np.pi/2,2.0*np.pi/2,n)
      d = p.hip[:,0:1]
      z0 = -d*np.sin(th) + lamb*np.cos(th + psi)
      G = {'a,l': dict(color='b'),
           'l,g': dict(color='b',dashes=[10,5]),
           'a,r': dict(color='g',dashes=[10,5]),
           'r,g': dict(color='g')}
      if False:
        ax.plot(th,z0[0],'b',lw=lw)
        ax.plot(th,z0[1],'g',lw=lw,dashes=[10,5])
        ax.fill(th,z0[0],hatch='\\',
                ec='none',zorder=1,fc=np.array([li,li,1.0]),label='$D_l$')
        ax.fill(th,z0[1],hatch='/',
                ec='none',zorder=1,fc=np.array([li,1.0,li]),label='$D_r$')
        ax.fill(np.hstack((th[0],th,th[-1])),
                np.hstack((2.,z0[0,:n/2],z0[1,n/2:],2.)),
                ec='none',zorder=1,fc=np.array([1.0,1.0,1.0]),label='$D_a$')
        ax.fill(th,np.hstack((z0[1,:n/2],z0[0,n/2:])),hatch='x',
                ec='none',zorder=1,fc=np.array([1.0,1.0,0.0]),label='$D_g$')
      else:
        #ax.plot(th,z0[0],'b',lw=lw)
        #ax.plot(th,z0[1],'g',lw=lw,dashes=[10,5])
        ax.plot(th[:n/2],z0[0][:n/2],lw=lw,**G['a,l'])
        ax.plot(th[:n/2],z0[1][:n/2],lw=lw,**G['l,g'])
        ax.plot(th[n/2:],z0[1][n/2:],lw=lw,**G['a,r'])
        ax.plot(th[n/2:],z0[0][n/2:],lw=lw,**G['r,g'])
        ax.fill(th,z0[0],
                ec='none',zorder=1,fc=np.array([li,li,1.0]))
        ax.fill(th,z0[1],
                ec='none',zorder=1,fc=np.array([li,1.0,li]))
        ax.fill(np.hstack((th[0],th,th[-1])),
                np.hstack((2.,z0[0,:n/2],z0[1,n/2:],2.)),
                ec='none',zorder=1,fc=np.array([1.0,1.0,1.0]))
        ax.fill(th,np.hstack((z0[1,:n/2],z0[0,n/2:])),
                ec='none',zorder=1,fc=np.array([1.0,1.0,0.0]))
        pad = 10
        bbox = dict(facecolor='w',pad=pad,lw=0.)
        thl = -1.25*np.pi/16
        zl = -d[0]*np.sin(thl) + lamb*np.cos(thl + psi)
        thr = +1.25*np.pi/16
        zr = -d[1]*np.sin(thr) + lamb*np.cos(thr + psi)
        ax.text(-1*np.pi/6,0.5,'$D_g$',ha='center',va='center',bbox=bbox,**fd)
        ax.text(-1*np.pi/6,zl,'$D_l$',ha='center',va='center',bbox=bbox,**fd)
        ax.text(+1*np.pi/6,0.5,'$D_r$',ha='center',va='center',bbox=bbox,**fd)
        ax.text(+1*np.pi/6,zl,'$D_a$',ha='center',va='center',bbox=bbox,**fd)
        ax.text(thl,zl,'$G_{(a,l)}$',ha='center',va='center',
                bbox=dict(facecolor='w',edgecolor=G['a,l']['color'],lw=3.,boxstyle='square'),**fd)
        ax.text(thl,0.8,'$G_{(l,g)}$',ha='center',va='center',
                bbox=dict(facecolor='w',edgecolor=G['l,g']['color'],lw=3.,boxstyle='square',ls='dashed'),**fd)
        ax.text(thr,zr,'$G_{(a,r)}$',ha='center',va='center',
                bbox=dict(facecolor='w',edgecolor=G['a,r']['color'],lw=3.,boxstyle='square',ls='dashed'),**fd)
        ax.text(thr,0.5,'$G_{(r,g)}$',ha='center',va='center',
                bbox=dict(facecolor='w',edgecolor=G['r,g']['color'],lw=3.,boxstyle='square'),**fd)

      #ax.arrow(0,ylim[1],0,-.6*dyl,lw=lw,fc='k',head_width=.05,zorder=10)

      #ax.legend(ncol=3,prop={'size':18},columnspacing=1,handletextpad=.2)
      ax.axis('equal')
      ax.set_xticks(np.arange(-2,3)*np.pi/8)
      ax.set_xticklabels(['$-\\pi/4$','$-\\pi/8$','$0$','$+\\pi/8$','$+\\pi/4$'])
      ax.set_yticks([0.,.5,1.,1.5])
      ax.set_xlim(xlim); ax.set_xlabel('$\\theta$',**fd);#ax.set_xticks([])
      ax.set_ylim(ylim); ax.set_ylabel('$z$',rotation='horizontal',**fd);#ax.set_yticks([])
      ax.tick_params(labelsize=fs)
      #ax.set_title('$\\dot{z} < 0,\\ \\dot{\\theta} = 0,\\ \\ell = 1,\\ \\psi = \\pi/%d$' % int(psi**-1 * np.pi),**fd)

      if not fn == '':
        for ext in exts:
          fig.savefig(fn+'_grd.'+ext,bbox_inches='tight',pad_inches=0.1,dpi=dpi)

      figs['grd'] = fig; axs['grd'] = ax

    return figs,axs

def step(op,plot=True):
  """
  step(op) 
  """
  p = Struct(**op.p)

  st = time.clock()
  trjs = rx.Euler((p.t0,p.tf), p.x0, p.p0, p.poly, p.h, p.eps, p.n, p.debug)
  print '%0.2f sec for sim' % (time.clock() - st)

  trj0 = trjs[0]
  trj1 = trjs[-1]

  obs = None
  if plot:
    obs = p.poly.plot(trjs)

  return trjs,obs


if __name__ == "__main__":

  import sys
  args = sys.argv

  assert len(args) > 1, 'specify .cfg file'

  op = opt.Opt()

  if '--profile' in args:
    args.pop(args.index('--profile'))
    fi = args[1]
    op.pars(fi=fi,Poly=Poly)
    profile.run("trjs = step(op,plot=False)", fi+'.pstats')
    cmd = 'gprof2dot.py -f pstats {0}.pstats | dot -Tsvg -o {0}.svg'.format(fi)
    print cmd
    os.system(cmd)
  else:
    fi = args[1]
    op.pars(fi=fi,Poly=Poly)
    trjs,obs = step(op,plot=True)

