# system
import os
import time
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
import poly
# util
from util import Struct
import opt

import sys
args = sys.argv

assert len(args) > 1, 'specify .cfg file'

sig = 0.
if len(args) > 2:
  sig = 2e-1
  seed = int(args[2])
  np.random.seed(int(seed))
  print 'seed = %d' % seed

offv = 400

op = opt.Opt()
op.pars(fi=args[1],Poly=poly.Poly)
p = Struct(**op.p)

def step(x0_,p=p):
  x_,z_,th_,dx_,dz_,dth_ = x0_

  t0 = p.t0
  x0 = p.poly.unsaggital([x_,z_,th_,dx_,dz_,dth_])
  p0 = p.p0.copy()
  trjs = rx.Euler((t0,np.infty), x0, p0, p.poly, p.h, p.eps, p.n, p.debug)

  x1 = trjs[-1].x[-1]
  x,y,z,thx,thy,thz,dx,dy,dz,dthx,dthy,dthz = x1
  assert np.allclose([y,thx,thz,dy,dthx,dthz],np.zeros(6)), "saggital dynamics"

  x1_ = p.poly.saggital(x1)

  return x1_,trjs

x0_ = p.poly.saggital(p.x0)

debug = True
#debug = False

savefig = True
#savefig = False


st = time.clock()
x_ = x0_
_,z0_,_,dx0_,_,_ = x0_
#(z_,dx_) = rx.fp( f, [z0_,dx0_], modes=[3], eps=1e-2, debug=debug )
dx_ = dx0_
f = lambda dx : step([0.,z0_,0.,dx,0.,0.],p=p)[0][3]
#def armijo(f,x,d,fx=None,Dfx=None):
#  return rx.armijo(f,x,d,gamma=1e-0)
#dx_ = rx.descent1(lambda x : np.sum((f(x) - x)**2),[dx0_],
#                     a=armijo,eps=1e-2,debug=debug)[0][0]
x_ = np.asarray([0.,z0_,0.,dx_,0.,0.])
if savefig:
  x_ += [0,0.05,np.pi/2**4,0,0,0] 
else:
  x_ += sig*np.random.randn(6)*[1,1,.6,1,0,1]
x__,trjs = step(x_,p=p)
fn = ''
#if savefig:
#  fn = 'fig/pronk_trj.pdf'
p.poly.plot(trjs,offv=offv,fn=fn,plots=['pos','eng'])
print '%0.2f sec; dx0 = %0.2f, dx1 = %0.2f' % (time.clock() - st,x_[3],x__[3])

text = True
#text = False

fps = 60.
if text:
  fps = 1e-10
n_  = 0.
fn = ''
if savefig:
  fn = 'fig/pronk'
exts = ['pdf','eps']
#exts = ['pdf']
#exts = ['eps']
#exts = []


if text:
  figs,axs = p.poly.sch(trjs[0].t[0],trjs[0].x[0],trjs[0].p,offv=offv,
                        text=text,fn=fn,plots=['sch','grd'],exts=exts)

  def tweak_grd(th0,dth0,ls='r',dashes=None,lw=4.,zorder=5,label=''):
    st = time.clock()
    x_ = np.asarray([0.,z0_,th0,dx_,0.,dth0])
    x__,trjs = step(x_,p=p)
    print '%0.2f sec; dx0 = %0.2f, dx1 = %0.2f' % (time.clock() - st,x_[3],x__[3])

    T,X,O,P = rx.stack(trjs)
    vars = O[0].keys()

    t = np.hstack(T)
    # observations until max compression
    o = Struct(**dict([( v, np.vstack([oo[v] for oo,pp in zip(O,P) if (pp.j is not None) and (pp.j[0] < 0)]) ) for v in vars]))

    z,th = o.q[:,[2,4]].T
    if dashes is None:
      axs['grd'].plot(th,z,ls,lw=lw,zorder=zorder,label=label)
    else:
      axs['grd'].plot(th,z,ls,lw=lw,zorder=zorder,label=label,dashes=dashes)
    axs['grd'].arrow(th[-20],z[-20],np.diff(th[-20:])[0],np.diff(z[-20:])[0],
                     lw=lw,fc='k',head_width=.04,zorder=10)

  d = .4 
  tweak_grd(0.,0.,ls='k',label='$\\dot{\\theta}_0=0$')
  tweak_grd(0.,-d,ls='k--',label='$\\dot{\\theta}_0=-%0.1f$'%d,dashes=[10,5])
  #tweak_grd(0.,+np.pi/d,ls='k--',label='$\\dot{\\theta}=+\\pi/%d$'%d,dashes=[3,3])
  #tweak_grd(-np.pi/d,ls='k--',label='$+\\pi/%d$'%d)
  #tweak_grd(+np.pi/d,ls='k--',label='$-\\pi/%d$'%d)

  def key((h,l)):
    if 'G' in l:
      return 0
    elif 'D_l' in l or 'D_r' in l:
      return 2
    elif 'D_g' in l or 'D_a' in l:
      return 4
    else:
      if '0' in l:
        return 1
      elif '+' in l:
        return 3
      elif '-' in l:
        return 5

  import operator
  h_,l_ = axs['grd'].get_legend_handles_labels()
  hl = sorted(zip(h_,l_),key=key)
  h,l = zip(*hl) #key=operator.itemgetter(1)))

  #axs['grd'].legend(h,l,ncol=2,prop={'size':16},columnspacing=1,handletextpad=.2,handleheight=1.5)
  #axs['grd'].legend(h,l,ncol=1,prop={'size':16},columnspacing=1,handletextpad=.2,handleheight=1.5)

  fs = 16.
  fd = {'family':'serif','size':fs}
  #bbox = dict(facecolor='w',edgecolor='k',pad=8)
  bbox = dict(facecolor='w',edgecolor='k',lw=2.,boxstyle='square')
  axs['grd'].text(0.,0.5,'$\\dot{\\theta}_0=0$',ha='center',va='center',bbox=bbox,zorder=6,**fd)
  #bbox = dict(facecolor='w',edgecolor='k',pad=8)#,dashes=[10,5])
  bbox = dict(facecolor='w',edgecolor='k',lw=2.,boxstyle='square',ls='dashed')
  axs['grd'].text(-np.pi/20,0.625,'$\\dot{\\theta}_0=-%0.1f$'%d,ha='center',va='center',bbox=bbox,zorder=6,**fd)

  if not fn == '':
    for ext in exts:
      figs['grd'].savefig(fn+'_grd.'+ext,bbox_inches='tight',pad_inches=0.1)


sys.exit(0)

for trj in trjs:
  for t,x in zip(trj.t,trj.x):
    if t >= p.t0 + n_/fps:
      p.poly.sch(t,x,trj.p,offv=offv,text=text,fn='vid/pronk_%04d'%n_,
                 plots=['sch'],exts=['png'],position=[0,0,1,1],pad_inches=0.)
      plt.draw()
      n_ += 1
  #time.sleep(1.)

os.system('ffmpeg -r %d -sameq -i vid/pronk_%%04d_sch.png vid/pronk.mp4' % fps)
#plt.figure(3); plt.clf();

#fx = f(x_)
#e = np.arange(-9,-2)
#D = [rx.jac(f, x_, fx=fx, d=d) for d in np.power(10.,e)]
#E = [np.linalg.eigvals(d) for d in D]

#plt.clf()
#plt.plot(e,E,'-.')

#v = np.array([  4.91e-01,   5.05e-06,  -5.31e-03,  -1.98e-01,   8.46e-01, 2.87e-06,  -5.58e-02])

