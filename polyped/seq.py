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
dx_ = dx0_
f = lambda dx : step([0.,z0_,0.,dx,0.,0.],p=p)[0][3]
x_ = np.asarray([0.,z0_,0.,dx_,0.,0.])
fn = ''

text = True
text = False

di = 'seq/'
exts = ['pdf','eps']
#exts = ['eps']
#exts = []

position = [0,0,1.00,1.00]
position = [0.2,0.2,0.6,0.6]

fill = False
fill = True

def seq(x_,p=p,fn='',xlim=None,
        props=dict(position=position,offv=550,text=text,pad_inches=0.,
                   light=0.,time=False,fn='',plots=['sch'],exts=exts,
                   alpha=0,offz=0,clf=True,ylim=(-0.1,1.45),fill=fill) ):
  x__,trjs = step(x_,p=p)

  props = props.copy()
  
  T,X,O,P = rx.stack(trjs)
  vars = O[0].keys()
  o = Struct(**dict([( v, np.vstack([oo[v] for oo in O]) ) for v in vars]))

  if len(trjs) == 7:
    trjs.pop(5); trjs.pop(1)

  n = len(trjs) / 2
  offx = -trjs[n].p.foot[:,0].mean()

  # enclose whole trajectory in figure
  #props.update(xlim=.7*trjs[0].p.lamb[0]*np.asarray([-1.,+1.])+[trjs[0].x[0,0],trjs[-1].x[0,0]])
  if xlim is None:
    xlim = (.7*trjs[0].p.lamb[0]*np.asarray([-1.,+1.])
            +[trjs[0].x[0,0],trjs[-1].x[0,0]]+offx)
  props.update(xlim=xlim,offx=offx)
  figs,axs = p.poly.sch(trjs[n].t[0],trjs[n].x[0],trjs[n].p,**props)
  props.update(figs=figs,axs=axs,clf=False,alpha=1,offz=0,lw=6,time=True)

  for k,trj in enumerate(trjs):
    #props.update(alpha=float(k+1)/len(trjs))
    if k+1 == len(trjs):
      props.update(fn=di+fn)
    props.update(offz=10*k,light=.9*(1.-float(k+1)/len(trjs)),offt=+1 - 2*(k == n))
    p.poly.sch(trj.t[0]-trjs[n].t[0],trj.x[0],trj.p,**props)
    #axs['sch'].plot(trj.x[:,0]+offx,trj.x[:,2],color=[.6,.6,.6],lw=6,zorder=1)

  axs['sch'].plot(o.q[:,0]+offx,o.q[:,2],color=.3*np.ones(3),lw=4,zorder=10*(k+1))
  axs['sch'].arrow(o.q[-1,0]+offx,o.q[-1,2],0.10,0.,head_width=0.06,lw=4,color=.3*np.ones(3),zorder=10*(k+1))
  plt.draw()
  return trjs,figs,axs,props

#seq(x_)

#d = np.asarray([0,0,np.pi,0,0,0]) / 2**5
d = np.asarray([0,0,0,0,0,1])

trjs,figs,axs,props = seq(x_ + 0.0*d,fn='pronk')
#trjs,figs,axs,props = seq(x_ + 0.4*d,fn='plus',xlim=props['xlim'])
#trjs,figs,axs,props = seq(x_ - 0.2*d,fn='minus',xlim=props['xlim'])
trjs,figs,axs,props = seq(x_ - 0.4*d,fn='minus')

#fn = 'plus'
#x__,trjs = step(x_+d,p=p)
#for k,trj in enumerate(trjs):
#      p.poly.sch(trj.t[0],trj.x[0],trj.p,offv=offv,text=text,fn=di+fn+'%01d' % k,
#                 plots=['sch'],exts=exts,position=position,pad_inches=0.)
#      plt.draw()

#fn = 'minus'
#x__,trjs = step(x_-d,p=p)
#for k,trj in enumerate(trjs):
#      p.poly.sch(trj.t[0],trj.x[0],trj.p,offv=offv,text=text,fn=di+fn+'%01d' % k,
#                 plots=['sch'],exts=exts,position=position,pad_inches=0.)
#      plt.draw()

#os.system('ffmpeg -r %d -sameq -i vid/pronk_%%04d_sch.png vid/pronk.mp4' % fps)

