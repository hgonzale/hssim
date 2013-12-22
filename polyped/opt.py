#!/usr/bin/python

#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public
#  License as published by the Free Software Foundation; either
#  version 3.0 of the License, or (at your option) any later version.
#
#  The library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
# (c) Sam Burden, Berkeley 2012 

import numpy as np
import pylab as plt

import os
from glob import glob

np.set_printoptions(precision=2)

dbg = False

# optimization libraries
try:
  import scipy as sp
except ImportError:
  if dbg:
    print 'WARNING: scipy not found'
  sp = None

try:
  import nlopt
except ImportError:
  if dbg:
    print 'WARNING: nlopt not found'
  nlopt = None

try:
  import openopt as oo
except ImportError:
  if dbg:
    print 'WARNING: openopt not found'
  oo = None

try:
  import optutil as ou
except ImportError:
  if dbg:
    print 'WARNING: optutil not found'
  ou = None

# default environment for evaluating parameter files
env = {'np':np,'array':np.array,
       'scipy':sp,
       'nlopt':nlopt,
       'openopt':oo,
       'optutil':ou,
       '__builtins__':__builtins__}

def bd(x,xm,xM,dbg=True):
  """
  Project x to keep it within bounds xm,xM

  Inputs:
    x - array - state
    xm,xM - array - elementwise min,max on x

  Outputs:
    x_ - array - xm < x < xM

  by Sam Burden 2012
  """
  x_ = np.asarray(x,dtype=np.float).copy(); 
  xm = np.asarray(xm,dtype=np.float); 
  xM = np.asarray(xM,dtype=np.float);
  jm = (x_ < xm).nonzero()[0]
  if jm.size > 0:
    if dbg:
      print 'WARNING: x%s < xm%s' % (list(jm),list(jm))
    x_[jm] = xm[jm]
  jM = (x_ > xM).nonzero()[0]
  if jM.size > 0:
    if dbg:
      print 'WARNING: x%s > xM%s' % (list(jM),list(jM))
    x_[jM] = xM[jM]
  return x_

class Opt(object):
  """
  Optimization superclass

  Usage:
    opt = Opt('pars.py') # pars.py specifies optimization parameters
    while opt.f(opt.x) > 1e-6:
      # update opt.x
      opt.unpack(x)
      print opt
  """
  def __init__( self,pth='' ):
    """
    Opt(pth)  creates optimization instance in directory / file given by pth

    Inputs:
      pth - str - parameter .py file or directory containing only one .py file
    """
    self.p = {}
    if pth and os.path.exists(pth):
      # sanitize filename
      if os.path.isdir(pth):
        fis = glob( os.path.join(pth,'*.py') )
        assert len(fis) == 1 # only one .py in directory ?
        fi = fis[0]
      else:
        fi = pth
      self.di,self.fi = os.path.split(fi)
      self.pars(os.path.join(self.fi))

  def pars( self,fi='',p0=env,**p ):
    """
    pars(fi)  reads parameters from .py file fi
    pars(**p) reads parameters from keyword arguments

    Inputs:
      (optional)
      fi - str - parameter .py file
      p0 - dict - initial parameter environment to use when exec'ing fi 
      **p - keywords - additional parameters specified as keyword arguments 

    Outputs:
      p - dict - current optimization params

    Effects:
      - exec's fi in self.p environment
      - runs self.sanitize on self.p
      - runs self.pack on self.p
    """
    if fi or p:
      self.p.update(p0)
      self.p.update(p)
      if fi:
        code = compile(open(fi).read(),fi,'exec')
        exec code in self.p
      self.p = self.sanitize(self.p,rm=p0.keys())
      if 'opt' in self.p.keys() and self.p['opt']:
        self.vars = self.p['opt']['vars']
        if 'cost' in self.p['opt']:
          self.cost = self.p['opt']['cost']
        self.pack(self.p)
    if hasattr(self, 'x') and hasattr(self, 'vars'):
      return dict([(var,(self.x[self.j[i]:self.j[i+1]])) 
                   for i,var in enumerate(self.vars)])
      
  def sanitize( self,p,rm={} ):
    """
    sanitize(p)  adds and/or resizes _m,_M,_s fields for each var

    Inputs:
      p - parameter dict
      (optional)
      rm - parameters to remove

    Outputs:
      p - sanitized parameter dict
    """
    # remove some keys
    for key in p.keys():
      if key in rm:
        p.pop(key)
    if not 'opt' in p.keys():
      return p
    # sanitize optimization bounds, scale
    for var in p['opt']['vars']:
      # min
      if var+'_m' not in p.keys():
        p[var+'_m'] = -np.inf*np.ones_like(p[var])
      elif not np.asarray(p[var+'_m']).size == np.asarray(p[var]).size:
        p[var+'_m'] = p[var+'_m']*np.ones_like(p[var])
      # max
      if var+'_M' not in p.keys():
        p[var+'_M'] = np.inf*np.ones_like(p[var])
      elif not np.asarray(p[var+'_M']).size == np.asarray(p[var]).size:
        p[var+'_M'] = p[var+'_M']*np.ones_like(p[var])
      # scale
      if var+'_s' not in p.keys():
        # default to unity
        if np.any(np.isinf(np.r_[p[var+'_M'],p[var+'_m']])):
          p[var+'_s'] = 1.*np.ones_like(p[var]) 
        # default to .1*(max - min)
        else:
          p[var+'_s'] = .1*(p[var+'_M']-p[var+'_m'])
      elif not np.asarray(p[var+'_s']).size == np.asarray(p[var]).size:
        p[var+'_s'] = p[var+'_s']*np.ones_like(p[var])
    # sanitized parameters
    return p

  def pack(self, p=None):
    """
    pack(p)  updates optimization variables given parameter dict 

    Inputs:
      p - updated parameter dict

    Outputs:
      x - array - optimization variables
      n - int - x.size
      j - int array - optimization variable indices
      m - array - min, x > m
      M - array - max, x < M
      s - array - scale, x ~ s

    Effects:
      updates self.x,n,j,m,M,s
    """
    if p is None:
      p = self.p
    # unpack optimization params
    a = [np.asarray(p[var]) for var in p['opt']['vars']]
    x = np.hstack(a) # optimization variables
    n = x.size # number of optimization variables
    j = np.hstack( (0,np.cumsum([aa.size for aa in a])) ) # indices
    m = np.hstack([np.asarray(p[var+'_m']) for var in self.vars]) # min
    M = np.hstack([np.asarray(p[var+'_M']) for var in self.vars]) # max
    s = np.hstack([np.asarray(p[var+'_s']) for var in self.vars]) # scale
    self.x = x; self.n = n; self.j = j;
    self.m = m; self.M = M; self.s = s;
    return x,n,j,m,M,s

  def unpack(self, x):
    """
    unpack(x)  updates parameter dict given optimization variables

    Inputs:
      x - array - optimization variables

    Outputs:
      self.p - updated parameter dict

    Effects:
      updates self.p
    """
    q = []
    if np.asarray(x).shape == ():
      x = np.asarray([x])
    for i,var in enumerate(self.p['opt']['vars']):
      if self.j[i+1] == self.j[i]+1:
        q.append((var,x[self.j[i]]))
      else:
        q.append((var,x[self.j[i]:self.j[i+1]]))
    self.p.update(q)
    return self.p

  def cost( self,x,*args ):
    """
    fx = cost(x,*args)  cost function

    Inputs:
      x - array - optimization variables
      (optional)
      args - list - extra parameters

    Outputs:
      fx - scalar - cost function value
    """
    return np.nan

  def __repr__( self ):
    """
    __repr__()  string representation of optimization state
    """
    if hasattr(self,'cost') and hasattr(self,'x'):
      return 'cost(x) = {:0.2f}; x = {!s}'.format( self.cost(self.x), self.x )
    else:
      return 'p = {}'.format(self.p)

class SPOpt(Opt):
  """
  scipy.optimize interface class

  Usage:
    op = SPOpt('pars.py') # pars.py specifies optimization parameters
    op.init()
  """
  dbg = True
  vars_tol = 1e-3
  cost_tol = 1e-3
  solver = None
  res = None
  ftol = 1e-3
  xtol = 1e-3

  def __init__( self,pth='',opt=None ):
    """
    SPOpt(pth)  creates optimization instance

    Inputs:
      pth - str - parameter .py file 
      op  - Opt - object with base class of Opt
    """
    Opt.__init__( self,pth )
    if opt is not None and isinstance(opt,Opt):
      self.pars(**opt.p)
    if 'opt' in self.p:
      self.vars_tol = self.p['opt'].get('vars_tol',self.vars_tol)
      self.cost_tol = self.p['opt'].get('cost_tol',self.cost_tol)
      self.solver   = self.p['opt'].get('solver',self.solver)
      self.ftol     = self.p['opt'].get('ftol',self.ftol)
      self.xtol     = self.p['opt'].get('xtol',self.xtol)

  def cost( self,x,*args ):
    """
    fx = cost(x,*args)  cost function

    Inputs:
      x - N-array - optimization variables
      (optional)
      args - list - extra parameters

    Outputs:
      fx - M-array - cost function values
    """
    return [np.nan]


  def solve( self ):
    """
    solve()  runs self.solver on optimization problem

    Effects:
      - assigns self.res
    """
    bounds = []
    for m,M in zip(self.m,self.M):
      bd = [m,M]
      if np.isinf(m):
        bd[0] = None
      if np.isinf(M):
        bd[1] = None
      bounds.append(bd)

    diag = self.s**-1

    if self.solver == 'leastsq':
      self.res = sp.optimize.leastsq(self.cost,self.x,full_output=True,
                                     xtol=self.xtol,ftol=self.ftol,
                                     diag=diag)
    #res = sp.optimize.fmin_l_bfgs_b(lambda opx : (np.sum(err(opx*op.s)[0]**2) / eta0.size),
    #                                opx1*op.s**-1,
    #                                pgtol=1e-3,bounds=bounds,
    #                                approx_grad=True,epsilon=1e-5)
    #res = sp.optimize.fmin_slsqp(lambda opx : (np.sum(err(opx*op.s)[0]**2) / eta0.size),
    #                             opx1*op.s**-1,
    #                             acc=1e-3,bounds=bounds,
    #                             epsilon=1e-5)
    return self.res

  def __repr__( self ):
    """
    __repr__()  string representation of optimization state
    """
    cost = self.cost(self.x)
    if len(cost) > 1:
      cost = np.sum(cost**2)
    if hasattr(self,'cost') and hasattr(self,'x'):
      return 'cost(x) = {:0.2f}; {!s}'.format( cost, self.pars() )

    else:
      return 'p = {}'.format(self.p)

class NM(Opt):
  """
  Nelder-Mead optimization class

  @article{NelderMead1965,
           Author = {Nelder, J. A. and Mead, R.},
           Journal = {The Computer Journal},
           Number = {4},
           Pages = {308-313},
           Title = {A Simplex Method for Function Minimization},
           Volume = {7},
           Year = {1965}}

  Usage:
    nm = NM('pars.py') # pars.py specifies optimization parameters
    nm.init()
    while nm.step():
      print nm
  """
  dbg = True
  vol = 1e-10
  vars_tol = 1e-3
  cost_tol = 1e-3

  def __init__( self,pth='' ):
    """
    NM(pth)  creates optimization instance in directory / file given by pth

    Inputs:
      pth - str - parameter .py file or directory containing only one .py file
    """
    Opt.__init__( self,pth )
    if 'vars_tol' in self.p['opt'].keys():
      self.vars_tol = self.p['opt']['vars_tol']
    if 'cost_tol' in self.p['opt'].keys():
      self.cost_tol = self.p['opt']['cost_tol']
    if 'vol' in self.p['opt'].keys():
      self.vol = self.p['opt']['vol']

  def init( self,x0=None,fx0=None ):
    """
    init(x0)  initializes Nelder-Mead optimization with simplex x0

    NOTE: if X.csv exists in self.di, simplex will be loaded

    Inputs:
      (optional)
      x0 - (n+1) x n - starting simplex
      fx0 - (n+1) x 1 - cost at starting simplex
    """
    self.k = 0
    x = self.x; s = self.s; p = self.p; n = self.n; di = self.di
    vars_tol = self.vars_tol; cost_tol = self.cost_tol
    if x0 is None:
      # initialize simplex
      x0 = np.vstack(( x, x + np.diag(s*(.5+np.random.rand(n))) ))
    assert x0.shape[1] == n # simplex consistent with parameter file ?
    assert x0.shape[0] >= n + 1 # full starting simplex ?
    if fx0 is not None:
      assert x0.shape[0] == fx0.shape[0] # cost evaluated on whole simplex ?
    X = []; F = []
    # load prior simplex
    pth = os.path.join(di,'X.csv')
    if os.path.exists(pth):
      if self.dbg:
        print 'loading {:s}'.format(pth)
      X = [list(a) for a in np.loadtxt(pth,delimiter=',')]
      assert a.size == n # X.csv compatible with *.py ?
      assert len(X) >= n + 1 # X.csv contains full starting simplex ?
      x0 = X[-(n+1):]
      F = [np.nan for _ in range(len(X))]
      # load prior cost
      pth = os.path.join(di,'F.csv')
      if os.path.exists(pth):
        F = list(np.loadtxt(pth,delimiter=','))
        assert len(X) == len(F) # X.csv compatible with F.csv ?
      fx0 = F[-(n+1):]
      #m = min(len(X),n+1)
      #x0 = np.vstack(( X[-m:], x + s*np.random.rand(n+1-m,n) ))
      #f0 = np.hstack(( F[-m:], np.nan*np.ones(n+1-m) ))
    self.X = X; self.F = F
    #print np.asarray(x0)
    #print np.asarray(fx0)
    self.op = ou.fminIter(x0,f0=fx0,xtol=vars_tol,ftol=cost_tol)
    #1/0
    self.fx = self.op.next()

  def step( self ):
    """
    step()  executes one step of Nelder-Mead optimization

    Effects:
      - saves X.csv,F.csv
    """
    # extract next element from NM iterator
    try:
      self.x = self.op.next()
    except StopIteration:
      return False
    self.k += 1
    # keep params within bounds
    self.x_ = bd( self.x, self.m, self.M, dbg=self.dbg )
    # evaluate cost function; penalize for leaving specified bounds
    fx = np.asarray([self.cost(**self.pars(x=self.x_))]).flatten()
    assert fx.size == 1 # cost function returns scalar ?
    self.fx[-1] = fx[0]
    self.fx[-1] += np.exp(np.sum(np.abs(self.x_ - self.x) / self.s)) - 1.
    # store history
    self.X.append(list(self.x))
    self.F.append(self.fx[-1])
    # check that simplex has volume
    n = self.n
    if len(self.X) >= n:
      _,s,_ = np.linalg.svd(self.X[-n:])
      r = np.sum(s > self.vol)
      if r < n and self.dbg:
        print 'WARNING: simplex has rank %d < %d = dim x' % (r,n)
    return True

  def __repr__( self ):
    """
    __repr__()  string representation of optimization state
    """
    return 'k = {:4d}; f = {:0.2e}; p = {:s}'.format( 
            self.k, self.fx[-1], self.pars())

  def save( self ):
    """
    save()  save nm progress

    Effects:
      - saves X.csv and F.csv in same directory as self.fi
    """
    np.savetxt(os.path.join(self.di,'X.csv'),self.X,delimiter=',')
    np.savetxt(os.path.join(self.di,'F.csv'),self.F,delimiter=',')

if __name__ == "__main__":

  pth = 'pars.py'

  import sys
  if len(sys.argv) > 1:
    pth = sys.argv[1]

  # random initial simplex
  nm = NM(pth)
  nm.init()
  while nm.step() and nm.k <= 100:
    print nm
  nm.save()
  print np.asarray(nm.X[-(nm.n+1):])
  print np.asarray(nm.F[-(nm.n+1):])

  # load initial simplex
  nm = NM(pth)
  nm.vars_tol *= .01
  nm.init()
  print np.asarray(nm.X[-(nm.n+1):])
  print np.asarray(nm.F[-(nm.n+1):])
  while nm.step() and nm.k <= 100:
    print nm

