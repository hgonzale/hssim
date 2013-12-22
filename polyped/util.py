import copy

class Lambda(object):
    """
    Lambda  lambda-like class

    Acts like a lambda function, but its string representation
    is Python code that yields the object when executed.

    >>> f = Lambda('x : x**2')
    >>> f(1.41)
    1.9880999999999998
    >>> g = eval(str(f))
    >>> g(2.82)
    7.952399999999999
    >>> f = Lambda('lambda x,y : x + y')
    >>> f(2,-0.1)
    1.8999999999999999
    """
    def __init__(self, lam):
        if not lam.strip().startswith('lambda'):
            lam = 'lambda '+lam
        self.lam = lam
        self.op  = eval(lam)
    def __call__(self, *args):
        return self.op(*args)
    def __repr__(self):
        return 'Lambda("'+self.lam+'")'

class Struct(object):
  """ 
  Struct  struct-like class

  Acts like a dict, but keys are members of the object.

  >>> a = Struct(foo='bar', goo=7)
  >>> b = a.copy()
  >>> b.hoo = [1,2,3]
  >>> print a
  {'goo': 7, 'foo': 'bar'}
  >>> print b
  {'hoo': [1, 2, 3], 'goo': 7, 'foo': 'bar'}
  >>> c = eval(str(a))
  >>> a.ioo = []
  >>> print a
  {'ioo': [], 'goo': 7, 'foo': 'bar'}
  >>> print c
  {'goo': 7, 'foo': 'bar'}
  """
  def __init__( self, file=None, **kwds ):
    if file:
      self.__dict__.update( eval( open(file).read() ) )
    self.__dict__.update( kwds )
  def copy( self ):
    return Struct(**copy.deepcopy(self.__dict__))
  def read( self, file, globals={"__builtins__":None}, locals={} ):
    #s = open(file).read()
    #for k in subs.keys():
    #  s = s.replace( k, subs[k] )
    self.__dict__.update( eval( open(file).read(), globals, locals ) )
  def write( self, file ):
    open(file,'w').write(str(self))
  def __repr__( self ):
    return str( self.__dict__ )
  def __getstate__(self):
    return self.__dict__
  def __setstate__( self, kwds ):
    self.__dict__.update( kwds )

if __name__ == "__main__":
    import doctest
    doctest.testmod()
