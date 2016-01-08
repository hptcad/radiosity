from .trace  import TraceCalls
import numpy as np
from collections import namedtuple
from numpy import pi, min, abs
eps = np.finfo(float).eps

Point = namedtuple('Point', ['z', 'r'], verbose=False)
Segment = namedtuple('Segment', ['a', 'b'], verbose=False)       

def Distance(one,two):
    return ((one.r-two.r)**2 +  (one.z-two.z)**2 )**0.5

zmin = lambda elem: min([elem.a.z,elem.b.z])
zmax = lambda elem: max([elem.a.z,elem.b.z])
rmin = lambda elem: min([elem.a.r,elem.b.r])
rmax = lambda elem: max([elem.a.r,elem.b.r])
length = lambda elem: Distance(elem.a,elem.b)
zlen = lambda elem: abs(elem.a.z-elem.b.z)
rlen = lambda elem: abs(elem.a.r-elem.b.r)
zcenter = lambda elem: 0.5*(elem.a.z+elem.b.z)
rcenter = lambda elem: 0.5*(elem.a.r+elem.b.r)

def Identical(one,two):
    """ test if elements are identical """
    a = np.array([one.a.z,one.b.z,one.a.r,one.b.r])
    b = np.array([two.a.z,two.b.z,two.a.r,two.b.r])
    if abs(sum(a-b))<eps:
        return True
    else:
        return False
def Overlapping(a,b):
    """ test if elements are overlapping """
    onright = bool(zmax(a) <= zmin(b))
    onleft = bool(zmin(a) >= zmax(b))    
    if not onright and not onleft:
        return True
    else:
        return False     
def Ring(segment):
    """ check if segment is a ring """
    if zlen(segment)<eps:
        return True
    else:
        return False
@TraceCalls()
def DiskVisibility(disk,elem):
    """ check if disk is completely visible from element """
    dz = (elem.a.z-elem.b.z) 
    dr = (elem.a.r-elem.b.r) 
    tanalpha = dr/dz
    if elem.a.r-tanalpha*(elem.a.z-disk.z)<disk.r-1e-11:
        raise ValueError('element not completely visible'+str(elem.a.r-tanalpha*(elem.a.z-disk.z)-disk.r))
def ElementVisibility(a,b):
    """ check if elements are completely visible (mutual)"""
    DiskVisibility(a.a,b)
    DiskVisibility(a.b,b)
    DiskVisibility(b.a,a)
    DiskVisibility(b.b,a)
    
def DiskArea(disk):
    return pi*disk.r**2
def Area(segment):
    """ area of a segment """
    a = segment.a
    b = segment.b
    dr = a.r-b.r
    dz = a.z-b.z
    s = (dr**2 + dz**2)**0.5
    return pi*(a.r+b.r)*s        

@TraceCalls()
def Radius2Radius(r1,r2):
    a1 = pi*r1**2
    a2 = pi*r2**2
    f12 = a2/a1 if r1>r2 else 1
    f21 = a1/a2 if r2>r1 else 1
    return f12, f21
@TraceCalls()
def Point2Disk(z,r):      
    sin_alpha = r/(z**2+r**2)**0.5
    f12 = (sin_alpha)**2
    f21 = 0
    a1 = 0
    a2 = pi*r**2
    return f12, f21
@TraceCalls()
def Disk2Point(r,z):
    a,b = Point2Disk(z,r)
    return b, a
@TraceCalls()
def Disk2Disk(a,b): 
    z=abs(a.z-b.z) 
    r1 = a.r
    r2 = b.r
    if r1<eps and r2>eps:
        return Point2Disk(z,r2)
    elif r2<eps and r1>eps:
        return Disk2Point(r1,z)
    elif r1<eps and r2<eps:
        return 0,0
    elif z<eps and r1>eps and r2>eps:
        return Radius2Radius(r1,r2)
    else:
        def R(r,z):
            return r/z   
        def X(R1,R2):
            return 1+(1+R2**2)/R1**2  
        def F(X,R1,R2):
            return  0.5 * (X - (X**2-4*(R2/R1)**2 )**0.5 )
        R1 = R(r1,z)
        R2 = R(r2,z)
        X = X(R1,R2)
        a2b = F(X,R1,R2)  
        b2a = DiskArea(a)/DiskArea(b) * a2b
        return a2b, b2a
@TraceCalls()
def Disk2Ring(disk,ring): 
    d2ra,_ = Disk2Disk(disk,ring.a)
    d2rb,_ = Disk2Disk(disk,ring.b)    
    d2r = abs(d2ra-d2rb)
    r2d = DiskArea(disk)/Area(ring) * d2r 
    return d2r, r2d
@TraceCalls()
def Ring2Disk(ring,disk):
    a, b = Disk2Ring(disk,ring)
    return b,a
@TraceCalls()
def Disk2Element(disk,elem): 
    DiskVisibility(disk,elem)
    d2ea,_ = Disk2Disk(disk,elem.a)
    d2eb,_ = Disk2Disk(disk,elem.b)
    d2e = abs(d2ea-d2eb)
    e2d = DiskArea(disk)/Area(elem) * d2e
    return d2e, e2d
@TraceCalls()
def Element2Disk(elem,disk):
    d2e, e2d = Disk2Element(disk,elem)
    return e2d, d2e
@TraceCalls()
def Ring2Ring(one,two):
    one2a,_ = Ring2Disk(one,two.a)
    one2b,_ = Ring2Disk(one,two.b)
    one2two = abs(one2a-one2b)
    two2one = Area(one)/Area(two) * one2two
    return one2two, two2one
@TraceCalls()
def Ring2Element(ring,elem):
    e2ra,_ = Element2Disk(elem,ring.a)
    e2rb,_ = Element2Disk(elem,ring.b)
    e2r = abs(e2ra-e2rb)
    r2e = Area(elem)/Area(ring) * e2r
    return r2e, e2r
@TraceCalls()
def Element2Ring(elem,ring):
    r2e, e2r = Ring2Element(ring,elem)
    return e2r, r2e
@TraceCalls()
def Element2Self(elem):
    aside,_ = Element2Disk(elem,elem.a)
    bside,_ = Element2Disk(elem,elem.b)
    self = 1-aside-bside
    return self,self
@TraceCalls()
def Ring2Self(ring):
    return 0, 0
@TraceCalls()
def Element2Element(one,two):
    ElementVisibility(one,two) 
    a2cd,cd2a = Disk2Element(one.a,two)
    b2cd,cd2b = Disk2Element(one.b,two)  
    cd2ab = abs(cd2b-cd2a)
    ab2cd = Area(two)/Area(one) * cd2ab
    return ab2cd, cd2ab
@TraceCalls() 
def ViewFactor(one,two):
    """ view factor between two segments"""
    if length(one)<eps or length(two)<eps:
        raise ValueError('Element with zero length')      
    if Overlapping(one,two):
        if Identical(one,two):
            if Ring(one):
                return Ring2Self(one)
            else:
                return Element2Self(one)  
        else:
            raise ValueError('Elements overlap')         
    elif Ring(one) and Ring(two) and abs(zmin(one)-zmin(two))<eps:
        return 0,0       
    elif Ring(one) and Ring(two) and abs(zmin(one)-zmin(two))>eps:
        return Ring2Ring(one,two)      
    elif Ring(one) and not Ring(two):
        return Ring2Element(one,two)
    elif Ring(two) and not Ring(one):
        return Element2Ring(one,two)  
    else:
        return Element2Element(one,two)
    