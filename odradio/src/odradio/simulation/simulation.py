import odradio.viewfactors.viewfactors as vf
from odradio.viewfactors.viewfactors import Point,Segment
import numpy as np
from numpy import linalg as LA
from collections import namedtuple

eps = np.finfo(float).eps

Regiondata = namedtuple('Regiondata', ['name', 'start', 'end','area', 'source', 'issource'], verbose=False)

def jacobi(A,b,Nmax,x=None):
    if x is None:
        x = np.zeros(len(A[0]))
    D = np.diag(A)
    R = A - np.diagflat(D)
    nb = LA.norm(b)
    for i in range(Nmax):
        res = b - np.dot(A,x)
        nres = LA.norm(res)
        if nres/nb < 1e-13:
            return x, i
        x = (b - np.dot(R,x)) / D
    return x, Nmax

class Simulator:
    def __init__(self, name):
        self.name = name
        self.regiondata = []
        self.region = []
        self.source = []
        self.stick = []
        self.n = None
        self.F = None
        self.rho = None
        self.src = None
        self.s = None
        self.RHO = None
        self.RHOF_T = None
        self.A = None
        self.rhs = None
        self.area = None
        self.abs = None
        self.absn = None
        self.absnsrc = None
        self.rc = None
        self.zc = None
        self.dirn = None
        self.dirnsrc = None
    def AddRegion(self, name, start, end, sticking, source, nelem):
        # create Segment
        areatotal = vf.Area(Segment(start,end))
        self.regiondata.append(Regiondata(name=name,start=start,end=end,area=areatotal,source=source,issource=True if source>0 else False))
        r = np.linspace(start.r,end.r,nelem+1)
        z = np.linspace(start.z,end.z,nelem+1)
        reg = []
        for i in range(nelem):
            reg.append(Segment(Point(z[i],r[i]),Point(z[i+1],r[i+1])) )
        src = []
        for s in reg:
            area = vf.Area(s)
            src.append(area/areatotal*source)
        self.source.append(src)
        self.region.append(reg)            
        self.stick.append(list(np.ones(nelem)*sticking))
    def SetupLSE(self):
        # problem size
        n = sum(len(r) for r in self.region)
        self.n = n
        # compute F
        F = np.zeros(shape=(n,n)) 
        area = np.zeros(shape=(n)) 
        rc = np.zeros(shape=(n)) 
        zc = np.zeros(shape=(n)) 
        for ri, ra in enumerate(self.region):
            ioff = sum(len(r) for r in self.region[0:ri])
            for i, _ in enumerate(ra):
                area[ioff+i] = vf.Area(ra[i])
                rc[ioff+i] = vf.rcenter(ra[i])
                zc[ioff+i] = vf.zcenter(ra[i])
                for rj, rb in enumerate(self.region):
                    joff = sum(len(r) for r in self.region[0:rj]) 
                    for j, _ in enumerate(rb):
                        if joff+j>ioff+i:
                            continue
                        else:
                            F[ioff+i,joff+j],F[joff+j,ioff+i] = vf.ViewFactor(ra[i],rb[j])
        self.F = F
        self.area = area
        self.rc = rc
        self.zc = zc
        # compute source 
        src = []
        for s in self.source:
            for e in s:
                src.append(e)
        self.src = np.array(src)
        # compute stick 
        s = []
        for st in self.stick:
            for e in st:
                s.append(e)
        self.s = np.array(s)
        self.rho = (1-self.s)
        # compute A
        self.RHO = np.diag(self.rho)
        I = np.identity(n)
        self.RHOF_T = np.dot(self.RHO,F).T
        self.A = np.asarray(I-self.RHOF_T)
        # compute RHS
        self.rhs = np.dot(F.T,src)
        # check if geometry is closed\
        if np.abs(np.sum(np.sum(F,axis=1))-n)>n*1e-8:
            raise ValueError('Geometry is not closed')
        
    def Simulate(self):
        #self.rec = LA.solve(self.A, self.rhs)
        #nmax = 5
        self.rec, nmax = jacobi(self.A,self.rhs,1000)
        print(nmax)
        self.abs = self.rec * self.s
        if abs(sum(self.abs)-sum(self.src))>1e-6*sum(self.src):
            raise ValueError('absorbed enegery does not match source energy')
        self.absn = self.abs/self.area
        srcregion = next(r for r in self.regiondata if r.issource)
        srcnorm = srcregion.source/srcregion.area
        self.absnsrc = self.absn/(srcnorm * self.s)
        
        self.dirn = self.rhs/self.area
        self.dirnsrc = (self.dirn * self.s)/(srcnorm * self.s)
        return nmax
    def ExtractResults(self,name):
        ri = [i for i,x in enumerate(self.regiondata) if x.name == name]
        if len(ri)>1:
            raise ValueError('more than one region with than name')
        st = sum(len(r) for r in self.region[0:ri[0]])
        en = st + len(self.region[ri[0]])
        return st, en, ri[0]