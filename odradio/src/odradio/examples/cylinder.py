from odradio.simulation.simulation import Simulator
from odradio.viewfactors.viewfactors import Point,Segment
from odradio.extractdata.viennats import Extract
import timeit
import matplotlib.pyplot as plt


def cylidner(diameter,depth,s):

    cyl = Simulator('cylinder')
    #diameter = 100
    #depth  = 4500
    r = diameter/2
    #s = 0.2
    z = depth
    cyl.AddRegion('top',start=Point(0,0 ),end=Point(0,r) ,sticking=1.00,source=1,nelem=1)
    cyl.AddRegion('wall',start=Point(0,r ),end=Point(z,r),sticking=s,source=0,nelem=200)
    cyl.AddRegion('bottom',start=Point(z,r),end=Point(z,0),sticking=1,source=0,nelem=20)
    
    #cone.SetupLSE()
    t = timeit.Timer(lambda: cyl.SetupLSE())
    print(t.timeit(number=1))
    
    size = 30
    ar = depth/diameter
    fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(size,size/ar))
    ax1.plot(cyl.zc,cyl.rc,'o')
    
    #cone.Simulate()
    t = timeit.Timer(lambda: cyl.Simulate())
    print(t.timeit(number=1))
    
    #basepath='/home/manstetten/hptcad_svn/plasma_etching/fitting/raw_data/'
    basepath=''
    fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(8,8))
    a,b,rid = cyl.ExtractResults('wall')
    ax1.semilogy(cyl.zc[a:b]/depth,cyl.absnsrc[a:b],'-b')
    ax1.semilogy(cyl.zc[a:b]/depth,cyl.dirnsrc[a:b],'--r')
    _, _, r, rval, _, _, z, zavg, zmin, zmax = Extract(diameter,depth,s,basepath)
    ax1.semilogy(z[::1],zavg[::1],'--')
    ax1.semilogy(z[::1],zmin[::1],'-')
    ax1.semilogy(z[::1],zmax[::1],'-')
    fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(8,8))
    
    a,b,rid = cyl.ExtractResults('bottom')
    ax1.semilogy(cyl.rc[a:b]/(diameter/2),cyl.absnsrc[a:b],'-ob')
    ax1.semilogy(cyl.rc[a:b]/(diameter/2),cyl.dirnsrc[a:b],'-r')
    _, _, r, rval, _, _, z, zavg, zmin, zmax = Extract(diameter,depth,s,basepath)
    ax1.semilogy(r[::1],rval[::1],'+')
    
    plt.show()
