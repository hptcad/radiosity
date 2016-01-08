from odradio.simulation.simulation import Simulator
from odradio.viewfactors.viewfactors import Point
from odradio.extractdata.viennats import Extract
import timeit
import matplotlib.pyplot as plt

def cone(diameter,depth,s,rfactor=1):

    cone = Simulator('cone')
    #diameter = 100
    #depth  = 4500
    r = diameter/2
    #s = 0.2
    z = depth
    #rfactor = 1
    cone.AddRegion('top',start=Point(0,0 ),end=Point(0,r) ,sticking=1.00,source=1,nelem=4)
    cone.AddRegion('wall',start=Point(0,r ),end=Point(z,r*rfactor),sticking=s,source=0,nelem=200)
    cone.AddRegion('bottom',start=Point(z,r*rfactor),end=Point(z,0),sticking=1,source=0,nelem=20)
    
    #cone.SetupLSE()
    t = timeit.Timer(lambda: cone.SetupLSE())
    print(t.timeit(number=1))
    
    size = 30
    ar = depth/diameter
    fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(size,size/ar))
    ax1.plot(cone.zc,cone.rc,'o')
    
    #cone.Simulate()
    t = timeit.Timer(lambda: cone.Simulate())
    print(t.timeit(number=1))
    
    #basepath='/home/manstetten/hptcad_svn/plasma_etching/fitting/raw_data/'
    basepath=''
    fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(8,8))
    a,b,rid = cone.ExtractResults('wall')
    ax1.semilogy(cone.zc[a:b]/depth,cone.absnsrc[a:b],'-b')
    ax1.semilogy(cone.zc[a:b]/depth,cone.dirnsrc[a:b],'--r')
    _, _, r, rval, _, _, z, zavg, zmin, zmax = Extract(diameter,depth,s,basepath)
    ax1.semilogy(z[::1],zavg[::1],'--')
    ax1.semilogy(z[::1],zmin[::1],'-')
    ax1.semilogy(z[::1],zmax[::1],'-')
    
    fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(8,8))
    a,b,rid = cone.ExtractResults('bottom')
    ax1.semilogy(cone.rc[a:b]/(diameter/2),cone.absnsrc[a:b],'-ob')
    ax1.semilogy(cone.rc[a:b]/(diameter/2),cone.dirnsrc[a:b],'-r')
    _, _, r, rval, _, _, z, zavg, zmin, zmax = Extract(diameter,depth,s,basepath)
    ax1.semilogy(r[::1],rval[::1],'+')
    
    plt.show()