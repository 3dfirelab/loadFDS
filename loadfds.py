import numpy as np 
import matplotlib.pyplot as plt
import fdsreader
from scipy.interpolate import RegularGridInterpolator, griddata
import sys
import pdb 



###########################
def load_slice(sim, id_var):
    slc = sim.slices[id_var] # to load temperature. see order of &SLCF in fds config file 
    data, grid = slc.to_global(return_coordinates=True, masked=True)  # need version 1.9.13

    times = slc.times

    return data, grid, times



###########################
def load_3dsmoke(sim, id_var, dx, dy, dz):
    
    data = sim.smoke_3d[id_var] #0, for smoke, 1 for HRRPUV

    datapermesh = list(data.subsmokes)
    nmesh = len(datapermesh)

    grids = []
    data = []

    xmin = 1.e6; xmax = -1.e6
    ymin = 1.e6; ymax = -1.e6
    zmin = 1.e6; zmax = -1.e6

    print('loading data ...')
    for imesh in range(nmesh):
        data.append(datapermesh[imesh][1].data)
        grids.append(datapermesh[imesh][1].mesh)

        xmin = min([xmin, grids[-1].coordinates['x'].min()])
        xmax = max([xmax, grids[-1].coordinates['x'].max()])
        ymin = min([ymin, grids[-1].coordinates['y'].min()])
        ymax = max([ymax, grids[-1].coordinates['y'].max()])
        zmin = min([zmin, grids[-1].coordinates['z'].min()])
        zmax = max([zmax, grids[-1].coordinates['z'].max()])

        if imesh == 0: 
            times = datapermesh[imesh][1].times
        
    xg = np.round(np.arange(xmin,xmax,dx), 4)
    yg = np.round(np.arange(ymin,ymax,dy), 4)
    zg = np.round(np.arange(zmin,zmax,dz), 4)

    ngx = xg.shape[0]
    ngy = yg.shape[0]
    ngz = zg.shape[0]

    ntime = times.shape[0]

    dataf = np.zeros([ntime,ngx,ngy,ngz]) - 999

    yyg, xxg, zzg = np.meshgrid(yg,xg,zg)
    jjg, iig, kkg = np.meshgrid(np.arange(ngy),np.arange(ngx),np.arange(ngz))

    print('running interpolation ...')
    for it in range(ntime)[::100]: 
        print('{:4.1f} % '.format(100.*it/ntime), end='\r')
        sys.stdout.flush()
        interps = []
        
        if it == 0 : 
            x,y,z = [], [], []
            xx,yy,zz = [], [], []
            bounds = []
            for im in range(nmesh):
                x.append(np.round(grids[im].coordinates['x'][:],4)[1:-1])
                y.append(np.round(grids[im].coordinates['y'][:],4)[1:-1])
                z.append(np.round(grids[im].coordinates['z'][:],4)[1:-1])
                
                x[-1][0] -= 0.001
                x[-1][-1] += 0.001
                y[-1][0] -= 0.001
                y[-1][-1] += 0.001
                z[-1][0] -= 0.001
                z[-1][-1] += 0.001

                bounds.append([x[-1].min(),x[-1].max(),y[-1].min(),y[-1].max(),z[-1].min(),z[-1].max()])

                yy_, xx_, zz_ = np.meshgrid(y[-1],x[-1],z[-1])
                
                xx.append(xx_)
                yy.append(yy_)
                zz.append(zz_)
        
        points = []; values = []
        for im in range(nmesh):

            [points.append( (x_,y_,z_) ) for x_,y_,z_ in zip(xx[im].flatten(),yy[im].flatten(),zz[im].flatten())]
            [values.append( data_) for data_ in data[im][it,1:-1,1:-1,1:-1].flatten() ]

            
        dataf[it,:,:,:] = griddata(points, values, (xxg,yyg,zzg), method='nearest') # any other method is very slow in 3D

    print('done          ')
    return dataf, (xg,yg,zg), times




###############
if __name__ == '__main__':
###############

    if False: 
        '''
        example for loading temperature from simulation NIST_PoolFire
        center vertical slice is shown at last time, t=2.5s
        '''
        sim = fdsreader.Simulation('./data/NIST_PoolFire/') #load simulation
        dx = 1.e-2
        dy = 1.e-2
        dz = 1.e-2

        temp, grid, times = load_slice(sim, 3)
        
        ax = plt.subplot(111)
        plt.imshow(temp[-1,:,30,:].T, origin='lower')
        plt.title('center vertical slice at t = {:.1f}s '.format(times[-1]))

        plt.show()

    if False: 
        '''
        example for loading 3d smoke from simulation NIST_PoolFire
        center vertical slice is shown at last time
        '''    
        sim = fdsreader.Simulation('./data/NIST_PoolFire/') #load simulation
        dx = 1.e-2  # also tested to dx=dy=dx=2cm
        dy = 1.e-2
        dz = 1.e-2

        smoke, grid, times = load_3dsmoke(sim, 0, dx, dy, dz)

        plt.imshow(smoke[800,:,int(smoke.shape[1]/2),:].T, origin='lower')
        plt.show()
        sys.exit()
    
    if True: 
        '''
        example for loading 3d smoke from simulation Entrepinos HF LC 7.6.5
        center vertical slice is shown at last time
        '''    
        sim = fdsreader.Simulation('./data/Entrepinos HF LC 7.6.5/') #load simulation
        dx = 20.e-2  
        dy = 20.e-2
        dz = 20.e-2

        smoke, grid, times = load_3dsmoke(sim, 0, dx, dy, dz)

        plt.imshow(smoke[200,:,int(smoke.shape[2]/2),:].T, origin='lower')
        plt.show()
        sys.exit()
