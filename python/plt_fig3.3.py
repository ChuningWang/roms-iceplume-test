import numpy as np
from scipy.io import loadmat
import pandas as pd
import netCDF4 as nc
# import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from cmocean import cm
import pyroms

# matplotlib.rcParams.update({'font.size': 10})

grid_file = '/home/chuning/projects/fjord/roms_files/roms_grd.nc'
vel_file = '/home/chuning/data/fjord_ks/KS_sect1_14transects_2014_v3.mat'
ctd_file = '/home/chuning/data/fjord_ks/KS2014_KS1_TS.mat'
# out_dir = '/mnt/c/Users/chuning/Documents/JBox/ms1detrainment/'
out_dir   = '/mnt/c/Users/chuning/Desktop/zview_uv/'

rcase1 = 'hmass'
rcase2 = 'vmass'
rcase3 = 'mix'

rcase = [rcase1, rcase2]

roms_file = ['/home/chuning/projects/fjord/data_mpdata/roms_avg_' + i + '.nc'
             for i in rcase]

t0 = 100
z0 = 50
x0 = 3
tni = 6
ang = 25

# ----------------------------------------------------------------------

vel_data = loadmat(vel_file)
tn = vel_data['timeN'].squeeze()
xn = vel_data['xN'].squeeze()
zn = vel_data['zN'].squeeze()*(-1)
un = -1*vel_data['u_rot_N']
vn = -1*vel_data['v_rot_N']

dtn = []
for i, ti in enumerate(tn):
    dtni = pd.to_datetime(ti-719529, unit='D')
    dtn.append(dtni.round('S'))

xn = xn - 2.65
xvel = np.ones(len(xn))*3
yvel = xn

un = np.ma.masked_invalid(un)
vn = np.ma.masked_invalid(vn)
un = un[:, ::-1]
vn = vn[:, ::-1]

ctd_data = loadmat(ctd_file)
sc = ctd_data['july2829'][0, 0][0]
tc = ctd_data['july2829'][0, 0][1]
zc = -1*ctd_data['july2829'][0, 0][2].squeeze()

sc = np.ma.masked_invalid(sc)
tc = np.ma.masked_invalid(tc)
sc = np.mean(sc, axis=1)
tc = np.mean(tc, axis=1)

# ----------------------------------------------------------------------

fh = nc.Dataset(grid_file)
xmask = fh.variables['mask_rho'][0, :]
ymask = fh.variables['mask_rho'][:, 2]
x = fh.variables['x_rho'][0, :]
fh.close()

idx = np.where(ymask == 1)[0]
eta1 = idx.min()
eta2 = idx.max()+1
idx = np.where((x > 0) & (xmask == 0))[0]
xi1 = idx.min()
xi2 = idx.max()+1

xi2 = xi2+5

grd = pyroms.grid.get_ROMS_grid(grid_file)
xr = grd.hgrid.x_rho[0, xi1:xi2]/1000
yr = grd.hgrid.y_rho[eta1:eta2, 0]/1000
xv = grd.hgrid.x_vert[0, xi1:xi2+1]/1000
yv = grd.hgrid.y_vert[eta1:eta2+1, 0]/1000
xp = grd.hgrid.x_psi[0, xi1:xi2-1]/1000
yp = grd.hgrid.y_psi[eta1:eta2-1, 0]/1000
zr = grd.vgrid.z_r[:][:, eta1, xi1]
zv = grd.vgrid.z_w[:][:, eta1, xi1]
h = grd.vgrid.h[eta1:eta2, xi1:xi2]

zi = np.argmin(np.abs(zr+abs(z0)))
xi = np.argmin(np.abs(xr-x0))

fh = nc.Dataset(roms_file[0])
time = fh.variables['ocean_time'][t0:]
temp = fh.variables['temp'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
salt = fh.variables['salt'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
uu = fh.variables['u'][t0:, :, eta1:eta2, xi1-1:xi2].mean(axis=0)
vv = fh.variables['v'][t0:, :, eta1-1:eta2, xi1:xi2].mean(axis=0)
vve = fh.variables['v'][t0:, :, eta1-1:eta2, xi1-1:xi2+1].mean(axis=0)
fh.close()

uu[uu.mask] = 0
vv[vv.mask] = 0
vve[vve.mask] = 0

ur = 0.5*(uu[:, :, 1:] + uu[:, :, :-1])
vr = 0.5*(vv[:, 1:, :] + vv[:, :-1, :])
vu = 0.25*(vve[:, 1:, 1:]+vve[:, 1:, :-1]+vve[:, :-1, :1]+vve[:, :-1, :-1])

fh = nc.Dataset(roms_file[1])
temp2 = fh.variables['temp'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
salt2 = fh.variables['salt'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
uu2 = fh.variables['u'][t0:, :, eta1:eta2, xi1-1:xi2].mean(axis=0)
vv2 = fh.variables['v'][t0:, :, eta1-1:eta2, xi1:xi2].mean(axis=0)
vve2 = fh.variables['v'][t0:, :, eta1-1:eta2, xi1-1:xi2+1].mean(axis=0)
fh.close()

uu2[uu2.mask] = 0
vv2[vv2.mask] = 0
vve2[vve2.mask] = 0

ur2 = 0.5*(uu2[:, :, 1:] + uu2[:, :, :-1])
vr2 = 0.5*(vv2[:, 1:, :] + vv2[:, :-1, :])
vu2 = 0.25*(vve2[:, 1:, 1:]+vve2[:, 1:, :-1] +
            vve2[:, :-1, :1]+vve2[:, :-1, :-1])

# ----------------------------------------------------------------------

fig = plt.figure()
gs = GridSpec(100, 100, figure=fig,
              left=0.06, right=0.95, bottom=0.06, top=0.98)

vel = ur[:, :, xi] + 1j*vr[:, :, xi]
theta = np.arctan(np.imag(vel.max())/np.real(vel.max()))
vel = vel*np.exp(-1j*theta)

axc = fig.add_subplot(gs[42:, :30])
axc.set_xlim(-2, 2)
axc.set_ylim(-150, 0)
axc.set_xticks([-2, -1, 0, 1, 2])
axc.set_yticks([-150, -100, -50, 0])
axc.set_yticklabels(['150', '100', '50', '0'])
axc.set_ylabel('Z [m]')
axc.set_xlabel('Y [km]')
axc.tick_params(axis='both', direction='in')

ctf = axc.contourf(yr, zr, np.real(vel), np.linspace(-0.2, 0.2, 41),
                   extend='both', cmap=cm.balance)
axc.contour(xn, zn, un[:, :, tni], [0], colors='k', linewidths=1.5)
axc.contour(xn, zn, un[:, :, tni], np.linspace(-0.3, 0.3, 13),
            colors='k', linewidths=0.5)

vel2 = ur2[:, :, xi] + 1j*vr2[:, :, xi]
theta2 = np.arctan(np.imag(vel2.max())/np.real(vel2.max()))
vel2 = vel2*np.exp(-1j*theta2)

axc2 = fig.add_subplot(gs[42:, 32:62])
axc2.set_xlim(-2, 2)
axc2.set_ylim(-150, 0)
axc2.set_xticks([-2, -1, 0, 1, 2])
axc2.set_yticks([-150, -100, -50, 0])
axc2.set_yticklabels([])
axc2.set_xlabel('Y [km]')
axc2.tick_params(axis='both', direction='in')

ctf = axc2.contourf(yr, zr, np.real(vel2), np.linspace(-0.2, 0.2, 41),
                    extend='both', cmap=cm.balance)
axc2.contour(xn, zn, un[:, :, tni], [0], colors='k', linewidths=1.5)
axc2.contour(xn, zn, un[:, :, tni], np.linspace(-0.3, 0.3, 13),
             colors='k', linewidths=0.5)

axc3 = fig.add_subplot(gs[42:, 64:94])
axc3.set_xlim(-2, 2)
axc3.set_ylim(-150, 0)
axc3.set_xticks([-2, -1, 0, 1, 2])
axc3.set_yticks([-150, -100, -50, 0])
axc3.set_yticklabels([])
axc3.set_xlabel('Y [km]')
axc3.tick_params(axis='both', direction='in')

ctf = axc3.contourf(xn, zn, un[:, :, tni], np.linspace(-0.2, 0.2, 41),
                    extend='both', cmap=cm.balance)

cbar_x0 = axc3.get_position().x1 + 0.01
cbar_y0 = axc3.get_position().y0
cbar_l0 = axc3.get_position().y1 - cbar_y0
cbar_ax = fig.add_axes([cbar_x0, cbar_y0, 0.01, cbar_l0])
fig.colorbar(ctf, cax=cbar_ax, ticks=np.linspace(-0.2, 0.2, 5),
             extendfrac=0, label='Velocity [m/s]')
cbar_ax.xaxis.tick_top()
cbar_ax.xaxis.set_label_position('top')

axz = fig.add_subplot(gs[:31, :46])
axz.set_xlim(0, 10)
axz.set_ylim(-2, 2)
axz.set_xlabel('X [km]')
axz.set_ylabel('Y [km]')
axz.plot([3, 3], [-2, 2], '--k', lw=0.5)
axz.tick_params(axis='both', direction='in')

dz = np.diff(zv)
dz = dz[zi:]
uuz = uu[zi:]
vuz = vu[zi:]
uuavg = np.zeros((eta2-eta1, xi2-xi1+1))
vuavg = np.zeros((eta2-eta1, xi2-xi1+1))
for i, dzi in enumerate(dz):
    uuavg = uuavg + uuz[i]*dzi
    vuavg = vuavg + vuz[i]*dzi
uuavg = uuavg/dz.sum()
vuavg = vuavg/dz.sum()
vel = np.sqrt(uuavg**2+vuavg**2)

zni = np.argmin(np.abs(zn+abs(z0)))
unr = un[:zni, :, tni]*np.cos(np.pi*ang/180) + \
    vn[:zni, :, tni]*np.sin(np.pi*ang/180)
vnr = -1 * un[:zni, :, tni]*np.sin(np.pi*ang/180) + \
    vn[:zni, :, tni]*np.cos(np.pi*ang/180)
unr = unr.mean(axis=0)
vnr = vnr.mean(axis=0)

ctf = axz.contourf(xv, yr, vel, np.linspace(0, 0.2, 101),
                   extend='both', cmap=cm.speed)
axz.quiver(xvel[::5], yvel[::5], unr[::5], vnr[::5],
           scale=2, alpha=0.75, color='r', angles='uv',
           width=0.005, minlength=0.01)
axz.quiver(xv[::2], yr, uuavg[:, ::2], vuavg[:, ::2],
           scale=2, alpha=0.75, angles='uv', minlength=0.01)
axz.plot(0, 0, '>r', markersize=15)

axz2 = fig.add_subplot(gs[:31, 48:94])
axz2.set_xlim(0, 10)
axz2.set_ylim(-2, 2)
axz2.set_xlabel('X [km]')
axz2.set_yticklabels([])
axz2.plot([3, 3], [-2, 2], '--k', lw=0.5)
axz2.tick_params(axis='both', direction='in')

dz = np.diff(zv)
dz = dz[zi:]
uuz2 = uu2[zi:]
vuz2 = vu2[zi:]
uuavg2 = np.zeros((eta2-eta1, xi2-xi1+1))
vuavg2 = np.zeros((eta2-eta1, xi2-xi1+1))
for i, dzi in enumerate(dz):
    uuavg2 = uuavg2 + uuz2[i]*dzi
    vuavg2 = vuavg2 + vuz2[i]*dzi
uuavg2 = uuavg2/dz.sum()
vuavg2 = vuavg2/dz.sum()
vel2 = np.sqrt(uuavg2**2+vuavg2**2)

zni = np.argmin(np.abs(zn+abs(z0)))
unr = un[:zni, :, tni]*np.cos(np.pi*ang/180) + \
    vn[:zni, :, tni]*np.sin(np.pi*ang/180)
vnr = -1 * un[:zni, :, tni]*np.sin(np.pi*ang/180) + \
    vn[:zni, :, tni]*np.cos(np.pi*ang/180)
unr = unr.mean(axis=0)
vnr = vnr.mean(axis=0)

ctf = axz2.contourf(xv, yr, vel2, np.linspace(0, 0.2, 101),
                    extend='both', cmap=cm.speed)
axz2.quiver(xvel[::5], yvel[::5], unr[::5], vnr[::5],
            scale=2, alpha=0.75, color='r', angles='uv',
            width=0.005, minlength=0.01)
axz2.quiver(xv[::2], yr, uuavg2[:, ::2], vuavg2[:, ::2],
            scale=2, alpha=0.75, angles='uv', minlength=0.01)
axz2.plot(0, 0, '>r', markersize=15)

cbar_x0 = axz2.get_position().x1 + 0.01
cbar_y0 = axz2.get_position().y0
cbar_l0 = axz2.get_position().y1 - cbar_y0
cbar_ax = fig.add_axes([cbar_x0, cbar_y0, 0.01, cbar_l0])
fig.colorbar(ctf, cax=cbar_ax, ticks=np.linspace(0, 0.2, 3),
             extendfrac=0, label='Flow Speed [m/s]')
cbar_ax.xaxis.tick_top()
cbar_ax.xaxis.set_label_position('top')

axz.text(0.02, 0.05, r'a) $\mathregular{H_{Mass}/OP}$',
         transform=axz.transAxes, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axc.text(0.03, 0.03, r'b) $\mathregular{H_{Mass}/OP}$',
         transform=axc.transAxes, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axz2.text(0.02, 0.05, 'c) Mix/NOP',
          transform=axz2.transAxes, ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axc2.text(0.03, 0.03, 'd) Mix/NOP',
          transform=axc2.transAxes, ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axc3.text(0.03, 0.03, 'e) Obs',
          transform=axc3.transAxes, ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

fig.savefig(out_dir+'fig3_3.png', dpi=300, bbox_inches='tight')
