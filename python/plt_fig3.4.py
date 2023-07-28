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

init_file = '/home/chuning/projects/fjord/roms_files/roms_rst_ini.nc'
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
x0 = 20
x1 = 40
tni = 6
ang = 25+180

# ----------------------------------------------------------------------

vel_data = loadmat(vel_file)
tn = vel_data['timeN'].squeeze()
xn = vel_data['xN'].squeeze()
zn = vel_data['zN'].squeeze()*(-1)
un = vel_data['u_rot_N']
vn = vel_data['v_rot_N']

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

fh = nc.Dataset(init_file)
temp0 = fh.variables['temp'][:, :, eta1:eta2, xi1:xi2].mean(axis=0)
rho0 = fh.variables['rho'][:, :, eta1:eta2, xi1:xi2].mean(axis=0)
fh.close()

fh = nc.Dataset(roms_file[0])
time = fh.variables['ocean_time'][t0:]
temp = fh.variables['temp'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
rho = fh.variables['rho'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
dye = fh.variables['dye_01'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
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
rho2 = fh.variables['rho'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
dye2 = fh.variables['dye_01'][t0:, :, eta1:eta2, xi1:xi2].mean(axis=0)
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

np.savez('fig3.4data.npz',
        xr=xr, xv=xv, yr=yr, yv=yv, zr=zr, zv=zv, ur=ur, vr=vr,
        ur2=ur2, vr2=vr2, xvel=xvel, yvel=yvel,
        xn=xn, zn=zn, un=un, vn=vn, uu=uu, vu=vu, uu2=uu2, vu2=vu2)

# ----------------------------------------------------------------------

data = np.load('fig3.3data.npz')
xr, xv, yr, yv = data['xr'], data['xv'], data['yr'], data['yv']
zr, zv = data['zr'], data['zv']
ur, vr, ur2, vr2 = data['ur'], data['vr'], data['ur2'], data['vr2']
xvel, yvel = data['xvel'], data['yvel']
xn, zn, un, vn = data['xn'], data['zn'], data['un'], data['un']
uu, vu, uu2, vu2 = data['uu'], data['vu'], data['uu2'], data['vu2']

zi = np.argmin(np.abs(zr+abs(z0)))
xi00 = np.argmin(np.abs(xr-x0))
xi01 = np.argmin(np.abs(xr-x1))

fig = plt.figure()
fig.set_size_inches(8, 7)
gs = GridSpec(100, 100, figure=fig, wspace=0.05, hspace=0.05,
              left=0.06, right=0.93, bottom=0.05, top=0.95)

axz = fig.add_subplot(gs[:21, :48])
axz.set_xlim(0, 60)
axz.set_ylim(-2, 2)
axz.set_xlabel('X [km]')
axz.set_ylabel('Y [km]')
axz.tick_params(axis='both', direction='in')

axz2 = fig.add_subplot(gs[:21, 50:98])
axz2.set_xlim(0, 60)
axz2.set_ylim(-2, 2)
axz2.set_xlabel('X [km]')
axz2.set_yticklabels([])
axz2.tick_params(axis='both', direction='in')

axc = fig.add_subplot(gs[30:60, 15:48])
axc.set_xlim(-2, 2)
axc.set_ylim(-150, 0)
axc.set_yticks([-150, -100, -50, 0])
axc.set_yticklabels([])
axc.set_xlabel('Y [km]')
axc.tick_params(axis='both', direction='in')

axc2 = fig.add_subplot(gs[30:60, 65:98])
axc2.set_xlim(-2, 2)
axc2.set_ylim(-150, 0)
axc2.set_yticks([-150, -100, -50, 0])
axc2.set_yticklabels([])
axc2.set_xlabel('Y [km]')
axc2.tick_params(axis='both', direction='in')

axp = fig.add_subplot(gs[30:60, :13])
axp.set_xlabel('U [cm/s]')
axp.set_xlim(-3, 6)
axp.set_ylim(-150, 0)
axp.set_xticks([0, 5])
axp.set_yticks([-150, -100, -50, 0])
axp.set_yticklabels(['150', '100', '50', '0'])
axp.set_yticks([-150, -100, -50, 0])
axp.set_ylabel('Z [m]')
axp.tick_params(axis='both', direction='in')

# axpv = axp.twiny()
# axpv.set_xlim(-0.3, 0.6)
# axpv.set_xlabel('V [cm/s]', color='r')
# axpv.tick_params(axis='x', colors='r')

axp2 = fig.add_subplot(gs[30:60, 50:63])
axp2.set_xlabel('U [cm/s]')
axp2.set_xlim(-3, 6)
axp2.set_ylim(-150, 0)
axp2.set_xticks([0, 5])
axp2.set_yticklabels([])
axp2.tick_params(axis='both', direction='in')

# axpv2 = axp2.twiny()
# axpv2.set_xlim(-0.3, 0.6)
# axpv2.set_xlabel('V [cm/s]', color='r')
# axpv2.tick_params(axis='x', colors='r')

axt = fig.add_subplot(gs[70:100, 15:48])
axt.set_xlim(-2, 2)
axt.set_ylim(-150, 0)
axt.set_yticks([-150, -100, -50, 0])
axt.set_yticklabels([])
axt.set_xlabel('Y [km]')
axt.tick_params(axis='both', direction='in')

axt2 = fig.add_subplot(gs[70:100, 65:98])
axt2.set_xlim(-2, 2)
axt2.set_ylim(-150, 0)
axt2.set_yticks([-150, -100, -50, 0])
axt2.set_yticklabels([])
axt2.set_xlabel('Y [km]')
axt2.tick_params(axis='both', direction='in')

axtp = fig.add_subplot(gs[70:100, :13])
axtp.set_xlabel(r'Temp [$^{\circ}$C]')
axtp.set_xlim(0, 3)
axtp.set_ylim(-150, 0)
axtp.set_xticks([0, 1, 2])
axtp.set_yticks([-150, -100, -50, 0])
axtp.set_yticklabels(['150', '100', '50', '0'])
axtp.set_ylabel('Z [m]')
axtp.tick_params(axis='both', direction='in')

axtp2 = fig.add_subplot(gs[70:100, 50:63])
axtp2.set_xlabel(r'Temp [$^{\circ}$C]')
axtp2.set_xlim(0, 3)
axtp2.set_ylim(-150, 0)
axtp2.set_xticks([0, 1, 2])
axtp2.set_yticks([-150, -100, -50, 0])
axtp2.set_yticklabels([])
axtp2.tick_params(axis='both', direction='in')

# ----------------------------------------------------------------------

axz.plot([x0, x0], [-2, 2], '--k')
axz.plot([x1, x1], [-2, 2], '--k')

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
vnr = - un[:zni, :, tni]*np.sin(np.pi*ang/180) + \
      vn[:zni, :, tni]*np.cos(np.pi*ang/180)
unr = unr.mean(axis=0)
vnr = vnr.mean(axis=0)

ctf = axz.contourf(xv, yr, vel, np.linspace(0, 0.2, 101),
                   extend='both', cmap=cm.speed)
axz.quiver(xvel[::5], yvel[::5], unr[::5], vnr[::5],
           scale=2, alpha=0.75, color='r', angles='uv',
           width=0.005, minlength=0.01)
axz.quiver(xv[::10], yr, uuavg[:, ::10], vuavg[:, ::10],
           scale=2, alpha=0.75, angles='uv', minlength=0.01)
axz.plot(0, 0, '>r', markersize=15)


axz2.plot([x0, x0], [-2, 2], '--k')
axz2.plot([x1, x1], [-2, 2], '--k')

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
vnr = - un[:zni, :, tni]*np.sin(np.pi*ang/180) + \
      vn[:zni, :, tni]*np.cos(np.pi*ang/180)
unr = unr.mean(axis=0)
vnr = vnr.mean(axis=0)

ctf = axz2.contourf(xv, yr, vel2, np.linspace(0, 0.2, 101),
                    extend='both', cmap=cm.speed)
axz2.quiver(xvel[::5], yvel[::5], unr[::5], vnr[::5],
            scale=2, alpha=0.75, color='r', angles='uv',
            width=0.005, minlength=0.01)
axz2.quiver(xv[::10], yr, uuavg2[:, ::10], vuavg2[:, ::10],
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

# ----------------------------------------------------------------------

vel = ur[:, :, xi00:xi01] + 1j*vr[:, :, xi00:xi01]
vel = vel.mean(axis=-1)

ctf = axc.contourf(yr, zr, np.real(vel), np.linspace(-0.2, 0.2, 41),
                   extend='both', cmap=cm.balance)

vel2 = ur2[:, :, xi00:xi01] + 1j*vr2[:, :, xi00:xi01]
vel2 = vel2.mean(axis=-1)

ctf = axc2.contourf(yr, zr, np.real(vel2), np.linspace(-0.2, 0.2, 41),
                    extend='both', cmap=cm.balance)

cbar_x0 = axc2.get_position().x1 + 0.01
cbar_y0 = axc2.get_position().y0
cbar_h0 = axc2.get_position().y1 - cbar_y0
cbar_ax = fig.add_axes([cbar_x0, cbar_y0, 0.01, cbar_h0])
fig.colorbar(ctf, cax=cbar_ax, ticks=np.linspace(-0.2, 0.2, 5),
             extendfrac=0, label='Velocity [m/s]')

axp.plot([0, 0], [-800, 0], '--', color='gray', lw=0.5)
axp.plot(np.real(vel).mean(axis=1)*100, zr, 'b', label='U')
axp.plot(0, -1000, 'r', label='V')
# axp.legend(fontsize=10)
# axpv.plot(np.imag(vel).mean(axis=1)*100, zr, 'r')


axp2.plot([0, 0], [-800, 0], '--', color='gray', lw=0.5)
axp2.plot(np.real(vel2).mean(axis=1)*100, zr, 'b', label='U')
axp2.plot(0, -1000, 'r', label='V')
# axp2.legend(fontsize=10)
# axpv2.plot(np.imag(vel2).mean(axis=1)*100, zr, 'r')

# ----------------------------------------------------------------------

clevs = [0.001, 1, 2, 3.0, 3.4]
lws = [1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]

tc0 = temp0[:, :, xi00:xi01].mean(axis=-1)

tc = temp[:, :, xi00:xi01].mean(axis=-1)
dc = dye[:, :, xi00:xi01].mean(axis=-1)*100

ctf = axt.contourf(yr, zr, tc-tc0, np.linspace(-1, 1, 41),
                   extend='both', cmap=cm.balance)
ct = axt.contour(yr, zr, dc, clevs,
                 colors='k', linewidths=lws)
ct.clabel(fontsize=8, fmt='%1.1f')

tc2 = temp2[:, :, xi00:xi01].mean(axis=-1)
dc2 = dye2[:, :, xi00:xi01].mean(axis=-1)*100

ctf = axt2.contourf(yr, zr, tc2-tc0, np.linspace(-1, 1, 41),
                    extend='both', cmap=cm.balance)
ct = axt2.contour(yr, zr, dc2, clevs,
                  colors='k', linewidths=lws)
ct.clabel(fontsize=8, fmt='%1.1f')

cbar_x0 = axt2.get_position().x1 + 0.01
cbar_y0 = axt2.get_position().y0
cbar_h0 = axt2.get_position().y1 - cbar_y0
cbar_ax = fig.add_axes([cbar_x0, cbar_y0, 0.01, cbar_h0])
fig.colorbar(ctf, cax=cbar_ax, ticks=np.linspace(-1, 1, 5),
             extendfrac=0, label=r'Temp Anomaly [$^{\circ}$C]')

axtp.plot(tc.mean(axis=1), zr, label='Model')
axtp.plot(tc0.mean(axis=1), zr, 'k', label='Init')
axtp.legend(fontsize=7)

axtp2.plot(tc2.mean(axis=1), zr, label='Model')
axtp2.plot(tc0.mean(axis=1), zr, 'k', label='Init')
axtp2.legend(fontsize=7)

# ----------------------------------------------------------------------

axz.text(0.98, 0.05, r'a) $\mathregular{H_{Mass}/OP}$',
         transform=axz.transAxes, ha='right', va='bottom',
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axc.text(0.03, 0.03, 'c)',
         transform=axc.transAxes, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axp.text(0.05, 0.03, 'b)',
         transform=axp.transAxes, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axt.text(0.03, 0.03, 'e)',
         transform=axt.transAxes, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axtp.text(0.05, 0.03, 'd)',
          transform=axtp.transAxes, ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axz2.text(0.98, 0.05, 'f) Mix/NOP',
          transform=axz2.transAxes, ha='right', va='bottom',
          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axc2.text(0.03, 0.03, 'h)',
          transform=axc2.transAxes, ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axp2.text(0.05, 0.03, 'g)',
          transform=axp2.transAxes, ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axt2.text(0.03, 0.03, 'j)',
          transform=axt2.transAxes, ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
axtp2.text(0.05, 0.03, 'i)',
           transform=axtp2.transAxes, ha='left', va='bottom',
           bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

# ----------------------------------------------------------------------

fig.savefig(out_dir+'fig3_4.png', dpi=300, bbox_inches='tight')
plt.close()
