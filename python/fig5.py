import pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import netCDF4 as nc
import pyroms
from cmocean import cm
import gsw

def get_zr(zeta, h, vgrid):
    """ get z at rho points from grid and zeta info. """

    ti = zeta.shape[0]
    zr = np.empty((ti, vgrid.N) + h.shape, 'd')
    if vgrid.Vtrans == 1:
        for k in range(vgrid.N):
            z0 = vgrid.hc * vgrid.s_rho[k] + (h - vgrid.hc) * vgrid.Cs_r[k]
            zr[:, k, :] = z0 + zeta * (1.0 + z0 / h)
    elif vgrid.Vtrans == 2 or vgrid.Vtrans == 4 or vgrid.Vtrans == 5:
        for k in range(vgrid.N):
            z0 = (vgrid.hc * vgrid.s_rho[k] + h * vgrid.Cs_r[k]) / (vgrid.hc + h)
            zr[:, k, :] = zeta + (zeta + h) * z0

    return zr


def get_zw(zeta, h, vgrid):
    """ get z at rho points from grid and zeta info. """

    ti = zeta.shape[0]
    zw = np.empty((ti, vgrid.Np) + h.shape, 'd')
    if vgrid.Vtrans == 1:
        for k in range(vgrid.Np):
            z0 = vgrid.hc * vgrid.s_w[k] + (h - vgrid.hc) * vgrid.Cs_w[k]
            zw[:, k, :] = z0 + zeta * (1.0 + z0 / h)
    elif vgrid.Vtrans == 2 or vgrid.Vtrans == 4 or vgrid.Vtrans == 5:
        for k in range(vgrid.Np):
            z0 = (vgrid.hc * vgrid.s_w[k] + h * vgrid.Cs_w[k]) / (vgrid.hc + h)
            zw[:, k, :] = zeta + (zeta + h) * z0

    return zw

z0raw_list = [0, 15, 25, 35, 50, 60, 70, 100, 150, 200]
t0 = 0
t1 = -1
dt = 1
qvstr = 1

rcase = 'iceplume_coupler'
app = '_mix_no_60s'
aver_file = '/Users/cw686/roms_archive/' + rcase + '/outputs' + app + '/fjord_avg.nc'
hist_file = '/Users/cw686/roms_archive/' + rcase + '/outputs' + app + '/fjord_his.nc'
out_dir   = '/Users/cw686/roms_archive/' + rcase + '/figs' + app + '/zview_uv/'
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

fh    = nc.Dataset(hist_file)
xmask = fh.variables['mask_rho'][0, :]
ymask = fh.variables['mask_rho'][:, 2]
x     = fh.variables['x_rho']   [0, :]
fh.close()
idx  = np.where(ymask==1)[0]
eta1 = idx.min()
eta2 = idx.max()+1

idx = np.where(xmask==1)[0]
xi1 = idx.min()
xi2 = len(x)-2

eta1 = eta1
eta2 = eta2

fh      = nc.Dataset(hist_file)
xr      = fh.variables['x_rho'] [0          , xi1:xi2  ]/1000
yr      = fh.variables['y_rho'] [eta1:eta2  , 0        ]/1000
xp      = fh.variables['x_psi'] [0          , xi1:xi2-1]/1000
yp      = fh.variables['y_psi'] [eta1:eta2-1, 0        ]/1000
h       = fh.variables['h']     [eta1:eta2  , xi1:xi2  ]
theta_b = fh.variables['theta_b'][0]
theta_s = fh.variables['theta_s'][0]
Tcline  = fh.variables['Tcline'][0]
N       = len(fh.dimensions['N'])
fh.close()
yoffset = (yr[0]+yr[-1])/2
yr = yr - yoffset
yp = yp - yoffset

dx = xr[1]-xr[0]
dy = yr[1]-yr[0]
xvert = np.arange(xr[0]-dx/2, xr[-1]+dx/2+0.01, dx)
yvert = np.arange(yr[0]-dy/2, yr[-1]+dy/2+0.01, dy)

vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=h)
zr   = vgrd.z_r[:][:, eta1, xi1]
zw   = vgrd.z_w[:][:, eta1, xi1]

zi_list = []
z0_list = []
for z0raw in z0raw_list:
    zi = np.argmin(np.abs(zr+abs(z0raw)))
    z0 = round(abs(zr[zi]))
    zi_list.append(zi)
    z0_list.append(z0)

fh   = nc.Dataset(aver_file)
if (t0 == 0) & (t1 == -1):
    time = fh.variables['ocean_time'][::dt]
    u = fh.variables['u'][::dt, :, eta1:eta2  , xi1:xi2-1]
    v = fh.variables['v'][::dt, :, eta1:eta2-1, xi1:xi2  ]
    w = fh.variables['w'][::dt, :, eta1:eta2  , xi1:xi2  ]
else:
    time = fh.variables['ocean_time'][t0:t1:dt]
    u = fh.variables['u'][t0:t1:dt, :, eta1:eta2  , xi1:xi2-1]
    v = fh.variables['v'][t0:t1:dt, :, eta1:eta2-1, xi1:xi2  ]
    w = fh.variables['w'][t0:t1:dt, :, eta1:eta2  , xi1:xi2  ]
fh.close()
u = 0.5*(u[:, :, 1:, :]+u[:, :, :-1, :])
v = 0.5*(v[:, :, :, 1:]+v[:, :, :, :-1])
w = 0.5*(w[:, 1:, :, :]+w[:, :-1, :, :])

for i, zi in enumerate(zi_list):
    fig, axs = plt.subplots(6, 4, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.02, hspace=0.05, left=0.05, right=0.90, bottom=0.05, top=0.94)
    axs[-1, 0].set_xlim(xvert[0], xvert[-1])
    axs[-1, 0].set_ylim(yvert[0], yvert[-1])
    axs[-1, 0].set_xlabel('X [km]')
    axs[-1, 0].set_ylabel('Y [km]')

    for ti, t in enumerate(time):
        ii = int(ti/4)
        jj = ti % 4
        pcm = axs[ii, jj].contourf(xr, yr, w[ti, zi]*1e2, np.linspace(-1, 1, 101), extend='both', cmap=cm.balance)
        qui = axs[ii, jj].quiver(xp[::qvstr], yp[::qvstr],
                u[ti, zi, ::qvstr, ::qvstr], v[ti, zi, ::qvstr, ::qvstr], scale=3, alpha=0.5, angles='uv')
        axs[ii, jj].text(0.9, 0.9, '%03d hr' % (t/3600), transform=axs[ii, jj].transAxes)
    y0 = axs[-1, 0].get_position().y0
    h0 = axs[ 0, 0].get_position().y1 - y0
    cbar_ax = fig.add_axes([0.905, y0, 0.012, h0])
    fig.colorbar(pcm, cax=cbar_ax, ticks=[0, 0.05], label=r'W [$10^{-2}$ m/s]', extendfrac='auto')
    rect = Rectangle((0.6, 0.1), 0.4, 0.35, transform=axs[-1, -1].transAxes, alpha=0.75, facecolor='w', edgecolor='w')
    axs[-1, -1].add_patch(rect)
    axs[-1, -1].quiverkey(qui, 0.8, 0.2, 0.2, '0.2 m/s')
    fig.savefig(out_dir + 'uv_%03d.png' % (z0_list[i]), dpi=300)
    plt.close()
