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

qvstrx = [6, 2, 1]
qvstry = 2

# z0raw_list = [0, 15, 25, 35, 50, 60, 70, 100, 150, 200]
z0raw_list = [50]

det_list = [r'$\mathregular{H_{Mass}/OP}$', 'Mix/NOP',
            r'$\mathregular{H_{Mass}/OP}$', 'Mix/NOP']
# dx_list = ['100 m', '200 m', '300 m', '400 m', '600 m']
dx_list = ['100 m', '300 m', '600 m']

# app_list = [['luv_det_100m', 'luv_det_200m', 'luv_det_300m', 'luv_det_400m', 'luv_det_600m'],
#             ['mix_no_100m',  'mix_no_200m',  'mix_no_300m',  'mix_no_400m',  'mix_no_600m' ]]
# blowup_list = [[False, False, False, False, False],
#                [False, False, False, False, False],]
app_list = [['luv_det_100m', 'luv_det_300m', 'luv_det_600m'],
            ['mix_no_100m',  'mix_no_300m',  'mix_no_600m' ]]
blowup_list = [[False, False, False],
               [False, False, False],]

rcase = 'iceplume_coupler'
app = app_list[0][0]
aver_file = '/home/chuning/projects/' + rcase + '/outputs_' + app + '/fjord_avg.nc'
hist_file = '/home/chuning/projects/' + rcase + '/outputs_' + app + '/fjord_his.nc'
# out_dir   = '/home/chuning/projects/' + rcase + '/figs/cmp/zview_uv/'
out_dir   = '/mnt/c/Users/chuning/Desktop/zview_uv/'
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size': 20})

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

fh = nc.Dataset(hist_file)
xr = fh.variables['x_rho'] [0          , xi1:xi2  ]/1000
xu = fh.variables['x_u']   [0          , xi1-1:xi2]/1000
yr = fh.variables['y_rho'] [eta1:eta2  , 0        ]/1000
fh.close()
xoffset = xu[0]
xr = xr - xoffset
yoffset = (yr[0]+yr[-1])/2
yr = yr - yoffset

dx = xr[1]-xr[0]
dy = yr[1]-yr[0]
xvert0 = np.arange(xr[0]-dx/2, xr[-1]+dx/2+0.01, dx)
yvert0 = np.arange(yr[0]-dy/2, yr[-1]+dy/2+0.01, dy)

for z0raw in z0raw_list:
    fig, axs = plt.subplots(2, len(dx_list), sharex=True, sharey=True, figsize=(12, 6))
    fig.subplots_adjust(wspace=0.02, hspace=0.05, left=0.05, right=0.88, bottom=0.1, top=0.94)
    axs[-1, 0].set_xlim(xvert0[0], xvert0[-1])
    axs[-1, 0].set_ylim(yvert0[0], yvert0[-1])
    axs[-1, 0].set_xlabel('X [km]')
    axs[-1, 0].set_ylabel('Y [km]')

    for ii in range(2):
        axs[ii, -1].set_ylabel(det_list[ii])
        axs[ii, -1].yaxis.set_label_position('right')
        # axs[ii+2, -1].set_ylabel(det_list[ii])
        # axs[ii+2, -1].yaxis.set_label_position('right')
    for jj in range(len(dx_list)):
        axs[0, jj].set_xlabel(dx_list[jj])
        axs[0, jj].xaxis.set_label_position('top')

    for ii in range(2):
        for jj in range(len(dx_list)):
            if blowup_list[ii][jj]:
                axs[ii, jj].text(0.5, 0.5, 'Blow Up', ha='center', va='center', transform=axs[ii, jj].transAxes)
            else:
                app = app_list[ii][jj]
                aver_file = '/home/chuning/projects/' + rcase + '/outputs_' + app + '/fjord_avg.nc'
                hist_file = '/home/chuning/projects/' + rcase + '/outputs_' + app + '/fjord_his.nc'

                fh      = nc.Dataset(hist_file)
                xr      = fh.variables['x_rho'] [0          , xi1:xi2  ]/1000
                yr      = fh.variables['y_rho'] [eta1:eta2  , 0        ]/1000
                xp      = fh.variables['x_psi'] [0          , xi1:xi2-1]/1000
                yp      = fh.variables['y_psi'] [eta1:eta2-1, 0        ]/1000
                xu      = fh.variables['x_u']   [0,           xi1-1:xi2]/1000
                h       = fh.variables['h']     [eta1:eta2  , xi1:xi2  ]
                theta_b = fh.variables['theta_b'][0]
                theta_s = fh.variables['theta_s'][0]
                Tcline  = fh.variables['Tcline'][0]
                N       = len(fh.dimensions['N'])
                fh.close()
                xoffset = xu[0]
                xr = xr - xoffset
                xp = xp - xoffset
                xu = xu - xoffset
                yoffset = (yr[0]+yr[-1])/2
                yr = yr - yoffset
                yp = yp - yoffset

                dx = xr[1]-xr[0]
                dy = yr[1]-yr[0]
                xvert = np.arange(xr[0]-dx/2, xr[-1]+dx/2+0.01, dx)
                yvert = np.arange(yr[0]-dy/2, yr[-1]+dy/2+0.01, dy)

                vgrd = pyroms.vgrid.SCoord(h, theta_b, theta_s, Tcline, N, hraw=h)
                zr   = vgrd.z_r[:][:, eta1, xi1]
                zw   = vgrd.z_w[:][:, eta1, xi1]
                zp   = zw[1:-1]

                zi = np.argmin(np.abs(zr+abs(z0raw)))
                z0 = round(abs(zr[zi]))

                fh   = nc.Dataset(aver_file)
                u = fh.variables['u'][-1, :, eta1:eta2  , xi1:xi2-1]
                v = fh.variables['v'][-1, :, eta1:eta2-1, xi1:xi2  ]
                w = fh.variables['w'][-1, :, eta1:eta2  , xi1:xi2  ]
                dye = fh.variables['dye_01'][-1, :, eta1:eta2, xi1:xi2]
                fh.close()
                u = 0.5*(u[:, 1:, :]+u[:, :-1, :])
                v = 0.5*(v[:, :, 1:]+v[:, :, :-1])
                w = 0.5*(w[1:, :, :]+w[:-1, :, :])

                pcmd = axs[ii, jj].contourf(xr, yr, dye[zi]*1e2, np.linspace(0, 5, 101), extend='both', cmap=cm.matter)
                axs[ii, jj].contour(xr, yr, dye[zi]*1e2, np.linspace(1, 3, 11), colors='w', linewidths=0.7)
                axs[ii, jj].contour(xr, yr, dye[zi]*1e2, [2], colors='k', linewidths=2)
                # pcm = axs[ii, jj].contourf(xr, yr, w[zi]*1e2, np.linspace(-1, 1, 101), extend='both', cmap=cm.balance)
                qui = axs[ii, jj].quiver(xp[::qvstrx[jj]], yp[::qvstry],
                        u[zi, ::qvstry, ::qvstrx[jj]], v[zi, ::qvstry, ::qvstrx[jj]],
                        width=0.0075, headwidth=2, scale=3, alpha=0.5, angles='uv')
    # y0 = axs[1, 0].get_position().y0
    # h0 = axs[0, 0].get_position().y1 - y0
    # cbar_ax = fig.add_axes([0.92, y0, 0.012, h0])
    # fig.colorbar(pcm, cax=cbar_ax, ticks=[-1, 0, 1], label=r'W [10$^{-2}$ m/s]', extendfrac='auto')
    rect = Rectangle((0.65, 0.0), 0.35, 0.225, transform=axs[1, -1].transAxes, alpha=0.75, facecolor='w', edgecolor='w')
    axs[1, -1].add_patch(rect)
    axs[1, -1].quiverkey(qui, 0.825, 0.05, 0.5, '0.5 m/s')
    y0 = axs[1, 0].get_position().y0
    h0 = axs[0, 0].get_position().y1 - y0
    cbar_ax = fig.add_axes([0.92, y0, 0.012, h0])
    fig.colorbar(pcmd, cax=cbar_ax, ticks=[0, 5], label=r'Dye [10$^{-2}$ kg/m$^3$]', extendfrac='auto')
    fig.savefig(out_dir + 'uv_dx_%03dm.png' % z0, dpi=300, bbox_inches='tight')
    plt.close()
