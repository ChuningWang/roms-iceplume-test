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
            z0 = (vgrid.hc * vgrid.s_rho[k] + h*vgrid.Cs_r[k]) / (vgrid.hc + h)
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


dt = 2

xstride = 10
zstride = 2

# s_list = [5, 10, 15, 20, 25, 30, 33, 34]
s_list = [5, 15, 25, 30, 33, 34]
q_list = [10, 25, 50, 100]

s_list = [33]
q_list = [50]

app2 = ''
# app2 = '_mix'

rcase = 'iceplume_2d'
for qi in q_list:
    for si in s_list:
        app = '_s%02d_q%03d' % (si, qi) + app2
        aver_file = '/home/chuning/projects/' + rcase + '/outputs' + app + '/fjord_avg.nc'
        hist_file = '/home/chuning/projects/' + rcase + '/outputs' + app + '/fjord_his.nc'
        out_dir = '/home/chuning/git/roms-iceplume-test/python/'
        # pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

        fh = nc.Dataset(hist_file)
        xmask = fh.variables['mask_rho'][0, :]
        ymask = fh.variables['mask_rho'][:, 2]
        x = fh.variables['x_rho'][0, :]
        fh.close()
        idx = np.where(ymask == 1)[0]
        eta1 = idx.min()+1
        eta2 = idx.max()

        idx = np.where(xmask == 1)[0]
        xi1 = idx.min()
        xi2 = len(x)-2

        eta1 = eta1
        eta2 = eta2

        fh = nc.Dataset(hist_file)
        xr = fh.variables['x_rho'][0, xi1:xi2]/1000
        yr = fh.variables['y_rho'][eta1:eta2, 0]/1000
        xp = fh.variables['x_psi'][0, xi1-1:xi2]/1000
        yp = fh.variables['y_psi'][eta1:eta2-1, 0]/1000
        h = fh.variables['h'][eta1:eta2, xi1:xi2]
        theta_b = fh.variables['theta_b'][0]
        theta_s = fh.variables['theta_s'][0]
        Tcline = fh.variables['Tcline'][0]
        N = len(fh.dimensions['N'])
        fh.close()
        yoffset = (yr[0]+yr[-1])/2
        yr = yr - yoffset
        yp = yp - yoffset

        dx = xr[1]-xr[0]
        dy = yr[1]-yr[0]
        xvert = np.arange(xr[0]-dx/2, xr[-1]+dx/2+0.01, dx)
        yvert = np.arange(yr[0]-dy/2, yr[-1]+dy/2+0.01, dy)

        # vgrd = pyroms.vgrid.s_coordinate_4(
        vgrd = pyroms.vgrid.SCoord(
            h, theta_b, theta_s, Tcline, N, hraw=h)
        zr = vgrd.z_r[:][:, eta1, xi1]
        zw = vgrd.z_w[:][:, eta1, xi1]
        zp = zw[1:-1]

        fh = nc.Dataset(aver_file)
        time = fh.variables['ocean_time'][::dt]
        dye1 = fh.variables['dye_01'][::dt, :, eta1:eta2, xi1:xi2]
        u = fh.variables['u'][::dt, :, eta1:eta2, xi1-1:xi2]
        v = fh.variables['v'][::dt, :, eta1:eta2-1, xi1:xi2]
        w = fh.variables['w'][::dt, :, eta1:eta2, xi1:xi2]
        rho = fh.variables['rho'][::dt, :, eta1:eta2, xi1:xi2]
        fh.close()
        fh = nc.Dataset(hist_file)
        rho0 = fh.variables['rho'][0, :, eta1:eta2, xi1:xi2]
        fh.close()
        u = 0.5*(u[:, :, 1:, :]+u[:, :, :-1, :])
        v = 0.5*(v[:, :, :, 1:]+v[:, :, :, :-1])
        w = 0.5*(w[:, 1:, :, :]+w[:, :-1, :, :])
        w = 0.25*(w[:, :, 1:, 1:]+w[:, :, 1:, :-1] +
                  w[:, :, :-1, 1:]+w[:, :, :-1, :-1])
        rhoD = 0.5*(rho0[0, 1, 1]+rho0[-1, 1, 1])

        zl = np.linspace(-200, 0, 9)
        rhol = np.interp(zl, zr, rho0[:, 1, 1])
        rhol = rhol[::-1]

        rhoa = rho-rho0

        rho2 = gsw.rho(35, 4, 125)
        rho1 = gsw.rho(si, 4, 25)
        gred = 9.8*(rho2-rho1)/1000
        cp = np.sqrt(gred*(50*150)/200)

        fig, axs = plt.subplots(4, 3, sharex=True, sharey=True)
        fig.subplots_adjust(wspace=0.02, hspace=0.05,
                            left=0.1, right=0.90, bottom=0.1, top=0.95)
        axs[-1, 0].set_xlim(xvert[0], xvert[-1])
        axs[-1, 0].set_ylim(zw[0], zw[-1])
        axs[-1, 0].set_xlabel('X [km]')
        axs[-1, 0].set_ylabel('Z [m]')
        axs[-1, 0].set_yticks([-150, -100, -50, 0])

        for ti, t in enumerate(time[:6]):
            ii = int(ti/3)
            jj = ti % 3
            hh = int(t/3600)
            mm = (t-int(t/3600)*3600)/60

            pcm = axs[ii, jj].contourf(
                xr, zr, dye1[ti].mean(axis=-2),
                np.linspace(0, 0.1, 51), extend='both', cmap=cm.matter)
            qui = axs[ii, jj].quiver(
                xp[::xstride], zr[::zstride],
                u[ti, ::zstride, :, ::xstride].mean(axis=-2),
                w[ti, ::zstride, :, ::xstride].mean(axis=-2),
                scale=5, alpha=0.5, angles='uv')
            axs[ii, jj].text(0.05, 0.05, '%02d hr %02d min' % (hh, mm),
                             transform=axs[ii, jj].transAxes, fontsize=15)
            axs[ii, jj].plot([xvert[0], xvert[-1]], [-50, -50], '--k', lw=0.7)
            axs[ii, jj].plot(cp*t/1000, 0, 'vk', ms=10)
            axs[ii, jj].plot([cp*t/1000, cp*t/1000], [0, -100], '--k')

        y0 = axs[1, 0].get_position().y0
        h0 = axs[0, 0].get_position().y1 - y0
        cbar_ax = fig.add_axes([0.905, y0, 0.012, h0])
        fig.colorbar(pcm, cax=cbar_ax, ticks=[0, 0.1],
                     label=r'Dye [kg/m$^3$]', extendfrac='auto')
        rect = Rectangle((0.7, 0.0), 0.3, 0.25,
                         transform=axs[1, -1].transAxes, alpha=0.75,
                         facecolor='w', edgecolor='w')
        axs[1, -1].add_patch(rect)
        axs[1, -1].quiverkey(qui, 0.85, 0.05, 0.5, '50 cm/s')

        for ti, t in enumerate(time[:6]):
            ii = int(ti/3)+2
            jj = ti % 3
            hh = int(t/3600)
            mm = (t-int(t/3600)*3600)/60

            pcm = axs[ii, jj].contourf(
                xr, zr, rhoa[ti].mean(axis=-2),
                np.linspace(-1, 1, 101), extend='both', cmap=cm.balance)
            axs[ii, jj].contour(
                xr, zr, rho[ti].mean(axis=-2),
                rhol, linewidths=0.3, colors='k')
            axs[ii, jj].text(0.05, 0.05, '%02d hr %02d min' % (hh, mm),
                             transform=axs[ii, jj].transAxes, fontsize=15)
            axs[ii, jj].plot(cp*t/1000, 0, 'vk', ms=10)
            axs[ii, jj].plot([cp*t/1000, cp*t/1000], [0, -100], '--k')

        y0 = axs[3, 0].get_position().y0
        h0 = axs[2, 0].get_position().y1 - y0
        cbar_ax = fig.add_axes([0.905, y0, 0.012, h0])
        fig.colorbar(pcm, cax=cbar_ax, ticks=[-1, 0, 1],
                     label=r'$\sigma$ [kg/m$^3$]', extendfrac='auto')

        fig.savefig(out_dir + 'fig3.png', dpi=300)
        plt.close()
