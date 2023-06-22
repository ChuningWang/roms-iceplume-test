import subprocess
from datetime import datetime
import numpy as np
import netCDF4 as nc
import pyroms
import pyroms_toolbox
from scipy.io import loadmat

# ----------------------------------------------------------------------
# Grid Construction

# grid name
grd_name = 'fjord'
grd_filename = '../roms_files/roms_grd.nc'

# basic paramters
add_sponge = True

# grid dimension
xsize = 60000.
ysize = 4180.
dx    = 300.
dy    = 220.
ddxs  = 50.
ddys  = 50.

xland  = np.arange(-2*dx, 1, dx)
xfjord = np.arange(dx, xsize+1, dx)
yfjord = np.arange(-1*ysize/2, ysize/2+1, dy)

xshelf = [xfjord[-1]]
for i in range(50):
    dx = dx + ddxs
    xshelf.append(xshelf[i]+dx)
xshelf = np.array(xshelf[1:])

yshelf = [yfjord[-1]]
for i in range(25):
    dy = dy + ddys
    yshelf.append(yshelf[i]+dy)
yshelf = np.array(yshelf[1:])

xvert = np.concatenate((xland, xfjord, xshelf))
yvert = np.concatenate((-yshelf[::-1], yfjord, yshelf))
Xvert = len(xvert)
Yvert = len(yvert)
Xrho  = Xvert-1
Yrho  = Yvert-1

# fjord depth and Coriolis Parameter
Dm   = 400.
f0   = 1.379*(1.e-4)
beta = 0.

# vertical grid specs
theta_s = 4.0
theta_b = 2.0
Tcline  = 100.
N       = 50

# horizontal grid construction
xrho = 0.5*(xvert[1:] + xvert[:-1])
yrho = 0.5*(yvert[1:] + yvert[:-1])
# meshgrid
xxvert, yyvert = np.meshgrid(xvert, yvert)
xxrho,  yyrho  = np.meshgrid(xrho, yrho)

# generate land mask
mask_rho        = np.ones((Yrho, Xrho))
mask_rho[:, :2] = 0.
land = (np.abs(yyrho)>yfjord[-1]) & (xxrho<xfjord[-1])
mask_rho[land]  = 0.

# write hgrid
hgrd          = pyroms.hgrid.CGrid(xxvert, yyvert)
# Coriolis Parameter
hgrd.f        = f0 + hgrd.y_rho*beta
hgrd.mask_rho = mask_rho

# vertical grid construction
h    = Dm*np.ones((Yrho, Xrho))

# write vertical grid
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=h)

# write grid
grd  = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
pyroms.grid.write_ROMS_grid(grd, filename=grd_filename)

if add_sponge:
    # add sponge layer
    len_sponge = 10000
    print('MAX dx : %08d' % grd.hgrid.dx.max())
    print('MAX dy : %08d' % grd.hgrid.dy.max())
    grid_size = grd.hgrid.dx*grd.hgrid.dy
    grid_size = grid_size/grid_size.min()
    x = hgrd.x_rho
    y = hgrd.y_rho
    dist = np.ones((Yrho, Xrho))
    visc_factor = np.ones((Yrho, Xrho))
    diff_factor = np.ones((Yrho, Xrho))
    for i in range(Yrho):
        for j in range(Xrho):
            dist[i, j] = min([np.abs(y[i, j]-y[ 0, 0]),
                              np.abs(y[i, j]-y[-1, 0])])
    visc_factor = (1-np.tanh(np.pi*(dist-len_sponge)/len_sponge))*0.5*9+1
    diff_factor = (1-np.tanh(np.pi*(dist-len_sponge)/len_sponge))*0.5*9+1
    visc_factor = visc_factor/grid_size
    diff_factor = diff_factor/grid_size

    fh = nc.Dataset(grd_filename, 'r+')

    fh.createVariable('visc_factor', 'f8', ('eta_rho', 'xi_rho'))
    fh.variables['visc_factor'].long_name     = 'horizontal viscosity sponge factor'
    fh.variables['visc_factor'].valid_min     = 0.
    fh.variables['visc_factor'].coordinates   = 'x_rho y_rho'
    fh.variables['visc_factor'].density_point = 1
    fh.variables['visc_factor'][:]            = visc_factor

    fh.createVariable('diff_factor', 'f8', ('eta_rho', 'xi_rho'))
    fh.variables['diff_factor'].long_name     = 'horizontal diffusivity sponge factor'
    fh.variables['diff_factor'].valid_min     = 0.
    fh.variables['diff_factor'].coordinates   = 'x_rho y_rho'
    fh.variables['diff_factor'].density_point = 1
    fh.variables['diff_factor'][:]            = diff_factor

    fh.close()

# ----------------------------------------------------------------------
# Make river file

# river file name
river_filename = '../roms_files/roms_rivers.nc'

# basic parameters
tracer1d   = False
river_time = np.arange(200., 500.001, 1)
ntime      = len(river_time)
sgtrs0     = 200.

ymask      = mask_rho[:, 2]
xpos       = np.where(ymask)[0]
xpos0      = float(round(0.5*(xpos[0]+xpos[-1])))
ypos0      = 2
ypos       = ypos0*np.ones(xpos.shape).astype(int)

nriver     = len(ypos)
river_dir  = np.zeros(nriver)
river_id   = np.ones(nriver)

if tracer1d:
    river_tracer_dim    = (ntime)
    sg_tracer_dim       = (ntime)
    nc_river_tracer_dim = ('river_time')
    nc_sg_tracer_dim    = ('river_time')
else:
    river_tracer_dim    = (ntime, N, nriver)
    sg_tracer_dim       = (ntime, nriver)
    nc_river_tracer_dim = ('river_time', 's_rho', 'river')
    nc_sg_tracer_dim    = ('river_time', 'river')

rtrs  = np.ones((ntime, nriver))*0.
rtemp = np.ones(river_tracer_dim)*0.
rsalt = np.ones(river_tracer_dim)*0.
rdye  = np.ones(river_tracer_dim)*0.

v_shape    = np.zeros((N, nriver))
v_shape[:] = 1./N

# For subglacial runoff
sgdepth  = np.ones(nriver)*(-260)
sgtype   = np.ones(nriver)*1
sglength = np.ones(nriver)*220/nriver
sggid    = np.ones(nriver)

sgtrs    = np.ones((ntime, nriver))*0.
sgtemp   = np.ones(sg_tracer_dim)*0.
sgsalt   = np.ones(sg_tracer_dim)*0.
sgdye    = np.ones(sg_tracer_dim)*1.

# Surface and subglacial discharge
for i in range(nriver):
    if xpos[i] == xpos0:
        sgtrs[:, i] = sgtrs0*np.tanh(0.1*(river_time-200))/nriver
        sgtype[i]   = 4  # sheet
sgtrs[sgtrs<0.] = 0.

# create file with all the objects
fh        = nc.Dataset(river_filename, 'w')
fh.type   = 'ROMS RIVERS file'
fh.title  = 'Fjord test'
fh.source = 'Analytical'

fh.createDimension('river_time', None)
fh.createDimension('river', nriver)
fh.createDimension('s_rho', N)
fh.createDimension('loc', 2)

fh.createVariable('river_time', 'f8', ('river_time'))
fh.variables['river_time'][:].units = 'days since 1900-01-01 00:00:00'
fh.variables['river_time'][:].long_name = 'river runoff time'
fh.variables['river_time'][:] = river_time

fh.createVariable('river', 'i4', ('river'))
fh.variables['river'].long_name = 'river runoff identification number'
fh.variables['river'][:] = river_id

fh.createVariable('river_Eposition', 'i4', ('river'))
fh.variables['river_Eposition'].long_name = 'river ETA-position at RHO-points'
fh.variables['river_Eposition'][:] = xpos

fh.createVariable('river_Xposition', 'i4', ('river'))
fh.variables['river_Xposition'].long_name = 'river XI-position at RHO-points'
fh.variables['river_Xposition'][:] = ypos

fh.createVariable('river_direction', 'i4', ('river'))
fh.variables['river_direction'].long_name = 'river runoff direction'
fh.variables['river_direction'][:] = river_dir

fh.createVariable('river_transport', 'f8', ('river_time', 'river'))
fh.variables['river_transport'].long_name = 'river runoff vertically integrated mass transport'
fh.variables['river_transport'].units = 'meter3 second-1'
fh.variables['river_transport'].time = 'river_time'
fh.variables['river_transport'][:] = rtrs

fh.createVariable('river_Vshape', 'f8', ('s_rho', 'river'))
fh.variables['river_Vshape'].long_name = 'river runoff mass transport vertical profile'
fh.variables['river_Vshape'][:] = v_shape

fh.createVariable('river_temp', 'f8', nc_river_tracer_dim)
fh.variables['river_temp'].long_name = 'river runoff potential temperature'
fh.variables['river_temp'].units = 'Celsius'
fh.variables['river_temp'].time = 'river_time'
fh.variables['river_temp'][:] = rtemp

fh.createVariable('river_salt', 'f8', nc_river_tracer_dim)
fh.variables['river_salt'].long_name = 'river runoff salinity'
fh.variables['river_salt'].time = 'river_time'
fh.variables['river_salt'][:] = rsalt

fh.createVariable('river_dye_01', 'f8', nc_river_tracer_dim)
fh.variables['river_dye_01'].long_name = 'river runoff dye 01 concentration'
fh.variables['river_dye_01'].time = 'river_time'
fh.variables['river_dye_01'][:] = rdye

fh.createVariable('river_dye_02', 'f8', nc_river_tracer_dim)
fh.variables['river_dye_02'].long_name = 'river runoff dye 02 concentration'
fh.variables['river_dye_02'].time = 'river_time'
fh.variables['river_dye_02'][:] = rdye

fh.createVariable('river_dye_03', 'f8', nc_river_tracer_dim)
fh.variables['river_dye_03'].long_name = 'river runoff dye 03 concentration'
fh.variables['river_dye_03'].time = 'river_time'
fh.variables['river_dye_03'][:] = rdye

# for subglacial runoff
fh.createVariable('subglacial_depth', 'f8', ('river'))
fh.variables['subglacial_depth'].units = 'meter'
fh.variables['subglacial_depth'][:] = sgdepth

fh.createVariable('subglacial_type', 'f8', ('river'))
fh.variables['subglacial_type'].units = 'nondimensional'
fh.variables['subglacial_type'][:] = sgtype

fh.createVariable('subglacial_length', 'f8', ('river'))
fh.variables['subglacial_length'].units = 'meter'
fh.variables['subglacial_length'][:] = sglength

fh.createVariable('subglacial_id', 'f8', ('river'))
fh.variables['subglacial_id'].units = 'nondimensional'
fh.variables['subglacial_id'][:] = sggid

fh.createVariable('subglacial_transport', 'f8', ('river_time', 'river'))
fh.variables['subglacial_transport'].long_name = 'subglacial runoff mass transport'
fh.variables['subglacial_transport'].units = 'meter3 second-1'
fh.variables['subglacial_transport'].time = 'river_time'
fh.variables['subglacial_transport'][:] = sgtrs

fh.createVariable('subglacial_temp', 'f8', nc_sg_tracer_dim)
fh.variables['subglacial_temp'].long_name = 'subglacial runoff potential temperature'
fh.variables['subglacial_temp'].units = 'Celsius'
fh.variables['subglacial_temp'].time = 'river_time'
fh.variables['subglacial_temp'][:] = sgtemp

fh.createVariable('subglacial_salt', 'f8', nc_sg_tracer_dim)
fh.variables['subglacial_salt'].long_name = 'subglacial runoff salinity'
fh.variables['subglacial_salt'].time = 'river_time'
fh.variables['subglacial_salt'][:] = sgsalt

fh.createVariable('subglacial_dye_01', 'f8', nc_sg_tracer_dim)
fh.variables['subglacial_dye_01'].long_name = 'subglacial runoff dye 01 concentration'
fh.variables['subglacial_dye_01'].time = 'river_time'
fh.variables['subglacial_dye_01'][:] = sgdye

fh.createVariable('subglacial_dye_02', 'f8', nc_sg_tracer_dim)
fh.variables['subglacial_dye_02'].long_name = 'subglacial runoff dye 02 concentration'
fh.variables['subglacial_dye_02'].time = 'river_time'
fh.variables['subglacial_dye_02'][:] = sgdye

fh.createVariable('subglacial_dye_03', 'f8', nc_sg_tracer_dim)
fh.variables['subglacial_dye_03'].long_name = 'subglacial runoff dye 03 concentration'
fh.variables['subglacial_dye_03'].time = 'river_time'
fh.variables['subglacial_dye_03'][:] = sgdye

fh.createVariable('subglacial_Erange', 'i4', ('loc', 'river'))
fh.variables['subglacial_Erange'].long_name = 'subglacial runoff average indices, ETA direction'
fh.variables['subglacial_Erange'][0, :] = xpos0-5
fh.variables['subglacial_Erange'][1, :] = xpos0+5

fh.createVariable('subglacial_Xrange', 'i4', ('loc', 'river'))
fh.variables['subglacial_Xrange'].long_name = 'subglacial runoff average indices, XI direction'
fh.variables['subglacial_Xrange'][0, :] = ypos0
fh.variables['subglacial_Xrange'][1, :] = ypos0+5

fh.close()

# ----------------------------------------------------------------------
# Make initial file

# initial file name
rcase_i = 'fjord_spinup'
rcase_o = 'fjord_melt'
ini_filename     = '../roms_files/roms_ini.nc'
rst_filename     = '../roms_archive/' + rcase_i + '/outputs/roms_rst.nc'
rst_ini_filename = '../roms_archive/' + rcase_o + '/roms_rst_ini.nc'

# basic parameter
use_rst = False

if use_rst:
    # generate initial file from restart file
    fh = nc.Dataset(rst_filename)
    rst_time = fh.variables['ocean_time'][:]
    rst_idx  = np.argmax(rst_time)
    fh.close()
    cmd = 'ncks -O -d ocean_time,%01d '  % rst_idx + rst_filename + ' ' + rst_ini_filename
    subprocess.call(cmd, shell=True)
else:
    # create initial file from CTD data
    ctddata = loadmat('KS2014_KS1_TS.mat')
    salt0 = np.nanmean(ctddata['july2829'][0, 0][0], axis=1)
    temp0 = np.nanmean(ctddata['july2829'][0, 0][1], axis=1)
    depth0 = -1*ctddata['july2829'][0, 0][2].squeeze()
    depth = grd.vgrid.z_r[:][:, 0, 0]
    spval = -1.0e20

    salt0 = salt0[::-1]
    temp0 = temp0[::-1]
    depth0 = depth0[::-1]
    salt0[np.isnan(salt0)] = salt0[498]
    temp0[np.isnan(temp0)] = temp0[498]

    salt = np.interp(depth, depth0, salt0)
    temp = np.interp(depth, depth0, temp0)
    salt = np.tile(salt, (Yrho, Xrho, 1, 1))
    salt = salt.transpose((2, 3, 0, 1))
    temp = np.tile(temp, (Yrho, Xrho, 1, 1))
    temp = temp.transpose((2, 3, 0, 1))

    # construct initial netCDF file
    class ocean_time_info(object):
        pass
    ocean_time = ocean_time_info()
    ocean_time.long_name = 'seconds since 00-00-00'
    ocean_time.units = 'second'
    pyroms_toolbox.nc_create_roms_file(ini_filename, grd, ocean_time, geogrid=False)
    fh = nc.Dataset(ini_filename, 'r+')
    fh.variables['ocean_time'][:] = 0

    fh.createVariable('zeta', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
    fh.variables['zeta'].long_name = 'free-surface'
    fh.variables['zeta'].units = 'meter'
    fh.variables['zeta'].field = 'free-surface, scalar, series'
    fh.variables['zeta'][:] = 0

    fh.createVariable('temp', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
    fh.variables['temp'].long_name = 'potential temperature'
    fh.variables['temp'].units = 'Celsius'
    fh.variables['temp'].field = 'temperature, scalar, series'
    fh.variables['temp'][:] = temp

    fh.createVariable('salt', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
    fh.variables['salt'].long_name = 'salinity'
    fh.variables['salt'].units = 'PSU'
    fh.variables['salt'].field = 'salinity, scalar, series'
    fh.variables['salt'][:] = salt

    fh.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
    fh.variables['u'].long_name = '3D u-momentum component'
    fh.variables['u'].units = 'meter second-1'
    fh.variables['u'].field = 'u-velocity, scalar, series'
    fh.variables['u'][:] = 0

    fh.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
    fh.variables['v'].long_name = '3D v-momentum component'
    fh.variables['v'].units = 'meter second-1'
    fh.variables['v'].field = 'v-velocity, scalar, series'
    fh.variables['v'][:] = 0

    fh.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
    fh.variables['ubar'].long_name = '2D u-momentum component'
    fh.variables['ubar'].units = 'meter second-1'
    fh.variables['ubar'].field = 'ubar-velocity, scalar, series'
    fh.variables['ubar'][:] = 0

    fh.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
    fh.variables['vbar'].long_name = '2D v-momentum component'
    fh.variables['vbar'].units = 'meter second-1'
    fh.variables['vbar'].field = 'vbar-velocity, scalar, series'
    fh.variables['vbar'][:] = 0

    fh.close()

# ----------------------------------------------------------------------
# Make boundary file

# boundary file name
bry_filename = 'roms_bry.nc'

# read in CTD data
ctddata = loadmat('KS2014_KS1_TS.mat')
salt0 = np.nanmean(ctddata['july2829'][0, 0][0], axis=1)
temp0 = np.nanmean(ctddata['july2829'][0, 0][1], axis=1)
depth0 = -1*ctddata['july2829'][0, 0][2].squeeze()
depth = grd.vgrid.z_r[:][:, 0, 0]
spval = -1.0e20

salt0 = salt0[::-1]
temp0 = temp0[::-1]
depth0 = depth0[::-1]
salt0[np.isnan(salt0)] = salt0[498]
temp0[np.isnan(temp0)] = temp0[498]

salt = np.interp(depth, depth0, salt0)
temp = np.interp(depth, depth0, temp0)
salt = np.tile(salt, (Yrho, Xrho, 1, 1))
salt = salt.transpose((2, 3, 0, 1))
temp = np.tile(temp, (Yrho, Xrho, 1, 1))
temp = temp.transpose((2, 3, 0, 1))

# construct boundary netCDF file
nc_ocean_time = np.array([0, 2000*24*60*60])
class ocean_time_info(object):
    pass
ocean_time = ocean_time_info()
ocean_time.long_name = 'seconds since 00-00-00'
ocean_time.units = 'second'
pyroms_toolbox.nc_create_roms_file(bry_filename, grd, ocean_time, geogrid=False)
fh = nc.Dataset(bry_filename, 'r+')
fh.variables['ocean_time'][:] = nc_ocean_time

# North
fh.createVariable('temp_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['temp_north'].long_name = 'potential temperature north boundary condition'
fh.variables['temp_north'].units = 'Celsius'
fh.variables['temp_north'].field = 'temp_north, scalar, series'
fh.variables['temp_north'][0] = temp[:, :, -1, :]
fh.variables['temp_north'][1] = temp[:, :, -1, :]

fh.createVariable('salt_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['salt_north'].long_name = 'salinity north boundary condition'
fh.variables['salt_north'].units = 'PSU'
fh.variables['salt_north'].field = 'salt_north, scalar, series'
fh.variables['salt_north'][0] = salt[:, :, -1, :]
fh.variables['salt_north'][1] = salt[:, :, -1, :]

fh.createVariable('dye_north_01', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['dye_north_01'].long_name = ' dye 01 north boundary condition'
fh.variables['dye_north_01'].units = ' '
fh.variables['dye_north_01'].field = 'dye_north_01, scalar, series'
fh.variables['dye_north_01'][:] = 0

fh.createVariable('dye_north_02', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['dye_north_02'].long_name = ' dye 02 north boundary condition'
fh.variables['dye_north_02'].units = ' '
fh.variables['dye_north_02'].field = 'dye_north_02, scalar, series'
fh.variables['dye_north_02'][:] = 0

fh.createVariable('dye_north_03', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['dye_north_03'].long_name = ' dye 03 north boundary condition'
fh.variables['dye_north_03'].units = ' '
fh.variables['dye_north_03'].field = 'dye_north_03, scalar, series'
fh.variables['dye_north_03'][:] = 0

# South
fh.createVariable('temp_south', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['temp_south'].long_name = 'potential temperature south boundary condition'
fh.variables['temp_south'].units = 'Celsius'
fh.variables['temp_south'].field = 'temp_south, scalar, series'
fh.variables['temp_south'][0] = temp[:, :, 0, :]
fh.variables['temp_south'][1] = temp[:, :, 0, :]

fh.createVariable('salt_south', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['salt_south'].long_name = 'salinity south boundary condition'
fh.variables['salt_south'].units = 'PSU'
fh.variables['salt_south'].field = 'salt_south, scalar, series'
fh.variables['salt_south'][0] = salt[:, :, 0, :]
fh.variables['salt_south'][1] = salt[:, :, 0, :]

fh.createVariable('dye_south_01', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['dye_south_01'].long_name = 'dye 01 south boundary condition'
fh.variables['dye_south_01'].units = ' '
fh.variables['dye_south_01'].field = 'dye_south_01, scalar, series'
fh.variables['dye_south_01'][:] = 0

fh.createVariable('dye_south_02', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['dye_south_02'].long_name = 'dye 02 south boundary condition'
fh.variables['dye_south_02'].units = ' '
fh.variables['dye_south_02'].field = 'dye_south_02, scalar, series'
fh.variables['dye_south_02'][:] = 0

fh.createVariable('dye_south_03', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
fh.variables['dye_south_03'].long_name = 'dye 03 south boundary condition'
fh.variables['dye_south_03'].units = ' '
fh.variables['dye_south_03'].field = 'dye_south_03, scalar, series'
fh.variables['dye_south_03'][:] = 0

# East
fh.createVariable('temp_east', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
fh.variables['temp_east'].long_name = 'potential temperature east boundary condition'
fh.variables['temp_east'].units = 'Celsius'
fh.variables['temp_east'].field = 'temp_east, scalar, series'
fh.variables['temp_east'][0] = temp[:, :, :, -1]
fh.variables['temp_east'][1] = temp[:, :, :, -1]

fh.createVariable('salt_east', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
fh.variables['salt_east'].long_name = 'salinity east boundary condition'
fh.variables['salt_east'].units = 'PSU'
fh.variables['salt_east'].field = 'salt_east, scalar, series'
fh.variables['salt_east'][0] = salt[:, :, :, -1]
fh.variables['salt_east'][1] = salt[:, :, :, -1]

fh.createVariable('dye_east_01', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
fh.variables['dye_east_01'].long_name = 'dye 01 east boundary condition'
fh.variables['dye_east_01'].units = ' '
fh.variables['dye_east_01'].field = 'dye_east_01, scalar, series'
fh.variables['dye_east_01'][:] = 0

fh.createVariable('dye_east_02', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
fh.variables['dye_east_02'].long_name = 'dye 02 east boundary condition'
fh.variables['dye_east_02'].units = ' '
fh.variables['dye_east_02'].field = 'dye_east_02, scalar, series'
fh.variables['dye_east_02'][:] = 0

fh.createVariable('dye_east_03', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
fh.variables['dye_east_03'].long_name = 'dye 03 east boundary condition'
fh.variables['dye_east_03'].units = ' '
fh.variables['dye_east_03'].field = 'dye_east_03, scalar, series'
fh.variables['dye_east_03'][:] = 0

fh.close()
