import numpy as np
import disc2radmc 
import matplotlib.pyplot as plt


def sigma(r, phi, r0, sigr):

    return np.exp(- 0.5* ((r-r0)/sigr)**2.)*(1.+np.sin(phi))

########################
###### PARAMETERS ######
########################

rmin=30.
rmax=300.
Mdust=0.5 # Mearth
par_sigma=(100.,20.)
Nr=30
Nphi=100
Nth=10
axisym=False # it can overwrite Nphi if True and set it to 1
mirror=True
dpc=50.
target='HDx'
Npix=512
dpix=0.05
inc=0.
PA=90.
wavelength=880. # um
fstar=-1.0
scattering_mode=0
h=0.05
alpha_dust=2.5
new_lambda=1300. # um

bg1=[Npix, dpix, 5.0e-3, -8., 8., 1., 0.5, 45.]
bg2=[Npix, dpix, 5.0e-3, -4., 4., 1., 0.5, 45.]
bg3=[Npix, dpix, 5.0e-3, -2., 2., 1., 0.5, 45.]


################################
###### CREATE INPUT FILES ######
################################

### PHYSICAL GRID
gridmodel=disc2radmc.model.physical_grid(rmin=rmin, rmax=rmax, Nr=Nr, Nphi=Nphi, Nth=Nth, mirror=mirror, logr=True, axisym=axisym)
gridmodel.save()

### WAVELENGTH GRID
lammodel=disc2radmc.model.wavelength_grid()
lammodel.save()

### STAR
starmodel=disc2radmc.model.star(lammodel, Tstar=5800., Rstar=1.0)
starmodel.save()

### DUST SIZE DISTRIBUTION AND OPACITY

path='/Users/Sebamarino/Astronomy/Opacities/Optical_constants/'
lnk_files=[path+'astrosilicate_ext.lnk' ,
           path+'ac_opct.lnk']#,
           # path+'ice_opct.lnk']
ddist=disc2radmc.model.dust_size_distribution(lammodel, Mdust=Mdust, lnk_file=lnk_files, densities=[4., 3.], N_species=1, N_per_bin=100, mass_weights=[100.0, 0.], tag='mix', compute_opct=True)
opac=ddist.compute_opacities()

### DUST DENSITY DISTRIBUTION
ddens=disc2radmc.model.dust_densities(grid=gridmodel, ddist=ddist, function_sigma=sigma, par_sigma=par_sigma, h=h)
ddens.write_density()

################################
###### RUN RADMC3D ###### ######
################################

sim=disc2radmc.model.simulation(scattering_mode=scattering_mode)

### RUN MCTHERM
sim.mctherm()

### MAKE IMAGE
sim.simimage(dpc=dpc, imagename='test', wavelength=880., Npix=Npix, dpix=dpix, inc=inc, PA=PA, tag=target, background_args=[bg1, bg2, bg3])

### MAKE 2nd image using a fixed disc spectral index
path_image='image_test_'+target+'.out'
path_fits='images/image_test_alpha_'+target+'.fits'
disc2radmc.model.convert_to_fits(path_image, path_fits, Npix, dpc,  background_args=[bg1, bg2, bg3], tag=target, alpha_dust=alpha_dust, new_lambda=new_lambda)


