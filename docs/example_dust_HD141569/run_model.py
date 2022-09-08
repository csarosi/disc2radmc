import numpy as np
import disc2radmc 

# import astropy.io.fits as pyfits
# import matplotlib.pyplot as plt
# import json
# from scipy import interpolate



#################################################################################################
####### SURFACE DENSITY DEFINITION ###############################################################
## Here you define the surface density distribution as a function of r and phi (azimuthal angle).
## It doesn't need to be normalized as the code will scale it to match the total dust mass specified.
#################################################################################################

def Sigma_dust(r, phi,  Ms, rcs, sigrs):
    ## Multiple radial Gaussians with masses Ms, central radii rcs, and standard deviations sigrs
    
    Ncomponents=len(Ms)
    S=np.zeros_like(r)
    for i in range(Ncomponents):
        S+=Ms[i]*np.exp( -0.5 * ((r-rcs[i])/sigrs[i])**2.)/(sigrs[i]*rcs[i]) # normalized so the relative mass of each ring is right
    return S



##############################
###### MODEL PARAMETERS ######
##############################

dpc=111.6  # [pc] distance to source 
target='HD141569'# name


## STAR (set to match HD141569) 
## not used as the file containing the stellar spectrum and location (stars.inp) has already been created
# Rstar=1.8 # [Solar radii]
# Tstar=8500 # [K] If this is negative, the code will consider the star as a blackbody. 
# g=4.0
# # For a realistic stellar model, you can download bt-settl models 
# # from http://svo2.cab.inta-csic.es/theory/newov2/index.php and indicate their directory below
# dir_stellar_templates='/Users/Sebamarino/Astronomy/Stellar_templates/BT-Settl/bt-settl/'
# # The code will then search for files named as 
# # dir_stellar_templates+'lte%03i-%1.1f-0.0a+0.0.BT-NextGen.7.dat.txt'%(Tstar//100,g) if Tstar>=2600
# # dir_stellar_templates+'lte%03i-%1.1f-0.0.BT-Settl.7.dat.txt'%(T//100,g)            else
# # and interpolate models with neighbouring temperatures if necessary

## DUST parameters
Mdusts=[5.0e-5, 1.0, 0.5, 0.5] # [Mearth] mass of each ring component
rcs=[0.5, 50., 200., 350.]     # [au] center of each Gaussian ring
sigrs=[0.2, 10., 10., 20.]     # [au] standard deviation of each ring
h=0.05                         # vertical aspect ratio =H/r, where H is the vertical standard deviation. This is a constant, but other parametrizations are possible. 
par_sigma=(Mdusts, rcs, sigrs) # list containing the parameters that define the dust surface density. They must be in the same order as in the definition of Sigma_dust
Mdust=np.sum(Mdusts) # total dust mass in Mearth


# grain size distribution parameters
amin=10.0      # [micron]  minimum grain size (1 by default)
amax=1.0e4     # [micron]  maximum grain size (1e4 by default)
N_species=3    #  (1 by default) Number of dust size bins to use for the radiative transfer calculations. 
slope = -3.5   #  (-3.5 by default) Slope of the size distribution. This is used for computing opacities and for the mass distribution of each size bin

## MODEL SPATIAL GRID
rmin=0.1 # [au] make sure it is small enough to sample the surface density
rmax=500.# [au] make sure it is large enough to sample the surface density
Nr=200    # number of radial cells (linearly or logspaced)
Nphi=100  # number of azimuthal cells
Nth=1     # number of polar angle cells (for only one emisphere)
thmax=np.arctan(h)*1 # maximum polar angle to sample as measured from the midplane.
axisym=False # Consider the disc to be axisymmetric to speed up calculations? it can overwrite Nphi if True and set it to 1
mirror=False  # Mirror the upper half to speed up calculations. This is incompatible with anisotropic scattering.
logr=True # Sample r logarithmically or linearly

# WAVELENGTH GRID (grid to sample the stellar flux in temperature calculations, see radmc3d manual)
lammin=0.09  # [mu]  minimum wavelength to consider in our model (important for temperature calculation)
lammax=1.0e5 # [mu] minimum wavelength to consider in our model (important for temperature calculation)
Nlam=150     # number of cells logarithmically spaced to sample the wavelength range.

# IMAGE PARAMETERS
Npix=400  # number of pixels
dpix=0.027491 # pixel size in arcsec
PA=-8.0    # [deg] position angle
inc=-55.0  # [deg] inclination (sign matters to define near or far side)

wavelength=[10.575, 11.3, 15.50] # [um] image wavelength. It can be a float or an array of wavelengths
scattering_mode=2 # scattering mode (0=no scattering, 1=isotropic, 2=anisotropic using H&G function)



#######################################################
#### SETUP THE MODEL USING DISC2RADMC #################
#######################################################

### PHYSICAL GRID (this only needs to be run once at the beginning to define the model)
gridmodel=disc2radmc.physical_grid(rmin=rmin, rmax=rmax, Nr=Nr, Nphi=Nphi, Nth=Nth, thmax=thmax, mirror=mirror, logr=logr, axisym=axisym)
gridmodel.save()

### WAVELENGTH GRID (this only needs to be run once at the beginning)
lammodel=disc2radmc.wavelength_grid(lammin=lammin, lammax=lammax, Nlam=Nlam)
lammodel.save()

### STAR (this only needs to be run once if you don't have the stars.inp file)
# starmodel=disc2radmc.star(lammodel,
#                           Tstar=Tstar,
#                           Rstar=Rstar,
#                           g=g,
#                           dir_stellar_templates=dir_stellar_templates # necessary line if Tstar>0
#                         )
# starmodel.save()


### DUST SIZE DISTRIBUTION AND OPACITY (this only needs to be run once at the beginning)
# path to optical constants that can be found at
# https://github.com/SebaMarino/disc2radmc/tree/main/opacities/dust_optical_constants
path_opct='/Users/Sebamarino/Astronomy/Codes/disc2radmc/opacities/dust_optical_constants/' 
lnk_files=[path_opct+'astrosilicate_ext.lnk',
           path_opct+'ac_opct.lnk',
           path_opct+'ice_opct.lnk']
densities=[4., 3., 1.] # densities in g/cm3
mass_weights=[1.0, 0., 0.] # mixing ratios by mass
dust=disc2radmc.dust(lammodel,
                     Mdust=Mdust,
                    lnk_file=lnk_files,
                    densities=densities,
                    N_species=N_species,
                    slope=slope,
                    N_per_bin=100, # number of species per size bin to have a good representation 
                    mass_weights=mass_weights,
                    tag='mix', # name to give to this new mixed species
                    compute_opct=False) # only set to true if you are changing the dust composition

# Compute dust opaicties with Mie theory only if you have changed the size distribution parameters or optical constants.
# It requires having the file makeopac in the same working directory. 
# makeopac can be found and compiled from files at https://github.com/SebaMarino/disc2radmc/tree/main/opacities/Mie
# dust.compute_opacities()

### DUST DENSITY DISTRIBUTION (this needs to be every time the spatial distribution or total masses are changed)
dust.dust_densities(grid=gridmodel,function_sigma=Sigma_dust, par_sigma=par_sigma, h=h)
dust.write_density()

## SET SOME BASIC PARAMETERS TO TELL RADMC3D HOW MANY PHOTONS AND WHAT KIND OF SCATTERING WE WANT TO CONSIDER (only necessary to run the first time)
sim=disc2radmc.simulation(nphot=10000000, # number of photon packages for thermal monte carlo
                        nphot_scat=1000000, # number of photon packages for image
                        nphot_spec=10000,   # number of photon packages for spectrum
                        scattering_mode=scattering_mode, 
                        modified_random_walk=0, # for very optically thick medium, this is a useful approximation
                        istar_sphere=0, # consider the star a sphere or a point.
                        setthreads=4,
                        verbose=True, 
                            )
### RUN MCTHERM to compute temperature (only necessary to run once, or every time the grid has changed or if the dust mass is high enough to affect the temperature field, i.e. optically thick)
sim.mctherm()

### MAKE IMAGE
sim.simimage(dpc=dpc, 
             imagename='MIRI_FQPM', 
             wavelength=wavelength,
             Npix=Npix,
             dpix=dpix,
             inc=inc,
             PA=PA,
             tag=target)
