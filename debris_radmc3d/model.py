import numpy as np
import cmath as cma
import os,sys
from debris_radmc3d.constants import *
from astropy.io import fits
from astropy.io.votable import parse


class simulation:
    """
    A class to run radmc3d and convert output files 
    """

    def __init__(self,  nphot=1000000, nphot_scat=1000000, nphot_spec=10000, nphot_mono=10000, scattering_mode=1, modified_random_walk=0, istar_sphere=0, tgas_eq_tdust=1, incl_lines=0, setthreads=4, rto_style=3, verbose=True):

        self.nphot=nphot
        self.nphot_scat=nphot_scat
        self.nphot_spec=nphot_spec
        self.nphot_mono=nphot_mono
        self.scattering_mode=scattering_mode
        self.modified_random_walk=modified_random_walk
        self.istar_sphere=istar_sphere
        self.tgas_eq_tdust=tgas_eq_tdust
        self.incl_lines=incl_lines
        self.setthreads=setthreads
        self.rto_style=rto_style
        self.verbose=verbose
        
        radmc_file=open('radmc3d.inp','w')
        radmc_file.write('nphot =       %1.0f \n'%self.nphot)
        radmc_file.write('nphot_scat=    %1.0f \n'%self.nphot_scat)
        radmc_file.write('nphot_spec=    %1.0f \n'%self.nphot_spec)
        radmc_file.write('nphot_mono=    %1.0f \n'%self.nphot_mono)
        radmc_file.write('scattering_mode_max = %1.0f \n'%self.scattering_mode)
        radmc_file.write('modified_random_walk = %1.0f \n'%self.modified_random_walk)
        radmc_file.write('istar_sphere= %1.0f \n'%self.istar_sphere)
        radmc_file.write('tgas_eq_tdust=%1.0f \n'%self.tgas_eq_tdust)
        radmc_file.write('incl_lines = %1.0f \n'%self.incl_lines)
        radmc_file.write('setthreads = %1.0f \n'%self.setthreads)
        radmc_file.write('rto_style = %1.0f'%self.rto_style )
        radmc_file.close()

        if not os.path.exists('./images'):
            os.makedirs('./images')

    ### methods
    # def mctherm
    # def make_image(self):
    # def make_cube(self):
    # def make_sed(self):

    def mctherm(self):
        if self.verbose:
            os.system('radmc3d mctherm')
        else:
            os.system('radmc3d mctherm > mctherm.log')

    def simimage(self, dpc=1., imagename='', wavelength=880., Npix=256, dpix=0.05, inc=0., PA=0., offx=0.0, offy=0.0, X0=0., Y0=0., tag='', omega=0.0, Npixf=-1, fstar=-1.0, background_args=[]):
        # X0, Y0, stellar position (e.g. useful if using a mosaic)
        # images: array of names for images produced at wavelengths
        # wavelgnths: wavelengths at which to produce image in um
        # fields: fields where to make images (=[0] unless observations are a mosaic)

        if Npixf==-1:
            Npixf=Npix

        sau=Npix*dpix*dpc
        print('image size = %f'%sau)
        if self.verbose:
            os.system('radmc3d image incl '+str(inc)+' phi '+str(omega)+' posang '+str(PA-90.0)+'  npix '+str(Npix)+' lambda '+str(wavelength)+' sizeau '+str(sau)+' secondorder')
        else:
            os.system('radmc3d image incl '+str(inc)+' phi '+str(omega)+' posang '+str(PA-90.0)+'  npix '+str(Npix)+' lambda '+str(wavelength)+' sizeau '+str(sau)+' secondorder  > simimgaes.log')

        pathin ='image_'+imagename+'_'+tag+'.out'
        pathout='image_'+imagename+'_'+tag+'.fits'
        os.system('mv image.out '+pathin)
        
        convert_to_fits(pathin, pathout, Npixf, dpc, mx=offx, my=offy, x0=X0, y0=Y0, omega=omega,  fstar=fstar, background_args=background_args, tag=tag)
        os.system('mv '+pathout+' ./images')

    
    
    

    
    
   

# class gas: # needs a quick way to redefine and save
#     """
#     A class used to define the gas densities and velocities
#     """
   

class dust_densities:
    
    """
    A class used to define the dust densities.
    
    """
    
    def __init__(self, grid=None, ddist=None, function_sigma=None, par_sigma=None, h=0.05):

        """
        funtion_sigma: function that defines the surface density. The first two arguments are two nd numpy arrays for rho and z
        """
        
        assert grid is not None, "grid object needed to define dust density distribution"
        assert ddist is not None, "dust size distribution needed to define dust density distribution"
        assert function_sigma is not None, "surface density profile needed to define dust density distribution"

        ## strange cases ntheta=1, mirror, etc
        self.grid=grid
        self.ddist=ddist
        self.rho_d=np.zeros((ddist.N_species,grid.Nth,grid.Nphi,grid.Nr)) # density field (only norther emisphere)

        thetam, phim, rm=np.meshgrid(grid.th, grid.phi, grid.r, indexing='ij' ) # so it goes from Northpole to equator. theta is still the angle from the equator.


        # if grid.Nphi>1:
        #     dThm, dPhim, dRm = np.meshgrid(Thedge[1:]-Thedge[:-1], Phiedge[1:]-Phiedge[:-1], Redge[1:]-Redge[:-1], indexing='ij' )
        # else:
        dthm, dphim, drm = np.meshgrid(grid.dth, grid.dphi, grid.dr, indexing='ij' )

        rhom=rm*np.cos(thetam) 
        zm=rm*np.sin(thetam)

        # Ms=ddist.Mgrid
        
        for ia in range(ddist.N_species):
            M_dust_temp= 0.0
        
            # nother emisphere
            if grid.Nth>1: # more than one cell per emisphere
                self.rho_d[ia,:,:,:]=self.rho_3d_dens(rhom, phim, zm, h, function_sigma, *par_sigma)

                # # now south emisphere is copy of nother emisphere
                # rho_d[ia,Nth-1:,:,:]=rho_d[ia,-Nth::-1,:,:]
                ### the line above works because this is how step and slicing work https://stackoverflow.com/questions/509211/understanding-slice-notation
                # a[::-1]    # all items in the array, reversed
                # a[1::-1]   # the first two items, reversed
                # a[:-3:-1]  # the last two items, reversed
                # a[-3::-1]  # everything except the last two items, reversed
            
            
            elif grid.Nth==1:# one cell

                self.rho_d[ia,:,:,:]=function_sigma(rhom, phim, *par_sigma)/(grid.dth[0]*rhom) # rho_3d_dens(rho, 0.0, 0.0, hs, sigmaf, *args )
                # rho_d[ia,1,:,:]=rho_d[ia,0,:,:]

            M_dust_temp=2.*np.sum(self.rho_d[ia,:,:,:]*(dphim*rhom)*(drm)*(dthm*rm))*au**3.0
            self.rho_d[ia,:,:,:]=self.rho_d[ia,:,:,:]*ddist.Mgrid[ia]/M_dust_temp*M_earth


    ###############
    ### methods ###
    ###############
    @staticmethod
    def rho_3d_dens(rho, phi, z, h, function_sigma, *arguments ):

        H=h*rho # au
        return function_sigma(rho,phi, *arguments)*np.exp(-z**2.0/(2.0*H**2.0))/(np.sqrt(2.0*np.pi)*H)
        
        
    def write_density(self):

        # Save 
        path='dust_density.inp'
        file_dust=open(path,'w')
    
        file_dust.write('1 \n') # iformat
        if self.grid.mirror:
            file_dust.write(str((self.grid.Nr)*(self.grid.Nth)*(self.grid.Nphi))+' \n') # iformat n cells
        else:
            file_dust.write(str((self.grid.Nr)*(2*self.grid.Nth)*(self.grid.Nphi))+' \n') # iformat n cells

        file_dust.write(str(self.ddist.N_species)+' \n') # n species

        for ai in range(self.ddist.N_species):
            for j in range(self.grid.Nphi):
                if self.grid.mirror:
                    for k in range(self.grid.Nth):
                        for i in range(self.grid.Nr):
                            file_dust.write(str(self.rho_d[ai,-(1+k),j,i])+' \n')
                else:
                    # northern emisphere
                    for k in range(self.grid.Nth):
                        for i in range(self.grid.Nr):
                            file_dust.write(str(self.rho_d[ai,-(1+k),j,i])+' \n')
                    # southern emisphere
                    for k in range(self.grid.Nth):
                        for i in range(self.grid.Nr):
                            file_dust.write(str(self.rho_d[ai,k,j,i])+' \n')
        file_dust.close()
        

    
class dust_size_distribution:
    """
    A class used to define the dust distribution and opacities
    """
    def __init__(self, wavelength_grid, Mdust=0.1, lnk_file=None,amin=1.0, amax=1.0e4, slope=-3.5, density=3.0, N_species=1, N_per_bin=50, densities=None, mass_weights=None, tag='i', compute_opct=True ):
        """
        Mdust: dust mass in earth masses
        amin: minimum grain size in um
        amax: maximum grain size in um
        """

        self.lnk_file=lnk_file if lnk_file is not None else sys.exit('invalid lnk_file')
        self.wavelength_grid=wavelength_grid
        self.amin=amin
        self.amax=amax
        self.slope=slope
        self.density=density
        self.N_species=N_species
        self.N_per_bin=N_per_bin
        self.tag=tag

        ### size grid
        self.Agrid_edges=np.logspace(np.log10(self.amin), np.log10(self.amax), self.N_species+1)
        self.Agrid=np.sqrt(self.Agrid_edges[1:]*self.Agrid_edges[:-1])
        self.Mgrid=self.Agrid**(self.slope+4.)
        self.Mgrid=Mdust*self.Mgrid/np.sum(self.Mgrid)
        
        if isinstance(lnk_file, str): # if one optical constant given
            self.lnk_file_p=lnk_file
            self.density=density
            
        elif isinstance(self.lnk_file, list):
            self.densities =np.array(densities) if densities is not None else sys.exit('error in densities array')
            self.mass_weights=np.array(mass_weights) if mass_weights is not None else sys.exit('error in mass weights array')

            if len(lnk_file)==len(mass_weights) and len(lnk_file)==len(densities):
                ### compute average and save
                self.lnk_file_p='opct_'+self.tag+'.lnk'
                if compute_opct:
                    print('Compute average optical constants')
                    Opct=self.mix_opct_bruggeman(pathout='opct_'+self.tag+'.lnk')
            else:
                sys.exit('mass_weights or densities do not have right length')
        else:
            sys.exit('lnk_file does not have the right format (list or string)')

    ###############
    ### methods ###
    ###############
            
    ### compute opacities
    def compute_opacities(self):
        path="./"
        if not os.path.exists('./Tempkappa'):
            os.makedirs('./Tempkappa')

        for j in range(self.N_species):

            Agrid_i=np.logspace(np.log10(self.Agrid_edges[j]), np.log10(self.Agrid_edges[j+1]), self.N_per_bin)
            weights=(Agrid_i)**(self.slope+4.) # w(a) propto n(a)*m(a)*da and da propto a
            weights=weights/np.sum(weights)

            os.system('rm '+path+'Tempkappa/*')

            opctj=np.zeros((self.wavelength_grid.Nlam,3))
            for i in range(self.N_per_bin):

                os.system('rm '+path+'param.inp')
                acm=Agrid_i[i]*1.0e-4
            
                file_inp=open(path+'param.inp','w')
                file_inp.write(self.lnk_file_p[:-4]+'\n')
                e=round(np.log10(acm))
                b=acm/(10.0**e)
                file_inp.write('%1.5fd%i \n' %(b,e))
                file_inp.write('%1.5f \n' %self.density) 
                file_inp.write('1')
                file_inp.close()
                
                os.system(path+'makeopac')
                os.system('mv '+path+'dustkappa_'+self.lnk_file_p[:-4]+'.inp '+path+'Tempkappa/dustkappa_temp_'+str(i+1)+'.inp ')

                # read opacity
                opcti=np.loadtxt('./Tempkappa/dustkappa_temp_'+str(i+1)+'.inp', skiprows=2)
                opctj=opctj+weights[i]*opcti[:,1:]

            ### write opacity file
            pathout='dustkappa_'+self.tag+'_'+str(j+1)+'.inp'
            file_opacity=open(pathout,'w')
            file_opacity.write('3 \n')
            file_opacity.write(str(self.wavelength_grid.Nlam)+'\n')
            for i in range(self.wavelength_grid.Nlam):
                file_opacity.write('%f \t %f \t %f \t %f\n' %(self.wavelength_grid.lams[i],opctj[i,0],opctj[i,1],opctj[i,2]))
            file_opacity.close()

        os.system('rm ./Tempkappa/*')
        path='dustopac.inp'
        file_list_opacities=open(path,'w')
        file_list_opacities.write("2               Format number of this file \n")
        file_list_opacities.write(str(self.N_species)+"              Nr of dust species \n")
        file_list_opacities.write("============================================================================ \n")
        for i in range(self.N_species):
            file_list_opacities.write("1               Way in which this dust species is read \n")
            file_list_opacities.write("0               0=Thermal grain \n")
            file_list_opacities.write(self.tag+"_"+str(i+1)+ " Extension of name of dustkappa_***.inp file \n")
            file_list_opacities.write("---------------------------------------------------------------------------- \n")
        file_list_opacities.close()
            
        
    def mix_opct_bruggeman(self, pathout='opct_mix.lnk'):

        # Mixing rule Bruggeman for max 3 species

        N_opct=len(self.lnk_file)
       

        self.volumes=self.mass_weights/self.densities
       
        self.density=np.sum(self.mass_weights)/np.sum(self.volumes) # final density

        self.volume_weights=self.volumes/np.sum(self.volumes)
      
        self.voli=self.volume_weights[1:]/self.volume_weights[0]
        
        print("final density = %1.1f g/cm3"%self.density)


        O1=np.loadtxt(self.lnk_file[0] )
        O2=np.loadtxt(self.lnk_file[1] )
   
        if N_opct==3:
            O3=np.loadtxt(self.lnk_file[2] )
            
        Opct1=np.zeros((self.wavelength_grid.Nlam,3))
        Opct1[:,0]=self.wavelength_grid.lams

        for i in range(self.wavelength_grid.Nlam):

            n1=Intextpol(O1[:,0],O1[:,1],Opct1[i,0])
            n2=Intextpol(O2[:,0],O2[:,1],Opct1[i,0])

            k1=Intextpol(O1[:,0],O1[:,2],Opct1[i,0])
            k2=Intextpol(O2[:,0],O2[:,2],Opct1[i,0])

            if N_opct==3:
                n3=Intextpol(O3[:,0],O3[:,1],Opct1[i,0])
                k3=Intextpol(O3[:,0],O3[:,2],Opct1[i,0])
                
                eff=effnk(n1,k1,n2,k2,n3,k3,self.voli[0],self.voli[1])
            else:
                eff=effnk(n1,k1,n2,k2,0.,0.,self.voli[0],0.)

                
            Opct1[i,1]=eff.real
            Opct1[i,2]=eff.imag

        np.savetxt(pathout,Opct1)



class star:
    """
    A class to define the star and companions
    """
    def __init__(self, lam_grid,
                 dir_stellar_templates='/Users/Sebamarino/Astronomy/Stellar_templates/BT-COND/bt-settl/',
                 Tstar=None,
                 Rstar=None,
                 Mstar=None,
                 g=4.
                 # companion=False
                 # separation=0.0,
                 # inc=0.0,
                 # plt=False, T_plt=0.0, R_plt=0.0, M_plt=0.0
                 ):

        default_Tstar=5750.
        default_Rstar=1.0
        default_Mstar=1.0

        if Tstar is not None:
            if (Tstar>=2000. and Tstar<=10000.):
                self.Tstar=Tstar
            elif Tstar<0.0:
                self.Tstar=Tstar
        else: self.Tstar= default_Tstar
           
        if Rstar is not None:
            self.Rstar=Rstar if Rstar>0. else default_Rstar
        else:
            self.Rstar=default_Rstar
        if Mstar is not None:
            self.Mstar=Mstar if Mstar>0. else default_Mstar
        else:
            self.Mstar=default_Mstar
            
        self.Nlam=lam_grid.Nlam
        self.lams=lam_grid.lams
        self.dlams=lam_grid.dlams
        self.model_directory=dir_stellar_templates
        #########################
        ###### load spectra
        #########################
        # Temperature resolution of NextGen grid

        if self.Tstar>0.0:

            if self.Tstar<7000.: 
                dT=100.
            else:
                dT=200.
            
            if self.Tstar%dT==0.: # no need to interpolate
                spectrum=self.get_spectrum(self.Tstar)
            else:
              
                    
                T1=dT*self.Tstar//dT
                T2=T1+dT

                spectrum1=self.get_spectrum(T1)
                spectrum2=self.get_spectrum(T2)

                w1=abs(T1-self.Tstar)
                w2=abs(T2-self.Tstar)
                spectrum=spectrum1*w1+spectrum2*w2/(w1+w2)

        ### compute flux at 1pc
        # spectrum is in units of erg/cm2/s/A
        # we want to convert Flambda to Fnu so we need to multiply by lambda**2/c
        spectrum_fnu=spectrum*(self.lams*1.0e4)**2.0/cc_a # now in units of erg/cm2/s/Hz

        ## flux at 1pc
        self.flux_1pc=spectrum_fnu*(self.Rstar*R_sun/pc)**2.

        


        
    def get_spectrum(self, T, g=4.):
        """
        returns model spectrum in units of erg/cm2/s/A 
        """

        path=self.model_directory+'lte%03i-%1.1f-0.0a+0.0.BT-NextGen.7.dat.txt'%(T//100,g)
        data=np.loadtxt(path)
        
        Nlam_model=data.shape[0]
    
        # down sample spectrum
        spectrum_lowr=np.zeros(self.Nlam)
        for i in range(self.Nlam):
            if self.lams[i]<1.0e3:
                mask=(data[:,0]/1.0e4>(self.lams[i]-self.dlams[i]/2.)) & (data[:,0]/1.0e4<(self.lams[i]+self.dlams[i]/2.)) # data in Armstrong and lams in um
                spectrum_lowr[i]=np.mean(data[mask,1])
            else: # Rayleigh Jeans extrapolation
                spectrum_lowr[i]=spectrum_lowr[i-1]*(self.lams[i]/self.lams[i-1])**(-4.)
        return spectrum_lowr
        
       
        
    def save(self):

        path='stars.inp'
        file_star=open(path,'w')
        file_star.write('2 \n')

        # if plt:
        # file_star.write('2  '+str(Nw)+'\n')
        # else:

        file_star.write('1  '+str(self.Nlam)+'\n')
        file_star.write(str(self.Rstar*R_sun)+'\t'+str(self.Mstar*M_sun)+'\t 0.0  0.0  0.0  \n')
        # if plt:
        #     file_star.write(str(R_plt)+'\t'+str(M_plt)+'\t'+'0.0 0.0 0.0   \n')

        ### wavelengths
        for i in range(self.Nlam):
            file_star.write(str(self.lams[i])+'\n')
        ### star
        if self.Tstar>0.0:
            for i in range(self.Nlam):
                file_star.write(str(self.flux_1pc[i])+'\n')
        else:
            file_star.write(str(self.Tstar)+'\n')
        # ### planet
        # if plt:
        #     file_star.write(str(-T_plt)+'\n')

        file_star.close()

class wavelength_grid:
    """
    A class used to define the wavelength grid (in um)
    """

    def __init__(self,lammin=None, lammax=None, Nlam=None):

        if lammin is not None:
            self.lammin=lammin if lammin>0. else 0.09 # um
        else: self.lammin=0.09
        if lammax is not None:
            self.lammax=lammax if lammax>self.lammin else 1.0e5 # um
        else: self.lammax=1.0e5
        if Nlam is not None:
            self.Nlam=int(Nlam) if Nlam>0 else 150
        else: self.Nlam=150

        self.lams=np.logspace(np.log10(self.lammin), np.log10(self.lammax), self.Nlam)
        P=(self.lams[1]-self.lams[0])/self.lams[0]
        self.dlams=P*self.lams
        
    def save(self):
        # ----- write wavelength_micron.inp

        path='wavelength_micron.inp'
        file_lams=open(path,'w')
        file_lams.write(str(self.Nlam)+'\n')
        for i in range(self.Nlam):
            file_lams.write(str(self.lams[i])+'\n')
        file_lams.close()



    
    
class physical_grid:
    """
    A class used to define and store grid parameters

    Theta is defined from midplane (theta=0) increasing towards the N pole (theta=pi/2), but in radmc3d theta=0 is the N pole and theta=90 is midplane.

    """
    

    
    def __init__(self, Nr=None, Nphi=None, Nth=None, rmin = None, rmax=None, thmin=None, thmax=None, logr=False, logtheta=False,  axisym=False, mirror=True, save=True, load=False):


        default_Nr=100
        default_Nphi=100
        default_Nth=20

        default_rmin=1.0
        default_rmax=300.
        default_thmin=0.03/10.
        default_thmax=0.1*5

        if rmin is not None:
            self.rmin=rmin if rmin>0. else default_rmin
        else:  self.rmin=default_rmin
            
        if rmax is not None:
            self.rmax=rmax if rmax>self.rmin else default_rmax
        else: self.rmax=default_rmax

        if thmin is not None:
            self.thmin=thmin if thmin>0. and thmin<np.pi/2.  else default_thmin
        else: self.thmin=default_thmin

        if thmax is not None:
            self.thmax=thmax if thmax>thmin and thmax<=np.pi/2. else default_thmax
        else: self.thmax=default_thmax
        
        if Nr is not None:
            self.Nr=int(Nr) if Nr>0 else default_Nr
        else: self.Nr=default_Nr
        
        if Nth is not None:
            self.Nth=int(Nth) if Nth>0 else default_Nth
        else: self.Nth=default_Nth
        
        self.axisym=axisym
        self.logtheta=logtheta
        self.logr=logr
        self.mirror=mirror
        
        if axisym:
            self.Nphi=1
        else:
            if Nphi is not None:
                self.Nphi=int(Nphi) if Nphi>0 else default_Nphi
            else: self.Nphi=default_Nphi


        ### R
        if self.logr: # log sampling
            self.redge=np.logspace(np.log10(self.rmin), np.log10(self.rmax), self.Nr+1)
        else:
            self.redge=np.linspace(self.rmin, self.rmax, self.Nr+1)
        self.dr=self.redge[1:]-self.redge[:-1]
        self.r=(self.redge[1:]+self.redge[:-1])/2.


        ### Theta 
        if self.logtheta and self.Nth>2: # log sampling

            self.thedge=np.zeros(self.Nth+1)
            self.thedge[0]=0.0
            self.thedge[1:]=np.logspace(np.log10(self.thmin), np.log10(self.thmax), self.Nth)
        else:
            self.thedge=np.linspace(0., self.thmax, self.Nth+1)

        self.dth=self.thedge[1:]-self.thedge[:-1]
        self.th=(self.thedge[1:]+self.thedge[:-1])/2
    
    
        ### Phi

        # ### linear sampling in phi
        self.phiedge=np.linspace(0.,2.0*np.pi, self.Nphi+1)
        self.phi=(self.phiedge[1:]+self.phiedge[:-1])/2.
        self.dphi=self.phiedge[1:]-self.phiedge[:-1]


    def save(self):
    
        path='amr_grid.inp' #'amr_grid.inp'

        gridfile=open(path,'w')
        gridfile.write('1 \n') # iformat: the format number, at present 1
        gridfile.write('0 \n') # Grid style (regular = 0)
        gridfile.write('101 \n') # coordsystem: If 100<=coordsystem<200 the coordinate system is spherical
        gridfile.write('0 \n') # gridinfo

        if self.axisym and self.Nphi==1:
            gridfile.write('1 \t 1 \t 0 \n') # incl x, incl y, incl z
        else:
            gridfile.write('1 \t 1 \t 1 \n') # incl x, incl y, incl z
            
        if self.mirror:
            gridfile.write(str(self.Nr)+ '\t'+ str(self.Nth)+'\t'+ str(self.Nphi)+'\n') 
        else:
            gridfile.write(str(self.Nr)+ '\t'+ str(self.Nth*2)+'\t'+ str(self.Nphi)+'\n') 

        for i in range(self.Nr+1):
            gridfile.write(str(self.redge[i]*au)+'\t')
        gridfile.write('\n')

        for i in range(self.Nth+1):
            gridfile.write(str(np.pi/2.0-self.thedge[self.Nth-i])+'\t')  # from northpole to equator
        if not self.mirror:
            for i in range(1,self.Nth+1):
                gridfile.write(str(np.pi/2.0+self.thedge[i])+'\t')       # from 0 to -pi/2
        gridfile.write('\n')

        for i in range(self.Nphi+1):
            gridfile.write(str(self.phiedge[i])+'\t')
        gridfile.write('\n')
            
        gridfile.close()


    def load(self):
        print('in progress')
            










#### Functions to mix optical constants following Bruggeman's mixing rule

def effnk(n1,k1,n2,k2,n3,k3,f2,f3): 

    # mixing rule Bruggeman http://en.wikipedia.org/wiki/Effective_medium_approximations
    # Sum fi* (epi-ep)/(epi+2ep) = 0, but normilizing by f1

 
    np1=n1+k1*1j # matrix
    np2=n2+k2*1j # inclusion 1
    np3=n3+k3*1j # inclusion 2

    e1=np1**2.0  # n = sqrt(epsilon_r x mu_r) and mu_r is aprox 1
    e2=np2**2.0
    e3=np3**2.0

    # polynomial of third order
    p=np.zeros(4, dtype=complex)

    p[3]=e1*e2*e3*(1.0+f2+f3) # 0 order
    p[2]=-e1*e3*f2 -e1*e2*f3 - e2*e3 + 2*(e1*e2*f2 + e1*e2 +e1*e3*f3+e1*e3+e2*e3*f2+e2*e3*f3) # 1st order
    p[1]= -2.0*(e1*f2+e1*f3+e3*f2+e2*f3+e2+e3)+4.0*(e1+e2*f2+e3*f3)# 2nd order
    p[0]= -4.0*(1.0+f2+f3)

    roots=np.roots(p) 

    # check roots
    for i in range(len(roots)):
        effi=roots[i]
        if effi.real>0.0 and effi.imag>0.0:
            return cma.sqrt(effi)
    ### if nothins satisfy the above condition
        
    return -1.0
        
def Intextpol(x,y,xi):

    Nx=len(x)
    if xi<=x[0]: return y[0] # extrapol                                                                                                                                                                         
    elif xi<=x[Nx-1]: #interpol                                                                                                                                                                                                              
        for l in range(1,Nx):
            if xi<=x[l]:
                return y[l-1]+(xi-x[l-1])*(y[l]-y[l-1])/(x[l]-x[l-1])

    elif xi>x[Nx-1]:    #extrapol                                                                                                                                                                                                            
        alpha=np.log(y[Nx-1]/y[Nx-2])/np.log(x[Nx-1]/x[Nx-2])
        return y[Nx-1]*(xi/x[Nx-1])**alpha


### functions to manipulate images

def convert_to_fits(path_image, path_fits, Npixf, dpc, mx=0.0, my=0.0, x0=0.0, y0=0.0, omega=0.0, fstar=-1.0, vel=False, continuum_subtraction=False, background_args=[], tag=''):

    ### load image
    image_in_jypix, nx, ny, nf, lam, pixdeg_x, pixdeg_y = load_image(path_image, dpc)

    ### manipulate central flux
    istar, jstar=star_pix(nx, omega)
    if fstar>=0.0: # change stellar flux given value of fstar.
         image_in_jypix[:, :, jstar,istar]=fstar
    print('Fstar=', image_in_jypix[0, 0, jstar,istar])

    ### shift image if necessary
    image_in_jypix_shifted= shift_image(image_in_jypix, mx, my, pixdeg_x, pixdeg_y, omega=omega)

    lam0=lam[0] # um   
    reffreq=cc/(lam0*1.0e-4) # Hz

    if nf==1: # single image
        flux = np.sum(image_in_jypix_shifted[0,0,:,:])
        print("flux [Jy] = ", flux)
    else: # image cube
        delta_freq= (lam[0] - lam[1])*cc*1.0e4/lam[nf//2]**2.0 # Hz
        delta_velocity = (lam[1] - lam[0])*cc*1e-5/lam0 # km/s

        if continuum_subtraction: # subtract continuum assuming it varies linearly with wavelength
            m=(image_in_jypix_shifted[0,-1,:,:]- image_in_jypix_shifted[0,0,:,:])/(lam[-1]-lam[0])
            I0=image_in_jypix_shifted[0,0,:,:]*1.
            for k in range(nf):
                Cont=I0+(lam[k]-lam[0])*m
                image_in_jypix_shifted[0,k,:,:]= image_in_jypix_shifted[0,k,:,:] - Cont 
        flux = np.sum(image_in_jypix_shifted[0,:,:,:])*delta_velocity
        print("flux [Jy km/s] = ", flux)
  
        
    ### BUILD HEADER


    # Make FITS header information:
    header = fits.Header()
    #header['SIMPLE']='T'
    header['BITPIX']=-32
    # all the NAXIS are created automatically header['NAXIS']=2
    header['OBJECT']=tag
    header['EPOCH']=2000.0
    # header['LONPOLE']=180.0

    header['EQUINOX']=2000.0
    header['SPECSYS']='LSRK'
    header['RESTFREQ']=reffreq
    header['VELREF']=0.0
    if nf==1:
        header['CTYPE3']='FREQ'
        header['CRPIX3'] = 1.0
        header['CDELT3']  = 1.0
        header['CRVAL3']= reffreq


    header['FLUX']=flux
    header['BTYPE'] = 'Intensity'
    header['BSCALE'] = 1
    header['BZERO'] = 0
    header['BUNIT'] = 'JY/PIXEL'#'erg/s/cm^2/Hz/ster'


    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    if Npixf%2==1: ## if odd number of pixels, central pixel coincides with center
        header['CRVAL1'] = x0
        header['CRVAL2'] = y0
        
    else:   ## if even number of pixels, there is no central pixel
        header['CRVAL1'] = x0+pixdeg_x
        header['CRVAL2'] = y0-pixdeg_y
     
    unit = 'DEG'
    multiplier = 1
    # RA
    header['CDELT1'] = -multiplier*pixdeg_x
    header['CUNIT1'] = unit
    # ...Zero point of coordinate system
    header['CRPIX1'] = (Npixf)//2+1

    # DEC
    header['CDELT2'] = multiplier*pixdeg_y
    header['CUNIT2'] = unit
    # ...Zero point of coordinate system
    header['CRPIX2'] = (Npixf)//2+1

    # FREQ
    if nf > 1:
        if vel==True:
            # multiple frequencies - set up the header keywords to define the
            #    third axis as frequency
            header['CTYPE3'] = 'VELOCITY'
            header['CUNIT3'] = 'km/s'
            header['CRPIX3'] = int(nf//2)+1
            if nf%2==0: # even
                header['CRVAL3'] = delta_velocity/2.
            else:
                header['CRVAL3'] = 0.0
            # Calculate the frequency step, assuming equal steps between all:
            header['CDELT3'] = delta_velocity
            header['RESTFRQ']=cc/(lam[int(nf//2)+1]*1.0e-4)

        else:
            header['CTYPE3'] = 'FREQ-LSR'
            header['CUNIT3'] = 'HZ'
            header['CRPIX3'] = int(nf//2)+1
            header['CRVAL3'] = cc*1.0e4/lam[int(nf//2)+1] # central channel if odd, 
            # Calculate the frequency step, assuming equal steps between all:
            header['CDELT3'] = delta_freq
            header['RESTFRQ']=cc/(lam[int(nf//2)+1]*1.0e-4)
    else:                # only one frequency
        header['RESTFRQ'] = cc/(lam[0]*1.0e-4)

    header['CUNIT3'] = 'Hz'

    # Make a FITS file!
    #

    # PAD IMAGE
    image_in_jypix_shifted=fpad_image(image_in_jypix_shifted, Npixf, Npixf, nx, ny)

    if len(background_args) != 0:
        for iback in background_args:
            image_in_jypix_shifted=image_in_jypix_shifted + background_object(*iback)

    image_in_jypix_float=image_in_jypix_shifted.astype(np.float32)
    fits.writeto(path_fits, image_in_jypix_float, header, output_verify='fix')







def load_image(path_image, dpc):

    f=open(path_image,'r')
    iformat=int(f.readline())

    if (iformat < 1) or (iformat > 4):
        sys.exit("ERROR: File format of image not recognized")

    nx, ny = tuple(np.array(f.readline().split(),dtype=int))
    nf = int(f.readline()) # number of wavelengths
    sizepix_x, sizepix_y = tuple(np.array(f.readline().split(),dtype=float))

    lam = np.empty(nf)
    for i in range(nf):
        lam[i] = float(f.readline())
    
    f.readline()  

    image = np.zeros((1,nf,ny,nx), dtype=float)

    for k in range(nf):
        for j in range(ny):
            for i in range(nx):

                image[0,k,j,i] = float(f.readline())

                if (j == ny-1) and (i == nx-1):
                    f.readline()

    f.close()

    # Compute the flux in this image as seen at dpc (pc)    
    pixdeg_x = 180.0*(sizepix_x/(dpc*pc))/np.pi
    pixdeg_y = 180.0*(sizepix_y/(dpc*pc))/np.pi

    # Compute the conversion factor from erg/cm^2/s/Hz/ster to erg/cm^2/s/Hz/ster at dpc
    pixsurf_ster = pixdeg_x*pixdeg_y * (np.pi/180.)**2
    factor = 1e+23 * pixsurf_ster
    # And scale the image array accordingly:
    image_in_jypix = factor * image

    return image_in_jypix, nx, ny, nf, lam, pixdeg_x, pixdeg_y


def star_pix(nx, omega):

    omega= omega%360.0

    if nx%2==0.0: # even number of pixels
        if omega>=0.0 and omega<=90.0:
            istar=nx//2 
            jstar=nx//2 
        elif omega>90.0 and omega<=180.0:
            istar=nx//2-1 
            jstar=nx//2 
        elif omega>180.0 and omega<=270.0:
            istar=nx//2-1 
            jstar=nx//2-1
        elif omega>270.0 and omega<360.0:
            istar=nx//2 
            jstar=nx//2 -1
    else:
        istar=nx//2
        jstar=nx//2
    return istar, jstar

def shift_image(image, mx, my, pixdeg_x, pixdeg_y, omega=0.0 ):

    if mx ==0.0 and my==0.0: return image

    mvx_pix=(mx/(pixdeg_x*3600.0))
    mvy_pix=(my/(pixdeg_y*3600.0))

    shiftVector=(0.0, 0.0, mvy_pix, -mvx_pix) # minus sign as left is positive 
    # cp star and remove it
    istar, jstar=star_pix(len(image[0,0,0,:]), omega)
    Fstar=image[0,0,jstar,istar]
    # print  image[0,0,jstar-1,istar-1  ]

    image[0,0,jstar,istar]=0.0
    
    # shift
    image_shifted=shift(image,shift=shiftVector, order=3)#,mode='wrap')
    # add star in new position
    image_shifted[0,0,jstar+int(mvy_pix),istar-int(mvx_pix)]=Fstar

    return image_shifted

def fpad_image(image_in, pad_x, pad_y, nx, ny):

    if image_in.shape[-2:] != (pad_x,pad_y):
        pad_image = np.zeros((1,1,pad_x,pad_y))
        if nx%2==0 and ny%2==0: # even number of pixels
            pad_image[0,0,
                      pad_y//2-ny//2:pad_y//2+ny//2,
                      pad_x//2-nx//2:pad_x//2+nx//2] = image_in[0,0,:,:]
        else:                  # odd number of pixels
            pad_image[0,0,
                      pad_y//2-(ny-1)//2:pad_y//2+(ny+1)//2,
                      pad_x//2-(nx-1)//2:pad_x//2+(nx+1)//2] = image_in[0,0,:,:]
        return pad_image

    else:                      # padding is not necessary as image is already the right size (potential bug if nx>pad_x)
        return image_in



def Gauss2d(xi , yi, x0,y0,sigx,sigy,theta):

        xp= (xi-x0)*np.cos(theta) + (yi-y0)*np.sin(theta)
        yp= -(xi-x0)*np.sin(theta) + (yi-y0)*np.cos(theta)

        a=1.0/(2.0*sigx**2.0)
        b=1.0/(2.0*sigy**2.0)

        return np.exp(- ( a*(xp)**2.0 + b*(yp)**2.0 ) )#/(2.0*np.pi*sigx*sigy)

def background_object(Ni, dpix, Flux, offx, offy, Rmaj, Rmin, Rpa):

    

    Xmax=(Ni-1)*dpix/2.
        
    xs=np.linspace(Xmax, -Xmax, Ni)
    ys=np.linspace(-Xmax, Xmax, Ni)

    Xs, Ys =np.meshgrid(xs, ys)
    rs=np.sqrt( (Xs-offx)**2. + (Ys-offy)**2. )

    F=Gauss2d(Xs , Ys, offx, offy, Rmaj, Rmin, -(Rpa+90.)*np.pi/180.)
    return F*Flux/np.sum(F)
