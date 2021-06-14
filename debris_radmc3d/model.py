import numpy as np
import cmath as cma
import os,sys
from debris_radmc3d.constants import *
from astropy.io import fits
from astropy.io.votable import parse


class simulation:
    """
    A class to create a complete radmc3d model 

    """

    # default_amin=1.0 # um
    # default_amax=1.0e4 # um
    # default_slope=-3.5
    # default_density=

    ### methods
    
    # def mctherm(self):
    # def make_image(self):
    # def make_cube(self):
    # def make_sed(self):
    # def convert_to_fits(self, image):
 
class dust_size_distribution:
    """
    A class used to define the dust distribution and opacities
    """
    def __init__(self, wavelength_grid, lnk_file=None,amin=1.0, amax=1.0e4, slope=-3.5, density=3.0, N_species=1, N_per_bin=50, densities=None, mass_weights=None, tag='i', compute_opct=True ):

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
            
        if isinstance(lnk_file, str): # if one optical constant given
            self.lnk_file_p=lnk_file
            self.density=density
            
        elif isinstance(self.lnk_file, list):
            self.densities =np.array(densities) if densities is not None else sys.exit('error in densities array')
            self.mass_weights=np.array(mass_weights) if mass_weights is not None else sys.exit('error in mass weights array')

            if len(lnk_file)==len(mass_weights) and len(lnk_file)==len(densities):
                print('Consider average optical constants')
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
                
            # calculate single opacity curve
            ##### THIS IS TOO SLOW, NEED TO GO BACK TO FORTRAN CODE
            # compute_opac_mie(self.lnk_file_p,
            #                      self.density,
            #                      Agrid_i/1.0e4, # in cm
            #                      self.wavelength_grid.lams/1.0e4, # in cm
            #                      wgt=weights,
            #                      # theta=None,logawidth=None,wfact=3.0,na=20,
            #                      # chopforward=0.0,errtol=0.01,verbose=False,
            #                      # extrapolate=False
            #                      verbose=True)

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



# class dust_density:
    
#     """
#     A class used to define the dust densities

#     """
# class gas_density: # needs a quick way to redefine and save
#     """
#     A class used to define the gas densities
#     """

    
# class gas_velocity:
#     """
#     A class used to define the gas velocities
#     """

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
        
       
        
    def save_spectrum(self):

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
            self.Nr==int(Nr) if Nr>0 else default_Nr
        else: self.Nr=default_Nr
        
        if Nth is not None:
            self.Nth==int(Nth) if Nth>0 else default_Nth
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
