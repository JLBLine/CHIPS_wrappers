import numpy as np
from copy import deepcopy
import os
import sys
import warnings
from astropy.cosmology import LambdaCDM
from numpy.ma import masked_array
from astropy import constants as const
from chips_wrappers.plotting.plot_2D import plot_wedge_cut

##Speed of light in m/s
SPEED_LIGHT = const.c.value

##21cm wavelength in m
WAVELENGTH_21CM = 0.2110611405

##Boltzmann constant
BOLTZMANN = const.k_B.value

##Lower k perp bins that always get ignored from outputs
KPERP_START = 2

##Upper k_parallel bins that always get ignored from outputs
KPARRA_END = 1

class ChipsDataProducts(object):
    """Class to read in a bunch of cosmological constants and obserational
    parameters, read in CHIPS data products, and apply them to the CHIPS
    outputs to format them into a 1D or 2D array for plotting"""

    def __init__(self, parser_args):
        """Setup CHIPS cosmological coords and obserational setup based on
        the arguments provided by the parser from argparse"""

        self.parser_args = parser_args

    def _create_k_coords(self, u_arr, v_arr, eta, z):
        """
        Convert u distance and eta into k_perp, k_parallel as per Morales et al. (2004).
        Constructs a LambdaCDM cosmology using the values from self.parser_args.

        Stores some cosmological constants for later calculations, being
        self.DM, self.Ez, self.hubble_distance.

        Parameters
        ----------
        u_arr : numpy array, float
            NDarray of u values. Should be in wavelengths.
        eta : numpy array, float
            1Darray of eta values.
        z : float
            Redshift at the central frequency of the band.

        Returns
        -------
        k_parra : numpy array, float
            NDarray of kparra values. Should be in units of h*Mpc^-1.
        k_perp : numpy array, float
            NDarray of kperp values. Should be in units of h*Mpc^-1.
        k_z : numpy array, float
            1Darray of kz values. Should be in units of h*Mpc^-1.
        """
        
        ##TODO allow the user to select a default astropy cosmology like Planck18
        # from astropy.cosmology import Planck18
        # cosmology = Planck18

        ##Create the cosmology and report the proprties
        cosmology = LambdaCDM(H0=self.parser_args.hubble,
                              Om0=self.parser_args.omega_matter,
                              Ode0=self.parser_args.omega_lambda,
                              Ob0=self.parser_args.omega_baryon)

        if self.parser_args.verbose:

            print('Cosmology being used has the following parameters:')
            print(f"\tH0 = {cosmology.H0:.2f}")
            print(f"\tOmega Matter = {cosmology.Om0:.4f}")
            print(f"\tOmega Lambda = {cosmology.Ode0:.4f}")
            print(f"\tOmega Baryon = {cosmology.Ob0:.4f}")

        nu_21 = (SPEED_LIGHT)/(WAVELENGTH_21CM) #[Hz]

        # Cosmological scaling parameter:
        h = cosmology.H(0).value/100 # Hubble parameter.

        E_z = cosmology.efunc(z) ## Scaling function, see (Hogg 2000)

        # Cosmological distances:
        Dm = cosmology.comoving_distance(z).value*h #[Mpc/h] Transverse co-moving distance.
        DH = 3000 # [Mpc/h] Hubble distance.

        # k parallel
        k_parra = eta * (2*np.pi*nu_21*E_z)/(DH*(1 + z)**2) # [Mpc^-1 h]

        # k perpendicular
        k_perp = u_arr * (2*np.pi/Dm) # [Mpc^-1 h]

        self.DM = Dm
        self.Ez = E_z
        self.h = h
        self.hubble_distance = DH

        return k_parra, k_perp

    def setup_chips_params(self):
        """When converting data from output CHIPS files into a 1D or 2D power
        spectrum, do all the the things that are common to the two operations,
        like setting up coordinates and DMs and things.

        There are many cosmological things here that I do not understand"""
        
        parser_args = self.parser_args

        self.central_freq = parser_args.lowerfreq + int(parser_args.N_chan / 2)*parser_args.chan_width

        ##Frequency bandwidth of the data
        bandwidth = float(parser_args.N_chan)*parser_args.chan_width

        ##The eta coords (FT pair of frequency)
        self.eta = np.zeros(int(parser_args.Neta))
        for i in range (0, int(parser_args.Neta)):
            self.eta[i] = (float(i)-0.5) / bandwidth
        self.eta[0] = self.eta[1]/2.

        #The bin length on the u,v plane the data were gridded to (in wavelengths)
        self.u_arr = np.zeros(parser_args.N_kperp)
        for i in range(0,parser_args.N_kperp):
            self.u_arr[i] = float(i)*parser_args.umax*1.1/float(parser_args.N_kperp)
        self.u_arr[0] = self.u_arr[1]/2.

        ##21cm radiation frequency in m/s
        f21 = SPEED_LIGHT / WAVELENGTH_21CM

        ##Redshift based on central frequency
        z = f21 / self.central_freq - 1.

        ##Calculate the k_parallel and k_perpendicular coords
        ##also sets:
        ##    self.DM ( Transverse co-moving distance, [Mpc/h])
        ##    self.Ez (scaling function from Hogg et al 2000)

        ## based on astropy cosmology
        k_parra, k_perp = self._create_k_coords(self.u_arr, self.u_arr, self.eta, z)

        ##There are parts of the 2D outputs that always get thrown away,
        ##so make the same cut on the output coords
        self.kpa = k_parra
        self.kper = k_perp

        ##There are parts of the 2D outputs that always get thrown away,
        ##so make the same cut on the output coords
        self.kper = self.kper[KPERP_START:]


        ##New way of doing it===================================================
        cent_wavelength = SPEED_LIGHT / self.central_freq
        beam_area_steradian = 0.07597
        ##Frequency bandwidth of the data
        bandwidth = float(parser_args.N_chan)*parser_args.chan_width

        norm_numer = cent_wavelength**4 * self.DM**2*self.hubble_distance * bandwidth * (1 + z)**2 * 1.e6
        norm_denom = parser_args.N_chan*(2*BOLTZMANN*1e26)**2*beam_area_steradian*f21*self.Ez

        self.normalisation = norm_numer / norm_denom

        ##Gridding causes decohence due to the grid points not being
        ##infite. Must multiply by this factor
        ##This is a density factor correction based on Barry et al 2019.
        
        if parser_args.density_correction:
            if parser_args.density_correction == 'use_fit':
                self.decoherence_factor = "use_fit"
            else:
                self.decoherence_factor = float(parser_args.density_correction)  
        else:
        
            self.decoherence_factor = 1.0
            
        if parser_args.verbose:

            print("Params either set or calculated:")
            print(f"\tCentral wavelength (m) {cent_wavelength:.3f}")
            print(f"\tDM {self.DM:.1f}", )
            print(f"\tHubble_distance {self.hubble_distance}")
            print(f"\tBandwidth {bandwidth:.5e}")
            print(f"\tredshift {z:.3f}")
            print(f"\tNum freq chan {parser_args.N_chan}")
            print(f"\tEz", self.Ez)
            print(f"\tGridding normalisation", self.decoherence_factor)

            print(f"OVERALL NORMALISATION APPLIED (without decoherence):  {self.normalisation:.8e}")

        if parser_args.wedge_factor >= 0:
            self.wedge_factor = parser_args.wedge_factor
        else:
            ##Used to do the wedge cut for 1D
            self.wedge_factor = self.DM * self.Ez / (self.hubble_distance * (z + 1))

    def _open_and_reshape(self, filename):
        """Read in a CHIPS binary output and reshape into a 2D array"""

        with open(filename, 'rb') as train_xf:
            # print(filename)
            ##Reads in as 1D
            data = np.fromfile(train_xf, dtype=np.float32)
            ##Make 2D - data were written out by spatial direction, then spectral,
            ##reshape into 2D using the parser arguments
            
            N_chans_present = int(len(data) / self.parser_args.N_kperp)
            
            ##Ok so this can happen if performing a masked FFT in newer versions of CHIPS
            ##You get get less eta channels out, but of the same resolution. The overal
            ##frequency bandwith should be the same though, so keep N_chan alive to
            ##propagate through the rest of the code
            if N_chans_present != self.parser_args.N_chan:
                print(f'Number of eta channels {N_chans_present} in file is not full N_chan {self.parser_args.N_chan} (assuming that number k_perp chans is {self.parser_args.N_kperp} as set by --N_kperp)')
                
                ##Known drop in number of channels from 384 to 336 when using nfft
                expec_ratio = 384 / 336
                
                print(f'If this is an nfft format file, would expect {self.parser_args.N_chan / expec_ratio}')
                
                self.parser_args.Neta = int(N_chans_present / 2)
            
            # data = np.reshape(data, (self.parser_args.N_kperp,self.parser_args.N_chan))
            data = np.reshape(data, (self.parser_args.N_kperp, N_chans_present))

            ##TODO through some useful error if the reshaping cannot be done

            ##python lists the 'y' axes first, 'x' second so swap the axes
            ##to play nicely with imshow and the like
            data = np.swapaxes((data),1,0)

            ##Limit data to only the positive spectral frequencies,
            ##including the DC term
            data = data[self.parser_args.Neta-1:,:]

            ##There are some bins that always get ignored so throw them away now
            data = data[:-KPARRA_END, KPERP_START:]
            
        return data
    
    def _8s_decoherence_factor(self, k_perp_mesh):
        """
        Calculates the decoherence factor for a given set of weights, using the 8s model
        as fit in Line et at. 2023 in prep
        
        Parameters
        -----------
        weights : np.ndarray
            An array of weights out of CHIPS (from file
            e.g. crosspower_xx_0.iter.{chips_tag}.dat )
            
        Returns
        --------
        correction : np.ndarray:
            An array of corrections for decoherence factors, one for each weight
            in the input array.
        """
        
        a_1 = 0.0
        delta = 0.01
        A = 0.377280447218
        a_2 = -0.161193559625
        x_b = 0.064354113566
        
        output = A*(k_perp_mesh/x_b)**(-a_1) * (0.5*(1+(k_perp_mesh/x_b)**(1/delta)))**((a_1-a_2)*delta)
        output[output > 1.0] = 1.0
        
        return 1 / output
    
    def _read_in_data_and_convert(self, polarisation, chips_tag=False, oneD=False):
        """Attempt to read in the data based on user provided paths and polarisation
        options. Will first check for files that have had kriging (have a zero
        in the title), if it can't find those, look for ones without kriging
        (have a one in the title).

        Convert to a 2D array for 2D plot by default, or if specified, a 1D array.
        1D array requires one extra CHIPS output file to have been downloaded"""

        ##For ratio plots, need to pick denominator or numerator
        if chips_tag:
            pass
        ##if not, should only need to use the chips_tag in self.parser_args.chips_tag
        else:
            chips_tag = self.parser_args.chips_tag

        # filename0 = f"{self.parser_args.basedir}crosspower_{polarisation}_0.iter.{chips_tag}.dat"
        
        file_found = False
        ##Try various running option numbers, and stop if we find valid files
        for run_opt in [0, 1, 20, 21, 22]:
            kriging = run_opt
            filename = f"{self.parser_args.basedir}/crosspower_{polarisation}_{kriging}.iter.{chips_tag}.dat"
            
            # filename = f"{self.parser_args.basedir}totpower_{polarisation}_{kriging}.iter.{chips_tag}.dat"
            
            if os.path.isfile(filename):
                file_found = True
                break
        
        if not file_found:
            msg = 'Could not open crosspower files based in input params.\n' \
            'Searched for files like :\n' \
                f'{self.parser_args.basedir}/crosspower_{polarisation}_*.iter.{chips_tag}.dat'

        crosspower = self._open_and_reshape(filename)
        ##Number of K_parra depends on shape of data, so setup the params after reading in
        ##the crosspower
        self.setup_chips_params()

        filename = f"{self.parser_args.basedir}/outputweights_{polarisation}_{kriging}.iter.{chips_tag}.dat"

        if not os.path.isfile(filename):
            msg = 'Could not open the outputweights file with the same kriging\n' \
            'number as the crosspower. Searched for the following file:\n' \
            f'{filename}'
            sys.exit(msg)

        weights = self._open_and_reshape(filename)

        crosspower = masked_array(crosspower, weights == 0)

        # print("Just read in", crosspower)
        crosspower = crosspower / weights

        ##32 is a CHIPS based number hard coded to make the weights sensible
        ##For anything else, stick to one??
        ##TODO used to be 36, pls why?
        ##TODO sort this max weight nightmare
        # max_weights = 2
        # weight_adjustment = 32 / 2

        if self.parser_args.RTS_outputs:
            weight_scheme = 32
        else:
            weight_scheme = 1

        self.weight_data = weights/(self.normalisation)**2*weight_scheme*np.sqrt(self.parser_args.N_chan)
        
        if self.decoherence_factor == 'use_fit':
            k_perp_mesh, k_parr_mesh = np.meshgrid(self.kper, self.kpa)
            self.decoherence_factor = self._8s_decoherence_factor(k_perp_mesh)
            
        self.weight_data /= self.decoherence_factor
        
        self.crosspower = crosspower*self.decoherence_factor*self.normalisation
        
        if self.parser_args.verbose:
            print(f"Max power in file {self.crosspower.max():.4e}")
        
        self.weights = weights

        np.savez("2D_coords_and_power.npz", k_perp=self.kper,
                                            k_parr=self.kpa,
                                            twoD_power=self.crosspower,
                                            weights=weights)
        
    def read_data_and_create_2Darray(self, polarisation, chips_tag=False):
        """Attempt to read in the data based on user provided paths and polarisation
        options. Will first check for files that have had kriging (have a zero
        in the title), if it can't find those, look for ones without kriging
        (have a one in the title).

        Convert data to a 2D array for a 2D plot"""

        self._read_in_data_and_convert(polarisation, chips_tag=chips_tag, oneD=False)

        ##make an 'extent' list for imshow, that details x,y coords for a 2D plot
        ##goes as [low x coord, high x coord, low y coord, high y coord]
        extent = [self.kper[0], self.kper[-1], self.kpa[0], self.kpa[-1]]

        twoD_ps_array = self.crosspower

        return twoD_ps_array, extent

    def _grid_2D_to_1D_k3_inside(self, k_perp, k_parra, twoD_data, twoD_weights, ktot_bin_edges, convert_to_delta=True):
        twoD_weights_sqrt = np.sqrt(twoD_weights)

        ##meshgrid them to find the length in kspace of all bins
        k_perp_mesh, k_parr_mesh = np.meshgrid(k_perp, k_parra)
        k_lengths_mesh = np.sqrt(k_perp_mesh**2 + k_parr_mesh**2)

        ##This does the wedge cut?
        wedge_cut = k_parr_mesh > (0.5*np.pi)*k_perp_mesh*self.wedge_factor
        
        ##Cuts off in k perpendicular, avoids small spatial scales in the 1D
        k_perp_cut = k_perp_mesh <= self.parser_args.kperp_max

        ##Cuts off in k perpendicular, avoids large spatial scales in the 1D
        k_perp_cut_min = k_perp_mesh > self.parser_args.kperp_min

        ##Cuts off in k parallel, avoiding things close to the wedge
        kparra_cut = k_parr_mesh > self.parser_args.kparra_min

        nozero_per = k_perp_mesh > 0.0
        nozero_par = k_parr_mesh > 0.0
        
        only_pos = twoD_data > 0.0

        ##Find the centre of all the bins as the gridding coords
        ktot_bins = (ktot_bin_edges[1:] + ktot_bin_edges[:-1])/2

        ##How many bins we have, and make a zero array for gridding 1D power
        num_ktot_bins = len(ktot_bins)
        oneD_power = np.zeros(int(num_ktot_bins))
        oneD_delta = np.zeros(int(num_ktot_bins))
        oneD_weights = np.zeros(int(num_ktot_bins))
        
        oneD_power_std = np.zeros(int(num_ktot_bins))
        self.bin_count = np.zeros(int(num_ktot_bins))
        ##Keep track of the locations of the bins on the 2D PS if we want to
        ##plot the wedge cut
        binning_array = np.ones(twoD_data.shape)*-1.0

        for k_tot_ind in range(num_ktot_bins):
            ##This finds all bins that sit inside the current annulus
            above_min = k_lengths_mesh > ktot_bin_edges[k_tot_ind]
            below_max = k_lengths_mesh <= ktot_bin_edges[k_tot_ind + 1]

            cut_inds = np.where(above_min & below_max & wedge_cut & k_perp_cut & nozero_par & nozero_per & kparra_cut & k_perp_cut_min)
            # cut_inds = np.where(above_min & below_max & wedge_cut & k_perp_cut & nozero_par & nozero_per & kparra_cut & k_perp_cut_min & only_pos)

            ##Always get annoying warnnging that a Masked element has been set
            ##to NaN here, so ignore them
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning,
                    message="Warning: converting a masked element to nan.")

                # oneD_delta[k_tot_ind] = np.nansum(twoD_data[cut_inds]*twoD_weights_sqrt[cut_inds]**2*k_lengths_mesh[cut_inds]**3)
                # oneD_weights[k_tot_ind] = np.nansum(twoD_weights_sqrt[cut_inds]**2)
                # oneD_power[k_tot_ind] = np.nansum(twoD_data[cut_inds]*twoD_weights_sqrt[cut_inds]**2)
                # oneD_power_std[k_tot_ind] = np.nanstd(twoD_data[cut_inds])
                
                oneD_delta[k_tot_ind] = np.nansum(twoD_data[cut_inds]*twoD_weights[cut_inds]*k_lengths_mesh[cut_inds]**3)
                oneD_weights[k_tot_ind] = np.nansum(twoD_weights[cut_inds])
                oneD_power[k_tot_ind] = np.nansum(twoD_data[cut_inds]*twoD_weights[cut_inds])
                
                weight_mean = oneD_power[k_tot_ind] / oneD_weights[k_tot_ind]
                
                
                # normed_weights = twoD_weights[cut_inds] / np.nansum(twoD_weights[cut_inds]) 
                
                normed_weights = twoD_weights[cut_inds]
                
                if np.nansum(normed_weights) == 0.0:
                    pass
                else:
                
                
                    ##weighted sample variance
                    oneD_power_std[k_tot_ind] = np.nansum(normed_weights*(twoD_data[cut_inds] - weight_mean)**2) / np.nansum(normed_weights)
                
                self.bin_count[k_tot_ind] = int(cut_inds[0].size)
                
            binning_array[cut_inds] = k_tot_ind + 1
            
        

        self.binning_array = binning_array
        # self.twoD_weights

        oneD_power = oneD_power / oneD_weights
        oneD_delta = oneD_delta / oneD_weights
        oneD_power_std = np.sqrt(oneD_power_std)
        
        ##Which sigma to report the noise to
        sigma = 2
        sqrt_weights = np.ones(len(oneD_weights))
        sqrt_weights[np.where(oneD_weights != 0)] = np.sqrt(oneD_weights[np.where(oneD_weights != 0)])
        oneD_noise = sigma / sqrt_weights

        ##Convert to delta
        oneD_delta = oneD_delta / (2*np.pi**2)
        oneD_noise = oneD_noise*ktot_bins**3 / (2*np.pi**2)
        
        # np.save("oneD_weights.npy", oneD_weights)
        
        self.oneD_power_std = oneD_power_std

        return ktot_bins, oneD_noise, oneD_power, oneD_delta

    def read_data_and_create_1Darray(self, polarisation, chips_tag=False):

        self._read_in_data_and_convert(polarisation, chips_tag=chips_tag, oneD=True)

        if self.parser_args.ktot_bin_edges:

            ktot_bin_edges = np.loadtxt(self.parser_args.ktot_bin_edges)

            # if os.path.isfile("./"+self.parser_args.ktot_bin_edges):
            #     ktot_bin_edges = np.loadtxt(self.parser_args.ktot_bin_edges)
            # else:
            #     exit(f"Cannot find --ktot_bin_edges={self.parser_args.ktot_bin_edges}, file doesn't exist")
        else:
            low_k_edge = self.parser_args.low_k_edge
            high_k_edge = self.parser_args.high_k_edge
            num_k_edges = self.parser_args.num_k_edges
            ktot_bin_edges = 10**np.linspace(np.log10(low_k_edge), np.log10(high_k_edge), num_k_edges)

        self.ktot_bin_edges = ktot_bin_edges

        # print("BOTTOM EDGE, UPPER EDGE", ktot_bin_edges[0],  ktot_bin_edges[-1])

        ktot_bins, oneD_noise, oneD_power, oneD_delta = self._grid_2D_to_1D_k3_inside(self.kper,
             self.kpa, self.crosspower, self.weight_data, ktot_bin_edges)

        ##If requested, make a 2D plot of the cuts applied
        if self.parser_args.plot_wedge_cut_2D:
            plot_wedge_cut(self)

        return ktot_bins, oneD_noise, oneD_power, oneD_delta
    
    def get_horizon_and_beam_lines(self):
        """Calculates the horizon line and primary beam line for plotting on"""
        
        grad_horiz = 0.5*np.pi # Horizon cut gradient.
        # grad_horiz = 1 # Horizon cut gradient.
        
        ##say a tile is like 4.5 metres across, FoV is lambda/D
        beam_width = (SPEED_LIGHT / self.central_freq) / 4.5
        grad_beam = beam_width # Beam FoV cut.
        # wedge_cut = grad*polySpectra.wedge_factor(self.z,self.cosmo)
        
        line_beam = grad_beam*self.wedge_factor*self.kper
        line_horiz = grad_horiz*self.wedge_factor*self.kper
        
        self.line_beam = line_beam
        self.line_horiz = line_horiz