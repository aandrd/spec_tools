# the plan here is to make a rebin a spectrum (x: wavelenght, y:flux ) decreasing the spectral resolution
from astropy.table import Table, hstack, vstack

def rebin_wavelength(wavelength, flux, dl): # dl is the spectral pixel size (\AA px^-1)

    # create a table with the flux arrays
    table = Table([wavelength, flux], names=['wavelength', 'flux'])

    # define the wavelenght array with the new spectral resolution (from lambda_0 to lambda_final at a step dl)
    wavelength_short = np.arange(
        wavelength[0], wavelength[len(wavelength) - 1], dl)

    new_flux = [] # define a void array 
    for i in range(len(wavelength_short)): # cicle in the rebinned wavelenght array
        if i < (len(wavelength_short) - 1): # skip the final value of the wavelenght array
            # select the those flux values between lambda_i and lambda_i+1
            
            left = wavelength_short[i] 
            right = wavelength_short[i + 1]

            flux = table['flux'][(wavelength > left) & (wavelength < right)]

            # take the mean value, this can be changes to sum, median or mean
            flux_binned = np.mean(flux) 
            new_flux.append(flux_binned) #append to the list

        else:
            last_flux = table['flux'][i]
            new_flux.append(last_flux)

    # return the spectra with the new binning 
    return wavelength_short+(dl/2), np.array(new_flux)
