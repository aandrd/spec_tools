# the plan here is to define a function that measures the emission lines by using a gaussian function


from specutils import Spectrum1D
from specutils.fitting import fit_lines, estimate_line_parameters
from specutils.fitting.continuum import fit_continuum
from specutils.analysis import line_flux, snr, snr_derived

import astropy.units as u
from astropy.modeling.polynomial import Chebyshev1D
from astropy.io import ascii
from astropy.table import Table, hstack, vstack


# the function inputs are: 
#   - lambda of the emisison line
#   - wavelenght array of the spectrum
#   - flux array of the spectrum
#   - initial guess of the width of the emission line in angstrom


def emission_line(center, wave, spec, std):

    # set the windows to measure the continuum flux to reduced, 
    # e.g. using a 10 angstrom window, 5 angstrom away the center of the emission line
    # the continuum window should be an array of tuples
    
    window_cont = [((center-15)*u.AA, (center-5)*u.AA),
              ((center+5)*u.AA, (center+15)*u.AA)]
  
    # set the spectrum as a Spectutils object in common optical spectrum units
    spectrum = Spectrum1D(flux=spec*(u.erg * u.cm**-2 *
                          u.s**-1), spectral_axis=wave*u.AA)

    # fit the spectral continuum by using 1d-chebyshev polynomials
    continuum = fit_continuum(
        spectrum, model=Chebyshev1D(1), window=window_cont)

    # convert the continuum fitting as an array 
    y_fit_cont = np.array(continuum(wave*u.AA))
    
    # apply the continuum subtraction 
    subtraction = spec-y_fit_cont

    # set the subtracted spectrum as an Spectutils object
    spectrum_subtracted = Spectrum1D(
        flux=spec*(u.erg * u.cm**-2 * u.s**-1), spectral_axis=wave*u.AA)

    # set the array index of the emission line in a 20A window
    line_idx = np.where((wave > (center-10)) &
                        (wave < (center+10)))[0]

    # take the flux value at the peak of the spectral line
    amp = np.max(spec)

    # model a 1d-gaussian function to the emission line in the subtracted spectrum
    # with the initial guesses: amplitude, center and sigma of the gaussian
    g_init = models.Gaussian1D(
        amplitude=amp*(u.erg * u.cm**-2 * u.s**-1), mean=center*u.AA, stddev=std*u.AA)

    # fit the modelled gaussian in the subtracted spectra
    g_fit = fit_lines(spectrum_subtracted, g_init)


    # set the spectral line as a Spectutils object
    linespec = Spectrum1D(
        flux=spec[line_idx]*(u.erg * u.cm**-2 * u.s**-1), spectral_axis=wave[line_idx]*u.AA)
    
    # measure the flux of the emission line 
    lineflux = line_flux(linespec)

    # estimate the error the emission line measurement
    try:
        error = np.sqrt(g_fit.cov_matrix[0, 0])
    # if the is no fittting, the measured error is 0 
    except TypeError:
        error = 0
    # set the emission line fitting as an 20A array 
    y_fit = g_fit(wave[line_idx]*u.AA)

    # return:
    #   - wavelenght array of fitting 
    #   - profile array of the gaussian fitting 
    #   - flux profile ofthe emission line 
    #   - gaussian fitting parameters (amp, center, std)
    #   - error of the gaussian fitting 
    return wave[line_idx], y_fit, spec[line_idx], g_fit, lineflux.value, error
