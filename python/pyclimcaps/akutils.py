"""
Utility functions for CLIMCAPS averaging kernels
"""

from datetime import datetime
import numpy as np

def climcaps_co2_prior(utc_time):
    """
    the climcaps prior: a global mean with secular increase,
    Page 2-59, CLIMCAPS Science User's guide

    input is assumed to be UTC time tuples, as is stored in the CLIMCAPS product.

    """

    b = 371.92429
    c = 1.8406818
    t0 = 2002.0

    # year fraction is the time length from input UTC to
    # Jan-01 00:00:00 in the same year (dt-dt1), divided by the length of
    # year we are in (dt2-dt1)
    dt1 = datetime(utc_time[0],   1, 1, 0, 0, 0)
    dt2 = datetime(utc_time[0]+1, 1, 1, 0, 0, 0)
    dt  = datetime(utc_time[0], utc_time[1], utc_time[2],
                   utc_time[3], utc_time[4], utc_time[5])
    year_fraction = (dt-dt1).total_seconds() / (dt2-dt1).total_seconds()

    fyear = utc_time[0] + year_fraction

    co2_prior = b + c*(fyear - t0)

    return co2_prior


def calc_BLmult(surf_pres, air_pres, pres_nsurf):
    """
    calculate the BL Mult factor
    """
    pdiff1 = surf_pres - air_pres[pres_nsurf-1]
    pdiff2 = air_pres[pres_nsurf] - air_pres[pres_nsurf-1]

    BLmult = pdiff1 / pdiff2

    return BLmult


def adjust_profile_bottom_layer(surf_pres, air_pres, pres_nsurf, prof, level_quantity=True):
    """
    adjust the surface containing layer with BLmult
    This either applies a linear extrapolation to the first below-surface pressure level,
    (should be used for temperature or CO2), if level_quantity is True; or scales the
    surface-containing layer by BLmult (all other retrieved molecular species that have
    layer-integrated molecular density.), if level_quantity is False.
    
    """
    BLmult = calc_BLmult(surf_pres, air_pres, pres_nsurf)
    if level_quantity:
        prof_delta = prof[pres_nsurf] - prof[pres_nsurf-1]
        prof[pres_nsurf-1] = BLmult * prof_delta
    else:
        prof[pres_nsurf-1] *= BLmult

    return prof


def adjust_ak_bottom_layer(air_pres, ak, ak_pidx, ak_nfunc, pres_nsurf, ak_peff):
    """
    adjust lower layer(s) of AK to account for surface pressure
    This will remove bottom layers from the AK (possibly reducing it's size).
    It also adjust the coarse mean layer pressures associated with the AK
    functions.
    """
    ak = ak[:ak_nfunc,:ak_nfunc].copy()
    ak_pidx = ak_pidx[:ak_nfunc+1].copy()
    ak_pidx[ak_nfunc] = pres_nsurf

    Pcoarse = ak_peff[:ak_nfunc].copy()
    bot_pidx = ak_pidx[ak_nfunc]
    top_pidx = ak_pidx[ak_nfunc-1]
    bot_pdiff = air_pres[bot_pidx-1] - air_pres[top_pidx-1]
    # double check here - should we be setting
    # ak_nfunc-1 or ak_nfunc.
    Pcoarse[ak_nfunc-1] = bot_pdiff/np.log(air_pres[bot_pidx-1]/air_pres[top_pidx-1])

    return ak, ak_pidx, Pcoarse


def calc_finv_mp(num_func, func_indx, ret_nlev, htop, hbot, air_pres):
    """
    calculates a scene-dependent transformation matrix (F_matrix) and
    its inverse the using Moore-Penrose pseudoinverse technique
    The transformation matrix consists of the retrieval variable
    trapezoid state functions [e.g., 21 for water vapor, 9 for ozone]
    on the standard 100 pressure level retrieval grid, e.g., dimensions [9 x 100]
    This matrix is scene-dependent because we use only those functions and
    pressure levels that are above Earth surface at a target scene.

    Parameters
    ----------
    num_func : int
        number of state functions above Earth surface [ave_kern/*_func_last_indx]
    func_index : ndarray of int
        trapezoid state function hinge-points as reported in [ave_kern/*_func_idxs]
    ret_nlev : int
        number of retrieval pressure levels (air_pres) above Earth surface
    htop : int
        value in [ave_kern/*_func_hbot]
    hbot : int
        value in [ave_kern/*_func_htop]
    air_pres : ndarray
        standard 100 level pressure grid [air_pres]/100 in hPa units

    Output:
    -------
    f_matrix : ndarray
        2D array shaped (L,j), where L=retrieval levels, j=trapezoid layers
    """

    # ------------------------------------------------
    # Step 1: calculate the transformation matrix: f_matrix
    # ------------------------------------------------
    	
    f_matrix = np.zeros((ret_nlev,num_func)) # [nL, nj]

    for ifunc in np.arange(num_func): 

        # Call slb2fin for one state function at a time setting the corresponding slbval = 1.0
        slbval = np.zeros(num_func)
       	slbval[ifunc] = 1.0
       	f_matrix[:,ifunc] = slb2fin(num_func, func_indx, ret_nlev,
                                    slbval, htop, hbot, air_pres)
       
    # ------------------------------------------------
    # Step 2: calculate the inverse of f_matrix using 
    #         the Moore-Penrose pseudoinverse method
    # ------------------------------------------------

    # This matrix, fftr, should have no zeros on diagonal, otherwise inversion fails
    fftr = f_matrix.T @ f_matrix # [nj x nj]
    finv1 = np.linalg.pinv(fftr)
    f_inv = finv1 @ f_matrix.T   # [nj x nL]
 
    return f_matrix, f_inv


def slb2fin(num_func, func_indx, ret_nlev, func_ampl, usehalftop, usehalfbot, air_pres,
            debug_print=False):
    """
    construct a trapezoid state function on 100 pressure levels
    one at a time, using [ave_kern/*_func_indxs] hinge points

    Parameters:
    -----------
    num_func : int
        number of trapezoid state functions above Earth surface at retrieval scene
    func_indx : ndarray of int
        index values for trapezoid hinge points, starting at 1
    ret_nlev : int
        number of fine grid levels
    func_ampl : ndarray
        amplitude of trapezoids, which is 1.0 in this case
        this is an ndarray with num_func elements, with the
        trapezoid amplitude set for one function, and zero elsewhere.
        for example, with 5 functions, to compute the second
        function, func_ampl would be [0,1,0,0,0].
    usehalftop : int
        0 is a trapezoid, 1 is a wedge
    usehalfbot : int
        0 is a trapezoid, 1 is a wedge
    air_pres : ndarray
        100-level retrieval pressure grid [air_pres]/100. in units [hPa]
    debug_print : bool
        If true, prints some information to console.

    returns:
    --------
    state_func_fine : ndarray
        the trapezoid state function, defined on the "fine" pressure level grid.
    """

    # allocate output ndarray - must be the largest index pressure level
    # referenced in the func_idx. (should match func_indx[-1]-1)
    state_func_fine = np.zeros(ret_nlev)

    # Construct the face of the trapezoid state function
    # the 'face' values are amplitudes associated with the hinge points,
    # therefore one additional value is needed.
    state_face = np.zeros(num_func + 1)

    # first hinge point value depends on whether the first function is a
    # wedge or trapezoid depending on usehalftop: If wedge, use full amplitude,
    # else use half.

    if usehalftop > 0:
        state_face[0] = 0.5 * func_ampl[0]
    else:
        state_face[0] = func_ampl[0]

    for n in np.arange(1, num_func):
        state_face[n] = 0.5 * (func_ampl[n] + func_ampl[n-1])

    if usehalfbot > 0:
        state_face[num_func] = 0.5 * func_ampl[num_func-1]
    else:
        state_face[num_func] = func_ampl[num_func-1]

    if debug_print:
        print('func_ampl: ', func_ampl)
        print('n iup idn state slope:')

    for n in range(num_func):
        idx_up = func_indx[n] - 1
        idx_down = func_indx[n+1] - 1

        state_up = state_face[n]
        state_down = state_face[n+1]

        pres_up = np.log(air_pres[idx_up])
        pres_down = np.log(air_pres[idx_down])

        slope = (state_down - state_up)/(pres_down - pres_up)
        for L in range(idx_up, idx_down):
            state_func_fine[L] = state_up + slope*(np.log(air_pres[L]) - pres_up)

        if debug_print:
            print('{:1d} {:3d} {:3d} {:4.1f} {:4.1f} {:12.6f}'.format(
                n, idx_up, idx_down, state_up, state_down, slope))

    state_func_fine[idx_down] = state_down

    return state_func_fine
