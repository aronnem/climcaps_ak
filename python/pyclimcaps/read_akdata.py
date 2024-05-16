from netCDF4 import Dataset
import numpy as np

from . import akutils

def read_climcaps_akdata(ncfile, ifoot, iscan, mol_name):
    """
    helper to read required data from a CLIMCAPS product file.
    This function only reads the needed fields out of the product file,
    and does no calculations (other than converting pressure units.)

    Parameters:
    -----------
    ncfile : str
        file + path to CLIMCAPS netCDF product file
    ifoot : int
        footprint (crosstrack) index to extract
    iscan : int
        along track scan index to extract
    mol_name : str
        string name specifying the desired averaging kernel
        options: air_temp, h2o_vap, ch4, co, co2, o3, hno3

    Returns
    -------
    air_pres : ndarray
        pressure profile, in hPa, of pressure levels.
        Field /air_pres in product file.
    surf_pres : float
        surface pressure, in hPa.
        aux/prior_surf_pres in product file
    htop : int
        ave_kern/mol_name_func_htop
    hbot : int
        ave_kern/mol_name_func_hbot
    ak_pidx : ndarray of int
        the index values, into the fine pressure level grid (air_pres),
        of the coarse trapezoid functions.
        ave_kern/mol_name_func_indxs in product
    ak : ndarray
        2D array, the coarse averaging kernel
        ave_kern/mol_name_ave_kern in product
    ak_nfunc : int
        the number of valid coarse AK layers; can be smaller than the
        typical number over high topography (low surface pressure)
        ave_kern/mol_name_func_last_indx in product
    ak_peff : ndarray
        unadjusted mean layer pressure for the coarse AK layers.
        ave_kern/mol_name_func_pres in product
    """

    with Dataset(ncfile, 'r') as nc:

        nc.set_auto_mask(False)

        air_pres = nc['air_pres'][:]
        surf_pres = nc['aux/prior_surf_pres'][iscan,ifoot]
        pres_nsurf = nc['air_pres_lay_nsurf'][iscan,ifoot]

        htop = nc['ave_kern/'+mol_name+'_func_htop'][0]
        hbot = nc['ave_kern/'+mol_name+'_func_hbot'][0]

        ak_pidx = nc['ave_kern/'+mol_name+'_func_indxs'][:]
        ak = nc['ave_kern/'+mol_name+'_ave_kern'][iscan,ifoot,:,:]
        ak_nfunc = nc['ave_kern/'+mol_name+'_func_last_indx'][iscan,ifoot]
        ak_peff = nc['ave_kern/'+mol_name+'_func_pres'][:]

    # Pa to hPa conversions
    air_pres /= 100.0
    surf_pres /= 100.0
    ak_peff /= 100.0
    
    return (air_pres, surf_pres, pres_nsurf,
            ak_pidx, htop, hbot, ak, ak_nfunc, ak_peff)


def calc_climcaps_ak(air_pres, surf_pres, pres_nsurf,
                     ak_pidx, htop, hbot, ak, ak_nfunc, ak_peff):
    """
    call utility functions to compute AK and other matrices
    from data stored in a CLIMCAPS data product file.

    Parameters:
    -----------
    identical to the return values from read_climcaps_akdata()

    Returns:
    --------
    Fmatrix : ndarray
        Matrix containing the trapezoid functions, shape is
        (L, j) elements, with L fine grid pressure levels and
        j coarse layers. L will be less than 100 (the full RT
        pressure level grid) because the truncation at the surface.
        j is equal to the number of trapezoids (coarse layers),
        or can be smaller if the surface pressure is low enough
        to require truncation of any coarse layers.
        (referred to as "F" in Maddy&Barnet 2008)
    Finv : ndarray
        Pseudoinverse of F matrix, shape (j, L)
        (referred to as "F+" in Maddy&Barnet 2008)
    AKcoarse : ndarray
        The coarse layer averaging kernel, shaped (j, j)
    Pcoarse : ndarray
        Pressure for AKcoarse, shaped (j)
    AKfine : ndarray
        The fine layer averaging kernel, shaped (L, L)
    Pfine : ndarray
        Pressure for AKfine, shaped (L)
    Skernel : ndarray
        Smoothing kernels, an (L, L) shaped matrix
        (FF+ in Maddy&Barnet 2008)    
    """

    # steps:
    # Adjust bottom layer
    AKcoarse, ak_pidx, Pcoarse = akutils.adjust_ak_bottom_layer(
        air_pres, ak, ak_pidx, ak_nfunc, pres_nsurf, ak_peff)
    
    # compute matrices
    ret_nlev = ak_pidx[-1]
    Fmatrix, Finv = akutils.calc_finv_mp(ak_nfunc, ak_pidx, ret_nlev,
                                         htop, hbot, air_pres)

    # create AK, etc.
    nL, nj = Fmatrix.shape
    AKfine = Fmatrix @ AKcoarse @ Finv
    Pfine = air_pres[0:nL]
    Skernel = Fmatrix @ Finv

    return Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel


def return_climcaps_akdata(ncfile, ifoot, iscan, mol_name):
    """
    function to read the required data from a CLIMCAPS netCDF file,
    and compute the F matrix (trapezoid functions) or the averaging
    kernels on coarse or fine pressure levels. The lower altitude
    (higher pressure) trapezoids are adjusted according to the surface
    pressure. This may involve removing one or more trapezoids if the
    surface pressure is low enough.

    Parameters:
    -----------
    ncfile : str
        file + path to CLIMCAPS netCDF product file
    ifoot : int
        footprint (crosstrack) index to extract
    iscan : int
        along track scan index to extract
    mol_name : str
        string name specifying the desired averaging kernel
        options: air_temp, h2o_vap, ch4, co, co2, o3, hno3

    Returns:
    --------    
    Fmatrix : ndarray
        Matrix containing the trapezoid functions, shape is
        (L, j) elements, with L fine grid pressure levels and
        j coarse layers. L will be less than 100 (the full RT
        pressure level grid) because the truncation at the surface.
        j is equal to the number of trapezoids (coarse layers),
        or can be smaller if the surface pressure is low enough
        to require truncation of any coarse layers.
        (referred to as "F" in Maddy&Barnet 2008)
    Finv : ndarray
        Pseudoinverse of F matrix, shape (j, L)
        (referred to as "F+" in Maddy&Barnet 2008)
    AKcoarse : ndarray
        The coarse layer averaging kernel, shaped (j, j)
    Pcoarse : ndarray
        Pressure for AKcoarse, shaped (j)
    AKfine : ndarray
        The fine layer averaging kernel, shaped (L, L)
    Pfine : ndarray
        Pressure for AKfine, shaped (L)
    Skernel : ndarray
        Smoothing kernels, an (L, L) shaped matrix
        (FF+ in Maddy&Barnet 2008)    

    """
    air_pres, surf_pres, pres_nsurf, ak_pidx, htop, hbot, ak, ak_nfunc, ak_peff = \
        read_climcaps_akdata(ncfile, ifoot, iscan, mol_name)

    Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel = \
        calc_climcaps_ak(air_pres, surf_pres, pres_nsurf,
                         ak_pidx, htop, hbot, ak, ak_nfunc, ak_peff)

    return Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel
