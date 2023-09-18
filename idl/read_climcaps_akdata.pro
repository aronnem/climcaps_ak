pro read_climcaps_akdata, ncfile, ifoot, iscan, mol_name, $
                          ret_pres, surf_pres, pres_nsurf, ak_pidx, $
                          htop, hbot, ak, ak_nfunc, ak_peff, $
                          adjust_surface=adjust_surface
  ;-------------------------------------------------
  ; read_climcaps_akdata:
  ;
  ; Helper procedure to read the required data from a CLIMCAPS netCDF file
  ; for computing the trapezoid functions, and apply a surface correction
  ; to adjust the near surface trapezoid function boundaries.
  ;
  ; INPUT:
  ;  Name              Description
  ; ---------    -----------------------------
  ; ncfile       string file + path to CLIMCAPS netCDF product file
  ; ifoot        footprint (crosstrack) index to extract
  ; iscan        along track scan index to extract
  ; mol_name     string name specifying the desired averaging kernel
  ;              options: air_temp, h2o_vap, ch4, co, co2, o3, hno3
  ;
  ; adjust_surface optional keyword specifying whether to apply
  ;                the surface adjustment. This is mainly for testing,
  ;                in most cases this should be enabled to adjust the
  ;                surface.
  ;
  ; OUTPUT:
  ;  Name              Description
  ; ---------   -----------------------------
  ; ret_pres    pressure levels in RT grid (100 elements), units hPa
  ; surf_pres   scalar surface pressure at requested footprint
  ; pres_nsurf  Scalar index specifying which element of ret_pres
  ;             is closest to surf_pres
  ; ak_pidx     integer index array of pressure level boundaries for the
  ;             trapezoid functions. For n coarse layers,
  ;             there will be n+1 elements in the pressure level index.
  ; htop, hbot  scalar integers specifying whether the upper and lower
  ;             functions are trapezoids or wedges
  ; ak          [n,n] array containing the averaging kernel matrix
  ; ak_nfunc    Scalar integer specifying the number of aks above surf_pres
  ; ak_peff     Effective pressure values of coarse ak layers, read from 
  ;             L2 file and calculated as the log-pressure value between
  ;             two pressure levels
  ;---------------------------------------------------------------


  ;--------
  ; STEP 1
  ;--------
  ; read the raw climcaps averaging kernel data, along with the other
  ; indexing values needed to produce the AK on the 'fine' levels

  fid = ncdf_open(ncfile)

  varid = ncdf_varid(fid, 'air_pres')
  ncdf_varget, fid, varid, ret_pres
  ; Pa to hPa
  ret_pres = ret_pres / 100.0

  varid = ncdf_varid(fid, 'air_pres_lay_nsurf')
  ncdf_varget, fid, varid, pres_nsurf
  pres_nsurf = pres_nsurf[ifoot,iscan]

  groupids = ncdf_groupsinq(fid)
  for i=0, n_elements(groupids)-1 do begin
     gname = ncdf_groupname(groupids[i])
     if strcmp(gname, 'aux') eq 1 then aux_grpid = groupids[i]
     if strcmp(gname, 'ave_kern') eq 1 then ak_grpid = groupids[i]
  end

  varid = ncdf_varid(aux_grpid, 'prior_surf_pres')
  ncdf_varget, aux_grpid, varid, surf_pres_all
  ; Pa to hPa conversion
  surf_pres = surf_pres_all[ifoot, iscan]/100

  varid = ncdf_varid(ak_grpid, mol_name + '_func_htop')
  ncdf_varget, ak_grpid, varid, htop
  varid = ncdf_varid(ak_grpid, mol_name + '_func_hbot')
  ncdf_varget, ak_grpid, varid, hbot

  varid = ncdf_varid(ak_grpid, mol_name + '_func_indxs')
  ncdf_varget, ak_grpid, varid, ak_pidx

  varid = ncdf_varid(ak_grpid, mol_name + '_ave_kern')
  ncdf_varget, ak_grpid, varid, ak_full

  ak = ak_full[*,*,ifoot,iscan]

  varid = ncdf_varid(ak_grpid, mol_name + '_func_last_indx')
  ncdf_varget, ak_grpid, varid, ak_nfunc
  ak_nfunc = ak_nfunc[ifoot,iscan]

  varid = ncdf_varid(ak_grpid, mol_name + '_func_pres')
  ncdf_varget, ak_grpid, varid, ak_peff
  ; Pa to hPa
  ak_peff = ak_peff/100.0
end



pro return_climcaps_akdata, ncfile, ifoot, iscan, mol_name, $
                            Fmatrix, Finv, AKcoarse, Pcoarse, $
                            AKfine, Pfine, Skernel

  ;-------------------------------------------------
  ; return_climcaps_akdata:
  ;
  ; procedure to read the required data from a CLIMCAPS netCDF file,
  ; and compute the F matrix (trapezoid functions) or the averaging
  ; kernels on coarse or fine pressure levels.
  ;
  ; INPUT:
  ;  Name              Description
  ; ---------    -----------------------------
  ; ncfile       string file + path to CLIMCAPS netCDF product file
  ; ifoot        footprint (crosstrack) index to extract
  ; iscan        along track scan index to extract
  ; mol_name     string name specifying the desired averaging kernel
  ;              options: air_temp, h2o_vap, ch4, co, co2, o3, hno3
  ;
  ; OUTPUT:
  ;  Name              Description
  ; ---------   -----------------------------
  ; Fmatrix     Matrix containing the trapezoid functions, shape is
  ;             (L, j) elements, with L fine grid pressure levels and
  ;             j coarse layers. L will be less than 100 (the full RT
  ;             pressure level grid) because the truncation at the surface.
  ;             j is equal to the number of trapezoids (coarse layers),
  ;             or can be smaller if the surface pressure is low enough
  ;             to require truncation of any coarse layers.
  ;             (referred to as "F" in Maddy&Barnet 2008)
  ; Finv        Pseudoinverse of F matrix, shape (j, L)
  ;             (referred to as "F+" in Maddy&Barnet 2008)
  ; AKcoarse    The coarse layer averaging kernel, shaped (j, j)
  ; Pcoarse     Pressure for AKcoarse, shaped (j)
  ; AKfine      The fine layer averaging kernel, shaped (L, L)
  ; Pfine       Pressure for AKfine, shaped (L)
  ; Skernel     Smoothing kernels, an (L, L) shaped matrix
  ;             (FF+ in Maddy&Barnet 2008)
  ;---------------------------------------------------------------


  ; call helper to get the various CLIMCAPS data, and apply
  ; surface correction
  read_climcaps_akdata, $
     ncfile, ifoot, iscan, mol_name, $
     ret_pres, surf_pres, pres_nsurf, ak_pidx, htop, hbot, ak, $
     ak_nfunc, ak_peff

  ; call helper to do all the AK calculations.
  calc_climcaps_akdata, ret_pres, surf_pres, pres_nsurf, ak_pidx, $
                        htop, hbot, ak, ak_nfunc, ak_peff, $
                        Fmatrix, Finv, AKcoarse, Pcoarse, $
                        AKfine, Pfine, Skernel, /adjust_surface
end
