pro read_raw_climcaps_akdata, ncfile, ifoot, iscan, mol_name, $
                         ret_pres, surf_pres, ak_pidx, htop, hbot, ak, $
                         adjust_surface=adjust_surface
  ;-------------------------------------------------
  ; read_raw_climcaps_akdata:
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
  ; ak_pidx     integer index array into pressure levels for the trapezoid
  ;             functions for the requested profile. For n coarse layers,
  ;             there will be n+1 elements in the pressure level index.
  ; htop, hbot  scalar integers specifying whether the upper and lower
  ;             functions are trapezoids or wedges
  ; ak          [n,n] array containing the averaging kernel matrix
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

  ; -------------------------
  ; STEP 2: adjust surface
  ; -------------------------

  ; Call changlev to determine number of pressure levels above surface
  ; pressure at retrieval scene
  ; ;  Name          Description
  ; ------    ----------------------
  ; surf_pres surface pressure
  ; ret_pres  Standard RTA pressure grid
  ; ret_nlev  number of levels (RTA grid)
  ; ak_nlev   number of coarse layers to be adjusted for topography
  ; ak_pres   coarse layers to be adjusted for topography 
  ;
  ; OUTPUT
  ; ak_nlev_scene    number of coarse levels at scene psurf
  ; ak_pidx_scene    coarse level boundaries at scene psurf

  if keyword_set(adjust_surface) then begin
     ret_nlev = n_elements(ret_pres)
     ak_nlev = n_elements(ak_pidx)
     changlev, surf_pres, ret_pres, ret_nlev, ak_nlev, ak_pidx, $
               ak_nlev_scene, ak_pidx_scene
     ak_pidx = ak_pidx_scene
     ; note: AK will be truncated if changlev removed a
     ; near-surface coarse layer. This happens if the surface is
     ; above one entire coarse layer.
     k = ak_nlev_scene-2
     ak = ak_full[0:k,0:k,ifoot,iscan]
  endif else begin
     ak = ak_full[*,*,ifoot,iscan]
  endelse

end



pro read_climcaps_akdata, ncfile, ifoot, iscan, mol_name, $
                          Fmatrix, Finv, AKcoarse, Pcoarse, $
                          AKfine, Pfine, Skernel

  ;-------------------------------------------------
  ; read_climcaps_akdata:
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
  read_raw_climcaps_akdata, $
     ncfile, ifoot, iscan, mol_name, $
     ret_pres, surf_pres, ak_pidx, htop, hbot, ak, /adjust_surface

  ; call helper to compute Func matrix and inverse (F and F+)
  num_func = n_elements(ak_pidx) - 1
  ret_nlev = n_elements(ret_pres)
  calc_finv_mp, num_func, ak_pidx, ret_nlev, htop, hbot, ret_pres, $
                Fmatrix, Finv

  s = size(Fmatrix, /dimensions)
  nL = s[0]
  nj = s[1]

  AKcoarse = ak
  Pcoarse = ret_pres[ak_pidx]

  AKfine = Fmatrix # ak # Finv
  Pfine = ret_pres[0:nL-1]

  Skernel = Fmatrix # Finv
  
end
