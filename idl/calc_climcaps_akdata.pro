pro calc_climcaps_akdata, ret_pres, surf_pres, ak_pidx, htop, hbot, ak, $
                          Fmatrix, Finv, AKcoarse, Pcoarse, $
                          AKfine, Pfine, Skernel, $
                          adjust_surface=adjust_surface

  ;---------------------------------------------------------------
  ; calc_climcaps_akdata:
  ;
  ; compute the full averaging kernel, trapezoids, etc. for a CLIMCAPS
  ; coarse compressed averaging kernel.
  ; generally when working with the original CLIMCAPS product file,
  ; one should use the read* functions in read_climcaps_akdata.
  ; This function would be used if the CLIMCAPS data was extracted
  ; already and you need to compute the full AK with the already
  ; known values.
  ;
  ; INPUT:
  ;  Name              Description
  ; ---------    -----------------------------
  ; ret_pres     the 100 element vertical pressure grid. [hPa]
  ;              The variable name in the CLIMCAPS L2 product is '/air_pres'
  ;              but should be converted from Pa to hPa.
  ; surf_pres    Scalar, surface pressure [hPa]
  ;              '/aux/prior_surf_pres' in L2 product, converted to hPa
  ; ak_pidx      The indices into the pressure grid for the trapezoid functions
  ;              'ave_kern/<mol_name>_func_indxs' in L2 product
  ; htop         Scalar integer specifying shape of TOA trapezoid function
  ;              'ave_kern/<mol_name>_func_htop' in L2 product
  ; hbot         Scalar integer specifying shape of near surface trapezoid function
  ;              'ave_kern/<mol_name>_func_hbot' in L2 product
  ; ak           Averaging kernel matrix for coarse trapezoid functions
  ;              'ave_kern/<mol_name>_ave_kern' in L2 product
  ;
  ; adjust_surface optional keyword specifying whether to apply
  ;                the surface adjustment. This is mainly for testing,
  ;                in most cases this should be enabled to adjust the
  ;                surface.
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

  ; -------------------------
  ; STEP 1: adjust surface
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
     AKcoarse = ak[0:k,0:k]
  endif else begin
     AKcoarse = ak
  endelse

  ; -------------------------
  ; STEP 2: call helper to compute Func matrix and inverse (F and F+)
  ; -------------------------
  num_func = n_elements(ak_pidx) - 1
  ret_nlev = n_elements(ret_pres)
  calc_finv_mp, num_func, ak_pidx, ret_nlev, htop, hbot, ret_pres, $
                Fmatrix, Finv

  s = size(Fmatrix, /dimensions)
  nL = s[0]
  nj = s[1]

  ; -------------------------
  ; STEP 3: compute AK, Smoothing kernels
  ; -------------------------
  AKcoarse = ak
  Pcoarse = ret_pres[ak_pidx]

  AKfine = Fmatrix # AKcoarse # Finv
  Pfine = ret_pres[0:nL-1]

  Skernel = Fmatrix # Finv

end
