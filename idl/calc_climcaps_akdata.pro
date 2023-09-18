pro calc_climcaps_akdata, ret_pres, surf_pres, pres_nsurf, ak_pidx, $
			  htop, hbot, ak, ak_nfunc, ak_peff, $
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
  ;              The variable name in the CLIMCAPS L2 product is 'air_pres'
  ;              but should be converted from Pa to hPa.
  ; surf_pres    Scalar, surface pressure [hPa]
  ;              '/aux/prior_surf_pres' in L2 product, converted to hPa
  ; pres_nsurf   Scalar index specifying which element of ret_pres
  ;              is closest to surf_pres
  ;              'air_pres_lay_nsurf' in L2 product
  ; ak_pidx      The indices into the pressure grid for the trapezoid functions
  ;              'ave_kern/<mol_name>_func_indxs' in L2 product
  ; htop         Scalar integer specifying shape of TOA trapezoid function
  ;              'ave_kern/<mol_name>_func_htop' in L2 product
  ; hbot         Scalar integer specifying shape of near surface trapezoid function
  ;              'ave_kern/<mol_name>_func_hbot' in L2 product
  ; ak           Averaging kernel matrix for coarse trapezoid functions
  ;              'ave_kern/<mol_name>_ave_kern' in L2 product
  ; ak_nfunc     Scalar integer specifying the number of aks above surf_pres
  ;              'ave_kern/<mol_name>_func_last_indx' in L2 product
  ; ak_peff      Pressure value of coarse trapezoid layer, calculated as the log-average
  ;              pressure between two pressure levels
  ;              'ave_kern/<mol_name>_func_pres' in L2 product
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
  ; Prepare arrays before computing the F matrices.
  ; this will remove coarse layer(s) that are fully below the surface.
  ; The mid-layer pressure of the nearest-surface coarse layer is also adjusted,
  ; as well as the number of pressure levels in the fine grid.

  if keyword_set(adjust_surface) then begin
     ret_nlev = pres_nsurf
     ak_nlev = ak_nfunc

     ; note: AK and ak_pidx will be truncated if
     ; ak_nfunc < n_elements(diag_matrix(ak))
     ; This happens if the surface is above one entire coarse layer.
     AKcoarse = ak[0:ak_nlev-1,0:ak_nlev-1]

     ; note: ak_pidx contain the pressure level indices that form
     ; the boundaries of the coarse AK layers (aka trapezoids). 
     ; This means that n_elements(ak_pidx) = ak_nlev + 1
     ak_pidx = ak_pidx[0:ak_nlev]
     ; Replace the bottom index of the trapezoid to surf_pres
     ; with pres_nsurf
     ak_pidx[ak_nlev] = pres_nsurf
 
     ; Adjust bottom coarse AK pressure layer to surf_pres
     Pcoarse = ak_peff[0:ak_nlev-1]
     bot_pidx = ak_pidx[ak_nlev]
     top_pidx = ak_pidx[ak_nlev-1]
     bot_pdiff = ret_pres[bot_pidx-1] - ret_pres[top_pidx-1]
     Pcoarse[ak_nlev-1] = bot_pdiff/alog(ret_pres[bot_pidx-1]/ret_pres[top_pidx-1])

  endif else begin
     ; since no surface adjustment is made, use the full set of trapezoids.
     ak_nlev = n_elements(ak_pidx) - 1
     ret_nlev = n_elements(ret_pres)
     AKcoarse = ak
  endelse

  ; -------------------------
  ; STEP 2: call helper to compute Func matrix and inverse (F and F+)
  ; -------------------------
  calc_finv_mp, ak_nlev, ak_pidx, ret_nlev, htop, hbot, ret_pres, $
                Fmatrix, Finv

  s = size(Fmatrix, /dimensions)
  nL = s[0]
  nj = s[1]

  ; -------------------------
  ; STEP 3: compute AK, Smoothing kernels
  ; -------------------------
  AKfine = Fmatrix # AKcoarse # Finv
  Pfine = ret_pres[0:nL-1]

  Skernel = Fmatrix # Finv

end
