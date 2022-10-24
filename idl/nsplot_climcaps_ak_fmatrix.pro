pro nsplot_climcaps_ak_fmatrix,type,func_matrix,func_inv,ak,ak100

; testing calc_finv.pro output using ozone
;fn='/peate_sci/nadias/data/L2/jpss/2020/SNDR.J1.CRIMSS.20201230T2248.m06.g229.L2_CLIMCAPS_RET.std.v02_28.G.210203073025.nc'
fn='../SNDR.J1.CRIMSS.20201230T2248.m06.g229.L2_CLIMCAPS_RET.std.v02_28.G.210203073025.nc'

; Read CLIMCAPS netcdf file for averaging kernels and CO column density profiles
ncread_climcaps_main,fn,clim

; Read AK group
nsread_nc_group,fn,'ave_kern',akgrp

; Read Aux gropu
nsread_nc_group,fn,'aux',auxgrp

fillval=9.96921e+36

htop = akgrp.o3_func_htop
hbot = akgrp.o3_func_hbot

nfoot=30; number of retrieval footprints per scanline
nscan=45; number of scanlines per granule

; Surface pressure
surf_pres=auxgrp.prior_surf_pres/100.; Pa to hPa

ifoot = 15
iscan = 20

; Retreival pressure layers
; loop over all retrieval scenes
	ret_nlev=100;clim.air_pres_nsurf(ifoot,iscan)
	psurf=surf_pres(ifoot,iscan)
	ret_pres=clim.air_pres/100.;(0:ret_nlev)/100.
	ak_pidx=akgrp.o3_func_indxs; 10 indices 
	ak_nlev=n_elements(ak_pidx)
   ak_plev=akgrp.o3_func_pres/100.

; -------------------------
; STEP 2: adjust surface
; -------------------------

; Call changlev to determine number of pressure levels above surface pressure at retrieval scene
; ;  Name          Description
; ------    ----------------------
; psurf     surface pressure
; ret_pres      Standard RTA pressure grid
; ret_nlev    number of levels (RTA grid)
; ak_nlev    number of coarse layers to be adjusted for topography
; ak_pres     coarse layers to be adjusted for topography 

; OUTPUT
; ak_nlev_scene    number of coarse levels at scene psurf
; ak_pidx_scene    coarse level boundaries at scene psurf
	changlev, psurf, ret_pres, ret_nlev, ak_nlev, ak_pidx, ak_nlev_scene, ak_pidx_scene

; ---------------------------------------------------------
; STEP 3: Calculate func_matrix and its inverse func_inv
; ---------------------------------------------------------
 
; Now compute the transformation matrix and its pseudo inverse 
; Equation 11 of Maddy and Barnet, IEEE TGRS, 2008 (MB08) 
; func_matrix is the [ak_nlev x ret_nlev] matrix of trapezoids
; func_inv is the inverse of func_matrix [ret_nlev x ak_nlev]
	calc_finv_mp, ak_nlev_scene, ak_pidx_scene, ret_nlev, htop, hbot, ret_pres,$
		func_matrix, func_inv
print,ak_pidx_scene
	s=size(func_matrix)
; Swtich to using L and j indices as in Maddy & Barnet 2008
	nL = s(1); ret_nlev
	nj = s(2); ak_nlev
; func_matrix [nL x nj]
; func_inv [nj x nL]
; ----------------------------------------------------------
; STEP 4: Plot averaging kernels f-matrix and its inverse
; ----------------------------------------------------------
	plev=ret_pres(0:nL-1)
	yrange=[1000.,10.]
	yticks=[1000,500,300,100,50,10]
	ytickname=['1000.','500.','300.','100.','50.','10.'];
	ct=colortable([[204,204,0],[0,144,255],[0,0,53]],ncolors=12)
	ct100=colortable([[204,204,0],[0,144,255],[0,0,53],[252,186,3],[252,3,173]],ncolors=100)
; 
; ****** 1) Type 1: Plot func_matrix and its inverse 
	if type eq 1 then begin
; Panel 1: func_matrix (referred to as "F" in Maddy&Barnet 2008)
		a=plot(func_matrix(*,0),plev,/ylog,layout=[2,1,1],/buffer,font_size=8,xrange=[0,1],$
		thick=1,color=reform(ct(0,*)),yrange=yrange,title='func_matrix',ytitle='Pressure [hPa]')
   		a.xtickfont_size=8
        	a.ytickfont_size=8
			a.ytickvalues=yticks
; loop over state functions
		for i=1,nj-1 do begin
			a=plot(func_matrix(*,i),plev,color=reform(ct100(i,*)),/overplot,/current)
		endfor

; Panel 2: func_inv (referred to as "F+" in Maddy&Barnet 2008)
		a=plot(func_inv(0,*),plev,/ylog,layout=[2,1,2],/buffer,font_size=8,/current,$
   	thick=1,color=reform(ct(0,*)),yrange=yrange,title='func_inv',ytitle='Pressure [hPa]')
   		a.xtickfont_size=8
        	a.ytickfont_size=8
                a.ysubgridstyle=0
         a.ytickvalues=yticks
; loop over state functions
		for i=1,nj-1 do begin
   		a=plot(func_inv(i,*),plev,color=reform(ct(i,*)),/overplot,/current)
		endfor
; save as jpg
		jname='test_fmatrix_finv.jpg'
		a.save,jname
; write arrays to csv files
		cname1='climcaps_f_matrix.csv'
		write_csv,cname1,func_matrix
		cname2='climcaps_finv.csv'
		write_csv,cname2,func_inv
		print,'saved '+jname+', '+cname1+' and '+cname2
	endif 
; ****** 2) Type 2: Plot averaging kernels on coarse layers, i.e., the AKs are reported in netcdf file
	if type eq 2 then begin
		ak_pres=ret_pres(ak_pidx_scene)
		print,ak_pres
      ak=reform(akgrp.o3_ave_kern(0:nj-1,0:nj-1,ifoot,iscan)); [nj x nj]
      a=plot(ak(0,*),ak_pres,/ylog,/buffer,font_size=8,$
       dimensions=[400,800],thick=1,color=reform(ct(0,*)),yrange=yrange,title='Coarse Averaging Kernels (AKs)',ytitle='Pressure [hPa]')
         a.xtickfont_size=8
         a.ytickfont_size=8
         a.ytickvalues=yticks
		for i=1,nj-1 do begin
         a=plot(ak(i,*),ak_pres,color=reform(ct(i,*)),/overplot,/current)
      endfor
; save as jpg
      jname='test_o3_coarse_ak.jpg'
      a.save,jname
; write array to csv file
		cname='climcaps_o3_coarse_aks.csv'
		write_csv,cname,ak
		print,'saved '+jname+' and '+cname
	endif
; ****** 3) Type 3: Plot averaging kernels on fine layers (referred to as the "effective AKs"
; in Maddy&Barnet 2008 and calculated as F*AK*F+) [nL x nL]
	if type eq 3 then begin
		ak=reform(akgrp.o3_ave_kern(0:nj-1,0:nj-1,ifoot,iscan)); [nj x nj]
; # denotes a matrix multiplication
; alternatively use IDL function "matrix_multiply"
		ak100=func_matrix#ak#func_inv; [nL x nj]*[nj x nj]*[nj x nL] = [nL x nL]

		a=plot(ak100(*,0),plev,/ylog,/buffer,font_size=8,$
  		 dimensions=[400,800],thick=1,color=reform(ct100(0,*)),yrange=yrange,title='Effective AKs (F#AK#F+)',ytitle='Pressure [hPa]')
         a.xtickfont_size=8
         a.ytickfont_size=8
         a.ytickvalues=yticks
		for i=1,nL-1 do begin
			a=plot(ak100(*,i),plev,color=reform(ct100(i,*)),/overplot,/current)
		endfor
; save as jpg
		jname='test_o3_fine_ak.jpg'
		a.save,jname
; write array to cvs file
		cname='climcaps_o3_effective_aks.csv'
		write_csv,cname,ak100
		print,'saved '+jname+' and '+cname
	endif 
; ****** 4) Type 4: Plot smoothing kernels (calculated as "FF+" in Maddy&Barnet 2008) 
	if type eq 4 then begin
		sfunc = func_matrix#func_inv; [nL x nj]*[nj x nL] = [nL x nL]
		a=plot(sfunc(0,*),plev,/ylog,/buffer,font_size=8,$
       dimensions=[400,800],thick=1,color=reform(ct100(0,*)),yrange=yrange,title='Smoothing Kernels (F#F+)',ytitle='Pressure [hPa]')
         a.xtickfont_size=8
         a.ytickfont_size=8
         a.ytickvalues=yticks
		for i=1,nL-1 do begin
         a=plot(sfunc(i,*),plev,color=reform(ct100(i,*)),/overplot,/current)
      endfor
; save as jpg 
      jname='test_o3_smooth_kernels.jpg'
      a.save,jname
; write to cvs file
		write_csv,'climcaps_o3_smooth_kernels.csv',sfunc
      print,'saved '+jname,' and climcaps_o3_smooth_kernels.csv'
	endif

end

