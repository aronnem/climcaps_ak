pro sample_ak_run

;
; batch run, over all profile: air temp, and five molecules.
; skipping hno3 since that does not have valid data in both sample obs.
;

; create AKs for each of the molecules and air temp.
molnames = ['air_temp', 'h2o_vap', 'o3', 'co', 'co2', 'ch4']

; for each of these two files, there is a hand-picked FOR.
; these two files should be downloaded with the data_download.sh file
; at the top level.

ncfiles = [ $
          '../SNDR.J1.CRIMSS.20190901T0336.m06.g037.L2_CLIMCAPS_RET.std.v02_28.G.200214174949.nc', $
          '../SNDR.J1.CRIMSS.20190901T2248.m06.g229.L2_CLIMCAPS_RET.std.v02_28.G.200214190447.nc' ]
ifoot = [ 4, 15]
iscan = [32, 20]

for m=0, n_elements(molnames)-1 do begin
   for i=0, n_elements(ncfiles)-1 do begin

      read_climcaps_akdata, $
         ncfiles[i], ifoot[i], iscan[i], molnames[m], $
         Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel

      outname = 'for'+strtrim(ifoot[i],1)+'_scan'+strtrim(iscan[i],1)+'_'+molnames[m]

      plot_fmatrix, Fmatrix, Finv, Pfine, 'example_plot_Fmatrix_'+outname+'.png'
      plot_ak, AKcoarse, Pcoarse, 'example_plot_AKcoarse_'+outname+'.png'
      plot_ak, AKfine, Pfine, 'example_plot_AKfine_'+outname+'.png'
      plot_sfunc, Skernel, Pfine, 'example_plot_Sfunc_'+outname+'.png'

      write_csv, 'example_Fmatrix_'+outname+'.csv', Fmatrix
      write_csv, 'example_Finv_'+outname+'.csv', Finv
      write_csv, 'example_AKfine_'+outname+'.csv', AKfine
      write_csv, 'example_AKcoarse_'+outname+'.csv', AKcoarse
      write_csv, 'example_Skernel_'+outname+'.csv', Skernel

   endfor
endfor

end
