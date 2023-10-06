pro write_idl_array, fid, x, varname
  ; helper to write a single idl array
  ; following https://www.nv5geospatialsoftware.com/docs/HDF5_Overview.html
  datatype_id = h5t_idl_create(x)
  dataspace_id = h5s_create_simple(size(x,/dimensions))
  dataset_id = h5d_create(fid, varname, datatype_id, dataspace_id)
  h5d_write, dataset_id, x
  h5d_close, dataset_id
  h5s_close, dataspace_id
  h5t_close, datatype_id
end

pro write_ak_run_output, h5_output_file, Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel
  ; helper to write the test data into h5 file
  ; following https://www.nv5geospatialsoftware.com/docs/HDF5_Overview.html
  fid = h5f_create(h5_output_file)
  write_idl_array, fid, Fmatrix, 'Fmatrix'
  write_idl_array, fid, Finv, 'Finv'
  write_idl_array, fid, AKcoarse, 'AKcoarse'
  write_idl_array, fid, Pcoarse, 'Pcoarse'
  write_idl_array, fid, AKfine, 'AKfine'
  write_idl_array, fid, Pfine, 'Pfine'
  write_idl_array, fid, Skernel, 'Skernel'
  h5f_close, fid
end

function read_idl_array, fid, varname
  ; helper to read a single idl array
  ; following https://www.nv5geospatialsoftware.com/docs/HDF5_Overview.html
  dataset_id = h5d_open(fid, varname)
  x = h5d_read(dataset_id)
  h5d_close, dataset_id
  return, x
end

pro test_ak_run_output, h5_output_file, Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel

  fid = h5f_open(h5_output_file)

  Fmatrix_stored = read_idl_array(fid, 'Fmatrix')
  print, 'Fmatrix max difference:   ', max(abs(Fmatrix - Fmatrix_stored))
  Finv_stored = read_idl_array(fid, 'Finv')
  print, 'Finv max difference:      ', max(abs(Finv - Finv_stored))

  AKcoarse_stored = read_idl_array(fid, 'AKcoarse')
  print, 'AKcoarse max difference:  ', max(abs(AKcoarse - AKcoarse_stored))
  Pcoarse_stored = read_idl_array(fid, 'Pcoarse')
  print, 'Pcoarse max difference:   ', max(abs(Pcoarse - Pcoarse_stored))

  AKfine_stored = read_idl_array(fid, 'AKfine')
  print, 'AKfine max difference:    ', max(abs(AKfine - AKfine_stored))
  Pfine_stored = read_idl_array(fid, 'Pfine')
  print, 'Pfine max difference:     ', max(abs(Pfine - Pfine_stored))

  Skernel_stored = read_idl_array(fid, 'Skernel')
  print, 'Skernel max difference:   ', max(abs(Skernel - Skernel_stored))

  h5f_close, fid

end

pro sample_ak_run, molnames, write_h5_output=write_h5_output, test_h5_output=test_h5_output

; for each of these two files, there is a hand-picked FOR.
; these two files should be downloaded with the data_download.sh file
; at the top level.

ncfiles = [ $
          '../test_data/SNDR.J1.CRIMSS.20190901T0336.m06.g037.L2_CLIMCAPS_RET.std.v02_28.G.200214174949.nc', $
          '../test_data/SNDR.J1.CRIMSS.20190901T2248.m06.g229.L2_CLIMCAPS_RET.std.v02_28.G.200214190447.nc' ]
ifoot = [3, 4]
iscan = [14, 32]

for m=0, n_elements(molnames)-1 do begin
   for i=0, n_elements(ncfiles)-1 do begin

      return_climcaps_akdata, $
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

      h5_output_file = string(i+1, molnames[m], format='("../test_data/test_case_",(1I1),"_",(A),".h5")')
      if keyword_set(write_h5_output) then begin
         print, 'writing output file: ', h5_output_file
         write_ak_run_output, $
            h5_output_file, $
            Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel
      endif
      if keyword_set(test_h5_output) then begin
         print, 'testing output relative to: ', h5_output_file
         test_ak_run_output, $
            h5_output_file, $
            Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel
      endif

   endfor
endfor

end
