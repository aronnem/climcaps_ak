pro plot_fmatrix, Fmatrix, Finv, Plev, plot_filename
  ;
  ; plot F matrix (and pseudoinverse) in a 2x1 subplot
  ;
  nj = (size(Fmatrix))[2]
  yrange = [1200.0, 0.1]
  colors = colortable([[204,204,0],[0,144,255],[0,0,53]],ncolors=nj)

  a = plot(Fmatrix[*,0], Plev, $
           /ylog, xrange=[-0.1,1.1], yrange=yrange, title='Func matrix', $
           ytitle='Pressure [hPa]', layout=[2,1,1], /buffer, $
           color=reform(colors[0,*]))
  a.ytickfont_size = 8
  a.xtickfont_size = 8

  for j = 1,nj-1 do begin
     a = plot(Fmatrix[*,j], Plev, /overplot, /current, color=reform(colors[j,*]))
  endfor

  a = plot(Finv[0,*], Plev, $
           /ylog, xrange=[-1.1,1.1], yrange=yrange, title='Func Inv', $
           ytitle='Pressure [hPa]', layout=[2,1,2], /buffer, $
           /current, color=reform(colors[0,*]))
  a.ytickfont_size = 8
  a.xtickfont_size = 8

  for j = 1,nj-1 do begin
     a = plot(Finv[j,*], Plev, /overplot, /current, color=reform(colors[j,*]))
  endfor

  a.save, plot_filename

end


pro plot_ak, ak, plev, plot_filename
  ;
  ; plot averaging kernels
  ;
  nj = (size(ak,/dimension))[1]
  yrange = [1200.0, 0.1]
  colors = colortable([[204,204,0],[0,144,255],[0,0,53]],ncolors=nj)
  ; round to nearest 0.1
  amax = (max(ak*10) + 1)/10.0
  amin = (min(ak*10) + 1)/10.0

  a = plot(ak[*,0], Plev, $
           /ylog, xrange=[amin,amax], yrange=yrange, title='Averaging Kernel', $
           ytitle='Pressure [hPa]', /buffer, $
           color=reform(colors[0,*]), dimensions=[400,800])
  a.ytickfont_size = 8
  a.xtickfont_size = 8

  for j = 1,nj-1 do begin
     a = plot(ak[*,j], Plev, /overplot, /current, color=reform(colors[j,*]))
  endfor

  a.save, plot_filename

end



pro plot_sfunc, sfunc, plev, plot_filename
  ;
  ; plot smoothing functions
  ;
  nj = (size(sfunc,/dimension))[0]
  yrange = [1200.0, 0.1]
  colors = colortable([[204,204,0],[0,144,255],[0,0,53]],ncolors=nj)
  ; round to nearest 0.1
  smax = (max(sfunc*10) + 1)/10.0

  a = plot(sfunc[0,*], Plev, $
           /ylog, xrange=[-0.1,smax], yrange=yrange, title='Smoothing Functions', $
           ytitle='Pressure [hPa]', /buffer, $
           color=reform(colors[0,*]), dimensions=[400,800])
  a.ytickfont_size = 8
  a.xtickfont_size = 8

  for j = 1,nj-1 do begin
     a = plot(sfunc[j,*], Plev, /overplot, /current, color=reform(colors[j,*]))
  endfor

  a.save, plot_filename

end
