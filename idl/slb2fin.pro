pro slb2fin, numfunc, func_indx, func_ampl, usehalftop, usehalfbot,$
    presbot, air_pres, func_fine
; -------------------------------------------------
; PURPOSE: construct a trapezoid state function on 100 pressure levels
;          one at a time, using [ave_kern/*_func_indxs] hinge points
; -------------------------------------------------
;
; INPUT: 
;  NAME            DESCRIPTION
; ----------    -------------------------
; numfunc        number of trapezoid state functions above Earth surface at retrieval scene 
; func_indx      index values for trapezoid hinge points, starting at 1
; func_ampl      amplitude of trapezoids, which is 1.0 in this case
;                this is an array with numfunc elements, with the
;                trapezoid amplitude set for one function, and zero
;                elsewhere.
;                for example, with 5 functions, to compute the second
;                function, func_ampl would be [0,1,0,0,0].
; usehalftop     0 is a trapezoid, 1 is a wedge
; usehalfbot     0 is a trapezoid, 1 is a wedge
; presbot        pressure [hPa] of bottom retrieval level = 1100.0       
; air_pres       100-level retrieval pressure grid [air_pres]/100. in units [hPa]

; OUTPUT:
; func_fine      trapezoid state function on standard pressure level grid

; CALLED BY: calc_finv.pro
; --------------------------------------------------
; Original code by Chris D. Barnet and Eric S. Maddy
; chrisdbarnet@gmail.com
;
; Modified by Nadia Smith
; nadias@stcnet.com
; Science and Technology Corporation
; 6 July 2020
; ---------------------------------------------------
; 
; ---------------------------------------------------
; STEP 1: Construct the face of the trapezoid state function
; the 'face' values are amplitudes associated with the hinge points,
; therefore one additional value is needed.
; ---------------------------------------------------
  state_face = fltarr(numfunc+1)

  ; first hinge point value depends on whether the first function is a
  ; wedge or trapezoid depending on usehalftop: If wedge, use full amplitude,
  ; else use half.
  IF(usehalftop GT 0) THEN BEGIN
     state_face(0) = 0.5*func_ampl(0)
  ENDIF ELSE BEGIN
     state_face(0) = func_ampl(0)
  ENDELSE

  ; interior hinge point amplitudes are averaged
  FOR n = 1, numfunc - 1 DO BEGIN
     state_face(n) = 0.5*(func_ampl(n) + func_ampl(n - 1))
  ENDFOR

  ; last value depends on usehalfbot (analagous to first.)
  IF (usehalfbot gt 0) THEN BEGIN
     state_face(numfunc) = 0.5*func_ampl(numfunc - 1)
  ENDIF ELSE BEGIN
     state_face(numfunc) = func_ampl(numfunc - 1)
  ENDELSE

; ---------------------------------------------------
; STEP 2: Calculate the state function on standard 100 level grid
; the state function is a piecewise linear interpolation of the
; hinge point values (the state_face array) in log-pressure space.
; ---------------------------------------------------
  FOR n = 0, numfunc - 1 DO BEGIN
     idx_up = func_indx(n) - 1
     idx_down = func_indx(n+1) - 1

     state_up = state_face(n)
     state_down = state_face(n+1)

     pres_up = alog(air_pres(idx_up))
     pres_down = alog(air_pres(idx_down))

     slope = (state_down - state_up)/(pres_down - pres_up)

     FOR L = idx_up, idx_down - 1 DO BEGIN
        func_fine(L) = state_up + slope * (alog(air_pres(L)) - pres_up)
     ENDFOR
  ENDFOR
  func_fine(idx_down) = state_down

END
