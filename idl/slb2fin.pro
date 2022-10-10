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
; usehalftop     0 is a trapezoid, 1 is a wedge
; usehalfbot     0 is a trapezoid, 1 is a wedge
; presbot        pressure [hPa] of bottom retrieval level = 1100.0       
; air_pres        100-level retrieval pressure grid [air_pres]/100. in units [hPa]

; OUTPUT:
; func_fine    trapezoid state function on standard pressure level grid

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
; STEP 1: Construct the face of the trapezoide state function
; ---------------------------------------------------
     	state_face = fltarr(numfunc)

      IF(usehalftop GT 0) THEN BEGIN
        	state_face(0) = 0.5*func_ampl(0)
      ENDIF ELSE BEGIN
        	state_face(0) = func_ampl(0)
      ENDELSE

      FOR n = 1, numfunc - 2 DO BEGIN
        	state_face(n) = 0.5*(func_ampl(n) + func_ampl(n - 1))
      ENDFOR

      IF (usehalfbot gt 0) THEN BEGIN
        	state_face(numfunc - 1) = 0.5*func_ampl(numfunc - 2)
      ENDIF ELSE BEGIN
      	state_face(numfunc - 1) = func_ampl(numfunc - 2)
      ENDELSE
; ---------------------------------------------------
; STEP 2: Calculate the state function on standard 100 level grid
; ---------------------------------------------------
		FOR n = 0, numfunc - 2 DO BEGIN
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
