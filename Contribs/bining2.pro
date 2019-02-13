;+
;
;  PRINCIPE DU BINING
;
;  - etant donnee une liste de densite (sigma) on cherche les  profondeurs
;    correspondants a ces valeurs. Pour cela on utilise une interpolation lineaire.
;
;  - puis, connaissant les profondeurs des isopycnes, on interpole les valeurs du contenu
;    vertical de x a ces profondeurs. Il est preferable d' interpoler le contenu de x plutot
;    que x afin d'assurer la conservation de x.
;
;  - on deduit la valeur moyenne de x dans chaque couche [i.e. entre deux isopycnes). Il ne
;    s'agit donc pas de la projection de x selon les isopycnes mais de la valeur moyenne entre
;    deux isopycnes. Les proprietes sont meilleures (conservation de la quantite totale de x).
;
; REMARQUES
;
;  - il n'y a aucun appel a un fichier de common (d'ou la liste des mots clef un peu longue)
;
;  - le calcul n'est rigoureux que pour les points T. Le bining de la vitesse doit etre adapte...
;;
;  - si l'utilisateur donne N interfaces, cela defini N+1 couches, d'ou les decalages dans les indices
;    et l'apparition d'un +1 dans la declaration des tableaux de sortie
;    la couche 0   contient tous les points dont la densite est inferieure a la densite minimale du bining
;    la couche N+1  "                   "                     " superieure      "       maximale
;
; @param DENSITY {in}
; density at T points (3D array)
;
; @param X1 {in}
; field at T points (e.g. Salinity or Temperature)
;
; @param SOBWLMAX {in}{type=2D array}
; bowl depth array
;
; @param SIG_BOWL {in}{type=integer}
; switch for bowl overlay
;
; @param DEPTH_BIN {out}
; depth of layers (3D array)
;
; @param THICK_BIN {out}
; thickness of layers (3D array)
;
; @param X1_BIN {out}
; X averaged for each sigma-layer (3D array)
;
; @param BOWL_BIN {out}
; bowl depth binned on density
;
; @keyword SIGMA {optional}
; bining values
;
; @keyword DEPTH_T {required}
; depth of T level
;
; @keyword DEPTH_W {required}
; depth of W level
;
; @keyword E3T
;
; @keyword E3W
;
; @keyword TMASK {required}
; tmask (3D array)
;
; @history
; - fplod 20100119T160644Z aedon.locean-ipsl.upmc.fr (Darwin)
;
;   * check parameters
;
; - fplod 20091208T102329Z aedon.locean-ipsl.upmc.fr (Darwin)
;
;   * syntax of array
;
; - 19/10/99 G. Roullet
;
; @version
; $Id: bining2.pro 206 2010-01-26 10:33:28Z pinsard $
;
;-
PRO bining2 $
    , density $
    , x1 $
    , sobwlmax $
    , sig_bowl $
    , depth_bin $
    , thick_bin $
    , x1_bin $
    , bowl_bin $
    , SIGMA=sigma $
    , DEPTH_T=depth_t $
    , DEPTH_W=depth_w $
    , E3T=e3t $
    , E3W=e3w $
    , TMASK=tmask
;
  compile_opt idl2, strictarrsubs
;

 usage='bining2' $
    + ', density' $
    + ', x1' $
    + ', sobwlmax' $
    + ', sig_bowl' $
    + ', depth_bin' $
    + ', thick_bin' $
    + ', x1_bin' $
    + ', bowl_bin' $
    + ', SIGMA=sigma' $
    + ', DEPTH_T=depth_t' $
    + ', DEPTH_W=depth_w' $
    + ', E3T=e3t' $
    + ', E3W=e3w' $
    + ', TMASK=tmask'
;
 nparam = N_PARAMS()
 IF (nparam LT 8) THEN BEGIN
    ras = report(['Incorrect number of arguments.' $
          + '!C' $
          + 'Usage : ' + usage])
    stop
 ENDIF

 arg_type = size(density,/type)
 IF (arg_type NE 4) THEN BEGIN
   ras = report(['Incorrect arg type density' $
          + '!C' $
          + 'Usage : ' + usage])
   stop
 ENDIF

 arg_type = size(sobwlmax,/type)
 IF (arg_type NE 4) THEN BEGIN
   ras = report(['Incorrect arg type sobwlmax' $
          + '!C' $
          + 'Usage : ' + usage])
   stop
 ENDIF

  size3d = size(density)
   jpi = size3d[1]
   jpj = size3d[2]
   jpk = size3d[3]
;print, 'size3d=', size3d
;print, 'size(tmask)=', size(tmask)

   indm = where(tmask eq 0)

   x1[indm] = !VALUES.F_NAN
   density[indm] = !VALUES.F_NAN

   N_z = jpk

; adjust profiles to avoid static instabilities (there are almost none IF using pac_40)
IF jpi EQ 180 THEN BEGIN
   density = npc(density)
ENDIF

;  N_s, defini ci dessous est le nombre d'interface

   IF keyword_set(sigma) THEN BEGIN
      s_s = sigma
      N_s = n_elements(s_s)
   ENDIF ELSE BEGIN
;
;  Definition d'un bining arbitraire (si density est un sigma potentiel, cela repere
;   les eaux intermediaires)
;
      N_s = 20
      s_s = 27.6+0.2*findgen(N_s)/(N_s-1)
      N_s = 3
      s_s = [26.8, 27.6, 27.8]
      s_s = [22.8, 27.6, 27.8]
   ENDELSE

;   print, 'Density bining :', s_s
;;
;; Definition des variables
;;
;   Profils selon les niveaux du modele (suffixe _z)
;
   c1_z = fltarr(N_z)    ; profil du contenu vertical de x
   s_z = fltarr(N_z)    ; profil de la densite
   z_zt = depth_t  ; profondeur au point T  (k=0 -> 5m)
   z_zw = depth_w  ;                     W  (k=0 -> 10m)

;   Profils selon les couches de densite (suffixe _s)
;
   z_s = fltarr(N_s+1)
   c1_s = fltarr(N_s+1)
   x1_s = fltarr(N_s+1)
   bowl_s = 0.

; Tableaux de sorties
;
   depth_bin = fltarr(jpi, jpj, N_s+1)
   thick_bin = fltarr(jpi, jpj, N_s+1)
   x1_bin = fltarr(jpi, jpj, N_s+1)
   bowl_bin = fltarr(jpi, jpj)

   x1_bin[*, *, *] = !VALUES.F_NAN
   depth_bin[*, *, *] = !VALUES.F_NAN
   thick_bin[*, *, *] = !VALUES.F_NAN
   bowl_bin[*, *] = !VALUES.F_NAN


;  Calcul du contenu de x sur la verticale (permet d assurer la conservation integrale)
;
   x1_content = fltarr(jpi, jpj, jpk)


;   x1_content[*, *, 0] = e3t[0] * x1[*, *, 0] * tmask [*, *, 0]
;   FOR k = 1, jpk-1 DO x1_content[*, *, k] = x1_content[*, *, k-1] $
;                        + e3t[k] * x1[*, *, k]* tmask [*, *, k]
   x1_content = x1

;
;  Boucle sur le domaine 2D
;

   FOR i = 0, (jpi-1) DO BEGIN
;  FOR i = 20, 20 DO BEGIN
;   FOR i = 162, 162 DO BEGIN
      FOR j = 0, (jpj-1) DO BEGIN
;     FOR j = 20, 20 DO BEGIN
;      FOR j = 116, 116 DO BEGIN
;
;  Indices des points T dans l ocean
;print, '    '
;print, ' ==============================='
;print, '    POINT I,J =', i, j
;
;stop
         i_ocean = where(tmask[i, j, *] EQ 1)
         z_s = z_s*0.
         c1_s[*] = !VALUES.F_NAN
         bowl_s = !VALUES.F_NAN
;         x1_s[*] = !VALUES.F_NAN


         IF i_ocean[0] NE -1 THEN BEGIN
; on n entre que si il y a des points ocean
;            print, 'ocean point'
            i_bottom = i_ocean[n_elements[i_ocean]-1]
;print, 'i_bottom=', i_bottom
            z_s[N_s] = z_zw[i_bottom]
            c1_s[N_s] = x1_content[i, j, jpk-1]

            s_z[*] = density[i, j, *]
            c1_z[*] = x1_content[i, j, *]
;            print, 'density profile s_z', s_z
;            print, 'field profile c1_z', c1_z

;            print, 'bins s_s', s_s

; extraction d'un sous profil strictement croissant
;            print, 's_z[i_ocean]=', s_z[i_ocean]
            mini = min(s_z[i_ocean])
            maxi = max(s_z[i_ocean])
            i_min = where(s_z[i_ocean] EQ mini)
            i_max = where(s_z[i_ocean] EQ maxi)
;            print, 'mini, maxi', mini, maxi
;            print, 'i_min, i_max', i_min, i_max
;   on prend le plus grand des indices min
;         et le plus petit des indices max
            i_min = i_min(n_elements(i_min)-1)
            i_max = i_max[0]
;            print, 'i_min, i_max', i_min, i_max
;IF i_min GT i_max THEN BEGIN
; print, 'WARNING: i_min, i_max=', i_min, i_max, ' at i,j=,', i, j
;ENDIF
;   on prend le plus grand des indices min

;            IF i_max GE jpk-1 THEN BEGIN
;             print, i, j, i_max
;            ENDIF
;            IF i_min LE 1 THEN BEGIN
;             print, i, j, i_min
;            ENDIF

; Si la valeur du niveau (s_s) est plus faible que la densite de surface,
; l isopycne est mise en surface (z_s=0)
;
;            print, 's_z[i_min]', s_z[i_min]
;print, 's_s=', s_s
            ind = where(s_s LT s_z[i_min])
            IF ind[0] NE -1 THEN BEGIN
;               IF i_min GT i_max THEN BEGIN
;                print, 'min reached at sigma indices', ind
;               ENDIF
               z_s[ind] = 0
               c1_s[ind] = !VALUES.F_NAN
            ENDIF

; Si la valeur du niveau (s_s) est plus elevee que la densite du fond,
; l isopycne est mise au fond (z_s=z_zw[i_bottom])
;
;            print, 's_z[i_max]', s_z[i_max]
            ind = where(s_s GT s_z[i_max])

            IF ind[0] NE -1 THEN BEGIN
;               IF i_min GT i_max THEN BEGIN
;                print, 'max reached at sigma indices', ind
;               ENDIF
               z_s[ind] = z_s[N_s]
               c1_s[ind] = c1_s[N_s]
               c1_s[ind] = !VALUES.F_NAN
            ENDIF
; cas general
            ind = where( (s_s GE s_z[i_min]) AND (s_s LE s_z[i_max]) )
;            IF i_min GT i_max THEN BEGIN
;             print, 's_s indices inside min/max', ind
;            ENDIF
            IF ind[0] NE -1 THEN BEGIN
               i_profil = i_ocean[i_min:i_max] ; problem line

               z_s[ind] = interpol(z_zt[i_profil], s_z[i_profil], s_s[ind]) ; original
;               z_s[ind] = interpol(z_zw[i_profil], s_z[i_profil], s_s[ind]) ; changed to z_zw 1/7/04

;
; j'utilise la fonction spline pour interpoler le contenu
; l'interpolation lineaire marche aussi. Elle donne des resultats differents. Elle
; me semble moins propre. Je la donne en remarque.
;
               IF n_elements[i_profil] GT 2 THEN BEGIN
;                  c1_s[ind] = spline(z_zw[i_profil], c1_z[i_profil], z_s[ind], 1)     ; cubic spline
                  c1_s[ind] = interpol(c1_z[i_profil], z_zt[i_profil], z_s[ind])      ; linear interpolation
;                   SPLINE_P, z_zw[i_profil], c1_z[i_profil], z_s[ind], c1_s[ind]     ; generalization of spline(*)
;print, z_zw[i_profil], c1_z[i_profil], z_s[ind], c1_s[ind]
               ENDIF ELSE BEGIN
                  c1_s[ind] = interpol(c1_z[i_profil], z_zt[i_profil], z_s[ind])
               ENDELSE
;
;     bowl depth binning

               IF sig_bowl EQ 1 THEN BEGIN
                  bowl_s = interpol(s_z[i_profil], z_zt[i_profil], sobwlmax[i, j])
               ENDIF ELSE BEGIN
                  bowl_s = 0
               ENDELSE

;stop
;
;
;               x1_s[ind] = (c1_s[ind+1]-c1_s[ind])/(z_s[ind+1]-z_s[ind])
;                     x1_s = c1_s

            ENDIF
;ELSE print, 'ind[0] = -1', ind, i_bottom
         ENDIF
;ELSE print, ' land ', tmask[i, j, 1]
         depth_bin [i, j, *] = z_s
         thick_bin [i, j, 0] = z_s[0]
         thick_bin [i, j, 1:N_s] = z_s[1:N_s]-z_s[0:N_s-1]
         IF total(thick_bin [i, j, 1:N_s] LT 0) GT 0 THEN BEGIN
            ;print, 'WARNING: negative thick_bin [i, j, 1:N_s] at i,j= ', i, j
            ;print, 'depth_bin [i, j, 1:N_s]= ', depth_bin [i, j, 1:N_s]
            ;print, 'thick_bin [i, j, 1:N_s]= ', thick_bin [i, j, 1:N_s]
            ;print, 's_z= ', s_z
            ;stop
         ENDIF
         x1_bin    [i, j, *] = c1_s
         bowl_bin  [i, j] = bowl_s
;         print, 'c1_s', c1_s
;
      ENDFOR
   ENDFOR
; mask depth_bin with undefined values of x_bin
   depth_bin[where(finite(x1_bin, /nan) EQ 1)] = !VALUES.F_NAN
   thick_bin[where(finite(x1_bin, /nan) EQ 1)] = !VALUES.F_NAN

; OPTIONAL: compute spiciness by removing isopycnal mean for each sigma (T or S)
   IF 0 THEN BEGIN
      print, '    Computing spiciness...'
      FOR i = 0, (size(x1_bin))[3]-1 DO BEGIN
         x1_bin[*,*,i] = x1_bin[*,*,i] - mean(x1_bin[*,*,i],/NAN)
      ENDFOR
   ENDIF

; OPTIONAL: remove domain mean (S)
   IF 0 THEN BEGIN
      print, '    Removing domain mean...'
      x1_bin[*,*,*] = x1_bin[*,*,*] - mean(x1_bin[*,*,*],/NAN)
   ENDIF

; OPTIONAL: remove 34.6psu (S)
   IF 1 THEN BEGIN
      print, '    Removing 34.6psu...'
      x1_bin[*,*,*] = x1_bin[*,*,*] - 34.6
   ENDIF

END
