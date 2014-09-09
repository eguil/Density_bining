#!/usr/local/uvcdat/latest/bin/cdat
#
#
# -------------------------------------------------------------------------------------------------
# Procedure to compute water mass transformation from surface buoyancy fluxes in density space
#
#  Input fields:
#    - sst, sss, E-P, Qnet (3D, time,j,i)
#    - density grid sigrid (1D)
#    - target grid, including ocean basins
#  Output fields (on target grid):
#    - density flux (total, heat, fresh water) (2D rho,time, per basin or specific region)
#    - transformation (2D rho, time, per basin or specific region)
#
# Following Walin (1982) and Speer and Tziperman (1992)
# -------------------------------------------------------------------------------------------------
#   E. Guilyardi Sept 2014
#
#
# inits
# -----
#
def alpha (t, s):
    # compute alpha=-1/rho (d rho / d T)
    dt = 0.05
    siga = eso_neutral(t, s)
    sigb = eso_neutral(t+0.05, s)
    alpha = -0.001*(sigb-siga)/dt/(1.+1.e-3*siga)
    return alpha
def betar (t, s):
    # compute beta= 1/rho (d rho / d S)
    ds = 0.01
    siga = eos_neutral(t, s)-1000.
    sigb = eos_neutral(t, s+ds)-1000.    
    beta = 0.001*(sigb-siga)/ds/(1.+1.e-3*siga)
    return beta
def cpsw (t, s, p):
    # Specific heat of sea water (J/KG C)
    CP1 = 0.
    CP2 = 0.
    SR=SQRT(ABS(S))
    # SPECIFIC HEAT CP0 FOR P=0 (MILLERO ET AL. 1973)
    A = (-1.38E-3*T+0.10727)*T-7.644
    B = (5.35E-5*T-4.08E-3)*T+0.177
    C = (((2.093236E-5*T-2.654387E-3)*T+0.1412855)*T-3.720283)*T+4217.4
    CP0 = (B*SR + A) * S + C
    # CP1 PRESSURE AND TEMPERATURE TERMS FOR S = 0
    A = (((1.7168E-8*T+2.0357E-6)*T-3.13885E-4)*T+1.45747E-2)*T-0.49592
    B = (((2.2956E-11*T-4.0027E-9)*T+2.87533E-7)*T-1.08645E-5)*T+2.4931E-4
    C = ((6.136E-13*T-6.5637E-11)*T+2.6380E-9)*T-5.422E-8
    CP1 = ((C*P+B)*P+A)*P
    # CP2 PRESSURE AND TEMPERATURE TERMS FOR S > 0
    A = (((-2.9179E-10*T+2.5941E-8)*T+9.802E-7)*T-1.28315E-4)*T+4.9247E-3
    B = (3.122E-8*T-1.517E-6)*T-1.2331E-4
    A = (A+B*SR)*S
    B = ((1.8448E-11*T-2.3905E-9)*T+1.17054E-7)*T-2.9558E-6
    B = (B+9.971E-8*SR)*S
    C = (3.513E-13*T-1.7682E-11)*T+5.540E-10
    C = (C-1.4300E-12*T*SR)*S
    CP2 = ((C*P+B)*P+A)*P
    cp = CP0 + CP1 + CP2
    return cp    

def surface_transf(sst, sss, emo, qnet, area, sigrid, regrido, outgrid, masks):
    # Define dimensions
    N_i = int(sst.shape[2])
    N_j = int(sst.shape[1])
    N_t = int(sst.shape[0])
    N_s = len(sigrid)
    # Read masking value
    valmask = sst._FillValue
    # reorganise i,j dims in single dimension data
    sst  = npy.reshape(sst, (N_t, N_i*N_j))
    sss  = npy.reshape(sss, (N_t, N_i*N_j))
    emp  = npy.reshape(emp, (N_t, N_i*N_j))
    qnet = npy.reshape(qnet, (N_t, N_i*N_j))
    area = npy.reshape(area, (N_i*N_j))

    # Compute density
    rhon = eso_neutral(sst, sss)

    # Compute buoyancy flux as mass fluxes in kg/m2/s (SI unts)
    P = 0          # surface pressure
    conwf = 1.e-3  # kg/m2/s=mm/s -> m/s
    fheat = (-alpha(sst,sss)/cpsw(sst,sss,P))*qnet
    fwafl = (1000.+rhon)*beta(sst,sss)*sss*emp*convwf
    # find non-masked points
    maskin = mv.masked_values(sst.data[0], valmask).mask 
    nomask = npy.equal(maskin,0)
    # init arrays
    areabin = npy.ones((N_t,N_s))*valmask # surface of bin
    denflxh = npy.ones((N_t,N_s))*valmask # heat flux contrib
    denflxw = npy.ones((N_t,N_s))*valmask # E-P contrib

    # Bin on density grid
    for t in range(N_t):
        # bining loop
        for ks in range(N_s-1):
            # find indices of points in density bin
            idxbin = npy.argwhere( (rhon[t,:] >= sigrid[ks]) & (rhon[t,:] < sigrid[ks+1]) )
            denflxh[t,ks] = cdu.averager(fheat[t, idxbin] * area[idxbin], axis=1, action='sum')
            denflxw[t,ks] = cdu.averager(fwafl[t, idxbin] * area[idxbin], axis=1, action='sum')
            areabin[t,ks] = cdu.averager(area[idxbin], axis=1, action='sum')
        # last bin
        idxbin = npy.argwhere( (rhon[t,:] >= sigrid[N_s]))
        denflxh[t,N_s] = cdu.averager(fheat[t, idxbin] * area[idxbin], axis=1, action='sum')
        denflxw[t,N_s] = cdu.averager(fwafl[t, idxbin] * area[idxbin], axis=1, action='sum')
        areabin[t,N_s] = cdu.averager(area[idxbin], axis=1, action='sum')
        for ks in range(N_s-1):
        
        
        

   + create a basins variables (loop on n masks)

    return denflx, transf
