#
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ctypes import *
plt.rcParams["font.size"]=15

#
DX_tmp       = 1000
SFC_TEMP_tmp = 300
ATM_TEMP_tmp = SFC_TEMP_tmp
SFC_PRES_tmp = 1e5
ATM_PRES_tmp = SFC_PRES_tmp
ATM_QV_tmp   = 2e-2
SFC_Z0M_tmp  = 1e-4
SFC_Z0H_tmp  = 1e-4
SFC_Z0E_tmp  = 1e-4
ATM_U_tmp    = np.arange(-4,4,1).reshape((2,2,2))
ATM_V_tmp    = -np.arange(-4,4,1).reshape((2,2,2))
ATM_Z1_tmp   = 50
PBL_tmp      = 1000
SFC_DENS_tmp = 1.2

#
IA = int(ATM_U_tmp.shape[2])
IS = int(1)
IE = int(IA)
JA = int(ATM_U_tmp.shape[1])
JS = int(1)
JE = int(JA)
TA = int(ATM_U_tmp.shape[0])
TS = int(1)
TE = int(TA)
DX = np.array([DX_tmp]).astype(np.float64)
DY = DX
ATM_TEMP = (ATM_TEMP_tmp*np.ones((TA,JA,IA))).astype(np.float64)
SFC_TEMP = (SFC_TEMP_tmp*np.ones((TA,JA,IA))).astype(np.float64)
ATM_PRES = (ATM_PRES_tmp*np.ones((TA,JA,IA))).astype(np.float64)
SFC_PRES = (SFC_PRES_tmp*np.ones((TA,JA,IA))).astype(np.float64)
ATM_QV   = (ATM_QV_tmp*np.ones((TA,JA,IA))).astype(np.float64)
ATM_U    = ATM_U_tmp.astype(np.float64)
ATM_V    = ATM_V_tmp.astype(np.float64)
ATM_Z1   = (ATM_Z1_tmp    *np.ones((JA,IA))).astype(np.float64)
PBL      = (PBL_tmp*np.ones((TA,JA,IA))).astype(np.float64)
SFC_Z0M  = (SFC_Z0M_tmp*np.ones((JA,IA))).astype(np.float64)
SFC_Z0H  = (SFC_Z0H_tmp*np.ones((JA,IA))).astype(np.float64)
SFC_Z0E  = (SFC_Z0E_tmp*np.ones((JA,IA))).astype(np.float64)
SFC_DENS = (SFC_DENS_tmp*np.ones((TA,JA,IA))).astype(np.float64)
#
SFLX_MW  = np.empty((TA,JA,IA),dtype=np.float64)
SFLX_MU  = np.empty((TA,JA,IA),dtype=np.float64)
SFLX_MV  = np.empty((TA,JA,IA),dtype=np.float64)

#
IA = c_int32(IA)
IS = c_int32(IS)
IE = c_int32(IE)
JA = c_int32(JA)
JS = c_int32(JS)
JE = c_int32(JE)
TA = c_int32(TA)
TS = c_int32(TS)
TE = c_int32(TE)
f = np.ctypeslib.load_library("sample.so",".")
f.sample_.argtypes = [
    POINTER(c_int32),
    POINTER(c_int32),
    POINTER(c_int32),
    POINTER(c_int32),
    POINTER(c_int32),
    POINTER(c_int32),
    POINTER(c_int32),
    POINTER(c_int32),
    POINTER(c_int32),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64)]
f.sample_.restype = c_void_p
f.sample_(
    byref(IA),
    byref(IS),
    byref(IE),
    byref(JA),
    byref(JS),
    byref(JE),
    byref(TA),
    byref(TS),
    byref(TE),
    DX,
    DY,
    ATM_TEMP,
    SFC_TEMP,
    ATM_PRES,
    SFC_PRES,
    ATM_QV  ,
    ATM_U   ,
    ATM_V   ,
    ATM_Z1  ,
    PBL     ,
    SFC_Z0M ,
    SFC_Z0H ,
    SFC_Z0E ,
    SFC_DENS,
    SFLX_MW,
    SFLX_MU,
    SFLX_MV )
del f

#
print("SFLX_MU.shape",SFLX_MU.shape)
SFLX_MU = SFLX_MU.astype(np.float32)
SFLX_MV = SFLX_MV.astype(np.float32)
for t1 in range(TA.value):
    print("t=",t1)
    y = np.arange(JA.value+1)*DY
    x = np.arange(IA.value+1)*DX
    fig,axes = plt.subplots(1,4,figsize=(12,3))

    axes[0].pcolor(x*1e-3,y*1e-3, ATM_U[t1,:,:],cmap="bwr",vmin=-5,vmax=5)
    axes[1].pcolor(x*1e-3,y*1e-3, ATM_V[t1,:,:],cmap="bwr",vmin=-5,vmax=5)
    axes[2].pcolor(x*1e-3,y*1e-3, SFLX_MU[t1,:,:],cmap="bwr",vmin=-0.005,vmax=0.005)
    axes[3].pcolor(x*1e-3,y*1e-3, SFLX_MV[t1,:,:],cmap="bwr",vmin=-0.005,vmax=0.005)

    axes[0].set_title("ATM_U")
    axes[1].set_title("ATM_V")
    axes[2].set_title("SFLX_MU")
    axes[3].set_title("SFLX_MV")

    axes[0].set_aspect("equal")
    axes[0].set_xlabel("x (km)")
    axes[0].set_ylabel("y (km)")

    for ind in range(1,4):
        axes[ind].set_aspect("equal")
        axes[ind].tick_params(labelbottom=False,labelleft=False,labelright=False,labeltop=False)

    fig.tight_layout()
    fig.savefig("SFLX_MOM_t{:02d}.png".format(t1))
    plt.close("all")
