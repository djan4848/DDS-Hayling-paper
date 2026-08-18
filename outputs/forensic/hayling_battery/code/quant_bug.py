import numpy as np, pandas as pd
B='/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/'
K=pd.read_csv('keys.csv'); d=np.load('arrays.npz')
Yerp=d['Yerp'].astype(float); Ycor=d['Ycor'].astype(float); R2=d['R2']
dds=pd.read_csv(B+'dds_peak_aligned_n450/dds_n450_results.csv')
dds=dds[dds.r2>0].set_index(['subject','cond','trial','roi']); dds=dds[~dds.index.duplicated()]
P=dds.reindex(pd.MultiIndex.from_frame(K[['subject','cond','trial','roi']]))
tmin=P.window_tmin_ms.values/1000.0
off=np.arange(101)*0.004
Tabs=tmin[:,None]+off[None,:]          # absolute epoch time used by the residual builder
Twin=off[None,:]*np.ones((len(K),1))   # window-relative time used by the fit
def M(T):
    return (P.A1.values[:,None]*np.exp(-P.gamma1.values[:,None]*T)*np.sin(2*np.pi*P.f1.values[:,None]*T+P.phi1.values[:,None])
          + P.A2.values[:,None]*np.exp(-P.gamma2.values[:,None]*T)*np.sin(2*np.pi*P.f2.values[:,None]*T+P.phi2.values[:,None]))
fit_wrong=M(Tabs); fit_right=M(Twin)
def r2v(y,yh):
    return 1-((y-yh)**2).sum(1)/(((y-y.mean(1,keepdims=True))**2).sum(1)+1e-300)
rw=r2v(Yerp,fit_wrong); rr=r2v(Yerp,fit_right)
resid_wrong=Yerp-fit_wrong
vr=resid_wrong.var(1)/(Yerp.var(1)+1e-300)
print("n units:",len(K))
print("R2 of the *as-published* fit inside its own fit window:")
print(pd.Series(rw).describe(percentiles=[.05,.25,.5,.75,.95]).to_string())
print("\nR2 of the fit re-evaluated on the correct time base (should equal stored r2):")
print(pd.Series(rr).describe(percentiles=[.05,.5,.95]).to_string())
print("max |rr - stored_r2| =", np.nanmax(np.abs(rr-R2)))
print("\nvar(published residual)/var(ERP), in-window:")
print(pd.Series(vr).describe(percentiles=[.05,.5,.95]).to_string())
print("\n%% units where published fit has R2<0 in-window: %.1f%%"%(100*(rw<0).mean()))
print("%% units where var(resid)>var(ERP): %.1f%%"%(100*(vr>1).mean()))
# how much of AIS is linear
K2=pd.read_csv('unit_battery.csv')
for v in ['pub','cor']:
    print(f"\n[{v}] mean AIS={K2['ais_'+v].mean():.4f}  mean AIS_surr={K2['ais_surr_'+v].mean():.4f}  "
          f"dA_NL={K2['dnl_'+v].mean():.4f}  ({100*K2['dnl_'+v].mean()/K2['ais_'+v].mean():.1f}% of AIS)")
