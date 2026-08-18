import numpy as np, pandas as pd
B='/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/'
res=pd.read_csv(B+'trial_roi_timeseries_residual_r2pos/trial_roi_timeseries_residual_r2pos.csv')
SUBS=['P1','P10','P11','P12','P13','P2','P20','P3','P30']
res=res[res.subject.isin(SUBS)]
dds=pd.read_csv(B+'dds_peak_aligned_n450/dds_n450_results.csv')
dds=dds[(dds.r2>0)&(dds.subject.isin(SUBS))].set_index(['subject','cond','trial','roi'])
dds=dds[~dds.index.duplicated()]
def model(t,p):
    return (p.A1*np.exp(-p.gamma1*t)*np.sin(2*np.pi*p.f1*t+p.phi1)
          + p.A2*np.exp(-p.gamma2*t)*np.sin(2*np.pi*p.f2*t+p.phi2))
def r2(y,yh): return 1-np.sum((y-yh)**2)/(np.sum((y-y.mean())**2)+1e-30)
rows=[]
for kk,g in res.groupby(['subject','cond','trial','roi'],sort=False):
    if kk not in dds.index: continue
    p=dds.loc[kk]
    g=g.sort_values('time_ms')
    t_ms=g.time_ms.values; y=g.value.values; fs=g.dds_fit.values
    m=(t_ms>=p.window_tmin_ms)&(t_ms<=p.window_tmax_ms)
    tw=(t_ms[m]-p.window_tmin_ms)/1000.0
    rows.append(dict(stored_r2=p.r2, r2_storedfit_win=r2(y[m],fs[m]),
                     r2_storedfit_full=r2(y,fs),
                     r2_correcttb_win=r2(y[m],model(tw,p)),
                     resid_var_ratio_full=np.var(g.residual.values)/(np.var(y)+1e-30)))
d=pd.DataFrame(rows)
print("n units:",len(d))
print(d.describe(percentiles=[.05,.25,.5,.75,.95]).T.to_string())
print("\ncorr(stored_r2, r2_correcttb_win) =", np.corrcoef(d.stored_r2,d.r2_correcttb_win)[0,1])
print("corr(stored_r2, r2_storedfit_win)  =", np.corrcoef(d.stored_r2,d.r2_storedfit_win)[0,1])
