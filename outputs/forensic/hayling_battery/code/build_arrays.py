import numpy as np, pandas as pd, pickle
B='/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/'
OUT='/tmp/claude-1000/-media-neuraldyn-PortableSSD-DDS-HAYLING/cc25759a-5311-4073-b69e-85fd4d2079cd/scratchpad/'
usec=['subject','cond','trial','roi','time_ms','value','dds_fit','residual','dds_r2']
res=pd.read_csv(B+'trial_roi_timeseries_residual_r2pos/trial_roi_timeseries_residual_r2pos.csv',usecols=usec)
dds=pd.read_csv(B+'dds_peak_aligned_n450/dds_n450_results.csv')
dds=dds[dds.r2>0].set_index(['subject','cond','trial','roi'])
dds=dds[~dds.index.duplicated()]
def model(t,p):
    return (p.A1*np.exp(-p.gamma1*t)*np.sin(2*np.pi*p.f1*t+p.phi1)
          + p.A2*np.exp(-p.gamma2*t)*np.sin(2*np.pi*p.f2*t+p.phi2))
keys=[];Ypub=[];Ycor=[];Yerp=[];R2=[];NW=[]
for kk,g in res.groupby(['subject','cond','trial','roi'],sort=False):
    if kk not in dds.index: continue
    p=dds.loc[kk]; g=g.sort_values('time_ms')
    t_ms=g.time_ms.values
    if len(t_ms)!=301: continue
    m=(t_ms>=p.window_tmin_ms-0.5)&(t_ms<=p.window_tmax_ms+0.5)
    if m.sum()!=101: NW.append(m.sum()); continue
    tw=(t_ms[m]-p.window_tmin_ms)/1000.0
    keys.append(kk); Ypub.append(g.residual.values.astype(np.float32))
    Ycor.append((g.value.values[m]-model(tw,p)).astype(np.float32))
    Yerp.append(g.value.values[m].astype(np.float32)); R2.append(float(p.r2))
K=pd.DataFrame(keys,columns=['subject','cond','trial','roi'])
np.savez_compressed(OUT+'arrays.npz',Ypub=np.array(Ypub),Ycor=np.array(Ycor),
                    Yerp=np.array(Yerp),R2=np.array(R2))
K.to_csv(OUT+'keys.csv',index=False)
print('units',len(K),'dropped_wrong_window',len(NW))
print(K.groupby(['cond','roi']).size())
print('subjects',K.subject.nunique())
