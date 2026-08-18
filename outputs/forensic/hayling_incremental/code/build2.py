"""FASE 2: extrae, del residuo CORREGIDO, las series por unidad trial x ROI:
   residual, ERP crudo (misma ventana peak-aligned), amplitud N450, R2."""
import numpy as np, pandas as pd, pickle
B='/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/'
F=B+'trial_roi_timeseries_residual_r2pos_FIXED/trial_roi_timeseries_residual_r2pos_FIXED.csv'
df=pd.read_csv(F,usecols=['subject','cond','trial','trial_uid','roi','sample_idx',
                          't_win_s','value','dds_fit','residual','dds_r2'])
keys=[];RES=[];ERP=[];AMP=[];PKA=[];R2=[]
for kk,g in df.groupby(['subject','cond','trial','roi'],sort=False):
    g=g.sort_values('sample_idx')
    t=g.t_win_s.values; y=g.value.values; e=g.residual.values
    trel=t-0.200                                  # 0 = pico N450
    m=np.abs(trel)<=0.0501                        # pico +/- 50 ms
    keys.append(kk); RES.append(e.astype(np.float32)); ERP.append(y.astype(np.float32))
    AMP.append(float(y[m].mean()))                # amplitud media N450 peak-aligned
    PKA.append(float(y[np.abs(trel).argmin()]))   # amplitud en el pico
    R2.append(float(g.dds_r2.iloc[0]))
K=pd.DataFrame(keys,columns=['subject','cond','trial','roi'])
K['amp_n450']=AMP; K['amp_peak']=PKA; K['r2']=R2; K['n']= [len(x) for x in RES]
with open('/tmp/claude-1000/-media-neuraldyn-PortableSSD-DDS-HAYLING/cc25759a-5311-4073-b69e-85fd4d2079cd/scratchpad/fase2_series.pkl','wb') as f:
    pickle.dump({'RES':RES,'ERP':ERP},f,protocol=4)
K.to_csv('/tmp/claude-1000/-media-neuraldyn-PortableSSD-DDS-HAYLING/cc25759a-5311-4073-b69e-85fd4d2079cd/scratchpad/fase2_keys.csv',index=False)
print('units',len(K),'len dist',K.n.value_counts().to_dict())
print(K.groupby(['cond','roi']).size().to_dict())
