import numpy as np, pandas as pd, time
from scipy import stats
from core import calculate_ais_shannon
d=np.load('arrays.npz'); K=pd.read_csv('keys.csv')
Ypub, Ycor, R2 = d['Ypub'].astype(float), d['Ycor'].astype(float), d['R2']
t0=time.time(); a=[calculate_ais_shannon(y) for y in Ypub[:500]]
print(f"AIS speed: {(time.time()-t0)/500*1000:.2f} ms/series (301 pts)")
K['ais_pub']=[calculate_ais_shannon(y) for y in Ypub]
K['ais_cor']=[calculate_ais_shannon(y) for y in Ycor]
K['r2']=R2
# sanity: reproduce published AIS
pub=pd.read_csv('/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/ais_n450_residual_r2pos/ais_n450_residual_r2pos_results.csv')
m=K.merge(pub[['subject','cond','trial','roi','ais_bits']],on=['subject','cond','trial','roi'])
print("reproduce published AIS: max|diff| =", np.nanmax(np.abs(m.ais_pub-m.ais_bits)), " n=",len(m))
def paired(col):
    out=[]
    for roi in ['frontal','cacc','lh','rh']:
        s=K[K.roi==roi].groupby(['subject','cond'])[col].mean().unstack().dropna()
        t,p=stats.ttest_rel(s.INIT,s.INHIB)
        diff=s.INIT-s.INHIB; dz=diff.mean()/diff.std(ddof=1)
        out.append(dict(roi=roi,n=len(s),mean_init=s.INIT.mean(),mean_inhib=s.INHIB.mean(),t=t,p=p,dz=dz))
    return pd.DataFrame(out)
print("\n--- AIS as-published residual (301 pts, -200..1000 ms) ---")
print(paired('ais_pub').to_string(index=False))
print("\n--- AIS CORRECTED residual (101 pts, fit window, correct time base) ---")
print(paired('ais_cor').to_string(index=False))
print("\n--- R2 confound: DDS r2 INIT vs INHIB, within subject ---")
print(paired('r2').to_string(index=False))
K.to_csv('unit_ais.csv',index=False)
