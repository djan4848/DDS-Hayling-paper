import numpy as np, pandas as pd, pickle, warnings; warnings.filterwarnings('ignore')
from scipy import stats
from scipy.signal import welch
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests
import sys; sys.path.insert(0,'.')
from core import calculate_ais_shannon, iaaft_batch, N_SURR, SEED
rng=np.random.default_rng(2026); NPERM=10000
K=pd.read_csv('fase2_unit.csv'); D=pickle.load(open('fase2_series.pkl','rb')); ERP=D['ERP']
def perm_p(d,n=NPERM):
    obs=abs(d.mean()); sg=rng.choice([-1.,1.],size=(n,len(d)))
    return float(((np.abs((sg*d).mean(1))>=obs-1e-15).sum()+1)/(n+1))
def con(df,col,f):
    s=df[f].groupby(['subject','cond'])[col].mean().unstack().dropna(subset=['INIT','INHIB'])
    d=(s.INIT-s.INHIB).values
    return dict(n=len(s),mean_init=s.INIT.mean(),mean_inhib=s.INHIB.mean(),
                dz=d.mean()/d.std(ddof=1),p_perm=perm_p(d))
# ---- IAAFT sobre el EEG CRUDO -----------------------------------------------
def dnl(y):
    a=calculate_ais_shannon(y,8,1)
    if not np.isfinite(a): return np.nan,np.nan
    v=[calculate_ais_shannon(q,8,1) for q in iaaft_batch(y,n_surr=N_SURR,seed=SEED)]
    v=[q for q in v if np.isfinite(q)]
    return (a-np.mean(v),np.mean(v)) if v else (np.nan,np.nan)
r=np.array(Parallel(n_jobs=18,batch_size=20)(delayed(dnl)(y) for y in ERP))
K['dnl_raw'],K['ais_surr_raw']=r[:,0],r[:,1]
rows=[]
for roi in ['frontal','cacc','lh','rh']:
    f=K.roi==roi
    rows.append(dict(roi=roi,quantity='AIS crudo (real)',**con(K,'ais_raw',f)))
    rows.append(dict(roi=roi,quantity='AIS crudo en surrogados IAAFT',**con(K,'ais_surr_raw',f)))
    rows.append(dict(roi=roi,quantity='dA_NL crudo (excedente no lineal)',**con(K,'dnl_raw',f)))
R=pd.DataFrame(rows); R.to_csv('closing_raw_iaaft.csv',index=False)
print(R.round(4).to_string(index=False))
print('\nfraccion no lineal del AIS crudo: %.4f / %.4f = %.1f%%'%(
      K.dnl_raw.mean(),K.ais_raw.mean(),100*K.dnl_raw.mean()/K.ais_raw.mean()))
K.to_csv('fase2_unit.csv',index=False)
# ---- diagnostico de colinealidad de la escalera ------------------------------
CHI=pd.read_csv('fase2_chi.csv'); F=K[K.roi=='frontal']
sub=F.groupby(['subject','cond'])[['ais_resid','ais_raw','amp_n450','chi_ols','resid_var']].mean().unstack()
chi=CHI[CHI.roi=='frontal'].set_index(['subject','cond']).chi_sp.unstack()
DF=pd.DataFrame({c:(sub[c]['INIT']-sub[c]['INHIB']) for c in ['ais_resid','ais_raw','amp_n450','chi_ols','resid_var']})
DF['chi_sp']=chi['INIT']-chi['INHIB']; DF=DF.dropna()
print('\n--- correlaciones entre diferencias INIT-INHIB (frontal, n=%d) ---'%len(DF))
print(DF.corr().round(3).to_string())
DF.to_csv('ladder_diffs.csv')
# ---- calidad del ajuste specparam -------------------------------------------
from specparam import SpectralModel
P=np.array(Parallel(n_jobs=18,batch_size=200)(delayed(lambda y: welch(y,fs=250,nperseg=len(y))[1])(y) for y in D['RES']))
f0=welch(D['RES'][0],fs=250,nperseg=len(D['RES'][0]))[0]
r2s=[]
for k,ix in K.groupby(['subject','cond','roi']).indices.items():
    fm=SpectralModel(peak_width_limits=[2,12],max_n_peaks=3,min_peak_height=0.1,aperiodic_mode='fixed',verbose=False)
    try:
        fm.fit(f0,P[ix].mean(0),[2.,40.]); r2s.append(float(fm.results.metrics.results['gof_rsquared']))
    except Exception: pass
print('\nspecparam GOF R2: n=%d  mediana %.3f  p05 %.3f'%(len(r2s),np.median(r2s),np.percentile(r2s,5)) if r2s else 'specparam R2 no disponible')
print('CLOSING_DONE')
