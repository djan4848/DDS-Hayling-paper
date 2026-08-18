"""FASE 2 — validez incremental. Criterios en PREREG.md (sellado antes de correr)."""
import numpy as np, pandas as pd, pickle, warnings; warnings.filterwarnings('ignore')
from scipy import stats
from scipy.signal import welch
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests
import sys; S='/tmp/claude-1000/-media-neuraldyn-PortableSSD-DDS-HAYLING/cc25759a-5311-4073-b69e-85fd4d2079cd/scratchpad/'
sys.path.insert(0,S)
from core import calculate_ais_shannon
FS=250.0; NPERM=10000; NBOOT=10000
rng=np.random.default_rng(2026)
K=pd.read_csv(S+'fase2_keys.csv'); D=pickle.load(open(S+'fase2_series.pkl','rb'))
RES,ERP=D['RES'],D['ERP']

# ---------- 2c: AIS sobre residuo corregido y sobre EEG CRUDO -----------------
K['ais_resid']=Parallel(n_jobs=18,batch_size=200)(delayed(calculate_ais_shannon)(y,8,1) for y in RES)
K['ais_raw']  =Parallel(n_jobs=18,batch_size=200)(delayed(calculate_ais_shannon)(y,8,1) for y in ERP)
print('[ais] done',flush=True)

# ---------- pendiente 1/f por unidad (OLS log-log) y PSD para specparam -------
def psd_of(y):
    f,p=welch(y,fs=FS,nperseg=len(y)); return f,p
f0,_=psd_of(RES[0]); band=(f0>=2)&(f0<=40)
P=np.array(Parallel(n_jobs=18,batch_size=200)(delayed(lambda y: psd_of(y)[1])(y) for y in RES))
K['chi_ols']=-np.polyfit(np.log10(f0[band]),np.log10(P[:,band].T+1e-300),1)[0]
K['resid_var']=[float(np.var(y)) for y in RES]
print('[psd] done',flush=True)

# specparam sobre PSD promediada por sujeto x cond x roi
from specparam import SpectralModel
grp=K.groupby(['subject','cond','roi']).indices
rows=[]
for k,ix in grp.items():
    psd=P[ix].mean(0)
    fm=SpectralModel(peak_width_limits=[2,12],max_n_peaks=3,min_peak_height=0.1,
                     aperiodic_mode='fixed',verbose=False)
    try:
        fm.fit(f0,psd,[2.0,40.0])
        chi=float(np.atleast_1d(fm.get_params('aperiodic'))[-1])
        r2=float(getattr(fm,'r_squared_',np.nan))
    except Exception:
        chi,r2=np.nan,np.nan
    rows.append(dict(subject=k[0],cond=k[1],roi=k[2],chi_sp=chi,chi_sp_r2=r2))
CHI=pd.DataFrame(rows); CHI.to_csv(S+'fase2_chi.csv',index=False)
print('[specparam] done  median R2=%.3f'%CHI.chi_sp_r2.median(),flush=True)
K.to_csv(S+'fase2_unit.csv',index=False)

# ---------- utilidades estadisticas ------------------------------------------
def perm_p(d,n=NPERM):
    obs=abs(d.mean()); sg=rng.choice([-1.,1.],size=(n,len(d)))
    return float(((np.abs((sg*d).mean(1))>=obs-1e-15).sum()+1)/(n+1))
def dz_ci(d,n=NBOOT):
    dz=d.mean()/d.std(ddof=1)
    b=[(lambda s:s.mean()/s.std(ddof=1))(rng.choice(d,len(d))) for _ in range(n)]
    return dz,*np.percentile(b,[2.5,97.5])
def sub(df,col,filt=None):
    x=df if filt is None else df[filt]
    return x.groupby(['subject','cond'])[col].mean().unstack().dropna(subset=['INIT','INHIB'])
def contrast(df,col,filt=None):
    s=sub(df,col,filt); d=(s.INIT-s.INHIB).values
    dz,lo,hi=dz_ci(d)
    return dict(n=len(s),mean_init=s.INIT.mean(),mean_inhib=s.INHIB.mean(),
                dz=dz,dz_lo=lo,dz_hi=hi,p_perm=perm_p(d)),s

# ---------- FAMILIA P: AIS x4 + TE x2 (residuo corregido) --------------------
TE=pd.read_csv('/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/'
               'te_n450_residual_r2pos_FIXED/te_n450_residual_r2pos_FIXED_results.csv')
rowsP=[]
for roi in ['frontal','cacc','lh','rh']:
    r,_=contrast(K,'ais_resid',K.roi==roi); rowsP.append(dict(family='P',contrast=f'AIS {roi}',**r))
for src,dst in [('frontal','cacc'),('cacc','frontal')]:
    t=TE[(TE.source_roi==src)&(TE.target_roi==dst)]
    r,_=contrast(t,'te_bits'); rowsP.append(dict(family='P',contrast=f'TE {src}->{dst}',**r))
P=pd.DataFrame(rowsP); P['q_perm']=multipletests(P.p_perm,method='fdr_bh')[1]
P['q_perm_AISonly']=np.nan
m=P.contrast.str.startswith('AIS'); P.loc[m,'q_perm_AISonly']=multipletests(P.loc[m,'p_perm'],method='fdr_bh')[1]

# ---------- FAMILIA C: AIS crudo x4 ------------------------------------------
rowsC=[]
for roi in ['frontal','cacc','lh','rh']:
    r,_=contrast(K,'ais_raw',K.roi==roi); rowsC.append(dict(family='C',contrast=f'AIS crudo {roi}',**r))
C=pd.DataFrame(rowsC); C['q_perm']=multipletests(C.p_perm,method='fdr_bh')[1]

# ---------- FAMILIA B: parametros DDS ----------------------------------------
dds=pd.read_csv('/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/'
                'dds_peak_aligned_n450/dds_n450_results.csv')
dds=dds[dds.r2>0].merge(K[['subject','cond','trial','roi']],on=['subject','cond','trial','roi'])
dds['absA1']=dds.A1.abs(); dds['absA2']=dds.A2.abs()
rowsB=[]
for par in ['absA1','gamma1','f1','absA2','gamma2','f2']:
    for roi in ['frontal','cacc','lh','rh']:
        r,_=contrast(dds,par,dds.roi==roi); rowsB.append(dict(family='B',contrast=f'{par} {roi}',**r))
B=pd.DataFrame(rowsB); B['q_perm']=multipletests(B.p_perm,method='fdr_bh')[1]
# fases: circular, exploratorio
rowsPh=[]
for par in ['phi1','phi2']:
    for roi in ['frontal','cacc','lh','rh']:
        d=dds[dds.roi==roi]
        cm=d.groupby(['subject','cond']).apply(lambda x: np.angle(np.mean(np.exp(1j*x[par]))))
        s=cm.unstack().dropna(subset=['INIT','INHIB'])
        diff=np.angle(np.exp(1j*(s.INIT-s.INHIB))).astype(float)
        dz,lo,hi=dz_ci(diff)
        rowsPh.append(dict(family='B-exploratorio',contrast=f'{par} {roi} (circular)',n=len(s),
                           mean_init=np.nan,mean_inhib=np.nan,dz=dz,dz_lo=lo,dz_hi=hi,p_perm=perm_p(diff)))
Ph=pd.DataFrame(rowsPh)
pd.concat([P,C,B,Ph],ignore_index=True).to_csv(S+'fase2_families.csv',index=False)
print('\n=== FAMILIA P ==='); print(P.round(4).to_string(index=False))
print('\n=== FAMILIA C ==='); print(C.round(4).to_string(index=False))
print('\n=== FAMILIA B ==='); print(B.round(4).to_string(index=False))
print('\n=== FASES (exploratorio) ==='); print(Ph.round(4).to_string(index=False))
print('FASE2_CORE_DONE',flush=True)
