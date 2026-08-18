import numpy as np, pandas as pd, sys, time, warnings
warnings.filterwarnings('ignore')
from scipy import stats
from scipy.signal import welch
from joblib import Parallel, delayed
from core import (calculate_ais_shannon, transfer_entropy_discrete, compute_epr,
                  iaaft_batch, iaaft_batch_bivariate, N_SURR, SEED)

NJOBS=18; NPERM=10000
d=np.load('arrays.npz'); K=pd.read_csv('keys.csv')
Y={'pub':d['Ypub'].astype(float),'cor':d['Ycor'].astype(float)}
K['r2']=d['R2']
FS=250.0

# ---------- per-unit scalars ----------
def unit_scalars(y):
    x=y-y.mean(); sd=x.std()
    r1=float(np.dot(x[:-1],x[1:])/(len(x)*sd**2)) if sd>0 else np.nan
    f,p=welch(y,fs=FS,nperseg=len(y))
    m=(f>=1)&(f<=40)
    slope=float(np.polyfit(np.log10(f[m]),np.log10(p[m]+1e-300),1)[0]) if m.sum()>4 else np.nan
    return calculate_ais_shannon(y), compute_epr(y), r1, -slope

for v in ['pub','cor']:
    r=Parallel(n_jobs=NJOBS,batch_size=200)(delayed(unit_scalars)(y) for y in Y[v])
    r=np.array(r)
    K[f'ais_{v}'],K[f'epr_{v}'],K[f'r1_{v}'],K[f'chi_{v}']=r[:,0],r[:,1],r[:,2],r[:,3]
    print(f'[scalars {v}] done',flush=True)

# ---------- IAAFT delta A_NL ----------
def dnl(y):
    a=calculate_ais_shannon(y)
    if not np.isfinite(a): return np.nan,np.nan
    s=iaaft_batch(y,n_surr=N_SURR,seed=SEED)
    vals=[calculate_ais_shannon(q) for q in s]
    vals=[q for q in vals if np.isfinite(q)]
    return (a-float(np.mean(vals))) if vals else np.nan, float(np.mean(vals)) if vals else np.nan
for v in ['pub','cor']:
    t0=time.time()
    r=np.array(Parallel(n_jobs=NJOBS,batch_size=20)(delayed(dnl)(y) for y in Y[v]))
    K[f'dnl_{v}'],K[f'ais_surr_{v}']=r[:,0],r[:,1]
    print(f'[iaaft {v}] done in {time.time()-t0:.0f}s',flush=True)
K.to_csv('unit_battery.csv',index=False)

# ---------- TE (frontal<->cacc) real + bivariate surrogate null ----------
idx={(s,c,t,ro):i for i,(s,c,t,ro) in enumerate(K[['subject','cond','trial','roi']].itertuples(index=False))}
trials=K[K.roi=='frontal'][['subject','cond','trial']].drop_duplicates()
pairs=[(idx[(s,c,t,'frontal')],idx[(s,c,t,'cacc')],s,c,t)
       for s,c,t in trials.itertuples(index=False) if (s,c,t,'cacc') in idx]
print('TE pairs:',len(pairs),flush=True)

def te_unit(i,j,v):
    x,y=Y[v][i],Y[v][j]
    te_fc=transfer_entropy_discrete(x,y); te_cf=transfer_entropy_discrete(y,x)
    SX,SY=iaaft_batch_bivariate(x,y,n_surr=N_SURR,seed=SEED)
    sfc=[transfer_entropy_discrete(a,b) for a,b in zip(SX,SY)]
    scf=[transfer_entropy_discrete(b,a) for a,b in zip(SX,SY)]
    sfc=[q for q in sfc if np.isfinite(q)]; scf=[q for q in scf if np.isfinite(q)]
    return (te_fc,te_cf,np.mean(sfc) if sfc else np.nan,np.mean(scf) if scf else np.nan)
TE={}
for v in ['pub','cor']:
    t0=time.time()
    r=np.array(Parallel(n_jobs=NJOBS,batch_size=20)(delayed(te_unit)(i,j,v) for i,j,_,_,_ in pairs))
    TE[v]=r; print(f'[te {v}] done in {time.time()-t0:.0f}s',flush=True)
TEdf=pd.DataFrame([(s,c,t) for _,_,s,c,t in pairs],columns=['subject','cond','trial'])
for v in ['pub','cor']:
    TEdf[f'te_fc_{v}'],TEdf[f'te_cf_{v}'],TEdf[f'tes_fc_{v}'],TEdf[f'tes_cf_{v}']=TE[v].T
    TEdf[f'dte_fc_{v}']=TEdf[f'te_fc_{v}']-TEdf[f'tes_fc_{v}']
    TEdf[f'dte_cf_{v}']=TEdf[f'te_cf_{v}']-TEdf[f'tes_cf_{v}']
TEdf.to_csv('unit_te.csv',index=False)
print('ALL DONE',flush=True)
