"""Barridos de sensibilidad sobre el residuo CORREGIDO (estimadores congelados)."""
import numpy as np, pandas as pd, pickle, warnings; warnings.filterwarnings('ignore')
from scipy import stats
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests
import sys; sys.path.insert(0,'/tmp/claude-1000/-media-neuraldyn-PortableSSD-DDS-HAYLING/cc25759a-5311-4073-b69e-85fd4d2079cd/scratchpad')
from core import calculate_ais_shannon, transfer_entropy_discrete
S='/tmp/claude-1000/-media-neuraldyn-PortableSSD-DDS-HAYLING/cc25759a-5311-4073-b69e-85fd4d2079cd/scratchpad/'
K=pd.read_csv(S+'fase2_keys.csv'); D=pickle.load(open(S+'fase2_series.pkl','rb')); RES=D['RES']
rng=np.random.default_rng(2026)
def perm_p(diff,n=10000):
    obs=abs(diff.mean()); sg=rng.choice([-1.,1.],size=(n,len(diff)))
    return float(((np.abs((sg*diff).mean(1))>=obs-1e-15).sum()+1)/(n+1))
def paired(K,col,roi):
    s=K[K.roi==roi].groupby(['subject','cond'])[col].mean().unstack().dropna(subset=['INIT','INHIB'])
    d=(s.INIT-s.INHIB).values; t,p=stats.ttest_rel(s.INIT,s.INHIB)
    return dict(n=len(s),mean_init=s.INIT.mean(),mean_inhib=s.INHIB.mean(),
                t=t,p_param=p,dz=d.mean()/d.std(ddof=1),p_perm=perm_p(d))
# ---- A) AIS: barrido de bins (lag=1 congelado) ------------------------------
rows=[]
for nb in [4,6,8,10]:
    v=Parallel(n_jobs=18,batch_size=200)(delayed(calculate_ais_shannon)(y,nb,1) for y in RES)
    K[f'ais_b{nb}']=v
    for roi in ['frontal','cacc','lh','rh']:
        rows.append(dict(measure='AIS',param='bins',value=nb,lag=1,roi=roi,**paired(K,f'ais_b{nb}',roi)))
    print('bins',nb,'done',flush=True)
# ---- B) AIS: barrido de ventana (pico +/-200 vs +/-100 vs +/-50 ms) ---------
for half,lab in [(50,'+/-200ms (ajuste)'),(25,'+/-100ms'),(12,'+/-50ms')]:
    v=[]
    for y in RES:
        c=len(y)//2
        v.append(calculate_ais_shannon(y[max(0,c-half):c+half+1],8,1))
    K[f'ais_w{half}']=v
    for roi in ['frontal','cacc','lh','rh']:
        rows.append(dict(measure='AIS',param='window',value=lab,lag=1,roi=roi,**paired(K,f'ais_w{half}',roi)))
    print('window',lab,'done',flush=True)
A=pd.DataFrame(rows)
for (m,p_,v) in A[['measure','param','value']].drop_duplicates().itertuples(index=False):
    m_=(A.measure==m)&(A.param==p_)&(A.value==v)
    A.loc[m_,'q_perm']=multipletests(A.loc[m_,'p_perm'],method='fdr_bh')[1]
A.to_csv(S+'sweep_ais.csv',index=False)
print(A.to_string(index=False),flush=True)
# ---- C) TE: barrido de lag x bins -------------------------------------------
idx={(s,c,t,r):i for i,(s,c,t,r) in enumerate(K[['subject','cond','trial','roi']].itertuples(index=False))}
tr=K[K.roi=='frontal'][['subject','cond','trial']].drop_duplicates()
pairs=[(idx[(s,c,t,'frontal')],idx[(s,c,t,'cacc')],s,c,t) for s,c,t in tr.itertuples(index=False) if (s,c,t,'cacc') in idx]
def te_pair(i,j,lag,nb):
    return (transfer_entropy_discrete(RES[i],RES[j],lag,nb), transfer_entropy_discrete(RES[j],RES[i],lag,nb))
rows=[]
for lag in [1,2,4]:
    for nb in [4,8]:
        r=np.array(Parallel(n_jobs=18,batch_size=200)(delayed(te_pair)(i,j,lag,nb) for i,j,_,_,_ in pairs))
        T=pd.DataFrame([(s,c,t) for _,_,s,c,t in pairs],columns=['subject','cond','trial'])
        T['fc'],T['cf']=r[:,0],r[:,1]; T['roi']='x'
        for d,lab in [('fc','frontal->cacc'),('cf','cacc->frontal')]:
            rows.append(dict(measure='TE',direction=lab,lag=lag,bins=nb,**paired(T,d,'x')))
        print('TE lag',lag,'bins',nb,'done',flush=True)
B=pd.DataFrame(rows)
for (l,b) in B[['lag','bins']].drop_duplicates().itertuples(index=False):
    m_=(B.lag==l)&(B.bins==b); B.loc[m_,'q_perm']=multipletests(B.loc[m_,'p_perm'],method='fdr_bh')[1]
B.to_csv(S+'sweep_te.csv',index=False)
print(B.to_string(index=False))
K.to_csv(S+'fase2_keys_sweep.csv',index=False)
print('SWEEPS_DONE')
