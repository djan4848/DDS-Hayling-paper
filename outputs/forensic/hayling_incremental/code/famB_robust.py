import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from scipy import stats
from statsmodels.stats.multitest import multipletests
rng=np.random.default_rng(2026); NPERM=10000
K=pd.read_csv('fase2_keys.csv')
dds=pd.read_csv('/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/'
                'dds_peak_aligned_n450/dds_n450_results.csv')
dds=dds[dds.r2>0].merge(K[['subject','cond','trial','roi']],on=['subject','cond','trial','roi'])
dds['absA1']=dds.A1.abs(); dds['absA2']=dds.A2.abs()
dds['logA1']=np.log10(dds.absA1+1e-12); dds['logA2']=np.log10(dds.absA2+1e-12)
def perm_p(d,n=NPERM):
    obs=abs(d.mean()); sg=rng.choice([-1.,1.],size=(n,len(d)))
    return float(((np.abs((sg*d).mean(1))>=obs-1e-15).sum()+1)/(n+1))
def con(df,col,roi,agg):
    s=df[df.roi==roi].groupby(['subject','cond'])[col].agg(agg).unstack().dropna(subset=['INIT','INHIB'])
    d=(s.INIT-s.INHIB).values
    return dict(n=len(s),dz=d.mean()/d.std(ddof=1),p_perm=perm_p(d))
tm=lambda x: stats.trim_mean(x,0.1)
rows=[]
for agg,lab in [('mean','media (pre-declarado)'),('median','mediana'),(tm,'media recortada 10%')]:
    for par in ['absA1','gamma1','f1','absA2','gamma2','f2']:
        for roi in ['frontal','cacc','lh','rh']:
            rows.append(dict(aggreg=lab,param=par,roi=roi,**con(dds,par,roi,agg)))
for par in ['logA1','logA2']:
    for roi in ['frontal','cacc','lh','rh']:
        rows.append(dict(aggreg='media de log10|A|',param=par,roi=roi,**con(dds,par,roi,'mean')))
R=pd.DataFrame(rows)
for a in R.aggreg.unique():
    m=R.aggreg==a; R.loc[m,'q_bh']=multipletests(R.loc[m,'p_perm'],method='fdr_bh')[1]
R.to_csv('famB_robust.csv',index=False)
print(R[R.q_bh<0.2].round(4).to_string(index=False))
print('\n--- solo frontal, todas las agregaciones ---')
print(R[R.roi=='frontal'].round(4).to_string(index=False))
print('\n--- fraccion de trials con |A1|>1e-4 por condicion (frontal) ---')
f=dds[dds.roi=='frontal']
print(f.assign(big=f.absA1>1e-4).groupby(['subject','cond']).big.mean().unstack().describe().to_string())
s=f.assign(big=(f.absA1>1e-4).astype(float)).groupby(['subject','cond']).big.mean().unstack().dropna()
d=(s.INIT-s.INHIB).values
print('contraste fraccion-outliers frontal: dz=%.3f p_perm=%.4f'%(d.mean()/d.std(ddof=1),perm_p(d)))
