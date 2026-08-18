"""FASE 2a (escalera de modelos) + extras (AUC, variabilidad trial-a-trial).
   Criterios pre-declarados en PREREG.md."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import statsmodels.formula.api as smf
S='/tmp/claude-1000/-media-neuraldyn-PortableSSD-DDS-HAYLING/cc25759a-5311-4073-b69e-85fd4d2079cd/scratchpad/'
NPERM=10000; NBOOT=10000
rng=np.random.default_rng(2026)
K=pd.read_csv(S+'fase2_unit.csv'); CHI=pd.read_csv(S+'fase2_chi.csv')
F=K[K.roi=='frontal'].copy()

# ---------- nivel sujeto: diferencias INIT-INHIB ------------------------------
sub=F.groupby(['subject','cond'])[['ais_resid','ais_raw','amp_n450','amp_peak','chi_ols','r2','resid_var']].mean()
sub=sub.unstack()
chi=CHI[CHI.roi=='frontal'].set_index(['subject','cond']).chi_sp.unstack()
DF=pd.DataFrame({c:(sub[c]['INIT']-sub[c]['INHIB']) for c in
                 ['ais_resid','ais_raw','amp_n450','amp_peak','chi_ols','r2','resid_var']})
DF['chi_sp']=chi['INIT']-chi['INHIB']
DF=DF.dropna()
print('n sujetos escalera:',len(DF))

def ancova(dv,covs):
    y=DF[dv].values
    X=np.column_stack([np.ones(len(DF))]+[DF[c].values-DF[c].mean() for c in covs])
    beta,_,_,_=np.linalg.lstsq(X,y,rcond=None)
    res=y-X@beta; s=res.std(ddof=X.shape[1])
    dz=beta[0]/s
    obs=abs(beta[0]); cnt=0
    for _ in range(NPERM):
        sg=rng.choice([-1.,1.],len(y))
        Xp=np.column_stack([np.ones(len(y))]+[sg*(DF[c].values-DF[c].mean()) for c in covs])
        b,_,_,_=np.linalg.lstsq(Xp,sg*y,rcond=None); cnt+=abs(b[0])>=obs-1e-15
    # dz bootstrap
    bs=[]
    for _ in range(2000):
        i=rng.integers(0,len(y),len(y)); yy=y[i]
        XX=np.column_stack([np.ones(len(i))]+[DF[c].values[i]-DF[c].values[i].mean() for c in covs])
        try:
            bb,_,_,_=np.linalg.lstsq(XX,yy,rcond=None); rr=yy-XX@bb
            bs.append(bb[0]/(rr.std(ddof=XX.shape[1])+1e-30))
        except Exception: pass
    lo,hi=np.percentile(bs,[2.5,97.5])
    return dict(dz=dz,dz_lo=lo,dz_hi=hi,p_perm=float((cnt+1)/(NPERM+1)),
                beta_cov={c:round(b,4) for c,b in zip(covs,beta[1:])})

MODELS=[('M0  AIS_resid ~ 1','ais_resid',[]),
        ('M-i  + amplitud N450','ais_resid',['amp_n450']),
        ('M-ii + exponente 1/f (specparam)','ais_resid',['chi_sp']),
        ('M-ii\' + exponente 1/f (OLS trial)','ais_resid',['chi_ols']),
        ('M-iii + AIS crudo','ais_resid',['ais_raw']),
        ('M-full + amp + 1/f + AIS crudo','ais_resid',['amp_n450','chi_sp','ais_raw']),
        ('M-R2 + R2 del ajuste','ais_resid',['r2'])]
rows=[]
for lab,dv,cov in MODELS:
    r=ancova(dv,cov); rows.append(dict(model=lab,dv=dv,covariables=', '.join(cov) or '—',**r))
# reciprocos: cada candidato como DV, solo y controlado por AIS_resid
for dv,lab in [('amp_n450','amplitud N450'),('chi_sp','exponente 1/f'),('ais_raw','AIS crudo'),
               ('amp_peak','amplitud en el pico'),('resid_var','varianza del residuo')]:
    r=ancova(dv,[]);            rows.append(dict(model=f'REC {lab} ~ 1',dv=dv,covariables='—',**r))
    r=ancova(dv,['ais_resid']); rows.append(dict(model=f'REC {lab} ~ 1 + AIS_resid',dv=dv,covariables='ais_resid',**r))
L=pd.DataFrame(rows); L.to_csv(S+'ladder.csv',index=False)
print(L[['model','dz','dz_lo','dz_hi','p_perm','beta_cov']].round(4).to_string(index=False))

# ---------- robustez: modelo mixto a nivel de trial ---------------------------
T=F.dropna(subset=['ais_resid','ais_raw','amp_n450','chi_ols']).copy()
T['cond_i']=(T.cond=='INIT').astype(int)
for name,form in [('cond solo','ais_resid ~ cond_i'),
                  ('+ amplitud','ais_resid ~ cond_i + amp_n450'),
                  ('+ 1/f (OLS)','ais_resid ~ cond_i + chi_ols'),
                  ('+ AIS crudo','ais_resid ~ cond_i + ais_raw'),
                  ('+ todo','ais_resid ~ cond_i + amp_n450 + chi_ols + ais_raw')]:
    try:
        m=smf.mixedlm(form,T,groups=T.subject).fit(method='lbfgs')
        print(f"[LMM trial] {name:14s} beta_cond={m.params['cond_i']:+.5f}  p={m.pvalues['cond_i']:.4g}")
    except Exception as e: print('[LMM]',name,'fail',e)

# ---------- EXTRA: AUC --------------------------------------------------------
Z=F.dropna(subset=['ais_resid']).copy()
Z['z']=Z.groupby('subject').ais_resid.transform(lambda x:(x-x.mean())/(x.std(ddof=1)+1e-12))
y=(Z.cond=='INIT').astype(int).values; z=Z.z.values
auc_pool=roc_auc_score(y,z)
subs=Z.subject.values; uniq=np.unique(subs)
null=[]
for _ in range(5000):
    yp=y.copy()
    for s_ in uniq:
        m=subs==s_; yp[m]=rng.permutation(y[m])
    null.append(roc_auc_score(yp,z))
null=np.array(null)
p_auc=float(((np.abs(null-0.5)>=abs(auc_pool-0.5)).sum()+1)/(len(null)+1))
# LOSO
pred=np.full(len(Z),np.nan)
for s_ in uniq:
    te=subs==s_; tr=~te
    lr=LogisticRegression().fit(Z.ais_resid.values[tr].reshape(-1,1),y[tr])
    pred[te]=lr.predict_proba(Z.ais_resid.values[te].reshape(-1,1))[:,1]
auc_loso=roc_auc_score(y,pred)
bs=[]
for _ in range(NBOOT//5):
    ss=rng.choice(uniq,len(uniq)); m=np.concatenate([np.where(subs==s_)[0] for s_ in ss])
    try: bs.append(roc_auc_score(y[m],pred[m]))
    except Exception: pass
lo,hi=np.percentile(bs,[2.5,97.5])
# nivel sujeto: signo de la diferencia
acc=float((DF.ais_resid>0).mean()); pb=stats.binomtest(int((DF.ais_resid>0).sum()),len(DF),0.5).pvalue
print(f"\n[AUC] pooled within-subject z : {auc_pool:.4f}  p_perm={p_auc:.4f}  (nulo {null.mean():.4f})")
print(f"[AUC] LOSO                    : {auc_loso:.4f}  IC95 [{lo:.4f}, {hi:.4f}]")
print(f"[SUJETO] fraccion con INIT>INHIB: {acc:.3f}  (binomial p={pb:.4f})")
pd.DataFrame([dict(auc_pooled=auc_pool,p_perm=p_auc,auc_loso=auc_loso,loso_lo=lo,loso_hi=hi,
                   subj_sign_acc=acc,subj_binom_p=pb)]).to_csv(S+'auc.csv',index=False)

# ---------- EXTRA: variabilidad trial-a-trial ---------------------------------
dds=pd.read_csv('/media/neuraldyn/PortableSSD/DDS-HAYLING/DDS-Hayling/derivatives/'
                'dds_peak_aligned_n450/dds_n450_results.csv')
dds=dds[dds.r2>0]; dds['absA1']=dds.A1.abs(); dds['absA2']=dds.A2.abs()
M=K.merge(dds[['subject','cond','trial','roi','absA1','gamma1','f1','absA2','gamma2','f2']],
          on=['subject','cond','trial','roi'],how='left')
VARS=['ais_resid','resid_var','absA1','gamma1','f1','absA2','gamma2','f2']
rows=[]
for roi in ['frontal','cacc']:
    d=M[M.roi==roi]
    out={v:{} for v in VARS}
    for s_,g in d.groupby('subject'):
        gi=g[g.cond=='INIT']; gh=g[g.cond=='INHIB']
        n=min(len(gi),len(gh))
        if n<10: continue
        acc={v:[[],[]] for v in VARS}
        for _ in range(20):
            a=gi.sample(n,random_state=int(rng.integers(1e6)))
            b=gh.sample(n,random_state=int(rng.integers(1e6)))
            for v in VARS:
                acc[v][0].append(a[v].std(ddof=1)); acc[v][1].append(b[v].std(ddof=1))
        for v in VARS: out[v][s_]=(np.mean(acc[v][0]),np.mean(acc[v][1]))
    for v in VARS:
        vals=np.array([x for x in out[v].values()]); vals=vals[np.isfinite(vals).all(1)]
        if len(vals)<10: continue
        d_=vals[:,0]-vals[:,1]; dz=d_.mean()/d_.std(ddof=1)
        obs=abs(d_.mean()); sg=rng.choice([-1.,1.],size=(NPERM,len(d_)))
        p=float(((np.abs((sg*d_).mean(1))>=obs-1e-15).sum()+1)/(NPERM+1))
        rows.append(dict(roi=roi,measure=f'SD entre trials de {v}',n=len(vals),dz=dz,p_perm=p))
V=pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
V['q_bh_exploratorio']=multipletests(V.p_perm,method='fdr_bh')[1]
V.to_csv(S+'trialvar.csv',index=False)
print('\n=== VARIABILIDAD TRIAL-A-TRIAL (exploratorio) ===')
print(V.round(4).to_string(index=False))
print('LADDER_DONE')
