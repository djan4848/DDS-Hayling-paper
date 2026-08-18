import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from scipy import stats
import diptest
from statsmodels.stats.multitest import multipletests
rng=np.random.default_rng(2026); NPERM=10000; NBOOT=10000

K=pd.read_csv('unit_battery.csv'); TE=pd.read_csv('unit_te.csv')

def subj(df,val,filt=None):
    x=df if filt is None else df[filt]
    return x.groupby(['subject','cond'])[val].mean().unstack().dropna(subset=['INIT','INHIB'])

def dz_ci(diff,nboot=NBOOT):
    dz=diff.mean()/diff.std(ddof=1)
    b=[(lambda s: s.mean()/s.std(ddof=1))(rng.choice(diff,len(diff))) for _ in range(nboot)]
    lo,hi=np.percentile(b,[2.5,97.5]); return dz,lo,hi

def perm_p(diff,nperm=NPERM):
    obs=abs(diff.mean())
    sg=rng.choice([-1.,1.],size=(nperm,len(diff)))
    return float(((np.abs((sg*diff).mean(1))>=obs-1e-15).sum()+1)/(nperm+1))

def mean_ci(x,nboot=NBOOT):
    b=[rng.choice(x,len(x)).mean() for _ in range(nboot)]
    return float(np.mean(x)),*np.percentile(b,[2.5,97.5])

def r2_control(dv,dr2):
    """subject-level ANCOVA on within-subject differences: dv_diff ~ 1 + r2_diff.
       intercept = condition effect at zero fit-quality difference."""
    X=np.column_stack([np.ones(len(dr2)),dr2-dr2.mean()])
    beta,_,_,_=np.linalg.lstsq(X,dv,rcond=None)
    res=dv-X@beta; s=res.std(ddof=2)
    dz_adj=beta[0]/s
    obs=abs(beta[0]); cnt=0
    for _ in range(NPERM):
        sg=rng.choice([-1.,1.],len(dv))
        b,_,_,_=np.linalg.lstsq(np.column_stack([np.ones(len(dr2)),sg*(dr2-dr2.mean())]),sg*dv,rcond=None)
        cnt+= abs(b[0])>=obs-1e-15
    return dz_adj,float((cnt+1)/(NPERM+1))

rows=[]
CONTRASTS=[('AIS frontal','frontal','ais'),('AIS cACC','cacc','ais'),
           ('AIS LH','lh','ais'),('AIS RH','rh','ais')]
for v in ['pub','cor']:
    for label,roi,_ in CONTRASTS:
        f=K.roi==roi
        s=subj(K,f'ais_{v}',f); dv=(s.INIT-s.INHIB).values
        sr=subj(K,'r2',f); dr2=(sr.INIT-sr.INHIB).reindex(s.index).values
        dz,lo,hi=dz_ci(dv); p=perm_p(dv)
        dza,pa=r2_control(dv,dr2)
        # delta A_NL by condition + contrast
        sn=subj(K,f'dnl_{v}',f)
        m_i,li,hi_=mean_ci(sn.INIT.values); m_h,lh,hh=mean_ci(sn.INHIB.values)
        dn=(sn.INIT-sn.INHIB).values; dznl,lnl,hnl=dz_ci(dn); pnl=perm_p(dn)
        # linear-surrogate AIS contrast: does the null reproduce the effect?
        ss=subj(K,f'ais_surr_{v}',f); dsu=(ss.INIT-ss.INHIB).values
        dzs,_,_=dz_ci(dsu); ps=perm_p(dsu)
        # EPR, 1/f slope, lag-1 autocorr
        se=subj(K,f'epr_{v}',f); de=(se.INIT-se.INHIB).values
        sc=subj(K,f'chi_{v}',f); dc=(sc.INIT-sc.INHIB).values
        s1=subj(K,f'r1_{v}',f); d1=(s1.INIT-s1.INHIB).values
        # dip test on subject-level AIS per condition
        dip_i=diptest.diptest(s.INIT.values)[1]; dip_h=diptest.diptest(s.INHIB.values)[1]
        rows.append(dict(variant=v,contrast=label,n=len(s),
            mean_init=s.INIT.mean(),mean_inhib=s.INHIB.mean(),
            dz=dz,dz_lo=lo,dz_hi=hi,p_perm=p,dz_r2adj=dza,p_perm_r2adj=pa,
            dNL_init=m_i,dNL_init_lo=li,dNL_init_hi=hi_,
            dNL_inhib=m_h,dNL_inhib_lo=lh,dNL_inhib_hi=hh,
            dz_dNL=dznl,dz_dNL_lo=lnl,dz_dNL_hi=hnl,p_perm_dNL=pnl,
            dz_surr=dzs,p_perm_surr=ps,
            dz_epr=dz_ci(de)[0],p_epr=perm_p(de),
            dz_chi=dz_ci(dc)[0],p_chi=perm_p(dc),
            dz_r1=dz_ci(d1)[0],p_r1=perm_p(d1),
            dip_p_init=dip_i,dip_p_inhib=dip_h,
            corr_ais_r1=np.corrcoef(K.loc[f,f'ais_{v}'],K.loc[f,f'r1_{v}'])[0,1]))
    # TE contrasts
    for label,c in [('TE frontal->cACC','fc'),('TE cACC->frontal','cf')]:
        s=subj(TE,f'te_{c}_{v}'); dv=(s.INIT-s.INHIB).values
        r2s=K[K.roi.isin(['frontal','cacc'])].groupby(['subject','cond']).r2.mean().unstack()
        dr2=(r2s.INIT-r2s.INHIB).reindex(s.index).values
        dz,lo,hi=dz_ci(dv); p=perm_p(dv); dza,pa=r2_control(dv,dr2)
        sn=subj(TE,f'dte_{c}_{v}')
        m_i,li,hi_=mean_ci(sn.INIT.values); m_h,lh,hh=mean_ci(sn.INHIB.values)
        dn=(sn.INIT-sn.INHIB).values; dznl,lnl,hnl=dz_ci(dn); pnl=perm_p(dn)
        ss=subj(TE,f'tes_{c}_{v}'); dsu=(ss.INIT-ss.INHIB).values
        rows.append(dict(variant=v,contrast=label,n=len(s),
            mean_init=s.INIT.mean(),mean_inhib=s.INHIB.mean(),
            dz=dz,dz_lo=lo,dz_hi=hi,p_perm=p,dz_r2adj=dza,p_perm_r2adj=pa,
            dNL_init=m_i,dNL_init_lo=li,dNL_init_hi=hi_,
            dNL_inhib=m_h,dNL_inhib_lo=lh,dNL_inhib_hi=hh,
            dz_dNL=dznl,dz_dNL_lo=lnl,dz_dNL_hi=hnl,p_perm_dNL=pnl,
            dz_surr=dz_ci(dsu)[0],p_perm_surr=perm_p(dsu),
            dip_p_init=diptest.diptest(s.INIT.values)[1],
            dip_p_inhib=diptest.diptest(s.INHIB.values)[1]))
R=pd.DataFrame(rows)
for v in ['pub','cor']:
    m=R.variant==v
    for col,new in [('p_perm','q_perm'),('p_perm_r2adj','q_r2adj'),('p_perm_dNL','q_dNL')]:
        R.loc[m,new]=multipletests(R.loc[m,col],method='fdr_bh')[1]
R.to_csv('battery_results.csv',index=False)
pd.set_option('display.width',250)
print(R[['variant','contrast','n','dz','dz_lo','dz_hi','p_perm','q_perm','dz_r2adj','q_r2adj',
         'dNL_init','dNL_init_lo','dNL_init_hi','dNL_inhib','dNL_inhib_lo','dNL_inhib_hi',
         'dz_dNL','p_perm_dNL','q_dNL','dz_surr','p_perm_surr']].round(4).to_string(index=False))
print()
print(R[['variant','contrast','dz_epr','p_epr','dz_chi','p_chi','dz_r1','p_r1','corr_ais_r1','dip_p_init','dip_p_inhib']].round(4).to_string(index=False))
