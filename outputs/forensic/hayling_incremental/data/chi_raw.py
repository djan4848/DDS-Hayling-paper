"""El chi reportado como '0.405' se calculo sobre el RESIDUO, no sobre el EEG crudo.
   Aqui se calcula tambien sobre el crudo, con el mismo procedimiento, para poder
   afirmar en el cierre lo que corresponde."""
import numpy as np, pandas as pd, pickle, warnings; warnings.filterwarnings('ignore')
from scipy.signal import welch
from joblib import Parallel, delayed
from specparam import SpectralModel
from statsmodels.stats.multitest import multipletests
rng=np.random.default_rng(2026); NPERM=10000
K=pd.read_csv('fase2_unit.csv'); D=pickle.load(open('fase2_series.pkl','rb'))
def perm_p(d,n=NPERM):
    obs=abs(d.mean()); sg=rng.choice([-1.,1.],size=(n,len(d)))
    return float(((np.abs((sg*d).mean(1))>=obs-1e-15).sum()+1)/(n+1))
def chi_table(series,tag):
    f0=welch(series[0],fs=250,nperseg=len(series[0]))[0]
    P=np.array(Parallel(n_jobs=18,batch_size=200)(
        delayed(lambda y: welch(y,fs=250,nperseg=len(y))[1])(y) for y in series))
    rows=[]
    for k,ix in K.groupby(['subject','cond','roi']).indices.items():
        fm=SpectralModel(peak_width_limits=[2,12],max_n_peaks=3,min_peak_height=0.1,
                         aperiodic_mode='fixed',verbose=False)
        try:
            fm.fit(f0,P[ix].mean(0),[2.,40.])
            chi=float(np.atleast_1d(fm.get_params('aperiodic'))[-1])
            r2=float(fm.results.metrics.results['gof_rsquared'])
        except Exception: chi,r2=np.nan,np.nan
        rows.append(dict(subject=k[0],cond=k[1],roi=k[2],chi=chi,gof=r2,src=tag))
    return pd.DataFrame(rows)
out=[]
for series,tag in [(D['RES'],'residuo corregido'),(D['ERP'],'EEG crudo')]:
    C=chi_table(series,tag); out.append(C)
    rows=[]
    for roi in ['frontal','cacc','lh','rh']:
        s=C[C.roi==roi].pivot_table(index='subject',columns='cond',values='chi').dropna()
        d=(s.INIT-s.INHIB).values
        rows.append(dict(src=tag,roi=roi,n=len(s),chi_init=s.INIT.mean(),chi_inhib=s.INHIB.mean(),
                         dz=d.mean()/d.std(ddof=1),p_perm=perm_p(d)))
    R=pd.DataFrame(rows); R['q_bh4']=multipletests(R.p_perm,method='fdr_bh')[1]
    print(f"\n=== chi aperiodico ({tag}) — GOF R2 mediana {C.gof.median():.3f} ===")
    print(R.round(4).to_string(index=False))
    R.to_csv(f'chi_{tag.split()[0]}.csv',index=False)
pd.concat(out).to_csv('chi_both.csv',index=False)
print('CHI_DONE')
