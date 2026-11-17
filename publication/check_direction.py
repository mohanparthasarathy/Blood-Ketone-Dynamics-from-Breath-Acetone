"""Recalculate matched-observation direction scores from saved cohort data."""
from pathlib import Path
import numpy as np
import pandas as pd
root=Path(__file__).resolve().parents[1]
m=pd.read_csv(root/'synthetic_data/outputs/cohort/cohort_metadata.csv')
m=m[m.scenario.isin(['developing','resolving'])].sort_values(['patient_id','breath_time_h']).copy()
m['delta']=m.groupby('patient_id').C_exh_plateau_ppm.diff()
y=m.dBHB_dt_mM_h>0
tr=m.split=='train'
x=m.C_exh_plateau_ppm
thresholds=np.linspace(x[tr].min(),x[tr].max(),200)
high=np.array([((x[tr]>t)==y[tr]).mean() for t in thresholds])
low=np.array([((x[tr]<t)==y[tr]).mean() for t in thresholds])
sign=1 if high.max()>=low.max() else -1
t=thresholds[(high if sign==1 else low).argmax()]
v=m.delta.notna()
a=int(((sign*x[v]>sign*t)==y[v]).sum())
b=int(((m.delta[v]>0)==y[v]).sum())
assert (a,b,int(v.sum()))==(139,211,240), (a,b,int(v.sum()))
pd.DataFrame([dict(current_correct=a,history_correct=b,n=int(v.sum()),threshold_ppm=t)]).to_csv(root/'publication/direction_scores.csv',index=False)
print(f'Current plateau: {a}/240; preceding-breath change: {b}/240; threshold: {t:.12g} ppm')
