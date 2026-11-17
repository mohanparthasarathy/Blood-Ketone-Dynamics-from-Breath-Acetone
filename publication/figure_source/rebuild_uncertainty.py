"""Render saved profile results; no model fitting or alteration of results."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
p=pd.read_csv(root/'gamma_profiles.csv'); b=pd.read_csv(root/'state_uncertainty_bands.csv')
plt.rcParams.update({'font.family':'serif','font.size':10,'pdf.fonttype':42})
fig,ax=plt.subplots(2,2,figsize=(8.0,6.4),layout='constrained')
for col,pid in enumerate([5,14]):
 q=p[p.patient_id==pid].sort_values('Gamma_A'); t=b[b.patient_id==pid].sort_values('time_h');a=ax[0,col]
 a.plot(q.Gamma_A,q.delta_chi_square,'o-',color='#0072B2',mfc='white',ms=5,label='Profile refits')
 ok=q.accepted.astype(str).str.lower().isin(['1','true']);a.scatter(q.loc[ok,'Gamma_A'],q.loc[ok,'delta_chi_square'],s=25,color='#0072B2',zorder=3,label='Retained refits')
 a.axhline(3.8415,color='#666666',ls='--',lw=1,label='Retention cutoff')
 a.set(xlabel=r'$\Gamma_A$',ylabel=r'$\Delta\chi^2$',title=('Stable patient 5' if pid==5 else 'Developing patient 14'))
 if col==0:a.legend(fontsize=8,loc='upper left')
 a=ax[1,col];x=t.time_h.to_numpy();a.fill_between(x,t.BHB_lower_mM.to_numpy(),t.BHB_upper_mM.to_numpy(),color='#b9d9ee',label='Profile sensitivity envelope')
 a.plot(x,t.BHB_truth_mM,'k-',lw=1.7,label='Simulated BHB')
 a.plot(x,t.BHB_point_mM,color='#0072B2',lw=1.5,label='Minimum-misfit profile refit')
 h=t.temporal_holdout==1;a.scatter(t.loc[h,'time_h'],t.loc[h,'BHB_truth_mM'],s=24,color='black',zorder=5,label='Withheld truth')
 a.set(xlabel='Time (h)',ylabel='BHB (mmol/L)',xlim=(.5,8))
 if col==0:
  a.set_ylim(.14,.24);a.legend(fontsize=7,loc='upper left')
 else:
  a.set_ylim(0,3.15)
  for v in [1.5,3]:a.axhline(v,ls='--',lw=.9,color='#a66b43')
for a in ax.flat:a.grid(alpha=.2);a.spines[['top','right']].set_visible(False)
fig.savefig(root.parent/'reconstruction_uncertainty.pdf',bbox_inches='tight')
fig.savefig(root.parent/'reconstruction_uncertainty.png',dpi=180,bbox_inches='tight')
