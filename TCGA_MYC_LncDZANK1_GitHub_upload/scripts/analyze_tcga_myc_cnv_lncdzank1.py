"""Analyze MYC GISTIC discrete CNV versus Lnc-DZANK1 in TCGA CRC and ESCA."""
from pathlib import Path
import argparse, gzip, json, platform, hashlib
import numpy as np, pandas as pd, requests
from scipy import stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=argparse.ArgumentParser(); P.add_argument('--input-dir',type=Path,default=Path('/Users/darius57/Documents/Codex/2026-09-23/yo/work')); P.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent); A=P.parse_args(); W,O=A.input_dir,A.output_dir; O.mkdir(parents=True,exist_ok=True)
STUDIES={'CRC':{'study':'coadread_tcga_pan_can_atlas_2018','cancers':['COAD','READ']},'ESCA':{'study':'esca_tcga_pan_can_atlas_2018','cancers':['ESCA']}}
ENCODING={-2:'deep deletion',-1:'shallow deletion',0:'diploid',1:'gain',2:'high-level amplification'}
def norm(s): return '-'.join(str(s).split('-')[:4])
def patient(s): return '-'.join(str(s).split('-')[:3])
ann=pd.read_csv(W/'probemap.tsv',sep='\t',dtype=str); ann['stable']=ann.id.str.split('.').str[0]; target=ann[ann.stable.eq('ENSG00000230010')]; myc=ann[ann.gene.eq('MYC')]
if len(target)!=1 or len(myc)!=1: raise RuntimeError('Exact MYC or ENSG00000230010 annotation unavailable')
genes={'MYC':myc.iloc[0],'LncDZANK1':target.iloc[0]}
def expr(cancer):
    path=W/f'{cancer.lower()}.star_tpm.tsv.gz'; found={k:None for k in genes}
    with gzip.open(path,'rt') as f:
        ids=next(f).rstrip('\n').split('\t')[1:]
        for line in f:
            z=line.rstrip('\n').split('\t'); gid=z[0].split('.')[0]
            for k,g in genes.items():
                if gid==g.stable: found[k]=pd.to_numeric(pd.Series(z[1:]),errors='coerce').to_numpy(float)
    if any(v is None for v in found.values()): raise RuntimeError('Missing expression gene')
    d=pd.DataFrame({'TCGA_sample_ID':ids,'Cancer':cancer}); d=d[d.TCGA_sample_ID.str[13:15].eq('01')].copy()
    for k,v in found.items(): d[k+'_expression']=v[d.index.to_numpy()]
    d['TCGA_sample_ID']=d.TCGA_sample_ID.map(norm); d['TCGA_patient_ID']=d.TCGA_sample_ID.map(patient)
    d=d[d[['MYC_expression','LncDZANK1_expression']].notna().all(axis=1)].sort_values('TCGA_sample_ID').drop_duplicates('TCGA_patient_ID')
    return d[['TCGA_sample_ID','TCGA_patient_ID','Cancer','MYC_expression','LncDZANK1_expression']]
def cnv(study):
    profile=f'{study}_gistic'; sl=f'{study}_cna'; url=f'https://www.cbioportal.org/api/molecular-profiles/{profile}/molecular-data/fetch'
    r=requests.post(url,json={'entrezGeneIds':[4609],'sampleListId':sl},timeout=120); r.raise_for_status(); x=r.json(); vals={int(a['value']) for a in x}
    if not vals.issubset(set(ENCODING)): raise RuntimeError(f'Unexpected GISTIC values {vals}')
    (O/f'{study}_MYC_GISTIC_raw.json').write_text(json.dumps(x,indent=2))
    return pd.DataFrame({'CNV_sample_ID':[norm(a['sampleId']) for a in x],'TCGA_patient_ID':[patient(a['sampleId']) for a in x],'MYC_GISTIC':[int(a['value']) for a in x]})
def stat(d,comparison,g1,g2,col='LncDZANK1_expression'):
    a=d.loc[d.CNV_group.eq(g1),col].to_numpy(float); b=d.loc[d.CNV_group.eq(g2),col].to_numpy(float)
    out={'Comparison':comparison,'Group_1':g1,'Group_2':g2,'N_group_1':len(a),'N_group_2':len(b),'Median_group_1':np.nan,'Median_group_2':np.nan,'Mean_group_1':np.nan,'Mean_group_2':np.nan,'IQR_group_1':np.nan,'IQR_group_2':np.nan,'Effect_size_rank_biserial':np.nan,'Median_ratio_TPM_scale':np.nan,'log2FC_TPM_scale':np.nan,'P_value':np.nan}
    if not len(a) or not len(b): return out
    u=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic',use_continuity=True); m1,m2=float(np.median(a)),float(np.median(b)); t1,t2=np.expm1(np.log(2)*m1),np.expm1(np.log(2)*m2); ratio=t1/t2 if t2>0 else np.nan
    out.update({'Median_group_1':m1,'Median_group_2':m2,'Mean_group_1':float(a.mean()),'Mean_group_2':float(b.mean()),'IQR_group_1':float(np.quantile(a,.75)-np.quantile(a,.25)),'IQR_group_2':float(np.quantile(b,.75)-np.quantile(b,.25)),'Effect_size_rank_biserial':2*float(u.statistic)/(len(a)*len(b))-1,'Median_ratio_TPM_scale':ratio,'log2FC_TPM_scale':float(np.log2(ratio)) if ratio>0 else np.nan,'P_value':float(u.pvalue)})
    return out
def primary_plot(d,c):
    gs=['MYC-non-amplified','MYC-amplified']; arr=[d.loc[d.CNV_group.eq(g),'LncDZANK1_expression'].to_numpy() for g in gs]; fig,ax=plt.subplots(figsize=(5,4.6)); vp=ax.violinplot(arr,positions=[0,1],showextrema=False,widths=.75)
    for b,col in zip(vp['bodies'],['#9AABB5','#45677D']): b.set_facecolor(col); b.set_edgecolor(col); b.set_alpha(.32)
    ax.boxplot(arr,positions=[0,1],widths=.25,showfliers=False,patch_artist=True,boxprops={'facecolor':'white'},medianprops={'color':'#222','linewidth':1.3}); rng=np.random.default_rng(20260929)
    for i,(v,col) in enumerate(zip(arr,['#9AABB5','#45677D'])): ax.scatter(i+rng.uniform(-.15,.15,len(v)),v,s=11,color=col,alpha=.55,edgecolors='none')
    q=stat(d,'primary_amp_vs_nonamp','MYC-amplified','MYC-non-amplified'); ax.set_xticks([0,1],[f'MYC-non-amplified\n(n={len(arr[0])})',f'MYC-amplified\n(n={len(arr[1])})']); ax.set_ylabel('Lnc-DZANK1 expression\nlog2(TPM + 1)'); ax.set_title(f'TCGA-{c}',fontweight='bold',pad=22); ax.text(.5,1.055,f"Wilcoxon P = {q['P_value']:.3g}",transform=ax.transAxes,ha='center'); ax.spines[['top','right']].set_visible(False); fig.subplots_adjust(top=.80,left=.20,bottom=.18,right=.97); stem=O/f'TCGA_{c}_MYC_CNV_LncDZANK1_primary'; fig.savefig(str(stem)+'.pdf'); fig.savefig(str(stem)+'.png',dpi=300); plt.close(fig)
def dosage_plot(d,c):
    cats=[-2,-1,0,1,2]; valid=[(x,d.loc[d.MYC_GISTIC.eq(x),'LncDZANK1_expression'].to_numpy()) for x in cats if (d.MYC_GISTIC==x).any()]; fig,ax=plt.subplots(figsize=(6,4.6)); arr=[v for x,v in valid]; pos=[x for x,v in valid]; vp=ax.violinplot(arr,positions=pos,showextrema=False,widths=.75)
    for b in vp['bodies']: b.set_facecolor('#6D8796'); b.set_edgecolor('#45677D'); b.set_alpha(.32)
    ax.boxplot(arr,positions=pos,widths=.25,showfliers=False,patch_artist=True,boxprops={'facecolor':'white'},medianprops={'color':'#222','linewidth':1.3}); rng=np.random.default_rng(20260930)
    for x,v in valid: ax.scatter(x+rng.uniform(-.14,.14,len(v)),v,s=11,color='#45677D',alpha=.55,edgecolors='none')
    ax.set_xticks(cats,['-2','-1','0','+1','+2']); ax.set_xlabel('MYC GISTIC copy-number category'); ax.set_ylabel('Lnc-DZANK1 expression\nlog2(TPM + 1)'); ax.set_title(f'TCGA-{c}',fontweight='bold'); ax.spines[['top','right']].set_visible(False); fig.tight_layout(); stem=O/f'TCGA_{c}_MYC_GISTIC_dosage_LncDZANK1'; fig.savefig(str(stem)+'.pdf'); fig.savefig(str(stem)+'.png',dpi=300); plt.close(fig)
all_d=[]; sums=[]; qcs=[]; source={'source':'cBioPortal TCGA PanCancer Atlas GISTIC discrete molecular profile','encoding':ENCODING,'MYC_entrezGeneId':4609,'study_profiles':{}}
for cohort,cfg in STUDIES.items():
    e=pd.concat([expr(c) for c in cfg['cancers']],ignore_index=True); n=cnv(cfg['study']); d=e.merge(n[['TCGA_patient_ID','MYC_GISTIC']].drop_duplicates('TCGA_patient_ID'),on='TCGA_patient_ID',how='inner'); d['Cancer']=cohort; d['MYC_CNV_group']=d.MYC_GISTIC.map({2:'MYC-amplified',1:'MYC-gain',0:'MYC-diploid',-1:'MYC-shallow-deletion',-2:'MYC-deep-deletion'}); d.to_csv(O/f'TCGA_{cohort}_MYC_CNV_LncDZANK1_analysis.csv',index=False); all_d.append(d)
    d['CNV_group']=np.where(d.MYC_GISTIC.eq(2),'MYC-amplified','MYC-non-amplified'); rows=[stat(d,'primary_amp_vs_nonamp','MYC-amplified','MYC-non-amplified')]
    x=d.copy(); x['CNV_group']=np.where(x.MYC_GISTIC.isin([1,2]),'MYC-gain','MYC-non-gain'); rows.append(stat(x,'sensitivity_gain_vs_non_gain','MYC-gain','MYC-non-gain'))
    x=d[d.MYC_GISTIC.isin([0,2])].copy(); x['CNV_group']=np.where(x.MYC_GISTIC.eq(2),'MYC-amplified','MYC-diploid'); rows.append(stat(x,'sensitivity_amp_vs_diploid','MYC-amplified','MYC-diploid'))
    rows.append(stat(d,'primary_amp_vs_nonamp_MYC_RNA','MYC-amplified','MYC-non-amplified','MYC_expression'))
    for metric in ['MYC_expression','LncDZANK1_expression']:
        z=d[['MYC_GISTIC',metric]].dropna(); sp=stats.spearmanr(z.MYC_GISTIC,z[metric]); rows.append({'Comparison':'dosage_spearman_'+metric,'Group_1':'-2_to_+2','Group_2':'','N_group_1':len(z),'N_group_2':'','Median_group_1':'','Median_group_2':'','Mean_group_1':'','Mean_group_2':'','IQR_group_1':'','IQR_group_2':'','Effect_size_rank_biserial':float(sp.statistic),'Median_ratio_TPM_scale':'','log2FC_TPM_scale':'','P_value':float(sp.pvalue)})
    pd.DataFrame(rows).assign(Cancer=cohort).to_csv(O/f'TCGA_{cohort}_MYC_CNV_summary.csv',index=False); primary_plot(d,cohort); dosage_plot(d,cohort); all_d[-1]['MYC_CNV_group']=all_d[-1].MYC_GISTIC.map({2:'MYC-amplified',1:'MYC-gain',0:'MYC-diploid',-1:'MYC-shallow-deletion',-2:'MYC-deep-deletion'})
    qcs.append({'Cancer':cohort,'RNA_primary_samples_before_patient_dedup':len(e),'CNV_samples_returned':len(n),'RNA_CNV_matched_patients':len(d),'+2_amplification':int((d.MYC_GISTIC==2).sum()),'+1_gain':int((d.MYC_GISTIC==1).sum()),'diploid':int((d.MYC_GISTIC==0).sum()),'deletion':int(d.MYC_GISTIC.isin([-2,-1]).sum()),'CNV_value_counts':{str(k):int(v) for k,v in d.MYC_GISTIC.value_counts().sort_index().items()},'No_normal_samples':True,'No_duplicate_patient_IDs':not d.TCGA_patient_ID.duplicated().any()}); source['study_profiles'][cohort]={'study':cfg['study'],'profile':f"{cfg['study']}_gistic"}
pd.concat(all_d,ignore_index=True).to_csv(O/'TCGA_MYC_CNV_LncDZANK1_analysis.csv',index=False); pd.concat([pd.read_csv(O/f'TCGA_{c}_MYC_CNV_summary.csv') for c in STUDIES],ignore_index=True).to_csv(O/'TCGA_MYC_CNV_LncDZANK1_summary.csv',index=False)
source['software']={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}; (O/'TCGA_MYC_CNV_source.json').write_text(json.dumps(source,indent=2)); (O/'TCGA_MYC_CNV_QC.json').write_text(json.dumps(qcs,indent=2)); (O/'TCGA_MYC_CNV_QC_log.txt').write_text(json.dumps(qcs,indent=2)); print(json.dumps(qcs,indent=2))
