"""Unified sample-level TCGA TPM analysis. Run: python analyze_tcga_myc_lncdzank1.py
Input: existing Xena STAR matrices (already log2(TPM+1)); no counts or TMM.
Optional: --input-dir PATH --output-dir PATH. Dependencies in requirements.txt.
"""
from pathlib import Path
import argparse, gzip, json, hashlib, sys, platform
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser()
p.add_argument('--input-dir',type=Path,default=Path('/Users/darius57/Documents/Codex/2026-09-23/yo/work'))
p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
a=p.parse_args(); W=a.input_dir; O=a.output_dir; O.mkdir(parents=True,exist_ok=True)
log=[]
def say(s):
    print(s,flush=True); log.append(str(s))
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
ann=pd.read_csv(W/'probemap.tsv',sep='\t',dtype=str)
ann['stable_id']=ann.id.str.split('.').str[0]
target=ann[ann.stable_id=='ENSG00000230010']
myc=ann[ann.gene=='MYC']
if len(target)!=1:
    candidates=ann[ann.gene.str.contains('DZANK1',case=False,na=False)]
    say('STOP: target cannot be matched uniquely. DZANK1 annotation candidates:\n'+candidates.to_string(index=False))
    candidates.to_csv(O/'DZANK1_annotation_candidates.csv',index=False)
    (O/'QC_log.txt').write_text('\n'.join(log));sys.exit(1)
assert len(myc)==1,'MYC must match exactly once'
genes={'MYC':myc.iloc[0], 'LncDZANK1':target.iloc[0]}
provenance={'annotation':{'path':str(W/'probemap.tsv'),'sha256':sha(W/'probemap.tsv')},'matrices':[], 'scope':'CRC: primary tumors (01), matching prior CRC primary tumor table; ESCA: primary + metastatic (01,06), matching prior ESCA tumor sample set.'}
frames=[]; inventory=[]; audit=[]
for cancer in ['COAD','READ','ESCA']:
    meta=json.loads((W/f'{cancer.lower()}.metadata.json').read_text())
    assert meta['unit'].lower()=='log2(tpm+1)', 'Unexpected input units; refusing implicit conversion'
    path=W/f'{cancer.lower()}.star_tpm.tsv.gz'; found={k:[] for k in genes}
    with gzip.open(path,'rt') as f:
        ids=next(f).rstrip('\n').split('\t')[1:]
        assert len(ids)==len(set(ids)),f'Duplicated sample IDs in {cancer}'
        for line in f:
            parts=line.rstrip('\n').split('\t'); gid=parts[0]
            for k,g in genes.items():
                if gid.split('.')[0]==g.stable_id:
                    found[k].append((gid,pd.to_numeric(pd.Series(parts[1:]),errors='coerce').to_numpy()))
    d=pd.DataFrame({'Sample_ID':ids,'Cancer_Type':cancer})
    codes=d.Sample_ID.str[13:15]
    names={'01':'Primary Solid Tumor','02':'Recurrent Solid Tumor','06':'Metastatic','11':'Solid Tissue Normal'}
    assert codes.isin(names).all()
    d['Sample_Type']=codes.map(names)
    for k,g in genes.items():
        assert len(found[k])==1,(cancer,k,found[k])
        gid,values=found[k][0];assert gid==g.id, 'Matrix/annotation version mismatch'
        assert len(values)==len(d)
        if np.isinf(values).any() or (values[np.isfinite(values)]<0).any():raise ValueError('Invalid expression values')
        d[k+'_log2TPM1']=values
        d[k+'_TPM']=np.expm1(values*np.log(2))
    allowed=['01','06'] if cancer=='ESCA' else ['01']
    d['Included_tumor_scope']=codes.isin(allowed)
    d['Both_genes_available']=d[['MYC_log2TPM1','LncDZANK1_log2TPM1']].notna().all(axis=1)
    inventory.append(d[['Sample_ID','Cancer_Type','Sample_Type','Included_tumor_scope','Both_genes_available']])
    for typ,sub in d.groupby('Sample_Type'):
        audit.append({'Cancer_Type':cancer,'Sample_Type':typ,'Matrix_n':len(sub),'Eligible_tumor_n':int(sub.Included_tumor_scope.sum()),'Final_n':int((sub.Included_tumor_scope&sub.Both_genes_available).sum())})
    frames.append(d[d.Included_tumor_scope].copy())
    provenance['matrices'].append({'cancer':cancer,'path':str(path),'sha256':sha(path),'url':f'https://gdc.xenahubs.net/download/TCGA-{cancer}.star_tpm.tsv.gz','metadata':meta,'matched_genes':{k:found[k][0][0] for k in genes}})
all_d=pd.concat(frames,ignore_index=True)
pd.concat(inventory).to_csv(O/'sample_inclusion_audit.csv',index=False)
pd.DataFrame(audit).to_csv(O/'sample_type_counts.csv',index=False)
cohorts={'CRC':all_d[all_d.Cancer_Type.isin(['COAD','READ'])].copy(),'ESCA':all_d[all_d.Cancer_Type=='ESCA'].copy()}
# Validate ESCA against saved previous sample membership; never use previous groups.
prior=W.parent/'outputs/gobp_185_tumors/sample_groups.tsv'
if prior.exists():
    prev=pd.read_csv(prior,sep='\t').iloc[:,0]
    assert set(prev)==set(cohorts['ESCA'].Sample_ID),'ESCA differs from previous tumor scope'
    provenance['previous_ESCA_membership_verified']=True
complete=all_d[all_d.Both_genes_available]
say('FINAL ELIGIBLE SAMPLES (both genes available)\nCRC:')
for c in ['COAD','READ']:say(f'{c} tumor n = {sum(complete.Cancer_Type==c)}')
say(f'Combined CRC tumor n = {sum(complete.Cancer_Type.isin(["COAD","READ"]))}')
es=complete[complete.Cancer_Type=='ESCA'];say('ESCA:')
say(f'Primary tumor n = {sum(es.Sample_Type=="Primary Solid Tumor")}')
say(f'Metastatic tumor n = {sum(es.Sample_Type=="Metastatic")}');say(f'Total ESCA tumor n = {len(es)}')
say('\nGENE MATCHES (user target name Lnc-DZANK1; annotation symbol retained)')
for k,g in genes.items():say(f'{k}: gene symbol = {g.gene}; Ensembl ID = {g.id}')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8})
colors=['#9AABB5','#45677D']; summaries=[]; descriptions=[]; qc={}; ranges=[]; results=[]
ymax=float(complete.LncDZANK1_log2TPM1.max())*1.24
xlims=(float(complete.MYC_log2TPM1.min())-.3,float(complete.MYC_log2TPM1.max())+.3)
def savefig(fig,stem):
    fig.savefig(O/(stem+'.pdf'));fig.savefig(O/(stem+'.png'),dpi=300);plt.close(fig)
def fp(v):return f'{v:.3g}'
for cohort,eligible in cohorts.items():
    d=eligible[eligible.Both_genes_available].sort_values(['MYC_log2TPM1','Sample_ID'],kind='stable').reset_index(drop=True)
    assert not d.Sample_ID.duplicated().any()
    assert d.Sample_Type.isin(['Primary Solid Tumor','Metastatic']).all()
    assert np.isfinite(d[['MYC_log2TPM1','LncDZANK1_log2TPM1']].values).all()
    n=len(d); cut=n//2;assert cut>0
    d['MYC_Group']=np.where(np.arange(n)<cut,'MYC-low','MYC-high')
    cols=['Sample_ID','Cancer_Type','Sample_Type','MYC_TPM','MYC_log2TPM1','LncDZANK1_TPM','LncDZANK1_log2TPM1','MYC_Group']
    d[cols].to_csv(O/f'TCGA_{cohort}_MYC_LncDZANK1_analysis.csv',index=False)
    low=d.loc[d.MYC_Group=='MYC-low','LncDZANK1_log2TPM1'].to_numpy();high=d.loc[d.MYC_Group=='MYC-high','LncDZANK1_log2TPM1'].to_numpy()
    test=stats.mannwhitneyu(low,high,alternative='two-sided',method='asymptotic',use_continuity=True)
    corr=stats.spearmanr(d.MYC_log2TPM1,d.LncDZANK1_log2TPM1,alternative='two-sided')
    summary={'Cohort':cohort,'Total_n':n,'MYC_low_n':len(low),'MYC_high_n':len(high),'LncDZANK1_median_low':float(np.median(low)),'LncDZANK1_median_high':float(np.median(high)),'Wilcoxon_P':float(test.pvalue),'Spearman_rho':float(corr.statistic),'Spearman_P':float(corr.pvalue)}
    summaries.append(summary)
    for label,v in [('MYC-low',low),('MYC-high',high)]:
        q1,median,q3=np.quantile(v,[.25,.5,.75],method='linear')
        descriptions.append({'Cohort':cohort,'MYC_Group':label,'n':len(v),'Median':median,'Q1':q1,'Q3':q3,'IQR':q3-q1,'Mean':v.mean(),'SD':v.std(ddof=1)})
    med=float(d.MYC_log2TPM1.median());left=float(d.MYC_log2TPM1.iloc[cut-1]);right=float(d.MYC_log2TPM1.iloc[cut]); tied=left==right
    patient_counts=d.Sample_ID.str[:12].value_counts()
    qc[cohort]={'Eligible_tumor_n':len(eligible),'Total_tumor_samples':n,'MYC_low':len(low),'MYC_high':len(high),'Missing_MYC':int(eligible.MYC_log2TPM1.isna().sum()),'Missing_LncDZANK1':int(eligible.LncDZANK1_log2TPM1.isna().sum()),'Excluded_missing_either':len(eligible)-n,'MYC_median':med,'MYC_low_max':left,'MYC_high_min':right,'Tie_crosses_split_boundary':bool(tied),'Samples_equal_MYС_median':int((d.MYC_log2TPM1==med).sum()),'Boundary_value_frequency_low_end':int((d.MYC_log2TPM1==left).sum()),'Boundary_value_frequency_high_start':int((d.MYC_log2TPM1==right).sum()),'Unique_patients':len(patient_counts),'Patients_with_multiple_samples':patient_counts[patient_counts>1].to_dict(),'No_normal_samples':True,'No_duplicated_sample_IDs':True,'No_NA_in_tests':True}
    say(f'\n{cohort}: MYC-low n = {len(low)}; MYC-high n = {len(high)}')
    say(f'MYC median = {med:.9g}; low maximum = {left:.9g}; high minimum = {right:.9g}; boundary tie = {tied}; samples equal median = {qc[cohort]["Samples_equal_MYС_median"]}')
    for k,g in genes.items():
        rr={'Cohort':cohort,'Display_name':k,'Gene_symbol':g.gene,'Ensembl_ID':g.id,'TPM_min':d[k+'_TPM'].min(),'TPM_max':d[k+'_TPM'].max(),'log2TPM1_min':d[k+'_log2TPM1'].min(),'log2TPM1_max':d[k+'_log2TPM1'].max()};ranges.append(rr);say(str(rr))
    say(pd.DataFrame(descriptions[-2:]).to_string(index=False))
    say(f'Two-sided Wilcoxon rank-sum: U = {test.statistic}; P = {test.pvalue:.12g}')
    say(f'Spearman rho = {corr.statistic:.9g}; P = {corr.pvalue:.12g}; n = {n}')
    # Same plotting parameters and axis ranges in both cohorts. Violin KDE ends at observations.
    fig,ax=plt.subplots(figsize=(4.7,4.8));fig.subplots_adjust(left=.21,right=.97,bottom=.18,top=.86)
    viol=ax.violinplot([low,high],positions=[0,1],widths=.72,showextrema=False)
    for b,c in zip(viol['bodies'],colors):b.set_facecolor(c);b.set_edgecolor(c);b.set_alpha(.32)
    rng=np.random.default_rng(24092026)
    for pos,v,c in zip([0,1],[low,high],colors):ax.scatter(pos+rng.uniform(-.16,.16,len(v)),v,s=9,c=c,alpha=.52,edgecolors='none',zorder=2)
    ax.boxplot([low,high],positions=[0,1],widths=.25,showfliers=False,patch_artist=True,boxprops={'facecolor':'white','edgecolor':'#333333','linewidth':.9},medianprops={'color':'#222222','linewidth':1.3},whiskerprops={'linewidth':.9},capprops={'linewidth':.9},zorder=3)
    ax.set_xticks([0,1],[f'MYC-low\n(n = {len(low)})',f'MYC-high\n(n = {len(high)})'])
    ax.set_xlim(-.6,1.6);ax.set_ylim(-.035*ymax,ymax);ax.set_ylabel('Lnc-DZANK1 expression\nlog2(TPM + 1)');ax.set_title('TCGA-'+cohort,pad=26,fontweight='bold',fontsize=13)
    ax.text(.5,1.035,f'Two-sided Wilcoxon P = {fp(test.pvalue)}',transform=ax.transAxes,ha='center',fontsize=10)
    savefig(fig,f'TCGA_{cohort}_MYC_high_low_LncDZANK1')
    fig,ax=plt.subplots(figsize=(4.9,4.8));fig.subplots_adjust(left=.20,right=.97,bottom=.18,top=.87)
    x=d.MYC_log2TPM1.to_numpy();y=d.LncDZANK1_log2TPM1.to_numpy()
    ax.scatter(x,y,s=16,color=colors[1],alpha=.5,edgecolors='none')
    slope,intercept=np.polyfit(x,y,1);xx=np.linspace(x.min(),x.max(),150)
    ax.plot(xx,slope*xx+intercept,color='#333333',lw=1.3,label='Linear fit (visual guide)')
    ax.set(xlim=xlims,ylim=(-.035*ymax,ymax),xlabel='MYC expression\nlog2(TPM + 1)',ylabel='Lnc-DZANK1 expression\nlog2(TPM + 1)')
    ax.set_title('TCGA-'+cohort,pad=12,fontweight='bold',fontsize=13)
    ax.text(.03,.97,f'Spearman ρ = {corr.statistic:.3f}\nP = {fp(corr.pvalue)}; n = {n}',transform=ax.transAxes,va='top',fontsize=10)
    ax.legend(loc='upper right',frameon=False,fontsize=8)
    savefig(fig,f'TCGA_{cohort}_MYC_LncDZANK1_correlation')
    cancer_name='colorectal cancer' if cohort=='CRC' else 'esophageal cancer'
    if test.pvalue<.05 and np.median(high)>np.median(low):
        clause=f'Lnc-DZANK1 expression was significantly higher in MYC-high than in MYC-low tumors in the TCGA {cancer_name} cohort'
    elif test.pvalue<.05 and np.median(high)<np.median(low):
        clause=f'Lnc-DZANK1 expression was significantly lower in MYC-high than in MYC-low tumors in the TCGA {cancer_name} cohort'
    elif test.pvalue<.05:
        clause=f'Lnc-DZANK1 expression distributions differed significantly between MYC groups in the TCGA {cancer_name} cohort'
    else:
        clause=f'No statistically significant difference in Lnc-DZANK1 expression between MYC-high and MYC-low tumors was detected in the TCGA {cancer_name} cohort'
    results.append(f'{clause} (high n = {len(high)}, low n = {len(low)}; median log2(TPM + 1), {np.median(high):.4f} versus {np.median(low):.4f}; two-sided Wilcoxon rank-sum P = {fp(test.pvalue)}). Spearman analysis yielded rho = {corr.statistic:.3f} (P = {fp(corr.pvalue)}, n = {n}).')
pd.DataFrame(summaries).to_csv(O/'TCGA_MYC_LncDZANK1_summary.csv',index=False)
pd.DataFrame(descriptions).to_csv(O/'TCGA_MYC_LncDZANK1_descriptive_statistics.csv',index=False)
pd.DataFrame(ranges).to_csv(O/'gene_annotation_expression_ranges.csv',index=False)
provenance['software']={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__}
(O/'source_provenance.json').write_text(json.dumps(provenance,indent=2))
(O/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in provenance['software'].items() if k!='python')+'\n')
qc['shared_methods']={'preprocessing':'Input log2(TPM+1), TPM back-transformed only for export; retain zero expression; complete-case samples only.','grouping':'Ascending MYC then Sample_ID; first floor(n/2) low, remaining high.','wilcoxon':'Two-sided Mann-Whitney U / Wilcoxon rank-sum, asymptotic normal approximation, tie correction and continuity correction.','spearman':'Two-sided Spearman using average ranks for ties; asymptotic t-based P value.','quartiles':'Linear interpolation (R type 7).','P_values':'Nominal, unadjusted; significance threshold 0.05.','scatter_line':'OLS visual guide, not Spearman fit or causal inference.'}
(O/'QC.json').write_text(json.dumps(qc,indent=2))
(O/'Results_template.txt').write_text('\n\n'.join(results)+'\n\nP values are nominal and unadjusted. These sample-level associations do not establish causality. Multiple tumor samples from some patients were retained to preserve the requested cohort; observations are not all independent at patient level.\n')
say('\nFINAL QC')
for c in cohorts:
    q=qc[c];say(f'{c}:\nTotal tumor samples = {q["Total_tumor_samples"]}\nMYC-low = {q["MYC_low"]}\nMYC-high = {q["MYC_high"]}\nMissing MYC = {q["Missing_MYC"]}\nMissing Lnc-DZANK1 = {q["Missing_LncDZANK1"]}\nUnique patients = {q["Unique_patients"]}; patients with multiple samples = {len(q["Patients_with_multiple_samples"])}')
say('Confirmed: no normal samples included; no duplicated sample IDs; no NA values used in statistical tests; CRC and ESCA used the same preprocessing and statistical methods.')
say('\nRESULTS\n'+'\n\n'.join(results))
(O/'QC_log.txt').write_text('\n'.join(log)+'\n')
