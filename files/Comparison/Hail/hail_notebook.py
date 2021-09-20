#%%
import os
import time
os.chdir(os.path.dirname(os.path.abspath(__file__)))
timer = time.perf_counter()
import hail as hl
from hail.genetics import reference_genome
hl.init()
# %%
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()
# %%
prefix = "mouse_hs1940"
# %%
# %%
hl.import_vcf(f'data/{prefix}.vcf.bgz', reference_genome='GRCm38').write(f'data/{prefix}.mt', overwrite=True)
#%%
table = hl.import_table(f'data/{prefix}.tsv', impute=True).key_by('Sample')
#%%
mt = hl.read_matrix_table(f'data/{prefix}.mt')
#%%
mt = mt.annotate_cols(pheno = table[mt.s])
#%%
mt = hl.sample_qc(mt)
#%%
#mt = mt.filter_cols((mt.sample_qc.dp_stats.mean >= 4) & (mt.sample_qc.call_rate >= 0.97))
#mt.describe()
#%%
#ab = mt.AD[1] / hl.sum(mt.AD)
#filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
#                        (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
#                        (mt.GT.is_hom_var() & (ab >= 0.9)))
#%%
#mt = mt.filter_entries(filter_condition_ab)
#%%
mt = hl.variant_qc(mt)
#%%
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.01)
#%%
eigenvalues, pcs, _ = hl.hwe_normalized_pca(mt.GT)
#%%
mt = mt.annotate_cols(scores = pcs[mt.s].scores)
#%%
gwas1 = hl.linear_regression_rows(
    y=mt.pheno.CD8,
    x=mt.GT.n_alt_alleles(),
    covariates=[1.0, mt.scores[0], mt.scores[1], mt.scores[2]])

gwas2 = hl.linear_regression_rows(
    y=mt.pheno.MCH,
    x=mt.GT.n_alt_alleles(),
    covariates=[1.0, mt.scores[0], mt.scores[1], mt.scores[2]])
# %%
p1 = hl.plot.manhattan(gwas1.p_value)
show(p1)
p2 = hl.plot.manhattan(gwas2.p_value)
show(p2)

#%%
gwas1.export(f'{prefix}1_results.tsv')
gwas2.export(f'{prefix}2_results.tsv')
# %%
timer2 = time.perf_counter()
timer_total = round(timer2 - timer, 2)
print(f'Runtime: {timer_total}')
# %%
