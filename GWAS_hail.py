#!/usr/bin/python3.6

import os
memory = '256g'
pyspark_submit_args = ' --driver-memory ' + memory + ' pyspark-shell'
os.environ["PYSPARK_SUBMIT_ARGS"] = pyspark_submit_args

import hail as hl
hl.init(min_block_size=128, master='local[16]')


# In[43]:


from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()


# In[89]:


hl.import_vcf('ALL.chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.split.AF_0.05_v2.vcf.bgz').write('ALL.chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.split.AF_0.05_v2.mt', overwrite = True)
mt = hl.read_matrix_table('ALL.chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.split.AF_0.05_v2.mt')


# In[90]:


#mt = mt.filter_rows(mt.locus.contig == '22')
table = (hl.import_table('integrated_call_samples_v3.20130502.ALL.pheno', impute=True)
         .key_by('sample'))
mt = mt.annotate_cols(pheno = table[mt.s])


# In[91]:


snp_counts = mt.aggregate_rows(hl.agg.counter(hl.Struct(ref=mt.alleles[0], alt=mt.alleles[1])))
pprint(snp_counts)
from collections import Counter
counts = Counter(snp_counts)
counts.most_common()


# In[92]:


mt.col.describe()
mt = hl.sample_qc(mt)
mt.col.describe()
mt = hl.variant_qc(mt)
mt.row.describe()
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.01)
mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6)
print('Samples: %d  Variants: %d' % (mt.count_cols(), mt.count_rows()))


# In[74]:


mt.describe()


# In[93]:


gwas = hl.linear_regression_rows(y=mt.pheno.pheno,
                                 x=mt.GT.n_alt_alleles(),
                                 covariates=[1.0])
gwas.row.describe()
gwas.export('ALL.chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.split.AF_0.05_v2_noPCAs.gwas.csv')
p = hl.plot.manhattan(gwas.p_value)
show(p)
p = hl.plot.qq(gwas.p_value)
show(p)
eigenvalues, pcs, _ = hl.hwe_normalized_pca(mt.GT, k = 5)
pprint(eigenvalues)
pcs.show(3, width=100)
mt = mt.annotate_cols(scores = pcs[mt.s].scores)
p = hl.plot.scatter(mt.scores[0],
                    mt.scores[1],
                    title='PCA', xlabel='PC1', ylabel='PC2')
show(p)


# In[94]:


gwas = hl.linear_regression_rows(
    y=mt.pheno.pheno,
    x=mt.GT.n_alt_alleles(),
    covariates=[1.0, mt.scores[0], mt.scores[1], mt.scores[2]])
p = hl.plot.qq(gwas.p_value)
show(p)
p = hl.plot.manhattan(gwas.p_value)
show(p)
hl.plot.pdf(gwas.p_value)


# In[95]:


gwas.order_by('p_value')
gwas.show(5)


# In[99]:


gwas.row.describe()


# In[96]:


gwas.export('ALL.chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.split.AF_0.05_v2.gwas.csv')

