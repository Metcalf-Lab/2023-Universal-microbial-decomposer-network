import glob
import pandas as pd

emp_md = pd.read_csv('emp_qiime_mapping_release1.tsv', sep='\t', dtype=str).set_index('#SampleID')
emp_ids = set(emp_md.index)

asv_metadata = {
    f.split('/')[1]: pd.read_csv(f, sep='\t', dtype=str).set_index('#SampleID')
    for f in glob.glob('results/*/samples.metadata.tsv')
}

md_results = []
feat_md = []
for asv, md in asv_metadata.items():
    n_emp = len(emp_ids & set(md.index))

    agp = {i for i in md.index if i.startswith('10317.')}
    agp_md = md.loc[list(agp)]
    agp_md = agp_md[agp_md['env_package'].isin(['human-gut', 'human-oral', 'human-skin'])]
    agp_ids = set(agp_md.index)
    n_agp = len(agp_ids & set(md.index))

    feat_md.append((asv, n_emp, n_agp))

    md['is_emp'] = False
    md['is_agp'] = False
    md.loc[[i for i in md.index if i in agp_ids], 'is_agp'] = True
    md.loc[[i for i in md.index if i in emp_ids], 'is_emp'] = True

    md_results.append(md.loc[list((emp_ids | agp_ids) & set(md.index))])

md = pd.concat(md_results).reset_index()
md = md[~md['#SampleID'].duplicated()]
md.set_index('#SampleID', inplace=True)

md.loc[[i for i in md.index if i.startswith('2229.')], 'empo_1'] = 'Host-associated'
md.loc[[i for i in md.index if i.startswith('2229.')], 'empo_2'] = 'Plant'
md.loc[[i for i in md.index if i.startswith('2229.')], 'empo_3'] = 'Plant surface'

empo_1_fix = {'host-associated': 'Host-associated'}
md['empo_1'] = md['empo_1'].apply(lambda x: empo_1_fix.get(x, x))
empo_2_fix = {'animal': 'Animal'}
md['empo_2'] = md['empo_2'].apply(lambda x: empo_2_fix.get(x, x))
empo_3_fix = {'animal surface': 'Animal surface',
              'animal secretion': 'Animal secretion',
              'animal distal gut': 'Animal distal gut'}
md['empo_3'] = md['empo_3'].apply(lambda x: empo_3_fix.get(x, x))

md.to_csv('sample-metadata.tsv', sep='\t', index=True, header=True)

feat_md = pd.DataFrame(feat_md, columns=['feature-id', 'number_emp', 'number_agp'])
feat_md.to_csv('feature-metadata.tsv', sep='\t', index=False, header=True)
