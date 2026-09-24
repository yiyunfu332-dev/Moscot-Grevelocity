from pathlib import Path
import os,json,csv,shutil,hashlib,base64,subprocess,html,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/moscot-professor-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
import numpy as np
import pymupdf as fitz
import argparse
parser=argparse.ArgumentParser(description='Build the MOSCOT professor package from preserved local analysis outputs. Requires matplotlib, numpy, Pillow and PyMuPDF.')
parser.add_argument('--source',type=Path,default=Path.home()/'dynamo')
parser.add_argument('--repository',type=Path,default=Path(__file__).resolve().parents[1])
args=parser.parse_args()
ROOT=args.repository.resolve();SRC=args.source.expanduser().resolve();PACK=ROOT/'professor_package';PACK.mkdir(exist_ok=True)
( PACK/'executed_notebooks').mkdir(exist_ok=True);(PACK/'overview').mkdir(exist_ok=True);(PACK/'previews').mkdir(exist_ok=True)
entries=[]; notebooks=[]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(p,target,status,context,kind):
 target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
 assert digest(p)==digest(target)
 row=dict(source=str(p.relative_to(SRC)),path=str(target.relative_to(PACK)),status=status,context=context,kind=kind,bytes=target.stat().st_size,sha256=digest(target));entries.append(row);return row
specs=[('zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb','CURRENT SAVED RUN','41 executed nonempty cells; no saved errors. Original paths and outputs preserved; not a new execution.'),('hematopoiesis_GraphVelo_moscot.ipynb','HISTORICAL / SAVED ERROR','Earlier combined model; saved categorical-median error. Superseded by native M0/M1.'),('hematopoiesis_graphvelo_2.ipynb','HISTORICAL / PARTIAL RUN','Earlier GraphVelo stage; some nonempty cells have no execution count.'),('hematopoiesis_graphvelo.ipynb','HISTORICAL / SAVED ERROR','Earlier annotation approach; saved SUBTYPE_RULES NameError.'),('zebrafish_human_cross_species_validation.before_expanded_panel.ipynb','SUPERSEDED HUMAN VALIDATION','Executed older marker-panel analysis, not the current expanded-panel result.')]
for name,status,context in specs:
 p=SRC/name;n=json.loads(p.read_text());row=add(p,PACK/'executed_notebooks'/name,status,context,'notebook')
 cs=[c for c in n['cells'] if c['cell_type']=='code' and ''.join(c.get('source',[])).strip()]
 errors=[{'cell':i+1,'name':o.get('ename'),'message':o.get('evalue')} for i,c in enumerate(n['cells']) for o in c.get('outputs',[]) if o.get('output_type')=='error']
 notebooks.append({**row,'nonempty_code_cells':len(cs),'executed_cells':sum(c.get('execution_count') is not None for c in cs),'errors':errors})
 heading=''
 for i,c in enumerate(n['cells']):
  if c['cell_type']=='markdown':
   headings=[l.lstrip('# ').strip() for l in ''.join(c.get('source',[])).splitlines() if l.startswith('#')]
   if headings:heading=headings[-1]
  for j,o in enumerate(c.get('outputs',[])):
   data=o.get('data',{})
   if 'image/png' in data:
    b=data['image/png'];b=''.join(b) if isinstance(b,list) else b
    target=PACK/'notebook_figures'/p.stem/f'cell_{i+1:03d}_output_{j+1:02d}.png';target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(base64.b64decode(b))
    entries.append(dict(source=f'{name}#cell-{i+1}-output-{j+1}',path=str(target.relative_to(PACK)),status=status,context=f'{heading}. {context}',kind='notebook_image',bytes=target.stat().st_size,sha256=digest(target)))
folders=[('blood_growth_graphvelo_moscot_v2/results/figures','CURRENT / READ WITH CORRECTED AUDITS','Native M0/M1; cell-level marker plots are exploratory. Use corrected directional evidence and support gates for interpretation.'),('cross_species_validation/figures','CURRENT EXPRESSION-PROGRAM VALIDATION','Human comparisons test expression programs, not fate, velocity or absolute growth.'),('blood_graphvelo_official/results/figures','HISTORICAL / SUPERSEDED','Earlier combined transport/velocity and trajectory analyses. These do not define the current native M0/M1 conclusions.'),('archive','ARCHIVED / SUPERSEDED','Earlier panels, M2 reweighting or cell-level marker inference; preserved to document unsuccessful/superseded approaches.'),('cluster23_pilot_trajectories','EXPLORATORY CLUSTER-23 PILOT','Cluster is not a lineage identity. Mixed blood/non-blood gene trajectories are exploratory and do not validate blood fate.'),('figs_graphvelo','EXPLORATORY WHOLE-ATLAS GRAPHVELO','Whole-atlas velocity visualizations; not independent validation of the final blood cohort.')]
for folder,status,context in folders:
 for p in sorted((SRC/folder).rglob('*')):
  if p.is_file() and p.suffix.lower() in {'.png','.pdf','.svg','.jpg','.jpeg'}:add(p,PACK/'original_figures'/p.relative_to(SRC),status,context,'standalone_figure')
(PACK/'notebook_execution_manifest.json').write_text(json.dumps(notebooks,indent=2)+'\n')
# Newly drawn summaries use saved tables, with no new fitting or invented values.
def rows(name):return list(csv.DictReader((ROOT/'results_snapshot/2026-09-15/moscot'/name).open()))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
overviews=[]
def save(fig,name,title,note):
 fig.suptitle(title,fontsize=17,y=.97);fig.text(.06,.035,textwrap.fill(note,135),fontsize=9,va='bottom');fig.tight_layout(rect=(.035,.13,.97,.92));p=PACK/'overview'/f'{name}.png';fig.savefig(p,dpi=160);fig.savefig(p.with_suffix('.pdf'));plt.close(fig)
 overviews.append(dict(path=str(p.relative_to(PACK)),status='CURRENT SUMMARY / SAVED TABLES',context=note,title=title))
r=rows('transport_time_support.csv');fig,ax=plt.subplots(figsize=(11.7,8.3));x=np.arange(len(r));vals=[int(t['total_cells']) for t in r];ax.bar(x,vals,color=['#bf5548' if v<50 else '#237f9e' for v in vals]);ax.set_xticks(x,[str(int(float(t['hpf']))) for t in r]);ax.set_xlabel('Developmental time (hpf)');ax.set_ylabel('Modeled cells');ax.axhline(50,color='gray',ls='--',label='50-cell reporting screen');ax.legend()
for i,v in enumerate(vals):ax.text(i,v+20,str(v),ha='center')
save(fig,'01_cohort_support','1. Cohort and sampling support','5,426 cells; 36 different fish; four fish per stage. The 12 hpf stage has 19 cells and cannot support type-level 12–14 hpf claims. Fish at different stages are different animals.')
r=rows('growth_prior_effect_summary.csv');labels=[f"{float(t['source_hpf']):g}–{float(t['target_hpf']):g}" for t in r];x=np.arange(len(r));fig,axes=plt.subplots(2,1,figsize=(11.7,8.3));axes[0].bar(x,[float(t['aggregate_source_mass_tvd_M0_M1']) for t in r],color='#237f9e');axes[0].axhline(.05,color='#bf5548',ls='--');axes[0].set_ylabel('Aggregate source-mass TVD')
y=np.array([float(t['mean_fish_source_mass_tvd_M0_M1']) for t in r]);lo=np.array([float(t['fish_bootstrap_source_tvd_ci_low']) for t in r]);hi=np.array([float(t['fish_bootstrap_source_tvd_ci_high']) for t in r]);axes[1].errorbar(x,y,yerr=[y-lo,hi-y],fmt='o',capsize=4,color='#237f9e');axes[1].axhline(.05,color='#bf5548',ls='--');axes[1].set_ylabel('Mean within-fish source TVD')
for a in axes:a.set_xticks(x,labels)
axes[1].set_xlabel('Interval (hpf)');save(fig,'02_growth_effect','2. Growth priors change late transport most','Top: aggregate source-mass effect. Bottom: mean within-fish normalized effect with descriptive 95% fish-bootstrap intervals; these are different estimands. No bootstrap refitting. The 0.05 line is an operational reporting rule, not biological significance.')
r=rows('moscot_graphvelo_directional_evidence.csv');fig,ax=plt.subplots(figsize=(11.7,8.3));
for i,t in enumerate(r):
 y=float(t['fish_weighted_mean_cosine']);lo=float(t['fish_bootstrap_ci_low']);hi=float(t['fish_bootstrap_ci_high']);color='#23865c' if t['directional_evidence']=='positive_concordance' else '#777777';ax.errorbar(i,y,yerr=[[y-lo],[hi-y]],fmt='o',capsize=5,color=color)
ax.axhline(0,color='black',lw=1);ax.set_xticks(range(len(r)),labels);ax.set_ylabel('Fish-weighted cosine and 95% bootstrap interval');ax.set_xlabel('Interval (hpf)');save(fig,'03_directional_agreement','3. Directional agreement is limited to two intervals','Green: positive absolute direction, positive fish-bootstrap lower bound, and within-fish permutation FDR <0.05 (16–19 and 48–72 hpf). Gray: criteria not met. Late negative cosine is not positive concordance even when better than the shuffled null. Shared expression/PCA means this is complementary computational evidence.')
r=rows('native_moscot_parameter_sensitivity_summary.csv');fig,ax=plt.subplots(figsize=(11.7,8.3));ax.barh(range(len(r)),[float(t['max_type_joint_mass_l1'])/2 for t in r],color='#237f9e');ax.set_yticks(range(len(r)),[t['run'] for t in r]);ax.invert_yaxis();ax.set_xlabel('Maximum type-joint TVD versus primary fit across intervals');save(fig,'04_parameter_sensitivity','4. Solver choices materially affect the coupling','Plotted TVD is half the saved normalized type-joint L1 distance. Eleven one-factor configurations; all fits converged. Numerical convergence does not establish biological accuracy. Balanced tau=1 is a stress test. Source: native_moscot_parameter_sensitivity_summary.csv.')
r=rows('graphvelo_gene_qc_sensitivity_summary.csv');fig,ax=plt.subplots(figsize=(11.7,8.3));x=np.arange(len(r));ax.plot(x,[float(t['cell_cosine_vs_primary_median']) for t in r],'o-',label='Median');ax.plot(x,[float(t['cell_cosine_vs_primary_q05']) for t in r],'s--',label='5th percentile');ax.axhline(0,color='gray');ax.set_xticks(x,[t['gamma_r2_threshold'] for t in r]);ax.set_xlabel('Velocity gene gamma-R2 threshold');ax.set_ylabel('Cell cosine relative to saved primary GraphVelo');ax.legend();save(fig,'05_velocity_qc','5. Stringent velocity-gene filtering destabilizes some cells','Thresholds 0.05 and 0.10 largely preserve primary directions. At 0.20, the lower tail becomes negative. Similarity to the primary fit measures stability, not correctness. Source: graphvelo_gene_qc_sensitivity_summary.csv.')
fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off');text='''WHAT THE CURRENT WORK SUPPORTS\n\n• A reproducible saved native M0/M1 run, plus fish-aware audits.\n• A larger growth-prior source-mass effect at 120–240 hpf.\n• Positive complementary direction evidence at two intervals.\n• Modest human expression-program concordance, with negative specificity tests.\n\nWHAT REMAINS UNRESOLVED\n\n• Actual division/death rates and true parent–descendant relationships.\n• Negative absolute late-interval velocity agreement.\n• Strong epsilon/tau sensitivity and sparse early sampling.\n• Whether growth priors improve independent predictive performance.\n\nHISTORY IS INCLUDED, NOT HIDDEN\n\n• Earlier combined/reweighted models and smaller-marker-panel results.\n• Cluster-23 exploratory trajectories, including mixed lineage markers.\n• Partial notebooks and saved errors, explicitly labeled as historical.''';ax.text(.02,.97,text,va='top',fontsize=13,linespacing=1.65);save(fig,'06_interpretation','6. What succeeded, what did not, and what is still unknown','Historical plots are a record of work, not additional validation of the current model. The separate non-OT project and unapproved integrated proposal are outside this sharing package.')
# Render all PDF pages for a browsable gallery and illustrated PDF appendix.
gallery=list(overviews);seen=set();pdf_renderer=shutil.which('pdftoppm')
for k,e in enumerate(entries):
 p=PACK/e['path']
 if e['kind']=='notebook':continue
 if e['sha256'] in seen:continue
 seen.add(e['sha256'])
 if p.suffix.lower()=='.svg':
  if p.with_suffix('.png').exists():continue
  continue # original SVG remains linked in the full manifest/gallery index
 if p.suffix.lower()=='.pdf':
  doc=fitz.open(p)
  for page_number,page in enumerate(doc,1):
   preview=PACK/'previews'/f'figure_{k:03d}-{page_number:03d}.png'
   scale=1400/max(page.rect.width,page.rect.height)
   page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False).save(preview)
   gallery.append({**e,'path':str(preview.relative_to(PACK)),'title':p.stem.replace('_',' ')+' / page '+str(page_number),'original':e['path']})
  doc.close()
 elif p.suffix.lower() in {'.png','.jpg','.jpeg'}:gallery.append({**e,'title':p.stem.replace('_',' '),'original':e['path']})
print('Copied',len(entries),'artifacts; rendering',len(gallery),'gallery pages',flush=True)
with PdfPages(PACK/'MOSCOT_project_figures_for_professor.pdf') as pdf:
 for idx,e in enumerate(gallery):
  fig=plt.figure(figsize=(11.7,8.3));fig.text(.035,.965,f"{idx+1}. {e.get('title','Figure')[:100]}",fontsize=13,va='top');fig.text(.035,.925,e['status'],fontsize=10,color='#a54330' if any(w in e['status'] for w in ['HISTOR','ARCHIV','SUPERSEDED','EXPLOR']) else '#176b65');ax=fig.add_axes([.035,.17,.93,.72]);ax.imshow(Image.open(PACK/e['path']));ax.axis('off');fig.text(.035,.13,textwrap.fill(e['context'],150),fontsize=8,va='top');fig.text(.035,.035,textwrap.shorten(e.get('original',e['path']),width=170,placeholder='...'),fontsize=7);pdf.savefig(fig,dpi=130);plt.close(fig)
with (PACK/'artifact_manifest.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(entries[0]),lineterminator='\n');w.writeheader();w.writerows(entries)
(PACK/'gallery_manifest.json').write_text(json.dumps(gallery,indent=2)+'\n')
parts=['<!doctype html><html lang="en"><meta charset="utf-8"><title>MOSCOT project figure gallery</title><style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:20px;color:#233}article{border-top:1px solid #ccc;padding:25px 0}img{max-width:100%;height:auto}.status{font-weight:bold;color:#9c4836}</style><h1>MOSCOT project: current results and complete figure history</h1><p>Saved analyses, including negative, partial and superseded work. No scientific fits were rerun for this package. Non-OT work is excluded.</p><p><a href="MOSCOT_project_figures_for_professor.pdf">Download the complete figure report</a> · <a href="README.md">Reading guide</a> · <a href="artifact_manifest.csv">All original files and provenance</a></p>']
for i,e in enumerate(gallery):parts.append(f'<article><h2>{i+1}. {html.escape(e.get("title","Figure"))}</h2><p class="status">{html.escape(e["status"])}</p><p>{html.escape(e["context"])}</p><img loading="lazy" src="{html.escape(e["path"],quote=True)}" alt="{html.escape(e.get("title","Figure"),quote=True)}"><p><a href="{html.escape(e.get("original",e["path"]),quote=True)}">Original figure</a></p></article>')
parts.append('</html>');(PACK/'gallery.html').write_text('\n'.join(parts))
# professor_package/README.md is maintained editorially; rebuilding figures must preserve it.

print(json.dumps({'source_artifacts':len(entries),'notebooks':len(notebooks),'figure_report_pages':len(gallery),'report_bytes':(PACK/'MOSCOT_project_figures_for_professor.pdf').stat().st_size}),flush=True)

# Companion report, index and package integrity checks.
import pymupdf
r=ROOT
p=PACK
# Six-page concise companion for an initial professor discussion.
doc=pymupdf.open()
for f in sorted((p/'overview').glob('*.pdf')):
 with pymupdf.open(f) as part:doc.insert_pdf(part)
doc.save(p/'MOSCOT_summary_for_professor.pdf');doc.close()
gallery=json.loads((p/'gallery_manifest.json').read_text())
md='# Figure index: current results and project history\n\nAll figures are included regardless of whether they support the current model. Status labels identify exploratory and superseded work. [Six-page summary](MOSCOT_summary_for_professor.pdf) · [80-page full report](MOSCOT_project_figures_for_professor.pdf) · [Executed notebooks](README.md#saved-notebook-history).\n\n'
for i,e in enumerate(gallery[:6],1):md+=f"## {i}. {e['title']}\n\n{e['context']}\n\n![{e['title']}]({e['path']})\n\n"
md+='## Complete displayed figure index\n\n| Page | Figure | Status | Original |\n|---:|---|---|---|\n'
for i,e in enumerate(gallery,1):md+=f"| {i} | [{e['title']}]({e['path']}) | {e['status']} | [File]({e.get('original',e['path'])}) |\n"
(p/'FIGURE_INDEX.md').write_text(md)

# Verify source copies, image decoding, all PDF-page coverage, and notebook state.
entries=list(csv.DictReader((p/'artifact_manifest.csv').open()))
for e in entries:
 f=p/e['path'];assert f.stat().st_size==int(e['bytes']);assert hashlib.sha256(f.read_bytes()).hexdigest()==e['sha256']
for e in gallery:
 with Image.open(p/e['path']) as im:im.verify()
unique_pdfs={e['sha256']:e for e in entries if e['path'].endswith('.pdf')}
for e in unique_pdfs.values():
 with pymupdf.open(p/e['path']) as d:expected=len(d)
 shown=sum(x.get('original')==e['path'] for x in gallery)
 if shown==0:
  duplicate_paths={x['path'] for x in entries if x['sha256']==e['sha256']};shown=sum(x.get('original') in duplicate_paths for x in gallery)
 assert shown==expected,(e['path'],expected,shown)
with pymupdf.open(p/'MOSCOT_project_figures_for_professor.pdf') as d:assert len(d)==len(gallery)==80
with pymupdf.open(p/'MOSCOT_summary_for_professor.pdf') as d:assert len(d)==6
n=json.loads((p/'executed_notebooks/zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb').read_text());cells=[c for c in n['cells'] if c['cell_type']=='code' and ''.join(c.get('source',[])).strip()];assert len(cells)==41 and all(c.get('execution_count') is not None for c in cells);assert not any(o.get('output_type')=='error' for c in cells for o in c.get('outputs',[]))
(p/'VERIFICATION.md').write_text(f'''# Sharing package verification

- {len(entries)} source/notebook-image artifacts: sizes and SHA-256 verified.
- Original primary notebook: exact source copy; 41/41 executed nonempty code cells; zero saved errors.
- Five notebooks: saved execution/error evidence preserved and recorded, including historical failures.
- {len(gallery)} gallery images decoded successfully.
- Every page of every unique source PDF is represented in the gallery/report; SVG counterparts and byte-identical duplicates remain available as original files.
- Full PDF: 80 pages. Concise PDF: six pages.
- Overview charts visually inspected for layout/readability and drawn directly from saved tables; no scientific fits rerun.
- No separate non-OT analysis or integrated proposal included.
''')
print('Verified',len(entries),'artifacts;',len(gallery),'gallery pages; primary notebook exact and executed')
