"""Generate legacy LaTeX tables for the RESCUE paper."""
import pandas as pd
import json
import numpy as np
import os

ROOT = os.path.join(os.path.dirname(__file__), '..', '..', '..')
OUT_DIR = os.path.join(ROOT, 'paper', 'v4', 'tables', 'generated_legacy')
os.makedirs(OUT_DIR, exist_ok=True)

def fmt_ci(mean, lo, hi, ndigits=3):
    """Format mean with CI as LaTeX."""
    fmt = f'{{:.{ndigits}f}}'
    return f'${fmt.format(mean)}$ ({fmt.format(lo)}, {fmt.format(hi)})'

def bootstrap_ci_dict(csv_path=None):
    if csv_path is None:
        csv_path = os.path.join(ROOT, 'paper', 'v2', 'figures', 'out', 'table_exp1_benchmark.csv')
    """Build lookup from V2 bootstrap CI CSV: {(tool_lower, species, regime): {f1, lo, hi}}."""
    df = pd.read_csv(csv_path)
    ci = {}
    # Map V2 tool names to V3 TSV keys
    tool_map = {
        'miniprot': 'miniprot',
        'Spaln3': 'spaln3',
        'GTH': 'gth',
        'Scipio': 'scipio',
        'genBlastG': 'genblastg',
        'mini+miniSplice': 'minisplice_miniprot',
    }
    for _, row in df.iterrows():
        v2_tool = row['tool']
        tsv_key = tool_map.get(v2_tool, v2_tool)
        sp = row['sp_key']
        # 1:1
        ci[(tsv_key, sp, 'one_to_one')] = {
            'f1': row['o2o_f1'], 'lo': row['o2o_lo'], 'hi': row['o2o_hi']
        }
        # Paralocal primary
        ci[(tsv_key, sp, 'paralog_primary')] = {
            'f1': row['para_prim_f1'], 'lo': row['para_prim_lo'], 'hi': row['para_prim_hi']
        }
        # Paralocal best_overlap
        ci[(tsv_key, sp, 'paralog_bo')] = {
            'f1': row['para_bo_f1'], 'lo': row['para_bo_lo'], 'hi': row['para_bo_hi']
        }
    return ci

# ============================================================
# TABLE 1: Complete End-to-End Benchmark + TPRS
# ============================================================
def gen_table1():
    df = pd.read_csv(os.path.join(ROOT, 'results/exp1/summary.tsv'), sep='\t')
    with open(os.path.join(ROOT, 'results/tprs/phase5/ablation_summary_all_species.json')) as f:
        tprs = json.load(f)
    ci = bootstrap_ci_dict()

    SPECIES = ['mouse', 'dog', 'chicken', 'zebrafish']
    SP_DISPLAY = ['Mouse', 'Dog', 'Chicken', 'Zebrafish']
    TOOLS = [
        ('miniprot', 'miniprot'),
        ('spaln3', 'Spaln3'),
        ('gth', 'GTH'),
        ('scipio', 'Scipio'),
        ('genblastg', 'genBlastG'),
        ('minisplice_miniprot', 'mini+miniSplice'),
    ]

    # Best e2e per species (paralog primary)
    best_paralog_primary = {}
    best_paralog_bo = {}
    for sp in SPECIES:
        bp, bbo = 0, 0
        for tool_key, _ in TOOLS:
            rp = df[(df['species']==sp) & (df['tool']==tool_key) & (df['regime']=='paralog')]
            if len(rp) > 0:
                bp = max(bp, rp['exon_f1_mean'].values[0])
            key = (tool_key, sp, 'paralog_bo')
            if key in ci:
                bbo = max(bbo, ci[key]['f1'])
        best_paralog_primary[sp] = bp
        best_paralog_bo[sp] = bbo

    # TPRS values from ablation E1
    tprs_primary = {}
    tprs_bo = {}
    tprs_p_ci = {}
    tprs_bo_ci = {}
    for sp in SPECIES:
        tprs_primary[sp] = tprs[sp]['primary']['E1']['mean']
        tprs_bo[sp] = tprs[sp]['best_overlap']['E1']['mean']
        tprs_p_ci[sp] = (tprs[sp]['primary']['E1']['ci_lower'],
                         tprs[sp]['primary']['E1']['ci_upper'])
        tprs_bo_ci[sp] = (tprs[sp]['best_overlap']['E1']['ci_lower'],
                          tprs[sp]['best_overlap']['E1']['ci_upper'])

    lines = []
    lines.append(r'\begin{table*}[t]')
    lines.append(r'\centering')
    lines.append(r'\caption{Complete end-to-end benchmark results (Exp~1) and TPRS performance. '
                 r'Mean exon F1 with 95\% paired-bootstrap CI ($10^4$ resamples, seed 4242). '
                 r'Paralog columns report \textbf{primary} / \emph{best-overlap} F1 '
                 r'(evaluated against actually-matched ortholog). '
                 r'$n = 2000$ (1:1), 488--500 (paralog) per species.}')
    lines.append(r'\label{tab:benchmark}')
    lines.append(r'\footnotesize')
    lines.append(r'\setlength{\tabcolsep}{2.5pt}')
    col_spec = 'l' + ''.join(['cc' for _ in SPECIES])
    lines.append(r'\resizebox{\textwidth}{!}{%')
    lines.append(r'\begin{tabular}{' + col_spec + '}')
    lines.append(r'\toprule')
    header1 = ' & '.join([f'\\multicolumn{{2}}{{c}}{{{d}}}' for d in SP_DISPLAY])
    lines.append(f'Tool & {header1} \\\\')
    header2_parts = []
    for sp in SPECIES:
        header2_parts.append('1:1 & Paralog')
    lines.append(f' & {" & ".join(header2_parts)} \\\\')
    lines.append(r'\midrule')

    # Tool rows: 1:1 = mean only, Paralog = primary (CI) / best_overlap (CI)
    for tool_key, tool_label in TOOLS:
        row = [f'\\texttt{{{tool_label}}}']
        for sp in SPECIES:
            key_1 = (tool_key, sp, 'one_to_one')
            if key_1 in ci:
                c = ci[key_1]
                row.append(f'{c["f1"]:.3f}')
            else:
                row.append('--')
            key_pp = (tool_key, sp, 'paralog_primary')
            key_pbo = (tool_key, sp, 'paralog_bo')
            if key_pp in ci and key_pbo in ci:
                pp = ci[key_pp]
                pbo = ci[key_pbo]
                row.append(f'{pp["f1"]:.3f} ({pp["lo"]:.2f}, {pp["hi"]:.2f}) / '
                           f'\\emph{{{pbo["f1"]:.3f} ({pbo["lo"]:.2f}, {pbo["hi"]:.2f})}}')
            else:
                row.append('--')
        lines.append(' & '.join(row) + r' \\')

    lines.append(r'\midrule')
    # TPRS row
    tprs_row = [r'\textbf{TPRS}']
    for sp in SPECIES:
        tprs_row.append('--')
        plo, phi = tprs_p_ci[sp]
        blo, bhi = tprs_bo_ci[sp]
        tprs_row.append(f'\\textbf{{{tprs_primary[sp]:.3f} ({plo:.2f}, {phi:.2f})}} / '
                        f'\\emph{{{tprs_bo[sp]:.3f} ({blo:.2f}, {bhi:.2f})}}')
    lines.append(' & '.join(tprs_row) + r' \\')

    # Delta rows
    delta_p_row = [r'$\Delta_{\mathrm{primary}}$ vs best e2e']
    delta_bo_row = [r'$\Delta_{\mathrm{best\text{-}overlap}}$ vs best e2e']
    for sp in SPECIES:
        delta_p_row.append('')
        delta_p = tprs_primary[sp] - best_paralog_primary[sp]
        delta_p_row.append(f'\\textbf{{{delta_p:+.3f}}}')
        delta_bo_row.append('')
        delta_bo = tprs_bo[sp] - best_paralog_bo[sp]
        delta_bo_row.append(f'\\textbf{{{delta_bo:+.3f}}}')
    lines.append(' & '.join(delta_p_row) + r' \\')
    lines.append(' & '.join(delta_bo_row) + r' \\')
    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular}')
    lines.append(r'}')  # close resizebox
    lines.append(r'\end{table*}')

    with open(f'{OUT_DIR}/table1_benchmark.tex', 'w') as f:
        f.write('\n'.join(lines))
    print('Table 1 written.')

# ============================================================
# TABLE 2: Pipeline Component Ablation (two single-column tables)
# ============================================================
def gen_table2():
    with open(os.path.join(ROOT, 'results/tprs/phase5/ablation_summary_all_species.json')) as f:
        data = json.load(f)

    SPECIES = ['mouse', 'dog', 'chicken', 'zebrafish']
    SP_LABELS = ['Mouse', 'Dog', 'Chicken', 'Zebrafish']

    BO_EXPS = [
        ('E1: Full TPRS', 'E1', False),
        ('E6: --last', 'E6_last', True),
        ('E6: --mmseqs', 'E6_mmseqs', True),
        ('E3: --Router', 'E3', True),
        ('E6: --diamond', 'E6_diamond', True),
        ('E6: --miniprot', 'E6_miniprot', True),
        ('E2: --Sp. prot.', 'E2', True),
        ('E4: --Co-orth.', 'E4', True),
        ('E5: --C-stage\\textsuperscript{\\dag}', 'E5', True),
    ]

    P_EXPS = [
        ('E1: Full TPRS', 'E1', False),
        ('E3: --Router', 'E3', True),
        ('E6: --last', 'E6_last', True),
        ('E6: --miniprot', 'E6_miniprot', True),
        ('E6: --diamond', 'E6_diamond', True),
        ('E6: --mmseqs', 'E6_mmseqs', True),
        ('E2: --Sp. prot.', 'E2', True),
        ('E4: --Co-orth.', 'E4', True),
        ('E5: --C-stage\\textsuperscript{\\dag}', 'E5', True),
    ]

    e1_bo = {sp: data[sp]['best_overlap']['E1']['mean'] for sp in SPECIES}
    e1_p = {sp: data[sp]['primary']['E1']['mean'] for sp in SPECIES}

    def write_one_panel(exps, e1_vals, metric_name, panel_label, out_path):
        lines = []
        lines.append(r'\begin{table}[t]')
        lines.append(r'\centering')
        lines.append(r'\caption{Pipeline component ablation -- ' + panel_label + r' F1. '
                     r'Mean F1. $n = 445$-$489$ per species. '
                     r'\textbf{Bold} = mean $|\Delta| > 0.20$ (critical).}')
        lines.append(r'\label{tab:ablation_' + panel_label.lower().replace(' ','_') + '}')
        lines.append(r'\scriptsize')
        lines.append(r'\setlength{\tabcolsep}{2pt}')
        lines.append(r'\begin{tabular}{l' + 'c' * (len(SPECIES)) + 'c}')
        lines.append(r'\toprule')
        sp_header = ' & '.join(SP_LABELS) + ' & Mean $\\Delta$'
        lines.append(f'Experiment & {sp_header} \\\\')
        lines.append(r'\midrule')

        for label, key, is_ablation in exps:
            row = [label]
            deltas = []
            for sp in SPECIES:
                v = data[sp][metric_name][key]['mean']
                if not is_ablation:
                    row.append(f'{v:.2f}')
                else:
                    d = v - e1_vals[sp]
                    deltas.append(d)
                    row.append(f'{v:.2f}')
            if is_ablation:
                mean_d = np.mean(deltas)
                is_crit = abs(mean_d) > 0.20
                d_str = f'\\textbf{{{mean_d:+.2f}}}' if is_crit else f'{mean_d:+.2f}'
                row.append(d_str)
            else:
                row.append('--')
            lines.append(' & '.join(row) + r' \\')

        lines.append(r'\bottomrule')
        lines.append(r'\end{tabular}')
        lines.append(r'\end{table}')

        with open(out_path, 'w') as f:
            f.write('\n'.join(lines))

    write_one_panel(BO_EXPS, e1_bo, 'best_overlap', 'best-overlap',
                    os.path.join(OUT_DIR, 'table2a_ablation_bo.tex'))
    write_one_panel(P_EXPS, e1_p, 'primary', 'primary',
                    os.path.join(OUT_DIR, 'table2b_ablation_primary.tex'))
    print('Table 2a, 2b written.')

# ============================================================
# TABLE 3: Forced-Locus Alignment Ceiling
# ============================================================
def gen_table3():
    df = pd.read_csv(os.path.join(ROOT, 'results/exp5/cell/per_cell.tsv'), sep='\t')
    df = df[df['padding'] == 50000]

    SPECIES = ['mouse', 'dog', 'chicken', 'zebrafish']
    SP_DISPLAY = ['Mouse', 'Dog', 'Chicken', 'Zebrafish']
    TOOLS = ['miniprot', 'genblastg', 'scipio', 'spaln3', 'gth', 'genewise', 'exonerate']
    TOOL_LABELS = {'miniprot': 'miniprot', 'genblastg': 'genBlastG', 'scipio': 'Scipio',
                   'spaln3': 'Spaln3', 'gth': 'GTH', 'genewise': 'Genewise', 'exonerate': 'Exonerate'}

    lines = []
    lines.append(r'\begin{table*}[t]')
    lines.append(r'\centering')
    lines.append(r'\caption{Forced-locus alignment ceiling (Exp~5, $\pm$50~kb padding). '
                 r'Mean exon F1 with 95\% bootstrap CI. '
                 r'$n = 2000$ (1:1), 488--500 (paralog) per species.}')
    lines.append(r'\label{tab:forced_locus}')
    lines.append(r'\scriptsize')
    lines.append(r'\setlength{\tabcolsep}{3.5pt}')
    col_spec = 'l' + ''.join(['cc' for _ in SPECIES])
    lines.append(r'\begin{tabular}{' + col_spec + '}')
    lines.append(r'\toprule')

    header1 = ' & '.join([f'\\multicolumn{{2}}{{c}}{{{d}}}' for d in SP_DISPLAY])
    lines.append(f'Tool & {header1} \\\\')
    header2_parts = []
    for sp in SPECIES:
        header2_parts.append('1:1 & Paralog')
    lines.append(f' & {" & ".join(header2_parts)} \\\\')
    lines.append(r'\midrule')

    for tool in TOOLS:
        label = TOOL_LABELS.get(tool, tool)
        row = [f'\\texttt{{{label}}}']
        for sp in SPECIES:
            r1 = df[(df['species']==sp) & (df['tool']==tool) & (df['regime']=='one2one')]
            rp = df[(df['species']==sp) & (df['tool']==tool) & (df['regime']=='paralog_broad')]
            v1 = r1['mean_exon_f1'].values[0] if len(r1) > 0 else float('nan')
            vp = rp['mean_exon_f1'].values[0] if len(rp) > 0 else float('nan')
            row.append(f'{v1:.3f}' if not np.isnan(v1) else '--')
            row.append(f'{vp:.3f}' if not np.isnan(vp) else '--')
        lines.append(' & '.join(row) + r' \\')

    lines.append(r'\midrule')
    gap_row = [r'\textbf{Identity Gap} \\ (1:1 $-$ Paralog)']
    for sp in SPECIES:
        r1 = df[(df['species']==sp) & (df['tool']=='miniprot') & (df['regime']=='one2one')]
        rp = df[(df['species']==sp) & (df['tool']=='miniprot') & (df['regime']=='paralog_broad')]
        v1 = r1['mean_exon_f1'].values[0] if len(r1) > 0 else 0
        vp = rp['mean_exon_f1'].values[0] if len(rp) > 0 else 0
        gap = v1 - vp
        gap_row.append(f'\\multicolumn{{2}}{{c}}{{\\textbf{{{gap:.3f}}}}}')
    lines.append(' & '.join(gap_row) + r' \\')
    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular}')
    lines.append(r'\end{table*}')

    with open(f'{OUT_DIR}/table3_forced_locus.tex', 'w') as f:
        f.write('\n'.join(lines))
    print('Table 3 written.')


if __name__ == '__main__':
    gen_table1()
    gen_table2()
    gen_table3()
