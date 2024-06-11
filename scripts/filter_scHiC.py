#! /usr/bin/env python
 
import os
from argparse import ArgumentParser
 
from matplotlib import pyplot as plt
import numpy as np
 
 
def get_cells(fname):
    cells = {}
    fh = open(fname, encoding="utf-8")
    header = ''
    for line in fh:
        if line.startswith("#"):
            header += line
        else:
            break
    fh.seek(len(header))
 
    for line in fh:
        _, c1, b1, _, l1, _, _, c2, b2, _, l2, _, _, ct = line.split()
        if c1 == c2:
            dist = int(abs((int(b2) + int(l2) / 2) - (int(b1) + int(l1) / 2)))
        else:
            dist = -1
        cells.setdefault(ct, []).append(dist)
    return cells

def get_cells_cluster(fname, reso, too_close, cluster_contacts):
    cells = {}
    fh = open(fname, encoding="utf-8")
    header = ''
    for line in fh:
        if line.startswith("#"):
            header += line
        else:
            break
    fh.seek(len(header))
 
    for line in fh:
        _, c1, b1, _, l1, _, _, c2, b2, _, l2, _, _, ct = line.split()
        
        b1, l1, b2, l2 = map(int, [b1, l1, b2, l2])
        p1 = b1 + l1 / 2
        #print(p1)
        #print(p2)
        p2 = b2 + l2 / 2
        if c1 == c2:
            if abs(int(p1) - int(p2)) < too_close:
                continue
            p1 = int(int(p1) // reso)
            p2 = int(int(p2) // reso)
        else:
            p1 = int(int(p1) // reso*10)
            p2 = int(int(p2) // reso*10)
            
        feature = f"{c1}:{p1}_{c2}:{p2}"
        if feature in cluster_contacts:
            if c1 == c2:
                dist = int(abs((int(b2) + int(l2) / 2) - (int(b1) + int(l1) / 2)))
            else:
                dist = -1
            cells.setdefault(ct, []).append(dist)
    return cells
 
 
def cis_ratio(cells, ratio_cut=0.25, total_cut=2000, cis_far=10_000):
    cell_ratio = []
    cell_total = []
    good_cells = {}
    min_total = min(100, total_cut / 2)
    for cell, counts in cells.items():
        total = len(counts)
        if total < min_total:
            #print(cell, total)
            continue
        trans  = sum(v == -1 for v in counts)
        cis10k = sum(v > cis_far for v in counts)
        cis0k  = len(counts) - trans
        try:
            ratio10k = cis10k / (cis10k + trans)
            cell_ratio.append(ratio10k)
        except ZeroDivisionError:
            continue
        cell_total.append(total)
        if ratio10k > ratio_cut and total > total_cut:
            good_cells[cell] = {
                'cis': cis0k, 'cis 10kb': cis10k, 'trans': trans, 
                'cis ratio': cis0k / (cis0k + trans),
                'cis ratio 10kb': cis10k / (cis10k + trans)}
    return cell_ratio, cell_total, ratio_cut, total_cut, good_cells

def cis_ratio_cluster(cells, cluster_cells, cis_far=10_000):

    good_cells = {}
    #min_total = min(100, total_cut / 2)
    for cell, counts in cells.items():
        if cell not in cluster_cells:
            continue
        total = len(counts)
        """         if total < min_total:
            continue """
        trans  = sum(v == -1 for v in counts)
        cis10k = sum(v > cis_far for v in counts)
        cis0k  = len(counts) - trans

        good_cells[cell] = {
                'cis': cis0k, 'cis 10kb': cis10k, 'trans': trans, 
                'cis ratio': cis0k / (cis0k + trans),
                'cis ratio 10kb': cis10k / (cis10k + trans)}
    return  good_cells
 
 
def do_plot(cell_ratio, cell_total, ratio_cut=0.25, total_cut=2000, outdir=None):
    good_total = sum(r > ratio_cut and t > total_cut 
                     for r, t in zip(cell_ratio, cell_total))
    _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9, 9))
 
    ax1.set_position((0.1, 0.1, 0.7, 0.7))
    ax1.plot(cell_ratio, cell_total, 'k.', ms=7, alpha=0.4, 
             label=f'All cells ({len(cell_total):,.0f})')
    ax1.set_ylabel('Total Counts')
    ax1.set_xlabel('Cis (>10kb) ratio')
    ax1.set_yscale('log')
    ax1.grid()
 
    x0, x1 = ax1.get_xlim()
    # x0 = 0
    y0, y1 = ax1.get_ylim()
    y0 = 100
 
    ax1.plot([ratio_cut, ratio_cut, 1],
            [y1, total_cut, total_cut], color='tab:red', 
            label=f'Good cells ({good_total:,.0f})')
    ax1.set_xlim(x0, x1)
    ax1.set_ylim(y0, y1)
 
    xdiff = x1 - x0
    ydiff = y1 - y0
 
    ax1.text(x1 - xdiff * 0.01, total_cut + ydiff * 0.0002, f"Count > {total_cut}", 
            color='tab:red', ha='right')
    ax1.text(ratio_cut + xdiff * 0.01, y1 - ydiff * 0.1, 
            f"Ratio > {ratio_cut}", color='tab:red', rotation=90,
            ha = 'left', va='top')
    ax1.legend(bbox_to_anchor=(1,1), loc="lower left")
 
    ax2.set_position((0.8, 0.1, 0.1, 0.7))
    ax2.hist(np.log10(cell_total), bins=50, edgecolor='k', 
             facecolor='tab:grey', alpha=0.4, orientation='horizontal')
    ax2.axhline(np.log10(total_cut), color='tab:red', ls='-')
    ax2.set_xlabel('Count')
    ax2.set_ylim(np.log10(y0), np.log10(y1))
    ax2.set_yticklabels([])
    ax2.grid(axis='y')
 
    ax3.set_position((0.1, 0.8, 0.7, 0.1))
    ax3.hist(cell_ratio, bins=50, edgecolor='k',
             facecolor='tab:grey', alpha=0.4)
    ax3.axvline(ratio_cut, color='tab:red', ls='-')
    ax3.set_xlim(x0, x1)
    ax3.set_xticklabels([])
    ax3.grid(axis='x')
    ax3.set_ylabel('Count')
 
    if outdir is not None:
        plt.savefig(os.path.join(outdir, "sc_cis-ratio_distribution.svg"),format='svg')
 
 
def main():
    opts = get_options()
 
    all_ditags = opts.all_ditags
    ratio_cut  = opts.cis_ratio
    total_cut  = opts.count_cut
    outdir     = opts.outdir
    
    cis_far = 10_000
    
    os.system(f"mkdir -p {outdir}")
 
    cells = get_cells(all_ditags)
    cell_ratio, cell_total, ratio_cut, total_cut, good_cells = cis_ratio(
        cells, ratio_cut=ratio_cut, total_cut=total_cut)
 
    # save data
    out = open(os.path.join(outdir, 'cell_stats.tsv'), 'w', encoding='utf-8')
    out.write(f"## Cis-far defined as cis interactions with a distance above {cis_far}\n")
    out.write(f"## Selected cell with a cis-far ratio above {ratio_cut}\n")
    out.write(f"## and a total di-tag count above {total_cut}\n")
    out.write(f"# cell-tag\tcis count\tcis-far count\ttrans count\tcis ratio\tcis-far ratio\n")
    for cell, vals in good_cells.items():
        out.write(f"{cell}\t{vals['cis']}\t{vals['cis 10kb']}\t{vals['trans']}\t{vals['cis ratio']}\t{vals['cis ratio 10kb']}\n")
    out.close()
 
    do_plot(cell_ratio, cell_total, ratio_cut=ratio_cut, total_cut=total_cut, outdir=outdir)
 
 
def get_options():
    parser = ArgumentParser()
 
    parser.add_argument('-i', dest='all_ditags', metavar='PATH', required=True,
                        default=False, help='''Input TSV file containing 
                        di-tags from all cells (last column should be the cell
                        tag).''')
    parser.add_argument('-o', dest='outdir', metavar='PATH', required=True,
                        default=True, help='Output directory.')
    parser.add_argument('-c', dest='count_cut', metavar='INT', default=2000, type=int,
                        help='Minimum number of di-tags per cell.')
    parser.add_argument('-r', dest='cis_ratio', metavar='FLOAT', default=0.25, type=float,
                        help='Minimum cis-ratio (excluding cis interactions bellow 10kb) per cell.')
    opts = parser.parse_args()
    return opts
 
if __name__ == "__main__":
    exit(main())
