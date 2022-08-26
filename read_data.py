import numpy as np


def get_gene_expr(infile, cell_line):

    # read in normalized gene expression levels
    gene_expr_data = np.genfromtxt(infile, dtype=None, delimiter=',', names=True, encoding='UTF-8-sig')
    # print(gene_expr_data.dtype.names)
    # print(gene_expr_data['Cell_line'])
    idx = np.where(gene_expr_data['Cell_line'] == cell_line)[0][0]
    ratio = {}  # { gene : ratio }
    for gene in gene_expr_data.dtype.names[1:]:
        val = gene_expr_data[gene][idx]
        # values are in fold-change, so for negative values convert to -1/val
        ratio[gene] = val if val >= 0.0 else -1/val
    return ratio

def get_t_div(infile):

    # read in doubling times
    doubling_times_data = np.genfromtxt(infile, dtype=None, delimiter=',', names=True, encoding='UTF-8-sig')
    # print(doubling_times_data.dtype.names)

    cell_lines = np.unique(doubling_times_data['cell_line'])

    # extract doubling times and calculate division times (hr)
    t_div = {}  # { cell_line: doubling_time }
    for cell_line in cell_lines:
        t_div[cell_line] = 0
        n_samples = 0
        for d in doubling_times_data:
            if d['cell_line'] == cell_line:
                t_div[cell_line] += d['doubling_time_hr']
                n_samples += 1
        t_div[cell_line] /= n_samples  # doubling time
        t_div[cell_line] /= np.log(2)  # division time

    return t_div
