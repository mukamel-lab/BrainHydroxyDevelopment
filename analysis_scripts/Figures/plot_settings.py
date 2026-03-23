import pandas as pd
import matplotlib.pyplot as plt

def apply_plot_settings():
   
    # Matplotlib settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'pdf.fonttype': 42,
        'figure.dpi': 150,
        'figure.titlesize': 8,
        'axes.titlesize': 8,
        'axes.labelsize': 8,
        'axes.spines.right': False,
        'axes.spines.top': False,
        'axes.linewidth': 0.5,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 1.5,
        'ytick.major.size': 1.5,
        'legend.fontsize': 6,
        'grid.linewidth': 0.5,
        'grid.color': 'lightgrey',
        'axes.axisbelow': True
    })

    # Boxplot properties
    BOXPROPS = {
        'boxprops': {'edgecolor': 'black', 'linewidth': 0.5},
        'medianprops': {'color': 'black', 'linewidth': 0.5},
        'whiskerprops': {'color': 'black', 'linewidth': 0.5},
        'capprops': {'color': 'black', 'linewidth': 0.5}
    }

    # Color palettes
    period_palette={'Infancy': '#F7C98A',
                    'Early_Childhood': '#F4A154',
                    'Late_Childhood': '#EF7A3B',
                    'Adolescence': '#D2614D',
                    'Early_Adulthood': '#9D4255',
                    'Adulthood': '#774147'}

    sex_palette = {'Female': '#BD2136', 'Male': '#3953A4'}
    celltype_palette = {'GLU': '#4FAF48', 'GABA': '#397FBA'}
    
    period_order = ['Infancy', 'Early_Childhood', 'Late_Childhood', 'Adolescence', 'Early_Adulthood', 'Adulthood']

    return BOXPROPS, period_palette, sex_palette, celltype_palette, period_order

def read_meta():

    gene_meta = pd.read_csv('./data/gencode.v37.annotation.intragenic.bed.gz',
                            names=['chr', 'start', 'end', 'gene_id', 'strand', 'tmp', 'gene_name', 'fn'],
                            index_col='gene_id', sep='\t')
    gene_meta['len'] = gene_meta['end'] - gene_meta['start']
    genedict_nametoid = dict(zip(gene_meta.gene_name, gene_meta.index))
    genedict_idtoname = dict(zip(gene_meta.index, gene_meta.gene_name))

    PECmeta = pd.read_csv('./data/TableS1_Donor_metadata.csv', index_col=0)
    agedict = dict(zip(PECmeta.index, PECmeta['Age']))
    sexdict = dict(zip(PECmeta.index, PECmeta['Sex']))
    perioddict = dict(zip(PECmeta.index, PECmeta['Age group']))
    
    return gene_meta, genedict_nametoid, genedict_idtoname, PECmeta, agedict, perioddict