import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
from pysam import VariantFile
import matplotlib.pyplot as plt

quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)

plt.savefig(snakemake.output[0])