try:
    import umap
except:
    subprocess.run(["pip", "install", "umap-learn"])

from vcf2gwas.utils import *
from vcf2gwas.parsing import *

argvals = "-h"

def main():
    Parser(argvals)


main()