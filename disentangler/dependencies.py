import sys
sys.path.insert(0, '/liulab/alynch/projects/multiomics/BatchEffect/MIRA')

import mira
assert mira.__file__ == '/liulab/alynch/projects/multiomics/BatchEffect/MIRA/mira/__init__.py'

import scanpy as sc
import seaborn as sns
import numpy as np
import pandas as pd
import seaborn as sns
import anndata
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_samples