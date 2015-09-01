from earm.lopez_embedded import model
import numpy as np
t=np.linspace(0,20000,20000)
from pysb.tools.PySB_SA import run_SA
run_SA(model,t, sp_SA='aSmac',N=2)