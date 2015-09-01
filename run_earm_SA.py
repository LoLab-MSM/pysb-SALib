from earm.lopez_embedded import model
from pysb_SA import run_SA
import numpy as np
t=np.linspace(0,20000,20000)
run_SA(model,t, sp_SA='aSmac',N=2)