import numpy as np
import pandas as pd
import loompy

ds = loompy.connect("/home/thimmamp/Downloads/l5_all.loom")
cluster_names = set(ds.ca.ClusterName)
row = ds.ra.Accession
frame = pd.DataFrame(data=ds[:,:],columns=ds.ca.CellID,index=row)
frame.to_csv(f"Linnarson_raw_count_data.csv", sep=',')
