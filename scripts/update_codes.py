import pandas as pd
import os

codes_n_temps = pd.read_csv("codes_n_temps.csv")
codes_n_temps.set_index("dir", inplace=True)

data_list = os.listdir("../raw_data")

# updates codes_n_temps data
for prot in data_list:   
     prot_path = os.path.join("../raw_data", prot)   
     out = [i for i in os.listdir(prot_path) if "temp" in i]
     out = ["{}K".format(i[4:-1]) for i in out]

     # print directory name for protein
     print("== Protein directory: {} ==".format(prot))

     # print files in directory (for debugging)
     #print("== Files ==")
     #print(os.listdir(prot_path))

     # get all temperature directories
     print("== Found temperatures list ==")
     print(out, end="\n\n")

     if len(out) == 0:
         # do nothing
         pass
     else:
        # add row if it's not in the dataset
        if prot not in codes_n_temps.index:
            codes_n_temps.loc[prot] = False
        # update temperature data 
        for temp in out:
            codes_n_temps.loc[prot][temp] = True
             
# save updated data
codes_n_temps.to_csv("codes_n_temps.csv")
