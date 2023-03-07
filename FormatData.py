import pandas as pd

raw_data = pd.read_csv("TiO2-Anatasa_PD_datos-procesados_2.csv", names=["c1"])


for i in range(len(raw_data)):
    try:
        a = str(raw_data.iloc[i]["c1"]).split(" ")
        x = a[len(a)-2]
        y = a[2]
        print(y)
    except:
        print(a)
        #print(a[7])

# raw_data = pd.read_csv("TiO2-Anatasa_PD_datos-procesados_2.csv", names=["x", "y"])
