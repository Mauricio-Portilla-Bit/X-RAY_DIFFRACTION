import pandas as pd

raw_data = pd.read_csv("TiO2-Anatasa_PD_datos-procesados_2.csv", names=["c1"])

x = []
y = []

for i in range(len(raw_data)):
    try:
        a = str(raw_data.iloc[i]["c1"]).split(" ")
        y.append(a[len(a)-2])

        if a[2] != "":
            x.append(a[2])
        else:
            x.append(a[3])

    except:
        continue
        print(a)

data = pd.DataFrame({
    "x": x,
    "y": y
})

print(data)

# raw_data = pd.read_csv("TiO2-Anatasa_PD_datos-procesados_2.csv", names=["x", "y"])
