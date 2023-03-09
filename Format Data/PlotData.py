import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("difractograma-anatasa.csv", names=["x","y"])
plt.scatter(data["x"], data["y"])
plt.title("DIFRACTOGRAMA DE DATOS")
plt.show()