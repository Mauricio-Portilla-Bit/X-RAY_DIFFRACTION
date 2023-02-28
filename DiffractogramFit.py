import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# CLASE DE FIT PEAKS
class FitPeaks():

    # Parámetros del Método
    H = 0.54
    A = 664.2
    ETA = 0.744
    BKG = 23.28
    SUMCHI = 0

    # Data Frame

    # Instanciar una clase para fit peaks
    def __init__(self, df, H, A, ETA, BKG, SUMCHI) -> None:

        # Instanciar Parámetros
        self.df = df
        self.H = H
        self.A = A
        self.ETA = ETA
        self.BKG = BKG
        self.SUMCHI = 0

        # Desplazar los datos al cero
        y_peak = max(self.df["y"])
        x_peak = float(df.loc[df["y"] == y_peak]["x"])
        self.df["x - x0"] = self.df["x"] - x_peak

        # Instanciar Lagrangiana
        self.df["L"] = (1/(2*np.pi))*(self.H / ((self.H/2)**2 + self.df["x - x0"]**2))

        # Instanciar Gaussiana
        self.df["G"] = (1/self.H)*(np.sqrt((4*np.log(2))/np.pi))*np.exp(-4*np.log(2)*((self.df["x - x0"]/self.H)**2))

        # Instanciar PseudoVoight
        self.df["PV"] = self.A*((self.ETA*self.df["L"]) + (1 - self.ETA)*self.df["G"])

        # Agregar el Background
        self.df["PV_BKG"] = self.df["PV"] + self.BKG

        # Instanciar el Error Cuadrado
        self.SUMCHI = sum(((self.df["PV_BKG"] - self.df["y"])**2)/self.df["y"])

    # Optimizar Parámetos
    def optimize_params(self):


        # MÉTODO DE OPTIMIZACIÓN
        #
        #
        #

        return {"H": self.H, "A": self.A, "ETA": self.ETA, "BKG": self.BKG, "SUMCHI": self.SUMCHI}

    # Retornar Parámetros
    def get_params(self):
        return {"H": self.H, "A": self.A, "ETA": self.ETA, "BKG": self.BKG, "SUMCHI": self.SUMCHI}

    # Probar el modelo
    def test_params(self):
        self.generate_function()
        self.sum_chi2()

    # Crear PseudoVoight + BKG
    def generate_function(self):
        self.df["L"] = (1/(2*np.pi))*(self.H / ((self.H/2)**2 + self.df["x - x0"]**2))
        self.df["G"] = (1/self.H)*(np.sqrt((4*np.log(2))/np.pi))*np.exp(-4*np.log(2)*((self.df["x - x0"]/self.H)**2))
        self.df["PV"] = self.A*((self.ETA*self.df["L"]) + (1 - self.ETA)*self.df["G"])
        self.df["PV_BKG"] = self.df["PV"] + self.BKG

    # Calcular el error del modelo
    def sum_chi2(self):
        self.SUMCHI = sum(((self.df["PV_BKG"] - self.df["y"])**2)/self.df["y"])

    # Graficar los puntos y el modelo
    def graph(self):
        plt.figure(1)
        plt.plot(self.df["x - x0"], self.df["PV_BKG"], c="b")
        plt.scatter(self.df["x - x0"], self.df["y"], c="r")
        plt.grid()
        plt.xlabel("2*Theta")
        plt.ylabel("Intensidad")
        plt.title("DIFRACTOGRAMA")
        plt.show()

    # Devolver el data frame
    def get_df(self):
        return self.df


# DATOS EXPERIMENTALES
raw_data = pd.read_csv("muestra_experimental.csv", names=["x", "y"])

H_ = 0.54
A_ = 664.2
ETA_ = 0.744
BKG_ = 23.28
SUMCHI_ = 0

a = FitPeaks(raw_data, H_, A_, ETA_, BKG_, SUMCHI_)
v = a.get_df()
print(v)
d = a.optimize_params()
print(d)
a.graph()
