import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

# CLASE DE FIT PEAKS
class FitPeaks():


    # Instanciar una clase para fit peaks
    def __init__(self, df, H, A, ETA, BKG) -> None:

        # Instanciar Parámetros
        self.df = df
        self.H = H
        self.A = A
        self.ETA = ETA
        self.BKG = BKG

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
        self.SUMCHI2 = sum(((self.df["PV_BKG"] - self.df["y"])**2)/self.df["y"])

        self.df_test = self.df

    # Optimizar Parámetos
    def optimize_params(self):

        # MÉTODO DE OPTIMIZACIÓN : RECORRIDO SIMULADO

        # Constantes de la Optimización
        Temp = 400 # Temperatura Inicial
        alfa = 1 # Mecanismo de descenso
        L = 200 # Número de Iteraciones en cada nivel
        Tempf = 0 # Temperatura Final
        Delta = 0

        # Vecino Aleatorio para la Posición actual
        H_curr = random.random()
        ETA_curr = random.random()
        A_curr = random.random() * 1000
        BKG_curr = random.random() * 100

        # ITERACIONES
        while Temp > Tempf:

            # Iterar L veces el model o
            for i in range(L):

                # Vecino Aleatorio para una Posición posible
                H_pos = random.random()
                ETA_pos = random.random()
                A_pos = random.random()*1000
                BKG_pos = random.random()*100

                # Probar la Posición
                SUMCHI2_pos = self.evaluate_function(H_pos, ETA_pos, A_pos, BKG_pos)
                SUMCHI2_curr = self.evaluate_function(H_curr, ETA_curr, A_curr, BKG_curr)

                # Evaluar si el error aumentó o disminuyó
                Delta = SUMCHI2_pos - SUMCHI2_curr

                if Delta < 0 or random.random() < np.exp(-Delta/Temp):

                    H_curr = H_pos
                    ETA_curr = ETA_pos
                    A_curr = A_pos
                    BKG_curr = BKG_pos

                    if SUMCHI2_pos - self.SUMCHI2 < 0:
                        self.H = H_pos
                        self.ETA = ETA_pos
                        self.A = A_pos
                        self.BKG = BKG_pos
                        self.update_vales()

            Temp = Temp - alfa
            print(str(Temp), ")", str(self.SUMCHI2))
        return {"H": self.H, "A": self.A, "ETA": self.ETA, "BKG": self.BKG, "SUMCHI2": self.SUMCHI2}

    # Retornar Parámetros
    def get_params(self):
        return {"H": self.H, "A": self.A, "ETA": self.ETA, "BKG": self.BKG, "SUMCHI2": self.SUMCHI2}


    # Evaluar la función dado unos parámetros
    def evaluate_function(self, H_, ETA_, A_, BKG_):
        df_test = self.df
        df_test["L"] = (1/(2*np.pi))*(H_ / ((H_/2)**2 + df_test["x - x0"]**2))
        df_test["G"] = (1/H_)*(np.sqrt((4*np.log(2))/np.pi))*np.exp(-4*np.log(2)*((df_test["x - x0"]/H_)**2))
        df_test["PV"] = A_*((ETA_*df_test["L"]) + (1 - ETA_)*df_test["G"])
        df_test["PV_BKG"] = df_test["PV"] + BKG_
        SUMCHI2_test = sum(((df_test["PV_BKG"] - df_test["y"])**2)/df_test["y"])
        return SUMCHI2_test

    # Calcular el error del modelo
    def sum_chi2(self):
        self.SUMCHI2 = sum(((self.df["PV_BKG"] - self.df["y"])**2)/self.df["y"])

    def update_vales(self):
        self.df["L"] = (1/(2*np.pi))*(self.H / ((self.H/2)**2 + self.df["x - x0"]**2))
        self.df["G"] = (1/self.H)*(np.sqrt((4*np.log(2))/np.pi))*np.exp(-4*np.log(2)*((self.df["x - x0"]/self.H)**2))
        self.df["PV"] = self.A*((self.ETA*self.df["L"]) + (1 - self.ETA)*self.df["G"])
        self.df["PV_BKG"] = self.df["PV"] + self.BKG
        self.SUMCHI2 = sum(((self.df["PV_BKG"] - self.df["y"])**2)/self.df["y"])

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

#H_0 = 0.54
#A_0 = 664.2
#ETA_0 = 0.744
#BKG_0 = 23.28

#H_0 = 0.521486278631984
#A_0 = 642.0679663975903
#ETA_0 = 0.7670427175795259
#BKG_0 = 25.737969372925107

H_0 = 0.2
A_0 = 120
ETA_0 = 0.5
BKG_0 = 20

a = FitPeaks(raw_data, H_0, A_0, ETA_0, BKG_0)
#a.graph()
#v = a.get_df()
# print(v)
d = a.optimize_params()
print(d)


# OBLIGATORIO
# - Encontrar los máximos
# - Scherrer
# - WH
# - HW
# - Determinar el tamaño de partícula
# - Gráficas de los métodos

# NICE TO HAVE:
# - Warren Averbach (Convolusión)
# - Determinar el elemento presente en el programa
