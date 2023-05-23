import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('output.txt', delimiter=',')

chosen_arms = data[:,:2]
lossescomp = data[:,2]
lossesdef = data[:,3]

# Estilo de las gráficas
plt.style.use('seaborn-darkgrid')

# Grafico 1
plt.figure(figsize=(12, 8))
plt.plot(chosen_arms[:,0], label='GPMW', color='red', linewidth=2)
plt.plot(chosen_arms[:,1], label='cGPMW', color='blue', linewidth=2)
plt.title('Comparación de brazos elegidos', fontsize=18)
plt.xlabel('Rondas', fontsize=14)
plt.ylabel('Brazos elegidos', fontsize=14)
plt.grid(True)
plt.legend(fontsize=12)
plt.show()

# Grafico 2
plt.figure(figsize=(12, 8))
plt.plot(lossesdef, label='cGPMW', color='blue', linewidth=2)
plt.plot(lossescomp, label='GPMW', color='red', linewidth=2)
plt.title('Pérdidas por ronda', fontsize=18)
plt.xlabel('Pérdidas', fontsize=14)
plt.ylabel('Ronda', fontsize=14)
plt.grid(True)
plt.legend(fontsize=12)
plt.show()

# Grafico 3
plt.figure(figsize=(12, 8))
plt.plot(np.cumsum(lossesdef), label='cGPMW', color='blue', linewidth=2)
plt.plot(np.cumsum(lossescomp), label='GPMW', color='red', linewidth=2)
plt.title('Pérdida acumulada', fontsize=18)
plt.xlabel('Rondas', fontsize=14)
plt.ylabel('Pérdidas acumuladas', fontsize=14)
plt.grid(True)
plt.legend(fontsize=12)
plt.show()


# Cálculo de las pérdidas medias
average_loss_def = np.cumsum(lossesdef) / (np.arange(len(lossesdef)) + 1)
average_loss_comp = np.cumsum(lossescomp) / (np.arange(len(lossescomp)) + 1)

# Gráfico 4
plt.figure(figsize=(12, 8))
plt.plot(average_loss_def, label='cGPMW', color='blue', linewidth=2)
plt.plot(average_loss_comp, label='GPMW', color='red', linewidth=2)
plt.title('Pérdida media por ronda', fontsize=18)
plt.xlabel('Rondas', fontsize=14)
plt.ylabel('Pérdida media', fontsize=14)
plt.grid(True)
plt.legend(fontsize=12)
plt.show()
