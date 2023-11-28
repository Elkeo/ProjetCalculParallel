import numpy as np

def calculateSU_EF(filename, n):
   tab = np.zeros((4, n))
   count = 0
   with open(filename, "r+") as fichier:
      for ligne in fichier:
         value = ligne.strip().split()
         tab[0, count] = int(value[0])
         tab[1, count] = float(value[1])
         count +=1
      tab[2, :] = tab[1, 0]/tab[1, :]
      tab[3, :] = tab[2, :]/tab[0, :]
      np.savetxt(filename, np.transpose(tab))

calculateSU_EF("metis/time.dat", 8)
calculateSU_EF("scotch/time.dat", 8)