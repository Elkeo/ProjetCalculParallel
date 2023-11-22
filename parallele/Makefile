# Compilateur utilisé
CC=mpic++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG -std=c++11
# Options en mode debug - La variable est DEBUG est définie comme vraie
DEBUG_FLAG = -g -DDEBUG -std=c++11 -Wall -pedantic -O0 -Wextra

# On choisit comment on compile
CXX_FLAGS = $(OPTIM_FLAG)

# Le nom de l'exécutable
PROG = main

# Les fichiers source à compiler
SRC = main.cpp conjGrad.cpp prodMatVect_BC.cpp fonctions.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	find . \( -name \*.gif -o -name \*.dat -o -name \*.o -o -name \*.dSYM -o -name \$(PROG) -o -name \*~ \) -type f -delete