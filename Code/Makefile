# Compilateur utilisé
CXX=mpic++
OMPI_CXX=g++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG -std=c++11
# Options en mode debug - La variable est DEBUG est définie comme vraie
DEBUG_FLAG = -O0 -g -DDEBUG -std=c++11 -Wall -Wextra -pedantic -fsanitize=undefined,address
#DEBUG_FLAG = -O0 -g3 -Wall -Wextra -Wconversion -Wdouble-promotion -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion -fsanitize=undefined -fsanitize-trap
# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = main

# Les fichiers source à compiler
SRC = main.cpp conjGrad.cpp prodMatVect_BC.cpp fonctions.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CXX) $(SRC) $(CXX_FLAGS) -o $(PROG)
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	find . \( -name \*.gif -o -name \*.dat -o -name \*.o -o -name \*.dSYM -o -name \$(PROG) -o -name \*~ \) -type f -delete