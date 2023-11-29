La compilation s'effectue avec un terminal ouvert dans le dossier contenant le Makefile.
On lance la commande :

	make

la compilation s'exécute avec les flags optimisés.

On change les paramètres désirés dans le fichiers parameters.txt

Pour lancer le programme, on écrit :

	mpirun --mca btl_vader_backing_directory /tmp -np [nproc] ./main

avec :
 - [nproc] un entier étant le nombre de processus souhaités (nproc >= 2)
 - "--mca btl_vader_backing_directory /tmp" qui est une commande pour éviter un troncature du dossier local, qui fait crasher MPI sur Mac quand on ajoute beaucoup de processeurs.

Pour visualiser la solution, il faut simplement aller chercher le fichier GIF dans le dossier "solutioni" correspondant au cas test i (les fichiers GIF et gnuplot sont créés automatiquement).

Pour accéder aux erreurs, elles sont tracées en fonction de l'itération dans un fichier errorFile.dat.