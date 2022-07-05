# Consegna Esercizi LSN - Martino Zanetti  

![](img/gipeto_small.png)   

## Organizzazione dei file in ogni cartella

🔴 In ciascuna cartella numerata si trovano:

- **Jupyter Notebook** `my-<exercise_number>.ipynb` con la spiegazione degli esercizi e l'esposizione dei risultati
- **README**           con istruzioni sintetiche per l'esecuzione ed eventuali informazioni sull'organizzazione dei file nella cartella 

e una ulteriore cartella numerata in cui è racchiuso il resto:
- **Codice sorgente**  `main.cpp` (o nomi più specifici, sempre con estensione `.cpp`)
- **Output files**     coi risultati delle simulazioni: `*.out` o `output.*` o simili (solitamente organizzati in cartelle)
- **Makefile**         (nota: spesso il comando `make` esegue direttamente il programma. <!-- fare ordine su questa cosa-->)
- **Eseguibile**       `main.exe` (o nomi più specifici, sempre con estensione `.exe`)
- **Header files**     (eventuali) `main.h`(o nomi più specifici, sempre con estensione `.h`)
- **Input files**      (eventuali) coi parametri della simulazione: `input.dat` o `input.*`, o simili
- altri file e cartelle con nomi più possibile autoesplicativi o comprensibili leggendo il Notebook

🔴 Sono inoltre presenti tre librerie di funzioni, alle quali i programmi accedono all'occorrenza in fase di compilazione:

- **`mylib`**          una piccola libreria personale 
- **`random`**         libreria per la generazione di numeri casuali 
- **`tsplib`**         libreria per l'algoritmo genetico applicato al problema del commesso viaggiatore (esercizi 9 e 10)

🔴 E per non far mancare niente ho messo pure la [licenza (GPU3)](LICENCE)

