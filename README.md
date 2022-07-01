# Consegna Esercizi LSN - Martino Zanetti  

![](img/gipeto_small.png)   

## Organizzazione dei file in ogni cartella

🔴 In ciascuna cartella numerata si trovano:

- **Jupyter Notebook** `my-<exercise_number>.ipynb` contenente la spiegazione degli esercizi con l'esposizione dei risultati in forma di grafici
- **Codice sorgente**  `main.cpp` (o nomi più specifici, sempre con estensione `.cpp`)
- **Output files**     coi risultati delle simulazioni: `*.out` o `output.*` o simili (solitamente organizzati in cartelle)
- **Makefile**         (nota: spesso il comando `make` esegue direttamente il programma. <!-- fare ordine su questa cosa-->)
- **Eseguibile**       `main.exe` (o nomi più specifici, sempre con estensione `.exe`)

🔴 Sono inoltre presenti tre librerie di funzioni, alle quali i programmi accedono all'occorrenza in fase di compilazione:

- **`mylib`**          una piccola libreria personale 
- **`random`**         libreria per la generazione di numeri casuali 
- **`tsplib`**         libreria per l'algoritmo genetico applicato al problema del commesso viaggiatore (esercizi 9 e 10)

🔴 In alcune cartelle sono presenti anche:

- **README**           quando sono necessarie ulteriori spiegazioni sull'organizzazione della cartella o sull'esercitazione 
- **Header files**    `main.h`(o nomi più specifici, sempre con estensione `.h`)
- **Input files**      coi parametri della simulazione: `input.dat` o `input.*`, o simili
- **Config files**     `config.*` con le configurazioni del sistema (sia input, che output)
