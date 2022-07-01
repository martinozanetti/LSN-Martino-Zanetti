# Esercitazione 6

Per **compilare** il codice, entrare nella cartella Ising_1D e comandare
```bash
make
```
Per **eseguire** normalmente, nella stessa cartella comandare
```bash
make run
```

Per semplificare l'**esecuzione dell'intera simulazione**, ho scritto un file `bash` eseguibile con
```bash
./run_all_temp.sh
```
il quale esegue 10 simulazioni su 10 temperature diverse tra $T=0.1$ e $T=3$. Per riempire tutti i file relativi alla simulazione è necessario eseguire lo script:
- ad $h=0$
- ad $h=0.02$
- ad $h=0.2$

con algoritmo di Metropolis e di Gibbs, quindi in totale 6 volte (modificando di volta in volta il file `input.dat`). Nota: quando $h\neq 0$ i file per i plot di energia, suscettività e capacità termica non vengono scritti; viceversa, quando $h=0$, non vengono scritti solo quelli della magnetizzazione.

Per pulire i file dei plot, eseguire
```bash
./clean_plots.sh
```

**Nota**: i link interni nel Jupyter funzionano solo se il notebook è aperto con un editor adeguato (e.g. jupyter-notebook, visual studio, ...).
