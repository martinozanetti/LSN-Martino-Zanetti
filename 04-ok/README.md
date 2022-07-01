# Esercitazione 4

**Organizzazione file:**  
`frames`: scommentando la linea di codice `ConfXYZ()`, raccoglie le configurazioni istantanee in formato XYZ, per la generazione di eventuali animazioni
`input`: contiene file che servono per avviare la simulazione (eccetto `config.out` e `velocity.out`, necessari per il *restart* dopo l'equilibrazione), tra cui `input.<phase>` nei quali si settano le proprietà del sistema.  
`output`: contiene i file scritti durante la simulazione, raggruppati per fase. Per ogni fase c'è una ulteriore sottocartella dedicata all'equilibrazione, contenente i file `config.out` e `velocity.out` letti per il *restart*.

**Compilare:**
```bash
make
```

**Eseguire:**
```bash
./NVE.exe <phase>
```
con phase={solid, liquid, gas}

**...alternativa:**
```bash
./NVE.exe
```
poi, quando richiesto, digitare la fase.

**Pulizie:**  
`./clean.sh`: rimuove *forzatamente* tutti i file dalla cartella `frames`  
`./clean<phase>.sh`: rimuove *forzatamente* tutti i file `.dat` dalla directory `<phase>`  
`make clean`: rimuove tutti i file `.o` ed `.exe`  
