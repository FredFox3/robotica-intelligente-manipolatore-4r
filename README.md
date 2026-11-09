# Robot Manipolatore Planare 4R con Gravità

Progetto MATLAB per cinematica, dinamica di Eulero–Lagrange e controllo Computed
Torque di quattro link rigidi nel piano XY. Unità SI, angoli in radianti,
rotazioni positive antiorarie attorno a z uscente; gravità lungo −y.

## Esecuzione su MATLAB Online

Caricare i sette file `.m` nella stessa cartella di MATLAB Online e selezionarla
come cartella corrente. La simulazione usa MATLAB base; solo la derivazione
simbolica richiede Symbolic Math Toolbox, disponibile nell'ambiente previsto.

```matlab
main_4R
```

Il programma esegue gli assert preliminari, integra con `ode45`, controlla la
convergenza, stampa errori e picchi delle coppie, poi mostra grafici e animazione.
Per derivare e confrontare separatamente il modello simbolico:

```matlab
model = derive_dynamics_4R();
model.M  % Espressioni simboliche; disponibili anche C, G, kinetic, potential, tauEL
```

La derivazione non scrive file: confronta M, C, G, energie ed equazioni complete
con la funzione numerica in tre casi. Non è richiamata dentro `ode45`.

## Modello e scelte

La tabella DH standard è:

| i | a_i | alpha_i | d_i | theta_i |
|---|---|---|---|---|
| 1 | L_1 | 0 | 0 | q_1 |
| 2 | L_2 | 0 | 0 | q_2 |
| 3 | L_3 | 0 | 0 | q_3 |
| 4 | L_4 | 0 | 0 | q_4 |

`fk_4R` moltiplica le quattro matrici DH. `ik_4R(x,y,phi,L)` sottrae il link 4
e sceglie analiticamente l'orientamento del terzo link dall'intervallo che rende
il polso residuo raggiungibile dal 2R. Il quinto argomento opzionale `+1`/`-1`
sceglie il segno di q_2 (i due gomiti coincidono nelle singolarità).
L'IK restituisce un vettore colonna con angoli equivalenti in [−pi,pi]; verifica
la posa con FK e segnala le pose irraggiungibili con `ik_4R:Unreachable`.
Sono ammesse lunghezze strettamente positive, senza limiti articolari.
La scelta deterministica di una soluzione IK può essere discontinua tra pose.

`dynamics_4R` assembla M dai Jacobiani dei baricentri, C dalle loro derivate
analitiche e dai simboli di Christoffel, G dal gradiente di V = sum(m_i*g*y_ci).
`params` contiene `L`, `m`, `I` (quattro elementi ciascuno) e `g` scalare.
Nel main le inerzie sono quelle di barre omogenee, I_i = m_i*L_i²/12.
Il controllo usa tau = M*(ddq_d + Kd*de + Kp*e) + C*dq + G; il riferimento
quintico ha velocità e accelerazione nulle agli estremi e prosegue con un tratto
fermo. Il modello è ideale e il controllore usa gli stessi parametri del robot;
non sono inclusi attrito, saturazioni o robustezza a errori di modello.

`plot_robot_4R` disegna le cinque terne DH; `animate_4R` aggiorna link, assi,
target e tracce con `drawnow`. L'animazione non garantisce il tempo reale.

## Validazione su MATLAB Online

Eseguiti con successo `main_4R` e `derive_dynamics_4R`: completati i controlli
FK/IK, irraggiungibilità, energie, M, C e G e la simulazione con `ode45`.
La derivazione simbolica è stata eseguita correttamente, con confronto numerico
riuscito in tre casi.

- Errore articolare massimo: da 6.981e-02 a 9.992e-16 rad.
- Errori finali di posa: 4.578e-16 m in posizione e 3.331e-16 rad in orientamento.
- Picchi delle coppie [N m]: [25.98 10.5 3.119 0.4644].

Cinque figure e un'animazione MATLAB sono state controllate visivamente.
Le cinque figure esportate in PNG sono disponibili nella cartella [`results/`](results/).
