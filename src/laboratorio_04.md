# Laboratorio 4: Sistemi lineari triangolari e fattorizzazione LU


In questo laboratorio ci focalizzeremo sul problema di risolvere
un sistema lineare. Partiremo dal caso dei sistemi triangolari per
poi occuparci del caso generale tramite la fattorizzazione $PA = LU$ della
matrice dei coefficienti.


## Soluzione di sistemi triangolari

I sistemi triangolari giocano un ruolo fondamentale nei calcoli matriciali.
Molti metodi sono costruiti sull'idea di ridurre un problema alla soluzione di uno o più sistemi triangolari, questo include praticamente tutti i **metodi diretti** per la risoluzione di sistemi lineari.

Due funzioni utili per lavorare con i sistemi triangolari sono le funzioni
`triu` e `tril`. Dal loro `help`:
```
tril Extract lower triangular part.
   tril(X) is the lower triangular part of X.
   tril(X,K) is the elements on and below the K-th diagonal
   of X .  K = 0 is the main diagonal, K > 0 is above the
   main diagonal and K < 0 is below the main diagonal.
```
ed analogamente per `triu`:
```
triu Extract upper triangular part.
   triu(X) is the upper triangular part of X.
   triu(X,K) is the elements on and above the K-th diagonal of
   X.  K = 0 is the main diagonal, K > 0 is above the main
   diagonal and K < 0 is below the main diagonal.
```
Possiamo usarle per costruire sistemi triangolari a partire da sistemi densi.

(forwardandbacwardsolve)=
## Sostituzione in avanti e all'indietro

Sui **computer seriali** i sistemi triangolari sono universalmente risolti dagli algoritmi standard di sostituzione in avanti e all'indietro. Sulle **macchine
parallele** la situazione è sensibilmente più varia, ma questo trascende gli
obiettivi di questo corso.

Sia $L$ una matrice triangolare inferiore, vogliamo risolvere un sistema della forma:
```{math}
L \mathbf{x} = \mathbf{b},
```
che possiamo riscrivere in forma estesa come:
```{math}
\begin{bmatrix}
l_{1,1} & 0       & 0      & \cdots & 0 \\
l_{2,1} & l_{2,2} & 0      & \cdots & 0 \\
\vdots  & \ddots  & \ddots & \ddots & \vdots \\
l_{n-1,1} & l_{n-1,2} & \cdots & l_{n-1,n} & 0 \\
l_{n,1} & l_{n,2} & \cdots & \cdots & l_{n,n}  \\
\end{bmatrix} \begin{bmatrix}
x_1 \\
x_2 \\
\vdots \\
x_{n-1}\\
x_n
\end{bmatrix}
= \begin{bmatrix}
b_1 \\
b_2 \\
\vdots \\
b_{n-1} \\
b_n
\end{bmatrix}
```
Da cui è facile ricavare che la componente $i$esima della soluzione è
```{math}
x_i = \frac{1}{l_{i,i}} \left( b_i - \sum_{j=1}^{i-1} l_{i,j} x_j \right), \qquad i=1,\ldots,n.
```

:::{admonition} Esercizio 1
Si scriva una *function* che implementa il metodo di sostituzione in avanti con
il seguente prototipo.
```matlab
function [x] = forwardsolve(A,b)
% FORWARDSOLVE Se A è una matrice triangolare inferiore risolve
% il sistema Ax=b tramite metodo di sostituzione in avanti.
%   Input:
%   A = matrice triangolare inferiore
%   b = vettore del termine noto
%   Output:
%   x = vettore della soluzione

end
```
Si presti attenzione:
- alla verifica degli *input*, il sistema **deve** essere triangolare inferiore e
quadrato;
- ad utilizzare le **operazioni vettorizzate** di MATLAB, ovvero si ragioni su
come evitare di usare un ciclo `for` per effettuare la somma che compare
nella formula di soluzione.

Si può testare il codice con il seguente problema:
```matlab
%% Sostituzione in avanti
n = 10;
A = tril(ones(n,n));
b = (1:n).';
x = forwardsolve(A,b);

fprintf("Sostituzione in avanti: Errore sulla soluzione è: %1.2e\n",norm(A*x-b));
```
:::

Sia $U$ una matrice triangolare superiore, vogliamo risolvere un sistema della forma:
```{math}
U \mathbf{x} = \mathbf{b},
```
che possiamo riscrivere in forma estesa come:
```{math}
\begin{bmatrix}
u_{1,1} & u_{1,2} & \cdots & \cdots & u_{1,n}  \\
0 & u_{2,2} & u_{2,3} & \cdots & u_{2,n} \\
\vdots  & \ddots  & \ddots & \ddots & \vdots \\
0      & \cdots & 0 & u_{n-1,n-1} & u_{n-1,n}  \\
0       & 0      & \cdots & 0 & u_{n,n} \\
\end{bmatrix} \begin{bmatrix}
x_1 \\
x_2 \\
\vdots \\
x_{n-1}\\
x_n
\end{bmatrix}
= \begin{bmatrix}
b_1 \\
b_2 \\
\vdots \\
b_{n-1} \\
b_n
\end{bmatrix}
```
Da cui di nuovo ricaviamo facilmente la componente $i$ma della soluzione come:
```{math}
x_i = \frac{1}{u_{i,i}} \left( b_i - \sum_{j=i+1}^{n} u_{i,j} x_j \right), \qquad i=n,n-1,\ldots,1.
```

:::{admonition} Esercizio 2
Si scriva una *function* che implementa il metodo di sostituzione all'indietro con
il seguente prototipo.
```matlab
function [x] = backwardsolve(A,b)
% BACKWARDSOLVE Se A è una matrice triangolare superiore risolve
% il sistema Ax=b tramite metodo di sostituzione all'indietro.
%   Input:
%   A = matrice triangolare inferiore
%   b = vettore del termine noto
%   Output:
%   x = vettore della soluzione

end
```
Si presti attenzione:
- alla verifica degli *input*, il sistema **deve** essere triangolare superiore e
quadrato;
- ad utilizzare le **operazioni vettorizzate** di MATLAB, ovvero si ragioni su
come evitare di usare un ciclo `for` per effettuare la somma che compare
nella formula di soluzione.

Si può testare il codice con il seguente problema:
```matlab
%% Sostituzione all'indietro
n = 10;
A = triu(ones(n,n));
b = (1:n).';
x = backwardsolve(A,b);

fprintf("Indietro: Errore sulla soluzione è: %1.2e\n",norm(A*x-b));
```
:::

## Tempo di calcolo

Si è visto a lezione che il numero di operazioni necessario a risolvere un
sistema triangolare con l'algoritmo di sostituzione è un $O(n^2)$ con $n$ la
dimensione del sistema. Proviamo ad indagare come si comporta a questo riguardo
la nostra implementazione. Per farlo sfruttiamo i codici che abbiamo appena
scritto per la sostituzione _forward_ o _backward_ e le funzioni `tic` e `toc`.
```matlab
sizes = floor(logspace(1,4,10));
times = zeros(size(sizes));
h = waitbar(0,"Calcolo in corso...");
for i=1:length(sizes)
    n = sizes(i);
    A = triu(ones(n,n));
    b = (1:n).';
    tic;
    x = backwardsolve(A,b); % Allo stesso modo con forwardsolve
    times(i) = toc;
    waitbar(i/length(sizes));
end
```
Ora che abbiamo i tempi possiamo visualizzare quello che abbiamo ottenuto
sfruttando alcune funzioni di MATLAB. Poiché sappiamo che il numero di
operazioni è una funzione che scala come un $O(n^2)$ possiamo assumere in
prima approssimazione che anche il tempo si comporti allo stesso modo.
Cerchiamo dunque di *fittare* un polinomio di secondo grado ai nostri dati:
```matlab
[P,S] = polyfit(sizes,times,2);
```
Con la mia macchina, questa procedure mi restituisce il polinomio:
```{math}
p(x) = 0.000005682501880 x^2  -0.015569138451058 x + 3.230452264767748
```
insieme ad una struttura $S$ che contiene informazioni relative all'algoritmo
con cui questa procedura è stata eseguita. Possiamo usare tutto questo per
confrontare i valori del *fit* con quelli dei dati
```matlab
[y_fit,delta] = polyval(P,sizes,S);
plot(sizes,times,'o-',...
    sizes,y_fit,'--',...
    sizes,y_fit+2*delta,'m--',sizes,y_fit-2*delta,'m--',...
    'Linewidth',2);
xlabel('Dimensione del sistema')
ylabel('Tempo di calcolo (s)')
legend({'Tempo misurato','Stima asintotica','Intervallo di Confidenza'},...
    'FontSize',14,'Location','northwest');
```
dove abbiamo stampato sia il tempo misurato, sia il *fit* con il suo
*intervallo di confidenza* entro il 95%. Dalla {numref}`triangularsolvetime`
vediamo che la predizione quadratica è piuttosto efficace.
```{figure} ./images/triangulartime.png
:name: triangularsolvetime

Tempo di calcolo come un $O(n^2)$ della dimensione.
```

 <!---
## Esercizi

:::{admonition} Esercizio 3: Numeri di Bernoulli

Consideriamo le condizioni
```{math}
B(x+1) - B(x) = n x^{n-1}, \quad \int_{0}^{1} B(x)\,{\rm d}x=0, \quad B(x)\, \text{polinomio}.
```
Queste definiscono in maniera unica la funzione $B(x)$. Se assumiamo che il grado
del polinomio monico $B(x)$ sia $n$ abbiamo così costruito i **polinomi di Bernoulli**. Si possono calcolare esattamente i primi polinomi come:
```{math}
B_0(x) = 1, \quad B_1(x) = x - \frac{1}{2}, \quad B_2(x) = x(x-1)+\frac{1}{6}, \quad B_3(x) = x(x-\frac{1}{2})(x-1), \ldots
```
Si può dimostrare che i polinomi di Bernoulli definiscono i coefficienti della rappresentazione in serie di potenze di diverse funzioni, ad esempio:
```{math}
\frac{t e^{xt}}{e^t - 1} = \sum_{k=0}^{+\infty} \frac{B_n(x)}{n!}t^n.
```
In particolare, se scegliamo $x = 0$, abbiamo che
```{math}
:label: bernoulliexpansion
\frac{t}{e^t - 1} = -\frac{1}{2}t + \sum_{k=0}^{+\infty} \frac{B_{2k}(0)}{(2k)!}t^{2k}.
```
Facciamo ora la seguente manipolazione algebrica, moltiplichiamo l'equazione
a sinistra e a destra per $e^{t} - 1$, espandiamo $e^{t}$ con la sua serie
di potenze, e poniamo uguale a zero i coefficienti di $t^i$ sul lato destro.
Così otteniamo il seguente sistema di equazioni:
```{math}
- \frac{1}{2} j + \sum_{k=0}^{[\frac{j-1}{2}]} \binom{j}{2k} B_{2k}(0), \qquad j=2,3,4,\ldots
```
Da cui otteniamo, per gli $j$ pari, il seguente sistema di equazioni lineari
```{math}
\begin{bmatrix}
\binom{2}{0} \\
\binom{4}{0} & \binom{4}{2} \\
\binom{6}{0} & \binom{6}{2} & \binom{6}{4} \\
\binom{8}{0} & \binom{8}{2} & \binom{8}{4} & \binom{8}{6} \\
\vdots & \vdots & \vdots & \vdots & \ddots \\
\end{bmatrix}
\begin{bmatrix}
B_0(0)\\
B_2(0)\\
B_4(0)\\
B_6(0)\\
\vdots
\end{bmatrix}
= \begin{bmatrix}
1\\
2\\
3\\
4\\
\vdots
\end{bmatrix}.
```
1. Si implementi una funzione che dato un intero $k$ restituisca i numeri di
Bernoulli $\{B_{2j}(0)\}_{j=0}^{k}$ risolvendo questo sistema lineare.
2. Si usino i numeri così ottenuti per determinare lo sviluppo in {eq}`bernoulliexpansion` e disegnare
- le due funzioni una accanto all'altra,
- l'errore di approssimazione commesso in scala logaritmica.

Un prototipo della funzione è:
```matlab
function [Bk] = nbernoulli(k)
% NBERNOULLI Produce i numeri di Bernoulli B_{2j}(0) per j=0,...,k
% risolvendo un sistema triangolare inferiore.
%   INPUT:
%   k =  Numero di termini (-1) da calcolare
%   OUTPUT:
%   Bk = Vettore di lunghezza k+1 che contiene i numeri richiesti
end
```

```{tip}
La funzione binomiale in MATLAB è implementata come `nchoosek`. Potete inoltre
confrontare i risultati ottenuti con la funzione nativa di MATLAB che calcola
i numeri di Bernoulli dal medesimo nome (`help bernoulli` per le informazioni).
```
:::
 -->


 <!---
(condizionamento-spettrale)=
## Condizionamento di un sistema lineare

Ricordiamo che il **numero di condizionamento** associato all'equazione lineare
$A \mathbf{x} = \mathbf{b}$ fornisce un limite sulla precisione della soluzione
$\mathbf{x}$ dopo che un qualunque algoritmo sarà stato applicato per la
sua approssimazione ed in maniera indipendente dalla precisione *floating point*
usata.

:::{margin} Norma di una matrice
Ricordiamo che una norma sullo spazio vettoriale $K^{m\times n}$ delle matrici
a elementi nel campo $K$ è una funzione
```{math}
\| \cdot \|:K^{m \times n} \mapsto \mathbb{R}^{+}
```
tale che per ogni coppia di matrici $A$ e $B$ e per ogni scalare $\lambda \in K$ si verifica:

- ${\displaystyle \|A\|=0}$ se e solo se ${\displaystyle A=0}$,
- ${\displaystyle \|\lambda A\|=|\lambda |\|A\|}$.
- ${\displaystyle \|A+B\|\leq \|A\|+\|B\|}$.

Se è data una norma sullo spazio vettoriale $K^n$, si può costruire una **norma indotta**
sullo spazio delle matrici in $K^{n \times n}$ considerando
```{math}
\|A\| = \sup_{x \neq \mathbf{0}} \frac{\| A \mathbf{x} \|}{\|\mathbf{x}\|},
```
ovvero
```{math}
\|A\| = \sup_{\|x\|  = 1} \| A \mathbf{x} \|.
```
Nei casi speciali in cui la norma vettoriale sia la norma  ${\displaystyle p=1,2,\infty ,}$
la norma matriciale indotta può essere calcolata come
```{math}
{\displaystyle \|A\|_{1}=\max _{1\leq j\leq n}\sum _{i=1}^{m}|a_{ij}|,}
```
che è il massimo della somma dei valori assoluti delle colonne di $A$;
```{math}
{\displaystyle \|A\|_{\infty }=\max _{1\leq i\leq m}\sum _{j=1}^{n}|a_{ij}|,}
```
che è il massimo della somma dei valori assoluti delle righe di $A$;
```{math}
{\displaystyle \|A\|_{2}={\sqrt {\lambda _{\max }\left(A^{*}A\right)}}.}
```
che è detta **norma spettrale** della matrice ${\displaystyle A}$, ovver la
radice quadrata the più grande autovalore di ${\displaystyle A^{* }A}$,
dove $A^{* }$ denota la coniugata trasposta di $A$.
:::

Sia $\mathbf{e}$ l'errore commesso rispetto a $\mathbf{b}$, assumiamo che $A$ sia
non singolare, allora l'errore sulla soluzione $A^{-1} \mathbf{b}$ è dato da
$A^{-1} \mathbf{e}$. Allora il rapporto tra gli errori relativi è dato da
```{math}
\frac{\| A^{-1} \mathbf{e} \|}{\| A^{-1} \mathbf{b} \|} / \frac{\|\mathbf{e}\|}{\|\mathbf{b}\|} = \frac{\| A^{-1}\mathbf{e}\|}{\|\mathbf{e}\|} \frac{\|\mathbf{b}\|}{\|A^{-1}\mathbf{b}\|},
```
per cui possiamo esibire la limitazione per il termine destro come
```{math}
:label: bound
\frac{\| A^{-1}\mathbf{e}\|}{\|\mathbf{e}\|} \frac{\|\mathbf{b}\|}{\|A^{-1}\mathbf{b}\|} \leq & \max_{\mathbf{e},\mathbf{b} \neq 0} \left\lbrace \frac{\| A^{-1}\mathbf{e}\|}{\|\mathbf{e}\|}, \frac{\|\mathbf{b}\|}{\|A^{-1}\mathbf{b}\|} \right\rbrace \\ = & \max_{\mathbf{e}\neq 0} \frac{\| A^{-1}\mathbf{e}\|}{\|\mathbf{e}\|} \max_{\mathbf{b}\neq 0} \frac{\|\mathbf{b}\|}{\|A^{-1}\mathbf{b}\|} \\
= & \max_{\mathbf{e}\neq 0} \frac{\| A^{-1}\mathbf{e}\|}{\|\mathbf{e}\|} \max_{\mathbf{b}\neq 0} \frac{\|A\mathbf{x}\|}{\|\mathbf{x}\|} \\
= & \| A^{-1} \| \|A \| \equiv \kappa(A) \geq \| A^{-1} A \| = 1.
```
Consideriamo ora il seguente esempio artificiale
```{math}
A \mathbf{x} = \begin{bmatrix}
4.1 & 2.8 \\
9.7 & 6.6
\end{bmatrix}  \begin{bmatrix}
x_1 \\ x_2
\end{bmatrix} = \begin{bmatrix}
4.1 \\ 9.7
\end{bmatrix}
```
Proviamo a risolvere direttamente con MATLAB il sistema lineare
```matlab
A = [4.1 2.8; 9.7 6.6];
b = A(:,1);
x = A\b
```
da cui otteniamo come risposta:
```
x =

   1.000000000000002
  -0.000000000000003
```
Ritracciamo i passi del bound sul termine noto, per cui introduciamo una
perturbazione di $0.01$ sulla prima componente di $\mathbf{b}$ e paragoniamo le
due soluzioni
```matlab
b2 = [4.11; 9.7];
x2 = A\b2
```
da cui otteniamo
```
x2 =

   0.339999999999957
   0.970000000000064
```

Andiamo ora a verificare la limitazione che abbiamo mostrato in termini del
numero di condizionamento. A questo scopo facciamo uso del comando `cond`
```
cond   Condition number with respect to inversion.
   cond(X) returns the 2-norm condition number (the ratio of the
   largest singular value of X to the smallest).  Large condition
   numbers indicate a nearly singular matrix.

   cond(X,P) returns the condition number of X in P-norm:

      NORM(X,P) * NORM(INV(X),P).

   where P = 1, 2, inf, or 'fro'.
```
:::{danger}
Il calcolo del numero di condizionamento per una matrice è bene che sia eseguito
tramite la funzione `cond`, la forma `NORM(X,P) * NORM(INV(X),P)` nella
documentazione è identificativa di ciò che stiamo calcolando, ma *non* è
rappresentativo del modo in cui è realmente calcolato.
:::

Calcoliamo quindi
```matlab
kappa = cond(A);
bound = kappa*norm(b-b2)/norm(b);
errore_relativo = norm(x-x2)/norm(x);
```
e osserviamo che
```
bound =

   1.54117722269025
errore_relativo =

  1.173243367763135
```
che quindi ci mostra che l'errore che abbiamo commesso è abbastanza vicino
alla limitazione data dalla {eq}`bound`.

:::{warning}
Per matrici di grandi dimensioni la funzione `cond` può risultare estremamente
lenta o andare **out of memory**. Per questo motivo sono disponibili altre
due funzioni chiamate `rcond` e `condest`.

La funzione `rcond(A)` calcola il reciproco di `cond(A,1)`, per cui più questo è
vicino a $0$ più la matrice è *malcondizionata*, più questo è vicino a $1$ più
la matrice è *ben condizionata*. Questa funzione si basa sulla libreria [LAPACK](http://www.netlib.org/lapack/).
Questa è una libreria di software in Fortran 90 che produce le implementazioni
standard delle operazioni dell'algebra lineare per matrici dense.

La funzione `condest` calcola un limite inferiore per il numero di condizionamento
in norma $1$ di una matrice quadrata $A$. Informazioni sull'algoritmo applicato
per ottenere questa approssimazione si trovano in {cite:p}`MR740850` e {cite:p}`MR1780268`.
:::
-->


# Il metodo di eliminazione di Gauss e la fattorizzazione LU


Il metodo di **eliminazione di Gauss** (Carl Friedrich Gauss, 1777-1855), noto anche come
algoritmo di riduzione per le righe, è un algoritmo per risolvere sistemi
quadrati di equazioni lineari. Come avete visto a lezione si riduce ad una
sequenza di operazioni eseguite sulla corrispondente matrice di coefficienti.

Alcune estensioni di questo metodo possono essere utilizzate anche per calcolare
il **rango** di una matrice, il **determinante** di una matrice quadrata e
l`**inversa di una matrice** invertibile.

 <!---
Come abbiamo visto nel secondo Laboratorio sul {ref}`condizionamento-spettrale`
la precisione ottenuta in uscita da questa procedura. Ricordiamo che per
valutare il *numero di condizionamento* di una matrice $A$, MATLAB mette a
disposizione due comandi:

- `cond`: Numero di condizione rispetto all'inversione. `cond(A)` restituisce
il numero di condizionamento rispetto alla norma 2 (che ricordiamo essere il
    rapporto tra valore singolare più grande di $A$ e di quello più piccolo);
- `rcond`: stimatore del reciproco del numero di condizionamento in norma 1.
Utilizza lo stimatore implementato in LAPACK;
- `condest`: calcola un limite inferiore $C$ per il numero di condizionamento in
norma 1 di una matrice quadrata $A$. Accetta come seconda opzione un parametro
intero positivo $T$ che rapprenta il numero di colonne utilizzate per l'iterazione
dell'algoritmo di stima sottostante. L'aumentare del numero di colonne fornisce
*di solito* una stima migliore del condizionamento, al prezzo di aumentarne il
costo computazionale. Il valore predefinito è $T = 2$, che fornisce
*quasi sempre* una stima corretta entro un fattore 2.

Si provi ad eseguire il seguente comando:
```matlab
A = gallery('grcar',8);
[cond(A,1) 1/rcond(A) condest(A)]
[cond(A,1) 1/rcond(A) condest(A)]
```
Come illustra questo esempio, `condest` non restituisce necessariamente lo
stesso risultato ad ogni chiamata, poiché fa uso della funzione `rand` al suo
interno (si tratta di uno stimatore stocastico).
 -->
 
## Implementare la fattorizzazione LU

Se $A$ è una matrice $n$ per $n$ su un campo $\mathbb{K}$ ($A  \in M_n(\mathbb{K})$),
allora $A$ si dice "avere una fattorizzazione $LU$" se esiste una matrice
**triangolare inferiore** $L \in  M_n(\mathbb{K})$ e una matrice **triangolare
superiore** $U  \in M_n(\mathbb{K})$ tale che
```{math}
A = LU.
```
Questo è particolarmente utile, poiché, se $A = LU$, e se $L$ ed $U$ sono invertibile, allora ogni soluzione $\mathbf{y}$ di $L\mathbf{y} = \mathbf{b}$
produrrà una soluzione $\mathbf{x}$ di $A\mathbf{x} = \mathbf{b}$ via $U\mathbf{x} = \mathbf{y}$.

:::{danger}
In **generale** non tutte le matrici $A$ ammettono una fattorizzazione $LU$ di
questa forma, si consideri ad esempio la matrice:
```{math}
A = \begin{bmatrix}
0 & 1 \\
1 & 0
\end{bmatrix} = \begin{bmatrix}
l_{11} & 0 \\
l_{12} & l_{22}
\end{bmatrix} \begin{bmatrix}
u_{11} & u_{12} \\
0 & u_{22}
\end{bmatrix} = LU,
```
questa richiederebbe di avere gli elementi diagonali $u_{ii}$, $l_{ii}$, $i=1,2$,
diversi da $0$, ma questo contraddirebbe $0 = A_{11} = l_{11}u_{11}$.

Una condizione **necessaria e sufficiente** affinché una fattorizzazione
di questa forma esista per una matrice $A \in M_n(\mathbb{K})$,
$\mathbb{K} = \mathbb{R},\mathbb{C}$, è che
```{math}
\operatorname{rank}(A_{11}) + k \geq \operatorname{rank}([A_{11}\,A_{12}]) + \operatorname{rank}\left(\begin{bmatrix}A_{11}\\A_{21}\end{bmatrix}\right),
```
per $k = 1,\ldots,n-1$ e
```{math}
A = \begin{bmatrix}
A_{11} & A_{12} \\
A_{21} & A_{22}
\end{bmatrix}, \qquad A_{11} \in M_k(\mathbb{K}).
```
:::

Quando esiste la decomposizione $LU$ **non è unica** (le combinazioni di $L$ e $U$
per una data $A$ sono infinite). Per determinare un algoritmo univoco è
necessario porre determinati vincoli su $L$ o $U$ ({numref}`tipidilu`).
:::{margin} Varianti dell'algoritmo $LU$.
```{list-table} Varianti dell'algoritmo $LU$.
:header-rows: 1
:name: tipidilu


* - Nome
  - Vincoli
* - Doolittle
  - $l_{i,i} = 1$, $i=1,\ldots,n$
* - Crout
  - $u_{i,i} = 1$, $i=1,\ldots,n$
* - Choleski
  - $U = L^T$
```
:::
Come primo esercizio implementeremo la fattorizzazione $LU$ nella forma di Doolittle,
che è quella più strettamente legata all'algoritmo di eleminazione di Gauss che
conoscete sin dal corso di *algebra lineare*.

Se scriviamo l'uguaglianza $A = LU$ con il vincolo aggiuntivo per cui
 $l_{i,i} = 1$, per $i=1,\ldots,n,$ possiamo ricavare che i termini della matrice
 $U$ sono dati da
```{math}
:label: lu1
\forall\,j\quad & \begin{array}{ll}
i = 1, & u_{i,j} = a_{i,j},\\
i > 1, & u_{i,j} = \displaystyle a_{i,j} - \sum_{k=1}^{i-1} l_{i,k} u_{k,j},
\end{array}
```
mentre quelli della matrice $L$ sono dati da:
```{math}
:label: lu2
\forall\,i\quad & \begin{array}{ll}
j = 1, & l_{i,j} = \frac{a_{i,j}}{u_{j,j}}, \\
j > 1, & l_{i,j} = \displaystyle \frac{a_{i,j} - \sum_{k=1}^{j-1} l_{i,k} u_{k,j} }{u_{j,j}}.
\end{array}
```

:::::{admonition} Esercizio
Si usino le relazioni {eq}`lu1` e {eq}`lu2` per implementare la versione di
Doolittle dell'algoritmo $LU$ per una matrice quadrata. Un prototipo della
funzione è dato da:
```matlab
function [L,U] = doolittlelu(A)
%%DOOLITTLELU Calcola la fattorizzazione LU della matrice A con i vincolo di
%Doolittle.

% Controllo degli input:

% Pre-allocazione memoria:
L = eye(n,n);
U = zeros(n,n);
%!! Questa potrebbe essere anche:
% L = zeros(n,n);
% U = zeros(n,n);
% Dipende da come decidete di utilizzare le formule...

% Utilizzo delle formule per l_{ij} e u_{ij}.

end
```
- Si faccia un **controllo dell'input**: la matrice $A$ è quadrata?
- Si utilizzino operazioni vettorizzate per le somme che compaiono in {eq}`lu1`
e {eq}`lu2`

Per verificare la buona riuscita dell'implementazione si può usare il seguente
problema di test.
```matlab
%% Test della versione di Doolittle dell'algoritmo LU  

clear; clc;

A = [ 3 -1  4;
     -2  0  5;
      7  2 -2 ];
[L,U] = doolittlelu(A);
fprintf("|| A - LU || = %1.2e\n",norm(A - L*U,2));

% Vediamone gli ingressi:
format rat
disp(A)
disp(L)
disp(U)
disp(L*U)
format short
```

```{tip}
Guardando alla formulazione delle equazioni {eq}`lu1` e {eq}`lu2` potremmo
osservare che in realtà possiamo sovrascrivere le entrate della matrice $A$
con la fattorizzazione $LU$. Se avete concluso la prima parte pensate a come
è possibile realizzare questa versione.
```
::::{admonition} Un suggerimento più estensivo
:class: tip, dropdown
Possiamo sviluppare il ciclo di istruzioni sostanzialmente in due ordine, in un primo procediamo alla fattorizzazioni **colonna per colonna**, ovvero:
```matlab
for k=1:n
    L(k:n,k) = A(k:n,k) / A(k,k); % Dividiamo per il Pivot
    U(k,1:n) = A(k,1:n);          % Aggiorniamo le entrate della U
    A(k+1:n,1:n) = A(k+1:n,1:n) - L(k+1:n,k)*A(k,1:n);
end
U(:,end) = A(:,end);
```
Nel caso dell'**esempio** vediamo:
```{math}
k = 0, \qquad L = \left(\begin{array}{ccc} 0 & 0 & 0\\ 0 & 0 & 0\\ 0 & 0 & 0 \end{array}\right) \qquad U = \left(\begin{array}{ccc} 0 & 0 & 0\\ 0 & 0 & 0\\ 0 & 0 & 0 \end{array}\right), \qquad A = \left(\begin{array}{ccc} 3 & -1 & 4\\ -2 & 0 & 5\\ 7 & 2 & -2 \end{array}\right)
```
```{math}
k = 1, \qquad L^{(1)} = \left(\begin{array}{ccc} 1 & 0 & 0\\ -\frac{2}{3} & 0 & 0\\ \frac{7}{3} & 0 & 0 \end{array}\right), \qquad \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & 0 & 0\\ 0 & 0 & 0 \end{array}\right) \qquad A^{(1)} = \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & -\frac{2}{3} & \frac{23}{3}\\ 0 & \frac{13}{3} & -\frac{34}{3} \end{array}\right), \\
```
```{math}
k = 2, \qquad L^{(2)} = \left(\begin{array}{ccc} 1 & 0 & 0\\ -\frac{2}{3} & 1 & 0\\ \frac{7}{3} & -\frac{13}{2} & 0 \end{array}\right), \qquad U^{(2)} = \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & -\frac{2}{3} & \frac{23}{3}\\ 0 & 0 & 0 \end{array}\right), \qquad A^{(2)} = \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & -\frac{2}{3} & \frac{23}{3}\\ 0 & 0 & \frac{77}{2} \end{array}\right), \\
```
```{math}
k = 3, \qquad L^{(3)} = \left(\begin{array}{ccc} 1 & 0 & 0\\ -\frac{2}{3} & 1 & 0\\ \frac{7}{3} & -\frac{13}{2} & 1 \end{array}\right), \qquad U^{(3)} = \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & -\frac{2}{3} & \frac{23}{3}\\ 0 & 0 & \frac{77}{2} \end{array}\right), \qquad A^{(3)} = \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & -\frac{2}{3} & \frac{23}{3}\\ 0 & 0 & \frac{77}{2} \end{array}\right)
```

La **seconda alternativa** è quella di eseguire questo ciclo procedendo lungo le
righe della matrice
```matlab
for i = 1:1:n
    for j = 1:(i - 1)
        L(i,j) = (A(i,j) - L(i,1:(j - 1))*U(1:(j - 1),j)) / U(j,j);
    end
    j = i:n;
    U(i,j) = A(i,j) - L(i,1:(i - 1))*U(1:(i - 1),j);
end
```
Nel caso dell'**esempio** vediamo:
```{math}
i = 0, \qquad L = \left(\begin{array}{ccc} 1 & 0 & 0\\ 0 & 1 & 0\\ 0 & 0 & 1 \end{array}\right), \qquad U = \left(\begin{array}{ccc} 0 & 0 & 0\\ 0 & 0 & 0\\ 0 & 0 & 0 \end{array}\right), \qquad A = \left(\begin{array}{ccc} 3 & -1 & 4\\ -2 & 0 & 5\\ 7 & 2 & -2 \end{array}\right),
```
```{math}
i = 1, \qquad L^{(1)} = \left(\begin{array}{ccc} 1 & 0 & 0\\ 0 & 1 & 0\\ 0 & 0 & 1 \end{array}\right), \qquad U^{(1)} = \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & 0 & 0\\ 0 & 0 & 0 \end{array}\right),
```
```{math}
i = 2, \qquad L^{(2)} = \left(\begin{array}{ccc} 1 & 0 & 0\\ -\frac{2}{3} & 1 & 0\\ 0 & 0 & 1 \end{array}\right), \qquad U^{(2)} = \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & -\frac{2}{3} & \frac{23}{3}\\ 0 & 0 & 0 \end{array}\right)',
```
```{math}
i = 3, \qquad L^{(3)} = \left(\begin{array}{ccc} 1 & 0 & 0\\ -\frac{2}{3} & 1 & 0\\ \frac{7}{3} & -\frac{13}{2} & 1 \end{array}\right), \qquad U^{(3)} = \left(\begin{array}{ccc} 3 & -1 & 4\\ 0 & -\frac{2}{3} & \frac{23}{3}\\ 0 & 0 & \frac{77}{2} \end{array}\right)',
```
L'**unica differenza** tra le due implementazioni è l'ordine, per colonna e per riga rispettivamente, in cui facciamo le operazioni.
::::
:::::

(pivoting)=
### La fattorizzazione con *pivoting*

Per fare sì che la fattorizzazione
$LU$ si possa calcolare per una *generica matrice* (ovvero anche
  quando si incontrano *pivot* nulli) 
e per **migliorare la stabilità dell'algoritmo** ,
possiamo introdurre una strategia di **scambio delle righe** della matrice detta, per l'appunto, di *pivoting parziale*.

Formalmente, questo vuol dire che otterremo una fattorizzazione della forma
```{math}
P A = L U,
```
in cui $P$ è una **matrice di permutazione**. Strutturalmente una matrice di
permutazione non è nient'altro che una matrice identità di cui abbiamo permutato
le righe. Vediamo un *veloce esempio*:
```matlab
format rat
H = hilb(5)
P = eye(5)
P = P([5:-1:1],:)
P*H
```
che costruisce le matrici
```{math}
H = \left(\begin{array}{ccccc} 1 & \frac{1}{2} & \frac{1}{3} & \frac{1}{4} & \frac{1}{5}\\ \frac{1}{2} & \frac{1}{3} & \frac{1}{4} & \frac{1}{5} & \frac{1}{6}\\ \frac{1}{3} & \frac{1}{4} & \frac{1}{5} & \frac{1}{6} & \frac{1}{7}\\ \frac{1}{4} & \frac{1}{5} & \frac{1}{6} & \frac{1}{7} & \frac{1}{8}\\ \frac{1}{5} & \frac{1}{6} & \frac{1}{7} & \frac{1}{8} & \frac{1}{9} \end{array}\right), \quad P = \left(\begin{array}{ccccc} 0 & 0 & 0 & 0 & 1\\ 0 & 0 & 0 & 1 & 0\\ 0 & 0 & 1 & 0 & 0\\ 0 & 1 & 0 & 0 & 0\\ 1 & 0 & 0 & 0 & 0 \end{array}\right), \quad PH = \left(\begin{array}{ccccc} \frac{1}{5} & \frac{1}{6} & \frac{1}{7} & \frac{1}{8} & \frac{1}{9}\\ \frac{1}{4} & \frac{1}{5} & \frac{1}{6} & \frac{1}{7} & \frac{1}{8}\\ \frac{1}{3} & \frac{1}{4} & \frac{1}{5} & \frac{1}{6} & \frac{1}{7}\\ \frac{1}{2} & \frac{1}{3} & \frac{1}{4} & \frac{1}{5} & \frac{1}{6}\\ 1 & \frac{1}{2} & \frac{1}{3} & \frac{1}{4} & \frac{1}{5} \end{array}\right),
```
e ci mostra anche come possiamo **permutare** due (o più) righe di una matrice
MATLAB. Abbiamo ora tutti gli ingredienti necessari a costruire la nostra versione
fatta in casa della fattorizzazione LU con pivoting in MATLAB.

::::{admonition} Esercizio
Si costruisca una funzione `ludecomp` per costruire la fattorizzazione $LU$ con
pivoting parziale per una matrice $A$. Un prototipo della funzione è quindi:
```matlab
function [L,U,P] = ludecomp(A)
%%LUDECOMPP Fattorizzazione LU con pivoting parziale.

% Controllo dell'input:

% Allocazione della memoria:
L = zeros(n);
U = zeros(n);
P = eye(n);

% Fattorizzazione con Pivoting
end
```
- Si faccia un **controllo dell'input**: la matrice $A$ è quadrata?
- Si utilizzino operazioni vettorizzate per le somme che compaiono in {eq}`lu1`
e {eq}`lu2`.

:::{tip}
Se il ciclo esterno è per $k=1,\ldots,n$ possiamo **individuare il pivot** facendo:
```matlab
[~,r] = max(abs(A(k:end,k)));
r = n-(n-k+1)+r; % Trasliamo l'indice in questione
```
Qual è l'effetto della traslazione di indice?
:::

::::

<!---
Possiamo **verificare** l'implementazione della versione con *pivoting parziale*
facendone un paragone con la versione dell'esercizio precedente.


Confrontiamo gli errori assoluti e relativi,
```{math}
E_{\text{abs}} = \| A - LU\|_2, \qquad E_{\text{rel}} = \frac{\| A - LU\|_2}{\| A\|_2},
```
tra la matrice $A$ ed il prodotto dei termini $LU$ ottenuti con i diversi algoritmi (con e senza pivoting parziale).
```matlab
clear; clc; close all;
sizes = 100:100:1500;
N = length(sizes);
maxerr = zeros(N,1);
relerr = zeros(N,1);
maxerrm = zeros(N,1);
relerrm = zeros(N,1);
for i = 1:N
   A = rand(sizes(i),sizes(i));
   [L,U] = doolittlelu(A);
   [Lm,Um,P] = ludecomp(A);
   normA = norm(A,2);
   maxerr(i) = norm(A - L*U,2);
   relerr(i) = maxerr(i)/normA;
   maxerrm(i) = norm(P*A - Lm*Um,2);
   relerrm(i) = maxerrm(i)/normA;
   fprintf("Dimensione %d\n \t || A - LU || = %1.2e\n",sizes(i),maxerr(i));
   fprintf("\t || A - LU ||/||A|| = %1.2e\n",relerr(i));
   fprintf("\t || P*A - LU ||/||A|| = %1.2e\n",maxerrm(i));
   fprintf("\t || P*A - LU ||/||A|| = %1.2e\n",relerrm(i));
end

figure(2)
loglog(sizes,maxerr,'o-',sizes,relerr,'x-',...
    sizes,maxerrm,'o-',sizes,relerrm,'x-','LineWidth',2);
xlabel('Dimensione');
legend({'Errore Assoluto','Errore Relativo',...
    'Errore Assoluto (Pivoting)','Errore Relativo (Pivoting)'},...
    'Location','best','FontSize',14);
```

Quello che osserviamo nella figura è una crescita di più di un'ordine di
grandezza dell'errore assoluto nella fattorizzazione $LU$ senza pivoting, ed un
certo numero di oscillazioni {numref}`lusipivoting`.

Ci aspetteremmo che almeno l'**errore relativo** si comportasse meglio. Quello
che sta accadendo è che la scelta dei *pivot* affidata alla forma originale della
matrice sta introducendo dell'instabilità numerica che ci fa perdere in precisione
relativa più di quanto sarebbe lecito aspettarsi guardando ai numeri di
condizionamento delle operazioni in questione.

L'implementazione con pivoting, invece, permette di mantenere un errore
relativo che è all'incirca della precisione di macchina, migliorando
sensibilmente le performance rispetto alla versione che utilizzava
l'ordinamento di partenza.
```{figure} ./images/lusipivoting.png
:name:  lusipivoting

Errori assoluti e relativi tra una matrice $A$ ad entrate casuali
uniformemente distribuite in $[0,1]$ e la relativa fattorizzazione $LU$
calcolata **senza** e **con** pivoting.
```
 -->

### Calcolo del determinante

Dopo aver calcolato la fattorizzazione $LU$ di una matrice $A$ è possibile
utilizzarla anche per calcolare il determinante della matrice $A$. Si può fare
sfruttando il **teorema di Binet** che ricordiamo collega il prodotto fra
matrici quadrate con il determinante:
```{math}
\det(AB) = \det(A)\det(B), \qquad \forall\,A,B \in M_n(\mathbb{K}),
```
ed, in particolare,
```{math}
\det(A) = \det(LU) = \det(L)\det(U),
```
che, se usiamo la formulazione di Doolittle si riduce quindi ad essere:
```{math}
\det(A) = \det(LU) = \det(L)\det(U) = \prod_{i=1}^{n} l_{ii} \prod_{i=1}^n u_{ii} = 1 \prod_{i=1}^n u_{ii} = \prod_{i=1}^n u_{ii}.
```
Che possiamo riportare facilmente su MATLAB come:
```matlab
function res = ludet(A)
%%LUDET Calcola il determinante della matrice A sfruttando la sua fattorizzazione
%LU in forma di Doolittle e senza pivoting.
[~,U] = doolittlelu(A);
res = prod(diag(U));
end
```
e testare rapidamente con:
```matlab
A = pascal(5,2);
fprintf('Errore relativo: %1.2e\n',abs(ludet(A) - det(A))/abs(det(A)));
```
che ci dovrebbe restituire un errore relativo, rispetto alla funzione `det`
nativa di MATLAB dell'ordine di `3.33e-16`.

:::{warning}
Per costruire il calcolo del determinante anche nel caso in cui si utilizzi la fattorizzazione $LU$ con *pivoting* è necessario tenere traccia del numero di
scambi di righe effettuati.
:::


## Le funzioni di MATLAB

MATLAB ha implementata tra le sue routine una versione della fattorizzazione $LU$ ed un
algoritmo per la soluzione dei sistemi lineari.

Per la fattorizzazione LU il comando ha il medesimo nome `lu` e permette una certa
flessibilità.
- `[L,U] = lu(A)`, restituisce una matrice triangolare superiore in $U$ e
una matrice triangolare inferiore **permutata** in $L$, tale che $A = LU$.
La matrice di input $A$ può essere *densa* o *sparsa*.
- `[L,U,P] = lu(A)`, restituisce una matrice triangolare inferiore $L$ con
diagonale di uno, matrice triangolare superiore $U$ ed una matrice di
permutazione $P$ tale che $PA = LU.$ Con un argomento di input *opzionale* è
possibile scegliere il formato della $P$, ovvero se chiama `lu(A, outputForm)`
con `outputForm` uguale a `'matrix'` (default) $P$ è memorizzata come una matrice,
se si chiama con l'opzione `'vector'` allora $P$ è un vettore tale che `A(P,:) = L*U`.
- `[L,U,P,Q] = lu(A)` restituisce una matrice triangolare inferiore $L$ con
diagonale di uno, matrice triangolare superiore $U$ e due matrici di
permutazione $P$ e $Q$ tali che $PAQ = LU$. Per una matrice sparsa questa opzione
è significativamente più efficiente della versione con tre output.

```{prf:definition}
Una **matrice sparsa** o **array sparso** è una matrice in cui la maggior parte
degli elementi è zero. Per dare una definizione rigorosa riguardo alla
proporzione di elementi di valore zero affinché una matrice sia qualificata
come sparsa è tipicamente necessario considerare la sorgente del problema che la
ha generata, e.g., la discretizzazione di una PDE. Un **criterio comune** è che
il numero di elementi diversi da zero sia approssimativamente uguale al numero
di righe o colonne.

Al contrario, se la maggior parte degli elementi è diversa da zero, la matrice è
considerata densa.
```

Informazioni sulle matrici sparse possono essere recuperate dal manuale di
MATLAB (`help sparse`). Non scenderemo qui in ulteriori dettagli.

La funzione per **risolvere i sistemi lineari** è la funzione `mldivide` `\`.
Si tratta di una delle funzioni più elaborate di MATLAB, il comando `A\B` è la
*divisione matriciale* di `A` in `B`, moralmente è la stessa operazione di
`INV(A)*B`, ma è **calcolata in modo diverso**, **stabile** ed **efficiente**.

Se $A$ è una matrice $N$ per $N$ e $\mathbf{b}$ è un vettore colonna con $N$
componenti, o una matrice con più di queste colonne, allora `X = A\b` è la
soluzione dell'equazione `A*X = B`.

:::{warning}
Viene stampato un messaggio di avviso se A è scalata male o malcondizionata, ad
**esempio**
```matlab
x = hilb(15)\ones(15,1);
```
restituisce
```
Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
RCOND =  5.460912e-19.
```
:::

## Esercizi

Consideriamo alcune applicazioni *ingegneristiche* del problema calcolare le
soluzione di un sistema {cite}`kiusalaas2015`.

:::{admonition} Esercizio
La formulazione per gli spostamenti di una trave reticolare piana è simile a
quella di un sistema di masse e molle. Le differenze sono che
1. i termini di rigidità sono dati da $k_{i} = (E A/L)_ {i}$, dove $E$ è il modulo
di elasticità, $A$ rappresenta l'area della sezione di taglio e $L$ è la
lunghezza dell'elemento;
2. ci sono due componenti dello scostamento in ogni punto di giuntura.

Per la trave reticolare in figura:
```{figure} ./images/planetruss.png

Profilo di un elemento della trave reticolare piana.
```

si ottiene quindi il seguente sistema di equazioni lineari
```{math}
K \mathbf{u} = \mathbf{p}
```
dove
```matlab
K = [27.58  7.004 -7.004 0 0
     7.004 29.57  -5.253 0 -24.32
    -7.004 -5.253 29.57  0 0
    0 0 0 27.58 -7.004
    0 -24.32 0 -7.004 29.57];
p = [0 0 0 0 -45]'; % kN
```
si determinino tutti gli scostamenti $u_i$ per il sistema.
:::

:::{admonition} Esercizio
Riprendiamo l'esperimento fatto ne {ref}`pivoting` con la **matrice di Hilbert**
```{math}
(H_n)_{i,j} = \int_{0}^{1} x^{i+j-2}\,{\rm d}x = \frac{1}{i+j-1}, \qquad i,j=1,\ldots,n
```
che è un noto esempio di matrice *malcondizionata*. Si scriva un programma
specializzato nella risoluzione delle equazioni $H \mathbf{x} = \mathbf{b}$
mediante il metodo di decomposizione $LU$ *con* e *senza* pivoting, dove $H$ è
la matrice di Hilbert di dimensione arbitraria $n \times n$, e
```{math}
b_i = \sum_{j=1}^{n} (H_n)_{i,j},
```
L'unico input del programma deve essere $n$. Eseguendo il programma, si determinino
i più grandi valori di $n$ per i due algoritmi per cui la soluzione è entro
6 cifre significative della soluzione esatta.
:::

:::{margin} Sistema di molle
```{figure} ./images/springlsys.png
:name: springlsys

Sistema lineare di molle.
```
:::
:::{admonition} Esercizio
Il sistema mostrato in {numref}`springlsys` è costituito da $n$ molle lineari che supportano
$n$ masse. Le costanti elastiche delle molle sono indicate con $k_i$, i pesi
delle masse sono $W_i$, e $x_i$ sono gli spostamenti delle masse (misurati dalle
posizioni in cui le molle sono indeformate). La cosiddetta formulazione dello
spostamento si ottiene scrivendo l'equazione di equilibrio di ciascuna massa e
sostituendo le forze elastiche con $F_i = k_i (x_{i+1} - x_{i} )$.

Il risultato è il sistema simmetrico e tridiagonale di equazioni
```{math}
\begin{split}
(k_1 + k_2) x_1 - k_2 x_2 = & W_1,\\
- k_i x_{i-1} + (k_i + k_{i+1}) x_i - k_{i+1} x_{i+1} = & W_i, \qquad i=2,3,\ldots,n-1,\\
-k_n x_{n-1} + k_n x_n = & W_n,
\end{split}
```
si scriva un programma che risolve queste equazioni per dati valori di $n$, $\mathbf{k}$ e $\mathbf{W}$. Si risolva il problema con $n=5$ e
```{math}
\begin{split}
k_1 = k_2 = k_3 = 10\, \text{N}/\text{mm}\;\; k_4 = k_5 = 5\, \text{N}/\text{mm}\\ W_1 = W_3 = W_5 = 100\,\text{N},\qquad W_2 = W_4 = 50\,\text{N}.
\end{split}
```
:::


## Bibliografia

 ```{bibliography}
 :filter: docname in docnames
 ```
 
 <!---
 -->
