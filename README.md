# pacs-elasticityequation
PACS PROJECT : ELASTICITY EQUATION by Alessandra Arrigoni and Sara Francesca Pichierri.

TODO (in base ai commenti di Formaggia)
1) Cambiare typedef con using in Core.h
2) usare smart pointers e inizializzarli con make_unique
3) MAI mettere come return type di una funzione una reference a un oggetto creato all'interno; posso però farlo se la ref è a uno dei paramentri in input, passati anche loro per reference, e la funzione serve per modificarli.
4) per passare parametri alle funzioni, per distinguere quelli che voglio modificare da quelli che devono rimanere costanti uso le references o le const references così è chiaro fin da subito.
5) metodi che non cambiano lo stato della classe vanno dichiarati const (il termine si riferisce sempre al "nome" che è alla sua sinistra) anche per avere un codice ottimizzato!

REFACTORING (idee possibili)
1) fare una classe Problem astratta da cui ereditano le 2 classi laplaciano e elasticity 
- la funzione assembleMatrix risulterebbe virtuale e va ridefinita nelle due classi derivate
- la funzione enforceStrongBC va ridefinita tenendo conto di Qdim anche per il laplaciano, ma poi non dovrei avere problemi perchè in quel caso sarebbe 1 quindi i loop su j si fermano a j==0
- stessa cosa per la funzione enforceInterfaceJump, exactSolution e bulkLoad (operatorsBulk.cc), costruttore del FEM (riceve anche qdim, spero che poi set_qdim funzioni anche nel caso qdim ==1), stressRHS (operatorsBD) però poi qui devo cambiare anche BCNeum nel file BC... 
!! In tutti i casi è meglio tenere le funzioni della cartella elasticità (contengono già Qdim), invece per le BC è meglio tenere quelle del laplaciano e ridefinirle in modo che vadano bene anche per un vettore di stringhe come già avviene per la forzante! Devo anche cambiare il modo di definire le condizioni nel file e lasciare una sola variabile M_BCNeum e non X e Y.

2) anche per il BulkData (gestisce la forzante, la soluzione esatta e i parametri collegati al dominio) dovrei definire classi diverse o trovare un modo intelligente di gestire i parametri: nel laplaciano ne ho uno solo, nell'elasticità 2 però vanno gestiti nello stesso identico modo! --> posso fare una classe che si chiama semplicemente bulkDatum che può essere o un coefficiente o una funzone e poi definire 3,4,5 oggetti di questa classe nel problema derivato... questa potrebbe essere l'idea migliore!
NB: nel constructor del problem ho già l'inizializzazione di tutti i vettori dei coefficienti lambda, mu, diff in bulkdata, quindi poi quando costruisco le matrici in OperatorsBulk non è più necessario creare i vettori locali nodo per nodo, ma basta recuperare quelli!

3) Nel LinearSystem ho lasciato tutto con i puntatori (std::shared_pointers) come nel codice originario mentre nelle altre classi cerco di mettere references
