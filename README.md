# pacs-elasticityequation
PACS PROJECT : ELASTICITY EQUATION by Alessandra Arrigoni and Sara Francesca Pichierri.

TODO (in base ai commenti di Formaggia)
1) Cambiare typedef con using in Core.h
2) usare smart pointers e inizializzarli con make_unique
3) MAI mettere come return type di una funzione una reference a un oggetto creato all'interno; posso però farlo se la ref è a uno dei paramentri in input, passati anche loro per reference, e la funzione serve per modificarli.
4) per passare parametri alle funzioni, per distinguere quelli che voglio modificare da quelli che devono rimanere costanti uso le references o le const references così è chiaro fin da subito.
5) metodi che non cambiano lo stato della classe vanno dichiarati const (il termine si riferisce sempre al "nome" che è alla sua sinistra) anche per avere un codice ottimizzato!
