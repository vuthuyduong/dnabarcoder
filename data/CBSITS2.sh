#!/bin/bash
cd /data/ProgLang/Python/dnabarcoder/data/CBSITS2
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank species -o CBSITS2.species.fasta
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank genus -o CBSITS2.genus.fasta
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank family -o CBSITS2.family.fasta
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank order -o CBSITS2.order.fasta
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank class -o CBSITS2.class.fasta
../../dnabarcoder.py sim -i CBSITS2.fasta -ml 50
../../dnabarcoder.py variation -i CBSITS2.fasta -c ITS_20211006.classification -rank class,order,family,genus,species -sim dnabarcoder/CBSITS2.sim
../../dnabarcoder.py remove -i CBSITS2.species.fasta -c ITS_20211006.classification -rank species -sim dnabarcoder/CBSITS2.sim
../../dnabarcoder.py predict -i dnabarcoder/CBSITS2.species.diff.fasta -c ITS_20211006.classification -st 0.9 -et 1 -s 0.001 -rank species -prefix CBSITS2 -sim dnabarcoder/CBSITS2.sim
../../dnabarcoder.py predict -i dnabarcoder/CBSITS2.species.diff.fasta -c ITS_20211006.classification -st 0.7 -et 1 -s 0.001 -rank species -prefix CBSITS2 -higherrank genus,family,order,class,phylum -sim dnabarcoder/CBSITS2.sim -minSeqNo 30
../../dnabarcoder.py predict -i CBSITS2.genus.fasta -c ITS_20211006.classification -st 0.7 -et 1 -s 0.001 -rank genus -prefix CBSITS2
../../dnabarcoder.py predict -i CBSITS2.genus.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -rank genus -higherrank family,order,class,phylum -prefix CBSITS2 -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.family.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -rank family -prefix CBSITS2
../../dnabarcoder.py predict -i CBSITS2.family.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -rank family -higherrank order,class,phylum -prefix CBSITS2 -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.order.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -rank order -prefix CBSITS2
../../dnabarcoder.py predict -i CBSITS2.order.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -rank order -higherrank class,phylum -prefix CBSITS2 -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.class.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix CBSITS2
../../dnabarcoder.py predict -i CBSITS2.class.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -rank class  -higherrank phylum -prefix CBSITS2 -minseqno 30
#visualization 
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -rank species,genus,family,order,class 
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -rank species -higherrank genus -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -rank genus -higherrank family -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -rank family -higherrank order -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -rank order -higherrank class -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -rank class -higherrank phylum  -minseqno 30
../../dnabarcoder.py visualize -i CBSITS2.fasta  -c ITS_20211006.classification -p 3 -n 10 -size 0.3
../../dnabarcoder.py best -i dnabarcoder/CBSITS2.cutoffs.json -c ITS_20211006.classification

