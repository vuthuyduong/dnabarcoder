#!/bin/bash
cd /data/ProgLang/Python/dnabarcoder/data/CBSITS2
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank species -o CBSITS2.species.fasta
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank genus -o CBSITS2.genus.fasta
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank family -o CBSITS2.family.fasta
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank order -o CBSITS2.order.fasta
../../aidscripts/selectsequences.py -i CBSITS2.fasta -c ITS_20211006.classification -rank class -o CBSITS2.class.fasta
../../dnabarcoder.py sim -i CBSITS2.fasta -ml 50
../../dnabarcoder.py variation -i CBSITS2.fasta -c ITS_20211006.classification -ranks class,order,family,genus,species -sim dnabarcoder/CBSITS2.sim
../../dnabarcoder.py remove -i CBSITS2.species.fasta -c ITS_20211006.classification -rank species -sim dnabarcoder/CBSITS2.sim
../../dnabarcoder.py predict -i dnabarcoder/CBSITS2.species.diff.fasta -c ITS_20211006.classification -st 0.9 -et 1 -s 0.001 -ranks species -prefix CBSITS2 -sim dnabarcoder/CBSITS2.sim
../../dnabarcoder.py predict -i dnabarcoder/CBSITS2.species.diff.fasta -c ITS_20211006.classification -st 0.7 -et 1 -s 0.001 -ranks species -prefix CBSITS2 -higherranks genus,family,order,class,phylum -sim dnabarcoder/CBSITS2.sim -minSeqNo 30
../../dnabarcoder.py predict -i CBSITS2.genus.fasta -c ITS_20211006.classification -st 0.7 -et 1 -s 0.001 -ranks genus -prefix CBSITS2
../../dnabarcoder.py predict -i CBSITS2.genus.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -ranks genus -higherranks family,order,class,phylum -prefix CBSITS2 -minSeqNo 30
../../dnabarcoder.py predict -i CBSITS2.family.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -ranks family -prefix CBSITS2
../../dnabarcoder.py predict -i CBSITS2.family.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -ranks family -higherranks order,class,phylum -prefix CBSITS2 -minSeqNo 30
../../dnabarcoder.py predict -i CBSITS2.order.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -ranks order -prefix CBSITS2
../../dnabarcoder.py predict -i CBSITS2.order.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -ranks order -higherranks class,phylum -prefix CBSITS2 -minSeqNo 30
../../dnabarcoder.py predict -i CBSITS2.class.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -ranks class -prefix CBSITS2
../../dnabarcoder.py predict -i CBSITS2.class.fasta -c ITS_20211006.classification -st 0.5 -et 1 -s 0.001 -ranks class  -higherranks phylum -prefix CBSITS2 -minSeqNo 30
#visualization 
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -ranks species,genus,family,order,class 
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -ranks species -higherranks genus -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -ranks genus -higherranks family -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -ranks family -higherranks order -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -ranks order -higherranks class -minseqno 30
../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -ranks class -higherranks phylum  -minseqno 30
../../dnabarcoder.py visualize -i CBSITS2.fasta  -c ITS_20211006.classification -p 3 -n 10 -size 0.3
../../dnabarcoder.py best -i dnabarcoder/CBSITS2.cutoffs.json -c ITS_20211006.classification

