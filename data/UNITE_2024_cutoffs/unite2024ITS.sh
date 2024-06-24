#!/bin/bash
cd /home/dvu/data/volume_2/cephstorage/ProgLang/Python/dnabarcoder/alldata/UNITE_ITS_2024

#select unique sequences
../../aidscripts/selectsequences.py -i unite2024ITS.fasta -unique yes -c unite2024ITS.classification -o unite2024ITS.unique.fasta
../../dnabarcoder.py length -i unite2024ITS.unique.fasta -l 100
../../dnabarcoder.py overview -i unite2024ITS.unique.fasta -c unite2024ITS.unique.classification

#species 
../../aidscripts/selectsequences.py -i unite2024ITS.unique.fasta -c unite2024ITS.unique.classification -rank species -o unite2024ITS.unique.species.fasta 
../../dnabarcoder.py predict -i unite2024ITS.unique.species.fasta -c unite2024ITS.unique.species.classification -st 0.7 -et 1 -s 0.001 -rank species -prefix unite2024ITS.unique -higherrank genus -maxproportion 0.9 -removecomplexes yes
../../dnabarcoder.py predict -i unite2024ITS.unique.species.fasta -c unite2024ITS.unique.species.classification -st 0.9 -et 1 -s 0.001 -rank species -prefix unite2024ITS.unique -removecomplexes yes
rm unite2024ITS.unique.species.*

#genus
../../aidscripts/selectsequences.py -i unite2024ITS.unique.fasta -c unite2024ITS.unique.classification -rank genus -maxseqnopergroup 1000 -o unite2024ITS.unique.genus.fasta  
../../dnabarcoder.py predict -i unite2024ITS.unique.genus.fasta -c unite2024ITS.unique.genus.classification -st 0.7 -et 1 -s 0.001 -rank genus -higherrank family -prefix unite2024ITS.unique -maxproportion 0.75 
../../dnabarcoder.py predict -i unite2024ITS.unique.genus.fasta -c unite2024ITS.unique.genus.classification -st 0.7 -et 1 -s 0.001 -rank genus -prefix unite2024ITS.unique 
rm unite2024ITS.unique.genus.*

#family
../../aidscripts/selectsequences.py -i unite2024ITS.unique.fasta -c unite2024ITS.unique.classification -rank family -maxseqnopergroup 1000 -o unite2024ITS.unique.family.fasta
../../dnabarcoder.py predict -i unite2024ITS.unique.family.fasta -c unite2024ITS.unique.family.classification -st 0.5 -et 1 -s 0.001 -rank family -higherrank order -prefix unite2024ITS.unique -maxproportion 0.75 
../../dnabarcoder.py predict -i unite2024ITS.unique.family.fasta -c unite2024ITS.unique.family.classification -st 0.5 -et 1 -s 0.001 -rank family -prefix unite2024ITS.unique
rm unite2024ITS.unique.family.*

#order
../../aidscripts/selectsequences.py -i unite2024ITS.unique.fasta -c unite2024ITS.unique.classification -rank order -maxseqnopergroup 1000 -o unite2024ITS.unique.order.fasta 
../../dnabarcoder.py predict -i unite2024ITS.unique.order.fasta -c unite2024ITS.unique.order.classification -st 0.5 -et 1 -s 0.001 -rank order -higherrank class -prefix unite2024ITS.unique -maxproportion 0.75
../../dnabarcoder.py predict -i unite2024ITS.unique.order.fasta -c unite2024ITS.unique.order.classification -st 0.5 -et 1 -s 0.001 -rank order -prefix unite2024ITS.unique
rm unite2024ITS.unique.order.*

#class
../../aidscripts/selectsequences.py -i unite2024ITS.unique.fasta -c unite2024ITS.unique.classification -rank class -maxseqnopergroup 1000 -o unite2024ITS.unique.class.fasta 
../../dnabarcoder.py predict -i unite2024ITS.unique.class.fasta -c unite2024ITS.unique.class.classification -st 0.5 -et 1 -s 0.001 -rank class -higherrank phylum -prefix unite2024ITS.unique -maxproportion 0.75
../../dnabarcoder.py predict -i unite2024ITS.unique.class.fasta -c unite2024ITS.unique.class.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix unite2024ITS.unique
rm unite2024ITS.unique.class.*

#visualization
../../dnabarcoder.py predict -i unite2024ITS1.unique.fasta -c unite2024ITS1.unique.classification -rank species,genus,family,order,class 

#compute best cutoffs
../../dnabarcoder.py best -i dnabarcoder/unite2024ITS.unique.cutoffs.json -c unite2024ITS.unique.classification -mincutoff 0.71
