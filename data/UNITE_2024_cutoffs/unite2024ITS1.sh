#!/bin/bash
#change the path to dnabarcoder and the UNITE data if needed
cd /home/dvu/data/volume_2/cephstorage/ProgLang/Python/dnabarcoder/alldata/UNITE_ITS1_2024
../../aidscripts/selectsequences.py -i unite2024ITS1.fasta -unique yes -c unite2024ITS1.classification -o unite2024ITS1.unique.fasta
../../dnabarcoder.py length -i unite2024ITS1.unique.fasta -l 50
../../dnabarcoder.py overview -i unite2024ITS1.unique.fasta -c unite2024ITS1.unique.classification

#species 
../../aidscripts/selectsequences.py -i unite2024ITS1.unique.fasta -c unite2024ITS1.unique.classification -rank species -o unite2024ITS1.unique.species.fasta 
../../dnabarcoder.py predict -i unite2024ITS1.unique.species.fasta -c unite2024ITS1.unique.species.classification -st 0.7 -et 1 -s 0.001 -rank species -prefix unite2024ITS1.unique -higherrank genus -maxproportion 0.9 -removecomplexes yes -ml 50
../../dnabarcoder.py predict -i unite2024ITS1.unique.species.fasta -c unite2024ITS1.unique.species.classification -st 0.9 -et 1 -s 0.001 -rank species -prefix unite2024ITS1.unique -removecomplexes yes -ml 50
rm unite2024ITS1.unique.species.*

#genus
../../aidscripts/selectsequences.py -i unite2024ITS1.unique.fasta -c unite2024ITS1.unique.classification -rank genus -maxseqnopergroup 1000 -o unite2024ITS1.unique.genus.fasta  
../../dnabarcoder.py predict -i unite2024ITS1.unique.genus.fasta -c unite2024ITS1.unique.genus.classification -st 0.7 -et 1 -s 0.001 -rank genus -higherrank family -prefix unite2024ITS1.unique -maxproportion 0.75 -ml 50
../../dnabarcoder.py predict -i unite2024ITS1.unique.genus.fasta -c unite2024ITS1.unique.genus.classification -st 0.7 -et 1 -s 0.001 -rank genus -prefix unite2024ITS1.unique -ml 50
rm unite2024ITS1.unique.genus.*

#family
../../aidscripts/selectsequences.py -i unite2024ITS1.unique.fasta -c unite2024ITS1.unique.classification -rank family -maxseqnopergroup 1000 -o unite2024ITS1.unique.family.fasta
../../dnabarcoder.py predict -i unite2024ITS1.unique.family.fasta -c unite2024ITS1.unique.family.classification -st 0.5 -et 1 -s 0.001 -rank family -higherrank order -prefix unite2024ITS1.unique -maxproportion 0.75 -ml 50 
../../dnabarcoder.py predict -i unite2024ITS1.unique.family.fasta -c unite2024ITS1.unique.family.classification -st 0.5 -et 1 -s 0.001 -rank family -prefix unite2024ITS1.unique -ml 50
rm unite2024ITS1.unique.family.*

#order
../../aidscripts/selectsequences.py -i unite2024ITS1.unique.fasta -c unite2024ITS1.unique.classification -rank order -maxseqnopergroup 1000 -o unite2024ITS1.unique.order.fasta 
../../dnabarcoder.py predict -i unite2024ITS1.unique.order.fasta -c unite2024ITS1.unique.order.classification -st 0.5 -et 1 -s 0.001 -rank order -higherrank class -prefix unite2024ITS1.unique -maxproportion 0.75 -ml 50
../../dnabarcoder.py predict -i unite2024ITS1.unique.order.fasta -c unite2024ITS1.unique.order.classification -st 0.5 -et 1 -s 0.001 -rank order -prefix unite2024ITS1.unique -ml 50
rm unite2024ITS1.unique.order.*

#class
../../aidscripts/selectsequences.py -i unite2024ITS1.unique.fasta -c unite2024ITS1.unique.classification -rank class -maxseqnopergroup 1000 -o unite2024ITS1.unique.class.fasta 
../../dnabarcoder.py predict -i unite2024ITS1.unique.class.fasta -c unite2024ITS1.unique.class.classification -st 0.5 -et 1 -s 0.001 -rank class -higherrank phylum -prefix unite2024ITS1.unique -maxproportion 0.75 -ml 50
../../dnabarcoder.py predict -i unite2024ITS1.unique.class.fasta -c unite2024ITS1.unique.class.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix unite2024ITS1.unique -ml 50
rm unite2024ITS1.unique.class.*

#visualization
../../dnabarcoder.py predict -i unite2024ITS1.unique.fasta -c unite2024ITS1.unique.classification -rank species,genus,family,order,class 

#compute the best cutoffs
../../dnabarcoder.py best -i dnabarcoder/unite2024ITS1.unique.cutoffs.json -c unite2024ITS1.unique.classification -mincutoff 0.71
