# wolPredictor
predict distribution of endosymbiotic bacteria among insects

CURRENTLY BEING EDITED

Methodological outline of wolPredictor.py 

Outline of Problem
It has been regularly suggested that highly prevalent Wolbachia induced reproductive isolation among arthropods is randomly distributed. If so, this implies that much arthropod biodiversity is a result of stochastically determined diversification events rather than process driven outcomes. For most arthropods we have limited knowledge about ecological contact which provides direct opportunity for horizontal exchange of microbes or genetic material between species. Here, I use Python programming to model our proposed mechanism that incorporates ecological contact in a phylogenetic context. The program runs quite quickly but is not yet optimised to take advantage of ‘numpy’ vectorisation potential.

Data requirements
(i) A CSV file with data columns featuring individual arthropod sample names, empirically derived Wolbachia strains (or absence of strain) for each sample, and host community name (here as host fig species).
(ii) A phylogenetic tree in nexus format. NB The program calls R (cophen4py.R – using the ‘ape’ library) using the Python library ‘subprocess’. If you cannot configure R (instructions for this will the added in the near future) to interact with Python you can comment out line 52 and use the ‘waspRaxml.tre_cophen.csv’ file for a test run (for your own analyses you would have to create a distance matrix of co-phenetic phylogenetic distances as formatted in the CSV file – see cophen4py.R for details – ordered alphabetically according to sample).

Functionality of wolPredictorXtra2_Kmeans_relaxed_numpy_mutliStrain.py
1. Initial step (main function)
Set required parameters that will determine the range of thresholds over which the program will iterate and delimit species, and to load in data files. All parameters have comments explaining their significance.

2. Species delimitation
In order to create putative species groupings, the program reformats the pairwise co-phenetic distance matrix (generated by the cophenetic R function in the ‘ape’ library) according to incremental thresholds set at the beginning of the run. For example, if the maximum pairwise co-phenetic distance from our data is 0.7 we may set ‘lower’ and ‘upper’ variables to 0 and 70, respectively, with the ‘increment’ variable set to 100 to give percentile increments between these values (alternatively, we may call 0, 700 & 1000 for these to increase the resolution that the program iterates over). If the maximum pairwise co-phenetic distance is unknown we can set lower and upper to 95 and 100, respectively.

	At a 1% threshold the program will reformulate the pairwise-distance matrix into 0s and 1s depending whether distances are less than or greater than 1%. At each incremental iteration (e.g. 2%, 3% etc) the matrix will be reformulated accordingly and use K-means clustering (functions: calcMat & calcK) to delimit species groupings. It should be emphasised that this method does not represent a systematic investigation of species delimitation space across different thresholds, but rather is a method to induce different random species delimitation hypotheses, contingent on phylogenetic relationships, that potentially represent either contemporary or historical clusterings of lineages through evolutionary time. The parameters ‘min_nSpp’ and ‘max_nSpp’ are set initially to place intuitively sensible limits on the number of species groups at each iteration. Other parameters are explained in comments.

3. Predict Wolbachia strain associations (function: addPredict)

The first step is to assign Wolbachia infection statuses. Under our working hypothesis, we expect that closely related, or contemporarily diverging, species in close ecological contact will host different Wolbachia strains that mediate reproductive isolation. Thus, wolPredictor assigns arbitrarily named strains to putative species groupings that co-occur on the same fig host. For example, if community #1 has been allocated three putative species groupings, the program will assign the strain names ‘com1_w1’, ‘com1_w2’ and ‘com1_w3’ to all individuals in the three putative clusters. However, if community #2 has only been allocated a single putative species grouping it will be assigned the strain name ‘noWol’ (i.e. no Wolbachia) as there is no additional putative species in the community from which it may require Wolbachia to mediate reproductive isolation. 

4. Remove Wolbachia strains according to co-phenetic distance between species clusters (function: wolPurger)

As Wolbachia has been shown to typically drop-out from host lineages after approx. 7 million years, our working hypothesis includes the possibility that species clusters separated by sufficient evolutionary time will have undergone Wolbachia purging. wolPredictor uses pairwise co-phenetic data to incrementally remove Wolbachia according to cut-off thresholds between putative species clusters within a community. Species clusters that are separated by a distance greater that the incremental threshold will have their assigned Wolbachia strain converted to ‘noWol’. This is done conservatively, if there are three species within a community and the evolutionary distances between two of them are greater than the threshold, Wolbachia will not be purged if their evolutionary distances from the third species are below the threshold. Calculations are made for all threshold purging cut-offs (using the parameter ‘purge’, which sets the maximum purging distance iterated to from zero, at intervals set by ‘pge_incr’). At this point, all Wolbachia strain predictions at every species delimitation iteration are written to a CSV file with a root name ‘wolPreds_threshClades_incr’; this includes calculations for species delimitation iterations without any purging having been performed.

5. Match arbitrarily assigned Wolbachia strain names to empirically derived strain names (function: matchStrains)

In order to calculate the accuracy of our model at each iteration, we need to match arbitrarily assigned Wolbachia strain names to the empirically derived strain names, as much as is possible. Again, this is done conservatively, if a community features two putative species delimitation clusters whose individuals all match up to a single Wolbachia strain, only one species will have their arbitrarily assigned Wolbachia strain names converted to match the empirical predictions. For example, if community #1 has two species clusters (e.g. ‘com1_w1’ and ‘com1_w2’) whose individuals have been shown to have the empirical strain ‘wspC6’ (see our data), then only the species cluster with the most individuals will be reassigned to the strain name ‘wspC6’. The second putative species will retain the arbitrarily assigned strain name (i.e. ‘com1_w1’ or ‘com1_w2’). 

6. outputInvestigator script

The program ‘outputInvestigator.py’ is a short script to investigate performance from the output files of the main program. 

run as: python outputInvestigator.py filename int1 int2

The first argument is the full name of the output file with the root name: correctedWolPreds; the 2nd and 3rd arguments are the minimum number of negative and positive Wolbachia strains that you require highlighting. The output will indicate which runs of the ‘wolPredictor’ program scored higher than the thresholds set here. There maybe a trade-off between negative and positive strain assignment by wolPredictor so it advisable to play with these two thresholds – starting low then increasing to highlight the best runs. From the output, individual columns of the output CSVs can be examined to see how predictions were made.

7. Output results (makePDFv1.py & makePDFv2.py)

run as: python makePDFvX.py filename[with “correctedWolPreds” root] 

We use the Python package matplotlib to output a graphical representation of predictive accuracy performed at each iteration. This will include the strains assigned as ‘noWol’ indicating no Wolbachia infection. For our dataset, ‘noWol’ accounts for about 70% of the dataset. Thus, by uniformly predicting ‘noWol’ for all individuals it is possible to gain an impressively predictive score of 70%. To guard against this outcome, the program makePDFv2.py (alongside makePDFv1.py which assesses negative strains) only assesses predictive power for individuals with positive empirical Wolbachia associations. wolPredictor may over-estimate ‘noWol’ assessments under certain conditions. For example, if virtually all individuals from the dataset are assigned to the same putative species cluster (this may occur at certain threshold cut-offs) – under such a scenario Wolbachia would not be required as multiple species are not co-occurring within communities. The ‘min_nSpp’ parameter is useful for guarding against unrealistically low species diversity assessments. The programs make a graphical PDF output and write CSVs of all runs’ predictive outputs.
