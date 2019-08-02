# wolPredictor
Predict distribution of endosymbiotic bacteria among insects

Methodological outline of wolPredictor.py 

Outline of Problem:
It has been regularly suggested that highly prevalent Wolbachia induced reproductive isolation among arthropods is randomly distributed. If so, this implies that much arthropod biodiversity is a result of stochastically determined diversification events rather than process driven outcomes. For most arthropods we have limited knowledge about ecological contact which provides direct opportunity for horizontal exchange of microbes or genetic material between species. Here, I use Python programming to model our proposed mechanism that incorporates ecological contact in a phylogenetic context. 

Providing all correct Python and R libraries are installed (see wolPredictor_xmeansDelim.py, makePDF.py & cophen4pyOut.R), the program should run directly from the unzipped download bundle. The program runs quite quickly (a few minutes at the featured parameter settings) but has not yet been optimised for vectorisation or parallelisation.

Data requirements:
(i) A CSV file with data columns featuring individual arthropod sample names, empirically derived Wolbachia strains (or absence of strain) for each sample (IMPORTANT: absence of Wolbachia should be coded as 'noWol'), and host community name (here labelled "community" - IMPORTANT: the first three characters of each community's name must be a unique combination of characters although the actual name can be longer than three characters).
(ii) A phylogenetic tree in nexus format. NB The program calls R (cophen4py.R – using the ‘ape’ library) using the Python library ‘subprocess’. If you cannot configure R (instructions for this will the added in the near future) to interact with Python you must use [NB VERSION to be added] with the ‘testPhylo.tre_cophen.csv’ file for a test run (for your own analyses without configuring R you will simply have to create a distance matrix of co-phenetic phylogenetic distances as formatted in the CSV file – see cophen4py.R for the code in R create the file, NB replace the last line of R code with:  `write.csv(as.matrix(phydist), quote = F, row.names = F) `).

How to run “wolPredictor_xmeansDelim.py”
The program can be run from a Python3 console (e.g. Spyder/Anaconda). Typing:
`python wolPredictor_xmeansDelim.py -h`
will bring up a menu of parameter options. Most are straightforward and relate to data input option (filenames & directories etc). However, three key considerations are:
`-m` min_nSpp: The minimum number of putative species clusters that Xmeans will divide the phylogenetic data into [default = 1]
`-M` max_nSpp: The maximum number of putative species clusters that Xmeans will divide the phylogenetic data into [default = 50] 
`-i` nSpp_incr: The incremental increase between min_nSpp and max_nSpp. NB must be ≥2 [default = 2]
-m & -M should be reasonable guestimates of minimum and maximum species richness from your tree. -i should be a reasonable increment between the two values (so that the runtime is not overlong: progress print-outs are outputted so you can judge if a run will take too long – change -i if necessary).
So, to run the program you might type:
`python wolPredictor_xmeansDelim.py -d data_directory -o out_directory -m 5 -M 100 -i 5`
You can also run the program from a Python IDE and alter the above parameters in the code.

Explanation of functions
1. R_cophen: call R (‘ape’ library) to create cophenetic distance matrix
2. calcX: X-means species delimitation:
In order to create putative species groupings, the program uses the pairwise co-phenetic distance matrix and X-means clustering to delimit species groupings. At each loop iteration, WOLPREDICTOR will calculate clusters of species for each value of species richness (i.e. between -m & -M). It should be emphasised that X-means will not always arrive at the same groupings for a specific level of species richness. Thus, you may want to run WOLPREDICTOR several times to assess optimum performance. NB optimum performance should occur whenever species clustering best reflects your empirical data’s Wolbachia distribution, so it is worth running the program a few times.

3. Predict Wolbachia strain associations (function: addPredict)
The program then assigns Wolbachia infection statuses. Under our working hypothesis, we expect that closely related, or contemporarily diverging, species in close ecological contact will host different Wolbachia strains that mediate reproductive isolation. Thus, WOLPREDICTOR assigns arbitrarily named strains to putative species groupings that co-occur in the same community. For example, if community #1 (labelled 'co1' in our data) has been allocated three putative species groupings, the program will assign the strain names ‘co1_w1’, ‘co1_w2’ and ‘co1_w3’ to all individuals in the three putative clusters. However, if community #2 has only been allocated a single putative species grouping it will be assigned the strain name ‘noWol’ (i.e. no Wolbachia) as there is no additional putative species in the community from which Wolbachia should mediate reproductive isolation. 

4. Remove Wolbachia strains according to co-phenetic distance between species clusters (function: wolPurger)
As Wolbachia has been shown to typically drop-out from host lineages after approx. 7 million years, our working hypothesis includes the possibility that species clusters separated by sufficient evolutionary time will have undergone Wolbachia purging. WOLPREDICTOR uses the co-phenetic distance matrix to incrementally remove Wolbachia according to cut-off thresholds between putative species clusters within a community. Species clusters that are separated by a distance greater that the incremental threshold will have their assigned Wolbachia strain converted to ‘noWol’. This is done conservatively, if there are three species within a community and the evolutionary distances between two of them are greater than the threshold, Wolbachia will not be purged if their evolutionary distances from the third species are below the threshold. Calculations are made for all threshold purging cut-offs (using the parameter ‘purge’, which sets the maximum purging distance iterated to from zero, at intervals set by ‘pge_incr’). At this point, all Wolbachia strain predictions at every species delimitation iteration are written to a CSV file with a root name ‘wolPreds_threshClades_incr’; this includes calculations for species delimitation iterations without any purging having been performed.

5. Match arbitrarily assigned Wolbachia strain names to empirically derived strain names (function: matchStrains)
In order to calculate the accuracy of our model at each iteration, we need to match arbitrarily assigned Wolbachia strain names to the empirically derived strain names, as much as is possible. Again, this is done conservatively, if a community features two putative species delimitation clusters whose individuals all match up to a single Wolbachia strain, only one species will have their arbitrarily assigned Wolbachia strain names converted to match the empirical predictions. For example, if community #1 has two species clusters (e.g. ‘co1_w1’ and ‘co1_w2’) whose individuals have been shown to have the empirical strain ‘wspC6’ (see our data), then only the species cluster with the most individuals will be reassigned to the strain name ‘wspC6’. The second putative species will retain the arbitrarily assigned strain name (i.e. ‘co1_w1’ or ‘co1_w2’) and will therefore not contribute to the subsequent assessment of prediction accuracy. NB if a species cluster in a separate community should also have the strain ‘wspC6’, WOLPREDICTOR will allow it to be named as such (as is seen to occur in natural populations).

6. makePDF – output figures (.pdf & .png) of prediction accuracy
We use the Python package matplotlib to output a graphical representation of predictive accuracy performed at each iteration. The flag: `-g`: indicates the gap between ticks on figure x-axis [default = 10; use multiples of 10 e.g. 100, only]. The flag: `-q`: can be used to switch off makePDF with `-q 0` [default: `-q 1`].
	
7. The program ‘outputInvestigator.py’ is a short script to investigate performance from the output files of the main program. 

	run as: `python outputInvestigator.py filename arg2 arg3`

The first argument is the full name of the output file with the root name "correctedWolPreds"; the 2nd and 3rd arguments are the minimum number of negative and positive Wolbachia strains that you require highlighting. The output will indicate which runs of the ‘wolPredictor’ program scored higher than the thresholds set here. There maybe a trade-off between negative and positive strain assignment by wolPredictor so it advisable to play with these two thresholds – starting low then increasing, to highlight the best runs. From the output, individual columns of the output CSVs (whose column names will be outputted) can be examined to see how predictions were made.

