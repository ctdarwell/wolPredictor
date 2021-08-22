# wolPredictor
**Predict distribution of endosymbiotic bacteria among insects**

*Methodological outline of wolPredictor* 

Outline of Problem:
It has been regularly suggested that highly prevalent Wolbachia induced reproductive isolation among arthropods appears randomly distributed among closely related host species. If so, this implies that much arthropod biodiversity is a result of stochastically determined diversification events rather than process driven outcomes. For most arthropods we have limited knowledge about ecological contact which provides direct opportunity for horizontal exchange of microbes or genetic material between species. Here, I use Python programming to model our proposed mechanism that incorporates ecological contact in a phylogenetic context. 

Providing all correct Python and R libraries are installed the program should run directly from the unzipped download bundle. The program has not yet been optimised for vectorisation or parallelisation (NB *wolPredictor_xmeans.py* runs very quickly).

PYTHON AND R LIBRARIES:

(i) The following library dependencies are required for the Python scripts to run:

*numpy* - numpy.org

*pandas* - pandas.pydata.org

*pyclustering* - pyclustering.github.io/docs/0.8.2/html/index.html (NB for wolPredictor_xmeans.py only)

*matplotlib* - matplotlib.org

(ii) R requires the *ape* library installing - `install.packages(ape)`


DATA REQUIREMENTS:

(i) A CSV file with data columns featuring individual arthropod sample names, empirically derived Wolbachia strains (or absence of strain) for each sample (IMPORTANT: absence of Wolbachia should be coded as 'noWol'), and host community name (here labelled "*sp.complex*" - IMPORTANT: the first three characters of each community's name must be a unique combination of characters although the actual name can be longer than three characters). For using ecoCladeGenerator.py you require a column of ecological categories (our data is elevation - column name "elevation" in our data)

(ii) A phylogenetic tree in nexus format. NB The program calls R (`cophen4pyOut.R`) using the Python library ‘subprocess’. If you cannot configure R (should be straightforward by adding R to your PATH in Windows: e.g. *C:\Program Files\R\R-3.6.2\bin* on my machine) to interact with Python you must use `wolPredictor_xmeansDelim_reduci.py` (see below) with the ‘exaData_faoReview.tre_cophen.csv’ file for a test run (for your own analyses without configuring R you will simply have to create a distance matrix of co-phenetic phylogenetic distances as formatted in the CSV file – see cophen4pyOut.R for the code in R create the file, NB replace the last line of R code with:  `write.csv(as.matrix(phydist), file='yourFile.csv', quote = F, row.names = F)` ).


**Initial consideration**

*wolPredictor* basically comes in two flavours according to method employed to delineate species richness clusters:

OPTION #1: *wolPredictor_MANUAL.py* requires an initial file to be generated that features different combinations of species designations among host insect samples in the dataset. For this purpose, the program *ecoCladeGenerator.py* divides samples into species clusters according to associated ecological categories within individual communities. For example, for our data it is clear that our wasp lineages are strongly correlated with population elevation. So, for a community featuring samples collected at 100m, 200m and 300m elevation, *ecoCladeGenerator.py* creates all combinations of species clusters: [[100, 200, 300]], [[100, 200], [300]], [[100], [200, 300]] and [[100], [200], [300]]. Whether it will also include species clusters featuring disjunct elevations (i.e. [[100, 300], [200]]) is decided by user input.

Workflow1: 1. *ecoCladeGenerator.py*; 2. *wolPredictor_MANUAL.py*; 3.*taxdegMatcher.py*; 4. *wolTabber.py*

OPTION #2: *wolPredictor_xmeans.py* divides samples into species clusters according to an input phylogenetic tree. Iteratively, is divides species into a (user inputted) range of species richness values using X-means evaluation of pairwise branch length distances.

Workflow2: 1. *wolPredictor_xmeans.py*; 2. *taxdegMatcher.py*; 3. *wolTabber.py* (replace *taxdegMatcher.py* and *wolTabber.py* with *outputInvestigator.py* - if there is no a priori hypothesis of species richness featured in the main data file)

**NB** It is recommended to use *wolPredictor_xmeans.py* when numbers of communities gets high (e.g. above 8 or so) as the number of permutations outputted by *ecoCladeGenerator.py* grows exponentially; or, if there is no apriori species richness hypothesis. However, using *ecoCladeGenerator.py* and *wolPredictor_MANUAL.py* may allow a more systematic and detailed investigation of partterns.

**OPTION #1 *wolPredictor*:** How to run *ecoCladeGenerator.py* and *wolPredictor_MANUAL.py*: 

*ecoCladeGenerator.py* can be run from a Python3 console (e.g. Spyder/Anaconda). Typing:

`python ecoCladeGenerator.py`

The following flags can be added to alter the default setting variables:

`-s` swch [default: 1]: This determines whether disjunct populations can be considered when deciding species cluster designations. The option are [0 - allow disjunct combinations] and [1 - don't allow]. This obviously only makes sense if the ecological variable being considered has some kind of linear relationship (for our data, it makes sense not to allow species clusters of low and high altitude populations).

`-e` ecoVar [default: "elevation"]: The column name of the ecological variable being examined. Variables can be added as any type e.g., categorical or numerical as long as they form discrete categories (e.g. our elevation variables consist of category values such as 200, 700, 1200 etc that repesent mean population altitude rather than the exact elevation of each sample)

`-c` community [default: "community"]: The column name of host communities. For our example we are considering host fig species.

`-d` filename [default: "exaData_faoReview.csv"]: CSV file of data

`-p` prfx [default: "eco"]: Prefix of outputted file.

Example: `python ecoCladeGenerator.py -d myData.csv -c species` - runs the program for a data file where community has the column heading "species"


How to run *wolPredictor_MANUAL.py*: The program can be run from a Python3 console (e.g. Spyder/Anaconda). Typing:

`python wolPredictor_MANUAL.py -h`

will bring up a menu of parameter options. Most are straightforward and relate to data input option (filenames & directories etc). However, key considerations are:

`-p` prefix: Prefix of outputted file [default: "manual"].

`-m` filename: Main input filename [defaults: "exaData_faoReview.csv"]

`-q` pdf: Make figure (Off/On: 0/1). More suited to *wolPredictor_xmeans.py* but can be run for general overview of results [default: 0].

`-g` gap: Gap between ticks on figure x-axis  [default: 10].

`-s` shuffle: Shuffle species clusters (Off/On: 0/1). Perform null analyses [default: 0].

`-C` cntrl: Control strains in multiple communities (Off/On: 0/1) [default: '0']. Setting to "1" performs a more conservative analysis, where members of the same designated clade occupying different communities will only be permitted a single matched strain (i.e. in one community only)

`-N` nPges: No. of purges at each species delim iteration (max=10) [default: '4']. This actually gives a non-precise outcome as it is used to calculate an integer from a division sum. The default should give a nice spread of purging thresholds (either 25%, 50%, 75% and 100% OR 20%, 40%, 60%, 80% and 100% of phylogenetic tree maximum pairwise branch length) - setting to 10 would give a more fine-scaled evaluation but is unlikely to be informative relative to the default and will just increase runtime.

So, to run the program you might type:

`python wolPredictor_MANUAL.py -m newDataFile.csv -d data_directory -o out_directory -p EcoVariables -N 10`


**NB** as it stands *ecoCladeGenerator.py* will output a file: *ecoDelim_communityXelevation_constr.csv* featuring 8192 permuations of species delimitations. This will take around 80 minutes to run in *wolPredictor_MANUAL.py*. We provide a smaller test file - *randomDegs_n200.csv* - to run: `python wolPredictor_MANUAL.py -r randomDegs_n200.csv`

You can also run the program from a Python IDE and alter the above parameters in the code.

After running the program you should use *taxdegMatcher.py* and *wolTabber.py* for analysis and evaluation of results.


**OPTION #2 *wolPredictor*:** How to run *wolPredictor_xmeans.py*: The program can be run from a Python3 console (e.g. Spyder/Anaconda). Typing:

`python wolPredictor_xmeans.py -h`

will bring up a menu of parameter options. Most are straightforward and relate to data input option (filenames & directories etc). However, key considerations are:

`-m` min_nSpp: The minimum number of putative species clusters that X-means will divide the phylogenetic data into. NB must be ≥2 [default = 2].

`-M` max_nSpp: The maximum number of putative species clusters that X-means will divide the phylogenetic data into [default = 50].

`-i` nSpp_incr: The incremental increase between min_nSpp and max_nSpp [default = 1].

`-m` & `-M` should be reasonable guestimates of minimum and maximum species richness from your tree. `-i` should be a reasonable increment between the two values (so that the runtime is not overlong: progress print-outs are outputted so you can judge if a run will take too long – change `-i` if necessary).

See also above options with *wolPredictor_MANUAL.py*

**NB Description of function *calcX* in *wolPredictor_xmeans.py* via X-means species delimitation:** 
In order to create putative species groupings, the program uses the pairwise co-phenetic distance matrix and X-means clustering to delimit species groupings. At each loop iteration, `wolPredictor` will calculate clusters of species for each value of species richness (i.e. between `-m` & `-M`). It should be emphasised that X-means will not always arrive at the same groupings for a specific level of species richness. Thus, you may want to run `wolPredictor` several times to assess optimum performance. NB optimum performance should occur whenever species clustering best reflects your empirical data’s Wolbachia distribution, so it is worth running the program a few times.

So, to run the program you might type:

`python wolPredictor_xmeans.py -d data_directory -o out_directory -m 5 -M 100 -i 5`

You can also run the program from a Python IDE and alter the above parameters in the code.

Also use *taxdegMatcher.py* and *wolTabber.py* for analysis and evaluation of results.


Explanation of *wolPredictor* functions common to *wolPredictor_xmeans.py* and *wolPredictor_MANUAL.py*

1. R_cophen: call R (‘ape’ library) to create cophenetic distance matrix

2. Uses X-means clustering to divide the branch length pair-wise distance matrix into *n* clusters (*wolPredictor_xmeans.py* only)

3. Predict Wolbachia strain associations (function: addPredict). This relates to the **contact contingency hypothesis** in the paper.

	The program then assigns Wolbachia infection statuses. Under our working hypothesis, we expect that closely related, or contemporarily diverging, species in close ecological contact will host different Wolbachia strains that mediate reproductive isolation. Thus, `wolPredictor` assigns arbitrarily named strains to putative species groupings that co-occur in the same community. For example, if community #1 (labelled 'co1' in our data) has been allocated three putative species groupings, the program will assign the strain names ‘co1_w1’, ‘co1_w2’ and ‘co1_w3’ to all individuals in the three putative clusters. However, if community #2 (i.e. 'co2') has only been allocated a single putative species grouping its wasps will be assigned the strain name ‘noWol’ (i.e. no Wolbachia) as there is no additional putative species in the community from which Wolbachia should mediate reproductive isolation. If community #3 (i.e. 'co3') is allocated two putative species, the strain names 'co3_w1' and 'co3_w2' will be assigned. 

4. Remove Wolbachia strains according to co-phenetic distance between species clusters (function: wolPurger). This relates to the **adaptive decay hypothesis** in the paper.

	As Wolbachia has been shown to typically drop-out from host lineages after approx. 7 million years, our working hypothesis includes the possibility that species clusters separated by sufficient evolutionary time will have undergone Wolbachia purging. `wolPredictor` uses the co-phenetic distance matrix to incrementally remove Wolbachia according to cut-off thresholds between putative species clusters within a community. Species clusters that are separated by a distance greater that the incremental threshold will have their assigned Wolbachia strain converted to ‘noWol’. This is done conservatively, if there are three species within a community and the evolutionary distances between two of them are greater than the threshold, Wolbachia will not be purged if their evolutionary distances from the third species are below the threshold. Calculations are made for all threshold purging cut-offs (using the parameter ‘purge’, which sets the maximum purging distance iterated to from zero, at intervals set by ‘pge_incr’). At this point, all Wolbachia strain predictions at every species delimitation iteration are written to a CSV file with a root name ‘wolPreds’; this includes calculations for species delimitation iterations without any purging having been performed.

5. Match arbitrarily assigned Wolbachia strain names to empirically derived strain names (function: matchStrains)

	In order to calculate the accuracy of our model at each iteration, we need to match arbitrarily assigned Wolbachia strain names to the empirically derived strain names, as much as is possible. Again, this is done conservatively, if a community features two putative species delimitation clusters whose individuals all match up to a single Wolbachia strain, only one species will have their arbitrarily assigned Wolbachia strain names converted to match the empirical predictions. For example, if community #1 has two species clusters (e.g. ‘co1_w1’ and ‘co1_w2’) whose individuals have been shown to have the empirical strain ‘wspC6’ (see our data), then only the species cluster with the most matching individuals will be reassigned to the strain name ‘wspC6’. The second putative species will retain the arbitrarily assigned strain name (i.e. ‘co1_w1’ or ‘co1_w2’) and will therefore not contribute to the subsequent assessment of prediction accuracy. This is important, if this pattern of single Wolbachia strains across multiple species within a community was found in empirical data, it would suggest our model did not accurately represent natural patterns. Thus, our accuracy assessment is conservative and means that only empirical datasets that reflect our theoretical propositions will return high accuracy assessments. Thus alternatively, if the two species clusters’ (‘co1_w1’ and ‘co1_w2’) individuals were shown to have two distinct empirical strains (e.g. ‘wspC5’ and ‘wspC6’), both would contribute to higher accuracy assessments (both receiving reassigned strain names) and corroborate our model. NB if a species cluster in a separate community should also have the same strain (e.g. ‘wspC6’), WOLPREDICTOR will allow it to be named as such (as can be seen to occur in natural populations and is not in violation of our predictions).

6. makePDF3 – output figures (.pdf & .png) of prediction accuracy

	We use the Python package matplotlib to output a graphical representation of predictive accuracy performed at each iteration. The flag: `-g`: indicates the gap between ticks on figure x-axis [default = 10]. The flag: `-q`: can be used to switch off makePDF with `-q 0` [default: `-q 1`]. The two flavours of *wolPredictor* determine some of the operating parameters: *tix* (0/1) determines whether the figures gets x-axis ticks (innappropriate for *wolPredictor_MANUAL.py*); *spl* ('div'/'nSpp') encodes a string used to process the data.


OTHER PROGRAMS

**(1) *taxdegMatcher.py*** evaluates whether the species delimitation clusters that are iterated thru by *wolPredictor* and generated by either *ecoCladeGenerator.py* or *wolPredictor* itself. 

It can be run as: `python taxdegMatcher.py filename1 filename2 arg3 arg4 arg5`

The user-inputted variables are [1] the name of the data file containing the sample taxon names and a column indicating the samples true species designation; [2] a file containing the species clusters generated by *ecoCladeGenerator.py* or when running *wolPredictor*; [3] the column name of sample names in the data file; [4] the column name of sample names in the species clustering file; [5] the columns name of the true species designations of species int he data file. It outputs a single file of matching values.



**(2) *wolTabber.py*** will generate a results table at each considered species clustering iteration. 

It can be run as: `python wolTabber.py filename1 filename2 filename3 arg4`

These are: [1] - the "correctedWolPreds" file outputted by *wolPredictor*; [2] - the "taxonDesignations" file outputted by *wolPredictor*; [3] - the file outputted by *taxdegMatcher.py*; [4] - a prefix for the outputted file name

Outputted columns are [1] 'match' - species delimitation accuracy at this species clustering iteration; [2] 'pgeIncr' - the best 'noWol' (i.e. uninfected with Wolbachia) prediction improvement after the *wolPurger* has been applied at different thresholds at each species clustering iteration; [3] 'pgeIncr_wo_loss' - the best 'noWol' prediction improvement after the *wolPurger* has been applied at different thresholds at each species clustering iteration without compromising positive strain prediction abilities; [4] 'maxScore' - the best combined  prediction score of positive strain prediction + 'noWol' prediction including the effcet of *wolPurger*; [5] 'noPge_noWol' - the 'noWol' prediction accuracy without the action of the *wolPurger* function at each species clustering iteration; [6] 'noPge_posWol' - the  Wolbachia strain prediction without the action of the *wolPurger* function at each species clustering iteration; [7] - 'match_nSpp' species richness at each species clustering iteration.

The format of the outputted file means it is simple to perform regression analyses on outputted data: e.g., *noPge_posWol ~ match*


**(3) *outputInvestigator.py*** is a short script to investigate performance from the output files of the main program. The program *wolTabber.py* is much better for subsequent analyses/evaluation.

run as: `python outputInvestigator.py filename1 filename2 arg3 arg4 arg5`

   The first and second arguments are the full name of the output files with the root name "correctedWolPreds" and "taxonDesignations"; the 3rd, 4th and 5th arguments are the minimum number of negative, positive and combined Wolbachia strains that you require highlighting. The output will indicate which runs of the ‘wolPredictor’ program scored higher than the thresholds set here. There maybe a trade-off between negative and positive strain assignment by wolPredictor so it advisable to play with these two thresholds – starting low then increasing, to highlight the best runs. From the output, individual columns of the output CSVs (whose column names will be outputted) can be examined to see how predictions were made.



**(4) *wolPredictor_xmeansDelim_reduci.py*** i.e. If R cannot be configured to interact with Python.

Variables must be set manually in the script text under the heading: "SET SOME PARAMETERS!!!!!" at Line 15. Here you must set the variables `purge` and `pge_incr`. To set purge you need the maximium cophenetic distance for your tree from R - `max(cophenetic(tree))`. If max distance is 0.17, enter "20" for `purge`. For `pge_incr` choose a sensible number of increments e.g. if `purge` is "20" set `pge_incr` as "2".
	
NB this program cannot randomly shuffle the Wolbachia strain data (`-s`) or perform the more conservative analysis (`-C`)


**(5) *ciFitnessModel.py*** is a standalone program to evaluate relative fitness outcomes between diverging lineages where cytoplasmic incompatibility (CI) is or is not operating under the assumption that CI entails a reduction in fecundity. This relates to the **fecundity trade-off hypothesis** in the paper.

run as: `python ciFitnessModel.py`

It models how offspring fitness might benefit from preferential ovipositioning facilitated by reduced foundress egg load (or via selective ovipositioning) despite losses in fecundity resulting from the post-zygotic mechanism of CI. The model allocates each foundress 1000 eggs whose fitnesses are determined by three variables: (i) the proportion of conspecifics in the population that the foundress can breed with; (ii) fitness of offspring of conspecific matings (between 0-1); and (iii) fitness of offspring of heterospecific matings (between 0-1). Further, a fig syconium is divided into five layers - at the centre, all eggs survive, while zero survive at the syconium edge (with 0.75, 0.5, 0.25 survival in the intermediate layers). We compares a 'CI' and a ‘no CI' model across variation in our three variables.
	
For example, for the no CI model at 75% conspecifics mating opportunities: 75% of eggs get conspecific relative mating fitness (e.g. ω = 0.8) and 25% of eggs get heterospecific relative mating fitness (e.g. ω = 0.2). Each egg gets randomly assigned into a layer of the syconium and each egg’s relative mating fitness is multiplied by its oviposition site fitness (i.e. randomly assigned: 1.0, 0.75, 0.5, 0.25, 0). Finally, all 1000 eggs’ scores are summed to obtain inclusive fitness for the foundress. For the CI model (at 75% conspecifics mating opportunities) we lose 25% of the egg load but all eggs get the conspecific relative mating fitness (i.e. 0.8) - then they are randomly placed into the best remaining 750 oviposition site positions (mimicking that fig wasps preferentially lay in optimum sites). Scores are again multiplied and summed across 10 replicates. The key modeling assumption here is that unviable eggs are not oviposited and thus do not waste premium oviposition sites (this is also equitable to a scenario where egg oviposition order is prioritised). Results are obtained by subtracting realised CI versus no CI fitnesses across all combinations of conspecific-heterospecific relative mating fitness (i.e. 0-1 for both) surfaces. Individual heat maps for outcomes at different percentages of conspecific mating opportunities are outputted.

Some variables can be altered at the beginnng of the script but this is mostly unnecessary as the program is modeling all fitnesses between 0 and 1 (no other values make sense). Increasing foundress egg load (`fec`) to more than 1000 would add greater but unnecessary precision. The `fldr` variable is an output directory and can be modified accordingly. The variable `oviLayers` may be worth exploring if empirical data are available for larval survival percentages from layers within a syconia. The program takes around five minutes. I have vectorised where I see possible. There are nested loops in the `Main` section - but to vectorise further would require some nested vectorisation of the `numpy random.choice` method (if it's possible, please let me know).
