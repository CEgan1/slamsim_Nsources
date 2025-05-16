Required softwares and packages:
-----------------------
<p><strong><a href="https://messerlab.org/slim/">SLiM</a> </strong></p>
<p><strong>Python:</strong>
<li>numpy</li></p>
<p><strong>R:</strong>
<li>e1071</li>
<li>optparse</li></p>
<p><strong>Other:</strong>
<li><a href="https://github.com/szpiech/asd/blob/master/README">asd software</a></li>
<li><a href="https://www.gnu.org/software/parallel/">GNU parallel</a></li>
<li><a href="https://vcftools.sourceforge.net/">vcftools</a></li></p>

To launch script, run :
-----------------------
<code>./SLAMSim.py --p SLAMSim.parfile</code> in the slamsim_Nsources folder.

Arguments are passed through the parfile:
-----------------------
<p>Format of the parfile is illustrated in the example parfile SLAMSim.parfile. Put ': ' between parameter name and parameter value. Details:</p>
<li>SLiMSimulation is the path to the SLiM script</li>
<li>SimulationParameters is the path to the folder where SLiM inputs are</li>
<li>ParallelJobs is the number of jobs to run in parallel</li>
<li>SourcePopulation is the number of source populations in the admixture process</li>
<li>SampleSize is the number of individual to sample in each population to compute summary statistics</li>
<li>bp is the length of the autosomal genome (in bp)</li>
<li>vcf is the path to the vcf for autosomes</li>
<li>fAM is the range of fAM to try in simulations. Format is min-max</li>
<li>sigma is the range of sigma to try in simulations. Format is min-max</li>
<li>RecombinationMap is a recombination rate map (two columns, one containing breakpoints and the other containing recombination rate between this breakpoint and the previous one) for all autosomes and X chromosome if applicable</li>
<li>ASDpath is the path to <a href="https://github.com/szpiech/asd/blob/master/README">asd software</a> (on user's computer)</li>
<li>activateX is YES to simulate an X chromosome and NO (or anything else) otherwise</li>
<li>activateY is YES to simulate an Y chromosome and NO (or anything else) otherwise</li>
<li>bpX is the length of the X chromosome (in bp). Optional argument, delete the line if there is no X.</li>
<li>bpX is the length of the Y chromosome (in bp). Optional argument, delete the line if there is no Y.</li>

SLiM Inputs
---------------------
<p>The simulation inputs include population sizes and genetic contributions of the source populations to the admixed populations at each generation.</p>
<li>A column named g and contains the generation number.</li>
<li>Columns named N1, N2 etc. containing the population sizes of each source populations at each generations.</li>
<li>A column named Nadm, containing the admixed population's sizes at each generation.</li>
<li>Columns named c1, c2, etc. containing the source populations' contributions to the admixed population at each generation. </li>

Output
---------------------
The output will be summary statistic files .sumstat. For each simulations, they will be found in the same folders as the parameters files. For each simulations, the social_ID file recalls the tag used by slim simulations for all individuals at each generation in the admixed population. There will also be a fAM_file, recalling the value of fAM for each simulation and a sigma_file recalling the value of sigma for each simulation. 
(As of now, summary statistics are computed on autosomes only...)

Example with 3 populations:
---------------------
<p><em>You may need to modify authorizations for the scripts by running: <code>chmod +x *</code> in the slamsim folder.</em></p>
<p><strong>Extend parameters files to make them compatible (up to 4 theoretical sources)</strong></p>
<p><code>./extend_parfiles.py --path test/</code></p>
<p><strong>Example with toy source vcf and default recombination rate (10⁻⁸), autosomes only:</strong></p>
<p><em>Don't forget to add path to asd software in the parfiles. It should end with /asd/src/asd.</em></p>
<p>Then run the following in the slamsim folder:</p>
<p><code>./SLAMSim.py --p SLAMSim_autosomes.parfile</code></p>
<p><strong>Example with toy source vcf and default recombination rate (10⁻⁸), autosomes and X chromosome:</strong></p>
<p><code>./SLAMSim.py --p SLAMSim_noY.parfile</code></p>
<p><strong>Example with toy source vcf and default recombination rate (10⁻⁸), autosomes, X chromosome and Y chromosome:</strong></p>
<p><code>./SLAMSim.py --p SLAMSim.parfile</code></p>
