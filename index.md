
# MIPanalyzer

MIPanalyzer is an R package for working with molecular inversion probe (MIP) data to perform various analysis and visualisation functions. The motivation for developing this package was that `vcf` data can be difficult to work with due to the way the important information is locked away inside a very specific character string format. What we want most of the time is simply the number of reads of the major allele (which we call `counts`), and the total number of reads (which we call `coverage`). From these two numbers alone we can calculate within sample allele frequencies, perform filtering, calculate relatedness, etc.

The second motivation was to represent this information in an efficient way. Most of the time we will have this information for every locus and every sample (apart from missing data), meaning long format is not efficient because the same locus/sample information is repeated many times. We opt for a simple matrix for `counts` and `coverage`, with the sample and locus data being stored in separate objects that always retain the correct dimension.

Finally, this package is intended to evolve and develop as new MIP analyses are needed. For example, rather than developing a nice figure for a paper in standalone script, this package provides a location to write this as a general function so it is available to others in the future. If you work on MIPs and there is something you would like to see in this package, then please build it an contribute!

Current main functionalities include:

- Filtering samples or loci based on `counts`, `coverage`, missingness, or other criteria.
- Population structure analysis through principal components analysis (PCA) and visualisation of these results.
- Calculation of pairwise genetic distances, including a maximum-likelihood F estimate.
- Basic simulation of MIP data
