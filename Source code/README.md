
MTGI function implements multiscale statistical algorithm to predict genomic islands (GIs) from a single genome.

Choices are: 

1) Input file
2) Process-input: -> Process Input file
     (1)Predict-separately (default)   -  For each DNA sequence of the input file, MTGIpick will process it separately and predict its GIs.
     (2)Assemble-predict               -  The sequences from the input file are assembled into a sequence according to the order of these sequences in the input file, and MTGIpick will process it as a sequence and predict its GIs.                                     
3)Word-size: -> the length of k-mer (2-5), and default is 4.
4)Transformed-window: -> the total number of the windows used in genomic transformation (1-4), and default is 4.
5)teration-time: -> the periods of time that are repeated to select windows whose scores are large enough to be considered statistically significant (1-10), and default is 5.
6)Core-feature: -> the size of selected features by the proposed kurtosis (1-4^k), and default is 256.
7)Eye-window: -> the size of eye windows used in the proposed divergence measure based on two sample t-test (1-6), and default is 5.
8)Standard-error-IST: -> the standard deviation of the mean of the window scores to select windows associated with putative GIs (0-1), and default is 0.05.
9)Total-scale: -> the total number of the scales in the multiscale segmentation algorithm (1-50), and default is 30.
10)Standard-error-MSA: -> the standard deviation of the mean of enrichment scores to select segments associated with putative GIs (0-1), and default is 0.3.
11)Minimum-GI:-> the minimum size of predicted genomic islands, and default is 10000. A smaller value is recommended if you want to predict GIs with small size.
12)Boundary-detection: -> 1 selects boundary detection algorithm, and 0 does not selects boundary detection algorithm
13)Upstream-RGI: -> the length of sequences around ¡®raw¡¯ GIs to refine the boundaries of predicted GIs, and default is 5000.
14)Downstream-RGI: -> the length of sequences around ¡®raw¡¯ GIs to refine the boundaries of predicted GIs, and default is 5000.

For example,

MTGI('Example.fasta','Process-input','Assemble-predict','Word-size',4,'Transformed-window',1,'Iteration-time',5,'Core-feature',256,'Eye-window',5,'Standard-error-IST',0.05,'Total-scale',30,'Standard-error-MSA',0.3,'Minimum-GI',10000,'Boundary-detection',1,'Upstream-RGI',5000,'Downstream-RGI',5000);


MTGI('Example.fasta','Process-input','Predict-separately','Word-size',4,'Transformed-window',1,'Iteration-time',5,'Core-feature',256,'Eye-window',5,'Standard-error-IST',0.05,'Total-scale',30,'Standard-error-MSA',0.3,'Minimum-GI',10000,'Boundary-detection',1,'Upstream-RGI',5000,'Downstream-RGI',5000);

