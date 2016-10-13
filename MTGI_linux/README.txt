
MTGIpick - A software developed intending to identify Genomic Islands
----------------------------------------------------------------------
1. What is MTGIpick?

MTGIpick implements a multiscale statistical algorithm to predict GIs from a single genome. It uses a small-scale test with large-scale features to score small region deviating from the host and a large-scale statistical test with small-scale features to identify multi-window segments for identification of GIs. MTGIpick can identify GIs from a single genome, without annotated information of genomes or prior knowledge from other datasets.
2. Requirements

MTGIpick has been compiled and tested under Sun Java interpreter and Matlab. MTGIpick can be used in Windows- and Linux- platforms. Java Virtual Machine and MATLAB Compiler Runtime (MCR) are required for MTGIpick setup on your platform. However, we strongly advise the use of openjdk (JDK) instead of the Oracle version of java virtual machine when working in linux-based machines as the Oracle version may result in some exceptions during the analyses.
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
_ _ Software _ _ _ _ _ _ _ _    On window(x64)_ _ _ _ _ _ _ _ On Linux(x86_64)_ _ _ _ _ _ 
  Java Virtual Machine		  JDK 1.8	                 JDK 1.8                  
  MATLAB Compiler Runtime	  MCR 8.4	                 MCR 8.1             	  
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _


3. Download

There are two ways to download the MTGIpick:
   1) You can download the MTGIpick package with JDK and MCR from our web ( http://bioinfo.zstu.edu.cn/MTGI )
   A) Windows (MTGIpick.zip)
       - MCRInstaller.exe                           # MCR 8.4 for Windows
       - MTGI_setup.exe                             # Main program
       - Example.fasta                              # Two sequences in FASTA format
       - README.txt                                 # Documentation
       - help.htm (help.files)                      # Introduction
   B) Linux (MTGIpick.zip)
       - jdk-8u102-linux-x64.tar.gz                 # JDK 1.8 for Linux
       - MCRInstaller.zip                           # MCR 8.1 for Linux
       - MTGI_linux.jar                             # Main program
       - install_jdk-mcr_linux.sh                   # Install JDK 1.8 and MCR 8.1
       - run_MTGI_linux.sh                          # Run MTGIpick software
       - Example.fasta                              # Two sequences in FASTA format
       - README.txt                                 # Documentation
       - help.htm (help.files)                      # Introduction
       - jsonFile (d3)                              # Visualization files
	   
   2) If you download the MTGIpick package without JDK or MCR from our webs ( http://bioinfo.zstu.edu.cn/MTGI, https://github.com/bioinfo0706/MTGIpick ), download the JDK and MCR for your platform from the following Web:
   JDK 1.8:http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html#jdk-8u102-oth-JPR
   MCR: http://www.mathworks.com/products/compiler/mcr/?refresh=true
   Make sure that the JDK and MCR are saved into the folder of the MTGIpick, and please rename the MCR as MCRInstaller.

4. Install and Run
  
   1) Windows (Tested on win7 x64)
   Before installing MTGIpick, make sure that the MCR 8.4 for windows is saved into the folder of MTGIpick software. Install MCR 8.4 first, and run MTGIpick setup directly.
   2) Linux (Tested on Centos7)
   Before installing MTGIpick, make sure that the JDK 1.8 (jdk-8u102-linux-x64.tar.gz) and MCR 8.1 for Linux are saved into the folder of MTGIpick software. Please follow the following steps for installing and running MTGIpick:
   Step 1
   To install JDK 1.8 (jdk-8u102-linux-x64.tar.gz) and MCR 8.1, it requires a simple command line as follow:
      > bash install_jdk-mcr_linux.sh
   Step 2
   To run the MTGIpick software, just type a simple command line as follow (Once the first step has run, execute the second step to run MTGIpick):
      >bash run_MTGI_linux.sh

5. Inputs

MTGIpick accepts DNA sequences, and the input file has to be in fasta, fa, or fna format. For example, the file name is example.fasta or example.fa, its content looks like this:  
>cya_GI1_93960_99000_GI2_208020_223020_GI3_408480_410520
CCCCATTCCCCCCATTCCCTCCTTTTCCACCATACCCTCTTTTCCCCTCGTTGCCCCCAA
ATTTTTACGCATTTCCCCATTAATGCGATGATCCCAGCGCGAAAGCATCTGTGATTAAGA
CGTCTATCAATTTATACTCGTTAGGGTTTTTTCTTCGGTGGTACCATCTGGGCGCCTACG
>::

Example DNA input (.fna format), its content looks like this:
>cya_GI1_93960_99000_GI2_208020_223020_GI3_408480_410520
CCCCATTCCCCCCATTCCCTCCTTTTCCACCATACCCTCTTTTCCCCTCGTTGCCCCCAA
ATTTTTACGCATTTCCCCATTAATGCGATGATCCCAGCGCGAAAGCATCTGTGATTAAGA
......
>ctt_ GI2_208020_223020_GI1_93960_99000_GI3_408480_410520
CCCCATTCCCCCCATTCCCTCCTTTTCCACCATACCCTCTTTTCCCCTCGTTGCCCCCAA
ATTTTTACGCATTTCCCCATTAATGCGATGATCCCAGCGCGAAAGCATCTGTGATTAAGA
......

Please note GenBank (from NCBI) will not work with MTGIpick. If using sequences from NCBI be sure to save them as FASTA format first. Sequence format conversion tools are available at http://www.ebi.ac.uk/Tools/sfc/. 
The sequences can be uploaded to MTGIpick in a file. It is very important that each of the sequences has a unique name. If they do not, the software will fail. There must be no empty lines, white spaces or control characters between sequences or at the top of the file. This will also cause the software to fail.

6.Output

The outputs of the MTGIpick consist of genomic signatures, conserved scores of predicted GIs and predicted GIs. They are stored in the same directory where the input file is stored. Output of the genomic signatures is a Zip file whose name is created by the input file name and signature. If the input file contains at least two sequences, the Zip file contains more signature files for all the sequences. 
  1) *_Total_SPGIs_*.gff
  2) *_Each_SPGIs_*.txt
  3) *_signa_*.txt 
  4) *_html_*.html

If you have any problem, please contact:
     
     daiailiu04@yahoo.com

