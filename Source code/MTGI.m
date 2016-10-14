% Qi Dai
% College of Life Sciences, Zhejiang Sci-Tech University, Hangzhou 310018, China
% Department of Biological Sciences, Center for Systems Biology, University of Texas at Dallas,
% Richardson, TX 75080, USA
%
% Oct 2016
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warranty Disclaimer and Copyright Notice
% 
% Copyright (C) 2014-2015 Zhejiang Sci-Tech University, Hangzhou 310018, China
% 
% The Zhejiang Sci-Tech University and the authors make no representation about the suitability or accuracy of this software for any purpose, and makes no warranties, either express or implied, including merchantability and fitness for a particular purpose or that the use of this software will not infringe any third party patents, copyrights, trademarks, or other rights. The software is provided "as is". The Institute for Systems Biology and the authors disclaim any liability stemming from the use of this software. This software is provided to enhance knowledge and encourage progress in the scientific community. 
% 
% This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with this library; if not, write to the Free Software Foundation,
% College of Life Sciences, Zhejiang Sci-Tech University, Hangzhou 310018, China
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function label_temp = MTGI(infile,varargin)

% MTGI() implements multiscale statistical algorithm to predict genomic islands (GIs) from a single genome.
% Choices are: 
%
% Process-input: -> Process Input file
%              Predict-separately     -  For each DNA sequence of the input file, MTGIpick will process it
%                (default)               separately and predict its GIs.
%              Assemble-predict       -  The sequences from the input file are assembled into a sequence 
%                                        according to the order of these sequences in the input file, and 
%                                        MTGIpick will process it as a
%                                        sequence and predict its GIs.
% Word-size: -> the length of k-mer (2-5), and default is 4.
% Transformed-window: -> the total number of the windows used in genomic transformation (1-4), and default is 4.
% Iteration-time: -> the periods of time that are repeated to select windows whose scores are large enough to 
%                    be considered statistically significant (1-10), and default is 5.
% Core-feature: -> the size of selected features by the proposed kurtosis (1-4^k), and default is 256.
% Eye-window: -> the size of eye windows used in the proposed divergence measure based on two sample 
%                t-test (1-6), and default is 5.
% Standard-error-IST: -> the standard deviation of the mean of the window scores to select windows 
%                     associated with putative GIs (0-1), and default is 0.05.
% Total-scale: -> the total number of the scales in the multiscale segmentation algorithm (1-50), 
%                 and default is 30.
% Standard-error-MSA: -> the standard deviation of the mean of enrichment scores to select segments 
%                        associated with putative GIs (0-1), and default is 0.3.
% Minimum-GI: the minimum size of predicted genomic islands, and default is 10000. 
%             A smaller value is recommended if you want to predict GIs with small size.
% Boundary-detection: -> 1 selects boundary detection algorithm, and 0 does not selects boundary detection algorithm
% Upstream-RGI: -> the length of sequences around ¡®raw¡¯ GIs to refine the boundaries of predicted GIs, 
%                  and default is 5000.
% Downstream-RGI: -> the length of sequences around ¡®raw¡¯ GIs to refine the boundaries of predicted GIs, 
%                  and default is 5000.

% For example
%            MTGI('Example.fasta','Process-input','Predict-separately',...
%             'Word-size',4,'Transformed-window',1,'Iteration-time',5,'Core-feature',256,...
%            'Eye-window',5,'Standard-error-IST',0.05,'Total-scale',30,'Standard-error-MSA',...
%            0.3,'Minimum-GI',10000,'Boundary-detection',1,'Upstream-RGI',5000,'Downstream-RGI',5000);

%% set defaults
Process='Predict-separately';% Process-input
k=4; % Word-size 
ww=1; %  Transformed-window
interation_whole=4; % Iteration-time
feature_size=256; % Core-feature
eye_window=5; % Eye-window
sig=0.05; % Standard-error-IST
L=30; % Total-scale
sig2=0.3; % Standard-error-MSA
sele=10000; % Minimum-GI
isbound=1; % Boundary-detection
upstream=5000; % Upstream-RGI
downstream=5000; % Downstream-RGI

%% check input sequences
  [pathstr, file_name, ext]=fileparts(infile);
  if strcmp(ext,'.fasta')==1| strcmp(ext,'.fa')==1|strcmp(ext,'.txt')==1|strcmp(ext,'.fna')==1|strcmp(ext,'.ffn')==1    
  else
      error('Bioinfo:seqpdist:IncorrectInputType',...
          'Input file must be  in .fasta, .fa, .ffn or fna format.')
  end

%% identify input arguments
if nargin > 1
    if rem(nargin,2) == 0
        error('Bioinfo:seqpdist:IncorrectNumberOfArguments',...
              'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'Process-input','Word-size','Transformed-window','Iteration-time',...
              'Core-feature','Eye-window','Standard-error-IST','Total-scale',...
              'Standard-error-MSA','Minimum-GI','Boundary-detection','Upstream-RGI','Downstream-RGI'};
          
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        kk = find(strncmpi(pname, okargs, length(pname)));
        if isempty(kk)
            error('Bioinfo:seqpdist:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(kk)>1
            error('Bioinfo:seqpdist:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(kk)
                case 1 % select method to Process input file
                    oka = {'Assemble-predict','Predict-separately'};
                    kkk = find(strncmpi(pval, oka, length(pval)));
                    if isempty(kkk)
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Unknown parameter name: %s.',pval);
                    else
                        Pro_v=pval;
                    end
                    
                case 2 % select word size
                    if  pval<2 | pval>5
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Word size should be from 1 to 5');
                    else
                        k=pval;
                    end
                    
                case 3 % select windowed transform
                    if  pval<1 | pval>5
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Transformed window size should be from 1 to 5');
                    else
                        ww=pval;
                    end
                case 4 % select iteration time
                    if  pval<1 | pval>10
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Interation time should be from 1 to 10');
                    else
                        interation_whole=pval;
                    end
                case 5 % select core feature size
                    if  pval<1 | pval>4^k
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Core feature size should be from 1 to 4^k');
                    else
                        feature_size=pval;
                    end
                case 6 % select eye window size
                    if  pval<1 | pval>6
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Eye window size should be from 1 to 6');
                    else
                        eye_window=pval;
                    end
                case 7 % select time standard error in IST-LFS
                    if  pval<0 | pval>1
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Time standard error should be from 0 to 1');
                    else
                        sig=pval;
                    end
                case 8 % select total scale
                    if  pval<1 | pval>50
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Total scale should be from 1 to 50');
                    else
                        L=pval;
                    end
                case 9 % select time standard error in MSA
                    if  pval<0 | pval>1
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Time standard error should be from 0 to 1');
                    else
                        sig2=pval;
                    end
                case 10 % select minimum GI size
                    if  pval<1
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Minimum GI size should be positive');
                    else
                        sele=pval;
                    end
                case 11 % Boundary detection
                    if  pval==0 | pval==1
                        isbound=pval;
                    else
                        error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Boundary detection should be 0 or 1');
                    end
                case 12 % select upstream of ¡®raw¡¯ GIs
                    if  pval<1
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Upstream of ¡®raw¡¯ GIs should be positive');
                    else
                        upstream=pval;
                    end
                case 13 % select downstream of ¡®raw¡¯ GIs
                    if  pval<1
                       error('Bioinfo:seqpdist:UnknownParameterName',...
                       'Downstream of ¡®raw¡¯ GIs should be positive');
                    else
                        downstream=pval;
                    end
            end
        end
    end
end  

%% Add path 
addpath(genpath(pwd));
addpath(genpath([pwd '/MTGI/']));

%% load dataset
  format short g;
  index1=findstr(infile,'.');
  out=infile(1:index1-1);
  Name=[];Seq=[];
  [Namee, Seqq] = fastaread(infile);
  if  strcmp(Pro_v,'Predict-separately')==1
      fprintf(['Predict each sequence separately \n' ]);
      % For each DNA sequence of the input file, MTGIpick will process it separately and predict its GIs.
      if (iscell(Namee))
         xh=length(Namee);      
      else
         xh=1;
      end
      for xunhuan=1:xh
          out1=strcat(out,'_score','_sequence',num2str(xunhuan),'.txt');
          out3=strcat(out,'_predict','_sequence',num2str(xunhuan),'.txt');
	      out4=strcat(out,'_html','_sequence',num2str(xunhuan),'.html');
         if (iscell(Namee))
		     Name=Namee{xunhuan};
		     Seq=Seqq{xunhuan};
		     Seq=upper(Seq);
         else
	         Name=Namee;
		     Seq=Seqq;
	         Seq=upper(Seq);
         end
         %% Extracted genomic signature
         N=length(Seq);
         mer=cell(1,1);
         each_window_length=[];
        % Calculate C+G content as genomic signatural
         cg_content=cgcomposition(Seq,'WINDOW',100);
        % Calculate k-mer as genomic signatural
         window=1000;
         slidelen=window;
         fprintf(['compute ' num2str(k) '-mer \n']);
         [mer{1} each_window_length signa_str]=cmer(Seq,window,slidelen,k);
   
        %% Calculate the score of each window using the IST-LFS algorithm
         fprintf(['IST-LFS Algorithm' ]);
         score_nucletide=[];ks_position_whole=[];ks_ascending_whole=[];
         [score_nucletide ks_position_whole ks_ascending_whole]=ISTLFS(Seq,mer{1},window,slidelen,sig,ww,interation_whole,feature_size,eye_window,each_window_length);
  
        %% identify segments using multiscale segmentation algorithm
        %% use MJD to detect the boundary
         fprintf(['Identify segments using multiscale segmentation algorithm \n' ]);
         label_temp=MSA(Seq,score_nucletide,cg_content,L,sig2,sele,isbound,upstream,downstream);
   
        %% print out results that consists of genomic signatures,conserved score of each predicted GIs 
        %% predicted GIs and visualization files.
        Printresult(Seq,file_name,feature_size,ks_position_whole,ks_ascending_whole,signa_str,L,label_temp,xunhuan);
      end   
  elseif strcmp(Pro_v,'Assemble-predict')==1 
      % The sequences from the input file are assembled into a sequence according to the
      % order of these sequences in the input file, and MTGIpick will process it as a sequence and predict its GIs.
      fprintf(['Assemble and predict \n' ]);  
      if (iscell(Namee))
		    Name=Namee{1};	
            xh=length(Namee); 
       	    for zuhe=1:xh
                Seq=strcat(Seq,Seqq{zuhe});
            end	
            Seq=upper(Seq);       
        else
	        Name=Namee;
		    Seq=Seqq;
	        Seq=upper(Seq);
        end  
       %% Extracted genomic signature
        N=length(Seq);
        mer=cell(1,1);
        each_window_length=[];
        % Calculate C+G content as genomic signatural
        cg_content=cgcomposition(Seq,'WINDOW',100);
        % Calculate k-mer as genomic signatural
        window=1000;
        slidelen=window;
        %   k=4;
        fprintf(['compute ' num2str(k) '-mer \n']);
        [mer{1} each_window_length signa_str]=cmer(Seq,window,slidelen,k);
   
       %% Calculate the score of each window using the IST-LFS algorithm
        fprintf(['IST-LFS Algorithm' ]);
        score_nucletide=[];ks_position_whole=[];ks_ascending_whole=[];
        [score_nucletide ks_position_whole ks_ascending_whole]=ISTLFS(Seq,mer{1},window,slidelen,sig,ww,interation_whole,feature_size,eye_window,each_window_length);
  
       %% identify segments using multiscale segmentation algorithm
       %% use MJD to detect the boundary
        fprintf(['Identify segments using multiscale segmentation algorithm \n' ]);
        label_temp=MSA(Seq,score_nucletide,cg_content,L,sig2,sele,isbound,upstream,downstream);
   
       %% print out results that consists of genomic signatures,conserved score of each predicted GIs 
       %% predicted GIs and visualization files.
        Printresult(Seq,file_name,feature_size,ks_position_whole,ks_ascending_whole,signa_str,L,label_temp,0);
  end
end
 
