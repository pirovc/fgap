% FGAP: an automated gap closing tool
% Vitor C Piro, Helisson Faoro, Vinicius A Weiss, Maria BR Steffens, Fabio O Pedrosa, Emanuel M Souza and Roberto T Raittz
% BMC Research Notes 2014, 7:371  doi:10.1186/1756-0500-7-371
%
% The MIT License (MIT)
% 
% Copyright (c) 2014 UFPR - Universidade Federal do Paraná (Vitor C. Piro - vitorpiro@gmail.com)
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%

function [] = fgap(varargin)

% Status:
% s1 - Regions were not found in datasets
% s2 - BLAST returned no significant results
% s3 - Can not find significant alignments (pairs)
% s4 - There were no pre-candidates in dataset
% s5 - There were no candidates for closing
% s6 - There were no compatible blast results
%Inicia contador de tempo
tic ();
   
global minScore maxEValue minIdentity contigEndLength edgeTrimLength gapChar maxRemoveLength maxInsertLength ...
    blastAlignParam blastMaxResults threads outputPrefix moreOutput blastPath positiveGap zeroGap negativeGap ...
    plusBasesMoreOutput saveDatasetsPath loadDatasetsPath;

minScore = 25;
maxEValue = 1e-7;
minIdentity = 70;

contigEndLength = 300;
edgeTrimLength = 0;
maxRemoveLength = 500;
maxInsertLength = 500;

positiveGap = 1;
zeroGap = 0;
negativeGap = 0;

gapChar = 'N';
blastPath = '';
blastAlignParam = '1,1,1,-3,15';
blastMaxResults = 200;
threads = 1;

moreOutput = 0;
outputPrefix = 'output_fgap';

%%%%% VERSAO %%%%%%%
global version;
version = '1.8.1';
%%%%%%%%%%%%%%%%%%%%
disp([repmat('-',1,42) 10 9 9 'FGAP v' version 10 repmat('-',1,42) 10]);
%%%%%%%%%%%%%%%%%%%%

%% Declaração  
    %disp(num2str(nargin));
    %disp(char(varargin));
    if(nargin==1)
        showHelp();
        return;
    elseif (nargin<4)
        disp('Not enough input arguments');
        showHelp();
        return;
    elseif (rem(nargin,2)~=0)
        disp('Incorrect number of arguments');
        showHelp();
        return;
    else
        okargs = {'-d','--draft-file',...
                  '-a','--datasets-files',...
                  '-S','--save-datasets',...
                  '-L','--load-datasets',...
                  '-s','--min-score', ...
                  '-e','--max-evalue', ...
                  '-i','--min-identity', ...
                  '-C','--contig-end-length', ...
                  '-T','--edge-trim-length', ...
                  '-R','--max-remove-length', ...
                  '-I','--max-insert-length', ...
                  '-p','--positive-gap', ...
                  '-z','--zero-gap', ...
                  '-g','--negative-gap', ...
                  '-c','--gap-char', ...
                  '-b','--blast-path', ...
                  '-l','--blast-alignment-parameters', ...
                  '-r','--blast-max-results', ...
                  '-t','--threads', ...
                  '-m','--more-output', ...
                  '-o','--output-prefix',...
                  '-h','--help'};
        
        % Caso possua argumentos que não pertencam a lista
        if length(setdiff(varargin(1:2:end),okargs(1:1:end)))>0
            disp('Incorrect arguments');
            showHelp();
            return;
        else   
            for i=1:2:nargin
                if(strcmp(varargin{i},'-d') || strcmp(varargin{i},'--draft-file'))
                    draftFile = varargin{i+1};
                elseif(strcmp(varargin{i},'-a') || strcmp(varargin{i},'--datasets-files'))
                    datasetsFiles = varargin{i+1};
                elseif(strcmp(varargin{i},'-s') || strcmp(varargin{i},'--min-score'))
                    minScore = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-e') || strcmp(varargin{i},'--max-evalue'))
                    maxEValue = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-i') || strcmp(varargin{i},'--min-identity'))
                    minIdentity = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-C') || strcmp(varargin{i},'--contig-end-length'))
                    contigEndLength = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-T') || strcmp(varargin{i},'--edge-trim-length'))
                    edgeTrimLength = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-R') || strcmp(varargin{i},'--max-remove-length'))
                    maxRemoveLength = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-I') || strcmp(varargin{i},'--max-insert-length'))
                    maxInsertLength = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-p') || strcmp(varargin{i},'--positive-gap'))
                    positiveGap = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-z') || strcmp(varargin{i},'--zero-gap'))
                    zeroGap = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-g') || strcmp(varargin{i},'--negative-gap'))
                    negativeGap = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-c') || strcmp(varargin{i},'--gap-char'))
                    gapChar = varargin{i+1};
                elseif(strcmp(varargin{i},'-b') || strcmp(varargin{i},'--blast-path'))
                    blastPath = varargin{i+1};
                elseif(strcmp(varargin{i},'-l') || strcmp(varargin{i},'--blast-alignment-parameters'))
                    blastAlignParam = varargin{i+1};
                elseif(strcmp(varargin{i},'-r') || strcmp(varargin{i},'--blast-max-results'))
                    blastMaxResults = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-t') || strcmp(varargin{i},'--threads'))
                    threads = str2num(varargin{i+1});
                elseif(strcmp(varargin{i},'-m') || strcmp(varargin{i},'--more-output'))
                    moreOutput = str2num(varargin{i+1});    
                elseif(strcmp(varargin{i},'-o') || strcmp(varargin{i},'--output-prefix'))
                    outputPrefix = varargin{i+1};
                end
            end
        end
    end

    % Quantidade de bases dos contigs ends no outputmore    
    plusBasesMoreOutput = contigEndLength;
    
    % Salvar dataset para uso posterior
    saveDatasetsPath = '';
    loadDatasetsPath = '';
%% Verificações

    if ~exist('draftFile','var') || ~exist('datasetsFiles','var')
       disp('Not enough input arguments');
       return;
    elseif ~exist(draftFile,'file')
       disp(['File not found: ' draftFile]);
       return;
    elseif ~isFasta(draftFile)
       disp(['File is not in fasta format: ' draftFile]);
       return;
    else
        if(sum(regexp(datasetsFiles,','))==0 && isempty(datasetsFiles))
           disp(['Wrong input arguments (use commas between datasets)']); 
           return;
        end
        ds = regexp(datasetsFiles,',','split');
        for dt = 1:length(ds) 
            df = char(ds(dt));
            if (length(ds)==1 && strcmp(df(end-8:end),'.datasets'))
                if ~exist([df '.mat'],'file') || ~exist([df '.nhr'],'file') || ~exist([df '.nin'],'file')  || ~exist([df '.nsq'],'file')
                    disp(['File(s) not found: ' df ' (.mat, .nhr, .nin, .nsq)']);
                    return;
                else
                    loadDatasetsPath = df;
                end
            elseif (dt==length(ds) && strcmp(df(end-8:end),'.datasets'))
                [fo, ~, ~] = fileparts(df);
                if ~isempty(fo) 
                    if ~exist(fo,'dir')
                        disp(['Directory not found: ' fo]);
                    else
                        saveDatasetsPath = df;
                    end
                else
                    saveDatasetsPath = df;
                end
            else
                if ~exist(df,'file')
                	disp(['File not found: ' df]);
                    return;
                elseif ~isFasta(df)
                    disp(['File is not in fasta format: ' df]);
                    return; 
                end
            end
        end
    end
    
    % Verifica se pasta de saída existe
    [folder_outputPrefix,~,~] = fileparts(outputPrefix);
    if ~exist(folder_outputPrefix,'dir') && ~isempty(folder_outputPrefix)
       disp('Output folder does not exist');
       return;
    elseif ~isempty(folder_outputPrefix)
        folder_outputPrefix = [folder_outputPrefix '/'];
    end
    
    % Verifica caminho blast
    global makeblastdbAlgorithm;
    [makeblastdbAlgortihmStatus, makeblastdbAlgorithm] = system([blastPath 'makeblastdb -version']);
    if makeblastdbAlgortihmStatus>0
        disp(['Error MAKEBLASTDB path: ' makeblastdbAlgorithm]);
        return; 
    else
        makeblastdbAlgorithm = strtrim(makeblastdbAlgorithm);
        disp(makeblastdbAlgorithm);
    end
    global blastnAlgorithm;
    [blastnAlgortihmStatus, blastnAlgorithm] = system([blastPath 'blastn -version']);
    if blastnAlgortihmStatus>0
        disp(['Error BLASTN path: ' blastnAlgorithm]);
        return; 
    else
        blastnAlgorithm = strtrim(blastnAlgorithm);
        disp(blastnAlgorithm);
    end

%% Inicialização de componentes
    % Para notação cientifica na matriz
    %format('shortG');
    
    % Cria pasta temporária
    global tmp_folder;
    tmp_folder = [folder_outputPrefix 'tmp_fgap/'];
    if exist(tmp_folder,'dir') 
        rmdir(tmp_folder,'s'); 
    end
    [mkdirStatus, mkdirMsg] = mkdir(tmp_folder);
    if mkdirStatus==0
        disp(['Temporary folder ' tmp_folder ' could not be created:  ' mkdirMsg]);
        return; 
    end
    
    % Arquivo de stats
    global stats_file;
    stats_file = [outputPrefix '.stats'];
    if exist(stats_file,'file') 
        delete(stats_file); 
    end
    
    % Estatísticas - totalizadores
    global stats;
    stats = struct('before',[],'after',[],'removedBases',0,'insertedBases',0,'totalGapsClosed',0,'datasetCount',[],'typeCount',[0 0 0]);
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% %% %%%%%%%%%%%%%%%% INICIA FGAP %%%%%%%%%%%%%%%% %% %% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    warning off
    [draft_genome_fasta error_msg] = loadDraftFile(draftFile);
    if isempty(error_msg)
        stats.before = getStats(draft_genome_fasta,0);        
        [datasets_fasta error_msg] = loadDatasetsFiles(datasetsFiles);
        if isempty(error_msg)
            if isempty(buildDatabase(datasets_fasta))
                [draft_genome_fasta error_msg] = identifyGaps(draft_genome_fasta);
                if isempty(error_msg)
                    
                    disp('Starting gap closure ...');
                    
                    pre_fasta = draftFile;
                    out_more_before = [outputPrefix '.before.fasta' ];
                    out_more_after = [outputPrefix '.after.fasta' ];
                    
                    total_char = getTotalChars(draft_genome_fasta);
                    total_char_before = 0;
                    log_round = [];
                    more_before = cell(1); more_after = cell(1);
                    round = 0;
                    while total_char~=0 && total_char_before~=total_char
                        round=round+1;

                        % Função principal 
                        pre_draft_genome_fasta = draft_genome_fasta;
                        [draft_genome_fasta candidates error_msg status_msg] = closeGaps(draft_genome_fasta, datasets_fasta);
                        if isempty(error_msg) && isempty(status_msg)
                            
                            out_fasta = [outputPrefix '_' num2str(round) '.fasta'];
                            out_final = [outputPrefix '.final.fasta'];
                            out_log = [outputPrefix '_' num2str(round) '.log' ];
                        
                            % Gera log
                            [log] = generateLog(candidates,pre_draft_genome_fasta,datasets_fasta,pre_fasta);
                            
                            % Salva total de gaps anteriores e atual
                            total_char_before = total_char;
                            total_char = getTotalChars(draft_genome_fasta);

                            stats.totalGapsClosed = stats.totalGapsClosed + (total_char_before-total_char);

                            % Mostra informações na tela
                            gapsclosed_round = ['Round ' num2str(round) ': ' num2str(total_char_before-total_char) ' gaps'];
                            disp(gapsclosed_round);

                            % Salva dados da rodada no arquivo de stats
                            log_round(round).log = [gapsclosed_round 10 ' Output log: ' hidePath(out_log) 10 ' Output fasta: ' hidePath(out_fasta) 10 ' After round ' num2str(round) ':' 10 '  ' getStats(draft_genome_fasta,1) 10];

                            % Escreve Log da rodada
                            writeLog(out_log,log);
                            % Escreve fasta da rodada
                            writeFasta(out_fasta, draft_genome_fasta);
                            if(moreOutput)
                                [more_before{round} more_after{round}] = generateMoreOutput(candidates,pre_draft_genome_fasta);
                            end
                            
                            pre_fasta = [outputPrefix '_' num2str(round) '.fasta'];
                        elseif ~isempty(status_msg)
                            disp(status_msg);
                            total_char_before = total_char;
                        elseif ~isempty(error_msg)
                            break;
                        end      
                    end
                    stats.after = getStats(draft_genome_fasta,0);
                    if(moreOutput)
                        writeMoreOutput(more_before,out_more_before,more_after,out_more_after);
                    end
                    %Escreve arquivo final
                    if(exist('out_fasta','var'))
                        copyfile(out_fasta,out_final);
                    end
                    elapsed = ['Elapsed time is ' num2str(toc ()) ' seconds'];
                    writeStats(draftFile,datasets_fasta,log_round,elapsed);    
                end
            end
        end
    end
    if ~isempty(error_msg)
       disp(error_msg);
    end
    disp('Cleaning files ...');
    rmdir(tmp_folder,'s');
    if exist('elapsed','var')
        disp(elapsed);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% %% %%%%%%%%%%%%%%%%%% FUNCOES %%%%%%%%%%%%%%%%%% %% %% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [out_draft candidates error_msg status_msg] = closeGaps(draft_genome_fasta, datasets_fasta)
global contigEndLength edgeTrimLength minScore maxEValue minIdentity blastMaxResults blastAlignParam threads tmp_folder blastPath stats saveDatasetsPath loadDatasetsPath;
out_draft = draft_genome_fasta;
candidates = [];
error_msg = '';
status_msg = '';

%% Procura Chars válidos no Draft
    for i = 1:length(draft_genome_fasta)       
        if ~isempty(draft_genome_fasta(i).Pos)
            % Marca gaps que já tentaram ser fechados (1) com -1
            draft_genome_fasta(i).Pos(draft_genome_fasta(i).Pos(:,4)==1,4) = -1;
            % Gaps que ainda nao foram fechados (0)
            new_chars = draft_genome_fasta(i).Pos(draft_genome_fasta(i).Pos(:,4)==0,1:3);

            if ~isempty(new_chars)
                % Pega apenas chars que estão distantes um do outros (areas de contigEndLength nao se encontram)
                cnt=1;
                draft_genome_fasta(i).Pos(draft_genome_fasta(i).Pos(:,3)==new_chars(1,3),4) = 1;
                for j=2:length(new_chars(:,1))
                    % Caso tenha diferenca de dois lados do ultimo gap valido, marca como candidato (1)
                    if (new_chars(j,1) - new_chars(j-cnt,2)) >= (edgeTrimLength*2)+(contigEndLength*2)
                        draft_genome_fasta(i).Pos(draft_genome_fasta(i).Pos(:,3)==new_chars(j,3),4) = 1;
                        cnt = 1;
                    else
                        cnt = cnt + 1;
                    end
                end
            end
        end
    end
    
%% Gera regiões de ancoragem 
    cont = 0;
    draft_contig = struct('draft_char',struct());
    for i = 1:length(draft_genome_fasta)
       if ~isempty(draft_genome_fasta(i).Pos)
           % Verifica apenas chars candidatos (1)
           cand_chars = draft_genome_fasta(i).Pos(draft_genome_fasta(i).Pos(:,4)==1,1:3);
           if ~isempty(cand_chars)    
               for j = 1:length(cand_chars(:,1))

                   % Caso edges ultrapassem regiao possível (sem deixar área para o contigEndLength)
                   if edgeTrimLength >= cand_chars(j,1)-1 || (cand_chars(j,2) + edgeTrimLength) >= length(draft_genome_fasta(i).Sequence)
                      continue; 
                   end
                          
                   % Posição inicial e final de ancoragem
                   seq_start = cand_chars(j,1)-(contigEndLength+edgeTrimLength);
                   seq_end = cand_chars(j,2)+(contigEndLength+edgeTrimLength);
                   
                   % Caso não possua região de aconragem suficiente (excede arquivo)              
                   if seq_start<1
                       seq_start_len = cand_chars(j,1)-edgeTrimLength-1;
                       seq_start = 1;
                   else
                       seq_start_len = contigEndLength;
                   end 
                   if seq_end>length(draft_genome_fasta(i).Sequence) 
                       seq_end_len = (length(draft_genome_fasta(i).Sequence)-cand_chars(j,2))-edgeTrimLength;
                       seq_end = length(draft_genome_fasta(i).Sequence);
                   else
                       seq_end_len = contigEndLength;
                   end
                   
                   % Para recuperação de dados
                   draft_contig(i).draft_char(cand_chars(j,3)).Header = ['draft_' num2str(i) '_' num2str(cand_chars(j,3))];
                   draft_contig(i).draft_char(cand_chars(j,3)).Sequence = draft_genome_fasta(i).Sequence(seq_start:seq_end);

                   % Para arquivo fasta, grava lado A (montante) e B (jusante) em relação ao char
                   cont = cont + 1;
                   draft_regions(cont).Header = ['draft_' num2str(i) '_' num2str(cand_chars(j,3)) '_A'];
                   draft_regions(cont).Sequence = draft_genome_fasta(i).Sequence(seq_start:seq_start+seq_start_len-1);
                   cont = cont + 1;
                   draft_regions(cont).Header = ['draft_' num2str(i) '_' num2str(cand_chars(j,3)) '_B'];
                   draft_regions(cont).Sequence = draft_genome_fasta(i).Sequence(seq_end-seq_end_len+1:seq_end);
                   
               end
           end 
       end
    end
    if ~exist('draft_regions','var')
        status_msg = 'There are no more possible gaps to be closed (s1)';
        return; 
    end
    
%% BLAST
    % Salva arquivo
    fastawrite2([tmp_folder 'draft_regions.fasta'],draft_regions);
    
    % Verifica se está carregando de arquivo específico ou pasta temporária
    % (se está salvando também está em pasta diferente)
    if ~isempty(loadDatasetsPath)
        dboutputfolder = loadDatasetsPath;
    elseif ~isempty(saveDatasetsPath)
        dboutputfolder = saveDatasetsPath;
    else
        dboutputfolder = [tmp_folder 'datasets'];
    end
    
    blastAlignParamSplit = regexp(blastAlignParam,',','split');
    % Executa BLAST - blastn
    blastn = ['-query ' tmp_folder 'draft_regions.fasta' ...
                ' -db ' dboutputfolder ...
                ' -task blastn ' ... 
                ' -min_raw_gapped_score ' num2str(minScore) ...
                ' -evalue ' num2str(maxEValue) ...
                ' -perc_identity ' num2str(minIdentity) ...
                ' -word_size ' char(blastAlignParamSplit(5)) ...
                ' -gapopen ' char(blastAlignParamSplit(1)) ...
                ' -gapextend ' char(blastAlignParamSplit(2)) ...
                ' -penalty ' char(blastAlignParamSplit(4)) ...
                ' -reward ' char(blastAlignParamSplit(3)) ...
                ' -max_target_seqs ' num2str(blastMaxResults) ...
                ' -out ' [tmp_folder 'blastn.out'] ...
                ' -num_threads ' num2str(threads) ...%' -dust no ' % fechou mais mas piorou a validacao ...
                ' -outfmt "6 qseqid sseqid score bitscore evalue pident nident mismatch sstrand qstart qseq qend sstart sseq send qlen length qcovhsp"'
                ];
                %qseqid sseqid score bitscore evalue pident nident mismatch sstrand qstart qseq qend sstart sseq send qlen
                %1 Query Seq-id
                %2 Subject Seq-id
                %3 Raw score
                %4 Bit score
                %5 Expect value
                %6 Percentage of identical matches
                %7 Number of identical matches
                %8 Number of mismatches
                %9 Subject Strand
                %10 Start of alignment in query
                %11 Aligned part of query sequence
                %12 End of alignment in query
                %13 Start of alignment in subject
                %14 Aligned part of subject sequence
                %15 End of alignment in subject
                %16 Query sequence length
                %17 Alignment length
                %18 Query Coverage Per HSP

    [blastnStatus, ~] = system([blastPath 'blastn ' blastn]);
    if blastnStatus>0
        error_msg = ['Error BLASTN: ' blastnStatus];
        return; 
    end

    fid = fopen([tmp_folder 'blastn.out'],'r');  % Open text file
    blastnResult = textscan(fid,'%s','delimiter','\n','bufsize',(contigEndLength*2)+1000);
    fclose(fid);
    
    if isempty(blastnResult) 
        status_msg = 'There are no more possible gaps to be closed (s2)';
        return; 
    end
    
    cont = 0;
    old_draft_head = '';
    for i=1:length(blastnResult{1})
        parts = regexp(blastnResult{1}{i},'\t','split');
       
        if ~strcmp(old_draft_head,parts(1))
            cont = cont + 1;
            posi = 1;
        else
            posi = posi + 1;
        end
        blast_result{cont}{posi} = parts;
        old_draft_head=parts(1);
    end
    
    if ~exist('blast_result','var')
        status_msg = 'There are no more possible gaps to be closed (s6)';
        return; 
    end    
 %% Filtra resultados do BLAST (retira lados orfãos)  
    i=0;
    cont = 0;
    % Apenas mantem resultados que possuem outro lado com resultados também
    while i <= length(blast_result) - 2
        if i <= length(blast_result) - 2
            i = i + 1;
            A = blast_result{i};
            A_q = A{1}{1};
            i = i + 1;
            B = blast_result{i};
            B_q = B{1}{1};
            if ~strcmp(A_q(1:end-2),B_q(1:end-2))
                i = i - 1;
            else
                cont = cont + 1;
                blast_result_filter{cont} = A;
                cont = cont + 1;
                blast_result_filter{cont} = B;
            end
        end
    end
    if ~exist('blast_result_filter','var')
        status_msg = 'There are no more possible gaps to be closed (s3)';
        return; 
    end
    
%% Loop nos resultados para pegar candidatos
    i=0;
    cont = 0;
    while i < length(blast_result_filter)
        
        % Pega HSPs dos dois lados
        i = i + 1;
        A = blast_result_filter{i};
        A_hsps = getHSPs(A,i);

        i = i + 1;
        B = blast_result_filter{i};
        B_hsps = getHSPs(B,i);

        % Caso possua HSPs compativeis com os parametros
        if (~isempty(A_hsps) && ~isempty(B_hsps))

            % Loop nos resultados do A procurando no B
            for j=1:length(A_hsps(:,1))
               vA = []; vB = []; vB_list=[];
               vA = A_hsps(j,:);
               % Compara coluna 1 (file) 2 (ctg) e 9 (strand) - se bateu no
               % mesmo arquivo e contig do dataset e strand
               vB_list = B_hsps(A_hsps(j,1)==B_hsps(:,1) & A_hsps(j,2)==B_hsps(:,2) & A_hsps(j,9)==B_hsps(:,9),:);

               % Caso possua hits no B_hsps(vB_list) referente ao A_hsps
               if (~isempty(vB_list))
                   
                   % Adiciona todos os possiveis candidatos do lado B que tenham referencia no lado A
                   for k=1:length(vB_list(:,1))
                       
                       % Pega ocorrencia de vB
                       vB = vB_list(k,:);

                       cont = cont +1;

                       id = regexp(A{1}{1}, '_', 'split');
                       id_draft_scaf = str2num(id{2});
                       id_char = str2num(id{3});

                       % Grava pré-candidatos
                       pre_candidates(cont).id_DtsetFile = vA(1);
                       pre_candidates(cont).id_DtsetContig = vA(2);
                       pre_candidates(cont).id_DraftScaf = id_draft_scaf;
                       pre_candidates(cont).id_DraftChar = id_char;

                       pre_candidates(cont).rawScore = [vA(16) vB(16)];
                       pre_candidates(cont).bitScore = [vA(3) vB(3)];
                       pre_candidates(cont).eValue = [vA(4) vB(4)];
                       pre_candidates(cont).identitiesPercent = [vA(12) vB(12)];
                       pre_candidates(cont).identitiesPossible = [vA(13) vB(13)];
                       pre_candidates(cont).identitiesMatch = [vA(14) vB(14)];
                       pre_candidates(cont).queryIndices = [vA(5) vA(6) vB(5) vB(6)];
                       pre_candidates(cont).subjectIndices = [vA(7) vA(8) vB(7) vB(8)];
                       pre_candidates(cont).alignmentBlastQuery = [blast_result_filter{vA(10)}{vA(11)}(11) ; blast_result_filter{vB(10)}{vB(11)}(11)];
                       pre_candidates(cont).alignmentBlastSubject = [blast_result_filter{vA(10)}{vA(11)}(14) ; blast_result_filter{vB(10)}{vB(11)}(14)];
                       pre_candidates(cont).strand = [vA(9) vB(9)];
                       pre_candidates(cont).lengthcontigEndLengths = [vA(15) vB(15)];
                       pre_candidates(cont).queryCoverage = [vA(17) vB(17)];
                       
                   end
               end
            end
        end
    end
    if ~exist('pre_candidates','var')
        status_msg = 'There are no more possible gaps to be closed (s4)';
        return; 
    end
    
%% Verifica candidatos
    pre_candidates = verifyPreCandidates(pre_candidates,draft_genome_fasta, datasets_fasta);
    if ~isempty(pre_candidates)
        candidates = selectCandidates(pre_candidates,datasets_fasta);
    end
    if isempty(candidates)
        status_msg = 'There are no more possible gaps to be closed (s5)';
        return; 
    end
    
%% Adiciona informações ao candidato
    candidates = addCandidateData(candidates,draft_genome_fasta,datasets_fasta);

%% Gera novo draft com gap fechado
    new_draft = draft_genome_fasta;
    for i=1:length(candidates)

        % NÃO USAR DA STRUCT - Precisa gerar a cada rodada pois dados da tabela mudam
        [remove_start remove_end remove_seq] = getRemovedSeqDraft(candidates(i),new_draft);
        
        % Caso tenha fechado mais de um gap, verifica se foi antes ou
        % depois e quantos para cada lado (para atualizar a tabela de Pos)
        closed_chars = getChars(remove_seq);  
        closed_chars = length(closed_chars(:,1));

        if closed_chars>1
            [n_start n_end] = getPosN(candidates(i),new_draft);
            rem_before = n_start - remove_start;
            rem_after = remove_end - n_end;
            closed_chars_before = size(getChars(remove_seq(1:rem_before)),1);
            closed_chars_after = size(getChars(remove_seq(end-rem_after+1:end)),1);

            % Verifica apenas chars disponiveis para fechamento (~=2)
            open_chars = new_draft(candidates(i).id_DraftScaf).Pos(new_draft(candidates(i).id_DraftScaf).Pos(:,4)~=2,:);
            % Identifica posicao relativa
            pos_open_chars = find(open_chars(:,3)==candidates(i).id_DraftChar);
            % Encontra quais chars antes e depois foram fechados
            chars_bef_aft = open_chars(pos_open_chars-closed_chars_before:pos_open_chars+closed_chars_after,3);
            index_chars_closed = ismember(new_draft(candidates(i).id_DraftScaf).Pos(:,3),chars_bef_aft);
        else
            index_chars_closed = candidates(i).id_DraftChar;
        end     
        
        % Tamanho anterior para comparacao de bases
        old_length = length(new_draft(candidates(i).id_DraftScaf).Sequence);

        % Insere nova sequencia no draft
        new_draft(candidates(i).id_DraftScaf).Sequence = [new_draft(candidates(i).id_DraftScaf).Sequence(1:remove_start-1) lower(candidates(i).insertedSeqDtset) new_draft(candidates(i).id_DraftScaf).Sequence(remove_end+1:end)]; 

        new_length = length(new_draft(candidates(i).id_DraftScaf).Sequence);
        % Move os outros gaps do contig para a posição correta apos inserir nova sequencia
        % Calcula diferenca de area
        new_area_size = new_length - old_length;
        % Pega todos os chars abaixo do gap fechado
        pos_chars_below = new_draft(candidates(i).id_DraftScaf).Pos(:,3)>candidates(i).id_DraftChar;
        % Altera posicao dos chars abaixo
        new_draft(candidates(i).id_DraftScaf).Pos(pos_chars_below,1:2) = new_draft(candidates(i).id_DraftScaf).Pos(pos_chars_below,1:2)+new_area_size;
        % Marca gaps como fechado (2) no Pos (verifica se nao fechou mais de um por vez)
        new_draft(candidates(i).id_DraftScaf).Pos(index_chars_closed,4) = 2;

        stats.datasetCount(candidates(i).id_DtsetFile) = stats.datasetCount(candidates(i).id_DtsetFile) + closed_chars;
        stats.typeCount(candidates(i).gapType+2) = stats.typeCount(candidates(i).gapType+2) + closed_chars;
        stats.removedBases = stats.removedBases + (remove_end - remove_start + 1);
        stats.insertedBases = stats.insertedBases + length(candidates(i).insertedSeqDtset);
    end
    out_draft = new_draft;  
end
%%
function [cand] = addCandidateData(cand,draft_genome_fasta,datasets_fasta)
    for i=1:length(cand)
        % Salva dados de inserção
        [insert_start insert_end insert_seq] = getInsertedSeqDtset(cand(i),datasets_fasta(cand(i).id_DtsetFile).fasta);
        cand(i).insertedSeqDtset = insert_seq;
        cand(i).insertedPosDtset = [insert_start insert_end];

        % Salva dados de remoção
        [remove_start remove_end remove_seq] = getRemovedSeqDraft(cand(i),draft_genome_fasta);
        cand(i).removedSeqDraft = remove_seq;
        cand(i).removedPosDraft = [remove_start remove_end];
    end
end
%%
function [log_complete] = generateLog(candidates,pre_draft_genome_fasta,datasets_fasta,pre_fasta)

    for c=1:length(candidates)

        % Query(draft/Cends)  1
        % ||||||||||||||||||  2
        % Subject (datasets)  3
        vAl = (candidates(c).alignmentBlastQuery{1}==candidates(c).alignmentBlastSubject{1})*'|';
        vAl(vAl==0) = 32;
        blastA = [candidates(c).alignmentBlastQuery{1}; vAl ;candidates(c).alignmentBlastSubject{1}];
        vBl = (candidates(c).alignmentBlastQuery{2}==candidates(c).alignmentBlastSubject{2})*'|';
        vBl(vBl==0) = 32;
        blastB = [candidates(c).alignmentBlastQuery{2}; vBl ;candidates(c).alignmentBlastSubject{2}];

        dtset_min = min(candidates(c).subjectIndices);
        dtset_max = max(candidates(c).subjectIndices);

        openSubjectA = 0;
        if(candidates(c).gapType==1)
            % 'empurra' começo do B quando existe open gap no lado A
            openSubjectA = length(regexp(blastA(3,:),'-')) ;
        end
        if(candidates(c).strand(1)==1)
            A_st = candidates(c).subjectIndices(1) - dtset_min + 1;
            B_st = candidates(c).subjectIndices(3) - dtset_min + 1 + openSubjectA;
        elseif(candidates(c).strand(1)==2)
            A_st = dtset_max - candidates(c).subjectIndices(1) + 1;
            B_st = dtset_max - candidates(c).subjectIndices(3) + 1 + openSubjectA;
        end
        A_en = A_st + length(candidates(c).alignmentBlastSubject{1}) - 1;
        B_en = B_st + length(candidates(c).alignmentBlastSubject{2}) - 1;

        % ajusta visualização de gaps sobrepostos
        if candidates(c).gapType==-1 && B_en > A_st
            if A_st > B_st
                openStart = length(regexp(candidates(c).alignmentBlastSubject{2}(1:A_st-1),'-')) ;
                A_st = A_st + openStart;
                A_en = A_en + openStart;
            elseif B_st > A_st
                openStart = length(regexp(candidates(c).alignmentBlastSubject{1}(1:B_st-1),'-')) ;
                B_st = B_st + openStart;
                B_en = B_en + openStart;
            end
        end
        
        % Tamamanho geral do log
        dtset_len = max([A_en B_en]);
        
        % Contig End 3' Sequencia
        a1A = blanks(dtset_len);
        a1A(A_st:A_en) = blastA(1,:);

        % Contig End 5' Sequencia
        a1B = blanks(dtset_len);
        a1B(B_st:B_en) = blastB(1,:);

        % Contig End 3' Alinhamento
        a2A = blanks(dtset_len);
        a2A(A_st:A_en) = blastA(2,:);
        
        % Contig End 5' Alinhamento
        a2B = blanks(dtset_len);
        a2B(B_st:B_en) = blastB(2,:);

        % Subject - Dataset Sequencia
        a3A = blanks(dtset_len);
        a3A(A_st:A_en) = blastA(3,:);
        
        a3B = blanks(dtset_len);
        a3B(B_st:B_en) = blastB(3,:);

        % Caso seja gap positivo, região inserida   
        insert_seq = '';
        if(candidates(c).gapType==1)
            insert_seq = lower(candidates(c).insertedSeqDtset);
            a3A(A_en+1:B_st-1) = insert_seq;
            a3B(A_en+1:B_st-1) = insert_seq;
        end

        % Posições
        [queryA_start queryA_end queryB_start queryB_end] = getRealQueryPos(candidates(c),pre_draft_genome_fasta);
        subjectA_start = candidates(c).subjectIndices(1);
        subjectA_end = candidates(c).subjectIndices(2);
        subjectB_start = candidates(c).subjectIndices(3);
        subjectB_end = candidates(c).subjectIndices(4);
        if(candidates(c).strand(1)==1 && candidates(c).gapType==1)
            subjectA_end = subjectA_end + length(insert_seq);
            subjectB_start = subjectB_start - length(insert_seq);
        elseif(candidates(c).strand(1)==2 && candidates(c).gapType==1)
            subjectA_end = subjectA_end - length(insert_seq);
            subjectB_start = subjectB_start + length(insert_seq);
        end
        
        gap_id = getGapId(candidates(c));
                
        % INICIO - escreve log em cell
        log = cell(1);
        l = 1;

        log{l} = ['Gap ID: ' gap_id];
        l=l+1;
        log{l} = ['Gap Type: ' gapTypeName(candidates(c).gapType) 10];
        l=l+1;
        log{l} = ['Draft file: ' hidePath(pre_fasta) ' (' pre_draft_genome_fasta(candidates(c).id_DraftScaf).Header ')'];
        l=l+1;
        log{l} = ['DtSet file: ' hidePath(datasets_fasta(candidates(c).id_DtsetFile).file) ' (' datasets_fasta(candidates(c).id_DtsetFile).fasta(candidates(c).id_DtsetContig).Header ')' 10];
        l=l+1;
        log{l} = ['Contig end 3'':'];
        l=l+1;
        log{l} = [' Bit Score: ' num2str(candidates(c).bitScore(1)) ' bits (' num2str(candidates(c).rawScore(1)) ')'];
        l=l+1;
        log{l} = [' E-Value: '  num2str(candidates(c).eValue(1))];
        l=l+1;
        log{l} = [' Identity: ' num2str(candidates(c).identitiesMatch(1)) '/' num2str(candidates(c).identitiesPossible(1)) ' (' num2str(candidates(c).identitiesPercent(1)) '%)'];
        l=l+1;
        log{l} = [' Query Cov.: ' num2str(candidates(c).queryCoverage(1)) '%'];
        l=l+1;
        log{l} = ['Contig end 5'':'];
        l=l+1;
        log{l} = [' Bit Score: ' num2str(candidates(c).bitScore(2)) ' bits (' num2str(candidates(c).rawScore(2)) ')'];
        l=l+1;
        log{l} = [' E-Value: '  num2str(candidates(c).eValue(2))];
        l=l+1;
        log{l} = [' Identity: ' num2str(candidates(c).identitiesMatch(2)) '/' num2str(candidates(c).identitiesPossible(2)) ' (' num2str(candidates(c).identitiesPercent(2)) '%)'];
        l=l+1;
        log{l} = [' Query Cov.: ' num2str(candidates(c).queryCoverage(2)) '%'];
        l=l+1;
        log{l} = ['Strand (CEnds/Dtset): ' verifyStrandRev(candidates(c).strand(1)) 10];
        l=l+1;

        %len_res = dtset_len;
        len_res = 60;
        len_max_number = length(num2str(max([queryA_start queryB_start subjectA_start subjectB_start])));
        pos_qA = ''; pos_qB = '';pos_sA='';pos_sB='';
        for x=1:len_res:dtset_len
            if(x+len_res>dtset_len)
                e = dtset_len;
                pos_qA = [' ' num2str(queryA_end,'%.f')];
                pos_qB = [' ' num2str(queryB_end,'%.f')];
                pos_sA = [' ' num2str(subjectA_end,'%.f')];
                pos_sB = [' ' num2str(subjectB_end,'%.f')];
            else
                e = x+len_res-1;
            end

            if x==1
                pre_qA = sprintf(['%' num2str(len_max_number) 'd'],queryA_start);
                pre_qB = sprintf(['%' num2str(len_max_number) 'd'],queryB_start);
                pre_sA = sprintf(['%' num2str(len_max_number) 'd'],subjectA_start);
                pre_sB = sprintf(['%' num2str(len_max_number) 'd'],subjectB_start);
            else
                pre_qA = blanks(len_max_number);
                pre_qB = blanks(len_max_number);
                pre_sA = blanks(len_max_number);
                pre_sB = blanks(len_max_number);
            end
            pre_al = blanks(len_max_number);

            log{l} = ['CEnd3 ' pre_qA ' ' a1A(x:e) pos_qA];
            l=l+1;
            log{l} = ['      ' pre_al ' ' a2A(x:e)];
            l=l+1;

            log{l} = ['Dtset ' pre_sA ' ' a3A(x:e) pos_sA];
            l=l+1;
            log{l} = ['Dtset ' pre_sB ' ' a3B(x:e) pos_sB];
            l=l+1;
            
            log{l} = ['      ' pre_al ' ' a2B(x:e)];
            l=l+1;
            log{l} = ['CEnd5 ' pre_qB ' ' a1B(x:e) pos_qB];
            l=l+1;

            log{l} = [10];
            l=l+1;
        end
        
        log{l} = ['Removed sequence (' num2str(length(candidates(c).removedSeqDraft)) 'bp):'];
        l=l+1;
        len_max_number = length(num2str(candidates(c).removedPosDraft(1)));
        remove_seq_p = regexp(candidates(c).removedSeqDraft, ['\w{1,' num2str(len_res) '}'], 'match');
        log{l} = [num2str(candidates(c).removedPosDraft(1),'%.f') ' ' remove_seq_p{1}]; 
        for p=2:length(remove_seq_p)
            log{l} = [log{l} 10 blanks(len_max_number) ' ' remove_seq_p{p}];      
        end
        log{l} = [log{l} ' ' num2str(candidates(c).removedPosDraft(2),'%.f')];
        l=l+1;
        
        if(candidates(c).gapType==1)
            log{l} = ['Inserted sequence (' num2str(length(insert_seq)) 'bp):'];
            l=l+1;
            if(candidates(c).strand==1)
               is = candidates(c).insertedPosDtset(1);
               ie = candidates(c).insertedPosDtset(2);
            else
               is = candidates(c).insertedPosDtset(2);
               ie = candidates(c).insertedPosDtset(1);
            end
            len_max_number = length(num2str(is));
            insert_seq_p = regexp(insert_seq, ['\w{1,' num2str(len_res) '}'], 'match');
            log{l} = [num2str(is,'%.f') ' ' insert_seq_p{1}]; 
            for p=2:length(insert_seq_p)
                log{l} = [log{l} 10 blanks(len_max_number) ' ' insert_seq_p{p}];
            end
            log{l} = [log{l} ' ' num2str(ie,'%.f')];
            l=l+1;
        end
        log{l} = [10 repmat('-',1,80)];
        l=l+1;
        
        log_complete{c} = log;
    end
end
%%
function [hsps_out] = getHSPs(blast_result,id_blast)
global minScore;
% Output:
%1 id_file
%2 id_contig
%3 Score
%4 eValue
%5 QueryIndices start
%6 QueryIndices end
%7 SubjectIndices start
%8 SubjectIndices end
%9 Strand (1 - Plus/Plus, 2 - Plus/Minus)
%10 id blast_result_filter
%11 id HSP
%12 Identity Percent
%13 Identity Possible (Alignment length)
%14 Identity Match
%15 Query Length
%16 Bit Score
%17 Query Coverage Per HSP

    hsps_out = [];
    for j=1:length(blast_result)
        id = regexp(char(blast_result{j}(2)), '_', 'split');
        id_dataset = str2num(id{2});
        id_contig = str2num(id{3});

        %if (str2num(blast_result{j}{4}) >= minScore)
            hsps_out = [hsps_out; ... 
                id_dataset ...
                id_contig ...
                str2num(blast_result{j}{4}) ...
                str2num(blast_result{j}{5}) ...
                str2num(blast_result{j}{10}) ...
                str2num(blast_result{j}{12}) ...
                str2num(blast_result{j}{13}) ...
                str2num(blast_result{j}{15}) ...
                verifyStrand(blast_result{j}{9}) ...
                id_blast j ...
                str2num(blast_result{j}{6}) ...
                str2num(blast_result{j}{17}) ...
                str2num(blast_result{j}{7}) ...
                str2num(blast_result{j}{16}) ...
                str2num(blast_result{j}{3}) ...
                str2num(blast_result{j}{18}) ...
                ];
        %end

    end
end
%%
function [pre_cand] = verifyPreCandidates(pre_candidates,draft_genome_fasta,datasets_fasta)
global gapChar maxInsertLength maxRemoveLength positiveGap zeroGap negativeGap;

    cont = 0;
    for i=1:length(pre_candidates)
        
        [insert_start insert_end] = getInsertedPosDtset(pre_candidates(i));

        % Tipo do gap
        if (insert_end-insert_start+1)==0
            % ZERO GAP
            pre_candidates(i).gapType = 0;
            if(~zeroGap) continue; end
        elseif (insert_end-insert_start+1)<0
            % NEGATIVE GAP
            pre_candidates(i).gapType = -1;
            if(~negativeGap) continue; end
        else
            % POSITIVE GAP
            pre_candidates(i).gapType = 1;
            if(~positiveGap) continue; end
        end
        
        % Verifica se nova regiao inserida é menor ou igual ao permitido
        if (insert_end-insert_start+1)>maxInsertLength
            continue;
        end 
        
        [remove_start remove_end] = getRemovedPosDraft(pre_candidates(i),draft_genome_fasta);
        if (remove_end-remove_start+1) > maxRemoveLength 
           continue;
        end

        if(pre_candidates(i).gapType==1)
            [~,~,insert_seq] = getInsertedSeqDtset(pre_candidates(i),datasets_fasta(pre_candidates(i).id_DtsetFile).fasta);
            % Verifica se possui o char na nova região (fecha gap com outro gap)
            if sum(insert_seq==gapChar)>=1
                continue;
            end 
        end
        
        cont = cont + 1;
        pre_cand(cont) = pre_candidates(i);
    end
    if ~exist('pre_cand','var')
        pre_cand = [];
    end
end
%%
function [cand] = selectCandidates(pre_candidates,datasets_fasta)
    pc = zeros(length(pre_candidates),8);
    for i=1:length(pre_candidates)
        p = pre_candidates(i);
        pc(i,:) = [i ...
                  p.id_DraftScaf ...
                  p.id_DraftChar ...
                  sum(p.rawScore) ...
                  sum(p.queryCoverage) ...
                  sum(p.identitiesPercent) ...
                  sum(p.eValue) ...
                  length(datasets_fasta(p.id_DtsetFile).fasta(p.id_DtsetContig).Sequence)];
    end
    pc = sortrows(pc,[2 3 -4 -5 -6 7 -8]);
    [~, idx, ~] = unique(pc(:,2:3),'rows','first');
    pc = pc(idx,:);
    
    cand = pre_candidates(pc(:,1));
    if ~exist('cand','var')
        cand = [];
    end
end

%%
function [stats] = getStats(fasta,inline)
global gapChar;
    if(inline==1)
        separator = 9;
    else
        separator = 10;
    end    
    c=0;g=0;chars=0;ls=[];
    for i=1:length(fasta)
        seq = upper(fasta(i).Sequence);
        c = c + sum(seq=='C');
        g = g + sum(seq=='G');
        chars = chars + sum(seq==gapChar);
        ls = [ls; length(seq)];
    end
    
    lseq = sum(ls);
    gc = (c+g) / lseq;
    
    ls = sortrows(ls,-1);
    aux_sum = 0;n50 = 0;
    h = lseq/2;
    for x=1:length(ls)
        aux_sum = aux_sum + ls(x);
        if(aux_sum >= h)
            n50 = ls(x);
            break;
        end
    end
        
    stats = [' Gaps: ' num2str(getTotalChars(fasta),'%.f') separator ' Sequences: ' num2str(length(fasta),'%.f') separator ' Length: ' num2str(lseq,'%.f') 'bp ' separator ' GC: ' num2str(gc*100) '%' separator ' N50: ' num2str(n50,'%.f') separator ' Min: ' num2str(min(ls),'%.f') separator ' Max: ' num2str(max(ls),'%.f') separator ' ' gapChar 's: ' num2str(chars,'%.f')];

end
%%
function [s] = verifyStrand(strand)
    % Strand (1 - Plus/Plus, 2 - Plus/Minus)
    if strcmp(strand,'plus')
        s=1;
    elseif strcmp(strand,'minus')
        s=2;
    end
end
%%
function [s] = verifyStrandRev(strand)
    % Strand (1 - Plus/Plus, 2 - Plus/Minus)
    if strand==1
        s='Plus/Plus';
    elseif strand==2
        s='Plus/Minus';
    end
end
%%
function [sequence_fasta] = fastaUpper(sequence_fasta)
    for i=1:length(sequence_fasta)
        sequence_fasta(i).Sequence = upper(sequence_fasta(i).Sequence);
    end
end
%%
function [char_list] = getChars(sequence)
global gapChar;
    char_list = [];
    if ~isempty(sequence)
        [sn,en] = regexpi(sequence,[gapChar '[' gapChar ']*']);
        char_list = [sn' en'];
    end
end
%%
function [total] = getTotalChars(sequence_fasta)
    total = 0;
    for i=1:length(sequence_fasta)
        ch = getChars(sequence_fasta(i).Sequence);
        if ~isempty(ch)
            total = total + length(ch(:,1));
        end
    end
end
%%
function [n_start n_end] = getPosN(cand,seq)
    n_start=[];n_end=[];
    N = seq(cand.id_DraftScaf).Pos(cand.id_DraftChar,1:2);
    n_start = N(1);
    n_end = N(2);
end
%%
function [remove_start remove_end] = getRemovedPosDraft(cand,seq)
    [queryA_start queryA_end queryB_start queryB_end] = getRealQueryPos(cand,seq);
    remove_start = queryA_end + 1;
    remove_end = queryB_start - 1;
    
    % Ajusta caso possua alinhamento negativo
    if(cand.gapType==-1)
        openQueryA= 0;
        openSubjectA = 0;
        
        % Identifica posições no alinhamento em relação ao lado A
        dtsetA_min = min(cand.subjectIndices(1:2));
        dtsetA_max = max(cand.subjectIndices(1:2));
        if(cand.strand(1)==1)
            A_st = cand.subjectIndices(1) - dtsetA_min + 1;
            B_st = cand.subjectIndices(3) - dtsetA_min + 1;
        elseif(cand.strand(1)==2)
            A_st = dtsetA_max - cand.subjectIndices(1) + 1;
            B_st = dtsetA_max - cand.subjectIndices(3) + 1;
        end
        A_en = A_st + length(cand.alignmentBlastSubject{1}) -1;
        B_en = B_st + length(cand.alignmentBlastSubject{2}) -1;
        % Apenas quando há sobreposição
        if B_en > A_st
            % Conta número de opengaps na sobreposição
            p = sort([A_st A_en B_st B_en]);
            openQueryA = length(regexp(cand.alignmentBlastQuery{1}(p(2):p(3)),'-'));
            openSubjectA = length(regexp(cand.alignmentBlastSubject{1}(p(2):p(3)),'-'));
        end

        remove_start = remove_start - (abs(cand.subjectIndices(2) - cand.subjectIndices(3)) + 1) + openQueryA - openSubjectA;
        if remove_start<1 remove_start=1; end
    end
end
%%
function [remove_start remove_end remove_seq] = getRemovedSeqDraft(cand,seq)
    [remove_start remove_end] = getRemovedPosDraft(cand,seq);
    if remove_start>0 && remove_end>0 && remove_end>=remove_start
        remove_seq = seq(cand.id_DraftScaf).Sequence(remove_start:remove_end);
    else
        remove_seq = [];
    end
end
%%
function [insert_start insert_end] = getInsertedPosDtset(cand)
    if(cand.strand(1)==1)
        insert_start = cand.subjectIndices(2)+1;
        insert_end = cand.subjectIndices(3)-1;
    elseif(cand.strand(1)==2)
        insert_start = cand.subjectIndices(3)+1;
        insert_end = cand.subjectIndices(2)-1;
    end
end
%%
function [insert_start insert_end insert_seq] = getInsertedSeqDtset(cand,seq)
    [insert_start insert_end] = getInsertedPosDtset(cand);
    insert_seq = [];
    if insert_end >= insert_start
        if(cand.strand(1)==1)
            insert_seq = seq(cand.id_DtsetContig).Sequence(insert_start:insert_end);
        elseif(cand.strand(1)==2)
            insert_seq = seqrcomplement2(seq(cand.id_DtsetContig).Sequence(insert_start:insert_end));
        end
    end
end
%%
function [queryA_start queryA_end queryB_start queryB_end] = getRealQueryPos(cand,seq)
global edgeTrimLength;
    [n_start n_end] = getPosN(cand,seq);
    queryA_start = n_start - edgeTrimLength - cand.lengthcontigEndLengths(1) + cand.queryIndices(1) - 1;
    queryA_end = n_start - edgeTrimLength - cand.lengthcontigEndLengths(1) + cand.queryIndices(1)  + length(regexp(cand.alignmentBlastQuery{1},'[^-]')) - 2;
    queryB_start = n_end + edgeTrimLength + cand.queryIndices(3);
    queryB_end = n_end + edgeTrimLength + cand.queryIndices(3) + length(regexp(cand.alignmentBlastQuery{2},'[^-]')) - 1;
end
%%
function [name] = gapTypeName(gapType)
    if(gapType==-1)
        name = 'Negative gap';
    elseif(gapType==0)
        name = 'Zero gap';
    elseif(gapType==1)
        name = 'Positive gap';
    end
end
%%
function [filename] = hidePath(filepath)
    %[~,n,e] = fileparts(filepath);
    %filename = [n e];
    filename = filepath;
end
%%
function fastawrite2(filename, fastadata)
    len_line = 60;
    fid = fopen(filename,'w');
    for i=1:length(fastadata)
        
        fprintf(fid,'>%s\n',fastadata(i).Header);
        len_seq = length(fastadata(i).Sequence);

        for s=1:len_line:len_seq
            if(s+len_line>len_seq)
                e = len_seq;
            else
                e = s+len_line-1;
            end
            fprintf(fid,'%s\n',fastadata(i).Sequence(s:e));
        end
    end
    fclose(fid);
end
%%
function fastawritefast(filename, fastadata)
    text = '';
    for d=1:length(fastadata)
        text = [text '>' fastadata(d).Header 10 fastadata(d).Sequence 10];
    end
    fid = fopen(filename,'w');
    fprintf(fid,'%s',text);
    fclose(fid);
end
%%
function [fasta] = fastaread2(filename)
    i = 1;
    fid = fopen(filename);
    l = fgets(fid);
    while l > -1
        fasta(i).Header = l(2:end-1);

        l = fscanf(fid, '%[^>]s');
        l(find(l==char(13) | l==char(10) | l==' ')) = [];
        fasta(i).Sequence = l;

        i = i + 1;
        l = fgets(fid);
    end
    fclose(fid);
end
%%
function [bool] = isFasta(filename)
    fid = fopen(filename);
    l = fgets(fid);
    if (l(1)~='>')
        bool=0;
    else
        bool=1;
    end
    fclose(fid);
end
%%
function [draft_genome_fasta error_msg] = loadDraftFile(draftFile)
error_msg='';
    disp([10 'Reading draft: ' draftFile]);
    draft_genome_fasta = fastaUpper(fastaread2(draftFile)); 
    if(getTotalChars(draft_genome_fasta)==0)
        error_msg = ['Draft file ' draftFile ' does not have gaps'];
        return; 
    end
end
%%
function [datasets_fasta error_msg] = loadDatasetsFiles(datasetsFiles)
global stats saveDatasetsPath loadDatasetsPath;
error_msg='';
    ds = regexp(datasetsFiles,',','split');
    % Caso esteja carregando os datasets de arquivo salvo
    if ~isempty(loadDatasetsPath)
        disp(['Loading datasets: ' loadDatasetsPath]);
        if 1==1
        end
        load([loadDatasetsPath '.mat']);
        for d=1:length(datasets_fasta)
            disp([' - ' datasets_fasta(d).file ' loaded']);
        end
    else
        % Caso rodando pela primeira vez e salvando
        if ~isempty(saveDatasetsPath)
            ds = ds(1:end-1);
        end
        % Carrega arquivos
        for dt = 1:length(ds) 
            disp(['Reading dataset: ' char(ds(dt))]);
            datasets_fasta(dt).file = char(ds(dt));
            datasets_fasta(dt).fasta = fastaUpper(fastaread2(char(ds(dt)))); 
        end
        % Caso rodando pela primeira vez e salvando
        if ~isempty(saveDatasetsPath)
           disp(['Saving datasets: ' saveDatasetsPath]);
           pause(0.001);
           save([saveDatasetsPath '.mat'],'datasets_fasta','-v7.3');
        end
    end
    stats.datasetCount(length(datasets_fasta)) = 0;
end
%%
function [error_msg] = buildDatabase(datasets_fasta)
global tmp_folder blastPath saveDatasetsPath loadDatasetsPath;
error_msg='';

    if ~isempty(loadDatasetsPath)
        disp('Loading database ...');
    else
            
        % Salva em arquivo específico ou temporário
        if ~isempty(saveDatasetsPath)
            dboutputfolder = saveDatasetsPath;
        else
            dboutputfolder = [tmp_folder 'datasets'];
        end
        
        disp('Building database ...');
        % Monta fasta para banco
        cont = 0;
        for dt=1:length(datasets_fasta)
            for i=1:length(datasets_fasta(dt).fasta)
                if(~isempty(datasets_fasta(dt).fasta(i).Sequence))
                    % Para arquivo fasta
                    cont = cont + 1;
                    datasets(cont).Header = ['dataset_' num2str(dt) '_' num2str(i)];
                    datasets(cont).Sequence = datasets_fasta(dt).fasta(i).Sequence;
                end
            end 
        end
        fastawrite2([tmp_folder 'datasets.fasta'], datasets);
        makeblastdb = ['-dbtype nucl -in ' tmp_folder 'datasets.fasta -out ' dboutputfolder ' -title datasets -logfile ' tmp_folder 'makeblastdb.log'];
        [makeblastdbStatus, makeblastdbResult] = system([blastPath 'makeblastdb ' makeblastdb]);
        if makeblastdbStatus>0
            error_msg = ['Error MAKEBLASTDB: ' makeblastdbResult];
            return; 
        end
    end
end
%%
function [draft_genome_fasta error_msg] = identifyGaps(draft_genome_fasta)
% Identifica gaps (0-inicial, 1-candidato, 2-fechado, -1-não fechado/pular)
error_msg='';
    for i = 1:length(draft_genome_fasta)
        chars = getChars(draft_genome_fasta(i).Sequence); 
        if ~isempty(chars)
            lc = length(chars(:,1));
            c = 1:lc;
            chars = [chars c' zeros(lc,1)];
        end
        draft_genome_fasta(i).Pos = chars;
    end
end
%%
function writeStats(draftFile,datasets_fasta,log_round,elapsed)
global stats stats_file minScore maxEValue minIdentity contigEndLength edgeTrimLength ...
    gapChar maxRemoveLength maxInsertLength blastAlignParam blastMaxResults blastPath threads ...
    outputPrefix moreOutput version blastnAlgorithm makeblastdbAlgorithm positiveGap zeroGap negativeGap ...
    saveDatasetsPath loadDatasetsPath;

    fid = fopen(stats_file, 'a');
    fprintf(fid, '%s GENERAL STATS %s\n',repmat('-',1,20),repmat('-',1,20));
    
    fprintf(fid, '\nClosed gaps (%s): %s\n\n', gapChar, num2str(stats.totalGapsClosed) );
    fprintf(fid, 'Before FGAP: \n%s\n\n', stats.before );
    fprintf(fid, 'After FGAP: \n%s\n\n', stats.after);
    
    fprintf(fid, 'Inserted: %sbp\n', num2str(stats.insertedBases,'%.f'));
    fprintf(fid, 'Removed : %sbp\n\n', num2str(stats.removedBases,'%.f'));
    
    fprintf(fid, 'Closed gaps by each dataset:\n');
    for i=1:length(datasets_fasta)
        fprintf(fid, ' %s: %s gaps\n', hidePath(datasets_fasta(i).file), num2str(stats.datasetCount(i),'%.f'));
    end
    
    fprintf(fid, '\nClosed gaps by type:\n');
    for i=1:3
        fprintf(fid, ' %s: %s gaps\n', gapTypeName(i-2), num2str(stats.typeCount(i),'%.f'));
    end
    
    if(~isempty(log_round))
        fprintf(fid, '\n%s STATS PER ROUND %s\n\n',repmat('-',1,20),repmat('-',1,20));
        for j=1:length(log_round)
            fprintf(fid, '%s\n', log_round(j).log);
        end
    end
    
    fprintf(fid, '%s PARAMETERS %s\n\n',repmat('-',1,20),repmat('-',1,20));
    fprintf(fid, '\tdraftFile: %s\n', hidePath(draftFile));
    if ~isempty(loadDatasetsPath)
        fprintf(fid, '\tDatasets loaded: %s\n', loadDatasetsPath);
    end
    for i=1:length(datasets_fasta)
        fprintf(fid, '\tdatasetsFiles: %s\n', hidePath(datasets_fasta(i).file));    
    end
    if ~isempty(saveDatasetsPath)
        fprintf(fid, '\tDatasets saved: %s\n', saveDatasetsPath);
    end
    fprintf(fid, '\tminScore: %s\n', num2str(minScore));
    fprintf(fid, '\tmaxEValue: %s\n', num2str(maxEValue));
    fprintf(fid, '\tminIdentity: %s\n', num2str(minIdentity));
    fprintf(fid, '\tcontigEndLength: %s\n', num2str(contigEndLength));
    fprintf(fid, '\tedgeTrimLength: %s\n', num2str(edgeTrimLength));
    fprintf(fid, '\tmaxRemoveLength: %s\n', num2str(maxRemoveLength));
    fprintf(fid, '\tmaxInsertLength: %s\n', num2str(maxInsertLength));
    fprintf(fid, '\tpositiveGap: %s\n', num2str(positiveGap));
    fprintf(fid, '\tzeroGap: %s\n', num2str(zeroGap));
    fprintf(fid, '\tnegativeGap: %s\n', num2str(negativeGap));
    fprintf(fid, '\tgapChar: %s\n', gapChar);
    fprintf(fid, '\tblastPath: %s\n', num2str(blastPath));
    fprintf(fid, '\tblastAlignParam: %s\n', num2str(blastAlignParam));
    fprintf(fid, '\tblastMaxResults: %s\n', num2str(blastMaxResults));
    fprintf(fid, '\tthreads: %s\n', num2str(threads));
    fprintf(fid, '\tmoreOutput: %s\n', num2str(moreOutput));
    fprintf(fid, '\toutputPrefix: %s\n\n', hidePath(outputPrefix));

    
    fprintf(fid, '%s\n',repmat('-',1,50));
    fprintf(fid, '\n%s\n', elapsed);
    fprintf(fid, '%s\n', datestr(now));
    fprintf(fid, 'FGAP v%s\n', version);
    fprintf(fid, '%s\n', makeblastdbAlgorithm);
    fprintf(fid, '%s\n', blastnAlgorithm);
    fclose(fid);
end
%%
function [more_before more_after] = generateMoreOutput(candidates,pre_draft_genome_fasta)
    for c=1:length(candidates)
        gap_id = getGapId(candidates(c));
           
        chars = getChars(candidates(c).removedSeqDraft);
        total_chars = length(chars(:,1));
    
        more_id = [gap_id '|' num2str(candidates(c).gapType) '|' num2str(total_chars,'%.f')];
        
        [more_before(c).Sequence gap_start gap_end] = generateBeforeSequence(candidates(c),pre_draft_genome_fasta);
        more_before(c).Header = [more_id '|' num2str(gap_start,'%.f') '|' num2str(gap_end,'%.f')];
        
        [more_after(c).Sequence insert_start insert_end] = generateAfterSequence(candidates(c),pre_draft_genome_fasta); 
        more_after(c).Header = [more_id '|' num2str(insert_start,'%.f') '|' num2str(insert_end,'%.f')];
    end
end
%%
function writeMoreOutput(more_before,out_more_before,more_after,out_more_after)
    more_before_all = [];
    more_after_all = [];
    for i=1:length(more_before)
        more_before_all = [more_before_all more_before{i}];
        more_after_all = [more_after_all more_after{i}];
    end
    if ~isempty(out_more_before)
        fastawrite2(out_more_before,more_before_all);
        fastawrite2(out_more_after,more_after_all);
    end
end
%%
function writeLog(out_log,log)
    fid = fopen(out_log, 'wt');
    for i=1:length(log)
        fprintf(fid, '%s\n', log{i}{:});
    end
    fclose(fid);
end
%%
function writeFasta(out_fasta, draft_genome_fasta)
    fastawrite2(out_fasta, draft_genome_fasta);
end
%%
function [gap_id] = getGapId(cand)

    gap_id = [num2str(cand.id_DraftScaf) '_' num2str(cand.id_DraftChar)];
end
%%
function [after_seq insert_start insert_end] = generateAfterSequence(cand,pre_draft_genome_fasta)
global gapChar plusBasesMoreOutput;
    

    len_seq = length(pre_draft_genome_fasta(cand.id_DraftScaf).Sequence);
    if(cand.gapType==-1)
        insert_start = 0;
        insert_end = 0;
        
        st_plus = cand.removedPosDraft(1) - plusBasesMoreOutput;
        if st_plus < 1 st_plus=1; end
        en_plus = cand.removedPosDraft(2) + plusBasesMoreOutput;
        if en_plus > len_seq en_plus = len_seq; end
        after_seq = [pre_draft_genome_fasta(cand.id_DraftScaf).Sequence(st_plus:cand.removedPosDraft(1)-1) ...
                    pre_draft_genome_fasta(cand.id_DraftScaf).Sequence(cand.removedPosDraft(2)+1:en_plus)];
        %Condensa gaps das pontas para melhorar validação
        after_seq = regexprep(after_seq,[gapChar '[' gapChar ']*'],gapChar);    
        
    else
        len_seq = length(pre_draft_genome_fasta(cand.id_DraftScaf).Sequence);

        st_plus = cand.removedPosDraft(1) - plusBasesMoreOutput;
        if st_plus < 1 st_plus=1; end
        en_plus = cand.removedPosDraft(2) + plusBasesMoreOutput;
        if en_plus > len_seq en_plus = len_seq; end

        pre_seq = pre_draft_genome_fasta(cand.id_DraftScaf).Sequence(st_plus:cand.removedPosDraft(1)-1);
        pos_seq = pre_draft_genome_fasta(cand.id_DraftScaf).Sequence(cand.removedPosDraft(2)+1:en_plus);
        
        % Retira sequencia caso possua N nos contigsEnds para fazer validação
        [~,n_pre_en] = regexpi(pre_seq,[gapChar '[' gapChar ']*']); 
        if ~isempty(n_pre_en)
            pre_seq = pre_seq(n_pre_en(length(n_pre_en))+1:end);
        end
        [n_pos_st,~] = regexpi(pos_seq,[gapChar '[' gapChar ']*']); 
        if ~isempty(n_pos_st)
            pos_seq = pos_seq(1:n_pos_st(1)-1);
        end
        after_seq = [pre_seq lower(cand.insertedSeqDtset) pos_seq ];  
        insert_start = length(pre_seq)+1;
        insert_end = insert_start + length(cand.insertedSeqDtset) - 1;
        
        %Condensa gaps das pontas para melhorar validação
        %after_seq = regexprep(after_seq,[gapChar '[' gapChar ']*'],gapChar); 
    end 
end
%%
function [before_seq gap_start gap_end] = generateBeforeSequence(cand,pre_draft_genome_fasta)
global edgeTrimLength plusBasesMoreOutput;

    len_seq = length(pre_draft_genome_fasta(cand.id_DraftScaf).Sequence);
    [n_start n_end] = getPosN(cand,pre_draft_genome_fasta);

    st_plus = n_start - edgeTrimLength - plusBasesMoreOutput;
    if st_plus < 1 st_plus=1; end
    en_plus = n_end + edgeTrimLength + plusBasesMoreOutput;
    if en_plus > len_seq en_plus = len_seq; end
    
    before_seq = pre_draft_genome_fasta(cand.id_DraftScaf).Sequence(st_plus:en_plus); 
    gap_start = edgeTrimLength + plusBasesMoreOutput + 1;
    gap_end = gap_start + (n_end - n_start);
    
end
%%
function [seqrc] = seqrcomplement2(seq)
    %Inverte    
    seqr = seq(end:-1:1);
    seqrc = ones(1,length(seqr));
    l = {'A','C','G','T','R','Y','K','M','S','W','B','D','H','V','N','-','*'};
    t = {'T','G','C','A','Y','R','M','K','S','W','V','H','D','B','N','-','*'};
    for i=1:length(l)
        seqrc(seqr==l{i})=t{i};
    end
    seqrc = char(seqrc);
end
%%
function [] = showHelp()  
    disp([ 10 'Usage in command-line mode (compiled): ./run_fgap.sh <MCR installation folder> -d <draft file> -a "<dataset(s) file(s)>" [parameters]']);
    disp(['Usage in Matlab/Octave (source): fgap -d <draft file> -a ''<dataset(s) file(s)>'' [parameters]' 10]);
    
    disp(['-d /--draft-file' 9 'Draft genome file [fasta format - Ex: ''draft.fasta'']']);
    disp(['-a /--datasets-files' 9 'List of datasets files to close gaps [fasta format - Ex: ''dataset1.fasta,dataset2.fasta'']' 10]);
    %disp(['--save-datasets']);
    %disp(['--load-datasets']);
    
    disp(['-s /--min-score' 9 9 'Min Score (raw) to return results from BLAST (integer) - Default: 25']);
    disp(['-e /--max-evalue' 9 'Max E-Value to return results from BLAST (float) - Default: 1e-7']);
    disp(['-i /--min-identity' 9 'Min identity (%) to return results from BLAST (integer [0-100]) - Default: 70' 10]);

    disp(['-C /--contig-end-length' 9 'Length (bp) of contig ends to perform BLAST alignment (integer) - Default: 300']);
    disp(['-T /--edge-trim-length' 9 'Length of ignored bases (bp) upstream and downstrem of the gap (integer) - Default: 0']);
    disp(['-R /--max-remove-length' 9 'Max number of bases (bp) that can be removed (integer) - Default: 500']);
    disp(['-I /--max-insert-length' 9 'Max number of bases (bp) that can be inserted (integer) - Default: 500' 10]);
    
    disp(['-p /--positive-gap' 9 'Enable closing of positive gaps (with insertion) (integer [0-1]) - Default: 1']);
    disp(['-z /--zero-gap' 9 9 'Enable closing of zero gaps (without insert any base) (integer [0-1]) - Default: 0']);
    disp(['-g /--negative-gap' 9 'Enable closing of negative gaps (overlapping contig ends) (integer [0-1]) - Default: 0' 10]);
    
    disp(['-c /--gap-char' 9 9 9 9 'Base that represents the gap (char) - Default: ''N''']);
    disp(['-b /--blast-path' 9 9 9 'Blast+ package path (only makeblastdb and blastn are needed, version 2.2.28+ or higher) - Default: ''''']);
    disp(['-l /--blast-alignment-parameters' 9 'BLAST alignment parameters (opengap,extendgap,match,mismatch,wordsize) - Default: ''1,1,1,-3,15''']);
    disp(['-r /--blast-max-results' 9 9 9 'Max results from BLAST for each query (integer) - Default: 200']);
    disp(['-t /--threads' 9 9 9 9 'Number of threads (integer) - Default: 1' 10]);
    
    disp(['-m /--more-output' 9 'More output files with gap regions after and before gap closing (integer [0-1]) - Default: 0']);
    disp(['-o /--output-prefix' 9 'Output prefix [File or folder - Ex: ''out'' or ''out_folder/out'' ] - Default: ''output_fgap''']);
    disp(['-h /--help' 9 9 'This help message']);
end
