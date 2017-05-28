%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Copyright (c) 2017, TCALAB                      %
%           Department of Computer Engineering & Informatics           %
%                         University of Patras                         %
%                                                                      %
% Authors: Filippos Filippou <filippo@ceid.upatras.gr>                 %
%          Michail Mavropoulos <mavropoulo@ceid.upatras.gr>            %
%          Georgios Keramidas <gkeramidas@ceid.upatras.gr>             %
%                                                                      %
% Description:                                                         %
% This script produces fault maps for a BTB structure based on a range %
% of failure probabilities. The produced fault map defines if the used %
% BTB entry is faulty or not. This fault map can be used by a modified %
% architecture emulator like gem5.                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = fault_map_btb_permutation(faulty_maps)
    clc;
    global statistics;
    global debugging;
    global export_fmaps;
    statistics = 1; % 1/0 -> Yes/No
    debugging = 0; % 1/0 -> Yes/No
    export_fmaps = 0; % 1/0 -> Yes/No
    
    %------------------------ BTB Configurations -------------------------%
    % BTB entry bits = tag(2 bytes - 32bit || 4 bytes - 64bit) 
    %                + target (4 bytes - 32bit || 8 bytes - 64bit) 
    %                + threadid (2 bytes)
    blocksize = 112; % bits ( 8 bytes = 64 bits, 14 bytes = 112 bits )
    btb_entries = [512,1024,2048,4096];
    associativity = [2,4];
    
    %{
     ---------------------------------pfails-----------------------------------
      0: 1.0*10^-5,  1: 5.0*10^-5,  2: 1.0*10^-4,  3: 2.0*10^-4,  4: 5.0*10^-4,
      5: 8.0*10^-4,  6: 1.0*10^-3,  7: 1.5*10^-3,  8: 2.0*10^-3,  9: 2.5*10^-3,
     10: 3.0*10^-3, 11: 3.5*10^-3, 12: 4.0*10^-3, 13: 4.5*10^-3, 14: 5.0*10^-3
     --------------------------------------------------------------------------
    %}
    pfails = [10^-5, 5*10^-5, 10^-4, 2*10^-4, 5*10^-4, 8*10^-4, 10^-3, 1.5*10^-3, 2*10^-3, 2.5*10^-3, 3*10^-3, 3.5*10^-3, 4*10^-3, 4.5*10^-3, 5*10^-3];
    per = 0;
    
    btb_types{1} = 'faulty';
    btb_types{2} = '1CR';
    btb_types{3} = '2CRs';
    btb_types{4} = '4CRs';
    btb_types{5} = '8CRs';
    btb_types{6} = '16CRs';
    btb_types{7} = '32CRs';
    btb_types{8} = '64CRs';
    btb_types{9} = '128CRs';
    
    group_permutation = 0; % 1/0 -> Yes/No
    read_existing_fmap = 0; % 1/0 -> Yes/No
    %---------------------------------------------------------------------%
    
    input_parent_dir = strcat('');
    if read_existing_fmap == 1
        input_parent_dir = strcat('subitmaps_folder/');
        if getDebug
            disp(strcat('Input Fault-Map Parent Dir "',input_parent_dir,'"'))
        end
        if exist(input_parent_dir, 'dir') == 0
            disp(strcat('Input Fault-Map Parent Dir "',input_parent_dir,'"'))
            error('ERROR - This fault-map parent directory does not exists!')
        end
    end
    
    output_parent_dir = strcat('');
    if getFmaps
        format shortg
        file_date = clock;
        output_parent_dir = strcat('subitmaps_',num2str(file_date(1)),'y_',num2str(file_date(2)),'m_',num2str(file_date(3)),'d/');
        disp(strcat('Output Fault-Map Parent Dir "',output_parent_dir,'"'))
        if exist(output_parent_dir, 'dir') == 0
            [status,message,messageid] = mkdir(strcat(output_parent_dir));
            if status
                disp(strcat(message))
            end
            % Header File
            max_btb_i = max(btb_entries)/min(associativity);
            max_btb_j = max(associativity);
            [fname,ErrMsg] = fopen(strcat(output_parent_dir,'create-subitmap.hh'),'w');
            error(ErrMsg);
            fprintf(fname,'#include <stdlib.h>\n\n');
            fprintf(fname,'#ifndef CREATESUBITMAP_H\n');
            fprintf(fname,'#define CREATESUBITMAP_H\n\n');
            fprintf(fname,'void updatesublkbitmap(unsigned int btb_entries, unsigned int assoc);\n\n');
            fprintf(fname,'extern int sbitmap[%d][%d];\n\n',max_btb_i,max_btb_j);
            fprintf(fname,'#endif /* CREATESUBITMAP_H */');
            fclose(fname);
        else
            disp(strcat('Directory "',output_parent_dir,'" already exists.'))
        end
    end
    
    if getStats
        [statsfile,ErrMsg] = fopen('subitmap_stats.txt','w');
        error(ErrMsg);
        fprintf(statsfile,'CRs, FFSETS Before, FFSETS After, Common FFSETS\n');
        fclose(statsfile);
    end
    
    for output_subdir=1:faulty_maps
        
        output_fmap_dir = strcat(output_parent_dir,'subitmaps_',num2str(output_subdir),'/');
        
        if exist(output_fmap_dir, 'dir') == 0
            input_fmap_dir = strcat('');
            if read_existing_fmap == 1
                input_subdir = output_subdir;
                input_fmap_dir = strcat(input_parent_dir,'subitmaps_',num2str(input_subdir),'/');
                if getDebug
                    disp(strcat('Input Fault-Map Dir "',input_fmap_dir,'"'))
                end
                if exist(input_fmap_dir, 'dir') == 0
                    disp(strcat('Input Fault-Map Dir "',input_fmap_dir,'"'))
                    error('ERROR - This fault-map directory does not exists!')
                end
            end
            
            if getDebug
                disp(strcat('Output Fault-Map Dir "',output_fmap_dir,'"'))
            end
            
            sbitmaps(blocksize,btb_entries,associativity,pfails,btb_types,group_permutation,per,read_existing_fmap,input_fmap_dir,output_fmap_dir);
            
            disp(strcat('Fault Map [',num2str(output_subdir),'] @ :/',output_fmap_dir))
            
            if getStats
                [statsfile,ErrMsg] = fopen('subitmap_stats.txt','a');
                error(ErrMsg);
                fprintf(statsfile,'\n');
                fclose(statsfile);
            end
        end %exist(output_fmap_dir)
        
    end %output_subdir
end %end of function fault_map_btb_permutation()

function [debugging_display] = getDebug
    global debugging;
    debugging_display = debugging;
end

function [fmaps_out] = getFmaps
    global export_fmaps;
    fmaps_out = export_fmaps;
end

function [stats_display] = getStats
    global statistics;
    stats_display = statistics;
end

function [] = sbitmaps(blocksize,btb_entries,associativity,pfails,btb_types,group_permutation,per,read_existing_fmap,input_fmap_dir,output_fmap_dir)
    % print possitions y,x
    % faults inside the block
    % bitmaps
    
    if getFmaps
        [status,message,messageid] = mkdir(strcat(output_fmap_dir));
        if status
            disp(strcat(message))
        end
        
        max_btb_i = max(btb_entries)/min(associativity);
        max_btb_j = max(associativity);
        
        for i=1:numel(pfails)
            for j=1:numel(btb_types)
                [fname,ErrMsg] = fopen(strcat(output_fmap_dir,'create-subitmap',num2str(i-1),'_',num2str(j),'.cc'),'w');
                error(ErrMsg);
                if j == 1
                    fprintf(fname,'//Subitmap of Faulty BTB (F-BTB) - pfails=%f\n',pfails(i));
                elseif j == 2
                    fprintf(fname,'//Subitmap of Permutation BTB using 1 CR/way (1CR P-BTB) - pfail = %f\n',pfails(i));
                else
                    CRs = 2^(j-1);
                    fprintf(fname,'//Subitmap of Permutation BTB using %d CRs/way (%dCRs P-BTB) - pfail = %f\n',CRs,CRs,pfails(i));
                end
                fprintf(fname,'#include <stdio.h>\n');
                fprintf(fname,'#include "create-subitmap.hh"\n\n');
                fprintf(fname,'int sbitmap[%d][%d]={0};\n\n',max_btb_i,max_btb_j);
                fprintf(fname,'void updatesublkbitmap(unsigned int btb_entries, unsigned int assoc)\n');
                fprintf(fname,'{  //btb_entries-assoc\n');
                fclose(fname);
            end
        end
    end
    
    for size=1:numel(btb_entries)
        if getFmaps
            for i=1:numel(pfails)
                for j=1:numel(btb_types)
                    [fname,ErrMsg] = fopen(strcat(output_fmap_dir,'create-subitmap',num2str(i-1),'_',num2str(j),'.cc'),'a');
                    error(ErrMsg);
                    if size==1
                        fprintf(fname,'  if (btb_entries == %d){\n',btb_entries(size));
                    else
                        fprintf(fname,'  else if (btb_entries == %d){\n',btb_entries(size)); 
                    end
                    fclose(fname);
                end
            end
        end
        
        for assoc=1:numel(associativity)
            %BTB shape
            sets = btb_entries(size) / associativity(assoc);
            set_bits = blocksize * associativity(assoc); % bits per set
            
            nfaults = zeros(numel(pfails),1);
            Xd = [];
            Yd = [];
            
            if getFmaps
                for i=1:numel(pfails)
                    for j=1:numel(btb_types)
                        [fname,ErrMsg] = fopen(strcat(output_fmap_dir,'create-subitmap',num2str(i-1),'_',num2str(j),'.cc'),'a');
                        error(ErrMsg);
                        if assoc==1
                            fprintf(fname,'    if (assoc == %d){  //size:%d\n',associativity(assoc),btb_entries(size));
                        else
                            fprintf(fname,'    else if (assoc == %d){  //size:%d\n',associativity(assoc),btb_entries(size));
                        end
                        fclose(fname);
                    end
                end
            end
            
            for p_number=1:numel(pfails)
                disp(strcat('BTB Entries = [',num2str(btb_entries(size)),...
                         '], Assoc. = [',num2str(associativity(assoc)),...
                         '], Pfail[',num2str(p_number-1),'] = [',num2str(pfails(p_number),'%.2e'),']'))
                
                if read_existing_fmap == 0
                    [btb,nfaults,Xd,Yd] = insert_faults(sets,btb_entries(size),blocksize,set_bits,per,pfails(p_number),nfaults,Xd,Yd,p_number);
                else % read_existing_fmap == 1
                    input_fmap_file = strcat(input_fmap_dir,'create-subitmap',num2str(p_number-1),'_1.cc');
                    if getDebug
                        disp(strcat('Input Fault-Map File "',input_fmap_file,'"'))
                    end
                    if exist(input_fmap_file, 'file') == 0
                        disp(strcat('Input Fault-Map File "',input_fmap_file,'"'))
                        error('ERROR - This fault-map file does not exists!')
                    end
                    [btb,faulty_bits] = read_fmap_file(btb_entries(size), associativity(assoc), blocksize, input_fmap_file);
                    nfaults(p_number) = faulty_bits;
                end
                
                run_all(btb,btb_entries(size),sets,associativity(assoc),blocksize,btb_types,p_number,set_bits,nfaults,group_permutation,output_fmap_dir);
            end %p_number
            
        end %assoc
        
        if getFmaps
            for i=1:numel(pfails)
                for j=1:numel(btb_types)
                    [fname,ErrMsg] = fopen(strcat(output_fmap_dir,'create-subitmap',num2str(i-1),'_',num2str(j),'.cc'),'a');
                    error(ErrMsg);
                    fprintf(fname,'    else { //incorrect associativity\n');
                    fprintf(fname,'      printf("ERROR - BTB Fault Map: Incorrect associativity\\n");\n');
                    fprintf(fname,'      exit(1);\n');
                    fprintf(fname,'    }\n');
                    fprintf(fname,'  }\n'); % btb_entries end
                    fclose(fname);
                end
            end
        end
        
    end %size
    
    if getFmaps
        for i=1:numel(pfails)
            for j=1:numel(btb_types)
                [fname,ErrMsg] = fopen(strcat(output_fmap_dir,'create-subitmap',num2str(i-1),'_',num2str(j),'.cc'),'a');
                error(ErrMsg);
                fprintf(fname,'  else {  //incorrect btb_entries\n');
                fprintf(fname,'    printf("ERROR - BTB Fault Map: Incorrect number of BTB entries\\n");\n');
                fprintf(fname,'    exit(1);\n');
                fprintf(fname,'  }\n');
                fprintf(fname,'}\n'); % updatesublkbitmap() end
                fclose(fname);
            end
        end
    end
    
end %end of function sbitmaps()

function [btb,nfaults,Xd,Yd] = insert_faults(sets,btb_entries,blocksize,set_bits,per,pfail,nfaults,Xd,Yd,p_number)
    
    if getDebug
        disp(strcat('Insert Faults: ',num2str(p_number-1)))
    end
    
    numfaults=round((btb_entries*blocksize)*pfail); %faults total number
    dev=(numfaults*per)/100;
    nfaults(p_number)=round(numfaults+dev);
    
    if p_number>1
        Xd(p_number,:)=Xd(p_number-1,1:nfaults(p_number-1));
        Yd(p_number,:)=Yd(p_number-1,1:nfaults(p_number-1));
    else
        Xd(1,:)=zeros(1,nfaults(1));
        Yd(1,:)=zeros(1,nfaults(1));
    end
    
    %insert data faults
    k=1;
    flag=0;
    while (k==1)
        if (flag==0)
            rng('shuffle');
            if p_number>1
                x_temp=randi(set_bits,nfaults(p_number)-nfaults(p_number-1),1);
                Xd(p_number,nfaults(p_number-1)+1:nfaults(p_number)) = x_temp;
            else
                Xd(p_number,1:nfaults(p_number))=randi(set_bits,nfaults(p_number),1);
            end
            rng('shuffle');
            if p_number>1
                y_temp=randi(sets,nfaults(p_number)-nfaults(p_number-1),1);
                Yd(p_number,nfaults(p_number-1)+1:nfaults(p_number)) = y_temp;
            else
                Yd(p_number,1:nfaults(p_number))=randi(sets,nfaults(p_number),1);
            end
            %x = unifrnd(1, set_bits, nfaults, 1);
            %y = unifrnd(1, sets, nfaults, 1);
            %x1d=round(x);
            %y1d=round(y);
            flag=1;
        else
            %x1=x1(randperm(nfaults));
            %y = unifrnd(1, sets, nfaults, 1);
            % y1=ceil(y);
            %y1=y1(randperm(nfaults));
        end
        
        k=0;
        for p=1:nfaults(p_number)-1
            xk=Xd(p_number,p);
            yk=Yd(p_number,p);
            for q=p+1:nfaults(p_number)
                if (xk==Xd(p_number,q) && yk==Yd(p_number,q))
                    k=1;
                    Xd(p_number,q)=mod(Xd(p_number,q),set_bits) + 1;
                end
            end
        end
    end
    
    btb=zeros(sets,set_bits);
    for k=1:nfaults(p_number)
        btb(Yd(p_number,k),Xd(p_number,k))=1;
    end
    
end %end of function insert_faults()

function [output_btb,faulty_bits] = read_fmap_file(btb_entries, associativity, blocksize, input_fmap_file)
    
    % count how many rows must be skipped from old fault map file
    head_lines = 8;
    comment_lines = 5;
    error_lines = 4;
    
    file_btb_entries = [512,1024,2048,4096];
    file_btb_assoc = [2,4];
    
    entries_index = find(file_btb_entries == btb_entries);
    assoc_index = find(file_btb_assoc == associativity);
    
    if numel(entries_index) == 0
        error('read_fmap_file() - ERROR: Wrong entries!')
    elseif numel(assoc_index) == 0
        error('read_fmap_file() - ERROR: Wrong associativity!')
    end
    
    skip_entries = 0;
    % remove previous btb configs
    for i = 1:(entries_index - 1) % remove btb entries until entries_index
        entries = file_btb_entries(i);
        skip_assoc = 0;
        for j = 1:numel(file_btb_assoc) % remove all assocs
            assoc = file_btb_assoc(j);
            sets = entries/assoc;
            skip_assoc = skip_assoc + (sets * assoc) + comment_lines + 2;
        end
        skip_assoc = skip_assoc + error_lines;
        skip_entries = skip_entries + skip_assoc + 2;
    end
    
    skip_entries_lines = skip_entries + 1;
    
    entries = file_btb_entries(entries_index);
    skip_assoc = 0;
    for j = 1:(assoc_index - 1) % remove assocs until assoc_index
        assoc = file_btb_assoc(j);
        sets = entries/assoc;
        skip_assoc = skip_assoc + (sets * assoc) + comment_lines + 2;
    end
    
    skip_assoc_lines = skip_assoc + 1;
    
    skipped_lines = head_lines + skip_entries_lines + skip_assoc_lines;
    
    if getDebug
        disp(strcat('Skepped Head Lines:  [',num2str(head_lines),']'))
        disp(strcat('Skepped Entry Lines: [',num2str(skip_entries_lines),']'))
        disp(strcat('Skepped Assoc Lines: [',num2str(skip_assoc_lines),']'))
        disp(strcat('Total Skepped Lines: [',num2str(skipped_lines),']'))
    end
    
    fileID = fopen(strcat(input_fmap_file),'r');
    textscan(fileID,'%s',skipped_lines,'delimiter','\n'); % skip lines
    InputText = textscan(fileID,'sbitmap[%f][%f] = %f;'); % get input fault map
    fclose(fileID);
    
    input_fmap = InputText{3};

    % create cache fault map
    sets = btb_entries/associativity;
    output_btb=zeros(sets,associativity * blocksize);
    faulty_bits = 0;
    for set=1:sets
        for way=1:associativity
            for sbit=1:blocksize
                % mark all bits of a faulty block as faulty =>
                % faulty_bit = faulty_blocks * block_size
                if (input_fmap((set-1)*associativity + way) == 1)
                    output_btb(set,(way-1)*blocksize + sbit) = 1;
                    faulty_bits = faulty_bits + 1;
                else
                    output_btb(set,(way-1)*blocksize + sbit) = 0;
                end
            end
        end
    end
    
end %end of function read_fmap_file()

function [] = run_all(btb,btb_entries,sets,assoc,blk_size,btb_types,p_number,set_bits,nfaults,group_permutation,output_fmap_dir)
    
    %reordered_flag = 0;
    %[reordered_btb] = btb_reorder(btb,sets,assoc,blk_size,set_bits);
    %reordered_flag = 1;
    
    if getStats
        [statsfile,ErrMsg] = fopen('subitmap_stats.txt','a');
        error(ErrMsg);
    end
    
    %------Faulty BTB (F-BTB)-----%
    if getFmaps
        if getDebug
            disp('BTB Type: F-BTB')
        end
        print_sbitmap(btb,btb_entries,sets,assoc,blk_size,p_number,set_bits,nfaults,output_fmap_dir,1,0);
    end
    
    for btb_t=2:numel(btb_types)
        %------Permutation BTB k-CR(s)/way (k-CR(s) P-BTB)-----%
        CRsPerWay = 2^(btb_t-2);
        subsets = sets/CRsPerWay;
        
        if subsets > 1
            if group_permutation == 0
                if getDebug
                    disp(strcat('BTB Type:',num2str(CRsPerWay),'-CR(s) P-BTB | Step Permutation'))
                end
                [CR,ffsets_before,ffsets_after,common_ffsets] = step_way_permutation(btb,sets,assoc,blk_size,CRsPerWay);
            else
                if getDebug
                    disp(strcat('BTB Type:',num2str(CRsPerWay),'-CR(s) P-BTB | Group Permutation'))
                end
                [CR,ffsets_before,ffsets_after,common_ffsets] = group_way_permutation(btb,sets,assoc,blk_size,CRsPerWay);
            end
            
            if getStats
                fprintf(statsfile,'%dCR, %d, %d, %d, ',CRsPerWay,ffsets_before,ffsets_after,common_ffsets);
                disp('###############################################')
                disp(strcat(num2str(CRsPerWay),'-CR(s) P-BTB Statistics - Step Permutation'))
                disp('###############################################')
                disp(strcat('Sets:',num2str(sets),', Subsets:',num2str(subsets),', Assoc.:',num2str(assoc),', Entries:',num2str(sets*assoc)))
                disp(strcat('Fully Faulty Sets Before: [',num2str(ffsets_before),']'))
                disp(strcat('Fully Faulty Sets After:  [',num2str(ffsets_after),']'))
                disp(strcat('Common Fully Faulty Sets: [',num2str(common_ffsets),']'))
                disp('###############################################')
            end
            
            if getFmaps
                % create new BTB using permutaion array CR
                black_btb=zeros(sets,set_bits);
                subsets = sets/CRsPerWay;
                for subarray=1:CRsPerWay
                    subsets_start = (subarray-1)*subsets;
                    for set=1:subsets
                        for way=1:assoc
                            set_offset = bitxor((set-1),CR(subarray,way))+1;
                            for sbit=1:blk_size
                                if (btb(subsets_start+set_offset,(way-1)*blk_size+sbit) == 1)
                                    black_btb(subsets_start+set,(way-1)*blk_size+sbit) = 1;
                                else
                                    black_btb(subsets_start+set,(way-1)*blk_size+sbit) = 0;
                                end
                            end
                        end
                    end
                end
                print_sbitmap(black_btb,btb_entries,sets,assoc,blk_size,p_number,set_bits,nfaults,output_fmap_dir,btb_t,common_ffsets);
            end
        else %subsets == 1
            if getFmaps
                [fname,ErrMsg] = fopen(strcat(output_fmap_dir,'create-subitmap',num2str(p_number-1),'_',num2str(btb_t),'.cc'),'a');
                error(ErrMsg);
                fprintf(fname,'      printf("ERROR - BTB Fault Map: Incorrect number of CRs\\n");\n');
                fprintf(fname,'      exit(1);\n');
                fprintf(fname,'    }\n'); % assoc end
                fclose(fname);
            end
        end %subsets > 1
    end %btb_t
    
    if getStats
        fclose(statsfile);
    end
end %end of function run_all()

function [reordered_btb] = btb_reorder(btb,sets,assoc,blk_size,set_bits)
    
    reordered_btb = zeros(sets,set_bits);
    % faulty blocks per way
    for way=1:assoc
        switch way
            case 1, new_way = 1;
            case 2, new_way = 3;
            case 3, new_way = 4;
            case 4, new_way = 2;
        end
        for set=1:sets
            for sbit=1:blk_size
                bit_pos = (way-1)*blk_size+sbit;
                bit_pos_new = (new_way-1)*blk_size+sbit;
                reordered_btb(set,bit_pos_new) = btb(set,bit_pos);
            end
        end
    end

end %end of function btb_reorder()

function [CR,ffsets_num_before,ffsets_num_after,common_ffsets] = step_way_permutation(btb,sets,assoc,blk_size,CRsPerWay)
    
    CR=zeros(CRsPerWay,assoc);
    subsets = sets/CRsPerWay;
    S = 0:(subsets-1);
    common_ffsets = 0;
    ffsets_num_before = 0;
    ffsets_num_after = 0;
    
    for subarray=1:CRsPerWay
        subsets_start = (subarray-1)*subsets;
        
        % FBw(i,j): Faulty block i @ way j
        % FBw_num(1,j): Number of faulty blocks @ way j
        % FSb(i,j): Faulty set i with j faulty blocks
        % FSb_num(1,j): Number of faulty sets with j faulty blocks
        FBw = zeros(subsets,assoc);
        FBw_num = zeros(1,assoc); % max_val = subsets
        FSb = zeros(subsets,assoc);
        FSb_num = zeros(1,assoc); % max_val = subsets
        
        % find faulty physical addresses per way
        for way=1:assoc
            for set=1:subsets
                errors = 0;
                for sbit=1:blk_size
                    if (btb(subsets_start+set,(way-1)*blk_size+sbit) == 1)
                        errors = 1;
                    end
                end
                if errors == 1
                    FBw_num(1,way) = FBw_num(1,way) + 1;
                    FBw_pnt = FBw_num(1,way);
                    FBw(FBw_pnt,way) = set - 1;
                end
            end
            if getDebug
                disp(strcat('FBw',num2str(way),' = [',num2str(FBw(:,way)'),']'))
            end
        end
        
        % find the Fully Faulty Sets initially
        FFS = zeros(subsets,1);
        FFS_num = 0;
        for set_addr=0:subsets-1
            FB_count = 0;
            for way=1:assoc
                FBw_end = FBw_num(1,way);
                if any(set_addr == FBw(1:FBw_end,way))
                    % if value exists in FBw(way)
                    FB_count = FB_count + 1;
                end
            end
            if FB_count == assoc
                FFS_num = FFS_num + 1;
                FFS(FFS_num,1) = set_addr; % insert #set to FFS_init
            end
        end
        
        if FFS_num > 0
            FFS_before = FFS(1:FFS_num,1);
        else
            FFS_before = [];
        end
        
        % initially FSb1 == FBw1
        FSb(:,1) = FBw(:,1);
        FSb_num(1,1) = FBw_num(1,1);
        
        % Permutation Algorithm
        for way=2:assoc
            if getDebug
                disp(strcat('Calculation of CR(',num2str(subarray),',',num2str(way),')'))
            end
            K = FBw_num(1,way);
            if K > 0
                CCR = S; % Candidates for CR = S
                for set_fblks=way-1:-1:1 % sets with "#set_fblks" faulty blocks until current way
                    if getDebug
                        disp(strcat('XFB = FSb',num2str(set_fblks),' XOR FBw',num2str(way)))
                    end
                    % FSb(way-1)-> ... -> FSb(1)
                    if numel(CCR) > 1
                        J = FSb_num(1,set_fblks);
                        if J > 0 %Calculate CR
                            if getDebug
                                disp(strcat('FSb',num2str(set_fblks),' = [',num2str(FSb(1:J,set_fblks)'),']'))
                                disp(strcat('FBw',num2str(way),' = [',num2str(FBw(1:K,way)'),']'))
                            end
                            CNT_XFB = zeros(subsets,1); % Repeats of each possible result
                            for j=1:J % FSb(j,set_fblks)
                                for k=1:K % FBw(k,way)
                                    XFB = bitxor(FSb(j,set_fblks),FBw(k,way)); % XFB = FSb(way-1) XOR FBw(way)
                                    CNT_XFB(XFB+1,1) = CNT_XFB(XFB+1,1) + 1; % Repeats++
                                end
                            end
                            
                            % from previous CCR values choose the one that
                            % respect CCR(k) column in CNT_XFB has min value
                            
                            % EXAMPLE:
                            % CNT_XFB = {D0,D1,D2,D3,D4,D5,D6,D7} = {8,1,2,7,3,2,0,1}
                            % CCR = {0,2,4,5}
                            % CNT_CCR = {D0,D2,D4,D5} = {8,2,3,2} -> CNT_MIN = 2, CCR_NUM = 2
                            % CCR_new = {2,5}
                            
                            CNT_CCR = zeros(numel(CCR),1);
                            for cr_num=1:numel(CCR)
                                CNT_CCR(cr_num,1) = CNT_XFB(CCR(cr_num)+1,1);
                            end
                            CNT_MIN = min(CNT_CCR); % Least repeats
                            CCR_NUM = histc(CNT_CCR,CNT_MIN); % How many values have the same number of least repeats
                            
                            CCR_new = zeros(1,CCR_NUM); % Candidates for CR
                            cnt_pnt = 1;
                            for cr_num=1:numel(CCR)
                                if (CNT_CCR(cr_num,1) == CNT_MIN)
                                    CCR_new(cnt_pnt) = CCR(cr_num); % Least repeated values
                                    cnt_pnt = cnt_pnt + 1;
                                end
                            end
                            if (CCR_NUM == 0) || (numel(CCR_new) ~= CCR_NUM)
                                ErrMsg = 'Something went wrong with CR calculation';
                                error(ErrMsg)
                            end
                            
                            if getDebug
                                disp(strcat('New Set Faulty Blocks = [',num2str(set_fblks),']'))
                                disp(strcat('CCR_input = [',num2str(CCR(:)'),'] -> [',num2str(numel(CCR)),'] elements'))
                                disp(strcat('CNT_XFB = [',num2str(CNT_XFB(:)'),'] -> [',num2str(numel(CNT_XFB)),'] elements'))
                                disp(strcat('CNT_CCR = [',num2str(CNT_CCR(:)'),'] -> [',num2str(numel(CNT_CCR)),'] elements'))
                                disp(strcat('CNT_MIN = ',num2str(CNT_MIN),'], CCR_NUM = [',num2str(CCR_NUM),']'))
                                disp(strcat('CCR_output = [',num2str(CCR_new(:)'),'] -> [',num2str(numel(CCR_new)),'] elements'))
                            end
                            
                            CCR = CCR_new;
                        else
                            if getDebug
                                disp(strcat('FSb',num2str(set_fblks),' = []'))
                            end
                        end %J
                    end %numel(CCR)
                end %set_fblks
                
                CR(subarray,way) = CCR(1,1); % choose the first one
                %CR(subarray,way) = randsample(CCR,1); % choose randomly
                
                % FBw(way) = FBw(way) XOR CR(way)
                FBw(:,way) = bitxor(FBw(:,way),CR(subarray,way));
                
                % FSb sets after current step
                FSb = zeros(subsets,assoc);
                FSb_num = zeros(1,assoc);
                for set_addr=0:subsets-1
                    FSb_count = 0;
                    for attached_way=1:way
                        FBw_end = FBw_num(1,attached_way);
                        if any(set_addr == FBw(1:FBw_end,attached_way))
                            % if value exists in FBw(attached_way)
                            FSb_count = FSb_count + 1;
                        end
                    end
                    if FSb_count > 0
                        FSb_num(1,FSb_count) = FSb_num(1,FSb_count) + 1;
                        FSb_pnt = FSb_num(1,FSb_count);
                        FSb(FSb_pnt,FSb_count) = set_addr; % insert #set to FSb
                    end
                end %set
            else
                if getDebug
                    disp(strcat('FBw',num2str(way),' = []'))
                end
                CR(subarray,way) = 0; %CR = 0 is the default value
            end %K
            
            if getDebug
                disp('FSb(sets,num_of_fblks):')
                disp(FSb)
            end
            
            if getDebug
                disp(strcat('CR(',num2str(subarray),',',num2str(way),'):',num2str(CR(subarray,way))))
            end
            
        end %way
        
        if getDebug
            for way=1:assoc
                disp(strcat('CR(',num2str(subarray),',',num2str(way),'):',num2str(CR(subarray,way))))
                disp(strcat('FBw',num2str(way),' = [',num2str(FBw(:,way)'),']'))
            end
        end
        
        % find the Fully Faulty Sets after permutation
        FFS = zeros(subsets,1);
        FFS_num = 0;
        for set_addr=0:subsets-1
            FB_count = 0;
            for way=1:assoc
                FBw_end = FBw_num(1,way);
                if any(set_addr == FBw(1:FBw_end,way))
                    % if value exists in FBw(way)
                    FB_count = FB_count + 1;
                end
            end
            if FB_count == assoc
                FFS_num = FFS_num + 1;
                FFS(FFS_num,1) = set_addr; % insert #set to FFS_init
            end
        end
        
        if FFS_num > 0
            FFS_semifinal = FFS(1:FFS_num,1);
        else
            FFS_semifinal = [];
        end
        
        % fix new FFSets to match with initial FFSets
        P = numel(FFS_before);
        R = numel(FFS_semifinal);
        if (P > 0 && R > 0)
            CNT_XFFS = zeros(subsets,1); % Repeats of each possible result
            for p=1:P
                for r=1:R
                    XFFS = bitxor(FFS_before(p),FFS_semifinal(r)); % XFFS = FFS_before XOR FFS_after
                    CNT_XFFS(XFFS+1,1) = CNT_XFFS(XFFS+1,1) + 1; % Repeats++
                end
            end
            
            [~,xffs_min_pos] = max(CNT_XFFS); % CR-0 value
            CR_offset = S(xffs_min_pos);
            for way=1:assoc
                CR(subarray,way) = bitxor(CR(subarray,way),CR_offset); % CR-way = CR-way XOR CR-0
                FBw(:,way) = bitxor(FBw(:,way),CR_offset); % FBw(way) = FBw(way) XOR CR-0
            end
        end
        
        % find the Fully Faulty Sets after permutation & fix
        FFS = zeros(subsets,1);
        FFS_num = 0;
        for set_addr=0:subsets-1
            FB_count = 0;
            for way=1:assoc
                FBw_end = FBw_num(1,way);
                if any(set_addr == FBw(1:FBw_end,way))
                    % if value exists in FBw(way)
                    FB_count = FB_count + 1;
                end
            end
            if FB_count == assoc
                FFS_num = FFS_num + 1;
                FFS(FFS_num,1) = set_addr; % insert #set to FFS_init
            end
        end
        
        if FFS_num > 0
            FFS_final = FFS(1:FFS_num,1);
        else
            FFS_final = [];
        end
        
        if getDebug
            disp('FFSETS Before Permutation:')
            if numel(FFS_before)
                disp(FFS_before')
            else
                disp('None')
            end
            disp('FFSETS After Permutation:')
            if numel(FFS_semifinal)
                disp(FFS_semifinal')
            else
                disp('None')
            end
            disp('FFSETS After Permutation & Fix:')
            if numel(FFS_final)
                disp(FFS_final')
            else
                disp('None')
            end
        end
        
        % number of common Full Faulty Sets
        common_ffsets = common_ffsets + numel(intersect(FFS_before,FFS_final));
        ffsets_num_before = ffsets_num_before + numel(FFS_before);
        ffsets_num_after = ffsets_num_after + numel(FFS_final);
    end %subarray
    
    if getDebug
        disp('Final CR(subarray,way):')
        disp(CR)
        disp(strcat('Number of Fully Faulty Sets Before: [',num2str(ffsets_num_before),']'))
        disp(strcat('Number of Fully Faulty Sets After: [',num2str(ffsets_num_after),']'))
        disp(strcat('Number of Common Fully Faulty Sets Before/After: [',num2str(common_ffsets),']'))
    end
    
end %end of function way_permutation()

function [CR,ffsets_num_before,ffsets_num_after,common_ffsets] = group_way_permutation(btb,sets,assoc,blk_size,CRsPerWay)
    
    CR=zeros(CRsPerWay,assoc);
    subsets = sets/CRsPerWay;
    S = 0:(subsets-1);
    common_ffsets = 0;
    ffsets_num_before = 0;
    ffsets_num_after = 0;
    
    for subarray=1:CRsPerWay
        subsets_start = (subarray-1)*subsets;
        
        % FBw(i,j): Faulty block i @ way j
        % FBw_num(1,j): Number of faulty blocks @ way j
        FBw = zeros(subsets,assoc);
        FBw_num = zeros(1,assoc); % max_val = subsets
        
        % find faulty physical addresses per way
        for way=1:assoc
            for set=1:subsets
                errors = 0;
                for sbit=1:blk_size
                    if (btb(subsets_start+set,(way-1)*blk_size+sbit) == 1)
                        errors = 1;
                    end
                end
                if errors == 1
                    FBw_num(1,way) = FBw_num(1,way) + 1;
                    FBw_pnt = FBw_num(1,way);
                    FBw(FBw_pnt,way) = set - 1;
                end
            end
            if getDebug
                disp(strcat('FBw',num2str(way),' = [',num2str(FBw(:,way)'),']'))
            end
        end
        
        % find the Fully Faulty Sets initially
        FFS = zeros(subsets,1);
        FFS_num = 0;
        for set_addr=0:subsets-1
            FB_count = 0;
            for way=1:assoc
                FBw_end = FBw_num(1,way);
                if any(set_addr == FBw(1:FBw_end,way))
                    % if value exists in FBw(way)
                    FB_count = FB_count + 1;
                end
            end
            if FB_count == assoc
                FFS_num = FFS_num + 1;
                FFS(FFS_num,1) = set_addr; % insert #set to FFS_init
            end
        end
        
        if FFS_num > 0
            FFS_before = FFS(1:FFS_num,1);
        else
            FFS_before = [];
        end
        
        group_permut_steps = log2(assoc);
        groups = assoc;
        for permut_step=1:group_permut_steps
            % Permutation between pairs of groups
            group_ways = assoc/groups;
            pair_ways = group_ways * 2;
            
            for group_pair=1:2:groups
                % FGSb_left(i,j): Faulty set i with j faulty blocks @ left group
                % FGSb_num_left(1,j): Number of faulty sets with j faulty blocks @ left group
                % FGSb_left(i,j): Faulty set i with j faulty blocks @ right group
                % FGSb_num_left(1,j): Number of faulty sets with j faulty blocks @ right group
                FGSb_left = zeros(subsets,group_ways);
                FGSb_num_left = zeros(1,group_ways);
                FGSb_right = zeros(subsets,group_ways);
                FGSb_num_right = zeros(1,group_ways);
                
                group_leader = group_pair*group_ways + 1;
                if getDebug
                    disp(strcat('Permutation Step [',num2str(permut_step),']: Calculation of CR(',num2str(subarray),',',num2str(group_leader),')'))
                end
                
                for set_addr=0:subsets-1
                    FW_count_left = 0;
                    FW_count_right = 0;
                    for subway=1:group_ways
                        left_way = (group_pair-1)*group_ways + subway;
                        right_way = group_pair*group_ways + subway;
                        FBw_end_left = FBw_num(1,left_way);
                        FBw_end_right = FBw_num(1,right_way);
                        if any(set_addr == FBw(1:FBw_end_left,left_way))
                            FW_count_left = FW_count_left + 1;
                        end
                        if any(set_addr == FBw(1:FBw_end_right,right_way))
                            FW_count_right = FW_count_right + 1;
                        end
                    end %subway
                    if FW_count_left > 0
                        FGSb_num_left(1,FW_count_left) = FGSb_num_left(1,FW_count_left) + 1;
                        FGSb_pnt = FGSb_num_left(1,FW_count_left);
                        FGSb_left(FGSb_pnt,FW_count_left) = set_addr; % insert #set to FGSb_left
                    end
                    if FW_count_right > 0
                        FGSb_num_right(1,FW_count_right) = FGSb_num_right(1,FW_count_right) + 1;
                        FGSb_pnt = FGSb_num_right(1,FW_count_right);
                        FGSb_right(FGSb_pnt,FW_count_right) = set_addr; % insert #set to FGSb_right
                    end
                end %set_addr
                
                % CNT_XFP: Count repeats of XORing results for each pair
                CNT_XFP = zeros(subsets,pair_ways);
                % for each total number of faulty blocks in group sets
                for left_group_fblk=group_ways:-1:1
                    left_fblk_num = FGSb_num_left(1,left_group_fblk);
                    FP_left = FGSb_left(1:left_fblk_num,left_group_fblk);
                    
                    for right_group_fblk=group_ways:-1:1
                        right_fblk_num = FGSb_num_right(1,right_group_fblk);
                        FP_right = FGSb_right(1:right_fblk_num,right_group_fblk);
                        
                        pair_set_fblks = left_group_fblk + right_group_fblk;
                        
                        if getDebug
                            disp(strcat('FPsw',num2str(left_group_fblk),'_left = [',num2str(FP_left'),']'))
                            disp(strcat('FPsw',num2str(right_group_fblk),'_right = [',num2str(FP_right'),']'))
                            disp(strcat('XFP = FPsw',num2str(left_group_fblk),'_left XOR FPsw',num2str(right_group_fblk),'_right'))
                        end
                        
                        for j=1:left_fblk_num
                            for k=1:right_fblk_num
                                XFP = bitxor(FP_left(j),FP_right(k)); % XFP = FPw(left_group) XOR FPw(right_group)
                                CNT_XFP(XFP+1,pair_set_fblks) = CNT_XFP(XFP+1,pair_set_fblks) + 1; % Repeats++
                            end
                        end
                    end %right_group_fblk
                end %left_group_fblk
                
                CCR = S; % Candidates for CR = S
                
                for pair_set_fblks = pair_ways:-1:2
                    % from previous CCR values choose the one that
                    % respect CCR(k) column in CNT_XFP has min value
                    
                    % EXAMPLE:
                    % CNT_XFP = {D0,D1,D2,D3,D4,D5,D6,D7} = {8,1,2,7,3,2,0,1}
                    % CCR = {0,2,4,5}
                    % CNT_CCR = {D0,D2,D4,D5} = {8,2,3,2} -> CNT_MIN = 2, CCR_NUM = 2
                    % CCR_new = {2,5}
                    
                    CNT_CCR = zeros(numel(CCR),1);
                    for cr_num=1:numel(CCR)
                        CNT_CCR(cr_num,1) = CNT_XFP(CCR(cr_num)+1,pair_set_fblks);
                    end
                    CNT_MIN = min(CNT_CCR); % Least repeats
                    CCR_NUM = histc(CNT_CCR,CNT_MIN); % How many values have the same number of least repeats
                    
                    CCR_new = zeros(1,CCR_NUM); % Candidates for CR
                    cnt_pnt = 1;
                    for cr_num=1:numel(CCR)
                        if (CNT_CCR(cr_num,1) == CNT_MIN)
                            CCR_new(cnt_pnt) = CCR(cr_num); % Least repeated values
                            cnt_pnt = cnt_pnt + 1;
                        end
                    end
                    if (CCR_NUM == 0) || (numel(CCR_new) ~= CCR_NUM)
                        ErrMsg = 'Something went wrong with CR calculation';
                        error(ErrMsg)
                    end
                    
                    if getDebug
                        disp(strcat('Pair-set Faulty Blocks = [',num2str(pair_set_fblks),']'))
                        disp(strcat('CCR_input = [',num2str(CCR(:)'),'] -> [',num2str(numel(CCR)),'] elements'))
                        disp(strcat('CNT_XFP = [',num2str(CNT_XFP(:,pair_set_fblks)'),'] -> [',num2str(numel(CNT_XFP(:,pair_set_fblks))),'] elements'))
                        disp(strcat('CNT_CCR = [',num2str(CNT_CCR(:)'),'] -> [',num2str(numel(CNT_CCR)),'] elements'))
                        disp(strcat('CNT_MIN = [',num2str(CNT_MIN),'], CCR_NUM = [',num2str(CCR_NUM),']'))
                        disp(strcat('CCR_output = [',num2str(CCR_new(:)'),'] -> [',num2str(numel(CCR_new)),'] elements'))
                    end
                    
                    CCR = CCR_new;
                end %pair_fblk
                
                CR_value = CCR(1,1); % choose the first one
                %CR_value = randsample(CCR,1); % choose randomly
                
                if getDebug
                    disp(strcat('CR_value=',num2str(CR_value)))
                end
                
                % Update right group CRs
                for subway=1:group_ways
                    way = group_pair*group_ways + subway;
                    CR(subarray,way) = bitxor(CR(subarray,way),CR_value);
                    % FBw(way) = FBw(way) XOR CR(way)
                    FBw(:,way) = bitxor(FBw(:,way),CR_value);
                    if getDebug
                        disp(strcat('CR(',num2str(subarray),',',num2str(way),'):',num2str(CR(subarray,way))))
                    end
                end
                
            end %group_pair
            
            groups = groups/2;
        end %group_step
        
        if getDebug
            for way=1:assoc
                disp(strcat('CR(',num2str(subarray),',',num2str(way),'):',num2str(CR(subarray,way))))
                disp(strcat('FBw',num2str(way),' = [',num2str(FBw(:,way)'),']'))
            end
        end
        
        % find the Fully Faulty Sets after permutation
        FFS = zeros(subsets,1);
        FFS_num = 0;
        for set_addr=0:subsets-1
            FB_count = 0;
            for way=1:assoc
                FBw_end = FBw_num(1,way);
                if any(set_addr == FBw(1:FBw_end,way))
                    % if value exists in FBw(way)
                    FB_count = FB_count + 1;
                end
            end
            if FB_count == assoc
                FFS_num = FFS_num + 1;
                FFS(FFS_num,1) = set_addr; % insert #set to FFS_semifinal
            end
        end
        
        if FFS_num > 0
            FFS_semifinal = FFS(1:FFS_num,1);
        else
            FFS_semifinal = [];
        end
        
        % fix new FFSets to match with initial FFSets
        P = numel(FFS_before);
        R = numel(FFS_semifinal);
        if (P > 0 && R > 0)
            CNT_XFFS = zeros(subsets,1); % Repeats of each possible result
            for p=1:P
                for r=1:R
                    XFFS = bitxor(FFS_before(p),FFS_semifinal(r)); % XFFS = FFS_before XOR FFS_after
                    CNT_XFFS(XFFS+1,1) = CNT_XFFS(XFFS+1,1) + 1; % Repeats++
                end
            end
            
            [~,xffs_min_pos] = max(CNT_XFFS); % CR-0 value
            CR_offset = S(xffs_min_pos);
            for way=1:assoc
                CR(subarray,way) = bitxor(CR(subarray,way),CR_offset); % CR-way = CR-way XOR CR-0
                FBw(:,way) = bitxor(FBw(:,way),CR_offset); % FBw(way) = FBw(way) XOR CR-0
            end
        end
        
        % find the Fully Faulty Sets after permutation & fix
        FFS = zeros(subsets,1);
        FFS_num = 0;
        for set_addr=0:subsets-1
            FB_count = 0;
            for way=1:assoc
                FBw_end = FBw_num(1,way);
                if any(set_addr == FBw(1:FBw_end,way))
                    % if value exists in FBw(way)
                    FB_count = FB_count + 1;
                end
            end
            if FB_count == assoc
                FFS_num = FFS_num + 1;
                FFS(FFS_num,1) = set_addr; % insert #set to FFS_final
            end
        end
        
        if FFS_num > 0
            FFS_final = FFS(1:FFS_num,1);
        else
            FFS_final = [];
        end
        
        if getDebug
            disp('FFSETS Before Permutation:')
            if numel(FFS_before)
                disp(FFS_before')
            else
                disp('None')
            end
            disp('FFSETS After Permutation:')
            if numel(FFS_semifinal)
                disp(FFS_semifinal')
            else
                disp('None')
            end
            disp('FFSETS After Permutation & Fix:')
            if numel(FFS_final)
                disp(FFS_final')
            else
                disp('None')
            end
        end
        
        % number of common Full Faulty Sets
        common_ffsets = common_ffsets + numel(intersect(FFS_before,FFS_final));
        ffsets_num_before = ffsets_num_before + numel(FFS_before);
        ffsets_num_after = ffsets_num_after + numel(FFS_final);
    end %subarray
    
    if getDebug
        disp('Final CR(subarray,way):')
        disp(CR)
        disp(strcat('Number of Fully Faulty Sets Before: [',num2str(ffsets_num_before),']'))
        disp(strcat('Number of Fully Faulty Sets After: [',num2str(ffsets_num_after),']'))
        disp(strcat('Number of Common Fully Faulty Sets Before/After: [',num2str(common_ffsets),']'))
    end
end %end of function group_way_permutation()

function [] = print_sbitmap(btb,btb_entries,sets,assoc,blk_size,p_number,set_bits,nfaults,output_fmap_dir,btb_type,common_ffsets)
    
    total_er = 0;
    fblks = 0;
    fset_fblk = zeros(1,assoc);
    way_fblk = zeros(1,sets);
    
    [fname,ErrMsg] = fopen(strcat(output_fmap_dir,'create-subitmap',num2str(p_number-1),'_',num2str(btb_type),'.cc'),'a');
    error(ErrMsg);
    
    for set=1:sets
        fblk_set = 0;
        for way=1:assoc
            errors = 0;
            for sbit=1:blk_size
                if (btb(set,(way-1)*blk_size+sbit) == 1)
                    errors = 1;
                    total_er = total_er + 1;
                end
            end
            if (errors == 1)
                fprintf(fname,'        sbitmap[%d][%d] = 1;\n',(set-1),(way-1));
                fblks = fblks + 1;
                fblk_set = fblk_set + 1;
                way_fblk(way) = way_fblk(way) + 1; % number of faulty blocks per way
            else
                fprintf(fname,'        sbitmap[%d][%d] = 0;\n',(set-1),(way-1));
            end
        end
        if (fblk_set>0)
            fset_fblk(fblk_set) = fset_fblk(fblk_set) + 1; % number of sets with #fblk_set faulty blocks
        end
    end
    
    fprintf(fname,'    }\n'); % assoc end
    fprintf(fname,'      // btb_entries=%d, assoc=%d, sets=%d, blk_size=%d(bits), set_bits=%d\n',btb_entries, assoc, sets, blk_size, set_bits);
    fprintf(fname,'      // totaler=%d, data_faults=%d, fblks=%d\n',total_er,nfaults(p_number),fblks);
    fprintf(fname,'      // total_ffsets=%d, new_ffsets=%d, common_ffsets=%d\n',fset_fblk(assoc),fset_fblk(assoc)-common_ffsets,common_ffsets);

    for i=1:assoc
        if i==1,
            fprintf(fname,'      // ');
        end
        fprintf(fname,'fblks_way_%d=%d',i,way_fblk(i));
        if i<assoc
            fprintf(fname,', ');
        end
    end
    fprintf(fname,'\n');
    
    for i=1:assoc
        if i==1,
            fprintf(fname,'      // ');
        end
        fprintf(fname,'fset_%d_fblk=%d',i,fset_fblk(i));
        if i<assoc
            fprintf(fname,', ');
        end
    end
    fprintf(fname,'\n');
    
    fclose(fname);
    
end %end of function print_sbitmap()
