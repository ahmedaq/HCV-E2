function WF_E2_script_Ng_b_beta(Ng,b,beta_param,imm_pres_site,strain_num)

pobj = parpool(16);

% imm_pres_site:    1 to 352 [can also input multiple sites]
% strain_num:       strain number to initialize the pop (select index of random strains generated using randperm with rng(0)

load weight_sequences_E2.mat
load test886_hepc_e2_send_Ahmed.mat
H = reshape(J_mat,length(single_double_mutant_mat),length(single_double_mutant_mat));
H = convertRayFLParamsToJohnParamsFormat(H);

phi_curr = num_mutants_combine_array;
mutant_order = amino_single_combine_array;

phi_cumulative(1) = phi_curr(1);
for kk = 2:length(mutant_order)
    phi_cumulative(kk) = sum(phi_curr(1:kk));
end

phi_cum = [0 phi_cumulative];

%
% msa_aa_ex = msa_bin;
msa = msa_aa(:,ind_non_conserve);

% % Cseq issue
% Cseq_aa = seqconsensus_aa(msa);
% for kk = 1:size(msa,2)
%     Cseq_mut_order(kk) = mutant_order{kk}(1);
% end
% find(Cseq_mut_order~=Cseq_aa)%different because of patient weights
% %thus we use Cseq_mut_order as the consensus seq
% save Cseq_mut_order Cseq_mut_order
load Cseq_mut_order Cseq_mut_order
Cseq_aa = Cseq_mut_order;

%Modified the genetic code 23: same as standard (code = 1) except changes in code.B, code.X, and code.stop
% msa_nt = aa2nt_ahmed(msa, 23);
% save msa_nt msa_nt
load msa_nt msa_nt

func_working = ~sum(sum(msa~=nt2aa_ahmed(msa_nt, 23)));
if func_working == 1
    fprintf('\n---Modified aa2nt and nt2aa functions working as expected---\n')
end

%Generating bin_matrix (definition of binary encoding)
for kk=1:max(phi_curr)+1
    temp_matrix=[];
    temp_matrix = fliplr(eye(kk));
    temp_matrix = [zeros(1,size(temp_matrix,2)) ; temp_matrix ];
    bin_matrix{kk} =temp_matrix;
end


N = 2000;           %Number of protein sequences in each generation
M = size(msa_nt,2);            %Number of sites in the protein
mu = 1e-4;          %Mutation rate per genome per replication

% Ng = 500; %4e3;     %Number of generations the simulation will be run

imm_pres_gen = 0;
% b = 10;%12;
% beta_param = 0.1;
freq_escape = 0.5;
ITER = 96;

Cseq_nt = aa2nt_ahmed(Cseq_aa, 23);

%%

%sites in extended binary under immune pressure
imm_pres_site_bin = [];
for kk = 1:length(imm_pres_site)
    imm_pres_site_bin = [imm_pres_site_bin phi_cum(imm_pres_site(kk))+1:phi_cum(imm_pres_site(kk)+1)];
end
imm_pres_site_bin


%% find seqs with WT at a particular position

msa_unique = unique(msa,'rows');
msa_nt_unique = aa2nt_ahmed(msa_unique,23);

for kk = 1:size(msa_unique,2)
    seqs_wt{kk} = find(msa_unique(:,kk)==mutant_order{kk}(1)); %k-th cell contains unique MSA sequence indices with WT at k-th position
    no_seqs_wt(kk) = length(seqs_wt{kk});
end

rng(0) %initializing the random number generator

%selecting 100 random MSA sequences with WT at imm_pres_site

rand_perm = randperm(no_seqs_wt(imm_pres_site));
indx_seq = seqs_wt{imm_pres_site}(rand_perm(strain_num));

init_pop = repmat(msa_nt_unique(indx_seq,:),N,1);


%%

parfor iter = 1:ITER
    
    init_pop_aa = nt2aa_ahmed(init_pop, 23);
    [init_pop_aa_unique,count_init_pop_aa] = unique_seqs(init_pop_aa);
    
    init_pop_bin_unique = zeros(size(init_pop_aa_unique,1),length(H));
    energy_init_pop_bin_unique = zeros(size(init_pop_aa_unique,1),1);
    energy_init_pop_bin = [];
    for kk = 1:size(init_pop_aa_unique,1)
        init_pop_bin_unique(kk,:) = convertAAseq2Bin_ahmed(mutant_order,bin_matrix,init_pop_aa_unique(kk,:));
        energy_init_pop_bin_unique(kk) = diag(init_pop_bin_unique(kk,:)*triu(-H)*init_pop_bin_unique(kk,:).').';
        energy_init_pop_bin = [energy_init_pop_bin; repmat(energy_init_pop_bin_unique(kk),count_init_pop_aa,1)];
    end
    
    energy_pop = cell(1,Ng+1);
    Ebar = zeros(1,Ng+1);
    Ebar(1) = mean(energy_init_pop_bin);
    
    Ebar_without_imm_pres = zeros(1,Ng+1);
    Ebar_without_imm_pres(1) = mean(energy_init_pop_bin);
    
    freq_imm_pres_site = zeros(Ng+1,1);
    freq_imm_pres_site(1,:) = mean(init_pop_bin_unique(:,imm_pres_site_bin));
    
    pop = init_pop;
    
    for kk = 1:Ng
        
        %%%%%%%%%%%%% MUTATION AT NT LEVEL %%%%%%%%%%%%%%%%%%
        mut_step = binornd(1,mu,N,M);
        
        %selecting sequences in which mutation is to be made
        seqs_with_mut = find(sum(mut_step,2).'~=0);
        
        %making nt mutations in the selected sequences randomly from ACGT
        mut_pop = pop;
        for mm = 1:length(seqs_with_mut)
            
            seq_sel = seqs_with_mut(mm);
            pos_mut = find(mut_step(seq_sel,:)~=0); %there can be multiple mutations per sequence; hence, next loop!
            
            for nn = 1:length(pos_mut)
                
                %nt_curr = pop(seq_sel,pos_mut(nn));
                
                %forming pool of nt available at this position, i.e., setdiff('ACGT',nt_curr)
                pool_nt = setdiff('ACGT',pop(seq_sel,pos_mut(nn)));
                %nt_new = pool_nt(randi(3));
                
                %selecting a random nt from the available pool
                mut_pop(seq_sel,pos_mut(nn)) = pool_nt(randi(3));
            end
            
        end
        
        mut_pop_aa = nt2aa_ahmed(mut_pop, 23);
        [mut_pop_aa_unique,count_mut_pop_aa] = unique_seqs(mut_pop_aa);
        
        mut_pop_bin_unique = zeros(size(mut_pop_aa_unique,1),length(H));
        energy_mut_pop_bin_unique = zeros(size(mut_pop_aa_unique,1),1);
        energy_mut_pop_bin = [];
        
        if kk>imm_pres_gen
            bb = b;
        else
            bb = 0;
        end
        
        count_mut_pop = 0;
        for kkk = 1:size(mut_pop_aa_unique,1)
            mut_pop_bin_unique(kkk,:) = convertAAseq2Bin_ahmed(mutant_order,bin_matrix,mut_pop_aa_unique(kkk,:));
            
            %imm_pressure apply or not
            if sum(xor(mut_pop_bin_unique(kkk,imm_pres_site_bin),zeros(1,length(imm_pres_site_bin))),2)==0
                energy_mut_pop_bin_unique(kkk) = mut_pop_bin_unique(kkk,:)*triu(-H)*mut_pop_bin_unique(kkk,:).' + bb;
            else
                energy_mut_pop_bin_unique(kkk) = mut_pop_bin_unique(kkk,:)*triu(-H)*mut_pop_bin_unique(kkk,:).';
                count_mut_pop = count_mut_pop+count_mut_pop_aa(kkk);
            end
            
            energy_mut_pop_bin = [energy_mut_pop_bin; repmat(energy_mut_pop_bin_unique(kkk),count_mut_pop_aa(kkk),1)];
        end
        
        freq_imm_pres_site(kk+1,:) = count_mut_pop/N;
        
        %when freq_imm_pres_site becomes >0 for the first time
        if freq_imm_pres_site(kk)==0 && freq_imm_pres_site(kk+1)>0
            escape_mut_appearance = kk+1
        end
        
        Ebar(kk+1) = mean(energy_mut_pop_bin);
        
        %John's metric but does not sum up to 1
        prob_survival_strains = (exp(beta_param*(Ebar(kk+1)-energy_mut_pop_bin_unique)))./...
            (1+exp(beta_param*(Ebar(kk+1)-energy_mut_pop_bin_unique)));
        prob_survival_strains = prob_survival_strains./sum(prob_survival_strains); %for summing to 1
        
        %My modified metric
        %     prob_survival_strains = (exp(beta_param*(Ebar(kk+1)-energy_pop{kk+1})))/sum((exp(beta_param*(Ebar(kk+1)-energy_pop{kk+1}))));
        
        counts_next_gen = mnrnd(N,prob_survival_strains);
        
        pop = [];
        for mm = 1:size(mut_pop_aa_unique,1)
            pop = [ pop;  repmat(aa2nt_ahmed(mut_pop_aa_unique(mm,:), 23),counts_next_gen(mm),1)];
        end
        
        if ismember(kk,[1e3:1e3:1e4])
            fprintf('No. of generations = %d\n',kk)
        end
        
        if sum(freq_imm_pres_site(kk+1,:) > freq_escape) > 0 %any targeted site > freq_escape
            escape_gen = kk+1
            escape_Ebar = Ebar(kk+1)
            escape_f = freq_imm_pres_site(kk+1,:)
            break;
        end
        
        if kk == Ng %last generation
            escape_gen = kk+1;
            escape_Ebar = Ebar(kk+1);
            escape_f = freq_imm_pres_site(kk+1,:)
        end
        
        %     time_gen(kk) = toc;
        
    end
    
    escape_mut_appearance_iter(iter) = escape_mut_appearance
    escape_gen_iter(iter) = escape_gen
    escape_Ebar_iter(iter) = escape_Ebar
    escape_f_iter(iter) = escape_f
end

escape_gen_iter
mean(escape_gen_iter)
escape_f_iter

save(sprintf('N%d_M%d_Ng%d_b%d_beta%.2f_fEscape%.1f_ITER%d_SiteIp%d_Seq%d.mat',...
    N,M,Ng,b,beta_param,freq_escape,ITER,imm_pres_site,strain_num),...
    'escape_mut_appearance_iter','escape_gen_iter','escape_Ebar_iter','escape_f_iter','indx_seq')

delete(pobj)