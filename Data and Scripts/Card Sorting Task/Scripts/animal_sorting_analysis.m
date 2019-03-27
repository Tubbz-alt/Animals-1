
%% Last edited 03/12/19
% Judy Sein Kim 
% Analyses for animal card sorting data 

% 1-objects, 2-habitat, 3-food, 4-shape, 5-skin, 6-color
tasks={'objects' 'habitat' 'food' 'shape' 'skin' 'color'}; 
load('animal_keys.mat')

%% PART 1: INDIVIDUAL-TO-GROUP CORRELATIONS %%

group = 'CB'; load('CB_subjects.mat'); filename = 'animals_sorting_CB_allData.mat';
%group = 'S'; load('S_subjects.mat'); filename = 'animals_sorting_S_allData.mat'; 
 
for i = 1:numel(tasks) 
    itemN = 30; 
    names = animalNames;
    if strcmp(tasks{i},'objects')
       itemN = 29;
       names = toolNames; 
    end
    subs = subsAll{i}; 
    group_matrix = zeros(itemN,itemN); 
    clear ind_matrix_all 
    for j  = 1:numel(subs) 
       % Load individual subject data 
       fid = fopen(sprintf('%s_%s_%s.csv',group,subs{j},tasks{i})); 
       data = textscan(fid,'%s %s %s','Delimiter',','); 
       org_vector = str2num(cell2mat(data{2}(2:end))); 
       ind_vector = org_vector; 
           
       % Turn vector into matrix 
       ind_matrix = tril(ones(itemN)) - eye([itemN itemN]);
       ind_matrix(logical(ind_matrix)) = ind_vector;
       ind_matrix = ind_matrix + ind_matrix' + eye([itemN itemN]);
           
       ind_matrix_all{j} = ind_matrix; 
       group_matrix = group_matrix + ind_matrix; 
    end
    group_matrix = group_matrix/numel(subs); 
    group_matrix_all{i} = group_matrix; 
    ind_matrix_all_all{i} = ind_matrix_all; 
       
%% A. Within-group, within-dimension, individual-to-group, leave-one-out correlations (B-B and S-S)

    spearRho = zeros(numel(subs),1); spearP = zeros(numel(subs),1); spearRho_Z = zeros(numel(subs),1); 
    for k = 1:numel(subs)
       leave_matrix = ind_matrix_all{k};
           
       group_LOO = ind_matrix_all;
       group_LOO{k} = zeros([itemN itemN]); 
       group_LOO_new = reshape(cell2mat(group_LOO),itemN,[],numel(subs));
       group_LOO_new = sum(group_LOO_new,3)/numel(subs);

       [I] = itril(size(leave_matrix,1),-1); 
       [spearRho(k),spearP(k)] = corr(group_LOO_new(I),leave_matrix(I),'type','Spearman');
       spearRho_Z(k) = 0.5*log((1+spearRho(k))/(1-spearRho(k))); % Fisher's Z-transform
    end   
    spearRho_all{i} = spearRho; 
    spearRhoZ_all{i} = spearRho_Z; 
     
end

save(filename,'group_matrix_all','spearRho_all','spearRhoZ_all','ind_matrix_all_all') 

%% B. Within-dimension, blind-individual-to-sighted-group correlations (B-S) 

load('animals_sorting_S_allData.mat') 
load('S_subjects.mat')
sightedSubs = subsAll;
sighted_group_all = group_matrix_all; 
sighted_ind_all = ind_matrix_all_all; 

load('animals_sorting_CB_allData.mat') 
load('CB_subjects.mat')
blindSubs = subsAll;
blind_group_all = group_matrix_all; 
blind_ind_all = ind_matrix_all_all; 

for i = 1:numel(tasks) 
    subs = blindSubs{i};  
    spearRho = zeros(numel(subs),1); spearP = zeros(numel(subs),1); spearRho_Z = zeros(numel(subs),1); 
    for j = 1:numel(subs) 
        [I] = itril(size(blind_ind_all{i}{j},1),-1); 
        [spearRho(j),spearP(j)] = corr(sighted_group_all{i}(I),blind_ind_all{i}{j}(I),'type','Spearman');
        spearRho_Z(j) = 0.5*log((1+spearRho(j))/(1-spearRho(j)));
    end
    BtoS_rho{i} = spearRho; 
    BtoS_rhoZ{i} = spearRho_Z; 
end

%% C. Within-dimension, individual-subject-to-taxonomy correlations (reported in Supplementals) 

load('taxonomy_matrix.mat')

ind_mats_both_groups{1} = sighted_ind_all; 
ind_mats_both_groups{2} = blind_ind_all; 
subs_both_groups{1} = sightedSubs;
subs_both_groups{2} = blindSubs;

for i = 1:2
    ind_mats = ind_mats_both_groups{i}; 
    subsAll = subs_both_groups{i};
    for j = 4:6 % only need to do this for shape, texture, color 
        subs = subsAll{j};  
        spearRho = zeros(numel(subs),1); spearP = zeros(numel(subs),1); spearRho_Z = zeros(numel(subs),1); 
        for k = 1:numel(subs) 
            [I] = itril(size(ind_mats{j}{k},1),-1); 
            [spearRho(k),spearP(k)] = corr(tax(I),ind_mats{j}{k}(I),'type','Spearman');
            spearRho_Z(k) = 0.5*log((1+spearRho(k))/(1-spearRho(k)));
        end
    tax_rho{i}{j} = spearRho; 
    tax_rhoZ{i}{j} = spearRho_Z; 
    end
end

%% PART 2: GROUP-TO-GROUP CORRELATIONS %%

% A. Across-group, within-dimension correlations
% diagnoals should be 0? not 1? 

zM_group = zeros(3,1); 
zN_group = zeros(3,1); 
rho_group = zeros(3,1); 
p_group = zeros(3,1);

for i = 4:6 
    [zM_group(i),zN_group(i),rho_group(i),p_group(i)] = mantel_test(sighted_group_all{i},blind_group_all{i},1000); 
end

% B. Within-group, across-dimension correlations 
% 1-shape vs. texture, 2-texture vs. color, 3-shape vs. color 

zM_S = zeros(3,1); zM_CB = zeros(3,1); 
zN_S = zeros(3,1); zN_CB = zeros(3,1);
rho_S = zeros(3,1); rho_CB = zeros(3,1); 
p_S = zeros(3,1); p_CB = zeros(3,1); 

[zM_S(1),zN_S(1),rho_S(1),p_S(1)] = mantel_test(sighted_group_all{4},sighted_group_all{5},1000); 
[zM_S(2),zN_S(2),rho_S(2),p_S(2)] = mantel_test(sighted_group_all{4},sighted_group_all{6},1000); 
[zM_S(3),zN_S(3),rho_S(3),p_S(3)] = mantel_test(sighted_group_all{5},sighted_group_all{6},1000);

[zM_CB(1),zN_CB(1),rho_CB(1),p_CB(1)] = mantel_test(blind_group_all{4},blind_group_all{5},1000); 
[zM_CB(2),zN_CB(2),rho_CB(2),p_CB(2)] = mantel_test(blind_group_all{4},blind_group_all{6},1000); 
[zM_CB(3),zN_CB(3),rho_CB(3),p_CB(3)] = mantel_test(blind_group_all{5},blind_group_all{6},1000);

% C. Correlations with taxonomy 

zM_S_tax = zeros(3,1); zM_CB_tax = zeros(3,1); 
zN_S_tax = zeros(3,1); zN_CB_tax = zeros(3,1);
rho_S_tax = zeros(3,1); rho_CB_tax = zeros(3,1); 
p_S_tax = zeros(3,1); p_CB_tax = zeros(3,1); 

for i = 4:6 
    %[zM_S_tax(i),zN_S_tax(i),rho_S_tax(i),p_S_tax(i)] = mantel_test(sighted_group_all{i},tax,1000); 
    [zM_CB_tax(i),zN_CB_tax(i),rho_CB_tax(i),p_CB_tax(i)] = mantel_test(blind_group_all{i},tax,1000); 
end
       