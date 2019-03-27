
%load('animals_sighted_labels_121018.mat') 
load('animals_blind_labels_031719.mat') 
rounds={'color' 'texture' 'shape'}; 
itemN=30; 
subNs=[19 19 19]; 
dataAll = {color texture shape};

for i = 1:numel(rounds)
    data = dataAll{i}; 
    subN = subNs(i); 
    for j = 1:itemN
        indItem = data(subN*(j-1)+1:subN*j,:); 
        uniqLabels = unique(indItem); 
        counts = zeros(1,numel(uniqLabels)-1); 
        for k = 2:numel(uniqLabels)
            counts(k) = numel(find(strcmp(indItem,uniqLabels{k}))); 
        end    
        counts = counts(2:end);
        countsAll{j} = counts;
        simpsonD(j) =  sum(counts.*(counts-1))/(sum(counts)*(sum(counts)-1));
    end
    simpsonDAll{i} = simpsonD; 
end