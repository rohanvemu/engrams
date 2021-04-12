%% Finding Engrams
load imaging_data_training_46.mat
neuron= permute(neuron_network_imaging,[2 1 3]); 
iter = floor(size(neuron, 2)./100);
sum_neuron = zeros(size(neuron, 1), 100, size(neuron, 3));
%sum across blocks to find more accurate local maxima 
for r = 1:iter
    sum_neuron = sum_neuron + neuron(:, (r-1)*100 + 1 : r*100, :);
end
sum_neuron = (sum_neuron>=2560); %convert to maxes find the first one 
%reorder the array so each row is a neuron
fire_seq = zeros(size(sum_neuron, 1), size(sum_neuron, 2), size(sum_neuron, 3));
% storage matrix for when the firing of each neuron occurs
for k = 1:size(sum_neuron, 3)
    for i = 1:size(sum_neuron,1)
        locs = islocalmax(sum_neuron(i, :, k), 'FlatSelection', 'first'); %find the local maxima
        fire_seq(i, :, k) = locs; %convert a fire to a 1
    end
end 
fire_seq = fire_seq(:, 1:100, :); %pick only 1 second = 1 engram
store_sum = sum(fire_seq, 2); %200x1x95, sum of each neuron for each char
%% Decoding Engrams
load imaging_data_text_46_7.mat
A = importdata('char_train_46.txt', '');
B = char(A);
chunk = neuron_network_imaging'; %orient data so row=neuron
num = size(chunk, 2)/100; %find how many 1 sec chunks exist
char_store = zeros(1, num); %preallocated the storing of characters
max_store = zeros(size(chunk, 1), size(chunk, 2));

for j = 1:size(chunk, 1)
    locs = islocalmax(chunk(j, :)); %assign 1 to a peak
    max_store(j, :) = locs; %store the peaks
end

for i=1:num
    look = max_store(:, (i-1)*100 + 1 : i*100); %obtain the chunk 
    sumlook = (sum(look, 2)); %sum each neuron up 
    [r, v] = max(sum(store_sum == sumlook)); %max trainining char matches
    char_store(i) = v; %store max match index into array
end
%% Displaying Final Sequence
disp("The sequence corresponds to characters:")
disp(B(char_store)) %take index and find corresponding character
%% Cumulative Deletion of Electrodes
runiters = 1;
cum_delete = zeros(1, 200);
char_store2 = zeros(1, num);
delete_store = zeros(runiters, 200);


for j = 1:size(chunk, 1)
    locs = islocalmax(chunk(j, :)); %assign 1 to a peak
    max_store(j, :) = locs; %store the peaks
end

max_store_change = max_store;
store_sum_change = store_sum;

randint = 1:200;
count = 1;

for r = randint
%     
idx = find(store_sum(r, :, :) == 1);

%     max_store_change(idx, :) = 10.*ones(length(idx), size(max_store, 2));
%     store_sum_change(idx, :, :) = 5.*ones(length(idx), size(store_sum, 2), size(store_sum, 3));
    
    for i=1:num
        
        look = max_store_change(:, (i-1)*100 + 1 : i*100); %obtain the chunk 
        idx = find(look(r, :) == 1);
        max_store_change(idx, :) = 10.*ones(length(idx), size(max_store, 2));
        sumlook = (sum(look, 2)); %sum each neuron up 
        [b, v] = max(sum(store_sum == sumlook)); %max trainining char matches
        char_store2(i) = v; %store max match index into array
    end 
% disp("The sequence corresponds to characters:")
disp(B(char_store2)) %take index and find corresponding character
cum_delete(count) = 100 - (sum(char_store2 == char_store) / length(char_store) .*100);
count = count+1;
end 

%%
figure(1) 
mean_delete = mean(delete_store);
plot(1:200, cum_delete(1:end))
xlabel("Years after Onset of Alzheimer's")
ylabel("Percentage of Original Message Lost")
title("Electrode Deletion")