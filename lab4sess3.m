%% Finding Engrams
load prelab_session3_training_chars.mat
neuron= permute(neuron_network_imaging,[2 1 3]); 
%reorder the array so each row is a neuron
fire_seq = zeros(size(neuron, 1), size(neuron, 2), size(neuron, 3));
%storage matrix for when the firing of each neuron occurs
for k = 1:size(neuron, 3)
    for i = 1:size(neuron,1)
        locs = islocalmax(neuron(i, :, k)); %find the local maxima
        fire_seq(i, :, k) = locs; %convert a fire to a 1
    end
end 
fire_seq = fire_seq(:, 1:100, :); %pick only 1 second = 1 engram
store_sum = sum(fire_seq, 2); %200x1x95, sum of each neuron for each char
%% Decoding Engrams
load prelab_session3_sequence.mat
A = importdata('prelab_session3_training_chars.txt', '');
B = char(A);
chunk = neuron_network_imaging'; %orient data so row=neuron
num = size(chunk, 2)/100; %find how many 1 sec chunks exist
char_store = zeros(1, num); %preallocated the storing of characters
max_store = zeros(size(chunk, 1), size(chunk, 2));

for j = 1:size(chunk, 1)
    [pks,locs] = islocalmax(chunk(j, :)); %assign 1 to a peak
    max_store(j, :) = locs; %store the peaks
end

for i=1:num
    look = max_store(:, (i-1)*100 + 1 : i*100); %obtain the chunk 
    sumlook = (sum(look, 2)>=1); %sum each neuron up 
    [r, v] = max(sum(store_sum == sumlook)); %max trainining char matches
    disp([r, v])
    char_store(i) = v; %store max match index into array
end
%% Displaying Final Sequence
disp("The sequence corresponds to characters:")
disp(B(char_store)) %take index and find corresponding character