elements = 7
pages = 2
array_A = rand(3, elements, pages);
array_B = rand(1, elements, pages);
array_B = [101 98 95 92 91 88 85];
array_B(:,:,2) = [12 14 16 18 20 22 24];

[array_C, array_mag_I] = sort(array_B);

array_D = diff(array_C,1,2);
array_C
array_D

% append(string(array_mag_index),"   ", string(array_mag_sort))
% for i = (1:elements)
%     array_mag_sort(i),array_mag_index(i)
% end