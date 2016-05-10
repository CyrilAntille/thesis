function [ new ] = copyStruct( input )
%COPYSTRUCT Summary of this function goes here
%   Detailed explanation goes here
new = feval(class(input));
p = fieldnames(struct(input));
for i = 1:length(p)
    new.(p{i}) = input.(p{i});
end
end

