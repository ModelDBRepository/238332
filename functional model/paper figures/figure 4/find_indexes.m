
function idxs = find_indexes(num)

type_number = [0,126,104,136,384,236,338,58]; 
for i=1:8
    v(i)=sum(type_number(1:i));
end

idxs = (v(num)+1):v(num+1);  