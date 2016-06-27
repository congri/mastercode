function [id] = get_id(nc)
%put in equation number, get back global node number

[eqs,i] = sort(nc(3,:));

id = [eqs',i'];

init = find(eqs == 1);

id(1:(init-1),:) = [];

id = id(:,2);
%convert to integer for efficient memory use
id = int16(id);
end

