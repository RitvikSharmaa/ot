# Simplex 

clear;
clc;

A = [1 1 1 0;
     1 -1 0 1];

B = [4;
     2];

C = [3 2 0 0];

index = [3,4];

m = length(B);
n = length(C);

BS = [];
CB = [];

for i = 1:m
    CB(i) = C(index(i));
    BS = [BS A(:,index(i))];
end

y = inv(BS) * A;
XB = inv(BS) * B;

z = CB * y - C;

while min(z) < 0
    
    [~, EV] = min(z);
    
    ratio = [];
    
    for i = 1:m
        if y(i,EV) > 0
            ratio(i) = XB(i) / y(i,EV);
        else
            ratio(i) = 1e8;
        end
    end
    
    if min(ratio) == 1e8
        fprintf('Unbounded solution\n');
        return;
    else
        [~, LV] = min(ratio);
        index(LV) = EV;
    end
    
    BS = [];
    CB = [];
    
    for i = 1:m
        CB(i) = C(index(i));
        BS = [BS A(:,index(i))];
    end
    
    y = inv(BS) * A;
    XB = inv(BS) * B;
    z = CB * y - C;
    
    zmax = CB * XB;
end

fprintf('Optimal solution:\n');
disp(XB);

fprintf('Maximum value of Z:\n');
disp(zmax);
