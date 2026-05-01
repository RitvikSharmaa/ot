### 🔹 Simplex Method

```matlab
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
```

---

### 🔹 Big-M Method

```matlab
clear;
clc;

M = 1e8;

A = [6 8 -1 0 1 0;
     7 12 0 -1 0 1];

B = [100;
     120];

C = [-12 -20 0 0 -M -M];

index = [5 6];

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
    end
    
    [~, LV] = min(ratio);
    index(LV) = EV;
    
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

fprintf('Optimal Basic Variables:\n');
disp(XB);

fprintf('Maximum value:\n');
disp(zmax);

fprintf('Minimum value:\n');
disp(-zmax);
```

---

### 🔹 Least Cost Method

```matlab
clear;
clc;

c = [19 30 50 10;
     70 30 40 60;
     40 8  70 20];

a = [7 9 18];
b = [5 8 7 14];

[m,n] = size(c);

if sum(a) ~= sum(b)
    if sum(a) < sum(b)
        c(end+1,:) = zeros(1,n);
        a(end+1) = sum(b) - sum(a);
    else
        c(:,end+1) = zeros(m,1);
        b(end+1) = sum(a) - sum(b);
    end
end

x = zeros(size(c));
initial_c = c;

while any(a) && any(b)
    
    min_cost = min(c(:));
    [row, col] = find(c == min_cost);
    
    i = row(1);
    j = col(1);
    
    allocation = min(a(i), b(j));
    x(i,j) = allocation;
    
    a(i) = a(i) - allocation;
    b(j) = b(j) - allocation;
    
    if a(i) == 0
        c(i,:) = Inf;
    end
    
    if b(j) == 0
        c(:,j) = Inf;
    end
end

z = 0;
for i = 1:size(x,1)
    for j = 1:size(x,2)
        z = z + initial_c(i,j) * x(i,j);
    end
end

disp(array2table(x));
fprintf('Total Cost = %f\n', z);
```

---

### 🔹 Steepest Descent Method

```matlab
clc;
clear;

syms x1 x2

f = x1 - x2 + 2*x1^2 + 2*x1*x2 + x2^2;

fobj = matlabFunction(f, 'Vars', {[x1 x2]});

grad = gradient(f, [x1 x2]);
gradf = matlabFunction(grad, 'Vars', {[x1 x2]});

H = hessian(f, [x1 x2]);
Hf = matlabFunction(H, 'Vars', {[x1 x2]});

x0 = [0; 0];

tol = 1e-4;
maxiter = 50;
iter = 0;

while norm(gradf(x0)) > tol && iter < maxiter
    
    g = gradf(x0);
    s = -g;
    
    Hval = Hf(x0);
    
    alpha = (g'*g) / (g'*Hval*g);
    
    x0 = x0 + alpha * s;
    
    iter = iter + 1;
end

fprintf('Optimal point = [%f , %f]\n', x0(1), x0(2));
fprintf('Optimal value = %f\n', fobj(x0));
fprintf('Iterations = %d\n', iter);
```

---

## 📌 Note

* These implementations are for **academic learning**
* Least Cost Method gives **initial solution**, not optimal
* Big-M uses artificial variables to handle ≥ constraints

---

## ⭐ If this helped

Give the repo a ⭐ and you're good to go 🚀

```
```
