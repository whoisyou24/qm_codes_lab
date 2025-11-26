#bisection
f[x_] := x^3 - 5 x + 1;
a = 0;
b = 1;
nmax = 15;
Plot[f[x], {x, -1, 1}]
Fori = 1, i ≤ nmax, i++, c = (a + b)  2;
Print"Root after iteration ", i, " is c = ", N[c, 6], " f(c) = ", N[f[c], 6], "\n";
If[f[a] * f[c] < 0, b = c, a = c];

#table
f[x_] := x^3 - 5 x + 1;
a = 0;
b = 1;
nmax = 15;

Plot[f[x], {x, -1, 1}]

For[i = 1, i <= nmax, i++ ,  c = (a + b)/2;
   Print[TableForm[
      {{i, N[c, 4], N[f[c], 4]}},
      TableHeadings -> {None, {"Iteration", "c", "f(c)"}}]];
  If[f[a] * f[c] < 0, b = c, a = c];
]


#secant
f[x_] := x^3 - 5 x + 1;

x0 = 2;
x1 = 3;
nmax = 5;

Plot[f[x], {x, x0, x1}]

For[i = 1, i <= nmax, i++,
  x2 = N[x1 - ((x1 - x0)/(f[x1] - f[x0])) * f[x1]];
  Print["Root after iteration ", i, " is ", x2];
  x0 = x1;
  x1 = x2;
]


#newton Raphson
f[x_] := x^3 + 2 x^2 - 3 x - 1;

a = 1;
b = 2;
nmax = 6;

Plot[f[x], {x, a - 1, b + 1}]

For[i = 1, i <= nmax, i++,
  c = N[b - f[b]/f'[b]];
  Print["The ", i, " th iteration is: ", c];a = b; b = c;
]



#Table
f[x_] := x^3 - 5 x + 1;
x0=0
x1 = 1;     (* initial guess *)
nmax = 6;
Plot[f[x], {x, x0 - 2, x0 + 2}]
For[i = 1, i <= nmax, i++,  
  x2 = x1 - f[x1]/f'[x1];   (* NR formula *)
  Print[TableForm[{{i, N[x2, 5], N[f[x2], 5]}} ]];
  x0 = x1;
 x1=x2;
]



#gauss Seidel
x = 0;
y = 0;
z = 0;
nmax = 10;

For[i = 1, i <= nmax, i++,
  xnew = N[(1/8) (15 - 3 y - 4 z), 5];
  ynew = N[(1/5) (1 + 2 xnew + 2 z), 5];
  znew = N[(1/3) (5 - xnew - ynew), 5];

  Print["The ", i, " th iteration is:  x = ", xnew, "   y = ", ynew,"   z = ", znew];
  x = xnew;
  y = ynew;
  z = znew;
]





#Gauss jacobi matrix
A = {{8, 3, 4},
     {-2, 5, -2},
     {1, 1, 3}};

B = {{15}, {1}, {5}};

nmax = 10;

jacobi[x0_, y0_, z0_, iterations_] :=
 Module[{x, y, z, list, i},
  
  (* initial values *)
  x[0] = x0;
  y[0] = y0;
  z[0] = z0;

  (* table list *)
  list = {{"Iteration", "x", "y", "z"}, {0, x[0], y[0], z[0]}};

  (* Jacobi loop *)
  For[i = 0, i < iterations, i++,
   
   x[i + 1] = N[(B[[1, 1]] - A[[1, 2]]*y[i] - A[[1, 3]]*z[i]) / A[[1, 1]]];
   y[i + 1] = N[(B[[2, 1]] - A[[2, 1]]*x[i] - A[[2, 3]]*z[i]) / A[[2, 2]]];
   z[i + 1] = N[(B[[3, 1]] - A[[3, 2]]*y[i] - A[[3, 1]]*x[i]) / A[[3, 3]]];

   AppendTo[list, {i + 1, x[i + 1], y[i + 1], z[i + 1]}];
  ];

  TableForm[list]
 ]

jacobi[0, 0, 0, nmax]




#lagrange,interpolation
X = {0, 1, 3};
Y = {1, 3, 55};
n = 3;

For[k = 1, k <= n, k++,
  L[n, k, x_] :=
    Product[(x - X[[j]])/(X[[k]] - X[[j]]), {j, 1, k - 1}] *
    Product[(x - X[[j]])/(X[[k]] - X[[j]]), {j, k + 1, n}];  
];

a = Sum[Y[[k]]*L[n, k, x], {k, 1, n}];

Print[a];
Simplify[a]






X = {-2, -1, 0, 1, 3, 4};
Y = {9, 16, 17, 18, 44, 81};
n = 6;

For[k = 1, k <= n, k++,
  L[n, k, x_] :=
    Product[(x - X[[j]])/(X[[k]] - X[[j]]), {j, 1, k - 1}] *
    Product[(x - X[[j]])/(X[[k]] - X[[j]]), {j, k + 1, n}]];

a = Sum[Y[[k]] * L[n, k, x], {k, 1, n}];

Print[a];
Simplify[a]

#newton,interpolation
X = {0, 1, 3};
Y = {1, 3, 55};
n = 2;
d = Table["", {n + 1}, {n + 1}];
d[[All, 1]] = Y[[All]];

(* Construct divided difference table *)
For[j = 1, j <= n, j++,
  For[k = j, k <= n, k++,
    d[[k + 1, j + 1]] = (d[[k + 1, j]] - d[[k, j]]) / (X[[k + 1]] - X[[k + 1 - j]])]];
  

(* p[k+1, x] polynomial term *)
For[k = 0, k <= n, k++,
  p[k + 1, x_] := Product[(x - X[[j]]), {j, 1, k}];
];

(* Final polynomial *)
a = Sum[d[[k + 1, k + 1]] * p[k + 1, x], {k, 0, n}];

Print[a]
Simplify[a]




#integeration
f[x_] := x + 1;
lowerlimit = a = 0;
upperlimit = b = 1;

fa = f[a];
fb = f[b];

h = b - a;

c = h * (fa + fb)/2;

Print["Integral = ", c]







#trapezoidal

f[x_] := x^2 - x - 1;
lowerlimit = a = 3;
upperlimit = b = 5;
fa = f[a];
fb = f[b];
h = b - a;
(* Trapezoidal rule *)
c = h * (fa + fb)/2;
Print["The area approximate by the trapezoidal rule is ", c];

(* Exact integral *)
int = Integrate[f[x], {x, a, b}];

Print[ "Original value of the integral is ", N[int], " and the error is ", N[c] - N[int]];



#composite
f[x_] := x + 1;
lowerlimit = a = 0;
upperlimit = b = 1;
n = 10;
fa = f[a];
fb = f[b];
h = (b - a)/n;
sum = 0;

For[i = 1, i < n, i++,sum = sum + N[f[a + i*h]]];

sum = (h/2) * (fa + 2*sum + fb);

Print["Integral = ", sum];







f[x_] := x + 1;

lowerlimit = a = 0;
upperlimit = b = 1;
fa = f[a];
fb = f[b];
h = (b - a)/2;

c = (h/3) * (fa + 4*f[a + h] + fb);

Print["Integral = ", N[c]];



#SIMPSON!1/3
f[x_] := x^2 + 2 x + 7;
lowerlimit = a = 1;
upperlimit = b = 2;
fa = f[a];
fb = f[b];
h = (b - a)/2;
c = (h/3) * (fa + 4*f[a + h] + fb);
Print["The area approximate by the Simpson's rule is ", N[c]];

(* Exact Integral *)
int = Integrate[f[x], {x, a, b}];

Print["Original value of the integral is ", N[int]," and the error is ", N[c] - N[int]];

(* Trapezoidal Rule *)
byTrap = N[(b - a)/2 * (fa + fb)];

Print["The area approximate by the trapezoidal rule is ", N[byTrap]];











#3/8
f[x_] := Sin[x]/(x^2 + 1);

lowlimit = a = 1;
upplimit = b = 3;

fa = f[a];
fb = f[b];

h = (b - a)/3;

(* Simpson's 3/8 Rule *)
c = (3*h/8) * (fa + 3*f[a + h] + 3*f[a + 2*h] + fb);

Print["The area approximated by Simpson's rule is ", N[c]];




#composite1/3
f[x_] := Exp[x];
lowlimit = a = 0;
upplimit = b = 4;
n = 20;
fa = f[a];
fb = f[b];
h = (b - a)/n;
sumodd = 0;
sumeven = 0;
(* odd terms *)
For[i = 1, i < n, i += 2, sumodd = sumodd + N[f[a + i*h]]];

(* even terms *)
For[i = 2, i < n, i += 2,sumeven = sumeven + N[f[a + i*h]]];

(* Simpson 1/3 composite formula *)
sum = (h/3) * (fa + 4*sumodd + 2*sumeven + fb);

Print["Integral = ", sum];

(* exact integral *)
exact = N[Integrate[f[x], {x, a, b}]];
Print["Exact value of integral is ", exact];

(* error *)
error = N[sum] - N[exact];

Print["The error between exact and numerical value is: ", error];






#euler
f[x_, y_] := 1 + (y/x);

a = 1;
b = 6;
n = 10;
h = (b - a)/n;
y[1] = 1;
temp = y[1];

For[i = 0, i <= n - 1, i++,
  x[i] = a + i*h;
  y[i] = temp; 
  y[i + 1] = y[i] + h*f[x[i], y[i]];
  
  Print["The ", i + 1, " approximation is ", N[y[i + 1]]];
  
  temp = y[i + 1];
];
