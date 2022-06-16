/*
Test whether a list of exponents is 0-APN in dimension N.
Uses auxillary files "exponents" and "wall_representatives",
where wall_representatives should be of the form

W := [...];

with outputs from generate_wall_reps.m
*/

function coordinates2element(number,dim)
	coords := Intseq(number,2);
	coords := Reverse(coords);
	if #coords lt dim then
		coords := [ 0 : i in [1..dim - #coords] ] cat coords;
	end if;
	a := GF(2^dim).1;
	return &+ [ a^(i-1) * coords[i] : i in [1..dim] ];
end function;

N := 35;
load "exponents";
load "wall_representatives";
Walls := [ coordinates2element(i,N) : i in W];
delete W;

t := Cputime();
for e in exp do
	for w in Walls do
		if w^e + (w+1)^e eq 1 then
			printf "\n%o is not 0-APN\n\n", e;
			break;
		end if;
	end for;
end for;
printf "Running time : %o seconds", Cputime(t);
