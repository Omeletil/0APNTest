/* 
Generate wall representatives 
Where "cosets" should be an auxillary defining the Magma sequence
s := [...];
with outputs from cosets.c
*/

load "util.m";
//n := 35;
ff<a> := GF(2^n);

load "cosets";
Results := [];
for i in s do
	c := a^i;
	Wall := &join {{x^(2^k) : k in [0..n-1]} : x in [c, 1+c, c^(-1), 1 + c^(-1), (1 + c)^(-1), 1 + (1 + c)^(-1)]};
	representative := Min(Wall);
	integer_representation := &+ [ 2^(i-1) * (Integers() ! B[i]) : i in [1..#B] ] where B is Reverse(Flat(representative));
	Append(~Results, integer_representation);
end for;

Results;
