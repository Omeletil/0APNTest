/* Generate wall representatives given a list of elements from a finite field. */
function getWallRepresentatives(ELS)
	ELS := {x : x in ELS};
	Walls := [];
	while #ELS gt 0 do
		c := Random(ELS);
		FF := Parent(c);
		n := Degree(FF);
		Wall := &join {{x^(2^k) : k in [0..n-1]} : x in [c, 1+c, c^(-1), 1 + c^(-1), (1 + c)^(-1), 1 + (1 + c)^(-1)]};
		Append(~Walls, Minimum(Wall));
		ELS := ELS diff Wall;
	end while;
	return Walls;
end function;

function getSubfields(FF)
	n := Degree(FF);
	subfields := &join {{b : b in GF(2^k)} : k in Divisors(n) | k ne n};
	return subfields;
end function;

function isAPN(d, FF)
	n := Degree(FF);
	pr<z> := PolynomialRing(FF);
	derivative_image := {Evaluate(z^d, b) + Evaluate(z^d, b + 1) : b in FF};
	return (#derivative_image ne 2^(n - 1)) select false else true;
end function;


/* Generate the elements that satisfy x^(2^k + 1) + 1 = 0 in GF(2^n) */
function generateAB0(ff : ks := [])
	n := Degree(ff);
	a := ff.1;
	if IsEmpty(ks) then
		ks := {k : k in [1..n-1] | n mod (2*k) eq 0};
		ks := Setseq(ks);
	end if;
	phi := [&+ [EulerPhi(c) : c in Divisors(2^k + 1) | (2^n-1) mod c eq 0] : k in ks];
	exponents := [(2^n - 1) div p : p in phi];
	AB0 := {};
	for i in [1..#ks] do
		k := ks[i];
		e := exponents[i];
		generator := a^e;
		AB0 := AB0 join {generator^j : j in [1..phi[i]]};
	end for;
	return AB0;
end function;

/* Generate the elements that satisfy x^(2^k + 1) + x^(2^k) + x = 0 in GF(2^n) */
function generateBB1(ff)
	n := Degree(ff);
	Ks := Setseq({GCD(k, n) : k in [1..n-1]});
	BB1 := generateAB0(ff : ks := Ks);
	return {b + 1 : b in BB1};
end function;

/* Generate the elements that satisfy x^(2^k + 1) + x^(2^k) + 1 = 0 in GF(2^n) */
function generateAB1(ff)
	n := Degree(ff);
	a := ff.1;
	// If 3 does not divide n, then the Bergen set is trivial
	if (n mod 3 ne 0) then
		return {};
	end if;
	AB1 := {};
	k := n div 3;
	generator := a^(2^k - 1);
	i := 1;
	element := generator^i;
	while element ne 1 do
		if Trace( a^i, GF(2^k) ) eq 0 then
			AB1 := AB1 join { element, element+1 };
		end if;
		element := element * generator;
		i +:= 1;
	end while;
	return AB1;
end function;

/* Generate the non-trivial elements that satisfy x^(2^k) + x + 1 = 0 in GF(2^n) */
function generateAA1(ff)
	n := Degree(ff);
	a := ff.1;
	if n mod 2 eq 1 then
		return {};
	end if;
	k := n div 2;
	generator := a^(2^k + 1);
	half_field := {generator^i : i in [1..2^k - 1]} join {0};
	// Find element with Trace(x, GF(2^k)) = 1
	elm := ff ! 1;
	while Trace(elm, GF(2^k)) ne 1 do	
		elm := elm * a;
	end while;
	AA1 := {b + elm : b in half_field};
	return AA1;
end function;

function getGBS(FF)
	A := generateAB0(FF);
	B := generateAB1(FF);
	C := generateAA1(FF);
	D := generateBB1(FF);
	GBS := &join {A, B, C, D};
	return GBS diff {0, 1};
end function;
