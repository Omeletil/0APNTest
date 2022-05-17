function GoldExponents(n)
	return { 2^i + 1 : i in [1..n] | GCD(i,n) eq 1 };
end function;

function KasamiExponents(n)
	limit := Round((n mod 2 eq 0) select (n-2)/2 else (n-1)/2);
	return { 2^(2*i) - 2^i + 1 : i in [1..limit] | GCD(i,n) eq 1 };
end function;

function WelchExponents(n)
	if n mod 2 eq 0 then
		return {};
	else
		t := (n-1) div 2;
		return {2^t + 3};
	end if;
end function;

function NihoExponents(n)
	if n mod 2 eq 0 then
		return {};
	else
		t := (n-1) div 2;
		if t mod 2 eq 0 then
			return {2^t + 2^(t div 2) - 1 };
		else
			return {2^t + 2^( (3*t + 1) div 2 ) -1 };
		end if;
	end if;
end function;

function InverseExponents(n)
	if n mod 2 eq 0 then
		return {};
	else
		t := (n-1) div 2;
		return {2^(2*t) - 1 };
	end if;
end function;

function DobbertinExponents(n)
	if n mod 5 ne 0 then
		return {};
	else
		i := n div 5;
		return {2^(4*i) + 2^(3*i) + 2^(2*i) + 2^i - 1};
	end if;
end function;

function AllExponents(n)
	return &join [GoldExponents(n), KasamiExponents(n), WelchExponents(n), NihoExponents(n), InverseExponents(n), DobbertinExponents(n)];
end function;

function CyclotomicCoset(k,n)
	coset := { (k*2^i) mod (2^n-1) : i in [0..n-1] };
	return coset;
end function;

function KnownExponents(n)
	known_exponents := AllExponents(n);
	known_exponents := &join {CyclotomicCoset(x, n) : x in known_exponents};
	known_exponents := known_exponents join (&join {CyclotomicCoset(InverseMod(x, 2^n - 1), n) : x in known_exponents | GCD(x, 2^n - 1) eq 1}); 
	return known_exponents;
end function;

function GenerateExponents(n)
	if (IsEven(n)) then
		exp := {a : a in [3..2^n -2] | GCD(a, 2^n - 1) eq 3};
	else
		exp := {a : a in [3..2^n - 2] | GCD(a, 2^n - 1) eq 1}; 
	end if;
	repr := {};
	while not IsEmpty(exp) do
		a := Random(exp);
		coset := CyclotomicCoset(a, n);
		a := Minimum(coset);
		min_representative := a;
		if GCD(a, 2^n - 1) eq 1 then
			inv_coset := CyclotomicCoset(InverseMod(a, 2^n - 1), n);
			exp := exp diff inv_coset;
			min_representative := Minimum(a, Minimum(inv_coset));
		end if;
		exp := exp diff coset;
		Include(~repr, min_representative);
	end while;
	known_exponents := AllExponents(n);
	known_exponents := &join {CyclotomicCoset(x, n) : x in known_exponents};
	known_exponents := known_exponents join (&join {CyclotomicCoset(InverseMod(x, 2^n - 1), n) : x in known_exponents | GCD(x, 2^n - 1) eq 1});
	repr := repr diff known_exponents;
	return repr;
end function;

function GenerateAllExponents(n)
	if (IsEven(n)) then
		exp := {a : a in [3..2^n - 2] | GCD(a, 2^n - 1) eq 3};
	else
		exp := {a : a in [3..2^n - 2] | GCD(a, 2^n - 1) eq 1}; 
	end if;
	repr := {};
	while not IsEmpty(exp) do
		a := Random(exp);
		coset := CyclotomicCoset(a, n);
		a := Minimum(coset);
		min_representative := a;
		if GCD(a, 2^n - 1) eq 1 then
			inv_coset := CyclotomicCoset(InverseMod(a, 2^n - 1), n);
			exp := exp diff inv_coset;
			min_representative := Minimum(a, Minimum(inv_coset));
		end if;
		exp := exp diff coset;
		Include(~repr, min_representative);
	end while;
	return repr;
end function;
