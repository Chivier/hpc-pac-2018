# Computing Part

`run_SOBI.m` is the main function. It consists of two parts, the first part is a Matrix-reader to read data in `DATA.mat`. And then it call `SOBI.m` to calculate all we want.

So our main mission is to convert `SOBI.m` into C++ version and then optimize it.

# Verification Part

All we need to do is to run `result_testing.m` and then we are able to verify our program.

# Understanding SOBI

SOBI is abbreviation for 'Second=Order Blind Identification'.

The story comes. We are given a time-related function `X(t)`. And we wanna get the source function `S(t)`, which satisfies the following equation.

$$
X(t) = A \times S(t)
$$

A is a signal matrix. And there is we are given this matrix's inverse matrix `T`.

# Analysis on `SOBI.m`

