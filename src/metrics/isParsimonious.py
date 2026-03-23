"""
Let's say we have a tree with 3 species: A, B, C.
A and B are the ingroup.
C is the outgroup.

Let's say the sequence identity between A and B is idt_AB.
Let's say the sequence identity between A and C is idt_AC.
Let's say the sequence identity between B and C is idt_BC.

Then the sequence mismatch rate between A and B is 1 - idt_AB = x.
Then the sequence mismatch rate between A and C is 1 - idt_AC = y.
Then the sequence mismatch rate between B and C is 1 - idt_BC = z.

Let's say the probablity of mutating from the root to A is p_A.
Let's say the probablity of mutating from the root to B is p_B.
Let's say the probablity of mutating from the root to C is p_C.

We can assume:
x = p_A + p_B
y = p_A + p_C
z = p_B + p_C

Then,
p_A = (x+y-z)/2
p_B = (x+z-y)/2
p_C = (y+z-x)/2

Species A (ingroup):
The probability of Parsimony p_pars = (substitution occurred only once) is p_A/3 * (1-p_B) * (1-p_C)
The probability of non-Parsimony p_nonpars = (substitution occurred twice) is (1-p_A) * p_B/3 * p_C/3

Species B (ingroup):
The probability of Parsimony p_pars = (substitution occurred only once) is (1-p_A) * p_B/3 * (1-p_C)
The probability of non-Parsimony p_nonpars = (substitution occurred twice) is p_A/3 * (1-p_B) * p_C/3

Thus the ratio of two senarios for species A is
p_pars / p_nonpars = (p_A/3 * (1-p_B) * (1-p_C)) / ((1-p_A) * p_B/3 * p_C/3)

The ratio of two senarios for species B is
p_pars / p_nonpars = ((1-p_A) * p_B/3 * (1-p_C)) / (p_A/3 * (1-p_B) * p_C/3)

Input:
        - idt_AB
        - idt_CA
        - idt_CB
Output:
        - isParsimonious for species A
        - isParsimonious for species B
"""

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("idt_inin_AB", type=float, help="Sequence percent identity between ingroup A and ingroup B (e.g. 97)")
    parser.add_argument("idt_outin_CA", type=float, help="Sequence percent identity between outgroup C and ingroup A (e.g. 93)")
    parser.add_argument("idt_outin_CB", type=float, help="Sequence percent identity between outgroup C and ingroup B (e.g. 92)")
    args = parser.parse_args()

    x = 100 - args.idt_inin_AB
    y = 100 - args.idt_outin_CA
    z = 100 - args.idt_outin_CB

    p_A = (x + y - z) / 2
    p_B = (x + z - y) / 2
    p_C = (y + z - x) / 2

    print(f"p_A = {p_A}, p_B = {p_B}, p_C = {p_C}")
    print(f"p_pars / p_nonpars for species A = (p_A/3 * (100-p_B) * (100-p_C)) / ((100-p_A) * p_B/3 * p_C/3) = {((p_A/3 * (100-p_B) * (100-p_C)) / ((100-p_A) * p_B/3 * p_C/3))}")
    print(f"p_pars / p_nonpars for species B = ((100-p_A) * p_B/3 * (100-p_C)) / (p_A/3 * (100-p_B) * p_C/3) = {((100-p_A) * p_B/3 * (100-p_C)) / (p_A/3 * (100-p_B) * p_C/3)}")
