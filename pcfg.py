
from collections import defaultdict

###############################################################################################

# Here's the PCFG from (2.57) on the handout. 
# I encode the functions I and R by wrapping a python dictionary in a little function that 
# looks up its argument in the dictionary and returns zero if no entry for it is found.

g0 = (  ["VP","NP","PP","V","P"] , 
        ["watches","spies","telescopes","with"] , 
        lambda nt: {"VP": 1.0}.get(nt,0) ,
        lambda *x: {("VP", ("V","NP")):   0.4, 
                    ("VP", ("VP","PP")):  0.2, 
                    ("VP", "watches"):    0.3, 
                    ("VP", "spies"):      0.1, 
                    ("NP", ("NP","PP")):  0.2, 
                    ("NP", "watches"):    0.3, 
                    ("NP", "spies"):      0.2, 
                    ("NP", "telescopes"): 0.3, 
                    ("PP", ("P","NP")):   1.0, 
                    ("V",  "watches"):    1.0, 
                    ("P",  "with"):       1.0, 
                   }.get(x,0)
     )

###############################################################################################

# This function calculates inside values (given a grammar g, a string xs and a nonterminal symbol n) 
# following equation (2.69) in a completely direct way.
# It works, but gets very slow if you try strings longer than about four or five symbols, because 
# previously-computed values don't get re-used the way they do if we work with a chart like in (2.68).
# Example:
#   >>> inside(g0, ["watches","spies","with","telescopes"], "VP")
#   0.009600000000000003
def inside(g, xs, n):
    (nonterms, sigma, initprob, ruleprob) = g
    if len(xs) == 1:
        return ruleprob(n, xs[0])
    else:
        return sum([ruleprob(n,(l,r)) * inside(g, xs[:i], l) * inside(g, xs[i:], r) 
                        for i in range(1,len(xs)) 
                        for l in nonterms 
                        for r in nonterms
                   ])

# This function calculates outside values (given a grammar g, a pair of strings 
# and a nonterminal symbol n) following equation (2.74) directly.
# Example:
#   >>> outside(g0, (["watches"],["with","telescopes"]), "NP")
#   0.048
#   >>> outside(g0, (["watches","spies"],[]), "PP")
#   0.03200000000000001
def outside(g, stringpair, n):
    (ys,zs) = stringpair
    (nonterms, sigma, initprob, ruleprob) = g
    if ys == [] and zs == []:
        return initprob(n)
    else:
        return (
            sum([ruleprob(p,(n,r)) * inside(g, zs[:i], r) * outside(g, (ys,zs[i:]), p) 
                    for i in range(1,len(zs)+1) 
                    for p in nonterms 
                    for r in nonterms
                ])
            +
            sum([ruleprob(p,(l,n)) * inside(g, ys[i:], l) * outside(g, (ys[:i],zs), p) 
                    for i in range(0,len(ys)) 
                    for p in nonterms 
                    for l in nonterms
                ])
        )

###############################################################################################

# This function constructs a chart like the one shown in (2.68), given a PCFG and a 
# sequence of terminal symbols. 
# The result is a dictionary that maps (symbol-sequence, nonterminal) pairs to probabilities. 
# The symbol-sequence when doing a lookup has to be a tuple, not a list, because python's lists are mutable and not hashable.
# Example:
#   >>> c = inside_chart(g0, ["watches","spies","with","telescopes"])
#   >>> c[(("watches","spies","with","telescopes"),"VP")]
#   0.009600000000000003
#   >>> c[(("spies","with","telescopes"),"NP")]
#   0.012000000000000002

def inside_chart(g, string):

    (nonterms, sigma, initprob, ruleprob) = g
    chart = {}

    # Each diagonal of the CKY chart represents substrings of a certain length.
    # We loop over these diagonals by looping over lengths from 1 up to the length of the entire string.
    for length in range(1, len(string)+1):

        # Each cell on this diagonal represents a particular substring of the chosen length.
        # We loop over these substrings by looping over starting positions for substrings, i.e. by looping 
        # over rows of the chart (but limited to those rows that contain a cell on this diagonal).
        for startpos in range(0, len(string)-length+1):

            # This is the substring/infix represented by this cell
            xs = tuple(string[startpos:startpos+length])

            for n in nonterms:

                # Calculate inside(xs)(n) using equation (2.69)
                if length == 1:
                    p = ruleprob(n, xs[0])
                else:
                    p = sum([ruleprob(n,(l,r)) * chart[(xs[:i],l)] * chart[(xs[i:],r)] 
                                for i in range(1,length) 
                                for l in nonterms 
                                for r in nonterms
                            ])

                # Store the result in the chart
                chart[(xs,n)] = p

    # Just for convenience: make all the zero entries in the chart implicit before returning it
    tidy_chart = defaultdict(float, [(key,p) for (key,p) in chart.items() if p != 0])
    return tidy_chart

# Total probability of a string, using inside values from a chart, following equation (2.66).
# Example:
#   >>> stringprob_via_inside_chart(g0, ["watches","spies","with","telescopes"])
#   0.009600000000000003
def stringprob_via_inside_chart(g, xs):
    (nonterms, sigma, initprob, ruleprob) = g
    chart = inside_chart(g, xs)
    return sum([initprob(n) * chart[(tuple(xs),n)] for n in nonterms])

###############################################################################################

# This function constructs a chart of outside values, given a PCFG and a sequence of 
# terminal symbols. We didn't discuss this in class, but it's analogous to constructing 
# a chart of inside values, just in reverse; see p.353 and Figure 3.2 of Jelinek et al (1992). 
# The result is a dictionary that maps ((symbol-sequence,symbol-sequence), nonterminal) pairs 
# to probabilities. 
# Example:
#   >>> oc = outside_chart(g0, ["watches","spies","with","telescopes"])
#   >>> oc[((("watches",),("with","telescopes")),"NP")]
#   0.048
#   >>> oc[((("watches","spies"),()),"PP")]
#   0.03200000000000001

def outside_chart(g, string):

    (nonterms, sigma, initprob, ruleprob) = g
    ochart = {}

    # First get a complete chart of all relevant inside values, which we will need
    ichart = inside_chart(g, string)

    # For each diagonal of the chart, now working down and left from the top-right corner
    for length in range(len(string), 0, -1):

        # For each row of the chart for which there is a cell on this diagonal
        for startpos in range(0, len(string)-length+1):

            ys = tuple(string[:startpos])
            zs = tuple(string[startpos+length:])

            for n in nonterms:

                # Calculate outside(ys,zs)(n) using equation (2.74)
                if len(ys) == 0 and len(zs) == 0:
                    prob = initprob(n)
                else:
                    prob = (sum([ruleprob(p,(n,r)) * ichart[(zs[:i],r)] * ochart[((ys,zs[i:]),p)]
                                    for i in range(1,len(zs)+1) 
                                    for p in nonterms 
                                    for r in nonterms
                                ])
                            +
                            sum([ruleprob(p,(l,n)) * ichart[(ys[i:],l)] * ochart[((ys[:i],zs),p)] 
                                    for i in range(0,len(ys)) 
                                    for p in nonterms 
                                    for l in nonterms
                                ])
                           )

                # Store the result in the chart
                ochart[((ys,zs),n)] = prob

    # Just for convenience: make all the zero entries in the chart implicit before returning it
    tidy_ochart = defaultdict(float, [(key,p) for (key,p) in ochart.items() if p != 0])
    return tidy_ochart

