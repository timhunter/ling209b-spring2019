
###############################################################################################

# Here's the PFSA from (2.22) on the handout. 
# I encode the functions I, F and Delta by wrapping a python dictionary in a little function that 
# looks up its argument in the dictionary and returns zero if no entry for it is found.

m0 = (  [1, 2, 3], 
        ["C","V"], 
        lambda q: {1: 1.0}.get(q,0), 
        lambda q: {1: 0.2}.get(q,0), 
        lambda *x: {(1,"V",1): 0.18, 
                    (1,"C",2): 0.5, 
                    (1,"V",3): 0.12, 
                    (2,"V",3): 0.4, 
                    (2,"V",1): 0.6, 
                    (3,"C",1): 1.0, 
                   }.get(x,0)
     )

###############################################################################################

# This function calculates forward values (given an automaton m, a string xs and a state q) 
# following equation (2.31).
# Example:
#   >>> forward(m0, ["V","C","V"], 1)
#   0.0756
def forward(m, xs, q):
    (states, sigma, initprob, finprob, trprob) = m
    if xs == []:
        return initprob(q)
    else:
        return sum([forward(m, xs[:-1], qp) * trprob(qp, xs[-1], q) for qp in states])  # read `qp' as `q prime'

# Total probability of a string, using forward values, following equation (2.27).
# Example:
#   >>> stringprob_via_forward(m0, ["V","C","V"])
#   0.015120000000000001
def stringprob_via_forward(m, xs):
    (states, sigma, initprob, finprob, trprob) = m
    return sum([forward(m, xs, q) * finprob(q) for q in states])

###############################################################################################

# This function calculates backward values (given an automaton m, a string xs and a state q) 
# following equation (2.38).
# Example:
#   >>> backward(m0, ["C","V"], 3)
#   0.036
def backward(m, xs, q):
    (states, sigma, initprob, finprob, trprob) = m
    if xs == []:
        return finprob(q)
    else:
        return sum([trprob(q, xs[0], qp) * backward(m, xs[1:], qp) for qp in states])   # read `qp' as `q prime'

# Total probability of a string, using backward values, following equation (2.34).
# Example:
#   >>> stringprob_via_backward(m0, ["V","C","V"])
#   0.015119999999999998
def stringprob_via_backward(m, xs):
    (states, sigma, initprob, finprob, trprob) = m
    return sum([initprob(q) * backward(m, xs, q) for q in states])

