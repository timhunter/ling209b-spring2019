module PFSA where

import Data.Map as M

type State = Int
type Symbol = String
type Prob = Double
type PFSA = ([State], [Symbol], State -> Prob, State -> Prob, (State,Symbol,State) -> Prob)

--------------------------------------------------------------------------

-- Here's the PFSA from (2.22) on the handout.
m0 :: PFSA
m0 = (  [1,2,3], 
        ["C","V"], 
        \q -> M.findWithDefault 0 q (M.fromList [(1,1.0)]), 
        \q -> M.findWithDefault 0 q (M.fromList [(1,0.2)]), 
        \t -> M.findWithDefault 0 t (M.fromList [((1,"V",1), 0.18), 
                                                 ((1,"C",2), 0.5), 
                                                 ((1,"V",3), 0.12), 
                                                 ((2,"V",3), 0.4), 
                                                 ((2,"V",1), 0.6), 
                                                 ((3,"C",1), 1.0)
                                                ])
     )

--------------------------------------------------------------------------

-- This function calculates forward values following equation (2.31).
-- NB:   init [a,b,c,d,e,f] == [a,b,c,d,e]    last [a,b,c,d,e,f] == f
-- Example:
--  *PFSA> forward m0 ["V","C","V"] 1
--  7.56e-2
forward :: PFSA -> [Symbol] -> State -> Prob
forward m xs q =
    let (states, sigma, initprob, finprob, trprob) = m in
    if xs == [] then
        initprob q
    else
        sum [forward m (init xs) q' * trprob (q',(last xs),q) | q' <- states]

-- Total probability of a string, using forward values, following equation (2.27).
-- Example:
--  *PFSA> stringprobViaForward m0 ["V","C","V"]
--  1.5120000000000001e-2
stringprobViaForward :: PFSA -> [Symbol] -> Prob
stringprobViaForward m xs =
    let (states, sigma, initprob, finprob, trprob) = m in
    sum [forward m xs q * finprob q | q <- states]

--------------------------------------------------------------------------

-- This function calculcates backward values following equation (2.38).
-- NB:   head [a,b,c,d,e,f] == a     tail [a,b,c,d,e,f] == [b,c,d,e,f]
-- Example:
--  *PFSA> backward m0 ["C","V"] 3
--  3.6e-2
backward :: PFSA -> [Symbol] -> State -> Prob
backward m xs q =
    let (states, sigma, initprob, finprob, trprob) = m in
    if xs == [] then
        finprob q
    else
        sum [trprob (q,(head xs),q') * backward m (tail xs) q' | q' <- states]

-- Total probability of a string, using backward values, following equation (2.34).
-- Example:
--  *PFSA> stringprobViaBackward m0 ["V","C","V"]
--  1.5119999999999998e-2
stringprobViaBackward :: PFSA -> [Symbol] -> Prob
stringprobViaBackward m xs =
    let (states, sigma, initprob, finprob, trprob) = m in
    sum [initprob q * backward m xs q | q <- states]

