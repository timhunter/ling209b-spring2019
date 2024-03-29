module PCFG where

import Data.Map as M

type Nonterminal = String
type Terminal = String
type Prob = Double
data Rule = BranchRule Nonterminal Nonterminal Nonterminal | EndRule Nonterminal Terminal deriving (Eq,Ord,Show)
type PCFG = ([Nonterminal], [Terminal], Nonterminal -> Prob, Rule -> Prob)

-- Here's the PCFG from (2.57) on the handout. 
g0 :: PCFG
g0 = (  ["VP","NP","PP","V","P"], 
        ["watches","spies","telescopes","with"], 
        \n -> M.findWithDefault 0 n (M.fromList [("VP",1.0)]), 
        \r -> M.findWithDefault 0 r (M.fromList [(BranchRule "VP" "V" "NP",     0.4), 
                                                 (BranchRule "VP" "VP" "PP",    0.2), 
                                                 (EndRule    "VP" "watches",    0.3), 
                                                 (EndRule    "VP" "spies",      0.1), 
                                                 (BranchRule "NP" "NP" "PP",    0.2), 
                                                 (EndRule    "NP" "watches",    0.3), 
                                                 (EndRule    "NP" "spies",      0.2), 
                                                 (EndRule    "NP" "telescopes", 0.3), 
                                                 (BranchRule "PP" "P" "NP",     1.0), 
                                                 (EndRule    "V"  "watches",    1.0), 
                                                 (EndRule    "P"  "with",       1.0)
                                                ])
     )

--------------------------------------------------------------------------

-- This function calculates inside values following equation (2.69).
-- It works, but gets very slow if you try strings longer than about four or five symbols, because 
-- previously-computed values don't get re-used in the way they do if we work with a chart like in (2.68).
-- Example:
--  *PCFG> inside g0 ["watches","spies","with","telescopes"] "VP"
--  9.600000000000003e-3
inside :: PCFG -> [Terminal] -> Nonterminal -> Prob
inside g xs n =
    let (nonterms, sigma, initprob, ruleprob) = g in
    case xs of
    [] -> 0
    [x] -> ruleprob (EndRule n x)
    _ -> sum [ruleprob (BranchRule n l r) * inside g (take i xs) l * inside g (drop i xs) r | 
                    i <- [1 .. length xs - 1], 
                    l <- nonterms, 
                    r <- nonterms
             ]

-- This function calculates outside values following equation (2.74).
-- Example:
--  *PCFG> outside g0 (["watches"],["with","telescopes"]) "NP"
--  4.8e-2
--  *PCFG> outside g0 (["watches","spies"],[]) "PP"
--  3.200000000000001e-2
outside :: PCFG -> ([Terminal],[Terminal]) -> Nonterminal -> Prob
outside g (ys,zs) n =
    let (nonterms, sigma, initprob, ruleprob) = g in
    if ys == [] && zs == [] then
        initprob n
    else
        sum [ruleprob (BranchRule p n r) * inside g (take i zs) r * outside g (ys, drop i zs) p | 
                    i <- [1 .. length zs], 
                    p <- nonterms, 
                    r <- nonterms
            ]
        +
        sum [ruleprob (BranchRule p l n) * inside g (drop i ys) l * outside g (take i ys, zs) p | 
                    i <- [0 .. length ys - 1], 
                    p <- nonterms, 
                    l <- nonterms
            ]

--------------------------------------------------------------------------

-- This function constructs a chart like the one shown in (2.68), given a PCFG and a 
-- list of terminal symbols. 
-- The result is a Map that maps (symbol-list, nonterminal) pairs to probabilities. 
-- Example:
--  *PCFG> let c = insideChart g0 ["watches","spies","with","telescopes"] in (c ! (["watches","spies","with","telescopes"],"VP"))
--  9.600000000000003e-3
--  *PCFG> let c = insideChart g0 ["watches","spies","with","telescopes"] in (c ! (["spies","with","telescopes"],"NP"))
--  1.2000000000000002e-2

insideChart :: PCFG -> [Terminal] -> M.Map ([Terminal],Nonterminal) Prob
insideChart g string =

    let (nonterms, sigma, initprob, ruleprob) = g in

    -- Set up the list of cells to be filled, in order; 
    -- each cell is represented by a substring
    let cellsToFill = [take len (drop startpos string) | 
                            len <- [1 .. length string], 
                            startpos <- [0 .. length string - len]
                      ]
    in

    -- Define a function that will add the inside value for a single (string,nonterminal) pair 
    -- to an existing chart
    let oneInsideValue xs n chart =
            let prob = case xs of
                       [] -> 0
                       [x] -> ruleprob (EndRule n x)
                       _ -> sum [ruleprob (BranchRule n l r) * (M.findWithDefault 0 (take i xs, l) chart) * (M.findWithDefault 0 (drop i xs, r) chart) |
                                    i <- [1 .. length xs - 1], 
                                    l <- nonterms, 
                                    r <- nonterms
                                ]
            in
            if prob /= 0 then M.insert (xs,n) prob chart else chart
    in

    -- Define a function which takes a list of cells and a chart, and 
    -- fills in those cells in the given order
    let fillCells cells chart =
            case cells of
            [] -> chart
            xs:rest -> let newChart = Prelude.foldl (\ch -> \n -> oneInsideValue xs n ch) chart nonterms in
                       fillCells rest newChart
    in

    -- Now, fill in all the cells, starting with an empty chart
    fillCells cellsToFill M.empty

-- Total probability of a string, using inside values from a chart, following equation (2.66).
-- Example:
--  *PCFG> stringProbViaInsideChart g0 ["watches","spies","with","telescopes"]
--  9.600000000000003e-3
stringProbViaInsideChart :: PCFG -> [Terminal] -> Prob
stringProbViaInsideChart g xs =
    let (nonterms, sigma, initprob, ruleprob) = g in
    let chart = insideChart g xs in
    sum [initprob n * M.findWithDefault 0 (xs,n) chart | n <- nonterms]

--------------------------------------------------------------------------

-- This function constructs a chart of outside values, given a PCFG and a sequence of 
-- terminal symbols. We didn't discuss this in class, but it's analogous to constructing 
-- a chart of inside values, just in reverse; see p.353 and Figure 3.2 of Jelinek et al (1992). 
-- The result is a Map that maps ((symbol-list,symbol-list), nonterminal) pairs to probabilities. 
-- Example:
--  *PCFG> let c = outsideChart g0 ["watches","spies","with","telescopes"] in (c ! ((["watches"],["with","telescopes"]),"NP"))
--  4.8e-2
--  *PCFG> let c = outsideChart g0 ["watches","spies","with","telescopes"] in (c ! ((["watches","spies"],[]),"PP"))
--  3.200000000000001e-2

outsideChart :: PCFG -> [Terminal] -> M.Map (([Terminal],[Terminal]),Nonterminal) Prob
outsideChart g string =

    let (nonterms, sigma, initprob, ruleprob) = g in

    -- Get a completed chart of inside values
    let ichart = insideChart g string in

    -- Set up the list of cells to be filled, in order; 
    -- each cell is represented by a pair of strings (ys,zs)
    let cellsToFill = [((take startpos string, drop (startpos+len) string), n) |
                            len <- [length string, length string - 1 .. 1], 
                            startpos <- [0 .. length string - len], 
                            n <- nonterms
                      ]
    in

    -- Define a function that will add the outside value for a single ((string,string),nonterminal) pair 
    -- to an existing chart
    let oneOutsideValue (ys,zs) n chart =
            let prob =
                    if ys == [] && zs == [] then
                        initprob n
                    else
                        sum [ruleprob (BranchRule p n r) * (M.findWithDefault 0 (take i zs, r) ichart) * (M.findWithDefault 0 ((ys, drop i zs), p) chart) |
                                    i <- [1 .. length zs], 
                                    p <- nonterms, 
                                    r <- nonterms
                            ]
                        +
                        sum [ruleprob (BranchRule p l n) * (M.findWithDefault 0 (drop i ys, l) ichart) * (M.findWithDefault 0 ((take i ys, zs), p) chart) |
                                    i <- [0 .. length ys - 1], 
                                    p <- nonterms, 
                                    l <- nonterms
                            ]
            in
            if prob /= 0 then M.insert ((ys,zs),n) prob chart else chart
    in

    -- Define a function which takes a list of cells and a chart, and 
    -- fills in those cells in the given order
    let fillCells cells chart =
            case cells of
            [] -> chart
            ((ys,zs),n):rest -> let newChart = Prelude.foldl (\ch -> \n -> oneOutsideValue (ys,zs) n ch) chart nonterms in
                                fillCells rest newChart
    in

    -- Now, fill in all the cells, starting with an empty chart
    fillCells cellsToFill M.empty

