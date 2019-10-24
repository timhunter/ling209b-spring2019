{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

module ForwardBackward where

import qualified Data.Map as M
import qualified Data.List as List

import PFSA(PFSA, State, Symbol, Prob, forward, backward, stringprobViaBackward)

-- Here's the PFSA in (3.23) in the class handouts.
m1 :: PFSA
m1 = (  [1,2], 
        ["a","b"], 
        \q -> M.findWithDefault 0 q (M.fromList [(1,1.0)]), 
        \q -> M.findWithDefault 0 q (M.fromList [(1,0.250)]), 
        \t -> M.findWithDefault 0 t (M.fromList [((1,"b",1), 0.250), 
                                                 ((1,"a",2), 0.250), 
                                                 ((1,"b",2), 0.250), 
                                                 ((2,"b",2), 0.333), 
                                                 ((2,"a",1), 0.333), 
                                                 ((2,"b",1), 0.333)
                                                ])
     )

---------------------------------------------------------
-- SoftCounts is a lookup table storing a collection of 
-- (soft) counts of the three sorts listed in (3.14).

type SoftCounts = M.Map HiddenEvent Prob

data HiddenEvent = HStart State | HEnd State | HStep State Symbol State deriving (Show,Eq,Ord)

-- Adds up the entries across a collection of 
-- these lookup table.
sumCounts :: [SoftCounts] -> SoftCounts
sumCounts xs =
    let zeroCounts = M.fromList [] in
    let addCounts m1 m2 = M.unionWith (+) m1 m2 in
    foldl addCounts zeroCounts xs

---------------------------------------------------------
-- Observing a string provides a collection of ``observed events''. 
-- For example, observing the string "abc" amounts to observing 
--  - a starting event for "abc"
--  - an ending event for "abc"
--  - a transition event for "a", after "" and before "bc"
--  - a transition event for "b", after "a" and before "c"
--  - a transition event for "c", after "ab" and before ""
-- The total number of such observed events in a string (or corpus) 
-- will equal the total of the three kinds of soft counts that we calculate.

data ObservedEvent = Start [Symbol] | End [Symbol] | Mid [Symbol] Symbol [Symbol] deriving Show

observedEvents :: [Symbol] -> [ObservedEvent]
observedEvents str = [Start str, End str] ++ [Mid (take i str) (str !! i) (drop (i+1) str) | i <- [0 .. (length str - 1)]]

---------------------------------------------------------
-- This next function, expectationsFromObs, is the core of the E-step. 
-- It calculates soft counts from a single observed event. 
--  - If the observed event is a starting event, then this function constructs a lookup table 
-- mapping each state q to the corresponding conditional probability in equation (3.18).
--  - If the observed event is an ending event, then this function constructs a lookup table 
-- mapping each state q to the corresponding conditional probability in equation (3.19).
--  - If the observed event is a transition event, then this function constructs a lookup table 
-- mapping each pair of states q,q' to the corresponding conditional probability in equation (3.22).
-- I take advantage of the fact that each resulting lookup table should have probabilities summing 
-- to one, and use a normalizing function rather than dividing by the probability of the observed string.

expectationsFromObs :: PFSA -> ObservedEvent -> SoftCounts
expectationsFromObs m e =
    let (states, sigma, initprob, finprob, trprob) = m in
    case e of
    Start str ->
        let likelihoods = M.fromList [(HStart q, initprob q * backward m str q) | q <- states] in
        normalize likelihoods
    End str ->
        let likelihoods = M.fromList [(HEnd q, forward m str q * finprob q) | q <- states] in
        normalize likelihoods
    Mid before x after ->
        let likelihoods = M.fromList [(HStep q1 x q2, forward m before q1 * trprob (q1,x,q2) * backward m after q2) | q1 <- states, q2 <- states] in
        normalize likelihoods

-- This just tallies up all the expected/soft counts across all observed events 
-- from all string in the corpus.
expectationsFromCorpus :: PFSA -> [[Symbol]] -> SoftCounts
expectationsFromCorpus fsa strings =
    sumCounts [expectationsFromObs fsa o | s <- strings, o <- observedEvents s]

---------------------------------------------------------
-- Now the M step, pretty simple

-- Given some (soft) counts, do the simple relative frequency estimation.
estimate_from_counts :: ([State],[Symbol]) -> SoftCounts -> PFSA
estimate_from_counts (states,sigma) counts =
    let startDist = normalize (M.fromList [(q, c) | (HStart q, c) <- M.assocs counts]) in
    let stepDist q = normalize (M.fromList $ [(Just (x,q2), c) | (HStep q1 x q2, c) <- M.assocs counts, q1 == q] ++ 
                                             [(Nothing, c)     | (HEnd q1, c) <- M.assocs counts, q1 == q])
    in
    let newInitFn = \q -> M.findWithDefault 0 q startDist in
    let newFinFn = \q -> M.findWithDefault 0 Nothing (stepDist q) in
    let newTrFn = \(q1,x,q2) -> M.findWithDefault 0 (Just (x,q2)) (stepDist q1) in
    (states, sigma, newInitFn, newFinFn, newTrFn)

---------------------------------------------------------
-- Put it all together

-- Given a current PFSA and a list of strings, produce a new PFSA after 
-- one round of expectation-maximization.
-- Example:
--  *ForwardBackward> update m1 [["b","a","a","b","a","a","a","a","b"]]
--  {I(1) = 1.0, F(1) = 0.1707174231332357, D(1,b,1) = 0.17780501219270242, D(1,a,2) = 0.5121522693997071, D(1,b,2) = 0.13932529527435467, D(2,a,1) = 0.7242236024844722, D(2,b,1) = 0.19701692892833175, D(2,b,2) = 7.875946858719625e-2}
update :: PFSA -> [[Symbol]] -> PFSA
update m strings =
    let (states, sigma, initprob, finprob, trprob) = m in
    let expected_counts = expectationsFromCorpus m strings in
    estimate_from_counts (states,sigma) expected_counts

-- Same as the plain update function but reports new likelihood as well
update' :: PFSA -> [[Symbol]] -> (PFSA, Prob)
update' fsa strings =
    let new_fsa = update fsa strings in
    let new_likelihood = likelihood new_fsa strings in
    (new_fsa, new_likelihood)

likelihood :: PFSA -> [[Symbol]] -> Prob
likelihood m strings = product (map (stringprobViaBackward m) strings)

---------------------------------------------------------
-- Utility functions, feel free to ignore

normalize :: M.Map a Prob -> M.Map a Prob
normalize m =
    let total = sum (M.elems m) in
    M.map (/total) m

printMap :: (Show a, Show b) => M.Map a b -> IO ()
printMap = putStr . unlines . showMap
    where
        showMap m = [show k ++ "\t" ++ show v | (k,v) <- M.assocs m]

showPFSA :: PFSA -> String
showPFSA (states, sigma, initprob, finprob, trprob) =
    let initStrs = ["I(" ++ show q ++ ") = " ++ show (initprob q) | q <- states, initprob q /= 0] in
    let finStrs = ["F(" ++ show q ++ ") = " ++ show (finprob q) | q <- states, finprob q /= 0] in
    let trStrs = ["D(" ++ List.intercalate "," [show q1, x, show q2] ++ ") = " ++ show (trprob (q1,x,q2)) | q1 <- states, q2 <- states, x <- sigma, 
                                                                                                            trprob (q1,x,q2) /= 0] in
    "{" ++ List.intercalate ", " (initStrs ++ finStrs ++ trStrs) ++ "}"

instance {-# OVERLAPPING #-} Show PFSA where
    show = showPFSA

