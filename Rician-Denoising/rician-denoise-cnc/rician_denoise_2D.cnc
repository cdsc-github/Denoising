// This file implements the CnC graph representation for Rician Denoise 2D

// Declarations
[double u <pair>];
[double f <pair>];
[double g <pair>];
[double r <pair>];
[double ug <pair>];
[double unext <pair>];
[double difference <pair>];

// Tags
<pair Iter_start>;
<pair ug_tag>;
<pair unext_tag>;
<pair converged_tag>;

// Step Prescriptions
<Iter_start> :: (compute_g);
<Iter_start> :: (compute_r);
<ug_tag> :: (compute_ug);
<unext_tag> :: (compute_unext);
<converged_tag> :: (compute_convergence);

// Input from the environment
env -> [u], [f], <Iter_start>, <unext_tag>;

// Step executions
[u] -> (compute_g) -> [g], <ug_tag>;
[u], [f] -> (compute_r) -> [r];
[u], [g] -> (compute_ug) -> [ug];
[u], [f], [g], [r], [ug] -> (compute_unext) -> [unext], <converged_tag>;
[u], [unext] -> (compute_convergence) -> [difference];

// Output to the environment
[difference], [unext] -> env;

