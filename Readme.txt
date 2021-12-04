====== README =====
Matlab replication code for:

Author(s): Q. Wang and P.C.B. Phillips
Title: "Structural Nonparametric Cointegrating Regression"
Journal: Econometrica
DOI: 10.3982/ECTA7732

Author(s): Q. Wang and P.C.B. Phillips
Title: "Nonparametric Cointegrating Regression with Endogeneity and Long Memory"
Journal: Econometric Theory
DOI: 10.1017/S0266466614000917

===== FUNCTIONS =====
The important functions are: KernelEst.m, NonparaCI.m, and TnStat.m (all located in the function folder).
 - KernelEst.m computes the nonparametric kernel estimate at a user-defined input point.
 - NonparaCI.m provides the lower and upper bound of the nonparametric estimate.
 - Wang and Phillips (2016) propose a parametric model specification test in their Section 3. The test statistic Tn is implemented in TnStat.m.
The Matlab implementation of these function is further illustrated with three examples.


===== OTHER FILES ====
 - Example1.m replicates Figure 4 of Wang and Phillips (2009).
 - A part of Table 1 of Wang and Phillips (2009) is replicated in Example2.m.
 - Example3.m implements a small Monte Carlo study. It uses the same cubic cointegrating relation as Example1.m and Example2.m and reports the empirical size of Wang and Phillips' specification test.



