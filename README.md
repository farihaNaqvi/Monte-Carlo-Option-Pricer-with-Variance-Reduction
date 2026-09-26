# Monte Carlo Option Pricer with Variance Reduction

A Java implementation of a Monte Carlo option pricer that uses antithetic variates and a delta control variate to reduce variance and improve pricing accuracy for European call and put options.

## Overview

This project prices options by simulating many possible stock price paths under the Black-Scholes framework and comparing the Monte Carlo estimate against the analytical Black-Scholes benchmark.

The program includes:

- Monte Carlo simulation with antithetic sampling
- Delta control variate for variance reduction
- Black-Scholes closed-form price and delta formulas
- Normal CDF approximation using Abramowitz & Stegun
- Console output with Monte Carlo value, standard error, and percentage error versus Black-Scholes

## Files

- `MonteCarloPricer.java` — main implementation
- `MonteCarloPricer.class` — compiled Java class (generated output)
- `.gitignore` — standard Git ignore file

## Core pricing approach

The model simulates geometric Brownian motion for the underlying asset:

- `S` = current stock price
- `K` = strike price
- `T` = maturity
- `r` = risk-free rate
- `sig` = volatility
- `div` = dividend yield
- `N` = number of time steps
- `M` = number of Monte Carlo paths

The program uses:

- antithetic paths to reduce variance
- delta-based control variate to improve estimator efficiency
- a Black-Scholes benchmark for validation

## How to run

Requirements:

- Java JDK 8 or later

From the repository root, compile and run:

```bash
javac MonteCarloPricer.java
java MonteCarloPricer
```

## Configuration

The default parameters are defined in `main(String[] args)`:

```java
double S = 100;
double r = 0.06;
double sig = 0.2;
double T = 1;
double K = 100;
int N = 300;
int M = 100_000;
double div = 0.03;
String option = "Put";
```

You can modify these values to price different scenarios or switch to a call option by changing:

```java
String option = "Call";
```

## Example output

```text
=== Monte Carlo (Antithetic + Delta Control Variate) ===
Option type:          Put
Paths:                100,000
Time steps:           300
MC Value:             8.123456
Black-Scholes Value:  8.100000
Error vs BSM:         0.2876%
Standard Deviation:   0.123456
Standard Error:       0.000389
Compute time:         245 ms
```

## Notes

This project is intended as a numerical finance and variance reduction exercise. It demonstrates how control variates and antithetic sampling can improve the quality of Monte Carlo pricing estimates without needing a more advanced library.

## License

This repository does not currently include a license file. If you plan to distribute or reuse the code publicly, consider adding an open-source license such as MIT or Apache 2.0.
