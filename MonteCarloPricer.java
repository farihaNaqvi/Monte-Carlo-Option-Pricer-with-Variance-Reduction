import java.util.concurrent.ThreadLocalRandom;

public class MonteCarloPricer {

    // ============================================================
    // Black-Scholes closed-form PRICE (for benchmarking)
    // ============================================================
    static double blackScholesPrice(double S, double K, double T, double r,
                                    double sigma, double div, double t, String option) {
        double t2t = T - t;
        double dplus = (Math.log(S / K) + (r - div + 0.5 * sigma * sigma) * t2t)
                    / (sigma * Math.sqrt(t2t));
        double dminus = dplus - sigma * Math.sqrt(t2t);

        if (option.equals("Call")) {
            return S * Math.exp(-div * t2t) * normalCDF(dplus)
                - K * Math.exp(-r * t2t) * normalCDF(dminus);
        } else {
            return K * Math.exp(-r * t2t) * normalCDF(-dminus)
                - S * Math.exp(-div * t2t) * normalCDF(-dplus);
        }
    }

    // ============================================================
    // Black-Scholes DELTA (this is what the control variate needs)
    // ============================================================
    static double blackScholesDelta(double S, double K, double T, double r,
                                    double sigma, double div, double t, String option) {
        double t2t = T - t;
        double dplus = (Math.log(S / K) + (r - div + 0.5 * sigma * sigma) * t2t)
                    / (sigma * Math.sqrt(t2t));

        if (option.equals("Call")) {
            return Math.exp(-div * t2t) * normalCDF(dplus);
        } else {
            return Math.exp(-div * t2t) * (normalCDF(dplus) - 1.0);
        }
    }

    // ============================================================
    // Standard normal CDF (Abramowitz & Stegun approximation)
    // ============================================================
    static double normalCDF(double x) {
        double a1 =  0.254829592;
        double a2 = -0.284496736;
        double a3 =  1.421413741;
        double a4 = -1.453152027;
        double a5 =  1.061405429;
        double p  =  0.3275911;

        int sign = (x < 0) ? -1 : 1;
        x = Math.abs(x) / Math.sqrt(2);

        double t = 1.0 / (1.0 + p * x);
        double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1)
                         * t * Math.exp(-x * x);

        return 0.5 * (1.0 + sign * y);
    }

    // ============================================================
    // MAIN — Monte Carlo with Antithetic + Delta Control Variate
    // ============================================================
    public static void main(String[] args) {
        // Parameters
        double S = 100;
        double r = 0.06;
        double sig = 0.2;
        double T = 1;
        double K = 100;
        int N = 300;          // time steps (scaled up)
        int M = 100_000;      // number of simulations (scaled up)
        double div = 0.03;
        String option = "Put";

        // Precompute constants
        double dt = T / N;
        double nu = r - div - 0.5 * sig * sig;
        double nudt = nu * dt;
        double sigsdt = sig * Math.sqrt(dt);
        double erddt = Math.exp((r - div) * dt);
        double beta1 = -1;

        double sumCT = 0;
        double sumCT2 = 0;

        long startTime = System.currentTimeMillis();

        for (int j = 0; j < M; j++) {                // FIXED: 0 to M (was 1 to M)
            double St1 = S;
            double St2 = S;
            double cv1 = 0;
            double cv2 = 0;

            for (int i = 0; i < N; i++) {            // FIXED: 0 to N (was 1 to N)
                double t = i * dt;                    // FIXED: was (i-1)*dt

                double delta1 = blackScholesDelta(St1, K, T, r, sig, div, t, option);
                double delta2 = blackScholesDelta(St2, K, T, r, sig, div, t, option);

                double eps = ThreadLocalRandom.current().nextGaussian();

                double Stn1 = St1 * Math.exp(nudt + sigsdt * eps);
                double Stn2 = St2 * Math.exp(nudt + sigsdt * (-eps));

                cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);

                St1 = Stn1;
                St2 = Stn2;
            }

            double CT;
            if (option.equals("Call")) {
                CT = 0.5 * (Math.max(0, St1 - K) + beta1 * cv1
                        + Math.max(0, St2 - K) + beta1 * cv2);
            } else {
                CT = 0.5 * (Math.max(0, K - St1) + beta1 * cv1
                        + Math.max(0, K - St2) + beta1 * cv2);
            }

            sumCT += CT;
            sumCT2 += CT * CT;
        }

        long elapsed = System.currentTimeMillis() - startTime;

        double value = (sumCT / M) * Math.exp(-r * T);
        double sd = Math.sqrt((sumCT2 - sumCT * sumCT / M) * Math.exp(-2 * r * T) / (M - 1));
        double se = sd / Math.sqrt(M);

        double bsmPrice = blackScholesPrice(S, K, T, r, sig, div, 0, option);
        double errorPct = Math.abs(value - bsmPrice) / bsmPrice * 100;

        System.out.println("=== Monte Carlo (Antithetic + Delta Control Variate) ===");
        System.out.printf("Option type:          %s%n", option);
        System.out.printf("Paths:                %,d%n", M);
        System.out.printf("Time steps:           %,d%n", N);
        System.out.printf("MC Value:             %.6f%n", value);
        System.out.printf("Black-Scholes Value:  %.6f%n", bsmPrice);
        System.out.printf("Error vs BSM:         %.4f%%%n", errorPct);
        System.out.printf("Standard Deviation:   %.6f%n", sd);
        System.out.printf("Standard Error:       %.6f%n", se);
        System.out.printf("Compute time:         %,d ms%n", elapsed);
    }
}