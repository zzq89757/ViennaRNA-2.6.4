#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#define BIG 1.44115188075855872e+17
#define MAXLOG 7.09782712893383996843e+02
#define MACHEP 1.11022302462515654042e-16
#define MAXNUM 1.7976931348623158e+308
#define EUL 0.57721566490153286060
#define PI 3.14159265358979323846

double expn(int n, double x) {
    if (n < 0 || x < 0) {
        fprintf(stderr, "Invalid input: n and x must be non-negative\n");
        exit(EXIT_FAILURE);
    }

    if (x > MAXLOG) {
        return 0.0;
    }

    if (x == 0.0) {
        if (n < 2) {
            fprintf(stderr, "Singularity: E_n(x) is undefined for n < 2\n");
            exit(EXIT_FAILURE);
        } else {
            return 1.0 / (n - 1.0);
        }
    }

    if (n == 0) {
        return exp(-x) / x;
    }

    if (n > 5000) {
        double xk = x + n;
        double yk = 1.0 / (xk * xk);
        double t = n;
        double ans = yk * t * (6.0 * x * x - 8.0 * t * x + t * t);
        ans = yk * (ans + t * (t - 2.0 * x));
        ans = yk * (ans + t);
        ans = (ans + 1.0) * exp(-x) / xk;
        return ans;
    }

    if (x > 1.0) {
        double psi = -EUL - log(x);
        for (int i = 1; i < n; ++i) {
            psi += 1.0 / i;
        }

        double z = -x;
        double xk = 0.0;
        double yk = 1.0;
        double pk = 1.0 - n;
        double ans = n != 1 ? 1.0 / pk : 0.0;

        while (1) {
            xk += 1.0;
            yk *= z / xk;
            pk += 1.0;
            if (pk != 0.0) {
                ans += yk / pk;
            }
            double t = fabs(yk / ans);
            if (t <= MACHEP) {
                break;
            }
        }

        double t = (double)n;
        double r = n - 1;
        ans = (pow(z, r) * psi / tgamma(t)) - ans;
        return ans;
    }

    int k = 1;
    double pkm2 = 1.0;
    double qkm2 = x;
    double pkm1 = 1.0;
    double qkm1 = x + n;
    double ans = pkm1 / qkm1;

    while (1) {
        ++k;
        double yk = (k & 1) ? 1.0 : x;
        double xk = (k & 1) ? n + (k - 1) / 2 : k / 2;

        double pk = pkm1 * yk + pkm2 * xk;
        double qk = qkm1 * yk + qkm2 * xk;

        if (qk != 0) {
            double r = pk / qk;
            double t = fabs((ans - r) / r);
            ans = r;
            if (t <= MACHEP) {
                break;
            }
        } else {
            if (fabs(pk) > BIG) {
                pkm2 /= BIG;
                pkm1 /= BIG;
                qkm2 /= BIG;
                qkm1 /= BIG;
            }
        }
    }

    ans *= exp(-x);
    return ans;
}

double kn(int nn, double x) {
    int n = abs(nn);

    if (x <= 0.0) {
        if (x < 0.0) {
            fprintf(stderr, "kn: Domain error\n");
            exit(EXIT_FAILURE);
        } else {
            fprintf(stderr, "kn: Singularity\n");
            exit(EXIT_FAILURE);
        }
    }

    if (x > 9.55) {
        // Asymptotic expansion for large x
        double k = (double)n;
        double pn = 4.0 * k * k;
        double pk = 1.0;
        double z0 = 8.0 * x;
        double fn = 1.0;
        double t = 1.0;
        double s = t;
        double nkf = MAXNUM;
        int i = 0;

        while (1) {
            double z = pn - pk * pk;
            t *= z / (fn * z0);
            double nk1f = fabs(t);
            if (i >= n && nk1f > nkf) {
                break;
            }
            nkf = nk1f;
            s += t;
            fn += 1.0;
            pk += 2.0;
            ++i;
            if (fabs(t / s) <= MACHEP) {
                break;
            }
        }

        double ans = exp(-x) * sqrt(PI / (2.0 * x)) * s;
        return ans;
    }

    double ans = 0.0;
    double z0 = 0.25 * x * x;
    double fn = 1.0;
    double pn = 0.0;
    double zmn = 1.0;
    double tox = 2.0 / x;

    if (n > 0) {
        pn = -EUL;
        double k = 1.0;
        for (int i = 1; i < n; ++i) {
            pn += 1.0 / k;
            k += 1.0;
            fn *= k;
        }

        zmn = tox;

        if (n == 1) {
            ans = 1.0 / x;
        } else {
            double nk1f = fn / n;
            double kf = 1.0;
            double s = nk1f;
            double z = -z0;
            double zn = 1.0;
            for (int i = 1; i < n; ++i) {
                nk1f = nk1f / (n - i);
                kf *= i;
                zn *= z;
                double t = nk1f * zn / kf;
                s += t;
                if ((MAXNUM - fabs(t)) < fabs(s)) {
                    fprintf(stderr, "kn: Overflow\n");
                    exit(EXIT_FAILURE);
                }
                if (tox > 1.0 && (MAXNUM / tox) < zmn) {
                    fprintf(stderr, "kn: Overflow\n");
                    exit(EXIT_FAILURE);
                }
                zmn *= tox;
            }

            s *= 0.5;
            double t = fabs(s);
            if (zmn > 1.0 && (MAXNUM / zmn) < t) {
                fprintf(stderr, "kn: Overflow\n");
                exit(EXIT_FAILURE);
            }
            if (t > 1.0 && (MAXNUM / t) < zmn) {
                fprintf(stderr, "kn: Overflow\n");
                exit(EXIT_FAILURE);
            }
            ans = s * zmn;
        }
    }

    double tlg = 2.0 * log(0.5 * x);
    double pk = -EUL;
    double t;
    if (n == 0) {
        pn = pk;
        t = 1.0;
    } else {
        pn += 1.0 / n;
        t = 1.0 / fn;
    }

    double s = (pk + pn - tlg) * t;
    double k = 1.0;
    while (1) {
        t *= z0 / (k * (k + n));
        pk += 1.0 / k;
        pn += 1.0 / (k + n);
        s += (pk + pn - tlg) * t;
        k += 1.0;
        if (fabs(t / s) <= MACHEP) {
            break;
        }
    }

    s = 0.5 * s / zmn;
    if (n % 2 != 0) {
        s = -s;
    }
    ans += s;
    return ans;
}

typedef struct {
    int i;
    int j;
    double energy;
    char structure[100];  // Adjust size as needed
} MFE;

const int BP_REV_DEFAULT[] = { 0, 2, 1, 4, 3, 6, 5, 7 };
const int BP_ALIAS_DEFAULT[] = { 0, 1, 2, 3, 4, 3, 2, 0 };
const int BP_ENCODING_DEFAULT[8][8] = {
    {0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 5, 0, 0, 5},
    {0, 0, 0, 1, 0, 0, 0, 0},
    {0, 0, 2, 0, 3, 0, 0, 0},
    {0, 6, 0, 4, 0, 0, 0, 6},
    {0, 0, 0, 0, 0, 0, 2, 0},
    {0, 0, 0, 0, 0, 1, 0, 0},
    {0, 6, 0, 0, 5, 0, 0, 0}};



