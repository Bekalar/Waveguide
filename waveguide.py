import math
import matplotlib.pyplot as plt

x0 = 0.001
xE = 1.0
Kmax = 20.0
dK = Kmax / 100.0
tol = Kmax / 10000.0
n = 10000
h = (xE - x0) / n


def numerov_cowells2(a, b, ndiv, y0, y1, funk, kk):
    hh = (b - a) / ndiv
    xp = [a, a + hh]
    yp = [y0, y1]
    for l in range(1, ndiv):
        xp.append(xp[l] + hh)
        yp.append(
            (2. * yp[l] * (1. - 5.0 * hh * hh * funk(xp[l], kk) / 12.0) / (1. + hh * hh * funk(xp[l + 1], kk) / 12.0)
             - yp[l - 1] * (1. + hh * hh * funk(xp[l - 1], kk) / 12.0) / (
                     1. + hh * hh * funk(xp[l + 1], kk) / 12.0)))
    return xp, yp


def bisec(xl, xr, tolx, fun, n_, h_):
    if fun(xr, n_, h_) * fun(xl, n_, h_) > 0:
        return "no zero"
    while (xr - xl) / 2.0 > tolx:
        xs = (xr + xl) / 2.0
        if fun(xr, n_, h_) * fun(xs, n_, h_) < 0:
            xl = xs
        else:
            xr = xs
    return (xr + xl) / 2.0


def fiMXY(k):
    yp1 = math.sqrt(x0)
    yp2 = math.sqrt(x0 + h)
    xx, psi = numerov_cowells2(x0, xE, n, yp1, yp2, k2, k)
    return xx, psi


def fiNet(k, n_, h_):
    yp1 = math.sqrt(x0)
    yp2 = math.sqrt(x0 + h_)
    xx, psi = numerov_cowells2(x0, xE, n_, yp1, yp2, k2, k)
    return psi[-1]


def k2(r, k):
    return k * k + 1.0 / (4.0 * r * r)



def main():    
    Kmax = 20.0
    dK = Kmax / 100.0
    tol = Kmax / 10000.0
    K1 = 0.001 * Kmax
    K2 = K1 + dK
    zero_points = []

    while K2 < Kmax:
        if fiNet(K1, n, h) * fiNet(K2, n, h) < 0:
            zero_points.append(bisec(K1, K2, tol, fiNet, n, h))
        K1 = K2
        K2 = K1 + dK

    filename = "waveguide.data"
    with open(filename, "w+") as file:
        for points in zero_points:
            file.write(str(points))
            file.write("\n")

    for points in zero_points:
        # print(points)
        x, y = fiMXY(points)
        for i in range(len(y)):
            y[i] = y[i] / math.sqrt(x[i])
        plt.plot(x, y)
    plt.xlabel("eigen values")
    plt.ylabel("eigen vectors")
    plt.legend()

    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()