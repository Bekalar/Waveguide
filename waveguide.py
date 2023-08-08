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


def bisecSiatka(xl, xr, tolx, fun, n_, h_):
    if fun(xr, n_, h_) * fun(xl, n_, h_) > 0:
        return 'no zero'
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


def fiSiatka(k, n_, h_):
    yp1 = math.sqrt(x0)
    yp2 = math.sqrt(x0 + h_)
    xx, psi = numerov_cowells2(x0, xE, n_, yp1, yp2, k2, k)
    return psi[-1]


def zeroPoints(h_val):
    K1 = 0.001 * Kmax
    K2 = K1 + dK
    zero_points = []

    n_val = round((xE - x0) / h_val)

    while K2 < Kmax:
        if fiSiatka(K1, n_val, h_val) * fiSiatka(K2, n_val, h_val) < 0:
            zero_points.append(bisecSiatka(K1, K2, tol, fiSiatka, n_val, h_val))
        K1 = K2
        K2 = K1 + dK

    return zero_points


def k2(r, k):
    return k * k + 1.0 / (4.0 * r * r)


# zad1
'''
Kmax = 20.0
dK = Kmax / 100.0
tol = Kmax / 10000.0
K1 = 0.001 * Kmax
K2 = K1 + dK
zero_points = []

while K2 < Kmax:
    if fiSiatka(K1, n, h) * fiSiatka(K2, n, h) < 0:
        zero_points.append(bisecSiatka(K1, K2, tol, fi, n, h))
    K1 = K2
    K2 = K1 + dK
'''
# siatka

h_lst = [0.1, 0.01, 0.001, 0.0001, 0.00001]
points_anal = [2.404826, 5.520078, 8.653728, 11.791534]
points_y = []
results = []
log_lst = []

for i in range(len(h_lst)):
    points_y.append(zeroPoints(h_lst[i]))
    log_lst.append(-1 * math.log10(h_lst[i]))
print(points_y)

for j in range(len(points_y[0])):
    temp = []
    for i in range(len(points_y)):
        temp.append(points_y[i][j])
    results.append(temp)

print(results)

for i in range(len(h_lst)):
    if i == 0:
        plt.plot(log_lst, results[i], label="numerical")
    else:
        plt.plot(log_lst, results[i])

for f in points_anal:
    if f == points_anal[0]:
        plt.axhline(f, color="black", label="analytical")
    else:
        plt.axhline(f, color="black")

plt.xlabel("h values")
plt.ylabel("Eigen values")
plt.legend(loc="upper right")

'''
# zad2

filename = "waveguide.data"
with open(filename, 'w+') as file:
    for points in zero_points:
        file.write(str(points))
        file.write('\n')

# zad3

for points in zero_points:
    print(points)
    x, y = fiMXY(points)
    for i in range(len(y)):
        y[i] = y[i] / math.sqrt(x[i])
    plt.plot(x, y, label="Funkcje radialne amplitud modów normalnych")
plt.xlabel("eigen values")
plt.ylabel("eigen vectors")
plt.legend()
'''
plt.grid()
plt.show()
