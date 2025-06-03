#include<iostream>
#include <fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include <iomanip>
int iter = 0;
int culc_func = 0;
const double alpha1 = 1.0;
const double alpha2 = 5.0;
const int Constobn = 4;
int indf = 3;
int h1 = 0;
int h2 = 0;
int indogr = 1;
double delta = 0.00001;
double r,betta;
double AimFunc(const double x, const double y)
{
    culc_func++;
    if (indogr == 1)
    {
        if (indf == 1)
        {
            return 10 * x * x - 4 * x * y + 7 * y * y - 4 * sqrt(5) * (5 * x - y) - 16 + r * (1 / (x + 10) + 1 / y + 1 / (1 - x - y));
        }
        if (indf == 2)
        {
            return 1 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1) + r * (1 / (x + 10) + 1 / y + 1 / (1 - x - y));
        }
        if (indf == 3)
        {
            return 5 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1) + r * (1 / (x + 10) + 1 / y + 1 / (1 - x - y));
        }
    }
    else
    {
        if (indf == 1)
        {
            return 10 * x * x - 4 * x * y + 7 * y * y - 4 * sqrt(5) * (5 * x - y) - 16 + r * (1 / (1 - (x + 3) * (x + 3) / 16 - (y + 4) * (y + 4) / 9));
        }
        if (indf == 2)
        {
            return 1 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1) + r * (1 / (1 - (x + 3) * (x + 3) / 16 - (y + 4) * (y + 4) / 9));
        }
        if (indf == 3)
        {
            return 5 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1) + r * (1 / (1 - (x + 3) * (x + 3) / 16 - (y + 4) * (y + 4) / 9));
        }
    }
}
double Funcx(const double x, const double y)
{
    return (AimFunc(x+delta,y)- AimFunc(x, y)) / delta;
}
double Funcy(const double x, const double y)
{
    return (AimFunc(x , y + delta)-AimFunc(x , y)) / delta;
}
double Funcxx(const double x, const double y)
{
    return (Funcx(x + delta, y) - Funcx(x, y)) / delta;
}
double Funcyy(const double x, const double y)
{
    return (Funcy(x, y + delta) - Funcy(x, y)) / delta;
}
double Funcxy(const double x, const double y)
{
    return (Funcy(x + delta, y) - Funcy(x, y)) / delta;
}
double NormVec(double x, double y)
{
    return sqrt(x * x + y * y);
}
double Pathx(double  fxx, double fxy, double fyy, double fx, double fy)
{
    return (fx * fyy - fy * fxy) / (fxx * fyy - fxy * fxy);
}
double Pathy(double  fxx, double fxy, double fyy, double fx, double fy)
{
    return (fxx * fy - fxy * fx) / (fxx * fyy - fxy * fxy);
}
std::vector<double> Newton(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda = 1;
    double xt, yt, fx, fy, fxx,fxy, fyy, px, py;
    double omegat;
    double etta = 0;
    double dzetta = 0;
    int indh = 0;
    xt = xn;
    yt = yn;
    fx = -Funcx(xt, yt);
    fy = -Funcy(xt, yt);
    points.push_back(xt);
    points.push_back(yt);
    betta = 0.5;
    r = 32;
    while (NormVec(fx, fy) > eps)
    {
        indh = 0;
        iter++;
        fxx = Funcxx(xt,yt);
        fxy = Funcxy(xt, yt);
        fyy = Funcyy(xt, yt);
        if ((fxx > 0) && ((fxx * fyy - fxy * fxy) > 0))
        {
            indh = 0;
            px = Pathx(fxx, fxy, fyy, fx, fy);
            py = Pathy(fxx, fxy, fyy, fx, fy);
        }
        else
        {
            indh = 1;
            while (indh==1)
            {
                fxx = fxx + eps;
                fyy = fyy + eps;
                if ((fxx > 0) && ((fxx * fyy - fxy * fxy) > 0))
                {
                    h1++;
                    indh = 0;
                }
            }
            px = Pathx(fxx, fxy, fyy, fx, fy);
            py = Pathy(fxx, fxy, fyy, fx, fy);
        }
        std::cout << fx << ' ' << fy << std::endl;
        std::cout << fxx<<' '<<fxy << std::endl;
        std::cout << fxy << ' ' << fyy << std::endl;
        std::cout << NormVec(px, py) << std::endl;
        xt = xt + px;
        yt = yt + py;
        fx = -Funcx(xt, yt);
        fy = -Funcy(xt, yt);
        points.push_back(xt);
        points.push_back(yt);
        r = r * betta;
    }
    points.push_back(xt);
    points.push_back(yt);

    return points;
}
int main()
{
    const double eps1 = 0.01;
    const double eps2 = 0.00001;
    std::vector<double> points1, points2;
    std::ofstream file;
    std::cout << std::setprecision(7);
    file << std::setprecision(7);
    indf = 3;
    file.open("Newtoneps1f3p3.txt");
    points1 = Newton(eps1, 0.2, -0.2);
    file << iter << ' ' << culc_func << std::endl;
    file << std::setprecision(7);
    for (int i = 0; i < points1.size() - 1; i += 2)
    {
        file << points1[i] << ' ' << points1[i + 1] << std::endl;
    }
    file.close();
    std::cout << std::endl;
    std::cout << h1 << ' ' << h2 << std::endl;
    std::cout << "Hello World!\n";
    return 0;
}


