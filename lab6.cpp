#include<iostream>
#include <fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include <iomanip>
#include<string>
int iter = 0;
int culc_func = 0;
const double alpha1 = 1.0;
const double alpha2 = 5.0;
int Constobn = 2;
int indf = 3;
int indogr = 1;
int h1 = 0;
int h2 = 0;
int count = 0;
double p, l;
double betta = 0.5;
double r;
double delta = 0.00001;
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
            return alpha1 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1) + r * (1 / (x + 10) + 1 / y + 1 / (1 - x - y));
        }
        if (indf == 3)
        {
            return alpha2 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1) + r * (1 / (x + 10) + 1 / y + 1 / (1 - x - y));
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
            return alpha1 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1) + r * (1 / (1 - (x + 3) * (x + 3) / 16 - (y + 4) * (y + 4) / 9));
        }
        if (indf == 3)
        {
            return alpha2 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1) + r * (1 / (1 - (x + 3) * (x + 3) / 16 - (y + 4) * (y + 4) / 9));
        }
    }
}
double Funcx(double x, double y)
{
    return (AimFunc(x + delta, y) - AimFunc(x, y)) / delta;
}
double Funcy(double x, double y)
{
    return (AimFunc(x, y + delta) - AimFunc(x, y)) / delta;
}
int Crit(double x, double y)
{
    if (indogr == 1)
    {
        if ((x >= -10) && (y >= 0) && ((x + y) <= 1))
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
    else
    {
        if (((x + 3) * (x + 3)) / 16 + ((y + 4) * (y + 4)) / 9 <= 1)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
}
double FuncP(double x, double y)
{
    if (indogr == 1)
    {
        return r * (1 / (x + 10) + 1 / y + 1 / (1 - x - y));
    }
    else
    {
        return r * (1 / (1 - (x + 3) * (x + 3) / 16 - (y + 4) * (y + 4) / 9));
    }
    
}
double AimFunclymbda(const double lymbda, const double pointx, const double pointy, const double fx, const double fy)
{
    double x, y;
    x = pointx + lymbda * fx;
    y = pointy + lymbda * fy;
    return AimFunc(x, y);
}
double NormVec(double x, double y)
{
    return sqrt(x * x + y * y);
}
double NormVecGramm(double x, double y, double xt, double yt)
{
    return sqrt(NormVec(Funcx(xt, yt), Funcx(xt, yt)) * x * x + y * y);
}
double Gold(const double eps, const double pointx, const double pointy, const double fx, const double fy)
{
    const double tau1 = 2 / (sqrt(5) + 3);
    const double tau2 = 2 / (sqrt(5) + 1);
    double x1, x2, f1, f2;
    double ak = 0;
    double bk = p/ (1 / NormVec(fx, fy));
    int count = 0;
    x1 = ak + tau1 * (bk - ak);
    x2 = ak + tau2 * (bk - ak);
    f1 = AimFunclymbda(x1, pointx, pointy, fx, fy);
    f2 = AimFunclymbda(x2, pointx, pointy, fx, fy);
    while (bk - ak > eps / 1000)
    {
        count++;
        if (f1 < f2)
        {
            bk = x2;
            x2 = x1;
            f2 = f1;
            x1 = ak + tau1 * (bk - ak);
            f1 = AimFunclymbda(x1, pointx, pointy, fx, fy);
        }
        else
        {
            ak = x1;
            x1 = x2;
            f1 = f2;
            x2 = ak + tau2 * (bk - ak);
            f2 = AimFunclymbda(x2, pointx, pointy, fx, fy);
        }
    }
    return (bk + ak) / 2;
}
std::vector<double> SoprgradPR(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda = 0;
    double xt, yt, fx, fy, fxl, fyl, px, py,xtl,ytl;
    double omegat;
    double gamma = 0;
    int count = 0;
    xtl = xn;
    ytl = yn;
    xt = xtl;
    yt = ytl;
    fx = -Funcx(xt, yt);
    fy = -Funcy(xt, yt);
    px = fx;
    py = fy;
    points.push_back(xt);
    points.push_back(yt);
    omegat = NormVec(px, py);
    while (abs(FuncP(xt,yt)) > eps)
    {
        iter++;
        count++;
        xtl = xt;
        ytl = yt;
        lymbda = Gold(eps, xtl, ytl, px, py);
        xt = xtl + lymbda * px;
        yt = ytl + lymbda * py;
        if (Crit(xt, yt) == 1)
        {
            xt = xtl;
            yt = ytl;
            r = betta * r;
            p = p / 2;
        }
        std::cout << xt << ' ' << yt <<' ' <<r<< std::endl;
        points.push_back(xt);
        points.push_back(yt);
        fxl = fx;
        fyl = fy;
        fx = -Funcx(xt, yt);
        fy = -Funcy(xt, yt);
        gamma = ((-fx) * (fxl - fx) + (-fy) * (fyl - fy)) / (fxl * fxl + fyl * fyl);
        if (count == Constobn)
        {
            gamma = 0;
            count = 0;
        }
        px = gamma * px + fx;
        py = gamma * py + fy;
    }
    return points;
}


int main()
{
    const double eps[2] = { 0.01,0.00001 };
    std::vector<double> points1, points2;
    std::ofstream file;
    std::string fil, end;
    double x0[2] = { -5,-6 };
    double y0[2] = { -5,2 };
    p =1;
    l =1;
    r = 100;
    //Constobn = 4;  // Константа обновления
    indf = 1;   //функция
    int i = 1;   //точность
    int j = 2; //точка
    fil = "PReps";
    end = std::to_string(i);
    end += "f";
    end += std::to_string(indf);
    end += "p";
    end += std::to_string(j);
    end += ".txt";
    std::cout << std::setprecision(7);
    file << std::setprecision(7);
    file.open(fil + end);
    culc_func = 0;
    iter = 0;
    fil = "PReps";
    file.open(fil + end);
    points2 = SoprgradPR(eps[i - 1], x0[j - 1], y0[j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    for (int i = 0; i < points2.size(); i += 2)
    {
        file << points2[i] << ' ' << points2[i + 1] << std::endl;
    }
    file.close();
    std::cout << "PR:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    std::cout << "Hello World!\n";
   
    
    std::cout << "Hello World!\n";
    return 0;
}



