#include <iostream>
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
int Constobn;
int indf = 0;
double p, l;
double AimFunc(const double x, const double y)
{
    culc_func++;
    if (indf == 1)
    {
        return 10 * x * x - 4 * x * y + 7 * y * y - 4 * sqrt(5) * (5 * x - y) - 16;
    }
    if (indf == 2)
    {
        return alpha1 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1);
    }
    if (indf == 3)
    {
        return alpha2 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1);
    }
}
double Funcx(const double x, const double y)
{
    culc_func++;
    if (indf == 1)
    {
        return 20 * x - 4 * y - 20 * sqrt(5);
    }
    if (indf == 2)
    {
        return -2 + 2 * x + 4 * x * x * x * alpha1 - 4 * x * y * alpha1;
    }
    if (indf == 3)
    {
        return -2 + 2 * x + 4 * x * x * x * alpha2 - 4 * x * y * alpha2;
    }
}
double Funcy(const double x, const double y)
{
    culc_func++;
    if (indf == 1)
    {
        return  4 * sqrt(5) - 4 * x + 14 * y;
    }
    if (indf == 2)
    {
        return -2 * x * x * alpha1 + 2 * y * alpha1;
    }
    if (indf == 3)
    {
        return -2 * x * x * alpha2 + 2 * y * alpha2;
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
double NormVecGramm(double x, double y,double xt,double yt)
{
    return sqrt(NormVec(Funcx(xt,yt), Funcx(xt, yt))*x * x + y * y);
}
double Gold(const double eps, const double pointx, const double pointy,const double fx,const double fy)
{
    const double tau1 = 2 / (sqrt(5) + 3);
    const double tau2 = 2 / (sqrt(5) + 1);
    double x1, x2, f1, f2;
    double ak = 0;
    double bk=1;
    if (iter <= 3)
    {
        bk = p * (1 / NormVec(fx, fy));
    }
    else
    {
        bk =l * (1 / NormVec(fx, fy));
    }
    int count = 0;
    x1 = ak + tau1 * (bk - ak);
    x2 = ak + tau2 * (bk - ak);
    f1 = AimFunclymbda(x1, pointx, pointy, fx, fy);
    f2 = AimFunclymbda(x2, pointx, pointy, fx, fy);
    while (bk - ak > eps/1000)
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
std::vector<double> Soprgrad(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda = 0;
    double xt, yt, fx, fy, fxl, fyl,px,py;
    double omegat;
    double gamma = 0;
    int count=0;
    xt = xn;
    yt = yn;
    fx = -Funcx(xt, yt);
    fy = -Funcy(xt, yt);
    px = fx;
    py = fy;
    points.push_back(xt);
    points.push_back(yt);
    omegat = NormVec(px, py);
    while (omegat > eps)
    {
        iter++;
        count++;
        lymbda = Gold(eps, xt, yt, px, py);
        xt = xt + lymbda * px;
        yt = yt + lymbda * py;
        points.push_back(xt);
        points.push_back(yt);
        fxl = fx;
        fyl = fy;
        fx = -Funcx(xt, yt);
        fy = -Funcy(xt, yt);
        gamma = (fx * fx + fy * fy) / (fxl * px + fyl  * py);
        if (count == Constobn)
        {
            gamma = 0;
            count = 0;
        }
        px = gamma * px + fx;
        py = gamma * py + fy;
        omegat = NormVec(fx, fy);

    }
  return points;
}
std::vector<double> SoprgradFR(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda = 0;
    double xt, yt, fx, fy, fxl, fyl, px, py;
    double omegat;
    double gamma = 0;
    int count = 0;
    xt = xn;
    yt = yn;
    fx = -Funcx(xt, yt);
    fy = -Funcy(xt, yt);
    px = fx;
    py = fy;
    points.push_back(xt);
    points.push_back(yt);
    omegat = NormVec(px, py);
    while (omegat > eps)
    {
        iter++;
        count++;
        lymbda = Gold(eps, xt, yt, px, py);
        xt = xt + lymbda * px;
        yt = yt + lymbda * py;
        points.push_back(xt);
        points.push_back(yt);
        fxl = fx;
        fyl = fy;
        fx = -Funcx(xt, yt);
        fy = -Funcy(xt, yt);
        gamma = (fx*fx+fy*fy) / (fxl * fxl + fyl * fyl);
        if (count == Constobn)
        {
            gamma = 0;
            count = 0;
        }
        px = gamma * px + fx;
        py = gamma * py + fy;
        omegat = NormVec(fx, fy);

    }
   return points;
}
std::vector<double> SoprgradPR(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda = 0;
    double xt, yt, fx, fy, fxl, fyl, px, py;
    double omegat;
    double gamma = 0;
    int count = 0;
    xt = xn;
    yt = yn;
    fx = -Funcx(xt, yt);
    fy = -Funcy(xt, yt);
    px = fx;
    py = fy;
    points.push_back(xt);
    points.push_back(yt);
    omegat = NormVec(px, py);
    while (omegat > eps)
    {
        iter++;
        count++;
        lymbda = Gold(eps, xt, yt, px, py);
        xt = xt + lymbda * px;
        yt = yt + lymbda * py;
        points.push_back(xt);
        points.push_back(yt);
        fxl = fx;
        fyl = fy;
        fx = -Funcx(xt, yt);
        fy = -Funcy(xt, yt);
        gamma = (( - fx) * (fxl - fx) + ( - fy) * (fyl - fy)) / (fxl * fxl + fyl * fyl);
        if (count == Constobn)
        {
            gamma = 0;
            count = 0;
        }
        px = gamma * px + fx;
        py = gamma * py + fy;
        omegat = NormVec(fx, fy);

    }
   return points;
}
int main()
{
    const double eps[2] = {0.01,0.00001};
    std::vector<double> points1, points2;
    std::ofstream file;
    std::string fil,end;
    double x0[2] = { -2,1.5 };
    double y0[2] = {5,2};
    p = 3;
    l = 2;
    Constobn = 1000;  // Константа обновления
    indf = 3;   //функция
    int i=2;   //точность
    int j =1; //точка
    fil = "Soprgradeps";
    end = std::to_string(i);
    end += "f";
    end += std::to_string(indf);
    end += "p";
    end += std::to_string(j);
    end += ".txt";
    std::cout << std::setprecision(7);
    file << std::setprecision(7);
    file.open(fil+end);
    points1 = Soprgrad(eps[i-1], x0[j-1], y0[j-1]);
    file << iter << ' ' << culc_func << std::endl;
    file << std::setprecision(7);
    for (int i = 0; i < points1.size() - 1; i += 2)
    {
        file << points1[i] << ' ' << points1[i + 1] << std::endl;
    }
    file.close();
    std::cout <<"Just sopr:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    culc_func = 0;
    iter = 0;
    fil = "SoprgradFReps";
    file.open(fil + end);
    points2 = SoprgradFR(eps[i - 1], x0[j - 1], y0[j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    for (int i = 0; i < points2.size(); i += 2)
    {
        file << points2[i] << ' ' << points2[i + 1] << std::endl;
    }
    file.close();
    std::cout << "FR sopr:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    culc_func = 0;
    iter = 0;
    fil = "SoprgradPReps";
    file.open(fil + end);
    points2 = SoprgradPR(eps[i - 1], x0[j - 1], y0[j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    for (int i = 0; i < points2.size(); i += 2)
    {
        file << points2[i] << ' ' << points2[i + 1] << std::endl;
    }
    file.close();
    std::cout << "PR sopr:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    std::cout << "Hello World!\n";
    return 0;
}

